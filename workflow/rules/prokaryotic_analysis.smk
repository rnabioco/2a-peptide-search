"""
Prokaryotic 2A-like discovery - Analysis, Validation, and Reporting

Validates discovered motifs against known stalling peptides and generates
comprehensive reports comparing different discovery approaches.
"""

import glob as pyglob
import re


def get_bacteria_cluster_hmms(wildcards):
    """Discover cluster HMMs from checkpoint output for bacteria database."""
    checkpoint_output = checkpoints.identify_consensus_patterns.get(
        database="bacteria"
    ).output.alignments_dir
    sto_files = pyglob.glob(f"{checkpoint_output}/cluster_*.sto")
    cluster_ids = sorted(
        [int(re.search(r"cluster_(\d+)", f).group(1)) for f in sto_files]
    )
    return expand(
        RESULTS_DIR + "/prokaryotic/gp_analysis/bacteria/models/cluster_{cluster_id}.hmm",
        cluster_id=cluster_ids,
    )


rule validate_against_known_peptides:
    """Compare discovered motifs to known stalling peptides (SecM, TnaC, etc.)."""
    input:
        hmms=get_bacteria_cluster_hmms,
        known_peptides="resources/stalling-peptides/known_stalling_peptides.fasta",
    output:
        validation=RESULTS_DIR + "/prokaryotic/validation/known_peptide_hits.tsv",
        summary=RESULTS_DIR + "/prokaryotic/validation/validation_summary.txt",
    params:
        hmms_pattern=RESULTS_DIR + "/prokaryotic/gp_analysis/bacteria/models/cluster_*.hmm",
    log:
        LOGS_DIR + "/prokaryotic/validate_known_peptides.log",
    threads: 4
    shell:
        """
        python workflow/scripts/validate_stalling_peptides.py \
            --hmms '{params.hmms_pattern}' \
            --known-peptides "{input.known_peptides}" \
            --validation "{output.validation}" \
            --summary "{output.summary}" \
            --threads {threads} \
            > "{log}" 2>&1
        """


# ============================================================================
# Comparative Analysis: All GP vs Inter-domain GP
# ============================================================================


rule compare_approaches:
    """Compare results from all-GP vs inter-domain-GP vs arrest-motif approaches."""
    input:
        all_gp_motifs=RESULTS_DIR + "/prokaryotic/gp_analysis/bacteria/all_gp_motifs.tsv.gz",
        interdomain_motifs=RESULTS_DIR
        + "/prokaryotic/gp_analysis/bacteria/interdomain_gp_motifs.tsv.gz",
        clusters=RESULTS_DIR + "/prokaryotic/gp_analysis/bacteria/gp_clusters.tsv.gz",
        gp_validation=RESULTS_DIR + "/prokaryotic/validation/known_peptide_hits.tsv",
        # APPROACH 1: Seed-based searches (merged from all databases)
        seed_comprehensive=RESULTS_DIR
        + "/prokaryotic/seed_search_results/comprehensive_merged.sto.gz",
        # Arrest motif approach results
        arrest_motifs=RESULTS_DIR + "/prokaryotic/arrest_analysis/bacteria/arrest_motifs.tsv.gz",
        arrest_validation=RESULTS_DIR + "/prokaryotic/arrest_analysis/bacteria/known_peptide_validation.tsv",
        arrest_summary=RESULTS_DIR + "/prokaryotic/arrest_analysis/bacteria/motif_summary.tsv",
    output:
        comparison=RESULTS_DIR + "/prokaryotic/analysis/approach_comparison.tsv",
        plots=directory(RESULTS_DIR + "/prokaryotic/analysis/comparison_plots/"),
    log:
        LOGS_DIR + "/prokaryotic/compare_approaches.log",
    shell:
        """
        python workflow/scripts/compare_gp_approaches.py \
            --all-gp-motifs "{input.all_gp_motifs}" \
            --interdomain-motifs "{input.interdomain_motifs}" \
            --clusters "{input.clusters}" \
            --gp-validation "{input.gp_validation}" \
            --arrest-motifs "{input.arrest_motifs}" \
            --arrest-validation "{input.arrest_validation}" \
            --arrest-summary "{input.arrest_summary}" \
            --comparison "{output.comparison}" \
            --plots "{output.plots}" \
            > "{log}" 2>&1
        """


# ============================================================================
# Final Report
# ============================================================================


rule prokaryotic_discovery_report:
    """Generate comprehensive report on prokaryotic 2A-like discovery."""
    input:
        comparison=RESULTS_DIR + "/prokaryotic/analysis/approach_comparison.tsv",
        validation=RESULTS_DIR + "/prokaryotic/validation/validation_summary.txt",
        conservation=RESULTS_DIR + "/prokaryotic/gp_analysis/bacteria/cluster_conservation.tsv.gz",
        consensus=RESULTS_DIR + "/prokaryotic/gp_analysis/bacteria/consensus_patterns.tsv",
        # Arrest motif analysis inputs
        arrest_motif_summary=RESULTS_DIR + "/prokaryotic/arrest_analysis/bacteria/motif_summary.tsv",
        arrest_validation=RESULTS_DIR + "/prokaryotic/arrest_analysis/bacteria/known_peptide_validation.tsv",
        arrest_known_summary=RESULTS_DIR + "/prokaryotic/arrest_analysis/bacteria/known_peptide_summary.txt",
    output:
        report=RESULTS_DIR + "/prokaryotic/reports/prokaryotic_discovery.html",
    log:
        LOGS_DIR + "/prokaryotic/generate_report.log",
    script:
        "../scripts/prokaryotic_discovery_report.qmd"
