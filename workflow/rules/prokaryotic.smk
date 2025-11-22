"""
Rules for discovering prokaryotic 2A-like ribosomal stalling peptides.

Strategy:
1. Extract all GP-containing sequences from prokaryotic proteomes
2. Annotate with domain boundaries (Pfam/InterPro)
3. Focus on inter-domain GP motifs
4. Cluster by sequence context
5. Build initial consensus and HMMs
6. Iterative refinement

Known prokaryotic stalling peptides:
- SecM: FXXXXWIXXXXGIRAGP (stalls at RAGP)
- TnaC: WDPXXXX...GP-like
- MifM: Contains stalling sequence
"""

# ============================================================================
# Data Preparation
# ============================================================================

rule download_prokaryotic_proteomes:
    """Download prokaryotic reference proteomes."""
    output:
        fasta="data/prokaryotic/uniprot_bacteria.fasta.gz"
    params:
        url=config["prokaryotic_databases"]["bacteria"]["url"]
    log:
        "logs/download/prokaryotic_proteomes.log"
    shell:
        """
        mkdir -p data/prokaryotic
        wget -c -o {log} {params.url} -O {output.fasta}
        """


rule download_domain_annotations:
    """Download Pfam or pre-computed domain annotations."""
    output:
        annotations="data/prokaryotic/domain_annotations.tsv.gz"
    params:
        url=config["prokaryotic_databases"]["pfam"]["url"]
    log:
        "logs/download/domain_annotations.log"
    shell:
        """
        mkdir -p data/prokaryotic
        wget -c -o {log} {params.url} -O {output.annotations}
        """


# ============================================================================
# GP Motif Extraction
# ============================================================================

rule extract_gp_motifs:
    """Extract all sequences containing GP motifs with context."""
    input:
        fasta="data/prokaryotic/uniprot_bacteria.fasta.gz"
    output:
        motifs="results/prokaryotic/gp_motifs/all_gp_motifs.tsv.gz",
        sequences="results/prokaryotic/gp_motifs/all_gp_sequences.fasta.gz"
    params:
        upstream=30,  # residues upstream of GP
        downstream=15  # residues downstream of GP
    log:
        "logs/prokaryotic/extract_gp_motifs.log"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/extract_gp_motifs.py"


rule annotate_domains:
    """Run InterProScan or Pfam scan on GP-containing sequences."""
    input:
        sequences="results/prokaryotic/gp_motifs/all_gp_sequences.fasta.gz"
    output:
        annotations="results/prokaryotic/gp_motifs/domain_annotations.tsv.gz"
    log:
        "logs/prokaryotic/annotate_domains.log"
    conda:
        "../envs/python.yaml"
    threads: 8
    resources:
        runtime=480,
        mem_mb=16000
    script:
        "../scripts/annotate_gp_domains.py"


rule filter_interdomain_gp:
    """Filter GP motifs that occur between protein domains."""
    input:
        motifs="results/prokaryotic/gp_motifs/all_gp_motifs.tsv.gz",
        annotations="results/prokaryotic/gp_motifs/domain_annotations.tsv.gz"
    output:
        interdomain="results/prokaryotic/gp_motifs/interdomain_gp_motifs.tsv.gz",
        intradomain="results/prokaryotic/gp_motifs/intradomain_gp_motifs.tsv.gz",
        statistics="results/prokaryotic/gp_motifs/gp_motif_stats.tsv"
    log:
        "logs/prokaryotic/filter_interdomain_gp.log"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/filter_interdomain_gp.py"


# ============================================================================
# Pattern Discovery and Clustering
# ============================================================================

rule cluster_gp_motifs:
    """Cluster GP motifs by sequence context using CD-HIT or MMseqs2."""
    input:
        interdomain="results/prokaryotic/gp_motifs/interdomain_gp_motifs.tsv.gz"
    output:
        clusters="results/prokaryotic/clusters/gp_clusters.tsv.gz",
        representatives="results/prokaryotic/clusters/cluster_representatives.fasta"
    params:
        identity=config["prokaryotic"]["clustering_identity"],
        coverage=config["prokaryotic"]["clustering_coverage"]
    log:
        "logs/prokaryotic/cluster_gp_motifs.log"
    conda:
        "../envs/python.yaml"
    threads: 8
    resources:
        runtime=120,
        mem_mb=16000
    script:
        "../scripts/cluster_gp_motifs.py"


rule analyze_cluster_conservation:
    """Analyze conservation patterns within each cluster."""
    input:
        clusters="results/prokaryotic/clusters/gp_clusters.tsv.gz",
        motifs="results/prokaryotic/gp_motifs/interdomain_gp_motifs.tsv.gz"
    output:
        conservation="results/prokaryotic/clusters/cluster_conservation.tsv.gz",
        logos="results/prokaryotic/clusters/logos/cluster_{cluster_id}.png"
    params:
        min_cluster_size=config["prokaryotic"]["min_cluster_size"]
    log:
        "logs/prokaryotic/analyze_conservation_{cluster_id}.log"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/analyze_gp_conservation.py"


rule identify_consensus_patterns:
    """Identify consensus patterns for major clusters."""
    input:
        conservation="results/prokaryotic/clusters/cluster_conservation.tsv.gz",
        clusters="results/prokaryotic/clusters/gp_clusters.tsv.gz"
    output:
        consensus="results/prokaryotic/consensus/consensus_patterns.tsv",
        alignments=expand(
            "results/prokaryotic/consensus/cluster_{cluster_id}.sto",
            cluster_id=range(1, 21)  # Top 20 clusters
        )
    params:
        min_conservation=config["prokaryotic"]["min_conservation"]
    log:
        "logs/prokaryotic/identify_consensus.log"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/identify_consensus_patterns.py"


# ============================================================================
# HMM Building for Prokaryotic Motifs
# ============================================================================

rule build_prokaryotic_hmms:
    """Build HMMs from consensus patterns."""
    input:
        alignment="results/prokaryotic/consensus/cluster_{cluster_id}.sto"
    output:
        hmm="results/prokaryotic/models/initial/cluster_{cluster_id}.hmm"
    params:
        name=lambda w: f"prok-2A-cluster-{w.cluster_id}"
    log:
        "logs/prokaryotic/build_hmm_cluster_{cluster_id}.log"
    conda:
        "../envs/hmmer.yaml"
    shell:
        """
        hmmbuild -n {params.name} {output.hmm} {input.alignment} 2> {log}
        """


rule search_prokaryotic_proteomes:
    """Search prokaryotic proteomes with discovered HMMs."""
    input:
        hmm="results/prokaryotic/models/initial/cluster_{cluster_id}.hmm",
        db="data/prokaryotic/uniprot_bacteria.fasta.gz"
    output:
        hmmsearch="results/prokaryotic/searches/cluster_{cluster_id}.hmmsearch.gz",
        tblout="results/prokaryotic/searches/cluster_{cluster_id}.tblout.gz",
        alignment="results/prokaryotic/searches/cluster_{cluster_id}.sto.gz"
    log:
        "logs/prokaryotic/search_cluster_{cluster_id}.log"
    conda:
        "../envs/hmmer.yaml"
    threads: 12
    resources:
        runtime=480,
        mem_mb=16000
    shell:
        """
        hmmsearch --cpu {threads} \
            --tblout >(gzip > {output.tblout}) \
            -A >(gzip > {output.alignment}) \
            --noali \
            {input.hmm} {input.db} 2> {log} | gzip > {output.hmmsearch}
        """


# ============================================================================
# Validation Against Known Stalling Peptides
# ============================================================================

rule validate_against_known_peptides:
    """Compare discovered motifs to known stalling peptides (SecM, TnaC, etc.)."""
    input:
        hmms=expand(
            "results/prokaryotic/models/initial/cluster_{cluster_id}.hmm",
            cluster_id=range(1, 21)
        ),
        known_peptides="resources/known_stalling_peptides.fasta"
    output:
        validation="results/prokaryotic/validation/known_peptide_hits.tsv",
        summary="results/prokaryotic/validation/validation_summary.txt"
    log:
        "logs/prokaryotic/validate_known_peptides.log"
    conda:
        "../envs/hmmer.yaml"
    threads: 4
    script:
        "../scripts/validate_stalling_peptides.py"


# ============================================================================
# Comparative Analysis: All GP vs Inter-domain GP
# ============================================================================

rule compare_approaches:
    """Compare results from all-GP vs inter-domain-GP approaches."""
    input:
        all_gp_motifs="results/prokaryotic/gp_motifs/all_gp_motifs.tsv.gz",
        interdomain_motifs="results/prokaryotic/gp_motifs/interdomain_gp_motifs.tsv.gz",
        clusters="results/prokaryotic/clusters/gp_clusters.tsv.gz",
        validation="results/prokaryotic/validation/known_peptide_hits.tsv"
    output:
        comparison="results/prokaryotic/analysis/approach_comparison.tsv",
        plots=directory("results/prokaryotic/analysis/comparison_plots/")
    log:
        "logs/prokaryotic/compare_approaches.log"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/compare_gp_approaches.py"


# ============================================================================
# Final Report
# ============================================================================

rule prokaryotic_discovery_report:
    """Generate comprehensive report on prokaryotic 2A-like discovery."""
    input:
        comparison="results/prokaryotic/analysis/approach_comparison.tsv",
        validation="results/prokaryotic/validation/validation_summary.txt",
        conservation="results/prokaryotic/clusters/cluster_conservation.tsv.gz",
        consensus="results/prokaryotic/consensus/consensus_patterns.tsv"
    output:
        report="results/prokaryotic/reports/prokaryotic_discovery.html"
    log:
        "logs/prokaryotic/generate_report.log"
    conda:
        "../envs/r-quarto.yaml"
    script:
        "../scripts/prokaryotic_discovery_report.qmd"
