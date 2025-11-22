"""
Rules for discovering prokaryotic 2A-like ribosomal stalling peptides.

Strategy (Two complementary approaches):

APPROACH 1: Seed-based discovery (Targeted)
1. Start with known stalling peptides (PMID 38565864)
2. Align by motif type (RAGP, RAPG, QAPP, etc.)
3. Build seed HMMs from each motif family
4. Search prokaryotic proteomes with seed HMMs
5. Iterative refinement

APPROACH 2: Unbiased GP discovery (Comprehensive)
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
        fasta=DATA_DIR + "/prokaryotic/uniprot_bacteria.fasta.gz",
    params:
        url=config["prokaryotic_databases"]["bacteria"]["url"],
    log:
        LOGS_DIR + "/download/prokaryotic_proteomes.log",
    shell:
        """
        mkdir -p data/prokaryotic
        wget -c -o {log} {params.url} -O {output.fasta}
        """


rule download_domain_annotations:
    """Download Pfam or pre-computed domain annotations."""
    output:
        annotations=DATA_DIR + "/prokaryotic/domain_annotations.tsv.gz",
    params:
        url=config["prokaryotic_databases"]["pfam"]["url"],
    log:
        LOGS_DIR + "/download/domain_annotations.log",
    shell:
        """
        mkdir -p data/prokaryotic
        wget -c -o {log} {params.url} -O {output.annotations}
        """


# ============================================================================
# APPROACH 1: Seed-based Discovery with Known Stalling Peptides
# ============================================================================

rule split_known_peptides_by_motif:
    """Split known stalling peptides into separate files by motif type."""
    input:
        fasta="resources/stalling-peptides/known_stalling_peptides.fasta"
    output:
        motif_list=RESULTS_DIR + "/prokaryotic/seeds/motif_list.txt",
        fastas=directory(RESULTS_DIR + "/prokaryotic/seeds/by_motif/")
    log:
        LOGS_DIR + "/prokaryotic/split_peptides_by_motif.log"
    script:
        "../scripts/split_peptides_by_motif.py"


rule align_all_seed_peptides:
    """Create comprehensive alignment from all known stalling peptides."""
    input:
        fasta="resources/stalling-peptides/known_stalling_peptides.fasta"
    output:
        alignment=RESULTS_DIR + "/prokaryotic/seeds/alignments/comprehensive.sto"
    log:
        LOGS_DIR + "/prokaryotic/align_all_seeds.log"
    shell:
        """
        # Use MUSCLE for alignment, convert to Stockholm format
        muscle -align {input.fasta} -output {output.alignment}.afa 2> {log}
        # Convert to Stockholm format
        esl-reformat stockholm {output.alignment}.afa > {output.alignment} 2>> {log}
        rm {output.alignment}.afa
        """


rule build_comprehensive_hmm:
    """Build comprehensive 'pan-stalling' HMM from all known peptides."""
    input:
        alignment=RESULTS_DIR + "/prokaryotic/seeds/alignments/comprehensive.sto"
    output:
        hmm=RESULTS_DIR + "/prokaryotic/models/seed/comprehensive.hmm"
    params:
        name="stall-pan"
    log:
        LOGS_DIR + "/prokaryotic/build_comprehensive_hmm.log"
    shell:
        """
        hmmbuild -n {params.name} {output.hmm} {input.alignment} 2> {log}
        """


rule search_with_comprehensive_hmm:
    """Search prokaryotic proteomes with comprehensive pan-stalling HMM."""
    input:
        hmm=RESULTS_DIR + "/prokaryotic/models/seed/comprehensive.hmm",
        db=DATA_DIR + "/prokaryotic/uniprot_bacteria.fasta.gz"
    output:
        hmmsearch=RESULTS_DIR + "/prokaryotic/seed_searches/comprehensive.hmmsearch.gz",
        tblout=RESULTS_DIR + "/prokaryotic/seed_searches/comprehensive.tblout.gz",
        alignment=RESULTS_DIR + "/prokaryotic/seed_searches/comprehensive.sto.gz"
    log:
        LOGS_DIR + "/prokaryotic/search_comprehensive.log"
    shell:
        """
        hmmsearch --cpu {threads} \
            --tblout >(gzip > {output.tblout}) \
            -A >(gzip > {output.alignment}) \
            --noali \
            {input.hmm} {input.db} 2> {log} | gzip > {output.hmmsearch}
        """


rule align_seed_peptides:
    """Create multiple sequence alignment for each motif family."""
    input:
        fasta=RESULTS_DIR + "/prokaryotic/seeds/by_motif/{motif}.fasta"
    output:
        alignment=RESULTS_DIR + "/prokaryotic/seeds/alignments/{motif}.sto"
    log:
        LOGS_DIR + "/prokaryotic/align_seeds_{motif}.log"
    shell:
        """
        # Use MUSCLE for alignment, convert to Stockholm format
        muscle -align {input.fasta} -output {output.alignment}.afa 2> {log}
        # Convert to Stockholm format (hmmer accepts various formats)
        esl-reformat stockholm {output.alignment}.afa > {output.alignment} 2>> {log}
        rm {output.alignment}.afa
        """


rule build_seed_hmms:
    """Build HMMs from seed alignments of known stalling peptides."""
    input:
        alignment=RESULTS_DIR + "/prokaryotic/seeds/alignments/{motif}.sto"
    output:
        hmm=RESULTS_DIR + "/prokaryotic/models/seed/{motif}.hmm"
    params:
        name=lambda w: f"stall-{w.motif}"
    log:
        LOGS_DIR + "/prokaryotic/build_seed_hmm_{motif}.log"
    shell:
        """
        hmmbuild -n {params.name} {output.hmm} {input.alignment} 2> {log}
        """


rule search_with_seed_hmms:
    """Search prokaryotic proteomes with seed HMMs from known peptides."""
    input:
        hmm=RESULTS_DIR + "/prokaryotic/models/seed/{motif}.hmm",
        db=DATA_DIR + "/prokaryotic/uniprot_bacteria.fasta.gz"
    output:
        hmmsearch=RESULTS_DIR + "/prokaryotic/seed_searches/{motif}.hmmsearch.gz",
        tblout=RESULTS_DIR + "/prokaryotic/seed_searches/{motif}.tblout.gz",
        alignment=RESULTS_DIR + "/prokaryotic/seed_searches/{motif}.sto.gz"
    log:
        LOGS_DIR + "/prokaryotic/seed_search_{motif}.log"
    shell:
        """
        hmmsearch --cpu {threads} \
            --tblout >(gzip > {output.tblout}) \
            -A >(gzip > {output.alignment}) \
            --noali \
            {input.hmm} {input.db} 2> {log} | gzip > {output.hmmsearch}
        """


# ============================================================================
# APPROACH 2: Unbiased GP Motif Extraction
# ============================================================================


rule extract_gp_motifs:
    """Extract all sequences containing GP motifs with context."""
    input:
        fasta=DATA_DIR + "/prokaryotic/uniprot_bacteria.fasta.gz",
    output:
        motifs=RESULTS_DIR + "/prokaryotic/gp_motifs/all_gp_motifs.tsv.gz",
        sequences=RESULTS_DIR + "/prokaryotic/gp_motifs/all_gp_sequences.fasta.gz",
    params:
        upstream=30,  # residues upstream of GP
        downstream=15,  # residues downstream of GP
    log:
        LOGS_DIR + "/prokaryotic/extract_gp_motifs.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/extract_gp_motifs.py"


rule annotate_domains:
    """Run InterProScan or Pfam scan on GP-containing sequences."""
    input:
        sequences=RESULTS_DIR + "/prokaryotic/gp_motifs/all_gp_sequences.fasta.gz",
    output:
        annotations=RESULTS_DIR + "/prokaryotic/gp_motifs/domain_annotations.tsv.gz",
    log:
        LOGS_DIR + "/prokaryotic/annotate_domains.log",
    conda:
        "../envs/python.yaml"
    threads: 8
    resources:
        runtime=480,
        mem_mb=16000,
    script:
        "../scripts/annotate_gp_domains.py"


rule filter_interdomain_gp:
    """Filter GP motifs that occur between protein domains."""
    input:
        motifs=RESULTS_DIR + "/prokaryotic/gp_motifs/all_gp_motifs.tsv.gz",
        annotations=RESULTS_DIR + "/prokaryotic/gp_motifs/domain_annotations.tsv.gz",
    output:
        interdomain=RESULTS_DIR + "/prokaryotic/gp_motifs/interdomain_gp_motifs.tsv.gz",
        intradomain=RESULTS_DIR + "/prokaryotic/gp_motifs/intradomain_gp_motifs.tsv.gz",
        statistics=RESULTS_DIR + "/prokaryotic/gp_motifs/gp_motif_stats.tsv",
    log:
        LOGS_DIR + "/prokaryotic/filter_interdomain_gp.log",
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
        interdomain=RESULTS_DIR + "/prokaryotic/gp_motifs/interdomain_gp_motifs.tsv.gz",
    output:
        clusters=RESULTS_DIR + "/prokaryotic/clusters/gp_clusters.tsv.gz",
        representatives=RESULTS_DIR
        + "/prokaryotic/clusters/cluster_representatives.fasta",
    params:
        identity=config["prokaryotic"]["clustering_identity"],
        coverage=config["prokaryotic"]["clustering_coverage"],
    log:
        LOGS_DIR + "/prokaryotic/cluster_gp_motifs.log",
    conda:
        "../envs/python.yaml"
    threads: 8
    resources:
        runtime=120,
        mem_mb=16000,
    script:
        "../scripts/cluster_gp_motifs.py"


rule analyze_cluster_conservation:
    """Analyze conservation patterns within each cluster."""
    input:
        clusters=RESULTS_DIR + "/prokaryotic/clusters/gp_clusters.tsv.gz",
        motifs=RESULTS_DIR + "/prokaryotic/gp_motifs/interdomain_gp_motifs.tsv.gz",
    output:
        conservation=RESULTS_DIR + "/prokaryotic/clusters/cluster_conservation.tsv.gz",
        logos=RESULTS_DIR + "/prokaryotic/clusters/logos/cluster_{cluster_id}.png",
    params:
        min_cluster_size=config["prokaryotic"]["min_cluster_size"],
    log:
        LOGS_DIR + "/prokaryotic/analyze_conservation_{cluster_id}.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/analyze_gp_conservation.py"


rule identify_consensus_patterns:
    """Identify consensus patterns for major clusters."""
    input:
        conservation=RESULTS_DIR + "/prokaryotic/clusters/cluster_conservation.tsv.gz",
        clusters=RESULTS_DIR + "/prokaryotic/clusters/gp_clusters.tsv.gz",
    output:
        consensus=RESULTS_DIR + "/prokaryotic/consensus/consensus_patterns.tsv",
        alignments=expand(
            RESULTS_DIR + "/prokaryotic/consensus/cluster_{cluster_id}.sto",
            cluster_id=range(1, 21),  # Top 20 clusters
        ),
    params:
        min_conservation=config["prokaryotic"]["min_conservation"],
    log:
        LOGS_DIR + "/prokaryotic/identify_consensus.log",
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
        alignment=RESULTS_DIR + "/prokaryotic/consensus/cluster_{cluster_id}.sto",
    output:
        hmm=RESULTS_DIR + "/prokaryotic/models/initial/cluster_{cluster_id}.hmm",
    params:
        name=lambda w: f"prok-2A-cluster-{w.cluster_id}",
    log:
        LOGS_DIR + "/prokaryotic/build_hmm_cluster_{cluster_id}.log",
    conda:
        "../envs/hmmer.yaml"
    shell:
        """
        hmmbuild -n {params.name} {output.hmm} {input.alignment} 2> {log}
        """


rule search_prokaryotic_proteomes:
    """Search prokaryotic proteomes with discovered HMMs."""
    input:
        hmm=RESULTS_DIR + "/prokaryotic/models/initial/cluster_{cluster_id}.hmm",
        db=DATA_DIR + "/prokaryotic/uniprot_bacteria.fasta.gz",
    output:
        hmmsearch=RESULTS_DIR
        + "/prokaryotic/searches/cluster_{cluster_id}.hmmsearch.gz",
        tblout=RESULTS_DIR + "/prokaryotic/searches/cluster_{cluster_id}.tblout.gz",
        alignment=RESULTS_DIR + "/prokaryotic/searches/cluster_{cluster_id}.sto.gz",
    log:
        LOGS_DIR + "/prokaryotic/search_cluster_{cluster_id}.log",
    conda:
        "../envs/hmmer.yaml"
    threads: 12
    resources:
        runtime=480,
        mem_mb=16000,
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
            RESULTS_DIR + "/prokaryotic/models/initial/cluster_{cluster_id}.hmm",
            cluster_id=range(1, 21),
        ),
        known_peptides="resources/stalling-peptides/known_stalling_peptides.fasta"
    output:
        validation=RESULTS_DIR + "/prokaryotic/validation/known_peptide_hits.tsv",
        summary=RESULTS_DIR + "/prokaryotic/validation/validation_summary.txt",
    log:
        LOGS_DIR + "/prokaryotic/validate_known_peptides.log",
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
        all_gp_motifs=RESULTS_DIR + "/prokaryotic/gp_motifs/all_gp_motifs.tsv.gz",
        interdomain_motifs=RESULTS_DIR
        + "/prokaryotic/gp_motifs/interdomain_gp_motifs.tsv.gz",
        clusters=RESULTS_DIR + "/prokaryotic/clusters/gp_clusters.tsv.gz",
        validation=RESULTS_DIR + "/prokaryotic/validation/known_peptide_hits.tsv",
    output:
        comparison=RESULTS_DIR + "/prokaryotic/analysis/approach_comparison.tsv",
        plots=directory(RESULTS_DIR + "/prokaryotic/analysis/comparison_plots/"),
    log:
        LOGS_DIR + "/prokaryotic/compare_approaches.log",
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
        comparison=RESULTS_DIR + "/prokaryotic/analysis/approach_comparison.tsv",
        validation=RESULTS_DIR + "/prokaryotic/validation/validation_summary.txt",
        conservation=RESULTS_DIR + "/prokaryotic/clusters/cluster_conservation.tsv.gz",
        consensus=RESULTS_DIR + "/prokaryotic/consensus/consensus_patterns.tsv",
    output:
        report=RESULTS_DIR + "/prokaryotic/reports/prokaryotic_discovery.html",
    log:
        LOGS_DIR + "/prokaryotic/generate_report.log",
    conda:
        "../envs/r-quarto.yaml"
    script:
        "../scripts/prokaryotic_discovery_report.qmd"
