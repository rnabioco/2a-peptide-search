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
        mkdir -p $(dirname "{output.fasta}")
        wget -c -o "{log}" "{params.url}" -O "{output.fasta}"
        """


rule download_pfam_database:
    """Download and decompress Pfam-A HMM database."""
    output:
        hmm=DATA_DIR + "/pfam/Pfam-A.hmm",
    params:
        url=config["prokaryotic_databases"]["pfam"]["url"],
    log:
        LOGS_DIR + "/download/pfam.log",
    shell:
        """
        mkdir -p $(dirname "{output.hmm}")
        wget -c -o "{log}" "{params.url}" -O "{output.hmm}.gz"
        gunzip -f "{output.hmm}.gz"
        """


# ============================================================================
# APPROACH 1: Seed-based Discovery with Known Stalling Peptides
# ============================================================================


rule split_known_peptides_by_motif:
    """Split known stalling peptides into separate files by motif type."""
    input:
        fasta="resources/stalling-peptides/known_stalling_peptides.fasta",
    output:
        motif_list=RESULTS_DIR + "/prokaryotic/seeds/motif_list.txt",
        fastas=directory(RESULTS_DIR + "/prokaryotic/seeds/by_motif/"),
    log:
        LOGS_DIR + "/prokaryotic/split_peptides_by_motif.log",
    shell:
        """
        python workflow/scripts/split_peptides_by_motif.py \
            --fasta "{input.fasta}" \
            --output-dir "{output.fastas}" \
            --motif-list "{output.motif_list}" \
            > "{log}" 2>&1
        """


rule align_all_seed_peptides:
    """Create comprehensive alignment from all known stalling peptides."""
    input:
        fasta="resources/stalling-peptides/known_stalling_peptides.fasta",
    output:
        alignment=RESULTS_DIR + "/prokaryotic/seeds/alignments/comprehensive.sto",
    log:
        LOGS_DIR + "/prokaryotic/align_all_seeds.log",
    shell:
        """
        mkdir -p $(dirname "{output.alignment}")

        # Use MAFFT for alignment, convert to Stockholm format
        mafft --auto "{input.fasta}" > "{output.alignment}.afa" 2> "{log}"

        # Convert to Stockholm format
        esl-reformat stockholm "{output.alignment}.afa" > "{output.alignment}" 2>> "{log}"

        rm "{output.alignment}.afa"
        """


rule build_comprehensive_hmm:
    """Build comprehensive 'pan-stalling' HMM from all known peptides."""
    input:
        alignment=RESULTS_DIR + "/prokaryotic/seeds/alignments/comprehensive.sto",
    output:
        hmm=RESULTS_DIR + "/prokaryotic/models/seed/comprehensive.hmm",
    params:
        name="stall-pan",
    log:
        LOGS_DIR + "/prokaryotic/build_comprehensive_hmm.log",
    shell:
        """
        hmmbuild -n {params.name} "{output.hmm}" "{input.alignment}" 2> "{log}"
        """


rule search_with_comprehensive_hmm:
    """Search prokaryotic proteomes with comprehensive pan-stalling HMM."""
    input:
        hmm=RESULTS_DIR + "/prokaryotic/models/seed/comprehensive.hmm",
        db=DATA_DIR + "/prokaryotic/uniprot_bacteria.fasta.gz",
    output:
        hmmsearch=RESULTS_DIR + "/prokaryotic/seed_searches/comprehensive.hmmsearch.gz",
        tblout=RESULTS_DIR + "/prokaryotic/seed_searches/comprehensive.tblout.gz",
        alignment=RESULTS_DIR + "/prokaryotic/seed_searches/comprehensive.sto.gz",
    log:
        LOGS_DIR + "/prokaryotic/search_comprehensive.log",
    shell:
        """
        hmmsearch --cpu {threads} \
            --tblout >(gzip > "{output.tblout}") \
            -A >(gzip > "{output.alignment}") \
            --noali \
            "{input.hmm}" "{input.db}" 2> "{log}" | gzip > "{output.hmmsearch}"
        """


rule align_seed_peptides:
    """Create multiple sequence alignment for each motif family."""
    input:
        fasta=RESULTS_DIR + "/prokaryotic/seeds/by_motif/{motif}.fasta",
    output:
        alignment=RESULTS_DIR + "/prokaryotic/seeds/alignments/{motif}.sto",
    log:
        LOGS_DIR + "/prokaryotic/align_seeds_{motif}.log",
    shell:
        """
        mkdir -p $(dirname "{output.alignment}")

        # Use MAFFT for alignment, convert to Stockholm format
        mafft --auto "{input.fasta}" > "{output.alignment}.afa" 2> "{log}"

        # Convert to Stockholm format (hmmer accepts various formats)
        esl-reformat stockholm "{output.alignment}.afa" > "{output.alignment}" 2>> "{log}"

        rm "{output.alignment}.afa"
        """


rule build_seed_hmms:
    """Build HMMs from seed alignments of known stalling peptides."""
    input:
        alignment=RESULTS_DIR + "/prokaryotic/seeds/alignments/{motif}.sto",
    output:
        hmm=RESULTS_DIR + "/prokaryotic/models/seed/{motif}.hmm",
    params:
        name=lambda w: f"stall-{w.motif}",
    log:
        LOGS_DIR + "/prokaryotic/build_seed_hmm_{motif}.log",
    shell:
        """
        hmmbuild -n {params.name} "{output.hmm}" "{input.alignment}" 2> "{log}"
        """


rule search_with_seed_hmms:
    """Search prokaryotic proteomes with seed HMMs from known peptides."""
    input:
        hmm=RESULTS_DIR + "/prokaryotic/models/seed/{motif}.hmm",
        db=DATA_DIR + "/prokaryotic/uniprot_bacteria.fasta.gz",
    output:
        hmmsearch=RESULTS_DIR + "/prokaryotic/seed_searches/{motif}.hmmsearch.gz",
        tblout=RESULTS_DIR + "/prokaryotic/seed_searches/{motif}.tblout.gz",
        alignment=RESULTS_DIR + "/prokaryotic/seed_searches/{motif}.sto.gz",
    log:
        LOGS_DIR + "/prokaryotic/seed_search_{motif}.log",
    shell:
        """
        hmmsearch --cpu {threads} \
            --tblout >(gzip > "{output.tblout}") \
            -A >(gzip > "{output.alignment}") \
            --noali \
            "{input.hmm}" "{input.db}" 2> "{log}" | gzip > "{output.hmmsearch}"
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
    shell:
        """
        python workflow/scripts/extract_gp_motifs.py \
            --fasta "{input.fasta}" \
            --motifs-out "{output.motifs}" \
            --sequences-out "{output.sequences}" \
            --upstream {params.upstream} \
            --downstream {params.downstream} \
            > "{log}" 2>&1
        """


rule run_hmmscan:
    """Run hmmscan to annotate domains in GP-containing sequences."""
    input:
        sequences=RESULTS_DIR + "/prokaryotic/gp_motifs/all_gp_sequences.fasta.gz",
        pfam_db=DATA_DIR + "/pfam/Pfam-A.hmm",
    output:
        domtblout=RESULTS_DIR + "/prokaryotic/gp_motifs/domains.domtblout",
    log:
        LOGS_DIR + "/prokaryotic/hmmscan.log",
    threads: 8
    resources:
        runtime=480,
        mem_mb=16000,
    shell:
        """
        # Decompress sequences and run hmmscan
        zcat "{input.sequences}" | hmmscan \
            --cpu {threads} \
            --domtblout "{output.domtblout}" \
            --cut_ga \
            "{input.pfam_db}" \
            - \
            > "{log}" 2>&1
        """


rule parse_domain_annotations:
    """Parse hmmscan output into domain annotations table."""
    input:
        domtblout=RESULTS_DIR + "/prokaryotic/gp_motifs/domains.domtblout",
        motifs=RESULTS_DIR + "/prokaryotic/gp_motifs/all_gp_motifs.tsv.gz",
    output:
        annotations=RESULTS_DIR + "/prokaryotic/gp_motifs/domain_annotations.tsv.gz",
    log:
        LOGS_DIR + "/prokaryotic/parse_annotations.log",
    shell:
        """
        python workflow/scripts/parse_domain_annotations.py \
            --domtblout "{input.domtblout}" \
            --motifs "{input.motifs}" \
            --annotations "{output.annotations}" \
            > "{log}" 2>&1
        """


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
    shell:
        """
        python workflow/scripts/filter_interdomain_gp.py \
            --motifs "{input.motifs}" \
            --annotations "{input.annotations}" \
            --interdomain-out "{output.interdomain}" \
            --intradomain-out "{output.intradomain}" \
            --stats-out "{output.statistics}" \
            > "{log}" 2>&1
        """


# ============================================================================
# Pattern Discovery and Clustering
# ============================================================================


rule prepare_clustering_fasta:
    """Extract GP context sequences for clustering."""
    input:
        interdomain=RESULTS_DIR + "/prokaryotic/gp_motifs/interdomain_gp_motifs.tsv.gz",
    output:
        fasta=RESULTS_DIR + "/prokaryotic/clusters/gp_motifs.fasta",
    shell:
        """
        # Extract context sequences to FASTA
        zcat "{input.interdomain}" | awk -F'\t' 'NR>1 {{
            printf ">%s_GP%s_pos%s\\n%s\\n", $1, $4, $5, $6
        }}' > "{output.fasta}"
        """


rule run_mmseqs_clustering:
    """Cluster GP motifs using MMseqs2."""
    input:
        fasta=RESULTS_DIR + "/prokaryotic/clusters/gp_motifs.fasta",
    output:
        cluster_tsv=RESULTS_DIR + "/prokaryotic/clusters/mmseqs_clusters.tsv",
    params:
        identity=config["prokaryotic"]["clustering_identity"],
        coverage=config["prokaryotic"]["clustering_coverage"],
        prefix=RESULTS_DIR + "/prokaryotic/clusters/mmseqs",
        tmpdir=RESULTS_DIR + "/prokaryotic/clusters/tmp",
    log:
        LOGS_DIR + "/prokaryotic/mmseqs_cluster.log",
    threads: 8
    resources:
        runtime=120,
        mem_mb=16000,
    shell:
        """
        mkdir -p "{params.tmpdir}"

        # Create MMseqs2 database
        mmseqs createdb "{input.fasta}" "{params.prefix}_db" 2>> "{log}"

        # Cluster
        mmseqs cluster \
            "{params.prefix}_db" \
            "{params.prefix}_cluster" \
            "{params.tmpdir}" \
            --min-seq-id {params.identity} \
            -c {params.coverage} \
            --threads {threads} \
            2>> "{log}"

        # Create TSV output
        mmseqs createtsv \
            "{params.prefix}_db" \
            "{params.prefix}_db" \
            "{params.prefix}_cluster" \
            "{output.cluster_tsv}" \
            2>> "{log}"

        # Cleanup temp files
        rm -rf "{params.tmpdir}" "{params.prefix}_db"* "{params.prefix}_cluster"*
        """


rule parse_gp_clusters:
    """Parse MMseqs2 clustering results."""
    input:
        cluster_tsv=RESULTS_DIR + "/prokaryotic/clusters/mmseqs_clusters.tsv",
        interdomain=RESULTS_DIR + "/prokaryotic/gp_motifs/interdomain_gp_motifs.tsv.gz",
    output:
        clusters=RESULTS_DIR + "/prokaryotic/clusters/gp_clusters.tsv.gz",
        representatives=RESULTS_DIR + "/prokaryotic/clusters/cluster_representatives.fasta",
    log:
        LOGS_DIR + "/prokaryotic/parse_clusters.log",
    shell:
        """
        python workflow/scripts/parse_gp_clusters.py \
            --cluster-tsv "{input.cluster_tsv}" \
            --interdomain "{input.interdomain}" \
            --clusters-out "{output.clusters}" \
            --representatives-out "{output.representatives}" \
            > "{log}" 2>&1
        """


rule analyze_cluster_conservation:
    """Analyze conservation patterns within each cluster."""
    input:
        clusters=RESULTS_DIR + "/prokaryotic/clusters/gp_clusters.tsv.gz",
        motifs=RESULTS_DIR + "/prokaryotic/gp_motifs/interdomain_gp_motifs.tsv.gz",
    output:
        conservation=RESULTS_DIR + "/prokaryotic/clusters/cluster_conservation.tsv.gz",
        logos=expand(
            RESULTS_DIR + "/prokaryotic/clusters/logos/cluster_{cluster_id}.png",
            cluster_id=range(1, 21),  # Top 20 clusters
        ),
    params:
        min_cluster_size=config["prokaryotic"]["min_cluster_size"],
        logos_dir=RESULTS_DIR + "/prokaryotic/clusters/logos",
    log:
        LOGS_DIR + "/prokaryotic/analyze_conservation.log",
    conda:
        "../envs/python.yaml"
    shell:
        """
        python workflow/scripts/analyze_gp_conservation.py \
            --clusters "{input.clusters}" \
            --motifs "{input.motifs}" \
            --conservation "{output.conservation}" \
            --logos-dir "{params.logos_dir}" \
            --min-cluster-size {params.min_cluster_size} \
            --top-n 20 \
            > "{log}" 2>&1
        """


rule identify_consensus_patterns:
    """Identify consensus patterns for major clusters."""
    input:
        conservation=RESULTS_DIR + "/prokaryotic/clusters/cluster_conservation.tsv.gz",
        clusters=RESULTS_DIR + "/prokaryotic/clusters/gp_clusters.tsv.gz",
        motifs=RESULTS_DIR + "/prokaryotic/gp_motifs/interdomain_gp_motifs.tsv.gz",
    output:
        consensus=RESULTS_DIR + "/prokaryotic/consensus/consensus_patterns.tsv",
        alignments=expand(
            RESULTS_DIR + "/prokaryotic/consensus/cluster_{cluster_id}.sto",
            cluster_id=range(1, 21),  # Top 20 clusters
        ),
    params:
        min_conservation=config["prokaryotic"]["min_conservation"],
        min_cluster_size=config["prokaryotic"]["min_cluster_size"],
        alignments_dir=RESULTS_DIR + "/prokaryotic/consensus",
    log:
        LOGS_DIR + "/prokaryotic/identify_consensus.log",
    conda:
        "../envs/python.yaml"
    shell:
        """
        python workflow/scripts/identify_consensus_patterns.py \
            --conservation "{input.conservation}" \
            --clusters "{input.clusters}" \
            --motifs "{input.motifs}" \
            --consensus "{output.consensus}" \
            --alignments-dir "{params.alignments_dir}" \
            --min-conservation {params.min_conservation} \
            --min-cluster-size {params.min_cluster_size} \
            --top-n 20 \
            > "{log}" 2>&1
        """


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
        hmmbuild -n {params.name} "{output.hmm}" "{input.alignment}" 2> "{log}"
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
            --tblout >(gzip > "{output.tblout}") \
            -A >(gzip > "{output.alignment}") \
            --noali \
            "{input.hmm}" "{input.db}" 2> "{log}" | gzip > "{output.hmmsearch}"
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
        known_peptides="resources/stalling-peptides/known_stalling_peptides.fasta",
    output:
        validation=RESULTS_DIR + "/prokaryotic/validation/known_peptide_hits.tsv",
        summary=RESULTS_DIR + "/prokaryotic/validation/validation_summary.txt",
    params:
        hmms_pattern=RESULTS_DIR + "/prokaryotic/models/initial/cluster_*.hmm",
    log:
        LOGS_DIR + "/prokaryotic/validate_known_peptides.log",
    conda:
        "../envs/hmmer.yaml"
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
    """Compare results from all-GP vs inter-domain-GP approaches."""
    input:
        all_gp_motifs=RESULTS_DIR + "/prokaryotic/gp_motifs/all_gp_motifs.tsv.gz",
        interdomain_motifs=RESULTS_DIR
        + "/prokaryotic/gp_motifs/interdomain_gp_motifs.tsv.gz",
        clusters=RESULTS_DIR + "/prokaryotic/clusters/gp_clusters.tsv.gz",
        validation=RESULTS_DIR + "/prokaryotic/validation/known_peptide_hits.tsv",
        # APPROACH 1: Seed-based searches
        seed_comprehensive=RESULTS_DIR + "/prokaryotic/seed_searches/comprehensive.sto.gz",
    output:
        comparison=RESULTS_DIR + "/prokaryotic/analysis/approach_comparison.tsv",
        plots=directory(RESULTS_DIR + "/prokaryotic/analysis/comparison_plots/"),
    log:
        LOGS_DIR + "/prokaryotic/compare_approaches.log",
    conda:
        "../envs/python.yaml"
    shell:
        """
        python workflow/scripts/compare_gp_approaches.py \
            --all-gp-motifs "{input.all_gp_motifs}" \
            --interdomain-motifs "{input.interdomain_motifs}" \
            --clusters "{input.clusters}" \
            --validation "{input.validation}" \
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
