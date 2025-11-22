"""
Rules for phage genome ORF prediction and GP motif discovery.

Workflow:
1. Download phage genome databases
2. Predict ORFs using Prodigal
3. Extract GP motifs from predicted ORFs
4. Integrate with prokaryotic discovery pipeline
"""

# ============================================================================
# Phage Genome Data
# ============================================================================


rule download_phage_genomes:
    """Download phage genome databases."""
    output:
        genomes=DATA_DIR + "/phage/{db}/genomes.fasta.gz",
    params:
        url=lambda w: config["phage_databases"][w.db]["url"],
    log:
        LOGS_DIR + "/download/phage_{db}.log",
    shell:
        """
        mkdir -p data/phage/{wildcards.db}
        wget -c -o {log} {params.url} -O {output.genomes}
        """


rule download_imgvr_phages:
    """Download IMG/VR phage database."""
    output:
        genomes=DATA_DIR + "/phage/imgvr/IMGVR_all_nucleotides.fna.gz",
    params:
        url=config["phage_databases"]["imgvr"]["url"],
    log:
        LOGS_DIR + "/download/phage_imgvr.log",
    shell:
        """
        mkdir -p data/phage/imgvr
        wget -c -o {log} {params.url} -O {output.genomes}
        """


rule download_inphared:
    """Download INPHARED phage database (curated phage genomes)."""
    output:
        genomes=DATA_DIR + "/phage/inphared/genomes.fasta.gz",
        metadata=DATA_DIR + "/phage/inphared/metadata.tsv",
    params:
        genomes_url=config["phage_databases"]["inphared"]["genomes_url"],
        metadata_url=config["phage_databases"]["inphared"]["metadata_url"],
    log:
        LOGS_DIR + "/download/phage_inphared.log",
    shell:
        """
        mkdir -p data/phage/inphared
        wget -c {params.genomes_url} -O {output.genomes} 2>> {log}
        wget -c {params.metadata_url} -O {output.metadata} 2>> {log}
        """


# ============================================================================
# ORF Prediction
# ============================================================================


rule split_phage_genomes:
    """Split large phage genome files into chunks for parallel ORF prediction."""
    input:
        genomes=DATA_DIR + "/phage/{db}/genomes.fasta.gz",
    output:
        chunks=expand(
            DATA_DIR + "/phage/{{db}}/chunks/chunk_{chunk}.fasta",
            chunk=range(1, 101),  # 100 chunks
        ),
    params:
        seqs_per_chunk=1000,
    log:
        LOGS_DIR + "/phage/split_{db}.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/split_fasta.py"


rule predict_orfs_prodigal:
    """Predict ORFs from phage genomes using Prodigal (metagenomics mode)."""
    input:
        chunk=DATA_DIR + "/phage/{db}/chunks/chunk_{chunk}.fasta",
    output:
        proteins=RESULTS_DIR + "/phage/{db}/orfs/chunk_{chunk}.faa",
        genes=RESULTS_DIR + "/phage/{db}/orfs/chunk_{chunk}.fna",
        gff=RESULTS_DIR + "/phage/{db}/orfs/chunk_{chunk}.gff",
    log:
        LOGS_DIR + "/phage/prodigal_{db}_chunk_{chunk}.log",
    conda:
        "../envs/orfs.yaml"
    threads: 1
    resources:
        runtime=60,
        mem_mb=4000,
    shell:
        """
        prodigal -i {input.chunk} \
            -a {output.proteins} \
            -d {output.genes} \
            -f gff \
            -o {output.gff} \
            -p meta \
            2> {log}
        """


rule merge_predicted_orfs:
    """Merge predicted ORFs from all chunks."""
    input:
        proteins=expand(
            RESULTS_DIR + "/phage/{{db}}/orfs/chunk_{chunk}.faa", chunk=range(1, 101)
        ),
    output:
        merged=RESULTS_DIR + "/phage/{db}/predicted_orfs.faa.gz",
    log:
        LOGS_DIR + "/phage/merge_orfs_{db}.log",
    shell:
        """
        cat {input.proteins} | gzip > {output.merged}
        echo "Merged $(zcat {output.merged} | grep -c '^>') ORFs" > {log}
        """


# ============================================================================
# GP Motif Extraction from Phage ORFs
# ============================================================================


rule extract_phage_gp_motifs:
    """Extract GP motifs from predicted phage ORFs."""
    input:
        fasta=RESULTS_DIR + "/phage/{db}/predicted_orfs.faa.gz",
    output:
        motifs=RESULTS_DIR + "/phage/{db}/gp_motifs/all_gp_motifs.tsv.gz",
        sequences=RESULTS_DIR + "/phage/{db}/gp_motifs/all_gp_sequences.fasta.gz",
    params:
        upstream=30,
        downstream=15,
    log:
        LOGS_DIR + "/phage/extract_gp_{db}.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/extract_gp_motifs.py"


rule annotate_phage_domains:
    """Run domain annotation on phage ORFs with GP motifs."""
    input:
        sequences=RESULTS_DIR + "/phage/{db}/gp_motifs/all_gp_sequences.fasta.gz",
    output:
        annotations=RESULTS_DIR + "/phage/{db}/gp_motifs/domain_annotations.tsv.gz",
    log:
        LOGS_DIR + "/phage/annotate_domains_{db}.log",
    conda:
        "../envs/python.yaml"
    threads: 8
    resources:
        runtime=480,
        mem_mb=16000,
    script:
        "../scripts/annotate_gp_domains.py"


rule filter_phage_interdomain_gp:
    """Filter phage GP motifs by domain boundaries."""
    input:
        motifs=RESULTS_DIR + "/phage/{db}/gp_motifs/all_gp_motifs.tsv.gz",
        annotations=RESULTS_DIR + "/phage/{db}/gp_motifs/domain_annotations.tsv.gz",
    output:
        interdomain=RESULTS_DIR + "/phage/{db}/gp_motifs/interdomain_gp_motifs.tsv.gz",
        intradomain=RESULTS_DIR + "/phage/{db}/gp_motifs/intradomain_gp_motifs.tsv.gz",
        statistics=RESULTS_DIR + "/phage/{db}/gp_motifs/gp_motif_stats.tsv",
    log:
        LOGS_DIR + "/phage/filter_interdomain_{db}.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/filter_interdomain_gp.py"


# ============================================================================
# Merge Prokaryotic and Phage Results
# ============================================================================


rule merge_prokaryotic_and_phage_motifs:
    """Combine GP motifs from prokaryotic proteomes and phage ORFs."""
    input:
        prok_interdomain=RESULTS_DIR
        + "/prokaryotic/gp_motifs/interdomain_gp_motifs.tsv.gz",
        phage_interdomain=expand(
            RESULTS_DIR + "/phage/{db}/gp_motifs/interdomain_gp_motifs.tsv.gz",
            db=config["phage_databases_to_use"],
        ),
    output:
        merged=RESULTS_DIR + "/combined/interdomain_gp_motifs_all.tsv.gz",
        summary=RESULTS_DIR + "/combined/source_summary.tsv",
    log:
        LOGS_DIR + "/combined/merge_motifs.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/merge_prokaryotic_phage_motifs.py"


rule cluster_combined_motifs:
    """Cluster combined prokaryotic and phage GP motifs."""
    input:
        interdomain=RESULTS_DIR + "/combined/interdomain_gp_motifs_all.tsv.gz",
    output:
        clusters=RESULTS_DIR + "/combined/clusters/gp_clusters.tsv.gz",
        representatives=RESULTS_DIR + "/combined/clusters/cluster_representatives.fasta",
    params:
        identity=config["prokaryotic"]["clustering_identity"],
        coverage=config["prokaryotic"]["clustering_coverage"],
    log:
        LOGS_DIR + "/combined/cluster_motifs.log",
    conda:
        "../envs/python.yaml"
    threads: 12
    resources:
        runtime=240,
        mem_mb=32000,
    script:
        "../scripts/cluster_gp_motifs.py"


# ============================================================================
# Phage-Specific Analysis
# ============================================================================


rule analyze_phage_gp_distribution:
    """Analyze distribution of GP motifs in phage genomes."""
    input:
        motifs=expand(
            RESULTS_DIR + "/phage/{db}/gp_motifs/all_gp_motifs.tsv.gz",
            db=config["phage_databases_to_use"],
        ),
        metadata=DATA_DIR + "/phage/inphared/metadata.tsv",
    output:
        distribution=RESULTS_DIR + "/phage/analysis/gp_distribution.tsv",
        plots=directory(RESULTS_DIR + "/phage/analysis/distribution_plots/"),
    log:
        LOGS_DIR + "/phage/analyze_distribution.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/analyze_phage_gp_distribution.py"


rule compare_phage_vs_bacterial:
    """Compare GP motif characteristics between phages and bacteria."""
    input:
        phage_motifs=RESULTS_DIR + "/combined/interdomain_gp_motifs_all.tsv.gz",
        clusters=RESULTS_DIR + "/combined/clusters/gp_clusters.tsv.gz",
    output:
        comparison=RESULTS_DIR + "/combined/analysis/phage_vs_bacterial.tsv",
        plots=directory(RESULTS_DIR + "/combined/analysis/comparison_plots/"),
    log:
        LOGS_DIR + "/combined/compare_phage_bacterial.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/compare_phage_vs_bacterial.py"
