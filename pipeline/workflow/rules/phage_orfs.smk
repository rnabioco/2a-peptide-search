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
        genomes="data/phage/{db}/genomes.fasta.gz"
    params:
        url=lambda w: config["phage_databases"][w.db]["url"]
    log:
        "logs/download/phage_{db}.log"
    shell:
        """
        mkdir -p data/phage/{wildcards.db}
        wget -c -o {log} {params.url} -O {output.genomes}
        """


rule download_imgvr_phages:
    """Download IMG/VR phage database."""
    output:
        genomes="data/phage/imgvr/IMGVR_all_nucleotides.fna.gz"
    params:
        url=config["phage_databases"]["imgvr"]["url"]
    log:
        "logs/download/phage_imgvr.log"
    shell:
        """
        mkdir -p data/phage/imgvr
        wget -c -o {log} {params.url} -O {output.genomes}
        """


rule download_inphared:
    """Download INPHARED phage database (curated phage genomes)."""
    output:
        genomes="data/phage/inphared/genomes.fasta.gz",
        metadata="data/phage/inphared/metadata.tsv"
    params:
        genomes_url=config["phage_databases"]["inphared"]["genomes_url"],
        metadata_url=config["phage_databases"]["inphared"]["metadata_url"]
    log:
        "logs/download/phage_inphared.log"
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
        genomes="data/phage/{db}/genomes.fasta.gz"
    output:
        chunks=expand(
            "data/phage/{{db}}/chunks/chunk_{chunk}.fasta",
            chunk=range(1, 101)  # 100 chunks
        )
    params:
        seqs_per_chunk=1000
    log:
        "logs/phage/split_{db}.log"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/split_fasta.py"


rule predict_orfs_prodigal:
    """Predict ORFs from phage genomes using Prodigal (metagenomics mode)."""
    input:
        chunk="data/phage/{db}/chunks/chunk_{chunk}.fasta"
    output:
        proteins="results/phage/{db}/orfs/chunk_{chunk}.faa",
        genes="results/phage/{db}/orfs/chunk_{chunk}.fna",
        gff="results/phage/{db}/orfs/chunk_{chunk}.gff"
    log:
        "logs/phage/prodigal_{db}_chunk_{chunk}.log"
    conda:
        "../envs/orfs.yaml"
    threads: 1
    resources:
        runtime=60,
        mem_mb=4000
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
            "results/phage/{{db}}/orfs/chunk_{chunk}.faa",
            chunk=range(1, 101)
        )
    output:
        merged="results/phage/{db}/predicted_orfs.faa.gz"
    log:
        "logs/phage/merge_orfs_{db}.log"
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
        fasta="results/phage/{db}/predicted_orfs.faa.gz"
    output:
        motifs="results/phage/{db}/gp_motifs/all_gp_motifs.tsv.gz",
        sequences="results/phage/{db}/gp_motifs/all_gp_sequences.fasta.gz"
    params:
        upstream=30,
        downstream=15
    log:
        "logs/phage/extract_gp_{db}.log"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/extract_gp_motifs.py"


rule annotate_phage_domains:
    """Run domain annotation on phage ORFs with GP motifs."""
    input:
        sequences="results/phage/{db}/gp_motifs/all_gp_sequences.fasta.gz"
    output:
        annotations="results/phage/{db}/gp_motifs/domain_annotations.tsv.gz"
    log:
        "logs/phage/annotate_domains_{db}.log"
    conda:
        "../envs/python.yaml"
    threads: 8
    resources:
        runtime=480,
        mem_mb=16000
    script:
        "../scripts/annotate_gp_domains.py"


rule filter_phage_interdomain_gp:
    """Filter phage GP motifs by domain boundaries."""
    input:
        motifs="results/phage/{db}/gp_motifs/all_gp_motifs.tsv.gz",
        annotations="results/phage/{db}/gp_motifs/domain_annotations.tsv.gz"
    output:
        interdomain="results/phage/{db}/gp_motifs/interdomain_gp_motifs.tsv.gz",
        intradomain="results/phage/{db}/gp_motifs/intradomain_gp_motifs.tsv.gz",
        statistics="results/phage/{db}/gp_motifs/gp_motif_stats.tsv"
    log:
        "logs/phage/filter_interdomain_{db}.log"
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
        prok_interdomain="results/prokaryotic/gp_motifs/interdomain_gp_motifs.tsv.gz",
        phage_interdomain=expand(
            "results/phage/{db}/gp_motifs/interdomain_gp_motifs.tsv.gz",
            db=config["phage_databases_to_use"]
        )
    output:
        merged="results/combined/interdomain_gp_motifs_all.tsv.gz",
        summary="results/combined/source_summary.tsv"
    log:
        "logs/combined/merge_motifs.log"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/merge_prokaryotic_phage_motifs.py"


rule cluster_combined_motifs:
    """Cluster combined prokaryotic and phage GP motifs."""
    input:
        interdomain="results/combined/interdomain_gp_motifs_all.tsv.gz"
    output:
        clusters="results/combined/clusters/gp_clusters.tsv.gz",
        representatives="results/combined/clusters/cluster_representatives.fasta"
    params:
        identity=config["prokaryotic"]["clustering_identity"],
        coverage=config["prokaryotic"]["clustering_coverage"]
    log:
        "logs/combined/cluster_motifs.log"
    conda:
        "../envs/python.yaml"
    threads: 12
    resources:
        runtime=240,
        mem_mb=32000
    script:
        "../scripts/cluster_gp_motifs.py"


# ============================================================================
# Phage-Specific Analysis
# ============================================================================

rule analyze_phage_gp_distribution:
    """Analyze distribution of GP motifs in phage genomes."""
    input:
        motifs=expand(
            "results/phage/{db}/gp_motifs/all_gp_motifs.tsv.gz",
            db=config["phage_databases_to_use"]
        ),
        metadata="data/phage/inphared/metadata.tsv"
    output:
        distribution="results/phage/analysis/gp_distribution.tsv",
        plots=directory("results/phage/analysis/distribution_plots/")
    log:
        "logs/phage/analyze_distribution.log"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/analyze_phage_gp_distribution.py"


rule compare_phage_vs_bacterial:
    """Compare GP motif characteristics between phages and bacteria."""
    input:
        phage_motifs="results/combined/interdomain_gp_motifs_all.tsv.gz",
        clusters="results/combined/clusters/gp_clusters.tsv.gz"
    output:
        comparison="results/combined/analysis/phage_vs_bacterial.tsv",
        plots=directory("results/combined/analysis/comparison_plots/")
    log:
        "logs/combined/compare_phage_bacterial.log"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/compare_phage_vs_bacterial.py"
