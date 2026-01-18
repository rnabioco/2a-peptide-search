"""
Prokaryotic 2A-like discovery - APPROACH 2: Unbiased GP Motif Discovery

Comprehensive discovery of GP-containing sequences without prior assumptions
about motif structure. Focuses on inter-domain GP motifs as potential
ribosomal stalling sites.

Strategy:
1. Extract all GP-containing sequences from prokaryotic proteomes
2. Annotate with domain boundaries (Pfam/InterPro)
3. Focus on inter-domain GP motifs
4. Cluster by sequence context
5. Build initial consensus and HMMs
6. Iterative refinement
"""

import glob as pyglob
import re


# ============================================================================
# Aggregation functions for dynamic cluster discovery
# ============================================================================


def get_cluster_ids_from_logos(wildcards):
    """Discover cluster IDs from logo checkpoint output."""
    checkpoint_output = checkpoints.analyze_cluster_conservation.get(
        database=wildcards.database
    ).output.logos_dir
    logo_files = pyglob.glob(f"{checkpoint_output}/cluster_*.png")
    cluster_ids = sorted(
        [int(re.search(r"cluster_(\d+)", f).group(1)) for f in logo_files]
    )
    return cluster_ids


def get_cluster_ids_from_alignments(wildcards):
    """Discover cluster IDs from consensus alignment checkpoint output."""
    checkpoint_output = checkpoints.identify_consensus_patterns.get(
        database=wildcards.database
    ).output.alignments_dir
    sto_files = pyglob.glob(f"{checkpoint_output}/cluster_*.sto")
    cluster_ids = sorted(
        [int(re.search(r"cluster_(\d+)", f).group(1)) for f in sto_files]
    )
    return cluster_ids


def aggregate_cluster_hmms(wildcards):
    """Aggregate all HMM files for a database based on discovered clusters."""
    cluster_ids = get_cluster_ids_from_alignments(wildcards)
    return expand(
        RESULTS_DIR + "/prokaryotic/gp_analysis/{database}/models/cluster_{cluster_id}.hmm",
        database=wildcards.database,
        cluster_id=cluster_ids,
    )


def aggregate_cluster_searches(wildcards):
    """Aggregate all HMM search results for a database."""
    cluster_ids = get_cluster_ids_from_alignments(wildcards)
    return expand(
        RESULTS_DIR + "/prokaryotic/gp_analysis/{database}/hmm_searches/cluster_{cluster_id}.tblout.gz",
        database=wildcards.database,
        cluster_id=cluster_ids,
    )


rule extract_gp_motifs:
    """Extract all sequences containing GP motifs with context."""
    input:
        fasta=get_prokaryotic_database_file,
    output:
        motifs=RESULTS_DIR + "/prokaryotic/gp_analysis/{database}/all_gp_motifs.tsv.gz",
        sequences=RESULTS_DIR + "/prokaryotic/gp_analysis/{database}/all_gp_sequences.fasta.gz",
    params:
        upstream=30,  # residues upstream of GP
        downstream=15,  # residues downstream of GP
    log:
        LOGS_DIR + "/prokaryotic/extract_gp_motifs_{database}.log",
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


checkpoint split_sequences_for_hmmscan:
    """Split GP sequences into chunks for parallel hmmscan."""
    input:
        sequences=RESULTS_DIR + "/prokaryotic/gp_analysis/{database}/all_gp_sequences.fasta.gz",
    output:
        chunk_dir=directory(RESULTS_DIR + "/prokaryotic/gp_analysis/{database}/chunks/"),
    params:
        n_chunks=config.get("prokaryotic", {}).get("hmmscan_chunks", 4),
    log:
        LOGS_DIR + "/prokaryotic/split_sequences_{database}.log",
    shell:
        """
        mkdir -p "{output.chunk_dir}"

        # Split FASTA into N parts using seqkit (preserves record boundaries)
        seqkit split2 \
            --by-part {params.n_chunks} \
            --out-dir "{output.chunk_dir}" \
            "{input.sequences}" \
            2>> "{log}"

        # Rename to consistent chunk_N.fasta.gz format
        cd "{output.chunk_dir}"
        for f in *.fasta.gz *.fa.gz 2>/dev/null; do
            if [[ -f "$f" ]]; then
                # Extract part number from seqkit naming (e.g., .part_001.fasta.gz)
                part=$(echo "$f" | grep -oP 'part_\\K\\d+')
                if [[ -n "$part" ]]; then
                    mv "$f" "chunk_${{part}}.fasta.gz"
                fi
            fi
        done

        echo "Split into {params.n_chunks} chunks" >> "{log}"
        """


rule run_hmmscan_chunk:
    """Run hmmscan on a single sequence chunk."""
    input:
        sequences=RESULTS_DIR + "/prokaryotic/gp_analysis/{database}/chunks/chunk_{chunk_id}.fasta.gz",
        pfam_db=DATA_DIR + "/pfam/Pfam-A.hmm",
        h3f=DATA_DIR + "/pfam/Pfam-A.hmm.h3f",
        h3i=DATA_DIR + "/pfam/Pfam-A.hmm.h3i",
        h3m=DATA_DIR + "/pfam/Pfam-A.hmm.h3m",
        h3p=DATA_DIR + "/pfam/Pfam-A.hmm.h3p",
    output:
        domtblout=RESULTS_DIR + "/prokaryotic/gp_analysis/{database}/chunks/chunk_{chunk_id}.domtblout",
    log:
        LOGS_DIR + "/prokaryotic/hmmscan_{database}_chunk_{chunk_id}.log",
    threads: 4
    resources:
        runtime=240,
        mem_mb=8000,
    shell:
        """
        # Copy Pfam DB to local scratch if available (much faster on clusters)
        if [[ -n "${{TMPDIR:-}}" && -d "$TMPDIR" ]]; then
            echo "Copying Pfam database to local scratch..." >> "{log}"
            cp "{input.pfam_db}" "{input.h3f}" "{input.h3i}" "{input.h3m}" "{input.h3p}" "$TMPDIR/"
            PFAM_DB="$TMPDIR/Pfam-A.hmm"
        else
            PFAM_DB="{input.pfam_db}"
        fi

        # Run hmmscan with --noali (skip alignment output for speed)
        zcat "{input.sequences}" | hmmscan \
            --cpu {threads} \
            --domtblout "{output.domtblout}" \
            --noali \
            --cut_ga \
            "$PFAM_DB" \
            - \
            > /dev/null 2>> "{log}"
        """


rule merge_hmmscan_results:
    """Merge domtblout results from all chunks."""
    input:
        chunks=aggregate_hmmscan_chunks,
    output:
        domtblout=RESULTS_DIR + "/prokaryotic/gp_analysis/{database}/domains.domtblout",
    log:
        LOGS_DIR + "/prokaryotic/merge_hmmscan_{database}.log",
    shell:
        """
        # Combine all chunk domtblout files
        # Keep header from first file, skip headers (lines starting with #) from rest
        first=true
        for chunk in {input.chunks}; do
            if [ "$first" = true ]; then
                cat "$chunk"
                first=false
            else
                grep -v '^#' "$chunk" || true
            fi
        done > "{output.domtblout}"

        echo "Merged $(echo '{input.chunks}' | wc -w) chunk files" >> "{log}"
        """


rule parse_domain_annotations:
    """Parse hmmscan output into domain annotations table."""
    input:
        domtblout=RESULTS_DIR + "/prokaryotic/gp_analysis/{database}/domains.domtblout",
        motifs=RESULTS_DIR + "/prokaryotic/gp_analysis/{database}/all_gp_motifs.tsv.gz",
    output:
        annotations=RESULTS_DIR + "/prokaryotic/gp_analysis/{database}/domain_annotations.tsv.gz",
    log:
        LOGS_DIR + "/prokaryotic/parse_annotations_{database}.log",
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
        motifs=RESULTS_DIR + "/prokaryotic/gp_analysis/{database}/all_gp_motifs.tsv.gz",
        annotations=RESULTS_DIR + "/prokaryotic/gp_analysis/{database}/domain_annotations.tsv.gz",
    output:
        interdomain=RESULTS_DIR + "/prokaryotic/gp_analysis/{database}/interdomain_gp_motifs.tsv.gz",
        intradomain=RESULTS_DIR + "/prokaryotic/gp_analysis/{database}/intradomain_gp_motifs.tsv.gz",
        statistics=RESULTS_DIR + "/prokaryotic/gp_analysis/{database}/gp_motif_stats.tsv",
    log:
        LOGS_DIR + "/prokaryotic/filter_interdomain_gp_{database}.log",
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
    """Extract GP context sequences for clustering.

    Uses ALL GP motifs (not just inter-domain) for comprehensive discovery.
    Domain position info is used as annotation, not a filter.
    """
    input:
        all_gp=RESULTS_DIR + "/prokaryotic/gp_analysis/{database}/all_gp_motifs.tsv.gz",
    output:
        fasta=RESULTS_DIR + "/prokaryotic/gp_analysis/{database}/gp_motifs.fasta",
    log:
        LOGS_DIR + "/prokaryotic/prepare_clustering_{database}.log",
    shell:
        """
        # Extract context sequences to FASTA (using all GP motifs)
        zcat "{input.all_gp}" | awk -F'\t' 'NR>1 {{
            printf ">%s_GP%s_pos%s\\n%s\\n", $1, $4, $5, $6
        }}' > "{output.fasta}" 2> "{log}"
        """


rule run_mmseqs_clustering:
    """Cluster GP motifs using MMseqs2."""
    input:
        fasta=RESULTS_DIR + "/prokaryotic/gp_analysis/{database}/gp_motifs.fasta",
    output:
        cluster_tsv=RESULTS_DIR + "/prokaryotic/gp_analysis/{database}/mmseqs_clusters.tsv",
    params:
        identity=config["prokaryotic"]["clustering_identity"],
        coverage=config["prokaryotic"]["clustering_coverage"],
        prefix=lambda w: RESULTS_DIR + f"/prokaryotic/gp_analysis/{w.database}/mmseqs",
        tmpdir=lambda w: RESULTS_DIR + f"/prokaryotic/gp_analysis/{w.database}/tmp",
    log:
        LOGS_DIR + "/prokaryotic/mmseqs_cluster_{database}.log",
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

        # Cleanup temp files (careful not to delete the output TSV)
        rm -rf "{params.tmpdir}"
        rm -f "{params.prefix}_db" "{params.prefix}_db".* "{params.prefix}_db_h" "{params.prefix}_db_h".*
        rm -f "{params.prefix}_cluster" "{params.prefix}_cluster".*
        """


rule parse_gp_clusters:
    """Parse MMseqs2 clustering results."""
    input:
        cluster_tsv=RESULTS_DIR + "/prokaryotic/gp_analysis/{database}/mmseqs_clusters.tsv",
        all_gp=RESULTS_DIR + "/prokaryotic/gp_analysis/{database}/all_gp_motifs.tsv.gz",
    output:
        clusters=RESULTS_DIR + "/prokaryotic/gp_analysis/{database}/gp_clusters.tsv.gz",
        representatives=RESULTS_DIR
        + "/prokaryotic/gp_analysis/{database}/cluster_representatives.fasta",
    log:
        LOGS_DIR + "/prokaryotic/parse_clusters_{database}.log",
    shell:
        """
        python workflow/scripts/parse_gp_clusters.py \
            --cluster-tsv "{input.cluster_tsv}" \
            --motifs "{input.all_gp}" \
            --clusters-out "{output.clusters}" \
            --representatives-out "{output.representatives}" \
            > "{log}" 2>&1
        """


checkpoint analyze_cluster_conservation:
    """Analyze conservation patterns within each cluster.

    This is a checkpoint because the number of clusters is dynamic.
    Downstream rules use get_cluster_ids_from_logos() to discover outputs.
    """
    input:
        clusters=RESULTS_DIR + "/prokaryotic/gp_analysis/{database}/gp_clusters.tsv.gz",
        motifs=RESULTS_DIR + "/prokaryotic/gp_analysis/{database}/all_gp_motifs.tsv.gz",
    output:
        conservation=RESULTS_DIR + "/prokaryotic/gp_analysis/{database}/cluster_conservation.tsv.gz",
        logos_dir=directory(RESULTS_DIR + "/prokaryotic/gp_analysis/{database}/logos"),
    params:
        min_cluster_size=config["prokaryotic"]["min_cluster_size"],
    log:
        LOGS_DIR + "/prokaryotic/analyze_conservation_{database}.log",
    shell:
        """
        mkdir -p "{output.logos_dir}"
        python workflow/scripts/analyze_gp_conservation.py \
            --clusters "{input.clusters}" \
            --motifs "{input.motifs}" \
            --conservation "{output.conservation}" \
            --logos-dir "{output.logos_dir}" \
            --min-cluster-size {params.min_cluster_size} \
            > "{log}" 2>&1
        """


checkpoint identify_consensus_patterns:
    """Identify consensus patterns for major clusters.

    This is a checkpoint because the number of clusters meeting quality
    thresholds is dynamic. Downstream rules use get_cluster_ids_from_alignments()
    to discover which cluster alignments were produced.
    """
    input:
        conservation=RESULTS_DIR + "/prokaryotic/gp_analysis/{database}/cluster_conservation.tsv.gz",
        clusters=RESULTS_DIR + "/prokaryotic/gp_analysis/{database}/gp_clusters.tsv.gz",
        motifs=RESULTS_DIR + "/prokaryotic/gp_analysis/{database}/all_gp_motifs.tsv.gz",
    output:
        consensus=RESULTS_DIR + "/prokaryotic/gp_analysis/{database}/consensus_patterns.tsv",
        alignments_dir=directory(RESULTS_DIR + "/prokaryotic/gp_analysis/{database}/consensus"),
    params:
        min_conservation=config["prokaryotic"]["min_conservation"],
        min_cluster_size=config["prokaryotic"]["min_cluster_size"],
    log:
        LOGS_DIR + "/prokaryotic/identify_consensus_{database}.log",
    shell:
        """
        mkdir -p "{output.alignments_dir}"
        python workflow/scripts/identify_consensus_patterns.py \
            --conservation "{input.conservation}" \
            --clusters "{input.clusters}" \
            --motifs "{input.motifs}" \
            --consensus "{output.consensus}" \
            --alignments-dir "{output.alignments_dir}" \
            --min-conservation {params.min_conservation} \
            --min-cluster-size {params.min_cluster_size} \
            > "{log}" 2>&1
        """


# ============================================================================
# HMM Building for Prokaryotic Motifs
# ============================================================================


rule build_prokaryotic_hmms:
    """Build HMMs from consensus patterns."""
    input:
        alignment=RESULTS_DIR + "/prokaryotic/gp_analysis/{database}/consensus/cluster_{cluster_id}.sto",
    output:
        hmm=RESULTS_DIR + "/prokaryotic/gp_analysis/{database}/models/cluster_{cluster_id}.hmm",
    params:
        name=lambda w: f"prok-2A-{w.database}-cluster-{w.cluster_id}",
    log:
        LOGS_DIR + "/prokaryotic/build_hmm_{database}_cluster_{cluster_id}.log",
    shell:
        """
        hmmbuild -n {params.name} "{output.hmm}" "{input.alignment}" 2> "{log}"
        """


rule search_prokaryotic_proteomes:
    """Search prokaryotic proteomes with discovered HMMs."""
    input:
        hmm=RESULTS_DIR + "/prokaryotic/gp_analysis/{database}/models/cluster_{cluster_id}.hmm",
        db=get_prokaryotic_database_file,
    output:
        hmmsearch=RESULTS_DIR
        + "/prokaryotic/gp_analysis/{database}/hmm_searches/cluster_{cluster_id}.hmmsearch.gz",
        tblout=RESULTS_DIR + "/prokaryotic/gp_analysis/{database}/hmm_searches/cluster_{cluster_id}.tblout.gz",
        alignment=RESULTS_DIR + "/prokaryotic/gp_analysis/{database}/hmm_searches/cluster_{cluster_id}.sto.gz",
    log:
        LOGS_DIR + "/prokaryotic/search_{database}_cluster_{cluster_id}.log",
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
# Aggregate Target Rules (dynamic cluster discovery)
# ============================================================================


rule build_all_prokaryotic_hmms:
    """Build HMMs for all discovered clusters in a database.

    This rule aggregates over the dynamically discovered cluster alignments
    from the identify_consensus_patterns checkpoint.
    """
    input:
        aggregate_cluster_hmms,
    output:
        touch(RESULTS_DIR + "/prokaryotic/gp_analysis/{database}/.hmms_built"),


rule search_all_prokaryotic_clusters:
    """Run HMM searches for all discovered clusters in a database.

    This rule aggregates over all discovered clusters and runs hmmsearch
    for each one against the source database.
    """
    input:
        aggregate_cluster_searches,
    output:
        touch(RESULTS_DIR + "/prokaryotic/gp_analysis/{database}/.searches_complete"),


rule prokaryotic_gp_analysis:
    """Complete GP-based prokaryotic analysis for a database.

    This is the main target rule for the GP motif discovery approach.
    """
    input:
        hmms_built=RESULTS_DIR + "/prokaryotic/gp_analysis/{database}/.hmms_built",
        searches_complete=RESULTS_DIR + "/prokaryotic/gp_analysis/{database}/.searches_complete",
    output:
        touch(RESULTS_DIR + "/prokaryotic/gp_analysis/{database}/.analysis_complete"),
