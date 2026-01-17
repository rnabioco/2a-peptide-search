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
# Helper Functions
# ============================================================================


def get_prokaryotic_database_file(wildcards):
    """Map database name to its file path, handling different download mechanisms.

    For databases that may contain invalid IUPAC characters (like IMG/VR),
    automatically use the sanitized version to prevent HMMER parse errors.
    """
    db_config = config["prokaryotic_databases"][wildcards.database]

    # IMG/VR: Always use sanitized version (may contain invalid chars like '-')
    if wildcards.database == "imgvr":
        return DATA_DIR + "/prokaryotic/imgvr.sanitized.fasta.gz"

    # Check if database has a local_path configured
    if "local_path" in db_config and db_config["local_path"]:
        return db_config["local_path"]

    # Databases with direct URL downloads (go to downloads/ subdirectory)
    if wildcards.database in ["bacteria", "archaea", "uniprot_viruses"]:
        return DATA_DIR + f"/prokaryotic/downloads/{wildcards.database}.fasta.gz"

    # Databases from merged splits (go to prokaryotic/ directly)
    if wildcards.database == "ncbi_viral_refseq":
        return DATA_DIR + "/prokaryotic/ncbi_viral_refseq.fasta.gz"

    # Databases from ORF prediction (go to prokaryotic/ directly)
    if wildcards.database in ["inphared_proteins", "millardlab_proteins"]:
        return DATA_DIR + f"/prokaryotic/{wildcards.database}.fasta.gz"

    # Gut Phage Database (CyVerse/iVirus)
    if wildcards.database == "gpd":
        return DATA_DIR + "/prokaryotic/gpd.fasta.gz"

    # MGnify uses main pipeline download
    if wildcards.database == "mgnify":
        return f"{DATA_DIR}/mgnify/mgnify_proteins.fasta.gz"

    # Fallback to downloads subdirectory
    return DATA_DIR + f"/prokaryotic/downloads/{wildcards.database}.fasta.gz"


def aggregate_hmmscan_chunks(wildcards):
    """Dynamically get all chunk domtblout files after checkpoint completes."""
    checkpoint_output = checkpoints.split_sequences_for_hmmscan.get(**wildcards).output[0]
    # Find all chunk files created by the checkpoint
    chunk_ids = glob_wildcards(
        RESULTS_DIR + f"/prokaryotic/gp_analysis/{wildcards.database}/chunks/chunk_{{n}}.fasta.gz"
    ).n
    return expand(
        RESULTS_DIR + "/prokaryotic/gp_analysis/{database}/chunks/chunk_{chunk_id}.domtblout",
        database=wildcards.database,
        chunk_id=chunk_ids,
    )


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
        db=get_prokaryotic_database_file,
    output:
        hmmsearch=RESULTS_DIR + "/prokaryotic/seed_search_results/{database}/comprehensive_hmm/results.hmmsearch.gz",
        tblout=RESULTS_DIR + "/prokaryotic/seed_search_results/{database}/comprehensive_hmm/results.tblout.gz",
        alignment=RESULTS_DIR + "/prokaryotic/seed_search_results/{database}/comprehensive_hmm/results.sto.gz",
    log:
        LOGS_DIR + "/prokaryotic/search_comprehensive_{database}.log",
    wildcard_constraints:
        database="[^/]+",
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
        db=get_prokaryotic_database_file,
    output:
        hmmsearch=RESULTS_DIR + "/prokaryotic/seed_search_results/{database}/{motif}.hmmsearch.gz",
        tblout=RESULTS_DIR + "/prokaryotic/seed_search_results/{database}/{motif}.tblout.gz",
        alignment=RESULTS_DIR + "/prokaryotic/seed_search_results/{database}/{motif}.sto.gz",
    log:
        LOGS_DIR + "/prokaryotic/seed_search_{database}_{motif}.log",
    wildcard_constraints:
        database="[^/]+",
        motif="[^/]+",
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


rule merge_comprehensive_searches:
    """Merge comprehensive search results from all databases."""
    input:
        alignments=expand(
            RESULTS_DIR + "/prokaryotic/seed_search_results/{database}/comprehensive_hmm/results.sto.gz",
            database=config["prokaryotic_databases_to_search"],
        ),
    output:
        merged=RESULTS_DIR + "/prokaryotic/seed_search_results/comprehensive_merged.sto.gz",
    log:
        LOGS_DIR + "/prokaryotic/merge_comprehensive_searches.log",
    shell:
        """
        # Decompress all alignments
        for aln in {input.alignments}; do
            zcat "$aln"
        done | \
        # Merge with esl-alimerge and compress
        esl-alimerge --list - 2> "{log}" | gzip > "{output.merged}"
        """


# ============================================================================
# APPROACH 2: Unbiased GP Motif Extraction
# ============================================================================


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
    """Extract GP context sequences for clustering."""
    input:
        interdomain=RESULTS_DIR + "/prokaryotic/gp_analysis/{database}/interdomain_gp_motifs.tsv.gz",
    output:
        fasta=RESULTS_DIR + "/prokaryotic/gp_analysis/{database}/gp_motifs.fasta",
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

        # Cleanup temp files
        rm -rf "{params.tmpdir}" "{params.prefix}_db"* "{params.prefix}_cluster"*
        """


rule parse_gp_clusters:
    """Parse MMseqs2 clustering results."""
    input:
        cluster_tsv=RESULTS_DIR + "/prokaryotic/gp_analysis/{database}/mmseqs_clusters.tsv",
        interdomain=RESULTS_DIR + "/prokaryotic/gp_analysis/{database}/interdomain_gp_motifs.tsv.gz",
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
            --interdomain "{input.interdomain}" \
            --clusters-out "{output.clusters}" \
            --representatives-out "{output.representatives}" \
            > "{log}" 2>&1
        """


rule analyze_cluster_conservation:
    """Analyze conservation patterns within each cluster."""
    input:
        clusters=RESULTS_DIR + "/prokaryotic/gp_analysis/{database}/gp_clusters.tsv.gz",
        motifs=RESULTS_DIR + "/prokaryotic/gp_analysis/{database}/interdomain_gp_motifs.tsv.gz",
    output:
        conservation=RESULTS_DIR + "/prokaryotic/gp_analysis/{database}/cluster_conservation.tsv.gz",
        logos=expand(
            RESULTS_DIR + "/prokaryotic/gp_analysis/{{database}}/logos/cluster_{cluster_id}.png",
            cluster_id=range(1, 21),  # Top 20 clusters
        ),
    params:
        min_cluster_size=config["prokaryotic"]["min_cluster_size"],
        logos_dir=lambda w: RESULTS_DIR + f"/prokaryotic/gp_analysis/{w.database}/logos",
    log:
        LOGS_DIR + "/prokaryotic/analyze_conservation_{database}.log",
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
        conservation=RESULTS_DIR + "/prokaryotic/gp_analysis/{database}/cluster_conservation.tsv.gz",
        clusters=RESULTS_DIR + "/prokaryotic/gp_analysis/{database}/gp_clusters.tsv.gz",
        motifs=RESULTS_DIR + "/prokaryotic/gp_analysis/{database}/interdomain_gp_motifs.tsv.gz",
    output:
        consensus=RESULTS_DIR + "/prokaryotic/gp_analysis/{database}/consensus_patterns.tsv",
        alignments=expand(
            RESULTS_DIR + "/prokaryotic/gp_analysis/{{database}}/consensus/cluster_{cluster_id}.sto",
            cluster_id=range(1, 21),  # Top 20 clusters
        ),
    params:
        min_conservation=config["prokaryotic"]["min_conservation"],
        min_cluster_size=config["prokaryotic"]["min_cluster_size"],
        alignments_dir=lambda w: RESULTS_DIR + f"/prokaryotic/gp_analysis/{w.database}/consensus",
    log:
        LOGS_DIR + "/prokaryotic/identify_consensus_{database}.log",
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
# Validation Against Known Stalling Peptides
# ============================================================================


rule validate_against_known_peptides:
    """Compare discovered motifs to known stalling peptides (SecM, TnaC, etc.)."""
    input:
        hmms=expand(
            RESULTS_DIR + "/prokaryotic/gp_analysis/bacteria/models/cluster_{cluster_id}.hmm",
            cluster_id=range(1, 21),
        ),
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
    """Compare results from all-GP vs inter-domain-GP approaches."""
    input:
        all_gp_motifs=RESULTS_DIR + "/prokaryotic/gp_analysis/bacteria/all_gp_motifs.tsv.gz",
        interdomain_motifs=RESULTS_DIR
        + "/prokaryotic/gp_analysis/bacteria/interdomain_gp_motifs.tsv.gz",
        clusters=RESULTS_DIR + "/prokaryotic/gp_analysis/bacteria/gp_clusters.tsv.gz",
        validation=RESULTS_DIR + "/prokaryotic/validation/known_peptide_hits.tsv",
        # APPROACH 1: Seed-based searches (merged from all databases)
        seed_comprehensive=RESULTS_DIR
        + "/prokaryotic/seed_search_results/comprehensive_merged.sto.gz",
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
        conservation=RESULTS_DIR + "/prokaryotic/gp_analysis/bacteria/cluster_conservation.tsv.gz",
        consensus=RESULTS_DIR + "/prokaryotic/gp_analysis/bacteria/consensus_patterns.tsv",
    output:
        report=RESULTS_DIR + "/prokaryotic/reports/prokaryotic_discovery.html",
    log:
        LOGS_DIR + "/prokaryotic/generate_report.log",
    script:
        "../scripts/prokaryotic_discovery_report.qmd"
