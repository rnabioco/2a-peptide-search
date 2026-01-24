"""
Prokaryotic 2A-like discovery - C-terminal Arrest Motif Analysis

Extract and analyze sequences with known bacterial arrest motifs (RAGP, QAPP, etc.)
near their C-terminus. This approach directly targets known stalling motif patterns
for validation against characterized arrest peptides.

Strategy:
1. Extract sequences with C-terminal arrest motifs (RAGP, RAPG, QAPP, etc.)
2. Cluster by motif type
3. Build HMMs from each motif type cluster
4. Validate against known stalling peptides (SecM, TnaC, etc.)
"""


# ============================================================================
# Arrest Motif Extraction and Analysis
# ============================================================================


rule extract_arrest_motifs:
    """Extract sequences with C-terminal arrest motifs."""
    input:
        fasta=get_prokaryotic_database_file,
    output:
        motifs=RESULTS_DIR + "/prokaryotic/arrest_analysis/{database}/arrest_motifs.tsv.gz",
        sequences=RESULTS_DIR + "/prokaryotic/arrest_analysis/{database}/arrest_sequences.fasta.gz",
    params:
        upstream=40,  # residues upstream of motif
        max_cterm_distance=50,  # max distance from C-terminus
    log:
        LOGS_DIR + "/prokaryotic/extract_arrest_motifs_{database}.log",
    shell:
        """
        python workflow/scripts/extract_arrest_motifs.py \
            --fasta "{input.fasta}" \
            --motifs-out "{output.motifs}" \
            --sequences-out "{output.sequences}" \
            --upstream {params.upstream} \
            --max-cterm-distance {params.max_cterm_distance} \
            > "{log}" 2>&1
        """


rule cluster_arrest_by_motif_type:
    """Group arrest sequences by motif type and create alignments."""
    input:
        motifs=RESULTS_DIR + "/prokaryotic/arrest_analysis/{database}/arrest_motifs.tsv.gz",
    output:
        alignments_dir=directory(RESULTS_DIR + "/prokaryotic/arrest_analysis/{database}/alignments"),
        summary=RESULTS_DIR + "/prokaryotic/arrest_analysis/{database}/motif_summary.tsv",
    log:
        LOGS_DIR + "/prokaryotic/cluster_arrest_{database}.log",
    shell:
        """
        python workflow/scripts/cluster_arrest_motifs.py \
            --motifs "{input.motifs}" \
            --alignments-dir "{output.alignments_dir}" \
            --summary "{output.summary}" \
            > "{log}" 2>&1
        """


rule build_arrest_hmms:
    """Build HMMs from arrest motif alignments."""
    input:
        alignments_dir=RESULTS_DIR + "/prokaryotic/arrest_analysis/{database}/alignments",
    output:
        hmm_dir=directory(RESULTS_DIR + "/prokaryotic/arrest_analysis/{database}/models"),
        hmm_db=RESULTS_DIR + "/prokaryotic/arrest_analysis/{database}/arrest_motifs.hmm",
    log:
        LOGS_DIR + "/prokaryotic/build_arrest_hmms_{database}.log",
    shell:
        """
        mkdir -p "{output.hmm_dir}"

        # Build individual HMMs for each motif type
        for sto in "{input.alignments_dir}"/*.sto; do
            if [[ -f "$sto" ]]; then
                motif=$(basename "$sto" .sto)
                hmmbuild -n "$motif" "{output.hmm_dir}/$motif.hmm" "$sto" 2>> "{log}"
            fi
        done

        # Concatenate all HMMs into a single database
        cat "{output.hmm_dir}"/*.hmm > "{output.hmm_db}"

        # Press the database for hmmscan
        hmmpress -f "{output.hmm_db}" 2>> "{log}"

        echo "Built HMM database from $(ls "{output.hmm_dir}"/*.hmm 2>/dev/null | wc -l) motif types" >> "{log}"
        """


rule validate_arrest_against_known:
    """Validate arrest motif HMMs against known stalling peptides."""
    input:
        hmm_db=RESULTS_DIR + "/prokaryotic/arrest_analysis/{database}/arrest_motifs.hmm",
        known_peptides="resources/stalling-peptides/known_stalling_peptides.fasta",
    output:
        validation=RESULTS_DIR + "/prokaryotic/arrest_analysis/{database}/known_peptide_validation.tsv",
        summary=RESULTS_DIR + "/prokaryotic/arrest_analysis/{database}/known_peptide_summary.txt",
    log:
        LOGS_DIR + "/prokaryotic/validate_arrest_{database}.log",
    shell:
        """
        python workflow/scripts/validate_arrest_motifs.py \
            --hmm-db "{input.hmm_db}" \
            --known-peptides "{input.known_peptides}" \
            --output "{output.validation}" \
            --summary "{output.summary}" \
            > "{log}" 2>&1
        """


rule arrest_motif_analysis:
    """Complete arrest motif analysis target."""
    input:
        validation=RESULTS_DIR + "/prokaryotic/arrest_analysis/{database}/known_peptide_validation.tsv",
        summary=RESULTS_DIR + "/prokaryotic/arrest_analysis/{database}/known_peptide_summary.txt",
    output:
        touch(RESULTS_DIR + "/prokaryotic/arrest_analysis/{database}/.analysis_complete"),


# ============================================================================
# Large Database Search with Arrest HMMs
# ============================================================================
#
# These rules search larger databases (UniParc, MGnify) using the arrest motif
# HMMs built from the bacteria database. This enables discovery of novel arrest
# peptides in metagenomic and comprehensive protein sequence archives.
#
# The arrest HMMs are located at:
#   - Combined: results/prokaryotic/arrest_analysis/bacteria/arrest_motifs.hmm
#   - Individual: results/prokaryotic/arrest_analysis/bacteria/models/{motif_type}.hmm
#
# Available motif types: HAPP, HGPP, QAPP, QGPP, RAGP, RAPG, RAPP, RPPP
# ============================================================================


# Define the arrest motif types discovered from bacteria
ARREST_MOTIF_TYPES = ["HAPP", "HGPP", "QAPP", "QGPP", "RAGP", "RAPG", "RAPP", "RPPP"]

# Large databases to search with arrest HMMs
ARREST_SEARCH_DATABASES = {
    "uniparc": DATA_DIR + "/uniparc/uniparc_active.fasta.gz",
    "mgnify": DATA_DIR + "/mgnify/mgnify_proteins.fasta.gz",
}


def get_arrest_search_database(wildcards):
    """Map database name to its file path for arrest HMM searches."""
    if wildcards.target_db in ARREST_SEARCH_DATABASES:
        return ARREST_SEARCH_DATABASES[wildcards.target_db]
    # Fallback to standard prokaryotic database lookup
    return get_prokaryotic_database_file(wildcards)


rule search_with_arrest_hmms:
    """Search a target database with a single arrest motif HMM.

    This rule searches large databases (UniParc, MGnify) using individual
    arrest motif HMMs built from bacterial proteins. Each motif type
    (RAGP, QAPP, etc.) is searched separately for fine-grained results.

    Input HMMs are from the bacteria arrest analysis pipeline.
    """
    input:
        hmm=RESULTS_DIR + "/prokaryotic/arrest_analysis/bacteria/models/{motif_type}.hmm",
        db=get_arrest_search_database,
    output:
        hmmsearch=RESULTS_DIR + "/prokaryotic/arrest_search/{target_db}/{motif_type}.hmmsearch.gz",
        tblout=RESULTS_DIR + "/prokaryotic/arrest_search/{target_db}/{motif_type}.tblout.gz",
        alignment=RESULTS_DIR + "/prokaryotic/arrest_search/{target_db}/{motif_type}.sto.gz",
    log:
        LOGS_DIR + "/prokaryotic/arrest_search_{target_db}_{motif_type}.log",
    threads: 12
    resources:
        runtime=1440,  # 24 hours for very large databases
        mem_mb=32000,
    shell:
        """
        hmmsearch --cpu {threads} \
            --tblout >(gzip > "{output.tblout}") \
            -A >(gzip > "{output.alignment}") \
            --noali \
            "{input.hmm}" "{input.db}" 2> "{log}" | gzip > "{output.hmmsearch}"
        """


rule search_with_all_arrest_hmms:
    """Search a target database with all arrest motif HMMs at once.

    Uses the combined HMM database for efficiency when searching very large
    databases. This is faster than running individual searches but provides
    less granular output.
    """
    input:
        hmm_db=RESULTS_DIR + "/prokaryotic/arrest_analysis/bacteria/arrest_motifs.hmm",
        db=get_arrest_search_database,
    output:
        hmmsearch=RESULTS_DIR + "/prokaryotic/arrest_search/{target_db}/all_motifs.hmmsearch.gz",
        tblout=RESULTS_DIR + "/prokaryotic/arrest_search/{target_db}/all_motifs.tblout.gz",
        domtblout=RESULTS_DIR + "/prokaryotic/arrest_search/{target_db}/all_motifs.domtblout.gz",
        alignment=RESULTS_DIR + "/prokaryotic/arrest_search/{target_db}/all_motifs.sto.gz",
    log:
        LOGS_DIR + "/prokaryotic/arrest_search_{target_db}_all.log",
    threads: 12
    resources:
        runtime=1440,  # 24 hours for very large databases
        mem_mb=32000,
    shell:
        """
        hmmsearch --cpu {threads} \
            --tblout >(gzip > "{output.tblout}") \
            --domtblout >(gzip > "{output.domtblout}") \
            -A >(gzip > "{output.alignment}") \
            --noali \
            "{input.hmm_db}" "{input.db}" 2> "{log}" | gzip > "{output.hmmsearch}"
        """


# ============================================================================
# MGnify Split Search (for parallel processing on clusters)
# ============================================================================
#
# MGnify is split into 25 files (~270GB total). For cluster environments,
# searching each split separately enables massive parallelization.
# ============================================================================


rule search_mgnify_split_with_arrest_hmms:
    """Search a single MGnify split with all arrest HMMs.

    MGnify is divided into 25 split files. This rule enables parallel
    searching across all splits on cluster environments.
    """
    input:
        hmm_db=RESULTS_DIR + "/prokaryotic/arrest_analysis/bacteria/arrest_motifs.hmm",
        split=DATA_DIR + "/mgnify/splits/mgy_proteins_{split_num}.fa.gz",
    output:
        hmmsearch=RESULTS_DIR + "/prokaryotic/arrest_search/mgnify_splits/split_{split_num}.hmmsearch.gz",
        tblout=RESULTS_DIR + "/prokaryotic/arrest_search/mgnify_splits/split_{split_num}.tblout.gz",
        alignment=RESULTS_DIR + "/prokaryotic/arrest_search/mgnify_splits/split_{split_num}.sto.gz",
    log:
        LOGS_DIR + "/prokaryotic/arrest_search_mgnify_split_{split_num}.log",
    threads: 12
    resources:
        runtime=720,  # 12 hours per split
        mem_mb=16000,
    shell:
        """
        hmmsearch --cpu {threads} \
            --tblout >(gzip > "{output.tblout}") \
            -A >(gzip > "{output.alignment}") \
            --noali \
            "{input.hmm_db}" "{input.split}" 2> "{log}" | gzip > "{output.hmmsearch}"
        """


rule merge_mgnify_arrest_search_results:
    """Merge arrest HMM search results from all MGnify splits."""
    input:
        tblouts=expand(
            RESULTS_DIR + "/prokaryotic/arrest_search/mgnify_splits/split_{split_num}.tblout.gz",
            split_num=range(1, config["databases"]["mgnify"]["num_splits"] + 1),
        ),
        alignments=expand(
            RESULTS_DIR + "/prokaryotic/arrest_search/mgnify_splits/split_{split_num}.sto.gz",
            split_num=range(1, config["databases"]["mgnify"]["num_splits"] + 1),
        ),
    output:
        tblout=RESULTS_DIR + "/prokaryotic/arrest_search/mgnify/merged.tblout.gz",
        alignment=RESULTS_DIR + "/prokaryotic/arrest_search/mgnify/merged.sto.gz",
    log:
        LOGS_DIR + "/prokaryotic/merge_mgnify_arrest_search.log",
    shell:
        """
        # Merge tblout files (keep header from first, skip headers from rest)
        (
            zcat "{input.tblouts[0]}" | head -3
            for f in {input.tblouts}; do
                zcat "$f" | grep -v '^#' || true
            done
        ) | gzip > "{output.tblout}" 2> "{log}"

        # Merge Stockholm alignments using esl-alimerge
        # First decompress all alignments to temp files
        tmpdir=$(mktemp -d)
        trap "rm -rf $tmpdir" EXIT

        for f in {input.alignments}; do
            base=$(basename "$f" .sto.gz)
            zcat "$f" > "$tmpdir/$base.sto"
        done

        # Merge alignments (esl-alimerge handles Stockholm format)
        esl-alimerge --outformat stockholm "$tmpdir"/*.sto 2>> "{log}" | gzip > "{output.alignment}"
        """


# ============================================================================
# Aggregation Rules
# ============================================================================


rule aggregate_arrest_searches_by_motif:
    """Aggregate all search results for a single motif type across databases."""
    input:
        tblouts=expand(
            RESULTS_DIR + "/prokaryotic/arrest_search/{target_db}/{{motif_type}}.tblout.gz",
            target_db=ARREST_SEARCH_DATABASES.keys(),
        ),
        alignments=expand(
            RESULTS_DIR + "/prokaryotic/arrest_search/{target_db}/{{motif_type}}.sto.gz",
            target_db=ARREST_SEARCH_DATABASES.keys(),
        ),
    output:
        tblout=RESULTS_DIR + "/prokaryotic/arrest_search/aggregated/{motif_type}.tblout.gz",
        alignment=RESULTS_DIR + "/prokaryotic/arrest_search/aggregated/{motif_type}.sto.gz",
        summary=RESULTS_DIR + "/prokaryotic/arrest_search/aggregated/{motif_type}.summary.tsv",
    log:
        LOGS_DIR + "/prokaryotic/aggregate_arrest_{motif_type}.log",
    shell:
        """
        # Merge tblout files
        (
            zcat "{input.tblouts[0]}" | head -3
            for f in {input.tblouts}; do
                zcat "$f" | grep -v '^#' || true
            done
        ) | gzip > "{output.tblout}" 2> "{log}"

        # Merge Stockholm alignments
        tmpdir=$(mktemp -d)
        trap "rm -rf $tmpdir" EXIT

        for f in {input.alignments}; do
            base=$(basename "$f" .sto.gz)
            zcat "$f" > "$tmpdir/$base.sto"
        done

        esl-alimerge --outformat stockholm "$tmpdir"/*.sto 2>> "{log}" | gzip > "{output.alignment}"

        # Generate summary statistics
        echo -e "motif_type\\tdatabase\\tn_hits\\tn_unique_sequences" > "{output.summary}"
        for f in {input.tblouts}; do
            db=$(basename $(dirname "$f"))
            n_hits=$(zcat "$f" | grep -v '^#' | wc -l || echo 0)
            n_unique=$(zcat "$f" | grep -v '^#' | awk '{{print $1}}' | sort -u | wc -l || echo 0)
            echo -e "{wildcards.motif_type}\\t$db\\t$n_hits\\t$n_unique"
        done >> "{output.summary}" 2>> "{log}"
        """


rule aggregate_all_arrest_searches:
    """Aggregate search results across all motif types and databases."""
    input:
        summaries=expand(
            RESULTS_DIR + "/prokaryotic/arrest_search/aggregated/{motif_type}.summary.tsv",
            motif_type=ARREST_MOTIF_TYPES,
        ),
        alignments=expand(
            RESULTS_DIR + "/prokaryotic/arrest_search/aggregated/{motif_type}.sto.gz",
            motif_type=ARREST_MOTIF_TYPES,
        ),
    output:
        summary=RESULTS_DIR + "/prokaryotic/arrest_search/arrest_search_summary.tsv",
        flag=touch(RESULTS_DIR + "/prokaryotic/arrest_search/.searches_complete"),
    log:
        LOGS_DIR + "/prokaryotic/aggregate_all_arrest_searches.log",
    shell:
        """
        # Combine all summary files
        head -1 "{input.summaries[0]}" > "{output.summary}"
        for f in {input.summaries}; do
            tail -n +2 "$f" >> "{output.summary}"
        done

        echo "Aggregated results from {input.summaries}" >> "{log}"
        """


# ============================================================================
# Target Rules for Arrest HMM Searches
# ============================================================================


rule search_uniparc_with_arrest_hmms:
    """Search UniParc with all arrest motif HMMs."""
    input:
        expand(
            RESULTS_DIR + "/prokaryotic/arrest_search/uniparc/{motif_type}.tblout.gz",
            motif_type=ARREST_MOTIF_TYPES,
        ),
    output:
        touch(RESULTS_DIR + "/prokaryotic/arrest_search/uniparc/.search_complete"),


rule search_mgnify_with_arrest_hmms:
    """Search MGnify with all arrest motif HMMs (via split files)."""
    input:
        RESULTS_DIR + "/prokaryotic/arrest_search/mgnify/merged.tblout.gz",
        RESULTS_DIR + "/prokaryotic/arrest_search/mgnify/merged.sto.gz",
    output:
        touch(RESULTS_DIR + "/prokaryotic/arrest_search/mgnify/.search_complete"),


rule search_all_large_databases_with_arrest_hmms:
    """Search all large databases with arrest HMMs.

    This is the main target for comprehensive arrest peptide discovery.
    Searches both UniParc and MGnify with all arrest motif types.
    """
    input:
        RESULTS_DIR + "/prokaryotic/arrest_search/uniparc/.search_complete",
        RESULTS_DIR + "/prokaryotic/arrest_search/mgnify/.search_complete",
        RESULTS_DIR + "/prokaryotic/arrest_search/.searches_complete",
    output:
        touch(RESULTS_DIR + "/prokaryotic/arrest_search/.all_searches_complete"),


# ============================================================================
# Cross-Database Arrest Motif Searches
# ============================================================================
#
# Search one database with HMMs built from another database to:
# 1. Find phage arrest peptides using bacterial HMMs
# 2. Identify horizontal transfer candidates (viral HMMs vs bacteria)
# 3. Discover divergent homologs with relaxed E-values
# ============================================================================


# Define cross-database search pairs
# Each tuple is (source_db, target_db)
CROSS_SEARCH_PAIRS = [
    # Bacterial HMMs searching viral databases
    ("bacteria", "uniprot_viruses"),
    ("bacteria", "ncbi_viral_refseq"),
    ("bacteria", "gpd"),
    ("bacteria", "imgvr"),
    # Viral HMMs searching bacterial database
    ("uniprot_viruses", "bacteria"),
    ("ncbi_viral_refseq", "bacteria"),
]


rule cross_database_arrest_search:
    """Search one database with HMMs from another database.

    This enables discovery of divergent arrest peptides by applying
    HMMs built from one domain of life to search another. For example,
    bacterial arrest motif HMMs can identify potential phage-encoded
    arrest peptides that manipulate host ribosomes.
    """
    input:
        hmm_db=RESULTS_DIR + "/prokaryotic/arrest_analysis/{source_db}/arrest_motifs.hmm",
        target_db=get_prokaryotic_database_file,
    output:
        tblout=RESULTS_DIR + "/prokaryotic/arrest_cross_search/{source_db}_vs_{database}.tblout.gz",
        alignment=RESULTS_DIR + "/prokaryotic/arrest_cross_search/{source_db}_vs_{database}.sto.gz",
    log:
        LOGS_DIR + "/prokaryotic/arrest_cross_search_{source_db}_vs_{database}.log",
    threads: 12
    resources:
        runtime=720,  # 12 hours
        mem_mb=16000,
    params:
        evalue=10,  # Relaxed E-value to capture divergent sequences
    shell:
        """
        # Check if HMM database exists and has content
        if [[ ! -s "{input.hmm_db}" ]]; then
            echo "Warning: HMM database {input.hmm_db} is empty or missing" >> "{log}"
            # Create empty output files
            echo "# No HMM models available" | gzip > "{output.tblout}"
            echo "# STOCKHOLM 1.0" | gzip > "{output.alignment}"
            exit 0
        fi

        hmmsearch --cpu {threads} \
            -E {params.evalue} \
            --tblout >(gzip > "{output.tblout}") \
            -A >(gzip > "{output.alignment}") \
            --noali \
            "{input.hmm_db}" "{input.target_db}" 2> "{log}" || true

        # Ensure output files exist even if no hits
        [[ -f "{output.tblout}" ]] || echo "# No hits" | gzip > "{output.tblout}"
        [[ -f "{output.alignment}" ]] || echo "# STOCKHOLM 1.0" | gzip > "{output.alignment}"
        """


rule cross_database_arrest_search_relaxed:
    """Search with very relaxed E-value to capture distant homologs.

    Use E-value of 100 to capture even very divergent sequences,
    then filter by domain architecture and conservation.
    """
    input:
        hmm_db=RESULTS_DIR + "/prokaryotic/arrest_analysis/{source_db}/arrest_motifs.hmm",
        target_db=get_prokaryotic_database_file,
    output:
        tblout=RESULTS_DIR + "/prokaryotic/arrest_cross_search/{source_db}_vs_{database}.relaxed.tblout.gz",
        alignment=RESULTS_DIR + "/prokaryotic/arrest_cross_search/{source_db}_vs_{database}.relaxed.sto.gz",
    log:
        LOGS_DIR + "/prokaryotic/arrest_cross_search_{source_db}_vs_{database}_relaxed.log",
    threads: 12
    resources:
        runtime=720,
        mem_mb=16000,
    params:
        evalue=100,  # Very relaxed E-value
    shell:
        """
        if [[ ! -s "{input.hmm_db}" ]]; then
            echo "Warning: HMM database {input.hmm_db} is empty or missing" >> "{log}"
            echo "# No HMM models available" | gzip > "{output.tblout}"
            echo "# STOCKHOLM 1.0" | gzip > "{output.alignment}"
            exit 0
        fi

        hmmsearch --cpu {threads} \
            -E {params.evalue} \
            --tblout >(gzip > "{output.tblout}") \
            -A >(gzip > "{output.alignment}") \
            --noali \
            "{input.hmm_db}" "{input.target_db}" 2> "{log}" || true

        [[ -f "{output.tblout}" ]] || echo "# No hits" | gzip > "{output.tblout}"
        [[ -f "{output.alignment}" ]] || echo "# STOCKHOLM 1.0" | gzip > "{output.alignment}"
        """


rule analyze_cross_search_results:
    """Analyze cross-database search results and generate summary."""
    input:
        tblout=RESULTS_DIR + "/prokaryotic/arrest_cross_search/{source_db}_vs_{database}.tblout.gz",
        alignment=RESULTS_DIR + "/prokaryotic/arrest_cross_search/{source_db}_vs_{database}.sto.gz",
    output:
        summary=RESULTS_DIR + "/prokaryotic/arrest_cross_search/{source_db}_vs_{database}.summary.tsv",
    log:
        LOGS_DIR + "/prokaryotic/analyze_cross_search_{source_db}_vs_{database}.log",
    shell:
        """
        # Generate basic summary statistics
        n_hits=$(zcat "{input.tblout}" | grep -v '^#' | wc -l || echo 0)
        n_unique=$(zcat "{input.tblout}" | grep -v '^#' | awk '{{print $1}}' | sort -u | wc -l || echo 0)

        echo -e "source_db\\ttarget_db\\tn_hits\\tn_unique_sequences" > "{output.summary}"
        echo -e "{wildcards.source_db}\\t{wildcards.database}\\t$n_hits\\t$n_unique" >> "{output.summary}"

        echo "Cross-search {wildcards.source_db} vs {wildcards.database}: $n_hits hits, $n_unique unique" >> "{log}"
        """


rule aggregate_cross_search_results:
    """Aggregate all cross-database search results."""
    input:
        summaries=[
            RESULTS_DIR + f"/prokaryotic/arrest_cross_search/{source}_vs_{target}.summary.tsv"
            for source, target in CROSS_SEARCH_PAIRS
        ],
    output:
        summary=RESULTS_DIR + "/prokaryotic/arrest_cross_search/cross_search_summary.tsv",
    log:
        LOGS_DIR + "/prokaryotic/aggregate_cross_search.log",
    run:
        import pandas as pd

        dfs = []
        for f in input.summaries:
            try:
                df = pd.read_csv(f, sep="\t")
                dfs.append(df)
            except Exception as e:
                print(f"Warning: Could not read {f}: {e}", file=open(log[0], "a"))

        if dfs:
            combined = pd.concat(dfs, ignore_index=True)
            combined.to_csv(output.summary, sep="\t", index=False)
        else:
            pd.DataFrame(columns=["source_db", "target_db", "n_hits", "n_unique_sequences"]).to_csv(
                output.summary, sep="\t", index=False
            )


# ============================================================================
# Iterative HMM Refinement
# ============================================================================


rule refine_arrest_hmms:
    """Refine arrest motif HMMs using hits from all databases.

    Collects hits from bacterial and viral databases, filters by quality,
    deduplicates, and builds refined HMMs with broader sequence diversity.
    """
    input:
        # Original bacterial arrest analysis
        bacteria_alignments=RESULTS_DIR + "/prokaryotic/arrest_analysis/bacteria/alignments",
        bacteria_summary=RESULTS_DIR + "/prokaryotic/arrest_analysis/bacteria/motif_summary.tsv",
        # Cross-search hits (bacteria HMMs vs viral DBs)
        cross_tblouts=expand(
            RESULTS_DIR + "/prokaryotic/arrest_cross_search/bacteria_vs_{target_db}.tblout.gz",
            target_db=["uniprot_viruses", "ncbi_viral_refseq"],
        ),
        cross_alignments=expand(
            RESULTS_DIR + "/prokaryotic/arrest_cross_search/bacteria_vs_{target_db}.sto.gz",
            target_db=["uniprot_viruses", "ncbi_viral_refseq"],
        ),
    output:
        refined_dir=directory(RESULTS_DIR + "/prokaryotic/arrest_analysis/refined/alignments"),
        refined_hmm=RESULTS_DIR + "/prokaryotic/arrest_analysis/refined/arrest_motifs.hmm",
        summary=RESULTS_DIR + "/prokaryotic/arrest_analysis/refined/refinement_summary.tsv",
    log:
        LOGS_DIR + "/prokaryotic/refine_arrest_hmms.log",
    params:
        cross_alignments=lambda w, input: " ".join(f"--cross-alignments {f}" for f in input.cross_alignments),
        cross_tblouts=lambda w, input: " ".join(f"--cross-tblouts {f}" for f in input.cross_tblouts),
    shell:
        """
        python workflow/scripts/refine_arrest_hmms.py \
            --bacteria-alignments "{input.bacteria_alignments}" \
            {params.cross_alignments} \
            {params.cross_tblouts} \
            --output-dir "{output.refined_dir}" \
            --output-hmm "{output.refined_hmm}" \
            --summary "{output.summary}" \
            > "{log}" 2>&1
        """


# ============================================================================
# Distribution Analysis
# ============================================================================


rule analyze_arrest_distribution:
    """Analyze the distribution of arrest motifs across databases.

    Creates summary statistics and visualizations showing:
    - Motif counts per database/organism
    - Enriched/depleted motif types
    - Host-phage comparison
    """
    input:
        summaries=expand(
            RESULTS_DIR + "/prokaryotic/arrest_analysis/{database}/motif_summary.tsv",
            database=["bacteria", "archaea", "uniprot_viruses", "ncbi_viral_refseq"],
        ),
        cross_summary=RESULTS_DIR + "/prokaryotic/arrest_cross_search/cross_search_summary.tsv",
    output:
        matrix=RESULTS_DIR + "/prokaryotic/analysis/arrest_distribution_matrix.tsv",
        host_phage=RESULTS_DIR + "/prokaryotic/analysis/host_phage_comparison.tsv",
        summary=RESULTS_DIR + "/prokaryotic/analysis/arrest_distribution_summary.txt",
    log:
        LOGS_DIR + "/prokaryotic/analyze_arrest_distribution.log",
    params:
        summaries=lambda w, input: " ".join(f"--summaries {f}" for f in input.summaries),
    shell:
        """
        python workflow/scripts/analyze_arrest_distribution.py \
            {params.summaries} \
            --cross-summary "{input.cross_summary}" \
            --matrix-out "{output.matrix}" \
            --host-phage-out "{output.host_phage}" \
            --summary-out "{output.summary}" \
            > "{log}" 2>&1
        """


# ============================================================================
# Target Rules for Cross-Database Analysis
# ============================================================================


rule run_all_cross_searches:
    """Run all cross-database arrest motif searches."""
    input:
        [
            RESULTS_DIR + f"/prokaryotic/arrest_cross_search/{source}_vs_{target}.summary.tsv"
            for source, target in CROSS_SEARCH_PAIRS
        ],
        RESULTS_DIR + "/prokaryotic/arrest_cross_search/cross_search_summary.tsv",
    output:
        touch(RESULTS_DIR + "/prokaryotic/arrest_cross_search/.all_cross_searches_complete"),


rule run_arrest_refinement:
    """Run iterative HMM refinement for arrest motifs."""
    input:
        RESULTS_DIR + "/prokaryotic/arrest_analysis/refined/arrest_motifs.hmm",
    output:
        touch(RESULTS_DIR + "/prokaryotic/arrest_analysis/refined/.refinement_complete"),


rule run_arrest_distribution_analysis:
    """Run distribution analysis for arrest motifs."""
    input:
        RESULTS_DIR + "/prokaryotic/analysis/arrest_distribution_matrix.tsv",
        RESULTS_DIR + "/prokaryotic/analysis/host_phage_comparison.tsv",
    output:
        touch(RESULTS_DIR + "/prokaryotic/analysis/.distribution_analysis_complete"),


rule complete_arrest_cross_analysis:
    """Complete arrest motif cross-database analysis pipeline."""
    input:
        RESULTS_DIR + "/prokaryotic/arrest_cross_search/.all_cross_searches_complete",
        RESULTS_DIR + "/prokaryotic/arrest_analysis/refined/.refinement_complete",
        RESULTS_DIR + "/prokaryotic/analysis/.distribution_analysis_complete",
    output:
        touch(RESULTS_DIR + "/prokaryotic/.arrest_cross_analysis_complete"),
