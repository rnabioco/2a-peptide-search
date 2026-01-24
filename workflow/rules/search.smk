"""
Rules for searching protein databases with HMM models.
"""


def get_database_file(wildcards):
    """Map database name to its file path."""
    # Check if database has a local_path configured (e.g., IMG/VR)
    if (
        wildcards.database in config["databases"]
        and "local_path" in config["databases"][wildcards.database]
    ):
        local_path = config["databases"][wildcards.database]["local_path"]
        if local_path:
            return local_path

    # Standard databases with download rules
    db_map = {
        "uniprot": "uniprot_sprot.fasta.gz",
        "reference_proteomes": "reference_proteomes.fasta.gz",
        "uniparc": "uniparc_active.fasta.gz",
        "mgnify": "mgnify_proteins.fasta.gz",
        "imgvr": "IMGVR_all_proteins.faa.gz",  # Fallback if no local_path set
    }
    return DATA_DIR + f"/{wildcards.database}/{db_map[wildcards.database]}"


def get_model_path(wildcards):
    """Map search iteration to model path.

    Iteration naming:
    - seed: initial search with seed models
    - iter1: search with models refined from seed results
    - iter2: search with models refined from iter1 results
    """
    if wildcards.iteration == "seed":
        model_dir = "seed"
    elif wildcards.iteration == "iter1":
        # iter1 searches use models built from seed results
        model_dir = "seed_refined"
    elif wildcards.iteration == "iter2":
        # iter2 searches use models built from iter1 results
        model_dir = "iter1_refined"
    else:
        # Default to iteration name
        model_dir = wildcards.iteration

    return RESULTS_DIR + f"/models/{model_dir}/2A-{wildcards.peptide_class}.hmm"


rule hmmsearch:
    """Search protein database with HMM model (excludes mgnify - handled separately)."""
    input:
        hmm=get_model_path,
        db=get_database_file,
    output:
        hmmsearch=SCRATCH_DIR
        + "/searches/{database}/{iteration}/2A-{peptide_class}.hmmsearch.gz",
        tblout=SCRATCH_DIR
        + "/searches/{database}/{iteration}/2A-{peptide_class}.tblout.gz",
        alignment=SCRATCH_DIR
        + "/searches/{database}/{iteration}/2A-{peptide_class}.sto.gz",
    log:
        LOGS_DIR + "/hmmsearch/{database}_{iteration}_{peptide_class}.log",
    wildcard_constraints:
        database="(?!mgnify).*",  # Exclude mgnify - handled by split rules below
    threads: 12
    resources:
        runtime=1440,  # 24 hours max
        mem_mb=16000,  # 16GB memory
    shell:
        """
        hmmsearch --cpu {threads} \
            --tblout >(gzip > {output.tblout}) \
            -A >(gzip > {output.alignment}) \
            --noali \
            {input.hmm} {input.db} 2> {log} | gzip > {output.hmmsearch}
        """


# =============================================================================
# MGnify Split-Based Search (for parallel processing of large database)
# =============================================================================


rule hmmsearch_mgnify_split:
    """Search one MGnify split file with HMM model."""
    input:
        hmm=get_model_path,
        db=DATA_DIR + "/mgnify/splits/mgy_proteins_{split_num}.fa.gz",
    output:
        tblout=SCRATCH_DIR
        + "/searches/mgnify_splits/{iteration}/split_{split_num}/2A-{peptide_class}.tblout.gz",
        alignment=SCRATCH_DIR
        + "/searches/mgnify_splits/{iteration}/split_{split_num}/2A-{peptide_class}.sto.gz",
    log:
        LOGS_DIR + "/hmmsearch/mgnify_split_{split_num}_{iteration}_{peptide_class}.log",
    threads: 4
    resources:
        runtime=480,  # 8 hours per split
        mem_mb=8000,
    shell:
        """
        hmmsearch --cpu {threads} \
            --tblout >(gzip > {output.tblout}) \
            -A >(gzip > {output.alignment}) \
            --noali \
            {input.hmm} {input.db} 2> {log} > /dev/null
        """


rule merge_mgnify_searches:
    """Merge MGnify split search results into combined output."""
    input:
        tblouts=expand(
            SCRATCH_DIR
            + "/searches/mgnify_splits/{{iteration}}/split_{split_num}/2A-{{peptide_class}}.tblout.gz",
            split_num=range(1, config["databases"]["mgnify"]["num_splits"] + 1),
        ),
        alignments=expand(
            SCRATCH_DIR
            + "/searches/mgnify_splits/{{iteration}}/split_{split_num}/2A-{{peptide_class}}.sto.gz",
            split_num=range(1, config["databases"]["mgnify"]["num_splits"] + 1),
        ),
    output:
        hmmsearch=SCRATCH_DIR
        + "/searches/mgnify/{iteration}/2A-{peptide_class}.hmmsearch.gz",
        tblout=SCRATCH_DIR
        + "/searches/mgnify/{iteration}/2A-{peptide_class}.tblout.gz",
        alignment=SCRATCH_DIR
        + "/searches/mgnify/{iteration}/2A-{peptide_class}.sto.gz",
    log:
        LOGS_DIR + "/hmmsearch/mgnify_merge_{iteration}_{peptide_class}.log",
    resources:
        runtime=60,
        mem_mb=16000,
    shell:
        """
        # Merge tblout files (concatenate, keeping single header)
        (zcat {input.tblouts[0]} | head -3; \
         for f in {input.tblouts}; do zcat "$f" | tail -n +4; done) | \
         gzip > {output.tblout} 2> {log}

        # Merge Stockholm alignments using esl-alimerge
        tmp_list=$(mktemp)
        trap "rm -f $tmp_list" EXIT

        for aln in {input.alignments}; do
            # Decompress to temp file for esl-alimerge
            tmp_aln=$(mktemp --suffix=.sto)
            zcat "$aln" > "$tmp_aln"
            echo "$tmp_aln" >> $tmp_list
        done

        esl-alimerge --list $tmp_list 2>> {log} | gzip > {output.alignment}

        # Clean up temp alignment files
        while read tmp_aln; do rm -f "$tmp_aln"; done < $tmp_list

        # Create empty hmmsearch output (not needed but keeps file structure consistent)
        echo "# MGnify search results merged from splits" | gzip > {output.hmmsearch}
        """


rule filter_alignment:
    """Filter alignment by E-value threshold."""
    input:
        alignment=SCRATCH_DIR
        + "/searches/{database}/{iteration}/2A-{peptide_class}.sto.gz",
        tblout=SCRATCH_DIR
        + "/searches/{database}/{iteration}/2A-{peptide_class}.tblout.gz",
    output:
        filtered=SCRATCH_DIR
        + "/alignments/{database}/{iteration}/2A-{peptide_class}.filtered.sto",
    params:
        evalue=config["thresholds"]["evalue"],
    log:
        LOGS_DIR + "/filter/{database}_{iteration}_{peptide_class}.log",
    resources:
        runtime=30,
        mem_mb=8000,
    shell:
        """
        python workflow/scripts/filter_alignment.py \
            <(zcat {input.alignment}) \
            <(zcat {input.tblout}) \
            {output.filtered} \
            --evalue {params.evalue} 2> {log}
        """


rule merge_database_alignments:
    """Merge alignments from all databases for a given iteration."""
    input:
        alignments=expand(
            SCRATCH_DIR
            + "/alignments/{database}/{{iteration}}/2A-{{peptide_class}}.filtered.sto",
            database=config["databases_to_search"],
        ),
    output:
        merged=SCRATCH_DIR + "/alignments/{iteration}/2A-{peptide_class}.merged.sto",
    log:
        LOGS_DIR + "/merge/{iteration}_{peptide_class}.log",
    resources:
        runtime=60,
        mem_mb=16000,
    shell:
        """
        # Count number of input alignments
        n_files=$(echo {input.alignments} | wc -w)

        if [[ "$n_files" -eq 1 ]]; then
            # Single file - just copy it (esl-alimerge requires RF annotation)
            cp {input.alignments} {output.merged}
            echo "Single alignment file, copied directly" > {log}
        else
            # Multiple files - merge with esl-alimerge
            # Create temp file with alignment list
            tmp_list=$(mktemp)
            trap "rm -f $tmp_list" EXIT

            # Write alignment paths to temp file
            for aln in {input.alignments}; do
                echo "$aln" >> $tmp_list
            done

            # Merge alignments
            esl-alimerge --list $tmp_list > {output.merged} 2> {log}
        fi
        """
