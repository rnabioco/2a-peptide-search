"""
Rules for searching protein databases with HMM models.
"""

def get_database_file(wildcards):
    """Map database name to its file path."""
    db_map = {
        "uniprot": "uniprot_sprot.fasta.gz",
        "reference_proteomes": "reference_proteomes.fasta.gz",
        "uniparc": "uniparc_active.fasta.gz",
        "mgnify": "mgnify_proteins.fasta.gz"
    }
    return DATA_DIR + f"/{wildcards.database}/{db_map[wildcards.database]}"


rule hmmsearch:
    """Search protein database with HMM model."""
    input:
        hmm=RESULTS_DIR + "/models/{iteration}/2A-{peptide_class}.hmm",
        db=get_database_file
    output:
        hmmsearch=SCRATCH_DIR + "/searches/{database}/{iteration}/2A-{peptide_class}.hmmsearch.gz",
        tblout=SCRATCH_DIR + "/searches/{database}/{iteration}/2A-{peptide_class}.tblout.gz",
        alignment=SCRATCH_DIR + "/searches/{database}/{iteration}/2A-{peptide_class}.sto.gz"
    log:
        "logs/hmmsearch/{database}_{iteration}_{peptide_class}.log"
    conda:
        "../envs/hmmer.yaml"
    threads: 12
    resources:
        runtime=1440,  # 24 hours max
        mem_mb=16000   # 16GB memory
    shell:
        """
        hmmsearch --cpu {threads} \
            --tblout >(gzip > {output.tblout}) \
            -A >(gzip > {output.alignment}) \
            --noali \
            {input.hmm} {input.db} 2> {log} | gzip > {output.hmmsearch}
        """


rule filter_alignment:
    """Filter alignment by E-value threshold."""
    input:
        alignment=SCRATCH_DIR + "/searches/{database}/{iteration}/2A-{peptide_class}.sto.gz",
        tblout=SCRATCH_DIR + "/searches/{database}/{iteration}/2A-{peptide_class}.tblout.gz"
    output:
        filtered=SCRATCH_DIR + "/alignments/{database}/{iteration}/2A-{peptide_class}.filtered.sto"
    params:
        evalue=config["thresholds"]["evalue"]
    log:
        "logs/filter/{database}_{iteration}_{peptide_class}.log"
    conda:
        "../envs/python.yaml"
    resources:
        runtime=30,
        mem_mb=8000
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
            SCRATCH_DIR + "/alignments/{database}/{{iteration}}/2A-{{peptide_class}}.filtered.sto",
            database=config["databases_to_search"]
        )
    output:
        merged=SCRATCH_DIR + "/alignments/{iteration}/2A-{peptide_class}.merged.sto"
    log:
        "logs/merge/{iteration}_{peptide_class}.log"
    conda:
        "../envs/python.yaml"
    resources:
        runtime=60,
        mem_mb=16000
    shell:
        """
        python workflow/scripts/merge_alignments.py \
            {input.alignments} \
            {output.merged} 2> {log}
        """
