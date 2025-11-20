"""
Rules for searching protein databases with HMM models.
"""

rule hmmsearch:
    """Search protein database with HMM model."""
    input:
        hmm="results/models/{iteration}/2A-{peptide_class}.hmm",
        db="data/{database}/{database_file}"
    output:
        hmmsearch="results/searches/{database}/{iteration}/2A-{peptide_class}.hmmsearch.gz",
        tblout="results/searches/{database}/{iteration}/2A-{peptide_class}.tblout.gz",
        alignment="results/searches/{database}/{iteration}/2A-{peptide_class}.sto.gz"
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
        alignment="results/searches/{database}/{iteration}/2A-{peptide_class}.sto.gz",
        tblout="results/searches/{database}/{iteration}/2A-{peptide_class}.tblout.gz"
    output:
        filtered="results/alignments/{database}/{iteration}/2A-{peptide_class}.filtered.sto"
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
            "results/alignments/{database}/{{iteration}}/2A-{{peptide_class}}.filtered.sto",
            database=config["databases_to_search"]
        )
    output:
        merged="results/alignments/{iteration}/2A-{peptide_class}.merged.sto"
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
