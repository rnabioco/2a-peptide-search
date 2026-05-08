"""
Rules for downloading protein databases.
"""


rule download_uniprot:
    """Download UniProt database."""
    output:
        fasta=DATA_DIR + "/uniprot/uniprot_sprot.fasta.gz",
    params:
        url=config["databases"]["uniprot"]["url"],
    log:
        LOGS_DIR + "/download/uniprot.log",
    shell:
        """
        mkdir -p $(dirname {output.fasta})
        wget -c -o {log} {params.url} -O {output.fasta}
        """


rule download_uniprot_relnotes:
    """Capture UniProt release version for provenance (e.g. Pfam submission)."""
    output:
        relnotes=DATA_DIR + "/uniprot_release.txt",
    log:
        LOGS_DIR + "/download/uniprot_relnotes.log",
    shell:
        """
        mkdir -p $(dirname {output.relnotes})
        wget -c -o {log} https://ftp.uniprot.org/pub/databases/uniprot/current_release/relnotes.txt -O {output.relnotes}
        """


rule download_reference_proteomes:
    """Download Reference Proteomes database."""
    output:
        tarball=DATA_DIR + "/reference_proteomes/Reference_Proteomes.tar.gz",
    params:
        url=config["databases"]["reference_proteomes"]["url"],
    log:
        LOGS_DIR + "/download/reference_proteomes.log",
    shell:
        """
        mkdir -p $(dirname {output.tarball})
        wget -c -o {log} {params.url} -O {output.tarball}
        """


rule extract_reference_proteomes:
    """Extract and concatenate reference proteomes."""
    input:
        tarball=DATA_DIR + "/reference_proteomes/Reference_Proteomes.tar.gz",
    output:
        fasta=DATA_DIR + "/reference_proteomes/reference_proteomes.fasta.gz",
    log:
        LOGS_DIR + "/download/extract_reference_proteomes.log",
    shell:
        r"""
        tar -xzf {input.tarball} -C $(dirname {input.tarball})/ 2> {log}
        find $(dirname {output.fasta}) -name "*.fasta.gz" -exec zcat {{}} \; | gzip > {output.fasta}
        """


rule download_uniparc:
    """Download UniParc database (warning: very large)."""
    output:
        fasta=DATA_DIR + "/uniparc/uniparc_active.fasta.gz",
    params:
        url=config["databases"]["uniparc"]["url"],
    log:
        LOGS_DIR + "/download/uniparc.log",
    threads: 1
    shell:
        """
        mkdir -p $(dirname {output.fasta})
        wget -c -o {log} {params.url} -O {output.fasta}
        """


rule download_mgnify_split:
    """Download one split file from MGnify protein database."""
    output:
        fasta=DATA_DIR + "/mgnify/splits/mgy_proteins_{split_num}.fa.gz",
    params:
        base_url=config["databases"]["mgnify"]["base_url"],
    log:
        LOGS_DIR + "/download/mgnify_split_{split_num}.log",
    shell:
        """
        mkdir -p $(dirname {output.fasta})
        wget -c -o {log} {params.base_url}/mgy_proteins_{wildcards.split_num}.fa.gz -O {output.fasta}
        """


rule download_mgnify:
    """Download all MGnify protein database splits."""
    input:
        splits=expand(
            DATA_DIR + "/mgnify/splits/mgy_proteins_{split_num}.fa.gz",
            split_num=range(1, config["databases"]["mgnify"]["num_splits"] + 1),
        ),
    output:
        flag=DATA_DIR + "/mgnify/download_complete.flag",
    shell:
        """
        touch {output.flag}
        """


rule merge_mgnify_splits:
    """Merge MGnify split files into single database (optional, for convenience)."""
    input:
        splits=expand(
            DATA_DIR + "/mgnify/splits/mgy_proteins_{split_num}.fa.gz",
            split_num=range(1, config["databases"]["mgnify"]["num_splits"] + 1),
        ),
    output:
        fasta=DATA_DIR + "/mgnify/mgnify_proteins.fasta.gz",
    log:
        LOGS_DIR + "/download/merge_mgnify.log",
    shell:
        """
        cat {input.splits} > {output.fasta} 2> {log}
        """
