"""
Rules for downloading protein databases.
"""

rule download_uniprot:
    """Download UniProt database."""
    output:
        fasta=DATA_DIR + "/uniprot/uniprot_sprot.fasta.gz"
    params:
        url=config["databases"]["uniprot"]["url"]
    log:
        "logs/download/uniprot.log"
    shell:
        """
        mkdir -p $(dirname {output.fasta})
        wget -c -o {log} {params.url} -O {output.fasta}
        """


rule download_reference_proteomes:
    """Download Reference Proteomes database."""
    output:
        tarball=DATA_DIR + "/reference_proteomes/Reference_Proteomes.tar.gz"
    params:
        url=config["databases"]["reference_proteomes"]["url"]
    log:
        "logs/download/reference_proteomes.log"
    shell:
        """
        mkdir -p $(dirname {output.tarball})
        wget -c -o {log} {params.url} -O {output.tarball}
        """


rule extract_reference_proteomes:
    """Extract and concatenate reference proteomes."""
    input:
        tarball=DATA_DIR + "/reference_proteomes/Reference_Proteomes.tar.gz"
    output:
        fasta=DATA_DIR + "/reference_proteomes/reference_proteomes.fasta.gz"
    log:
        "logs/download/extract_reference_proteomes.log"
    shell:
        r"""
        tar -xzf {input.tarball} -C $(dirname {input.tarball})/ 2> {log}
        find $(dirname {output.fasta}) -name "*.fasta.gz" -exec zcat {{}} \; | gzip > {output.fasta}
        """


rule download_uniparc:
    """Download UniParc database (warning: very large)."""
    output:
        fasta=DATA_DIR + "/uniparc/uniparc_active.fasta.gz"
    params:
        url=config["databases"]["uniparc"]["url"]
    log:
        "logs/download/uniparc.log"
    threads: 1
    shell:
        """
        mkdir -p $(dirname {output.fasta})
        wget -c -o {log} {params.url} -O {output.fasta}
        """


rule download_mgnify:
    """Download MGnify protein database."""
    output:
        fasta=DATA_DIR + "/mgnify/mgnify_proteins.fasta.gz"
    params:
        url=config["databases"]["mgnify"]["url"]
    log:
        "logs/download/mgnify.log"
    shell:
        """
        mkdir -p $(dirname {output.fasta})
        wget -c -o {log} {params.url} -O {output.fasta}
        """
