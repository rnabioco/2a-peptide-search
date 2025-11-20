"""
Rules for downloading protein databases.
"""

rule download_uniprot:
    """Download UniProt database."""
    output:
        fasta="data/uniprot/uniprot_sprot.fasta.gz"
    params:
        url=config["databases"]["uniprot"]["url"]
    log:
        "logs/download/uniprot.log"
    shell:
        """
        mkdir -p data/uniprot
        wget -c -o {log} {params.url} -O {output.fasta}
        """


rule download_reference_proteomes:
    """Download Reference Proteomes database."""
    output:
        tarball="data/reference_proteomes/Reference_Proteomes.tar.gz"
    params:
        url=config["databases"]["reference_proteomes"]["url"]
    log:
        "logs/download/reference_proteomes.log"
    shell:
        """
        mkdir -p data/reference_proteomes
        wget -c -o {log} {params.url} -O {output.tarball}
        """


rule extract_reference_proteomes:
    """Extract and concatenate reference proteomes."""
    input:
        tarball="data/reference_proteomes/Reference_Proteomes.tar.gz"
    output:
        fasta="data/reference_proteomes/reference_proteomes.fasta.gz"
    log:
        "logs/download/extract_reference_proteomes.log"
    shell:
        """
        tar -xzf {input.tarball} -C data/reference_proteomes/ 2> {log}
        find data/reference_proteomes -name "*.fasta.gz" -exec zcat {{}} \; | gzip > {output.fasta}
        """


rule download_uniparc:
    """Download UniParc database (warning: very large)."""
    output:
        fasta="data/uniparc/uniparc_active.fasta.gz"
    params:
        url=config["databases"]["uniparc"]["url"]
    log:
        "logs/download/uniparc.log"
    threads: 1
    shell:
        """
        mkdir -p data/uniparc
        wget -c -o {log} {params.url} -O {output.fasta}
        """


rule download_mgnify:
    """Download MGnify protein database."""
    output:
        fasta="data/mgnify/mgnify_proteins.fasta.gz"
    params:
        url=config["databases"]["mgnify"]["url"]
    log:
        "logs/download/mgnify.log"
    shell:
        """
        mkdir -p data/mgnify
        wget -c -o {log} {params.url} -O {output.fasta}
        """
