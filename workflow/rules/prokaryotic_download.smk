"""
Download rules for prokaryotic 2A-like discovery pipeline.

Handles downloading:
- Prokaryotic protein databases (bacteria, archaea)
- Viral protein databases (NCBI, UniProt, IMG/VR)
- Phage genome databases (INPHARED, Millard Lab)
- Pfam domain database
"""


rule download_prokaryotic_proteomes:
    """Download prokaryotic reference proteomes."""
    output:
        fasta=DATA_DIR + "/prokaryotic/downloads/{database}.fasta.gz",
    params:
        url=lambda w: config["prokaryotic_databases"][w.database]["url"],
    log:
        LOGS_DIR + "/download/prokaryotic_{database}.log",
    shell:
        """
        mkdir -p $(dirname "{output.fasta}")
        wget -c -o "{log}" "{params.url}" -O "{output.fasta}"
        """


rule download_ncbi_viral_split:
    """Download one NCBI Viral RefSeq protein file."""
    output:
        fasta=DATA_DIR
        + "/prokaryotic/ncbi_viral/splits/viral.{file_num}.protein.faa.gz",
    params:
        base_url=lambda w: config["prokaryotic_databases"]["ncbi_viral_refseq"][
            "base_url"
        ],
    log:
        LOGS_DIR + "/download/ncbi_viral_{file_num}.log",
    shell:
        """
        mkdir -p $(dirname "{output.fasta}")
        wget -c -o "{log}" "{params.base_url}/viral.{wildcards.file_num}.protein.faa.gz" -O "{output.fasta}"
        """


rule download_ncbi_viral_refseq:
    """Download all NCBI Viral RefSeq protein files."""
    input:
        splits=expand(
            DATA_DIR
            + "/prokaryotic/ncbi_viral/splits/viral.{file_num}.protein.faa.gz",
            file_num=range(
                1,
                config["prokaryotic_databases"]["ncbi_viral_refseq"]["num_files"] + 1,
            ),
        ),
    output:
        flag=DATA_DIR + "/prokaryotic/ncbi_viral/download_complete.flag",
    shell:
        """
        touch "{output.flag}"
        """


rule merge_ncbi_viral_splits:
    """Merge NCBI Viral RefSeq splits into single database."""
    input:
        splits=expand(
            DATA_DIR
            + "/prokaryotic/ncbi_viral/splits/viral.{file_num}.protein.faa.gz",
            file_num=range(
                1,
                config["prokaryotic_databases"]["ncbi_viral_refseq"]["num_files"] + 1,
            ),
        ),
        flag=DATA_DIR + "/prokaryotic/ncbi_viral/download_complete.flag",
    output:
        fasta=DATA_DIR + "/prokaryotic/ncbi_viral_refseq.fasta.gz",
    log:
        LOGS_DIR + "/download/merge_ncbi_viral.log",
    shell:
        """
        cat {input.splits} > "{output.fasta}" 2> "{log}"
        """


rule download_uniprot_viruses:
    """Download UniProt viral proteins via REST API."""
    output:
        fasta=DATA_DIR + "/prokaryotic/uniprot_viruses.fasta.gz",
    params:
        url=config["prokaryotic_databases"]["uniprot_viruses"]["url"],
    log:
        LOGS_DIR + "/download/uniprot_viruses.log",
    shell:
        """
        mkdir -p $(dirname "{output.fasta}")
        wget -c -o "{log}" "{params.url}" -O "{output.fasta}"
        """


rule download_inphared_genomes:
    """Download INPHARED phage genomes."""
    output:
        genomes=DATA_DIR + "/prokaryotic/phage_genomes/inphared_genomes.fasta",
    params:
        url=config["phage_databases"]["inphared"]["genomes_url"],
    log:
        LOGS_DIR + "/download/inphared_genomes.log",
    shell:
        """
        mkdir -p $(dirname "{output.genomes}")
        wget -c -o "{log}" "{params.url}" -O "{output.genomes}"
        """


rule download_millardlab_genomes:
    """Download Millard Lab phage genomes."""
    output:
        genomes=DATA_DIR + "/prokaryotic/phage_genomes/millardlab_genomes.fasta.gz",
    params:
        url=config["phage_databases"]["millardlab"]["url"],
    log:
        LOGS_DIR + "/download/millardlab_genomes.log",
    shell:
        """
        mkdir -p $(dirname "{output.genomes}")
        wget -c -o "{log}" "{params.url}" -O "{output.genomes}"
        """


rule download_pfam_database:
    """Download and decompress Pfam-A HMM database."""
    output:
        hmm=DATA_DIR + "/pfam/Pfam-A.hmm",
        h3f=DATA_DIR + "/pfam/Pfam-A.hmm.h3f",
        h3i=DATA_DIR + "/pfam/Pfam-A.hmm.h3i",
        h3m=DATA_DIR + "/pfam/Pfam-A.hmm.h3m",
        h3p=DATA_DIR + "/pfam/Pfam-A.hmm.h3p",
    params:
        url=config["prokaryotic_databases"]["pfam"]["url"],
    log:
        LOGS_DIR + "/download/pfam.log",
    shell:
        """
        mkdir -p $(dirname "{output.hmm}")
        wget -c -o "{log}" "{params.url}" -O "{output.hmm}.gz"
        gunzip -f "{output.hmm}.gz"

        # Press HMM database to create binary auxfiles
        hmmpress "{output.hmm}" 2>> "{log}"
        """
