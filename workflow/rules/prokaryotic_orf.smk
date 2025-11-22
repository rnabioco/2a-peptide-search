"""
ORF prediction rules for prokaryotic 2A-like discovery pipeline.

Uses Prodigal to predict protein-coding genes from phage genomes.
"""


def get_phage_genome_file(wildcards):
    """Get phage genome file path, handling compressed and uncompressed files."""
    phage_db = wildcards.phage_db

    # INPHARED is uncompressed, Millard Lab is compressed
    if phage_db == "inphared":
        return DATA_DIR + f"/prokaryotic/phage_genomes/{phage_db}_genomes.fasta"
    else:  # millardlab
        return DATA_DIR + f"/prokaryotic/phage_genomes/{phage_db}_genomes.fasta.gz"


rule run_prodigal:
    """Run Prodigal to predict ORFs from phage genomes.

    Uses metagenomic mode for diverse phage sequences.
    Handles both compressed and uncompressed input.
    """
    input:
        genomes=get_phage_genome_file,
    output:
        proteins=DATA_DIR + "/prokaryotic/{phage_db}_proteins.fasta.gz",
        genes=DATA_DIR + "/prokaryotic/{phage_db}_genes.gff",
    params:
        mode=config["orf_prediction"]["mode"],
        table=config["orf_prediction"]["translation_table"],
    log:
        LOGS_DIR + "/prodigal/{phage_db}.log",
    wildcard_constraints:
        phage_db="inphared|millardlab",
    threads: 1
    shell:
        """
        # Handle compressed or uncompressed input
        if [[ "{input.genomes}" == *.gz ]]; then
            # Decompress on the fly for compressed input
            zcat "{input.genomes}" | \
            prodigal -a /dev/stdout \
                -f gff \
                -o "{output.genes}" \
                -p {params.mode} \
                -g {params.table} \
                2> "{log}" | \
            gzip > "{output.proteins}"
        else
            # Direct input for uncompressed
            prodigal -i "{input.genomes}" \
                -a >(gzip > "{output.proteins}") \
                -f gff \
                -o "{output.genes}" \
                -p {params.mode} \
                -g {params.table} \
                2> "{log}"
        fi
        """
