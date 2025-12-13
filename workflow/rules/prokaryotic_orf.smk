"""
ORF prediction rules for prokaryotic 2A-like discovery pipeline.

Uses Prodigal to predict protein-coding genes from phage genomes.
"""


def get_phage_genome_file(wildcards):
    """Get phage genome file path (all compressed as of 2025)."""
    phage_db = wildcards.phage_db
    return DATA_DIR + f"/prokaryotic/phage_genomes/{phage_db}_genomes.fasta.gz"


rule run_prodigal:
    """Run Prodigal to predict ORFs from phage genomes.

    Uses metagenomic mode for diverse phage sequences.
    All inputs are gzipped as of 2025.
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
        phage_db="inphared",  # millardlab disabled - URL no longer available
    threads: 1
    shell:
        """
        zcat "{input.genomes}" | \
        prodigal -a /dev/stdout \
            -f gff \
            -o "{output.genes}" \
            -p {params.mode} \
            -g {params.table} \
            2> "{log}" | \
        gzip > "{output.proteins}"
        """
