"""
Rules for domain annotation using InterProScan.

InterProScan is optional - set interproscan.enabled: true in config.yaml
and ensure interproscan.sh is available in PATH or specify interproscan.path.
"""


def interproscan_enabled():
    """Check if InterProScan is enabled in config."""
    return config.get("interproscan", {}).get("enabled", False)


rule extract_hit_sequences:
    """Extract protein sequences from HMM search hits for domain annotation."""
    input:
        alignment=RESULTS_DIR + "/alignments/final/2A-{peptide_class}.curated.sto",
        tblout=expand(
            SCRATCH_DIR + "/searches/{database}/final/2A-{{peptide_class}}.tblout.gz",
            database=config["databases_to_search"],
        ),
    output:
        fasta=RESULTS_DIR + "/domains/hit_sequences.{peptide_class}.fasta",
        ids=RESULTS_DIR + "/domains/hit_ids.{peptide_class}.txt",
    log:
        LOGS_DIR + "/domains/extract_hits_{peptide_class}.log",
    resources:
        runtime=30,
        mem_mb=8000,
    script:
        "../scripts/extract_hit_sequences.py"


if interproscan_enabled():

    rule run_interproscan:
        """Run InterProScan on hit sequences for domain annotation."""
        input:
            fasta=RESULTS_DIR + "/domains/hit_sequences.{peptide_class}.fasta",
        output:
            tsv=RESULTS_DIR + "/domains/interproscan.{peptide_class}.tsv.gz",
        params:
            analyses=",".join(config.get("interproscan", {}).get("analyses", ["Pfam"])),
            interproscan_path=config.get("interproscan", {}).get("path") or "interproscan.sh",
        log:
            LOGS_DIR + "/domains/interproscan_{peptide_class}.log",
        threads: config.get("interproscan", {}).get("threads", 16)
        resources:
            runtime=1440,  # InterProScan can take a long time
            mem_mb=32000,
        shell:
            """
            {params.interproscan_path} \
                -i {input.fasta} \
                -f TSV \
                -appl {params.analyses} \
                --cpu {threads} \
                -o /dev/stdout 2> {log} | gzip > {output.tsv}
            """

    rule parse_interproscan:
        """Parse InterProScan TSV output into standardized domain format."""
        input:
            tsv=RESULTS_DIR + "/domains/interproscan.{peptide_class}.tsv.gz",
        output:
            domains=RESULTS_DIR + "/domains/parsed_domains.{peptide_class}.tsv",
        log:
            LOGS_DIR + "/domains/parse_interproscan_{peptide_class}.log",
        script:
            "../scripts/parse_interproscan.py"


rule generate_domain_tblout:
    """Generate domain hits table from HMM search results.

    This provides 2A peptide domain coordinates from our HMM models,
    complementing InterProScan's external database annotations.
    """
    input:
        tblouts=expand(
            SCRATCH_DIR + "/searches/{database}/final/2A-{{peptide_class}}.tblout.gz",
            database=config["databases_to_search"],
        ),
    output:
        domains=RESULTS_DIR + "/domains/2a_domains.{peptide_class}.tsv",
    log:
        LOGS_DIR + "/domains/generate_domain_tblout_{peptide_class}.log",
    script:
        "../scripts/generate_domain_tblout.py"
