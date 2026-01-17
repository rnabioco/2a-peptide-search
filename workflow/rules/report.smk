"""
Rules for generating final analysis report.
"""


rule collect_statistics:
    """Collect statistics from all searches."""
    input:
        tblouts=expand(
            SCRATCH_DIR + "/searches/{database}/final/2A-{peptide_class}.tblout.gz",
            database=config["databases_to_search"],
            peptide_class=["class-1", "class-2"],
        ),
    output:
        stats=RESULTS_DIR + "/reports/search_statistics.tsv",
    log:
        LOGS_DIR + "/report/collect_stats.log",
    script:
        "../scripts/collect_statistics.py"


rule generate_report:
    """Generate final Quarto report."""
    input:
        stats=RESULTS_DIR + "/reports/search_statistics.tsv",
        models=expand(
            RESULTS_DIR + "/models/final/2A-{peptide_class}.hmm",
            peptide_class=["class-1", "class-2"],
        ),
        alignments=expand(
            RESULTS_DIR + "/alignments/final/2A-{peptide_class}.curated.sto",
            peptide_class=["class-1", "class-2"],
        ),
    output:
        report=RESULTS_DIR + "/reports/2A-peptide-analysis.html",
    log:
        LOGS_DIR + "/report/generate_report.log",
    shell:
        """
        quarto render workflow/report.qmd -o {output.report} 2> {log}
        """


def get_domain_report_inputs(wildcards):
    """Get inputs for domain report based on InterProScan availability."""
    inputs = {
        "domains": RESULTS_DIR + f"/domains/2a_domains.{wildcards.peptide_class}.tsv",
    }
    # Add InterProScan if enabled
    if config.get("interproscan", {}).get("enabled", False):
        inputs["interproscan"] = RESULTS_DIR + f"/domains/interproscan.{wildcards.peptide_class}.tsv.gz"
    return inputs


rule generate_domain_report:
    """Generate domain co-occurrence analysis report."""
    input:
        unpack(get_domain_report_inputs),
    output:
        report=RESULTS_DIR + "/reports/domain-analysis-{peptide_class}.html",
    params:
        results_dir=RESULTS_DIR,
        interproscan_enabled=config.get("interproscan", {}).get("enabled", False),
    log:
        LOGS_DIR + "/report/domain_analysis_{peptide_class}.log",
    resources:
        runtime=30,
        mem_mb=8000,
    shell:
        """
        quarto render workflow/reports/domain-analysis.qmd \
            -P class={wildcards.peptide_class} \
            -P results_dir={params.results_dir} \
            -o {output.report} 2> {log}
        """
