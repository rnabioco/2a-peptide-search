"""
Rules for generating final analysis report.
"""

rule collect_statistics:
    """Collect statistics from all searches."""
    input:
        tblouts=expand(
            SCRATCH_DIR + "/searches/{database}/final/2A-{peptide_class}.tblout.gz",
            database=config["databases_to_search"],
            peptide_class=["class-1", "class-2"]
        )
    output:
        stats=RESULTS_DIR + "/reports/search_statistics.tsv"
    log:
        LOGS_DIR + "/report/collect_stats.log"
    script:
        "../scripts/collect_statistics.py"


rule generate_report:
    """Generate final Quarto report."""
    input:
        stats=RESULTS_DIR + "/reports/search_statistics.tsv",
        models=expand(
            RESULTS_DIR + "/models/final/2A-{peptide_class}.hmm",
            peptide_class=["class-1", "class-2"]
        ),
        alignments=expand(
            RESULTS_DIR + "/alignments/final/2A-{peptide_class}.curated.sto",
            peptide_class=["class-1", "class-2"]
        )
    output:
        report=RESULTS_DIR + "/reports/2A-peptide-analysis.html"
    log:
        LOGS_DIR + "/report/generate_report.log"
    shell:
        """
        quarto render workflow/report.qmd -o {output.report} 2> {log}
        """
