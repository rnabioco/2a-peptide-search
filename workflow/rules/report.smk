"""
Rules for generating final analysis report.
"""

rule collect_statistics:
    """Collect statistics from all searches."""
    input:
        tblouts=expand(
            "results/searches/{database}/final/2A-{peptide_class}.tblout.gz",
            database=config["databases_to_search"],
            peptide_class=["class-1", "class-2"]
        )
    output:
        stats="results/reports/search_statistics.tsv"
    log:
        "logs/report/collect_stats.log"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/collect_statistics.py"


rule generate_report:
    """Generate final Quarto report."""
    input:
        stats="results/reports/search_statistics.tsv",
        models=expand(
            "results/models/final/2A-{peptide_class}.hmm",
            peptide_class=["class-1", "class-2"]
        ),
        alignments=expand(
            "results/alignments/final/2A-{peptide_class}.curated.sto",
            peptide_class=["class-1", "class-2"]
        )
    output:
        report="results/reports/2A-peptide-analysis.html"
    log:
        "logs/report/generate_report.log"
    conda:
        "../envs/r-quarto.yaml"
    shell:
        """
        quarto render workflow/report.qmd -o {output.report} 2> {log}
        """
