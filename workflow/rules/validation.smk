"""
Rules for validating HMM models against known literature 2A sequences.

This module benchmarks the trained models against a curated set of
published 2A peptide sequences with known activity measurements.
"""


rule convert_literature_to_fasta:
    """Convert literature sequences TSV to FASTA format for hmmsearch."""
    input:
        tsv="resources/literature/known-2a-sequences.tsv",
    output:
        fasta=RESULTS_DIR + "/validation/literature_sequences.fasta",
    log:
        LOGS_DIR + "/validation/convert_literature.log",
    script:
        "../scripts/convert_literature_fasta.py"


rule validate_literature_sequences:
    """Search literature sequences with trained HMM models."""
    input:
        fasta=RESULTS_DIR + "/validation/literature_sequences.fasta",
        model=RESULTS_DIR + "/models/final/2A-{peptide_class}.hmm",
    output:
        hmmsearch=RESULTS_DIR + "/validation/lit.{peptide_class}.hmmsearch.gz",
        tblout=RESULTS_DIR + "/validation/lit.{peptide_class}.tblout",
        alignment=RESULTS_DIR + "/validation/lit.{peptide_class}.sto",
    log:
        LOGS_DIR + "/validation/hmmsearch_lit_{peptide_class}.log",
    threads: 2
    resources:
        runtime=10,
        mem_mb=2000,
    shell:
        """
        hmmsearch --cpu {threads} \
            --tblout {output.tblout} \
            -A {output.alignment} \
            --noali \
            {input.model} {input.fasta} 2> {log} | gzip > {output.hmmsearch}
        """


rule compare_literature_hits:
    """Compare HMM hits against known activity data."""
    input:
        tblout_c1=RESULTS_DIR + "/validation/lit.class-1.tblout",
        tblout_c2=RESULTS_DIR + "/validation/lit.class-2.tblout",
        sequences="resources/literature/known-2a-sequences.tsv",
    output:
        comparison=RESULTS_DIR + "/validation/literature_comparison.tsv",
    log:
        LOGS_DIR + "/validation/compare_literature.log",
    script:
        "../scripts/compare_literature.py"


rule generate_validation_report:
    """Generate validation summary report."""
    input:
        comparison=RESULTS_DIR + "/validation/literature_comparison.tsv",
        tblouts=expand(
            RESULTS_DIR + "/validation/lit.{peptide_class}.tblout",
            peptide_class=["class-1", "class-2"],
        ),
    output:
        report=RESULTS_DIR + "/reports/validation-report.html",
    params:
        results_dir=RESULTS_DIR,
    log:
        LOGS_DIR + "/report/validation_report.log",
    resources:
        runtime=15,
        mem_mb=4000,
    shell:
        """
        quarto render workflow/reports/validation-report.qmd \
            -P results_dir={params.results_dir} \
            -o {output.report} 2> {log}
        """
