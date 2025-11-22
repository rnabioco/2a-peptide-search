"""
Rules for building and refining HMM models.
"""


rule build_seed_models:
    """Build initial HMM models from seed alignments."""
    input:
        alignment="resources/seed-alignments/2A-{peptide_class}.sto.gz",
    output:
        hmm=RESULTS_DIR + "/models/seed/2A-{peptide_class}.hmm",
    params:
        name=lambda w: f"2A-{w.peptide_class}",
    log:
        LOGS_DIR + "/hmmbuild/seed_{peptide_class}.log",
    shell:
        """
        hmmbuild -n {params.name} {output.hmm} {input.alignment} 2> {log}
        """


rule build_refined_model:
    """Build refined HMM from merged alignment."""
    input:
        alignment=SCRATCH_DIR + "/alignments/{iteration}/2A-{peptide_class}.merged.sto",
    output:
        hmm=RESULTS_DIR + "/models/{iteration}_refined/2A-{peptide_class}.hmm",
    params:
        name=lambda w: f"2A-{w.peptide_class}",
    log:
        LOGS_DIR + "/hmmbuild/{iteration}_refined_{peptide_class}.log",
    shell:
        """
        hmmbuild -n {params.name} {output.hmm} {input.alignment} 2> {log}
        """


rule create_checkpoint:
    """Create checkpoint for manual curation."""
    input:
        alignments=expand(
            SCRATCH_DIR + "/alignments/{{iteration}}/2A-{peptide_class}.merged.sto",
            peptide_class=["class-1", "class-2"],
        ),
    output:
        checkpoint=touch(RESULTS_DIR + "/checkpoints/{iteration}.curated"),
    message:
        "Manual curation required for iteration {wildcards.iteration}. "
        "Review alignments in results/alignments/{wildcards.iteration}/ and "
        "touch {output.checkpoint} when done."


rule build_final_models:
    """Build final production HMMs after manual curation."""
    input:
        alignment=RESULTS_DIR + "/alignments/final/2A-{peptide_class}.curated.sto",
        checkpoint=RESULTS_DIR + "/checkpoints/iter2.curated",
    output:
        hmm=RESULTS_DIR + "/models/final/2A-{peptide_class}.hmm",
    params:
        name=lambda w: f"2A-{w.peptide_class}",
    log:
        LOGS_DIR + "/hmmbuild/final_{peptide_class}.log",
    shell:
        """
        hmmbuild -n {params.name} {output.hmm} {input.alignment} 2> {log}
        """
