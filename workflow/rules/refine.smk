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
    """Create checkpoint for manual curation (optional workflow)."""
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


rule auto_curate_alignment:
    """Automatically curate alignment for final model building.

    Applies quality filters to produce clean alignments without manual intervention:
    - E-value filtering
    - Minimum sequence length
    - Maximum gap percentage
    - C-terminal motif validation (PGP conservation)
    - Bit score outlier removal
    """
    input:
        alignment=SCRATCH_DIR + "/alignments/iter2/2A-{peptide_class}.merged.sto",
    output:
        curated=RESULTS_DIR + "/alignments/final/2A-{peptide_class}.auto-curated.sto",
    params:
        evalue=config["thresholds"]["evalue"],
        min_length=15,  # 2A peptides are ~15-20 residues
        max_gap_pct=0.5,
    log:
        LOGS_DIR + "/auto_curate/{peptide_class}.log",
    resources:
        runtime=30,
        mem_mb=4000,
    shell:
        """
        python workflow/scripts/auto_curate_alignment.py \
            {input.alignment} \
            {output.curated} \
            --evalue {params.evalue} \
            --min-length {params.min_length} \
            --max-gap-pct {params.max_gap_pct} \
            --require-motif \
            --bitscore-percentile 10 \
            2> {log}
        """


rule build_final_models:
    """Build final production HMMs from auto-curated alignments."""
    input:
        alignment=RESULTS_DIR + "/alignments/final/2A-{peptide_class}.auto-curated.sto",
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


rule build_manually_curated_models:
    """Build final HMMs from manually curated alignments (optional).

    Use this rule when you want maximum control over the final models.
    Requires manual curation of alignments in results/alignments/final/.
    """
    input:
        alignment=RESULTS_DIR + "/alignments/final/2A-{peptide_class}.curated.sto",
        checkpoint=RESULTS_DIR + "/checkpoints/iter2.curated",
    output:
        hmm=RESULTS_DIR + "/models/manual/2A-{peptide_class}.hmm",
    params:
        name=lambda w: f"2A-{w.peptide_class}",
    log:
        LOGS_DIR + "/hmmbuild/manual_{peptide_class}.log",
    shell:
        """
        hmmbuild -n {params.name} {output.hmm} {input.alignment} 2> {log}
        """
