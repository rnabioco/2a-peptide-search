"""
Rules for building and refining HMM models.
"""

rule build_seed_models:
    """Build initial HMM models from seed alignments."""
    input:
        alignment="resources/seed-alignments/2A-{peptide_class}.sto.gz"
    output:
        hmm="results/models/seed/2A-{peptide_class}.hmm"
    params:
        name=lambda w: f"2A-{w.peptide_class}"
    log:
        "logs/hmmbuild/seed_{peptide_class}.log"
    conda:
        "../envs/hmmer.yaml"
    shell:
        """
        hmmbuild -n {params.name} {output.hmm} {input.alignment} 2> {log}
        """


rule build_refined_model:
    """Build refined HMM from merged alignment."""
    input:
        alignment="results/alignments/{iteration}/2A-{peptide_class}.merged.sto"
    output:
        hmm="results/models/{iteration}_refined/2A-{peptide_class}.hmm"
    params:
        name=lambda w: f"2A-{w.peptide_class}"
    log:
        "logs/hmmbuild/{iteration}_refined_{peptide_class}.log"
    conda:
        "../envs/hmmer.yaml"
    shell:
        """
        hmmbuild -n {params.name} {output.hmm} {input.alignment} 2> {log}
        """


rule create_checkpoint:
    """Create checkpoint for manual curation."""
    input:
        alignments=expand(
            "results/alignments/{{iteration}}/2A-{peptide_class}.merged.sto",
            peptide_class=["class-1", "class-2"]
        )
    output:
        checkpoint=touch("results/checkpoints/{iteration}.curated")
    message:
        "Manual curation required for iteration {wildcards.iteration}. "
        "Review alignments in results/alignments/{wildcards.iteration}/ and "
        "touch {output.checkpoint} when done."


rule build_final_models:
    """Build final production HMMs after manual curation."""
    input:
        alignment="results/alignments/final/2A-{peptide_class}.curated.sto",
        checkpoint="results/checkpoints/iter2.curated"
    output:
        hmm="results/models/final/2A-{peptide_class}.hmm"
    params:
        name=lambda w: f"2A-{w.peptide_class}"
    log:
        "logs/hmmbuild/final_{peptide_class}.log"
    conda:
        "../envs/hmmer.yaml"
    shell:
        """
        hmmbuild -n {params.name} {output.hmm} {input.alignment} 2> {log}
        """
