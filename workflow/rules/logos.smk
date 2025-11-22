"""
Rules for generating sequence logos using Skylign API.
"""


rule generate_logo_seed:
    """Generate sequence logo from seed alignment."""
    input:
        alignment="resources/seed-alignments/2A-{peptide_class}.sto.gz",
    output:
        logo=RESULTS_DIR + "/logos/seed/2A-{peptide_class}.png",
    params:
        api_url=config.get("skylign", {}).get("api_url", "http://skylign.org"),
        processing=config.get("skylign", {}).get("processing", "hmm"),
        retry_attempts=config.get("skylign", {}).get("retry_attempts", 3),
        retry_delay=config.get("skylign", {}).get("retry_delay", 5),
    log:
        LOGS_DIR + "/skylign/seed_{peptide_class}.log",
    shell:
        """
        python workflow/scripts/generate_skylign_logos.py \
            {input.alignment} \
            {output.logo} \
            --api-url {params.api_url} \
            --processing {params.processing} \
            --retry-attempts {params.retry_attempts} \
            --retry-delay {params.retry_delay} \
            > {log} 2>&1
        """


rule generate_logo_iter1:
    """Generate sequence logo from iteration 1 merged alignment."""
    input:
        alignment=SCRATCH_DIR + "/alignments/iter1/2A-{peptide_class}.merged.sto",
    output:
        logo=RESULTS_DIR + "/logos/iter1/2A-{peptide_class}.png",
    params:
        api_url=config.get("skylign", {}).get("api_url", "http://skylign.org"),
        processing=config.get("skylign", {}).get("processing", "hmm"),
        retry_attempts=config.get("skylign", {}).get("retry_attempts", 3),
        retry_delay=config.get("skylign", {}).get("retry_delay", 5),
    log:
        LOGS_DIR + "/skylign/iter1_{peptide_class}.log",
    shell:
        """
        python workflow/scripts/generate_skylign_logos.py \
            {input.alignment} \
            {output.logo} \
            --api-url {params.api_url} \
            --processing {params.processing} \
            --retry-attempts {params.retry_attempts} \
            --retry-delay {params.retry_delay} \
            > {log} 2>&1
        """


rule generate_logo_iter2:
    """Generate sequence logo from iteration 2 merged alignment."""
    input:
        alignment=SCRATCH_DIR + "/alignments/iter2/2A-{peptide_class}.merged.sto",
    output:
        logo=RESULTS_DIR + "/logos/iter2/2A-{peptide_class}.png",
    params:
        api_url=config.get("skylign", {}).get("api_url", "http://skylign.org"),
        processing=config.get("skylign", {}).get("processing", "hmm"),
        retry_attempts=config.get("skylign", {}).get("retry_attempts", 3),
        retry_delay=config.get("skylign", {}).get("retry_delay", 5),
    log:
        LOGS_DIR + "/skylign/iter2_{peptide_class}.log",
    shell:
        """
        python workflow/scripts/generate_skylign_logos.py \
            {input.alignment} \
            {output.logo} \
            --api-url {params.api_url} \
            --processing {params.processing} \
            --retry-attempts {params.retry_attempts} \
            --retry-delay {params.retry_delay} \
            > {log} 2>&1
        """


rule generate_logo_final:
    """Generate sequence logo from final HMM model."""
    input:
        hmm=RESULTS_DIR + "/models/final/2A-{peptide_class}.hmm",
    output:
        logo=RESULTS_DIR + "/logos/final/2A-{peptide_class}.png",
    params:
        api_url=config.get("skylign", {}).get("api_url", "http://skylign.org"),
        processing=config.get("skylign", {}).get("processing", "hmm"),
        retry_attempts=config.get("skylign", {}).get("retry_attempts", 3),
        retry_delay=config.get("skylign", {}).get("retry_delay", 5),
    log:
        LOGS_DIR + "/skylign/final_{peptide_class}.log",
    shell:
        """
        python workflow/scripts/generate_skylign_logos.py \
            {input.hmm} \
            {output.logo} \
            --api-url {params.api_url} \
            --processing {params.processing} \
            --retry-attempts {params.retry_attempts} \
            --retry-delay {params.retry_delay} \
            > {log} 2>&1
        """


# Convenience rules for generating all logos at each stage
rule logos_seed:
    """Generate all seed alignment logos."""
    input:
        expand(
            RESULTS_DIR + "/logos/seed/2A-{peptide_class}.png",
            peptide_class=["class-1", "class-2"],
        ),


rule logos_iter1:
    """Generate all iteration 1 logos."""
    input:
        expand(
            RESULTS_DIR + "/logos/iter1/2A-{peptide_class}.png",
            peptide_class=["class-1", "class-2"],
        ),


rule logos_iter2:
    """Generate all iteration 2 logos."""
    input:
        expand(
            RESULTS_DIR + "/logos/iter2/2A-{peptide_class}.png",
            peptide_class=["class-1", "class-2"],
        ),


rule logos_final:
    """Generate all final model logos."""
    input:
        expand(
            RESULTS_DIR + "/logos/final/2A-{peptide_class}.png",
            peptide_class=["class-1", "class-2"],
        ),


rule logos_all:
    """Generate all logos for main 2A peptide pipeline."""
    input:
        rules.logos_seed.input,
        rules.logos_iter1.input,
        rules.logos_iter2.input,
        rules.logos_final.input,
