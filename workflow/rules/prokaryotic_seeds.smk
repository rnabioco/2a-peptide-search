"""
Prokaryotic 2A-like discovery - APPROACH 1: Seed-based Discovery

Uses known stalling peptides (SecM, TnaC, MifM, etc.) as seeds to build HMMs
and search prokaryotic proteomes for related sequences.

Strategy:
1. Start with known stalling peptides (PMID 38565864)
2. Align by motif type (RAGP, RAPG, QAPP, etc.)
3. Build seed HMMs from each motif family
4. Search prokaryotic proteomes with seed HMMs
5. Iterative refinement
"""


rule split_known_peptides_by_motif:
    """Split known stalling peptides into separate files by motif type."""
    input:
        fasta="resources/stalling-peptides/known_stalling_peptides.fasta",
    output:
        motif_list=RESULTS_DIR + "/prokaryotic/seeds/motif_list.txt",
        fastas=directory(RESULTS_DIR + "/prokaryotic/seeds/by_motif/"),
    log:
        LOGS_DIR + "/prokaryotic/split_peptides_by_motif.log",
    shell:
        """
        python workflow/scripts/split_peptides_by_motif.py \
            --fasta "{input.fasta}" \
            --output-dir "{output.fastas}" \
            --motif-list "{output.motif_list}" \
            > "{log}" 2>&1
        """


rule align_all_seed_peptides:
    """Create comprehensive alignment from all known stalling peptides."""
    input:
        fasta="resources/stalling-peptides/known_stalling_peptides.fasta",
    output:
        alignment=RESULTS_DIR + "/prokaryotic/seeds/alignments/comprehensive.sto",
    log:
        LOGS_DIR + "/prokaryotic/align_all_seeds.log",
    shell:
        """
        mkdir -p $(dirname "{output.alignment}")

        # Use MAFFT for alignment, convert to Stockholm format
        mafft --auto "{input.fasta}" > "{output.alignment}.afa" 2> "{log}"

        # Convert to Stockholm format
        esl-reformat stockholm "{output.alignment}.afa" > "{output.alignment}" 2>> "{log}"

        rm "{output.alignment}.afa"
        """


rule build_comprehensive_hmm:
    """Build comprehensive 'pan-stalling' HMM from all known peptides."""
    input:
        alignment=RESULTS_DIR + "/prokaryotic/seeds/alignments/comprehensive.sto",
    output:
        hmm=RESULTS_DIR + "/prokaryotic/models/seed/comprehensive.hmm",
    params:
        name="stall-pan",
    log:
        LOGS_DIR + "/prokaryotic/build_comprehensive_hmm.log",
    shell:
        """
        hmmbuild -n {params.name} "{output.hmm}" "{input.alignment}" 2> "{log}"
        """


rule search_with_comprehensive_hmm:
    """Search prokaryotic proteomes with comprehensive pan-stalling HMM."""
    input:
        hmm=RESULTS_DIR + "/prokaryotic/models/seed/comprehensive.hmm",
        db=get_prokaryotic_database_file,
    output:
        hmmsearch=RESULTS_DIR + "/prokaryotic/seed_search_results/{database}/comprehensive_hmm/results.hmmsearch.gz",
        tblout=RESULTS_DIR + "/prokaryotic/seed_search_results/{database}/comprehensive_hmm/results.tblout.gz",
        alignment=RESULTS_DIR + "/prokaryotic/seed_search_results/{database}/comprehensive_hmm/results.sto.gz",
    log:
        LOGS_DIR + "/prokaryotic/search_comprehensive_{database}.log",
    wildcard_constraints:
        database="[^/]+",
    threads: 12
    resources:
        runtime=480,
        mem_mb=16000,
    shell:
        """
        hmmsearch --cpu {threads} \
            --tblout >(gzip > "{output.tblout}") \
            -A >(gzip > "{output.alignment}") \
            --noali \
            "{input.hmm}" "{input.db}" 2> "{log}" | gzip > "{output.hmmsearch}"
        """


rule align_seed_peptides:
    """Create multiple sequence alignment for each motif family."""
    input:
        fasta=RESULTS_DIR + "/prokaryotic/seeds/by_motif/{motif}.fasta",
    output:
        alignment=RESULTS_DIR + "/prokaryotic/seeds/alignments/{motif}.sto",
    log:
        LOGS_DIR + "/prokaryotic/align_seeds_{motif}.log",
    shell:
        """
        mkdir -p $(dirname "{output.alignment}")

        # Use MAFFT for alignment, convert to Stockholm format
        mafft --auto "{input.fasta}" > "{output.alignment}.afa" 2> "{log}"

        # Convert to Stockholm format (hmmer accepts various formats)
        esl-reformat stockholm "{output.alignment}.afa" > "{output.alignment}" 2>> "{log}"

        rm "{output.alignment}.afa"
        """


rule build_seed_hmms:
    """Build HMMs from seed alignments of known stalling peptides."""
    input:
        alignment=RESULTS_DIR + "/prokaryotic/seeds/alignments/{motif}.sto",
    output:
        hmm=RESULTS_DIR + "/prokaryotic/models/seed/{motif}.hmm",
    params:
        name=lambda w: f"stall-{w.motif}",
    log:
        LOGS_DIR + "/prokaryotic/build_seed_hmm_{motif}.log",
    shell:
        """
        hmmbuild -n {params.name} "{output.hmm}" "{input.alignment}" 2> "{log}"
        """


rule search_with_seed_hmms:
    """Search prokaryotic proteomes with seed HMMs from known peptides."""
    input:
        hmm=RESULTS_DIR + "/prokaryotic/models/seed/{motif}.hmm",
        db=get_prokaryotic_database_file,
    output:
        hmmsearch=RESULTS_DIR + "/prokaryotic/seed_search_results/{database}/{motif}.hmmsearch.gz",
        tblout=RESULTS_DIR + "/prokaryotic/seed_search_results/{database}/{motif}.tblout.gz",
        alignment=RESULTS_DIR + "/prokaryotic/seed_search_results/{database}/{motif}.sto.gz",
    log:
        LOGS_DIR + "/prokaryotic/seed_search_{database}_{motif}.log",
    wildcard_constraints:
        database="[^/]+",
        motif="[^/]+",
    threads: 12
    resources:
        runtime=480,
        mem_mb=16000,
    shell:
        """
        hmmsearch --cpu {threads} \
            --tblout >(gzip > "{output.tblout}") \
            -A >(gzip > "{output.alignment}") \
            --noali \
            "{input.hmm}" "{input.db}" 2> "{log}" | gzip > "{output.hmmsearch}"
        """


rule merge_comprehensive_searches:
    """Merge comprehensive search results from all databases."""
    input:
        alignments=expand(
            RESULTS_DIR + "/prokaryotic/seed_search_results/{database}/comprehensive_hmm/results.sto.gz",
            database=config["prokaryotic_databases_to_search"],
        ),
    output:
        merged=RESULTS_DIR + "/prokaryotic/seed_search_results/comprehensive_merged.sto.gz",
    log:
        LOGS_DIR + "/prokaryotic/merge_comprehensive_searches.log",
    shell:
        """
        # Decompress all alignments
        for aln in {input.alignments}; do
            zcat "$aln"
        done | \
        # Merge with esl-alimerge and compress
        esl-alimerge --list - 2> "{log}" | gzip > "{output.merged}"
        """
