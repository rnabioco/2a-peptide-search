"""
Rules for discovering prokaryotic 2A-like ribosomal stalling peptides.

This module orchestrates two complementary discovery approaches:

APPROACH 1: Seed-based discovery (Targeted)
  - Uses known stalling peptides (SecM, TnaC, MifM) as seeds
  - Defined in: prokaryotic_seeds.smk

APPROACH 2: Unbiased GP discovery (Comprehensive)
  - Extracts all GP-containing sequences from prokaryotic proteomes
  - Focuses on inter-domain GP motifs
  - Defined in: prokaryotic_gp.smk

Analysis and Reporting:
  - Validation against known peptides
  - Comparative analysis of approaches
  - Defined in: prokaryotic_analysis.smk

Known prokaryotic stalling peptides:
- SecM: FXXXXWIXXXXGIRAGP (stalls at RAGP)
- TnaC: WDPXXXX...GP-like
- MifM: Contains stalling sequence
"""

# ============================================================================
# Helper Functions
# ============================================================================


def get_prokaryotic_database_file(wildcards):
    """Map database name to its file path, handling different download mechanisms.

    For databases that may contain invalid IUPAC characters (like IMG/VR),
    automatically use the sanitized version to prevent HMMER parse errors.
    """
    db_config = config["prokaryotic_databases"][wildcards.database]

    # IMG/VR: Always use sanitized version (may contain invalid chars like '-')
    if wildcards.database == "imgvr":
        return DATA_DIR + "/prokaryotic/imgvr.sanitized.fasta.gz"

    # Check if database has a local_path configured
    if "local_path" in db_config and db_config["local_path"]:
        return db_config["local_path"]

    # Databases with direct URL downloads (go to downloads/ subdirectory)
    if wildcards.database in ["bacteria", "archaea", "uniprot_viruses"]:
        return DATA_DIR + f"/prokaryotic/downloads/{wildcards.database}.fasta.gz"

    # Databases from merged splits (go to prokaryotic/ directly)
    if wildcards.database == "ncbi_viral_refseq":
        return DATA_DIR + "/prokaryotic/ncbi_viral_refseq.fasta.gz"

    # Databases from ORF prediction (go to prokaryotic/ directly)
    if wildcards.database in ["inphared_proteins", "millardlab_proteins"]:
        return DATA_DIR + f"/prokaryotic/{wildcards.database}.fasta.gz"

    # Gut Phage Database (CyVerse/iVirus)
    if wildcards.database == "gpd":
        return DATA_DIR + "/prokaryotic/gpd.fasta.gz"

    # GOV (Global Ocean Virome) predicted proteins
    if wildcards.database == "tara_proteins":
        return RESULTS_DIR + "/phage/tara/predicted_orfs.faa.gz"

    if wildcards.database == "malaspina_proteins":
        return RESULTS_DIR + "/phage/malaspina/predicted_orfs.faa.gz"

    # MGnify uses main pipeline download
    if wildcards.database == "mgnify":
        return f"{DATA_DIR}/mgnify/mgnify_proteins.fasta.gz"

    # Fallback to downloads subdirectory
    return DATA_DIR + f"/prokaryotic/downloads/{wildcards.database}.fasta.gz"


def aggregate_hmmscan_chunks(wildcards):
    """Dynamically get all chunk domtblout files after checkpoint completes."""
    checkpoint_output = checkpoints.split_sequences_for_hmmscan.get(**wildcards).output[0]
    # Find all chunk files created by the checkpoint
    chunk_ids = glob_wildcards(
        RESULTS_DIR + f"/prokaryotic/gp_analysis/{wildcards.database}/chunks/chunk_{{n}}.fasta.gz"
    ).n
    return expand(
        RESULTS_DIR + "/prokaryotic/gp_analysis/{database}/chunks/chunk_{chunk_id}.domtblout",
        database=wildcards.database,
        chunk_id=chunk_ids,
    )


# ============================================================================
# Include Modular Rule Files
# ============================================================================

# Approach 1: Seed-based discovery with known stalling peptides
include: "prokaryotic_seeds.smk"

# Approach 2: Unbiased GP motif extraction and clustering
include: "prokaryotic_gp.smk"

# Approach 3: C-terminal arrest motif discovery (RAGP, QAPP, etc.)
include: "prokaryotic_arrest.smk"

# Validation, comparative analysis, and reporting
include: "prokaryotic_analysis.smk"
