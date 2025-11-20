# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a bioinformatics research project for identifying 2A peptides in protein databases using profile hidden Markov models (HMMs). 2A peptides are short (~15 residue) cis-acting oligopeptides that cause ribosomal skipping at Gly-Pro dipeptides. The project identifies two distinct classes of 2A peptides (class 1 and class 2) with conserved sequence features.

The repository is organized as a Snakemake pipeline for reproducible analysis from data download through iterative HMM refinement to final reporting.

## Key Software Dependencies

- **Snakemake** - Workflow management system
- **Conda/Mamba** - Environment management (environments defined in `workflow/envs/`)
- **HMMER suite** (`hmmsearch`, `hmmbuild`, `hmmalign`) - Core tool for HMM-based sequence searches
- **Easel tools** (`esl-alimerge`) - For merging multiple sequence alignments
- **Python 3** with Biopython (`Bio.AlignIO`, `Bio.SeqIO`) - For sequence processing
- **R/Quarto** - For analysis and reporting

## Repository Structure

```
workflow/
├── Snakefile                    # Main pipeline orchestration
├── config.yaml                  # Configuration (databases, thresholds)
├── envs/                        # Conda environment specifications
├── scripts/                     # Python scripts for data processing
└── rules/                       # Modular Snakemake rules
    ├── download.smk            # Database download
    ├── search.smk              # HMM searches
    ├── refine.smk              # Model building
    └── report.smk              # Report generation

resources/seed-alignments/       # Curated starting alignments (version controlled)
results/                         # All pipeline outputs (gitignored)
data/                           # Downloaded databases (gitignored)
legacy/                         # Archived historical results and scripts
```

## Common Commands

### Running the full pipeline

```bash
# Dry run to see what will be executed
snakemake -n

# Run with conda environments and 12 cores
snakemake --use-conda --cores 12

# Generate workflow visualization
snakemake --dag | dot -Tpng > workflow.png
```

### Testing with UniProt only

```bash
# Quick test with just UniProt database
snakemake test --use-conda --cores 12
```

### Specific pipeline stages

```bash
# Download all configured databases
snakemake download_all --use-conda --cores 4

# Build seed models from curated alignments
snakemake build_seeds --use-conda

# Run iteration 1 (automated)
snakemake iter1 --use-conda --cores 12

# Run iteration 2 (requires manual curation checkpoint)
snakemake iter2 --use-conda --cores 12
```

### Configuration

Edit `workflow/config.yaml` to:
- Select which databases to search (comment out large databases for testing)
- Adjust E-value and identity thresholds
- Modify resource allocation (threads, retries)

### Manual Curation Checkpoints

The pipeline includes manual curation checkpoints:

1. Pipeline pauses after iteration 2
2. Review alignments in `results/alignments/iter2/`
3. Manually curate and save to `results/alignments/final/*.curated.sto`
4. Create checkpoint: `touch results/checkpoints/iter2.curated`
5. Continue pipeline to build final models

## Pipeline Workflow

1. **Download** → Fetch protein databases (UniProt, Reference Proteomes, UniParc, MGnify)
2. **Build Seed Models** → Create initial HMMs from curated seed alignments
3. **Search Iteration 1** → Search all databases with seed HMMs
4. **Refine Iteration 1** → Build refined HMMs from high-confidence hits
5. **Search Iteration 2** → Search with refined HMMs
6. **Manual Checkpoint** → User reviews and curates alignments
7. **Build Final Models** → Create production HMMs from curated alignments
8. **Generate Report** → Create Quarto document with analysis summary

## Key Architectural Notes

### Two-Class System

The project distinguishes two classes of 2A peptides:
- **Class 1**: N-terminal leucine stretch, central GDVE motif, C-terminal NPGP
- **Class 2**: N-terminal invariant tryptophan, central EEGIE motif, C-terminal PNPGP/PHPGP

Always search with both class 1 and class 2 models, as they identify distinct but related sequences.

### File Naming Conventions

- HMM models: `2A-class-{1,2}.hmm`
- Alignments: `2A-class-{1,2}.{filtered|merged|curated}.sto`
- Search results: `2A-class-{1,2}.{hmmsearch|tblout|sto}.gz`
- Always use gzip compression for large files

### Stockholm Alignment Format

Alignments use Stockholm format (`.sto` files), which includes both the alignment and metadata. These are the input for `hmmbuild` and output from `hmmsearch -A`. The coordinate information in sequence IDs follows the pattern: `accession/start-end`.

### Snakemake Best Practices

- Rules are modular and in separate files under `workflow/rules/`
- Each rule specifies its conda environment
- Wildcards enable parallel execution across databases and peptide classes
- Checkpoints allow for manual intervention in automated workflows

## Data Sources

Major protein databases searched:
- **UniProt** - Curated reference proteins (~500k sequences)
- **Reference Proteomes** - Subset of UniProt used for Pfam
- **UniParc** - Non-redundant uncurated proteins (millions of sequences)
- **MGnify** - Environmental/metagenomic sequences (very large)

## Development Notes

- Scripts in `workflow/scripts/` are called by Snakemake rules
- Each script should be standalone and use click for CLI
- Test rules individually: `snakemake <target> --use-conda --cores 1`
- Legacy code is archived in `legacy/` for reference
- The models can identify partial cross-matches between classes due to the conserved C-terminal PGP motif

## Troubleshooting

- **Out of memory**: Reduce `--cores` or disable large databases (uniparc, mgnify) in config
- **Download fails**: Check URLs in `workflow/config.yaml` are current
- **Conda issues**: Use `--conda-frontend mamba` for faster environment resolution
- **Missing checkpoint**: Manually create the checkpoint file to continue pipeline
