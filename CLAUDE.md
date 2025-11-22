# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a bioinformatics research project for identifying 2A peptides in protein databases using profile hidden Markov models (HMMs). 2A peptides are short (~15 residue) cis-acting oligopeptides that cause ribosomal skipping at Gly-Pro dipeptides. The project identifies two distinct classes of 2A peptides (class 1 and class 2) with conserved sequence features.

The repository is organized as a Snakemake pipeline for reproducible analysis from data download through iterative HMM refinement to final reporting.

## Key Software Dependencies

- **Snakemake** - Workflow management system
- **Pixi** - Package and environment management (all dependencies defined in `pixi.toml`)
- **HMMER suite** (`hmmsearch`, `hmmbuild`, `hmmalign`) - Core tool for HMM-based sequence searches
- **Easel tools** (`esl-alimerge`) - For merging multiple sequence alignments
- **Python 3** with Biopython (`Bio.AlignIO`, `Bio.SeqIO`) - For sequence processing
- **R/Quarto** - For analysis and reporting

**Note:** This project uses **Pixi only** for dependency management, not conda/mamba.

## Repository Structure

```
workflow/
├── Snakefile                    # Main pipeline orchestration
├── config.yaml                  # Configuration (databases, thresholds)
├── README.md                    # Detailed workflow documentation
├── scripts/                     # Python scripts for data processing
└── rules/                       # Modular Snakemake rules
    ├── download.smk            # Database download
    ├── search.smk              # HMM searches
    ├── refine.smk              # Model building
    ├── logos.smk               # Sequence logo generation
    └── report.smk              # Report generation

cluster/slurm/                   # SLURM cluster configuration
├── config.yaml                  # SLURM resource specifications
└── README.md                    # SLURM setup and usage guide

resources/seed-alignments/       # Curated starting alignments (version controlled)
results/                         # All pipeline outputs (gitignored)
data/                           # Downloaded databases (gitignored)
legacy/                         # Archived historical results and scripts

pixi.toml                        # Pixi environment and dependency specification
submit-slurm.sh                  # SLURM orchestrator submission script
submit-test.sh                   # Quick test submission script
```

## Common Commands

### Running the full pipeline

**Local/Workstation:**
```bash
# Dry run to see what will be executed
pixi run dry-run

# Run full pipeline with 12 cores
pixi run run

# Generate workflow visualization
pixi run dag

# Or run snakemake directly within pixi environment
pixi run snakemake --cores 12
```

**SLURM Cluster (Alpine):**
```bash
# Configure account in cluster/slurm/config.yaml first
# Then submit orchestrator job:
sbatch submit-slurm.sh

# Or test first:
sbatch submit-test.sh

# Using profile directly:
pixi run snakemake --profile cluster/slurm
```

### Testing with UniProt only

```bash
# Local
pixi run test

# SLURM
sbatch submit-test.sh
# or
pixi run snakemake --profile cluster/slurm test
```

### Specific pipeline stages

```bash
# Download all configured databases
pixi run download

# Build seed models from curated alignments
pixi run build-seeds

# Run iteration 1 (automated)
pixi run snakemake iter1 --cores 12

# Run iteration 2 (requires manual curation checkpoint)
pixi run snakemake iter2 --cores 12
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
- All dependencies are managed via Pixi (defined in `pixi.toml`)
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
- **Plotting**: Use `plotnine` (NOT matplotlib or seaborn) for all data visualization
- **Image manipulation**: Use `pillow` (PIL) for creating images
- **Code formatting**: Use `ruff` for Python formatting (available in dev environment)
- Test rules individually: `pixi run snakemake <target> --cores 1`
- Legacy code is archived in `legacy/` for reference
- The models can identify partial cross-matches between classes due to the conserved C-terminal PGP motif

## SLURM Cluster Usage

The pipeline is configured for the CU Boulder Alpine cluster with comprehensive SLURM support:

### Setup
1. Edit `cluster/slurm/config.yaml` and set your account:
   ```yaml
   slurm_account: amc-general  # Change to your allocation
   ```

2. Submit orchestrator job:
   ```bash
   sbatch submit-slurm.sh [target]
   ```

### Resource Specifications
- Small jobs (downloads, model building): 1-2 CPUs, 2-4GB, 10-60 min
- Medium jobs (alignment processing): 2-4 CPUs, 8-16GB, 30-60 min
- Large jobs (HMM searches): 12 CPUs, 16GB, up to 24 hours

### Monitoring
```bash
# Check your jobs
squeue -u $USER

# View orchestrator log
tail -f logs/orchestrator_*.out

# Check individual rule logs
ls logs/

# Cancel jobs
scancel -u $USER
```

See `cluster/slurm/README.md` for complete SLURM documentation.

## Troubleshooting

- **Out of memory**: Increase `mem_mb` in `cluster/slurm/config.yaml` or reduce `--cores` for local runs
- **Download fails**: Check URLs in `workflow/config.yaml` are current
- **Dependency issues**: Run `pixi install` to update environment, or `pixi update` to upgrade packages
- **Missing checkpoint**: Manually create the checkpoint file to continue pipeline
- **SLURM job fails**: Check `.snakemake/slurm_logs/` and increase resources in cluster config
- **Timeout on cluster**: Increase `runtime` in `cluster/slurm/config.yaml` for specific rules
