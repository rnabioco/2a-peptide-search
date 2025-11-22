# Using Pixi for Dependency Management

Pixi is now configured to manage the core workflow dependencies (Snakemake, Python, etc.). Individual tool environments (HMMER, Biopython, R, etc.) are still managed by Snakemake's `--use-conda` flag.

## Quick Start

```bash
# Activate the pixi environment (runs commands with pixi)
pixi shell

# Or run commands directly with pixi
pixi run snakemake --version
```

## Available Tasks

Pixi provides convenient shortcuts for common commands:

```bash
# Run full pipeline
pixi run run

# Dry run (preview what will execute)
pixi run dry-run

# Test with UniProt only
pixi run test

# Download databases
pixi run download

# Build seed models
pixi run build-seeds

# Generate workflow visualization
pixi run dag
```

## Environments

Two environments are available:

- **default**: Core Snakemake environment (activated by default)
- **dev**: Includes development tools (ruff, ipython)

```bash
# Activate dev environment
pixi shell -e dev
```

## How It Works

1. **Pixi manages**: Snakemake, Python, Graphviz (for DAG visualization)
2. **Snakemake manages**: Tool-specific environments via `--use-conda` flag
   - HMMER environment (`workflow/envs/hmmer.yaml`)
   - Python/Biopython environment (`workflow/envs/python.yaml`)
   - ORF prediction environment (`workflow/envs/orfs.yaml`)
   - R/Quarto environment (`workflow/envs/r-quarto.yaml`)

This hybrid approach leverages pixi's speed for the core workflow while maintaining compatibility with Snakemake's existing conda integration.

## Benefits

- **Fast**: Pixi is significantly faster than conda for environment resolution
- **Reproducible**: `pixi.lock` ensures exact versions across systems
- **Convenient**: Pre-defined tasks simplify common operations
- **Compatible**: Works seamlessly with existing Snakemake conda workflows

## Traditional Usage Still Works

You can still use the traditional approach:

```bash
# Activate pixi shell first
pixi shell

# Then run snakemake commands as usual
snakemake --use-conda --cores 12
snakemake test --use-conda --cores 4
```

## Updating Dependencies

```bash
# Update all dependencies to latest compatible versions
pixi update

# Add a new dependency
pixi add <package-name>

# Add to dev environment only
pixi add --feature dev <package-name>
```
