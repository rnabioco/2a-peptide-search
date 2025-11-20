# 2A Peptide Search Pipeline

Snakemake workflow for identifying 2A peptides in protein databases using profile HMMs and iterative refinement.

## Quick Start

```bash
# Install Snakemake (if not already installed)
conda install -n base -c conda-forge -c bioconda snakemake

# Dry run to see planned steps
snakemake -n

# Run full pipeline with 12 cores
snakemake --use-conda --cores 12

# Or test with UniProt only
snakemake test --use-conda --cores 12
```

## Pipeline Overview

This workflow performs:

1. **Download** - Fetch protein databases (UniProt, Reference Proteomes, optionally UniParc/MGnify)
2. **Seed Search** - Search databases with curated seed HMMs
3. **Iteration 1** - Build refined HMMs from high-confidence hits, re-search
4. **Iteration 2** - Further refinement with automated quality filtering
5. **Manual Curation** - Review and curate alignments (checkpoint)
6. **Final Models** - Build production HMMs from curated data
7. **Report** - Generate comprehensive analysis report

## Configuration

Edit `workflow/config.yaml` to customize:

### Select Databases

```yaml
databases_to_search:
  - uniprot              # ~500k sequences, quick
  - reference_proteomes  # Medium size
  # - uniparc           # Millions of sequences, slow
  # - mgnify            # Very large, environmental
```

**Recommendation**: Start with `uniprot` and `reference_proteomes` for testing. Enable `uniparc` and `mgnify` for comprehensive searches (requires significant compute time).

### Adjust Thresholds

```yaml
thresholds:
  evalue: 1e-5           # E-value cutoff for hits
  identity: 0.95         # Redundancy filtering threshold
  min_coverage: 0.8      # Minimum alignment coverage
```

### Resource Allocation

```yaml
resources:
  hmmsearch_threads: 12  # Threads per search
  download_retries: 3    # Network retry attempts
```

## Pipeline Targets

### Complete Pipeline
```bash
snakemake --use-conda --cores 12
```
Runs entire workflow through report generation (requires manual curation).

### Specific Stages

```bash
# Download all configured databases
snakemake download_all --use-conda --cores 4

# Build seed models only
snakemake build_seeds --use-conda

# Run iteration 1 (automated)
snakemake iter1 --use-conda --cores 12

# Run iteration 2 (automated, creates checkpoint)
snakemake iter2 --use-conda --cores 12

# Test with UniProt only
snakemake test --use-conda --cores 12
```

## Manual Curation Checkpoint

After iteration 2, the pipeline pauses for manual curation:

1. **Review alignments** in `results/alignments/iter2/`
   ```bash
   # View with Jalview or similar alignment viewer
   jalview results/alignments/iter2/2A-class-1.merged.sto
   ```

2. **Curate alignments** (remove low-quality sequences, adjust boundaries)

3. **Save curated alignments** to `results/alignments/final/`
   ```bash
   cp results/alignments/iter2/2A-class-1.merged.sto \
      results/alignments/final/2A-class-1.curated.sto
   ```

4. **Create checkpoint** to resume pipeline
   ```bash
   touch results/checkpoints/iter2.curated
   ```

5. **Continue pipeline**
   ```bash
   snakemake --use-conda --cores 12
   ```

## Output Structure

```
results/
├── models/
│   ├── seed/                      # Initial HMMs from seed alignments
│   ├── iter1_refined/             # Refined HMMs after iteration 1
│   ├── iter2_refined/             # Refined HMMs after iteration 2
│   └── final/                     # Production HMMs (curated)
├── searches/
│   ├── uniprot/
│   │   ├── seed/                  # Seed model searches
│   │   ├── iter1_refined/         # Iteration 1 searches
│   │   └── final/                 # Final searches
│   └── ...                        # Other databases
├── alignments/
│   ├── iter1/                     # Merged alignments from iteration 1
│   ├── iter2/                     # Merged alignments from iteration 2
│   └── final/                     # Manually curated alignments
├── checkpoints/
│   └── iter2.curated              # Manual curation checkpoint
└── reports/
    └── 2A-peptide-analysis.html   # Final comprehensive report
```

## Key Files

- **Snakefile** - Main workflow logic
- **config.yaml** - Configuration parameters
- **rules/*.smk** - Modular workflow rules
- **scripts/*.py** - Data processing scripts
- **envs/*.yaml** - Conda environment specifications

## Workflow Visualization

```bash
# Generate workflow DAG
snakemake --dag | dot -Tpng > workflow.png

# Generate rule graph
snakemake --rulegraph | dot -Tpng > rulegraph.png
```

## Troubleshooting

### Out of Memory
- Reduce `--cores` to limit parallel jobs
- Disable large databases (uniparc, mgnify) in config.yaml
- Use `--resources mem_mb=<limit>` to constrain memory

### Download Failures
- Check internet connection
- Verify URLs in config.yaml are current
- Increase `download_retries` in config.yaml

### Conda Environment Issues
```bash
# Use mamba for faster resolution
snakemake --use-conda --conda-frontend mamba --cores 12

# Clean conda cache if corrupted
conda clean --all
```

### Resume After Failure
Snakemake automatically resumes from where it stopped:
```bash
snakemake --use-conda --cores 12
```

## Advanced Usage

### Run on HPC Cluster

For LSF clusters:
```bash
snakemake --cluster "bsub -n {threads} -o logs/{rule}.out -e logs/{rule}.err" \
          --jobs 100 \
          --use-conda
```

For SLURM clusters:
```bash
snakemake --cluster "sbatch --cpus-per-task={threads} --output=logs/{rule}.out" \
          --jobs 100 \
          --use-conda
```

### Profile Creation

Create a Snakemake profile for your cluster in `~/.config/snakemake/`:
```yaml
# Example: ~/.config/snakemake/lsf/config.yaml
cluster: "bsub -n {threads} -o logs/{rule}.out"
jobs: 100
use-conda: true
```

Then run:
```bash
snakemake --profile lsf
```

### Custom Rule Execution

Run a specific rule with specific parameters:
```bash
snakemake results/models/seed/2A-class-1.hmm --use-conda
```

## Dependencies

All dependencies are managed via conda environments (see `workflow/envs/`):

- **hmmer** - HMMER suite (hmmsearch, hmmbuild, easel tools)
- **python-bio** - Python with Biopython, pyfaidx, click, pandas
- **r-quarto** - R with tidyverse, cowplot, glue, and Quarto

## Citation

If you use these 2A peptide models, please cite:

- Luke GA, et al. (2008) Occurrence, function and evolutionary origins of '2A-like' sequences in virus genomes. J Gen Virol.
- de Lima JGS, Lanza DCF (2021) 2A and 2A-like Sequences: Distribution in Different Virus Species and Applications in Biotechnology. Viruses.

## Support

For issues or questions:
- Check [GitHub Issues](https://github.com/rnabioco/2a-peptide-search/issues)
- Review `CLAUDE.md` for development documentation
- Consult Snakemake documentation: https://snakemake.readthedocs.io/
