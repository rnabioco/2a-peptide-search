# LSF Executor Profile for 2A Peptide Search Pipeline

This directory contains the Snakemake profile configuration for running the 2A peptide search pipeline on LSF clusters (e.g., Fiji) using the LSF executor.

## Setup

### 1. Install LSF Executor Plugin

```bash
# Using pixi (recommended)
pixi add snakemake-executor-plugin-lsf

# Or using pip
pip install snakemake-executor-plugin-lsf
```

### 2. Configure Your LSF Project

Edit `cluster/lsf/config.yaml` and set your LSF project:

```yaml
default-resources:
  - 'lsf_project=your-project-here'  # e.g., "rnabioco"
  - 'lsf_queue=rna'                  # Change queue if needed
```

**Important:** You must set the project before running the pipeline, or jobs will fail to submit.

## Usage

### Quick Test Run

```bash
# Dry run to see execution plan
pixi run snakemake --profile cluster/lsf -n

# Run test target (UniProt only)
pixi run snakemake --profile cluster/lsf test
```

### Production Runs

```bash
# Build seed models and run iteration 1
pixi run snakemake --profile cluster/lsf iter1

# After reviewing results, run iteration 2
pixi run snakemake --profile cluster/lsf iter2

# Continue with auto-curation and final models
pixi run snakemake --profile cluster/lsf
```

### Specific Targets

```bash
# Download all databases
pixi run snakemake --profile cluster/lsf download_all

# Build seed models
pixi run snakemake --profile cluster/lsf build_seeds
```

## Resource Configuration

Resources are configured at two levels:

1. **Default Resources** (`cluster/lsf/config.yaml`): Applied to all rules
   - `lsf_queue: "rna"` - Default queue
   - `runtime: 120` - Default 2-hour limit (minutes)
   - `mem_mb: 8` - Default 8GB memory (NOTE: values are in GB despite the name)

2. **Rule-Specific Overrides** (`set-resources`): Custom settings per rule
   - HMM searches get 24 hours and 16GB memory
   - Downloads get minimal resources
   - Reporting gets moderate resources

### Memory Note

The LSF executor uses `mem_mb` parameter but values are specified in **GB**, not MB. For example:
- `mem_mb=8` means 8GB
- `mem_mb=16` means 16GB

## Resource Specifications

### Small Jobs (downloads, model building)
- Runtime: 10-60 minutes
- Memory: 2-4GB

### Medium Jobs (alignment processing)
- Runtime: 30-60 minutes
- Memory: 8-16GB

### Large Jobs (HMM searches)
- Runtime: Up to 24 hours
- Memory: 16GB

## Monitoring Jobs

### Check LSF Queue

```bash
# Your jobs
bjobs

# All jobs with details
bjobs -l

# Job history
bhist -a
```

### Check Logs

```bash
# Snakemake log
tail -f .snakemake/log/*.log

# Individual rule logs
ls logs/

# LSF job output
ls .snakemake/lsf_logs/
```

### Cancel Jobs

```bash
# Cancel specific job
bkill <JOB_ID>

# Cancel all your jobs
bkill 0
```

## Troubleshooting

### Jobs Not Submitting

1. **Check project is set**: Verify `lsf_project` in `cluster/lsf/config.yaml`
2. **Check queue access**: `bqueues`
3. **Check your limits**: `blimits`

### Jobs Failing

1. **Check LSF logs**: `.snakemake/lsf_logs/`
2. **Check rule logs**: `logs/<rule>.log`
3. **Increase resources**: Edit `set-resources` in `cluster/lsf/config.yaml`

### Out of Memory

Increase `mem_mb` (in GB) for the specific rule:

```yaml
set-resources:
  - hmmsearch:mem_mb=32  # Increase to 32GB
```

### Timeout

Increase `runtime` (in minutes) for the specific rule:

```yaml
set-resources:
  - hmmsearch:runtime=2880  # Increase to 48 hours
```

## References

- [Snakemake LSF Executor Plugin](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/lsf.html)
- [Snakemake Profiles](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles)
