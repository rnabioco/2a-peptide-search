# SLURM Executor Profile for 2A Peptide Search Pipeline

This directory contains the Snakemake profile configuration for running the 2A peptide search pipeline on the CU Boulder Alpine cluster using the SLURM executor.

## Setup

### 1. Install SLURM Executor Plugin

```bash
# Activate your conda/mamba environment
conda activate snakemake

# Install the SLURM executor plugin
conda install -c conda-forge -c bioconda snakemake-executor-plugin-slurm

# Or using pip
pip install snakemake-executor-plugin-slurm
```

### 2. Configure Your Alpine Account

Edit `cluster/slurm/config.yaml` and set your Alpine allocation account:

```yaml
default-resources:
  slurm_account: "your-account-here"  # e.g., "amc-general" or "ucb123_asc1"
```

**Important:** You must set this before running the pipeline, or jobs will fail to submit.

### 3. Configure Output Directories

Edit `workflow/config.yaml` and set paths appropriate for Alpine:

```yaml
# Use scratch for temporary/intermediate files (fast I/O, 90-day purge)
# Use /projects for permanent storage

# Example for user jay
scratch_dir: "/scratch/alpine/jay@xsede.org/2a-peptide-search"
projects_dir: "/projects/jay@xsede.org/2a-peptide-search"
```

## Usage

### Quick Test Run

Test the pipeline with UniProt only (fast, ~1-2 hours):

```bash
# Dry run to see execution plan
snakemake --profile cluster/slurm -n

# Run test target (UniProt only)
snakemake --profile cluster/slurm test
```

### Production Runs

Run the full pipeline with all databases:

```bash
# Build seed models and run iteration 1
snakemake --profile cluster/slurm iter1

# After reviewing results, run iteration 2
snakemake --profile cluster/slurm iter2

# After manual curation checkpoint
# 1. Review alignments in results/alignments/iter2/
# 2. Curate and save to results/alignments/final/*.curated.sto
# 3. Create checkpoint: touch results/checkpoints/iter2.curated
# 4. Continue pipeline
snakemake --profile cluster/slurm
```

### Specific Targets

```bash
# Download all databases
snakemake --profile cluster/slurm download_all

# Build seed models
snakemake --profile cluster/slurm build_seeds

# Search specific database
snakemake --profile cluster/slurm results/searches/uniprot/seed/2A-class-1.hmmsearch.gz
```

## How It Works

### Orchestrator Pattern

The SLURM executor uses an "orchestrator" pattern:

1. **Orchestrator Job**: Submit initial Snakemake job to SLURM
2. **Worker Jobs**: Snakemake submits each rule as a separate SLURM job
3. **Job Management**: Orchestrator monitors workers and manages DAG execution

### Resource Management

Resources are configured at three levels:

1. **Default Resources** (`cluster/slurm/config.yaml`): Applied to all rules
   - `slurm_partition: "amilan"` - CPU partition
   - `runtime: 120` - Default 2-hour limit
   - `mem_mb: 8000` - Default 8GB memory
   - `cpus_per_task: 4` - Default 4 CPUs

2. **Rule-Specific Overrides** (`set-resources`): Custom settings per rule
   - HMM searches get 24 hours and 16GB memory
   - Downloads get minimal resources
   - Reporting gets moderate resources

3. **Workflow Config** (`workflow/config.yaml`): Pipeline parameters
   - Database selection
   - E-value thresholds
   - Number of threads for specific tools

## Resource Specifications

### Small Jobs (downloads, model building)
- Runtime: 10-60 minutes
- Memory: 2-4GB
- CPUs: 1-2

### Medium Jobs (alignment processing)
- Runtime: 30-60 minutes
- Memory: 8-16GB
- CPUs: 2-4

### Large Jobs (HMM searches)
- Runtime: Up to 24 hours
- Memory: 16GB
- CPUs: 12

**Note**: UniParc and MGnify searches can take many hours. Consider disabling in `workflow/config.yaml` for testing.

## Monitoring Jobs

### Check SLURM Queue

```bash
# Your jobs
squeue -u $USER

# Specific job details
scontrol show job <JOB_ID>

# Jobs for this pipeline
squeue -u $USER --name=2a-*
```

### Check Logs

```bash
# Snakemake orchestrator log
tail -f .snakemake/log/*.log

# Individual rule logs (from workflow)
ls logs/

# SLURM job logs
ls .snakemake/slurm_logs/
```

### Cancel Jobs

```bash
# Cancel specific job
scancel <JOB_ID>

# Cancel all your jobs
scancel -u $USER

# Cancel by name pattern
scancel --name=2a-
```

## Typical Workflow

### 1. Initial Setup
```bash
# Configure account in cluster/slurm/config.yaml
# Configure databases in workflow/config.yaml (start with UniProt only)
# Note: ALL outputs go to scratch by default
```

### 2. Dry Run
```bash
snakemake --profile cluster/slurm -n
```

### 3. Test Run
```bash
# Test with UniProt only (~1-2 hours)
snakemake --profile cluster/slurm test
```

### 4. First Iteration
```bash
# Run seed searches and build refined models
snakemake --profile cluster/slurm iter1

# This will:
# - Build seed HMMs from curated alignments
# - Search configured databases
# - Filter and merge alignments
# - Build refined models
```

### 5. Second Iteration
```bash
# Search with refined models
snakemake --profile cluster/slurm iter2

# This will:
# - Search with iter1 refined models
# - Create checkpoint for manual curation
```

### 6. Manual Curation
```bash
# Review alignments
jalview results/alignments/iter2/2A-class-1.merged.sto

# Curate and save
cp results/alignments/iter2/2A-class-1.merged.sto \
   results/alignments/final/2A-class-1.curated.sto

# Create checkpoint
touch results/checkpoints/iter2.curated
```

### 7. Final Models and Report
```bash
# Build final models and generate report
snakemake --profile cluster/slurm

# This will:
# - Build final production HMMs
# - Run final comprehensive searches
# - Generate HTML report
```

### 8. Archive Results to /projects
```bash
# After pipeline completion, archive important results
./scripts/archive_to_projects.sh

# This will:
# - Copy final models, reports, curated alignments to /projects
# - Compress large intermediate files (searches, alignments)
# - Create manifest of archived files
# - Leave originals on scratch (90-day purge)

# Manual archiving (alternative):
cd /scratch/alpine/$USER/2a-peptide-search
mkdir -p /projects/$USER/2a-peptide-search

# Copy essential files
cp -r results/models/final /projects/$USER/2a-peptide-search/
cp results/reports/*.html /projects/$USER/2a-peptide-search/

# Archive large files
tar -czf searches.tar.gz scratch/searches
mv searches.tar.gz /projects/$USER/2a-peptide-search/
```

## Troubleshooting

### Jobs Not Submitting

1. **Check account is set**: Verify `slurm_account` in `cluster/slurm/config.yaml`
2. **Check allocation**: `sacctmgr show assoc user=$USER format=account,qos`
3. **Check partition access**: `sinfo -p amilan`

### Jobs Failing

1. **Check SLURM logs**: `.snakemake/slurm_logs/rule-<rule>-<jobid>.out`
2. **Check rule logs**: `logs/<rule>.log`
3. **Increase resources**: Edit `set-resources` in `cluster/slurm/config.yaml`

### Out of Memory

Increase `mem_mb` for the specific rule:

```yaml
set-resources:
  hmmsearch:
    mem_mb: 32000  # Increase to 32GB
```

### Timeout

Increase `runtime` (in minutes) for the specific rule:

```yaml
set-resources:
  hmmsearch:
    runtime: 2880  # Increase to 48 hours
```

### Download Failures

Database URLs may change. Update in `workflow/config.yaml`:

```yaml
databases:
  uniprot:
    url: "https://ftp.uniprot.org/pub/databases/uniprot/current_release/..."
```

## Performance Tips

### Start Small
- Begin with `uniprot` and `reference_proteomes` only
- Disable `uniparc` and `mgnify` in `workflow/config.yaml`
- Test with dry run first: `snakemake --profile cluster/slurm -n`

### Use Scratch Storage
- Set `data/` directory to `/scratch/alpine/$USER/2a-peptide-search/data`
- Faster I/O than `/projects`
- Remember to move important results to `/projects` (90-day purge policy)

### Monitor Resource Usage
```bash
# After jobs complete, check efficiency
seff <JOB_ID>

# Adjust resources based on actual usage
```

### Parallel Execution
- Pipeline automatically parallelizes across databases and peptide classes
- SLURM profile allows up to 100 concurrent jobs
- Adjust in `cluster/slurm/config.yaml`: `jobs: 100`

## References

- [Snakemake SLURM Executor Plugin](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html)
- [Alpine Documentation](https://curc.readthedocs.io/en/latest/clusters/alpine/)
- [Alpine SLURM Guide](https://curc.readthedocs.io/en/latest/running-jobs/batch-jobs.html)
- [Snakemake Profiles](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles)
