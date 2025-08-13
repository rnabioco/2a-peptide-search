# Running the 2A Peptide Search Pipeline with LSF

This pipeline has been configured to work with IBM Spectrum LSF (Load Sharing Facility) batch scheduler.

## Prerequisites

1. **Snakemake**: Install via conda/mamba
   ```bash
   conda install -c bioconda snakemake
   ```

2. **Required software**: Ensure these are available on your cluster
   - HMMER suite (hmmbuild, hmmsearch)
   - Easel utilities (esl-alimerge)
   - Python 3 with PyYAML

3. **Database files**: Update paths in `config.yaml` to point to your database files

## Quick Start

### Method 1: Using the run script (recommended)
```bash
# Make the run script executable
chmod +x run-lsf.sh

# Submit the entire pipeline
./run-lsf.sh

# Or run specific targets
./run-lsf.sh --until merge_alignments
```

### Method 2: Using Snakemake profile
```bash
# Make the LSF submit script executable  
chmod +x lsf-submit.py

# Run with the LSF profile
snakemake --profile profiles/lsf
```

### Method 3: Direct Snakemake command
```bash
snakemake \
    --cluster "python3 lsf-submit.py {jobscript} --cluster-config cluster.yaml" \
    --cluster-config cluster.yaml \
    --jobs 50 \
    --latency-wait 60 \
    --keep-going \
    --rerun-incomplete
```

## Configuration

### Cluster Resources (`cluster.yaml`)
Modify resource requirements per rule:
- `queue`: LSF queue name
- `memory`: Memory requirement (e.g., "8GB") 
- `cores`: Number of CPU cores
- `walltime`: Runtime limit (format: "HH:MM")

### Pipeline Settings (`config.yaml`)
- Update database file paths
- Adjust HMMER parameters (E-value, CPU cores)

## Key LSF Features

- **Job names**: Automatically formatted with rule and wildcards
- **Resource requests**: Memory, cores, and walltime per job
- **Queue selection**: Different queues for different job types
- **Log files**: Stored in `logs/` directory with job ID
- **Module loading**: Easily add module load commands in `lsf-submit.py`

## Monitoring Jobs

```bash
# Check job status
bjobs

# Check specific job details
bjobs -l <job_id>

# View job output
bpeek <job_id>

# Check queue information
bqueues
```

## Troubleshooting

1. **Permission errors**: Make sure scripts are executable
   ```bash
   chmod +x run-lsf.sh lsf-submit.py
   ```

2. **Module loading**: Edit `lsf-submit.py` to add necessary module loads
3. **Queue names**: Update queue names in `cluster.yaml` to match your system
4. **File paths**: Ensure all paths in `config.yaml` are correct
5. **Log files**: Check `logs/` directory for job-specific error messages

## Example Workflow

1. Test locally first:
   ```bash
   snakemake -n  # Dry run
   ```

2. Submit to cluster:
   ```bash
   ./run-lsf.sh
   ```

3. Monitor progress:
   ```bash
   bjobs  # Check running jobs
   snakemake --summary  # Check completed rules
   ```

## Advanced Usage

### Run specific rules
```bash
./run-lsf.sh hmmsearch_single
```

### Force re-run of specific files
```bash
./run-lsf.sh --forcerun merge_alignments
```

### Dry run to see what would be executed
```bash
./run-lsf.sh -n
```