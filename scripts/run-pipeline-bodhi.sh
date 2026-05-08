#!/bin/bash
#SBATCH --job-name=2a-peptide
#SBATCH --partition=rna
#SBATCH --account=rbi
#SBATCH --time=3-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/pipeline_%j.out
#SBATCH --error=logs/pipeline_%j.err

# 2A Peptide Search Pipeline - SLURM Submission Script (amc-bodhi)
# ================================================================
# Identifies eukaryotic/viral 2A peptides via iterative HMM search
# against UniProt and reference proteomes.
#
# Usage:
#   sbatch scripts/run-pipeline-bodhi.sh [target] [snakemake_args...]
#
# Examples:
#   sbatch scripts/run-pipeline-bodhi.sh                       # Full pipeline (rule all)
#   sbatch scripts/run-pipeline-bodhi.sh test                  # Quick UniProt-only test
#   sbatch scripts/run-pipeline-bodhi.sh download_all          # Download databases
#   sbatch scripts/run-pipeline-bodhi.sh iter1                 # Up through iter1
#   sbatch scripts/run-pipeline-bodhi.sh iter2                 # Through iter2 checkpoint
#   sbatch scripts/run-pipeline-bodhi.sh all --rerun-triggers=mtime
#
# Configuration:
#   - Edit cluster/slurm-bodhi/config.yaml to adjust resource allocations
#   - Edit workflow/config.yaml for analysis parameters and database URLs
#
# Cluster info:
#   - CPU jobs: partition=rna, account=rbi (this pipeline is CPU-only)

set -euo pipefail

# Ensure pixi is on PATH (not inherited by compute nodes)
export PATH="$HOME/.pixi/bin:$PATH"

# Change to submission directory (project root)
cd "$SLURM_SUBMIT_DIR"

# ============================================================================
# Configuration
# ============================================================================

if [ $# -eq 0 ]; then
    SNAKEMAKE_ARGS=("all")
else
    SNAKEMAKE_ARGS=("$@")
fi

SNAKEFILE="workflow/Snakefile"
CONFIGFILE="workflow/config.yaml"
PROFILE="cluster/slurm-bodhi"

# ============================================================================
# Setup
# ============================================================================

mkdir -p logs logs/slurm

echo "=========================================="
echo "2A Peptide Search Pipeline (amc-bodhi)"
echo "=========================================="
echo "Started: $(date)"
echo "Snakemake args: ${SNAKEMAKE_ARGS[*]}"
echo "Snakefile: $SNAKEFILE"
echo "Config: $CONFIGFILE"
echo "Profile: $PROFILE"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURMD_NODENAME"
echo ""

# ============================================================================
# Dry Run
# ============================================================================

echo "Running dry-run to check workflow..."
echo ""

if pixi run -e default snakemake \
    --snakefile "$SNAKEFILE" \
    --configfile "$CONFIGFILE" \
    --profile "$PROFILE" \
    -n "${SNAKEMAKE_ARGS[@]}"; then
    echo ""
    echo "Dry-run successful. Proceeding with execution..."
    echo ""
else
    echo ""
    echo "ERROR: Dry-run failed. Check workflow configuration."
    echo ""
    exit 1
fi

# ============================================================================
# Execute Pipeline
# ============================================================================

echo "=========================================="
echo "Starting pipeline execution"
echo "=========================================="
echo ""

pixi run -e default snakemake \
    --snakefile "$SNAKEFILE" \
    --configfile "$CONFIGFILE" \
    --profile "$PROFILE" \
    "${SNAKEMAKE_ARGS[@]}"

EXIT_CODE=$?

# ============================================================================
# Summary
# ============================================================================

echo ""
echo "=========================================="
echo "Pipeline Completed"
echo "=========================================="
echo "Finished: $(date)"
echo "Exit code: $EXIT_CODE"
echo ""

if [ $EXIT_CODE -eq 0 ]; then
    echo "SUCCESS: Pipeline completed successfully"
    echo ""
    echo "Output locations:"
    echo "  Results: scratch/results/"
    echo "  Data: scratch/data/"
    echo "  Logs: logs/"
    echo ""
    echo "Manual curation checkpoint (between iter2 and final):"
    echo "  1. Review alignments in scratch/results/alignments/iter2/"
    echo "  2. Curate and save to scratch/results/alignments/final/"
    echo "  3. Touch checkpoint: touch scratch/results/checkpoints/iter2.curated"
    echo "  4. Resubmit: sbatch scripts/run-pipeline-bodhi.sh"
else
    echo "ERROR: Pipeline failed with exit code $EXIT_CODE"
    echo ""
    echo "Troubleshooting:"
    echo "  - Check logs/pipeline_${SLURM_JOB_ID}.err for errors"
    echo "  - Check .snakemake/log/ for detailed Snakemake logs"
    echo "  - View SLURM job status: squeue -u \$USER"
    echo "  - View failed job details: sacct -j <JOBID> --format=JobID,State,ExitCode,Reason"
fi

exit $EXIT_CODE
