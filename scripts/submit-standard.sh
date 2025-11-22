#!/bin/bash
#SBATCH --job-name=2a-standard
#SBATCH --comment="2a-standard"
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --account=amc-general
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=/scratch/alpine/%u/2a-peptide-search/logs/standard_%j.out
#SBATCH --error=/scratch/alpine/%u/2a-peptide-search/logs/standard_%j.err

# Standard 2A Peptide Search Pipeline
# ====================================
# Identifies eukaryotic/viral 2A peptides using curated seed alignments
# for class-1 and class-2 2A peptides.
#
# Usage:
#   sbatch scripts/submit-standard.sh [target]
#
# Examples:
#   sbatch scripts/submit-standard.sh                 # Run full pipeline
#   sbatch scripts/submit-standard.sh test            # Test with UniProt only
#   sbatch scripts/submit-standard.sh iter1           # Run iteration 1
#   sbatch scripts/submit-standard.sh download_all    # Download all databases
#
# Configuration:
#   - Edit cluster/slurm/config.yaml to set your SLURM account
#   - Edit workflow/config/config.yaml to select databases and set parameters

set -euo pipefail

# ============================================================================
# Configuration
# ============================================================================

# Target to run (default: all)
TARGET="${1:-all}"

# Snakemake configuration
SNAKEFILE="workflow/Snakefile"
CONFIGFILE="workflow/config/config.yaml"
PROFILE="cluster/slurm"

# ============================================================================
# Configure output directories for Alpine scratch filesystem
# ============================================================================

# Use Alpine's fast scratch filesystem for all outputs
# Format: /scratch/alpine/<username>/<project-name>
# Note: Scratch has 90-day purge policy - move important results to /projects after completion
export SCRATCH_DIR="/scratch/alpine/${USER}/2a-peptide-search"
export RESULTS_DIR="${SCRATCH_DIR}/results"
export DATA_DIR="${SCRATCH_DIR}/data"
export LOGS_DIR="${SCRATCH_DIR}/logs"

# Create base directories
mkdir -p "$SCRATCH_DIR"
mkdir -p "$LOGS_DIR"

# ============================================================================
# Environment Setup
# ============================================================================

echo "=========================================="
echo "Standard 2A Peptide Search Pipeline"
echo "=========================================="
echo "Started: $(date)"
echo "Target: $TARGET"
echo "Snakefile: $SNAKEFILE"
echo "Config: $CONFIGFILE"
echo "Profile: $PROFILE"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURMD_NODENAME"
echo "Scratch dir: $SCRATCH_DIR"
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
    -n "$TARGET"; then
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
    "$TARGET"

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
    echo "Next steps:"
    echo "  - Check results in results/ directory"
    echo "  - Review logs in logs/ directory"
    echo "  - For manual curation checkpoint:"
    echo "    1. Review alignments in results/alignments/iter2/"
    echo "    2. Curate and save to results/alignments/final/"
    echo "    3. Touch checkpoint: touch results/checkpoints/iter2.curated"
    echo "    4. Resubmit: sbatch submit-slurm.sh"
else
    echo "ERROR: Pipeline failed with exit code $EXIT_CODE"
    echo ""
    echo "Troubleshooting:"
    echo "  - Check logs/orchestrator_${SLURM_JOB_ID}.err for errors"
    echo "  - Check .snakemake/log/ for detailed Snakemake logs"
    echo "  - Check individual rule logs in logs/"
    echo "  - View SLURM job status: squeue -u \$USER"
    echo "  - View failed job details: sacct -j <JOBID> --format=JobID,State,ExitCode,Reason"
fi

exit $EXIT_CODE
