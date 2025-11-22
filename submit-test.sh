#!/bin/bash
#SBATCH --job-name=2a-test
#SBATCH --partition=amilan
#SBATCH --account=amc-general
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --qos=normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=/scratch/alpine/%u/2a-peptide-search/logs/test_%j.out
#SBATCH --error=/scratch/alpine/%u/2a-peptide-search/logs/test_%j.err

# 2A Peptide Search Pipeline - Quick Test
# ========================================
# This script runs a quick test of the pipeline with UniProt only.
# Good for verifying the setup before running the full pipeline.
#
# Usage:
#   sbatch submit-test.sh

set -euo pipefail

# ============================================================================
# Configure output directories for Alpine scratch filesystem
# ============================================================================

# Use Alpine's fast scratch filesystem for all outputs
# Format: /scratch/alpine/<username>/<project-name>
export SCRATCH_DIR="/scratch/alpine/${USER}/2a-peptide-search"
export RESULTS_DIR="${SCRATCH_DIR}/results"
export DATA_DIR="${SCRATCH_DIR}/data"
export LOGS_DIR="${SCRATCH_DIR}/logs"

# Create base directories
mkdir -p "$SCRATCH_DIR"
mkdir -p "$LOGS_DIR"

echo "=========================================="
echo "2A Peptide Search Pipeline - Test Run"
echo "=========================================="
echo "Started: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Scratch dir: $SCRATCH_DIR"
echo ""

echo "Running test with UniProt only..."
echo ""

# Dry run first
echo "Dry run:"
pixi run -e default snakemake \
    --snakefile pipeline/workflow/Snakefile \
    --configfile pipeline/workflow/config/config.yaml \
    --profile pipeline/cluster/slurm \
    -n test

echo ""
echo "Executing test target..."
echo ""

pixi run -e default snakemake \
    --snakefile pipeline/workflow/Snakefile \
    --configfile pipeline/workflow/config/config.yaml \
    --profile pipeline/cluster/slurm \
    test

EXIT_CODE=$?

echo ""
echo "=========================================="
if [ $EXIT_CODE -eq 0 ]; then
    echo "SUCCESS: Test completed"
    echo "Check results/searches/uniprot/ for output"
else
    echo "ERROR: Test failed with exit code $EXIT_CODE"
fi
echo "Finished: $(date)"
echo "=========================================="

exit $EXIT_CODE
