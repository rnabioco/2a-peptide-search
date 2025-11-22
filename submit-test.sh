#!/bin/bash
#SBATCH --job-name=2a-test
#SBATCH --partition=amilan
#SBATCH --account=amc-general
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --qos=normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/test_%j.out
#SBATCH --error=logs/test_%j.err

# 2A Peptide Search Pipeline - Quick Test
# ========================================
# This script runs a quick test of the pipeline with UniProt only.
# Good for verifying the setup before running the full pipeline.
#
# Usage:
#   sbatch submit-test.sh

set -euo pipefail

mkdir -p logs

echo "=========================================="
echo "2A Peptide Search Pipeline - Test Run"
echo "=========================================="
echo "Started: $(date)"
echo "Job ID: $SLURM_JOB_ID"
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
