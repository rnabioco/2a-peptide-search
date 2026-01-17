#!/bin/bash
#SBATCH --job-name=2a-prokaryotic
#SBATCH --comment="2a-prokaryotic"
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --account=amc-general
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=/scratch/alpine/%u/2a-peptide-search/logs/prokaryotic_%j.out
#SBATCH --error=/scratch/alpine/%u/2a-peptide-search/logs/prokaryotic_%j.err

# Prokaryotic 2A-like Peptide Discovery Pipeline (SLURM)
# =======================================================
# Discovers prokaryotic ribosomal stalling peptides using GP-motif
# extraction, domain annotation, and pattern clustering.
#
# Usage:
#   sbatch scripts/submit-prokaryotic-slurm.sh [target]
#
# Examples:
#   sbatch scripts/submit-prokaryotic-slurm.sh                      # Run full prokaryotic discovery
#   sbatch scripts/submit-prokaryotic-slurm.sh extract_gp_motifs    # Extract GP motifs
#   sbatch scripts/submit-prokaryotic-slurm.sh cluster_gp_motifs    # Cluster motifs
#   sbatch scripts/submit-prokaryotic-slurm.sh prokaryotic_report   # Generate report
#
# Configuration:
#   - Edit cluster/slurm/config.yaml to set your SLURM account
#   - Edit workflow/config/config-prokaryotic.yaml to configure phage databases and parameters
#
# See workflow/PROKARYOTIC-DISCOVERY.md for detailed documentation

set -euo pipefail

# Change to submission directory (project root)
cd "$SLURM_SUBMIT_DIR"

# ============================================================================
# Configuration
# ============================================================================

# Target to run (default: prokaryotic_discovery)
TARGET="${1:-prokaryotic_discovery}"

# Snakemake configuration
SNAKEFILE="workflow/Snakefile"
CONFIGFILE="workflow/config/config-prokaryotic.yaml"
PROFILE="cluster/slurm"

# ============================================================================
# Configure output directories for Alpine scratch filesystem
# ============================================================================

# Use Alpine's fast scratch filesystem for all outputs
# Format: /scratch/alpine/<username>/<project-name>
# Note: Scratch has 90-day purge policy - move important results to /projects after completion
# Use $USER as-is (contains @ on Alpine cluster)
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
echo "Prokaryotic 2A-like Discovery Pipeline"
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
    echo "    4. Resubmit: sbatch scripts/submit-prokaryotic-slurm.sh"
else
    echo "ERROR: Pipeline failed with exit code $EXIT_CODE"
    echo ""
    echo "Troubleshooting:"
    echo "  - Check logs/prokaryotic_${SLURM_JOB_ID}.err for errors"
    echo "  - Check .snakemake/log/ for detailed Snakemake logs"
    echo "  - Check individual rule logs in logs/"
    echo "  - View SLURM job status: squeue -u \$USER"
    echo "  - View failed job details: sacct -j <JOBID> --format=JobID,State,ExitCode,Reason"
fi

exit $EXIT_CODE
