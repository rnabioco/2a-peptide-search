#!/bin/bash
#BSUB -J 2a-prokaryotic
#BSUB -q rna
#BSUB -n 1
#BSUB -R "rusage[mem=4GB]"
#BSUB -W 24:00
#BSUB -o logs/prokaryotic_%J.out
#BSUB -e logs/prokaryotic_%J.err

# Prokaryotic 2A-like Peptide Discovery Pipeline (LSF)
# =====================================================
# Discovers prokaryotic ribosomal stalling peptides using GP-motif
# extraction, domain annotation, and pattern clustering.
#
# Usage:
#   bsub < scripts/submit-prokaryotic-lsf.sh
#   bsub -J 2a-extract scripts/submit-prokaryotic-lsf.sh extract_gp_motifs
#
# Examples:
#   bsub < scripts/submit-prokaryotic-lsf.sh                      # Run full prokaryotic discovery
#   bsub scripts/submit-prokaryotic-lsf.sh extract_gp_motifs      # Extract GP motifs
#   bsub scripts/submit-prokaryotic-lsf.sh cluster_gp_motifs      # Cluster motifs
#   bsub scripts/submit-prokaryotic-lsf.sh prokaryotic_report     # Generate report
#
# Configuration:
#   - Edit cluster/lsf/config.yaml to set your LSF project and queue
#   - Edit workflow/config/config-prokaryotic.yaml to configure phage databases and parameters
#
# See workflow/PROKARYOTIC-DISCOVERY.md for detailed documentation

set -euo pipefail

# Change to submission directory (project root)
# LSF uses LSB_SUBCWD for the directory where bsub was run
cd "${LSB_SUBCWD:-.}"

# ============================================================================
# Configuration
# ============================================================================

# Target to run (default: prokaryotic_discovery)
TARGET="${1:-prokaryotic_discovery}"

# Snakemake configuration
SNAKEFILE="workflow/Snakefile"
CONFIGFILE="workflow/config/config-prokaryotic.yaml"
PROFILE="cluster/lsf"

# ============================================================================
# Configure output directories (via direnv)
# ============================================================================

# Source direnv configuration for host-specific paths
# This sets RESULTS_DIR, DATA_DIR, LOGS_DIR based on hostname
if [[ -f .envrc ]]; then
    source .envrc
fi

# Ensure directories exist
mkdir -p "${RESULTS_DIR:-results}"
mkdir -p "${LOGS_DIR:-logs}"
mkdir -p "logs"

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
echo "Job ID: ${LSB_JOBID:-N/A}"
echo "Host: $(hostname)"
echo "Working dir: $(pwd)"
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
    echo "  - Check results in ${RESULTS_DIR:-results}/ directory"
    echo "  - Review logs in ${LOGS_DIR:-logs}/ directory"
    echo "  - For manual curation checkpoint:"
    echo "    1. Review alignments in ${RESULTS_DIR:-results}/alignments/iter2/"
    echo "    2. Curate and save to ${RESULTS_DIR:-results}/alignments/final/"
    echo "    3. Touch checkpoint: touch ${RESULTS_DIR:-results}/checkpoints/iter2.curated"
    echo "    4. Resubmit: bsub < scripts/submit-prokaryotic-lsf.sh"
else
    echo "ERROR: Pipeline failed with exit code $EXIT_CODE"
    echo ""
    echo "Troubleshooting:"
    echo "  - Check logs/prokaryotic_${LSB_JOBID:-unknown}.err for errors"
    echo "  - Check .snakemake/log/ for detailed Snakemake logs"
    echo "  - Check individual rule logs in ${LOGS_DIR:-logs}/"
    echo "  - View LSF job status: bjobs"
    echo "  - View job history: bhist -l ${LSB_JOBID:-<JOBID>}"
fi

exit $EXIT_CODE
