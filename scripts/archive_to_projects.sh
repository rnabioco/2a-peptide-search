#!/bin/bash
#
# Archive important results from scratch to /projects
#
# Usage:
#   ./scripts/archive_to_projects.sh [scratch_dir] [projects_dir]
#
# Example:
#   ./scripts/archive_to_projects.sh \
#     /scratch/alpine/$USER/2a-peptide-search \
#     /projects/$USER/2a-peptide-search

set -euo pipefail

# Get directories
SCRATCH_DIR="${1:-/scratch/alpine/$USER/2a-peptide-search}"
PROJECTS_DIR="${2:-/projects/$USER/2a-peptide-search}"

echo "=========================================="
echo "Archive Results to Projects"
echo "=========================================="
echo "Source (scratch): $SCRATCH_DIR"
echo "Destination (projects): $PROJECTS_DIR"
echo ""

# Create destination
mkdir -p "$PROJECTS_DIR"

# Archive final models
echo "Archiving final models..."
if [ -d "$SCRATCH_DIR/results/models/final" ]; then
    cp -r "$SCRATCH_DIR/results/models/final" "$PROJECTS_DIR/"
    echo "  ✓ Copied models/final/"
fi

# Archive all models (optional)
echo "Archiving all models..."
if [ -d "$SCRATCH_DIR/results/models" ]; then
    mkdir -p "$PROJECTS_DIR/models"
    cp -r "$SCRATCH_DIR/results/models"/* "$PROJECTS_DIR/models/"
    echo "  ✓ Copied all models/"
fi

# Archive reports
echo "Archiving reports..."
if [ -d "$SCRATCH_DIR/results/reports" ]; then
    mkdir -p "$PROJECTS_DIR/reports"
    cp "$SCRATCH_DIR/results/reports"/*.html "$PROJECTS_DIR/reports/" 2>/dev/null || true
    cp "$SCRATCH_DIR/results/reports"/*.tsv "$PROJECTS_DIR/reports/" 2>/dev/null || true
    echo "  ✓ Copied reports/"
fi

# Archive curated alignments
echo "Archiving curated alignments..."
if [ -d "$SCRATCH_DIR/results/alignments/final" ]; then
    mkdir -p "$PROJECTS_DIR/alignments"
    cp -r "$SCRATCH_DIR/results/alignments/final" "$PROJECTS_DIR/alignments/"
    echo "  ✓ Copied alignments/final/"
fi

# Archive checkpoints
echo "Archiving checkpoints..."
if [ -d "$SCRATCH_DIR/results/checkpoints" ]; then
    cp -r "$SCRATCH_DIR/results/checkpoints" "$PROJECTS_DIR/"
    echo "  ✓ Copied checkpoints/"
fi

# Compress and archive search results (optional - large)
echo "Compressing search results (this may take a while)..."
if [ -d "$SCRATCH_DIR/scratch/searches" ]; then
    cd "$SCRATCH_DIR/scratch"
    tar -czf "$PROJECTS_DIR/searches.tar.gz" searches/
    echo "  ✓ Compressed searches/ → searches.tar.gz"
fi

# Compress and archive alignments (optional - large)
echo "Compressing alignment results..."
if [ -d "$SCRATCH_DIR/scratch/alignments" ]; then
    cd "$SCRATCH_DIR/scratch"
    tar -czf "$PROJECTS_DIR/alignments.tar.gz" alignments/
    echo "  ✓ Compressed alignments/ → alignments.tar.gz"
fi

# Copy config for reproducibility
echo "Copying configuration..."
if [ -f "$SCRATCH_DIR/workflow/config.yaml" ]; then
    cp "$SCRATCH_DIR/workflow/config.yaml" "$PROJECTS_DIR/"
    echo "  ✓ Copied config.yaml"
fi

# Create manifest
echo "Creating manifest..."
cat > "$PROJECTS_DIR/MANIFEST.txt" <<EOF
2A Peptide Search - Archived Results
=====================================
Archived: $(date)
Source: $SCRATCH_DIR
Destination: $PROJECTS_DIR

Contents:
- models/final/          Final HMM models
- models/                All models (seed, refined, final)
- reports/               HTML reports and statistics
- alignments/final/      Curated alignments
- checkpoints/           Manual curation checkpoints
- searches.tar.gz        Compressed search results
- alignments.tar.gz      Compressed alignment intermediates
- config.yaml            Pipeline configuration

To extract compressed archives:
  tar -xzf searches.tar.gz
  tar -xzf alignments.tar.gz

Note: Original files remain on scratch until purge (90 days)
EOF

echo "  ✓ Created MANIFEST.txt"

# Summary
echo ""
echo "=========================================="
echo "Archive Complete"
echo "=========================================="
echo "Destination: $PROJECTS_DIR"
du -sh "$PROJECTS_DIR"
echo ""
echo "Important files archived:"
ls -lh "$PROJECTS_DIR"
echo ""
echo "Note: Original files remain on scratch for 90 days"
echo "      You can delete them manually if needed:"
echo "      rm -rf $SCRATCH_DIR"
