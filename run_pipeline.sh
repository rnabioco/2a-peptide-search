#!/usr/bin/env bash

# Wrapper script for running the 2A peptide search Snakemake pipeline

set -e

# Default parameters
CORES=1
DRYRUN=false
CLUSTER=false
UNLOCK=false
CLEAN=false
HELP=false

# Function to show help
show_help() {
    cat << EOF
Usage: $0 [OPTIONS]

Options:
    -c, --cores CORES       Number of cores to use (default: 1)
    -d, --dry-run          Perform a dry run (show what would be executed)
    -s, --cluster          Submit jobs to cluster using SLURM
    -u, --unlock           Unlock the working directory
    -x, --clean            Clean intermediate files
    -h, --help             Show this help message

Examples:
    # Run locally with 4 cores
    $0 -c 4
    
    # Dry run to see what would be executed
    $0 -d
    
    # Submit to SLURM cluster
    $0 -s
    
    # Clean intermediate files
    $0 -x
    
    # Unlock directory (useful if pipeline was interrupted)
    $0 -u

EOF
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -c|--cores)
            CORES="$2"
            shift 2
            ;;
        -d|--dry-run)
            DRYRUN=true
            shift
            ;;
        -s|--cluster)
            CLUSTER=true
            shift
            ;;
        -u|--unlock)
            UNLOCK=true
            shift
            ;;
        -x|--clean)
            CLEAN=true
            shift
            ;;
        -h|--help)
            HELP=true
            shift
            ;;
        *)
            echo "Unknown option: $1"
            show_help
            exit 1
            ;;
    esac
done

# Show help if requested
if [[ "$HELP" == "true" ]]; then
    show_help
    exit 0
fi

# Create necessary directories
mkdir -p logs results tmp

# Handle special operations
if [[ "$UNLOCK" == "true" ]]; then
    echo "Unlocking Snakemake working directory..."
    snakemake --unlock
    exit 0
fi

if [[ "$CLEAN" == "true" ]]; then
    echo "Cleaning intermediate files..."
    snakemake clean
    exit 0
fi

# Build the snakemake command
SNAKEMAKE_CMD="snakemake"

# Add cores
SNAKEMAKE_CMD="$SNAKEMAKE_CMD --cores $CORES"

# Add dry run if requested
if [[ "$DRYRUN" == "true" ]]; then
    SNAKEMAKE_CMD="$SNAKEMAKE_CMD --dry-run"
fi

# Add cluster configuration if requested
if [[ "$CLUSTER" == "true" ]]; then
    SNAKEMAKE_CMD="$SNAKEMAKE_CMD --cluster-config cluster.yaml"
    SNAKEMAKE_CMD="$SNAKEMAKE_CMD --cluster 'sbatch --partition={cluster.partition} --time={cluster.time} --mem={cluster.mem} --cpus-per-task={cluster.cpus} --job-name={cluster.name} --output={cluster.output} --error={cluster.error}'"
    SNAKEMAKE_CMD="$SNAKEMAKE_CMD --jobs 50"  # Allow up to 50 concurrent jobs
fi

# Add other useful flags
SNAKEMAKE_CMD="$SNAKEMAKE_CMD --printshellcmds --reason --stats snakemake_stats.txt"

# Run the pipeline
echo "Running: $SNAKEMAKE_CMD"
eval $SNAKEMAKE_CMD