#!/bin/bash
# Submit Snakemake pipeline to LSF cluster

# Create logs directory
mkdir -p logs

# Run Snakemake with LSF cluster execution
snakemake \
    --cluster "python3 lsf-submit.py {jobscript} --cluster-config cluster.yaml" \
    --cluster-config cluster.yaml \
    --jobs 50 \
    --latency-wait 60 \
    --keep-going \
    --rerun-incomplete \
    --printshellcmds \
    "$@"