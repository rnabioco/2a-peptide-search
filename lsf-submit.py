#!/usr/bin/env python3
"""
LSF submission script for Snakemake
Based on the Snakemake LSF executor pattern
"""

import os
import sys
import argparse
from pathlib import Path

def parse_jobscript(jobscript):
    """Parse the jobscript to extract job information"""
    job_properties = {}
    
    with open(jobscript, 'r') as f:
        for line in f:
            if line.startswith('# properties'):
                # Extract job properties from the comment
                import json
                props_str = line.split('=', 1)[1].strip()
                try:
                    job_properties = json.loads(props_str)
                except:
                    pass
                break
    
    return job_properties

def create_lsf_script(jobscript, job_properties, cluster_config):
    """Create LSF submission script"""
    
    # Get rule name
    rule = job_properties.get('rule', 'unknown')
    wildcards = job_properties.get('wildcards', {})
    
    # Get cluster configuration for this rule
    rule_config = cluster_config.get(rule, cluster_config.get('__default__', {}))
    
    # Format job name with wildcards
    jobname = rule_config.get('jobname', rule)
    if wildcards:
        wildcard_str = '.'.join([f"{k}-{v}" for k, v in wildcards.items()])
        jobname = f"{rule}.{wildcard_str}"
    
    # Create LSF script content
    lsf_script = f"""#!/bin/bash
#BSUB -J {jobname}
#BSUB -q {rule_config.get('queue', 'normal')}
#BSUB -n {rule_config.get('cores', 1)}
#BSUB -R "rusage[mem={rule_config.get('memory', '4GB')}]"
#BSUB -W {rule_config.get('walltime', '2:00')}
#BSUB -o logs/{jobname}.%J.out
#BSUB -e logs/{jobname}.%J.err

# Make sure logs directory exists
mkdir -p logs

# Source any necessary environment setup
# module load hmmer/3.3.2  # Uncomment and modify as needed
# module load easel/0.48   # Uncomment and modify as needed

# Execute the actual job
{jobscript}
"""
    
    return lsf_script

def main():
    parser = argparse.ArgumentParser(description='Submit Snakemake job to LSF')
    parser.add_argument('jobscript', help='Path to the jobscript')
    parser.add_argument('--cluster-config', help='Path to cluster configuration')
    
    args = parser.parse_args()
    
    # Load cluster configuration
    cluster_config = {}
    if args.cluster_config and os.path.exists(args.cluster_config):
        import yaml
        with open(args.cluster_config, 'r') as f:
            cluster_config = yaml.safe_load(f)
    
    # Parse job properties
    job_properties = parse_jobscript(args.jobscript)
    
    # Create LSF script
    lsf_script = create_lsf_script(args.jobscript, job_properties, cluster_config)
    
    # Write LSF script to temporary file
    lsf_script_path = f"tmp_lsf_{os.getpid()}.sh"
    with open(lsf_script_path, 'w') as f:
        f.write(lsf_script)
    
    # Submit job
    os.system(f"bsub < {lsf_script_path}")
    
    # Clean up
    os.remove(lsf_script_path)

if __name__ == '__main__':
    main()