#!/usr/bin/env python3
"""
Collect statistics from HMMER search results for the final report.

This script parses tblout files from HMM searches and generates a summary
table with hit counts, E-value distributions, and sequence statistics per
database and peptide class.

Called by Snakemake rule collect_statistics in report.smk.
"""

import gzip
import re
import sys
from pathlib import Path

import pandas as pd
from Bio import SearchIO


def parse_tblout(tblout_path):
    """Parse HMMER tblout file and extract hit statistics.

    Args:
        tblout_path: Path to gzipped tblout file

    Returns:
        List of dicts with hit information
    """
    hits = []

    # Handle gzipped files
    opener = gzip.open if str(tblout_path).endswith('.gz') else open

    with opener(tblout_path, 'rt') as f:
        for result in SearchIO.parse(f, 'hmmer3-tab'):
            for hit in result.hits:
                hits.append({
                    'target_id': hit.id,
                    'evalue': hit.evalue,
                    'score': hit.bitscore,
                    'bias': hit.bias,
                })

    return hits


def extract_metadata(tblout_path):
    """Extract database and peptide class from file path.

    Expected path pattern:
        .../searches/{database}/final/2A-{peptide_class}.tblout.gz

    Args:
        tblout_path: Path to tblout file

    Returns:
        Tuple of (database, peptide_class)
    """
    path = Path(tblout_path)

    # Extract peptide class from filename (e.g., "2A-class-1.tblout.gz")
    match = re.match(r'2A-(class-\d+)\.tblout\.gz', path.name)
    peptide_class = match.group(1) if match else 'unknown'

    # Extract database from path (parent of 'final' directory)
    parts = path.parts
    try:
        final_idx = parts.index('final')
        database = parts[final_idx - 1]
    except (ValueError, IndexError):
        database = 'unknown'

    return database, peptide_class


def compute_statistics(hits, database, peptide_class):
    """Compute summary statistics for a set of hits.

    Args:
        hits: List of hit dicts from parse_tblout
        database: Database name
        peptide_class: Peptide class (class-1 or class-2)

    Returns:
        Dict with summary statistics
    """
    if not hits:
        return {
            'database': database,
            'peptide_class': peptide_class,
            'hit_count': 0,
            'unique_proteins': 0,
            'mean_evalue': None,
            'median_evalue': None,
            'min_evalue': None,
            'max_evalue': None,
            'mean_score': None,
            'median_score': None,
        }

    df = pd.DataFrame(hits)
    unique_proteins = df['target_id'].str.split('/').str[0].nunique()

    return {
        'database': database,
        'peptide_class': peptide_class,
        'hit_count': len(hits),
        'unique_proteins': unique_proteins,
        'mean_evalue': df['evalue'].mean(),
        'median_evalue': df['evalue'].median(),
        'min_evalue': df['evalue'].min(),
        'max_evalue': df['evalue'].max(),
        'mean_score': df['score'].mean(),
        'median_score': df['score'].median(),
    }


def main(tblout_files, output_path, log_path=None):
    """Main function to collect statistics from all tblout files.

    Args:
        tblout_files: List of paths to tblout files
        output_path: Path for output TSV file
        log_path: Optional path for log file
    """
    # Set up logging
    log_file = open(log_path, 'w') if log_path else sys.stderr

    def log(msg):
        print(msg, file=log_file)

    log(f"Collecting statistics from {len(tblout_files)} tblout files")

    all_stats = []

    for tblout_path in tblout_files:
        log(f"Processing: {tblout_path}")

        database, peptide_class = extract_metadata(tblout_path)
        log(f"  Database: {database}, Class: {peptide_class}")

        try:
            hits = parse_tblout(tblout_path)
            log(f"  Found {len(hits)} hits")

            stats = compute_statistics(hits, database, peptide_class)
            all_stats.append(stats)

        except Exception as e:
            log(f"  ERROR: {e}")
            # Add empty stats on error
            all_stats.append({
                'database': database,
                'peptide_class': peptide_class,
                'hit_count': 0,
                'unique_proteins': 0,
                'mean_evalue': None,
                'median_evalue': None,
                'min_evalue': None,
                'max_evalue': None,
                'mean_score': None,
                'median_score': None,
            })

    # Create output DataFrame
    df = pd.DataFrame(all_stats)

    # Sort by database and peptide class
    df = df.sort_values(['database', 'peptide_class'])

    # Write output
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, sep='\t', index=False)

    log(f"Wrote statistics to {output_path}")
    log(f"Total databases: {df['database'].nunique()}")
    log(f"Total hits: {df['hit_count'].sum()}")

    if log_path:
        log_file.close()


# Snakemake script interface
if __name__ == '__main__':
    # When called from Snakemake, use snakemake object
    if 'snakemake' in dir():
        main(
            tblout_files=snakemake.input.tblouts,
            output_path=snakemake.output.stats,
            log_path=snakemake.log[0] if snakemake.log else None,
        )
    else:
        # CLI fallback for testing
        import click

        @click.command()
        @click.argument('tblout_files', nargs=-1, type=click.Path(exists=True))
        @click.option('--output', '-o', required=True, help='Output TSV file')
        @click.option('--log', '-l', help='Log file')
        def cli(tblout_files, output, log):
            """Collect statistics from HMMER tblout files."""
            main(tblout_files, output, log)

        cli()
