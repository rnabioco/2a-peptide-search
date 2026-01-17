#!/usr/bin/env python3
"""
Extract protein sequences from HMM search hits for domain annotation.

This script extracts sequence IDs from curated alignments and corresponding
tblout files, outputting them for downstream domain analysis.
"""

import gzip
import sys
from Bio import AlignIO, SearchIO
from pathlib import Path


def main():
    """Extract hit sequences from alignment and tblout files."""
    # Snakemake provides inputs/outputs via snakemake object
    alignment_file = snakemake.input.alignment
    tblout_files = snakemake.input.tblout
    output_fasta = snakemake.output.fasta
    output_ids = snakemake.output.ids

    # Collect all hit IDs and their E-values from tblout files
    hits = {}
    for tblout_file in tblout_files:
        opener = gzip.open if tblout_file.endswith('.gz') else open
        with opener(tblout_file, 'rt') as f:
            for result in SearchIO.parse(f, 'hmmer3-tab'):
                for hit in result.hits:
                    # Use full ID including coordinates
                    hits[hit.id] = hit.evalue

    # Read alignment and extract sequences
    aln = AlignIO.read(alignment_file, 'stockholm')

    with open(output_fasta, 'w') as fasta_out, open(output_ids, 'w') as ids_out:
        for record in aln:
            # Remove gaps for FASTA output
            seq_nogaps = str(record.seq).replace('-', '').replace('.', '')
            fasta_out.write(f">{record.id}\n{seq_nogaps}\n")
            ids_out.write(f"{record.id}\n")

    print(f"Extracted {len(aln)} sequences", file=sys.stderr)


if __name__ == '__main__':
    main()
