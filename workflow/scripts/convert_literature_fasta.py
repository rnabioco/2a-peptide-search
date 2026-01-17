#!/usr/bin/env python3
"""
Convert literature 2A sequences TSV to FASTA format for hmmsearch.

Creates sequence IDs that encode the row number and analysis ID for
later cross-referencing with the original data.
"""

import csv
import sys


def main():
    """Convert TSV to FASTA."""
    input_tsv = snakemake.input.tsv
    output_fasta = snakemake.output.fasta

    with open(input_tsv, 'r') as f_in, open(output_fasta, 'w') as f_out:
        reader = csv.DictReader(f_in, delimiter='\t')
        row_num = 0

        for row in reader:
            row_num += 1

            name = row.get('Name', f'seq_{row_num}')
            # Clean name for FASTA header (replace problematic chars)
            name_clean = name.replace(' ', '_').replace('/', '-').replace('|', '-')

            anal_id = row.get('Analysis ID', 'NA')
            if anal_id == 'N/A':
                anal_id = 'NA'

            seq = row.get('Sequence', '')
            if not seq:
                continue

            # Create FASTA ID that encodes row number for cross-referencing
            # Format: name|row_num|analysis_id
            seq_id = f"{name_clean}|{row_num}|{anal_id}"

            f_out.write(f">{seq_id}\n{seq}\n")

    print(f"Converted {row_num} sequences to FASTA", file=sys.stderr)


if __name__ == '__main__':
    main()
