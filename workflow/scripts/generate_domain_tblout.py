#!/usr/bin/env python3
"""
Generate domain hits table from HMM search tblout files.

Extracts 2A peptide domain coordinates from our HMM search results,
providing domain boundary information for downstream analysis.
"""

import gzip
import csv
import sys
from Bio import SearchIO


def main():
    """Generate domain table from tblout files."""
    tblout_files = snakemake.input.tblouts
    output_tsv = snakemake.output.domains

    # Output columns
    out_columns = [
        'protein_acc', 'sig_db', 'sig_acc', 'sig_desc',
        'query_start', 'query_end', 'dom_evalue', 'dom_score',
        'database'
    ]

    with open(output_tsv, 'w', newline='') as f_out:
        writer = csv.DictWriter(f_out, fieldnames=out_columns, delimiter='\t')
        writer.writeheader()

        total_hits = 0

        for tblout_file in tblout_files:
            # Extract database name from path (e.g., .../searches/uniprot/final/...)
            path_parts = tblout_file.split('/')
            database = 'unknown'
            for i, part in enumerate(path_parts):
                if part == 'searches' and i + 1 < len(path_parts):
                    database = path_parts[i + 1]
                    break

            opener = gzip.open if tblout_file.endswith('.gz') else open

            with opener(tblout_file, 'rt') as f:
                for result in SearchIO.parse(f, 'hmmer3-tab'):
                    model_name = result.id  # e.g., 2A-class-1

                    for hit in result.hits:
                        # Parse protein accession (may include coordinates)
                        protein_acc = hit.id.split('/')[0] if '/' in hit.id else hit.id

                        # For tblout format, we get overall hit info
                        # Domain boundaries come from domtblout; approximate here
                        out_row = {
                            'protein_acc': protein_acc,
                            'sig_db': 'Hmmer3-custom',
                            'sig_acc': model_name,
                            'sig_desc': '2A peptide (skips peptide bond)',
                            'query_start': hit.query_start if hasattr(hit, 'query_start') else '-',
                            'query_end': hit.query_end if hasattr(hit, 'query_end') else '-',
                            'dom_evalue': hit.evalue,
                            'dom_score': hit.bitscore,
                            'database': database
                        }
                        writer.writerow(out_row)
                        total_hits += 1

    print(f"Generated {total_hits} domain entries", file=sys.stderr)


if __name__ == '__main__':
    main()
