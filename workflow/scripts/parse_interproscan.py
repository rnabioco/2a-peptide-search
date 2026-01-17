#!/usr/bin/env python3
"""
Parse InterProScan TSV output into standardized domain format.

InterProScan TSV output columns:
1. protein_acc - Protein accession
2. seq_md5 - Sequence MD5 digest
3. seq_len - Sequence length
4. sig_db - Signature database
5. sig_acc - Signature accession
6. sig_desc - Signature description
7. start_loc - Start location
8. stop_loc - Stop location
9. score - E-value or score
10. status - Status (T for matches)
11. date - Run date
12. ipr_annot_acc - InterPro accession
13. ipr_annot_desc - InterPro description
14. go_annot - GO annotations (optional)
15. pathways - Pathway annotations (optional)
"""

import gzip
import csv
import sys


def main():
    """Parse InterProScan TSV to standardized format."""
    input_tsv = snakemake.input.tsv
    output_tsv = snakemake.output.domains

    # InterProScan TSV column names
    ipr_columns = [
        'protein_acc', 'seq_md5', 'seq_len', 'sig_db', 'sig_acc',
        'sig_desc', 'start_loc', 'stop_loc', 'score', 'status',
        'date', 'ipr_annot_acc', 'ipr_annot_desc', 'go_annot', 'pathways'
    ]

    # Output columns (subset of input, renamed for clarity)
    out_columns = [
        'protein_acc', 'sig_db', 'sig_acc', 'sig_desc',
        'query_start', 'query_end', 'dom_evalue',
        'ipr_acc', 'ipr_desc'
    ]

    opener = gzip.open if input_tsv.endswith('.gz') else open

    with opener(input_tsv, 'rt') as f_in, open(output_tsv, 'w', newline='') as f_out:
        writer = csv.DictWriter(f_out, fieldnames=out_columns, delimiter='\t')
        writer.writeheader()

        reader = csv.DictReader(f_in, fieldnames=ipr_columns, delimiter='\t')
        row_count = 0

        for row in reader:
            # Skip empty lines
            if not row.get('protein_acc'):
                continue

            out_row = {
                'protein_acc': row['protein_acc'],
                'sig_db': row['sig_db'],
                'sig_acc': row['sig_acc'],
                'sig_desc': row['sig_desc'],
                'query_start': row['start_loc'],
                'query_end': row['stop_loc'],
                'dom_evalue': row['score'],
                'ipr_acc': row.get('ipr_annot_acc', '-'),
                'ipr_desc': row.get('ipr_annot_desc', '-')
            }
            writer.writerow(out_row)
            row_count += 1

    print(f"Parsed {row_count} domain annotations", file=sys.stderr)


if __name__ == '__main__':
    main()
