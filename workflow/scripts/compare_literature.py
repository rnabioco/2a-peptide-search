#!/usr/bin/env python3
"""
Compare HMM search results against known literature 2A sequences.

This script cross-references HMM hits with the original literature data,
including activity measurements where available.
"""

import csv
import gzip
import sys
from collections import defaultdict
from Bio import SearchIO


def parse_tblout(tblout_file, class_name):
    """Parse hmmsearch tblout file and return hit IDs."""
    hits = {}
    for result in SearchIO.parse(tblout_file, 'hmmer3-tab'):
        for hit in result.hits:
            if hit.domain_included_num >= 1:
                hits[hit.id] = {
                    'evalue': hit.evalue,
                    'score': hit.bitscore,
                    'class': class_name
                }
    return hits


def main():
    """Compare literature sequences with HMM hits."""
    tblout_c1 = snakemake.input.tblout_c1
    tblout_c2 = snakemake.input.tblout_c2
    sequences_file = snakemake.input.sequences
    output_file = snakemake.output.comparison

    # Parse HMM search results
    hits_c1 = parse_tblout(tblout_c1, 'class-1')
    hits_c2 = parse_tblout(tblout_c2, 'class-2')

    # Load original sequence data
    ref_seqs = {}
    row_num = 0
    with open(sequences_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            row_num += 1
            name = row.get('Name', f'seq_{row_num}')
            name_clean = name.replace(' ', '_').replace('/', '-').replace('|', '-')
            anal_id = row.get('Analysis ID', 'NA')
            if anal_id == 'N/A':
                anal_id = 'NA'

            seq_id = f"{name_clean}|{row_num}|{anal_id}"
            ref_seqs[seq_id] = {
                'name': name,
                'sequence': row.get('Sequence', ''),
                'class_original': row.get('ClassA/B', ''),
                'source': row.get('Host/virus/synthetic', ''),
                'activity_tested': row.get('Tested activity yes/no', ''),
                'activity': row.get('Activity', ''),
                'old_new': row.get('Old/new', '')
            }

    # Generate comparison output
    output_cols = [
        'seq_id', 'name', 'sequence', 'class_original', 'source',
        'activity_tested', 'activity', 'old_new',
        'hmm_class1_hit', 'hmm_class1_evalue', 'hmm_class1_score',
        'hmm_class2_hit', 'hmm_class2_evalue', 'hmm_class2_score',
        'classification'
    ]

    with open(output_file, 'w', newline='') as f_out:
        writer = csv.DictWriter(f_out, fieldnames=output_cols, delimiter='\t')
        writer.writeheader()

        hits_found = 0
        for seq_id, ref_info in ref_seqs.items():
            c1_hit = hits_c1.get(seq_id)
            c2_hit = hits_c2.get(seq_id)

            # Determine classification
            if c1_hit and not c2_hit:
                classification = 'class-1'
            elif c2_hit and not c1_hit:
                classification = 'class-2'
            elif c1_hit and c2_hit:
                # Both hit - assign to better scoring model
                if c1_hit['score'] > c2_hit['score']:
                    classification = 'class-1 (dual-hit)'
                else:
                    classification = 'class-2 (dual-hit)'
            else:
                classification = 'no-hit'

            if c1_hit or c2_hit:
                hits_found += 1

            out_row = {
                'seq_id': seq_id,
                'name': ref_info['name'],
                'sequence': ref_info['sequence'],
                'class_original': ref_info['class_original'],
                'source': ref_info['source'],
                'activity_tested': ref_info['activity_tested'],
                'activity': ref_info['activity'],
                'old_new': ref_info['old_new'],
                'hmm_class1_hit': 'Yes' if c1_hit else 'No',
                'hmm_class1_evalue': c1_hit['evalue'] if c1_hit else '',
                'hmm_class1_score': c1_hit['score'] if c1_hit else '',
                'hmm_class2_hit': 'Yes' if c2_hit else 'No',
                'hmm_class2_evalue': c2_hit['evalue'] if c2_hit else '',
                'hmm_class2_score': c2_hit['score'] if c2_hit else '',
                'classification': classification
            }
            writer.writerow(out_row)

    total = len(ref_seqs)
    print(f"Compared {total} sequences: {hits_found} hits ({100*hits_found/total:.1f}%)", file=sys.stderr)


if __name__ == '__main__':
    main()
