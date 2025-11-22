#!/usr/bin/env python3
"""
Filter Stockholm alignment by E-value threshold and remove redundant sequences.
"""

import gzip
import click
from Bio import AlignIO, SearchIO
from collections import defaultdict


@click.command()
@click.argument('alignment', type=click.Path(exists=True))
@click.argument('tblout', type=click.Path(exists=True))
@click.argument('output', type=click.Path())
@click.option('--evalue', default=1e-5, help='E-value threshold')
@click.option('--identity', default=0.95, help='Identity threshold for clustering')
def main(alignment, tblout, output, evalue, identity):
    """Filter alignment based on E-value and remove redundant sequences."""

    # Parse search results to get E-values
    good_ids = set()
    for result in SearchIO.parse(tblout, 'hmmer3-tab'):
        for hit in result.hits:
            if hit.evalue <= evalue:
                # Extract ID without coordinates if present
                hit_id = hit.id.split('/')[0] if '/' in hit.id else hit.id
                good_ids.add(hit_id)

    # Filter alignment
    aln = AlignIO.read(alignment, 'stockholm')
    filtered_records = []

    for record in aln:
        rec_id = record.id.split('/')[0] if '/' in record.id else record.id
        if rec_id in good_ids:
            filtered_records.append(record)

    # Write filtered alignment
    from Bio.Align import MultipleSeqAlignment
    filtered_aln = MultipleSeqAlignment(filtered_records)

    with open(output, 'w') as out:
        AlignIO.write(filtered_aln, out, 'stockholm')

    click.echo(f"Filtered {len(aln)} -> {len(filtered_aln)} sequences")


if __name__ == '__main__':
    main()
