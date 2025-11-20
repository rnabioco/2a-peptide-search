#!/usr/bin/env python3
"""
Extract sequences from FASTA database with flanking regions based on
alignment coordinates from Stockholm format file.
"""

import gzip
import sys
from collections import namedtuple

import click
from Bio import AlignIO, SeqIO


Id = namedtuple('Id', ['fasta_id', 'start', 'end'])


def parse_id(seq_id, descrip):
    """Parse sequence ID with coordinates: tr|X4YD30|X4YD30_9PICO/692-706"""
    fasta_id, coords = seq_id.split('/')
    start, end = map(int, coords.split('-'))
    return Id(fasta_id, start, end)


@click.command()
@click.argument('aln_fname', type=click.Path(exists=True))
@click.argument('db_fname', type=click.Path(exists=True))
@click.option('--flank', default=20, help='Number of flanking residues to include')
def main(aln_fname, db_fname, flank):
    """Extract sequences with flanking regions from alignment hits."""

    # Parse alignment and extract IDs with coordinates
    align = AlignIO.read(aln_fname, "stockholm")
    ids = {}
    for record in align:
        seq_id = parse_id(record.id, record.description)
        ids[seq_id.fasta_id] = seq_id

    # Extract sequences from database with flanking regions
    with gzip.open(db_fname, 'rt') as fasta:
        for record in SeqIO.parse(fasta, 'fasta'):
            if record.id in ids:
                seq_id = ids[record.id]
                start = max(0, seq_id.start - flank)
                end = min(len(record.seq), seq_id.end + flank)
                seq = record.seq[start:end]
                print(f'>{seq_id.fasta_id}/{start}-{end}\n{seq}')


if __name__ == '__main__':
    main()
