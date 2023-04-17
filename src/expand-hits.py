#! /usr/bin/env python3

import pdb
import gzip

import click

from collections import namedtuple

from Bio import AlignIO, SeqIO

@click.command()
@click.argument('aln_fname')
@click.argument('db_fname')

def main(aln_fname, db_fname):

    align = AlignIO.read(aln_fname, "stockholm")

    ids = {}
    for record in align:
        id = parse_id(record.id, record.description)
        ids[id.fasta_id] = id
    
    with gzip.open(db_fname, 'rt') as fasta:
        for record in SeqIO.parse(fasta, 'fasta'):
            if record.id in ids:
                id = ids[record.id]
                seq = record.seq[id.start-20:id.end+20]
                print(f'>{id.orig_id}\n{seq}')

Id = namedtuple('Id', ['fasta_id', 'orig_id', 'uniprot_id', 'org', 'start', 'end'])
def parse_id(id, descrip):
    ''' tr|X4YD30|X4YD30_9PICO/692-706 '''
 
    fasta_id, coords = id.split('/')
    orig_id = fasta_id + ' ' + descrip
    start, end = map(int, coords.split('-'))
    uniprot_id, org = fasta_id.split('|')[2].split('_')

    return Id(fasta_id, orig_id, uniprot_id, org, start, end)

if __name__ == '__main__':
    main()
