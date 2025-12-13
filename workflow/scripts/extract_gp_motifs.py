#!/usr/bin/env python3
"""
Extract all GP-containing sequences from prokaryotic proteomes.

Extracts sequence context around GP motifs for downstream analysis.
"""

import gzip
import re
from collections import defaultdict

import click
import pandas as pd
from Bio import SeqIO


def find_gp_motifs(sequence, upstream=30, downstream=15):
    """Find all GP motif positions in a sequence with context."""
    motifs = []

    # Find all GP positions
    for match in re.finditer(r'GP', str(sequence)):
        gp_pos = match.start()

        # Extract context
        start = max(0, gp_pos - upstream)
        end = min(len(sequence), gp_pos + 2 + downstream)

        context_seq = str(sequence[start:end])

        # Relative position of GP in context
        gp_in_context = gp_pos - start

        motifs.append({
            'gp_position': gp_pos,
            'context_sequence': context_seq,
            'gp_relative_pos': gp_in_context,
            'upstream_length': gp_pos - start,
            'downstream_length': end - gp_pos - 2
        })

    return motifs


@click.command()
@click.option('--fasta', required=True, help='Input FASTA file (can be gzipped)')
@click.option('--motifs-out', required=True, help='Output TSV of GP motifs')
@click.option('--sequences-out', required=True, help='Output FASTA of sequences with GP')
@click.option('--upstream', default=30, help='Residues upstream of GP')
@click.option('--downstream', default=15, help='Residues downstream of GP')
def main(fasta, motifs_out, sequences_out, upstream, downstream):
    """Extract GP motifs from proteome."""

    click.echo(f"Extracting GP motifs from {fasta}")
    click.echo(f"Context: {upstream} upstream, {downstream} downstream")

    all_motifs = []
    gp_protein_ids = set()

    # Open input FASTA
    if fasta.endswith('.gz'):
        handle = gzip.open(fasta, 'rt')
    else:
        handle = open(fasta, 'r')

    total_proteins = 0
    proteins_with_gp = 0
    total_gp_motifs = 0

    for record in SeqIO.parse(handle, 'fasta'):
        total_proteins += 1

        # Find GP motifs
        motifs = find_gp_motifs(record.seq, upstream, downstream)

        if motifs:
            proteins_with_gp += 1
            total_gp_motifs += len(motifs)
            gp_protein_ids.add(record.id)

            for i, motif in enumerate(motifs):
                all_motifs.append({
                    'protein_id': record.id,
                    'protein_description': record.description,
                    'protein_length': len(record.seq),
                    'gp_index': i + 1,  # 1-indexed
                    'gp_position': motif['gp_position'],
                    'context_sequence': motif['context_sequence'],
                    'gp_relative_pos': motif['gp_relative_pos'],
                    'upstream_length': motif['upstream_length'],
                    'downstream_length': motif['downstream_length'],
                    # Additional features
                    'n_terminal_distance': motif['gp_position'],
                    'c_terminal_distance': len(record.seq) - motif['gp_position'] - 2,
                })

        if total_proteins % 10000 == 0:
            click.echo(f"Processed {total_proteins:,} proteins...")

    handle.close()

    # Save motifs table
    df = pd.DataFrame(all_motifs)
    df.to_csv(motifs_out, sep='\t', index=False, compression='gzip')

    # Save sequences with GP
    if fasta.endswith('.gz'):
        handle = gzip.open(fasta, 'rt')
    else:
        handle = open(fasta, 'r')

    with gzip.open(sequences_out, 'wt') as out:
        for record in SeqIO.parse(handle, 'fasta'):
            if record.id in gp_protein_ids:
                SeqIO.write(record, out, 'fasta')

    handle.close()

    # Print summary
    click.echo("\n" + "="*60)
    click.echo("GP Motif Extraction Summary")
    click.echo("="*60)
    click.echo(f"Total proteins: {total_proteins:,}")

    if total_proteins == 0:
        click.echo("ERROR: No proteins found in input file!", err=True)
        click.echo("Please check that the FASTA file exists and contains sequences.", err=True)
        raise ValueError("No proteins found in input file")

    click.echo(f"Proteins with GP: {proteins_with_gp:,} ({proteins_with_gp/total_proteins*100:.2f}%)")
    click.echo(f"Total GP motifs: {total_gp_motifs:,}")

    if proteins_with_gp > 0:
        click.echo(f"Average GP per protein with GP: {total_gp_motifs/proteins_with_gp:.2f}")

    click.echo(f"\nOutput files:")
    click.echo(f"  Motifs: {motifs_out}")
    click.echo(f"  Sequences: {sequences_out}")


if __name__ == '__main__':
    main()
