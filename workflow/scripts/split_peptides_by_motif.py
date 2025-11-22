#!/usr/bin/env python3
"""
Split known stalling peptides into separate FASTA files by motif type.
"""

from pathlib import Path
from collections import defaultdict
from Bio import SeqIO
import click


@click.command()
@click.option('--fasta', required=True, help='Input FASTA file')
@click.option('--output-dir', required=True, help='Output directory for motif FASTAs')
@click.option('--motif-list', required=True, help='Output file listing all motifs')
def main(fasta, output_dir, motif_list):
    """Split known stalling peptides into separate FASTA files by motif type."""

    input_fasta = Path(fasta)
    output_dir_path = Path(output_dir)
    motif_list_file = Path(motif_list)

    # Create output directory
    output_dir_path.mkdir(parents=True, exist_ok=True)

    # Group sequences by motif type
    motif_sequences = defaultdict(list)

    for record in SeqIO.parse(input_fasta, "fasta"):
        # Extract motif from header (format: MOTIF|organism|gene|...)
        motif = record.id.split("|")[0]
        motif_sequences[motif].append(record)

    # Write separate FASTA files for each motif
    motif_counts = {}
    for motif, sequences in sorted(motif_sequences.items()):
        output_fasta = output_dir_path / f"{motif}.fasta"
        SeqIO.write(sequences, output_fasta, "fasta")
        motif_counts[motif] = len(sequences)
        click.echo(f"Wrote {len(sequences)} sequences for motif {motif}")

    # Write motif list for downstream rules
    with open(motif_list_file, "w") as f:
        for motif in sorted(motif_sequences.keys()):
            f.write(f"{motif}\n")

    click.echo(f"\nSummary:")
    click.echo(f"Total motif types: {len(motif_sequences)}")
    click.echo(f"Motif distribution:")
    for motif, count in sorted(motif_counts.items(), key=lambda x: -x[1]):
        click.echo(f"  {motif}: {count} sequences")


if __name__ == "__main__":
    main()
