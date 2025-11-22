#!/usr/bin/env python3
"""
Split known stalling peptides into separate FASTA files by motif type.
"""

from pathlib import Path
from collections import defaultdict
from Bio import SeqIO
import sys


def main():
    input_fasta = Path(snakemake.input.fasta)
    output_dir = Path(snakemake.output.fastas)
    motif_list_file = Path(snakemake.output.motif_list)

    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)

    # Group sequences by motif type
    motif_sequences = defaultdict(list)

    for record in SeqIO.parse(input_fasta, "fasta"):
        # Extract motif from header (format: MOTIF|organism|gene|...)
        motif = record.id.split("|")[0]
        motif_sequences[motif].append(record)

    # Write separate FASTA files for each motif
    motif_counts = {}
    for motif, sequences in sorted(motif_sequences.items()):
        output_fasta = output_dir / f"{motif}.fasta"
        SeqIO.write(sequences, output_fasta, "fasta")
        motif_counts[motif] = len(sequences)
        print(f"Wrote {len(sequences)} sequences for motif {motif}", file=sys.stderr)

    # Write motif list for downstream rules
    with open(motif_list_file, "w") as f:
        for motif in sorted(motif_sequences.keys()):
            f.write(f"{motif}\n")

    print(f"\nSummary:", file=sys.stderr)
    print(f"Total motif types: {len(motif_sequences)}", file=sys.stderr)
    print(f"Motif distribution:", file=sys.stderr)
    for motif, count in sorted(motif_counts.items(), key=lambda x: -x[1]):
        print(f"  {motif}: {count} sequences", file=sys.stderr)


if __name__ == "__main__":
    main()
