#!/usr/bin/env python3
"""
Extract sequences with C-terminal arrest motifs from prokaryotic proteomes.

Focuses on known bacterial ribosomal arrest motif patterns:
- RAGP, RAPG (SecM-like)
- QAPP, QGPP
- HAPP, HGPP
- RAPP, RPPP

Extracts upstream context for building HMMs that can be validated
against known stalling peptides.
"""

import gzip
import re
from collections import defaultdict

import click
import pandas as pd
from Bio import SeqIO


# Known arrest motif patterns (4-mer)
ARREST_MOTIFS = [
    "RAGP",  # SecM-like (E. coli)
    "RAPG",  # SecM variant
    "QAPP",  # B. subtilis stalling
    "QGPP",  # stalling variant
    "HAPP",  # stalling variant
    "HGPP",  # stalling variant
    "RAPP",  # stalling variant
    "RPPP",  # stalling variant
]


def find_arrest_motifs(sequence, upstream=40, max_distance_from_cterm=50):
    """
    Find arrest motifs near the C-terminus of a sequence.

    Args:
        sequence: Protein sequence string
        upstream: Residues to extract upstream of motif
        max_distance_from_cterm: Maximum distance from C-terminus to consider

    Returns:
        List of dicts with motif info and context
    """
    motifs = []
    seq_str = str(sequence).upper()
    seq_len = len(seq_str)

    for motif in ARREST_MOTIFS:
        # Find all occurrences of this motif
        for match in re.finditer(motif, seq_str):
            motif_start = match.start()
            motif_end = match.end()

            # Calculate distance from C-terminus
            c_term_distance = seq_len - motif_end

            # Only consider motifs near C-terminus
            if c_term_distance <= max_distance_from_cterm:
                # Extract upstream context
                context_start = max(0, motif_start - upstream)
                # Include the motif itself plus a few residues downstream
                context_end = min(seq_len, motif_end + 10)

                context_seq = seq_str[context_start:context_end]

                motifs.append({
                    "motif_type": motif,
                    "motif_position": motif_start,
                    "c_terminal_distance": c_term_distance,
                    "context_sequence": context_seq,
                    "motif_in_context": motif_start - context_start,
                    "upstream_length": motif_start - context_start,
                    "downstream_length": context_end - motif_end,
                    "is_at_cterm": c_term_distance == 0,
                })

    return motifs


@click.command()
@click.option("--fasta", required=True, help="Input FASTA file (can be gzipped)")
@click.option("--motifs-out", required=True, help="Output TSV of arrest motifs")
@click.option("--sequences-out", required=True, help="Output FASTA of context sequences")
@click.option("--upstream", default=40, help="Residues upstream of motif to extract")
@click.option(
    "--max-cterm-distance",
    default=50,
    help="Maximum distance from C-terminus to consider",
)
def main(fasta, motifs_out, sequences_out, upstream, max_cterm_distance):
    """Extract arrest motifs from proteome."""

    click.echo(f"Extracting C-terminal arrest motifs from {fasta}")
    click.echo(f"Motifs searched: {', '.join(ARREST_MOTIFS)}")
    click.echo(f"Context: {upstream} upstream, max {max_cterm_distance} from C-terminus")

    all_motifs = []
    motif_counts = defaultdict(int)

    # Open input FASTA
    if fasta.endswith(".gz"):
        handle = gzip.open(fasta, "rt")
    else:
        handle = open(fasta)

    total_proteins = 0
    proteins_with_motifs = 0

    for record in SeqIO.parse(handle, "fasta"):
        total_proteins += 1

        # Find arrest motifs
        motifs = find_arrest_motifs(
            record.seq, upstream=upstream, max_distance_from_cterm=max_cterm_distance
        )

        if motifs:
            proteins_with_motifs += 1

            for i, motif in enumerate(motifs):
                motif_counts[motif["motif_type"]] += 1

                all_motifs.append({
                    "protein_id": record.id,
                    "protein_description": record.description,
                    "protein_length": len(record.seq),
                    "motif_index": i + 1,
                    "motif_type": motif["motif_type"],
                    "motif_position": motif["motif_position"],
                    "c_terminal_distance": motif["c_terminal_distance"],
                    "context_sequence": motif["context_sequence"],
                    "motif_in_context": motif["motif_in_context"],
                    "upstream_length": motif["upstream_length"],
                    "downstream_length": motif["downstream_length"],
                    "is_at_cterm": motif["is_at_cterm"],
                })

        if total_proteins % 100000 == 0:
            click.echo(f"Processed {total_proteins:,} proteins...")

    handle.close()

    if not all_motifs:
        click.echo("WARNING: No arrest motifs found!")
        # Create empty outputs
        pd.DataFrame().to_csv(motifs_out, sep="\t", index=False, compression="gzip")
        with gzip.open(sequences_out, "wt") as out:
            pass
        return

    # Save motifs table
    df = pd.DataFrame(all_motifs)
    df.to_csv(motifs_out, sep="\t", index=False, compression="gzip")

    # Save context sequences as FASTA (for clustering/alignment)
    with gzip.open(sequences_out, "wt") as out:
        for _, row in df.iterrows():
            # Create unique ID: protein_motiftype_position
            seq_id = f"{row['protein_id']}|{row['motif_type']}|pos{row['motif_position']}"
            out.write(f">{seq_id}\n{row['context_sequence']}\n")

    # Print summary
    click.echo("\n" + "=" * 60)
    click.echo("Arrest Motif Extraction Summary")
    click.echo("=" * 60)
    click.echo(f"Total proteins: {total_proteins:,}")
    click.echo(
        f"Proteins with arrest motifs: {proteins_with_motifs:,} "
        f"({100*proteins_with_motifs/total_proteins:.4f}%)"
    )
    click.echo(f"Total arrest motifs found: {len(all_motifs):,}")

    click.echo("\nMotifs by type:")
    for motif in ARREST_MOTIFS:
        count = motif_counts[motif]
        at_cterm = df[(df["motif_type"] == motif) & df["is_at_cterm"]].shape[0]
        click.echo(f"  {motif}: {count:,} ({at_cterm} at C-terminus)")

    click.echo(f"\nOutput files:")
    click.echo(f"  Motifs: {motifs_out}")
    click.echo(f"  Sequences: {sequences_out}")


if __name__ == "__main__":
    main()
