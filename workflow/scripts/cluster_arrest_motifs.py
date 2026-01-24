#!/usr/bin/env python3
"""
Cluster arrest motif sequences by motif type and create alignments.

Groups extracted arrest motifs (RAGP, QAPP, etc.) by their motif type
and creates Stockholm alignments suitable for HMM building.
"""

import gzip
from collections import defaultdict
from pathlib import Path

import click
import pandas as pd


def create_stockholm_alignment(sequences, motif_type, output_path):
    """
    Create a Stockholm alignment from sequences.

    Sequences are aligned by their motif position - the motif itself
    serves as the anchor point for alignment.
    """
    if not sequences:
        return 0

    # Find the maximum upstream and downstream lengths
    max_upstream = max(s["motif_in_context"] for s in sequences)
    max_downstream = max(
        len(s["context_sequence"]) - s["motif_in_context"] - 4 for s in sequences
    )

    # Pad sequences to align by motif position
    aligned_seqs = []
    for i, s in enumerate(sequences):
        seq = s["context_sequence"]
        motif_pos = s["motif_in_context"]

        # Pad upstream
        upstream_pad = "-" * (max_upstream - motif_pos)
        # Pad downstream
        downstream_len = len(seq) - motif_pos - 4
        downstream_pad = "-" * (max_downstream - downstream_len)

        aligned_seq = upstream_pad + seq + downstream_pad

        # Create informative ID: protein_id/motif_position
        # Use the actual protein position in the original sequence
        protein_motif_pos = s.get("motif_position", i)
        seq_id = f"{s['protein_id']}/{protein_motif_pos}"
        # Truncate if too long for Stockholm format (max ~50 chars for name)
        if len(seq_id) > 50:
            seq_id = seq_id[:47] + "..."
        aligned_seqs.append((seq_id, aligned_seq, s["protein_id"]))

    # Write Stockholm format
    with open(output_path, "w") as f:
        f.write("# STOCKHOLM 1.0\n")
        f.write(f"#=GF ID {motif_type}_arrest\n")
        f.write(f"#=GF DE Bacterial arrest motif {motif_type}\n")
        f.write(f"#=GF SQ {len(aligned_seqs)}\n")

        for seq_id, aligned_seq, protein_id in aligned_seqs:
            f.write(f"{seq_id} {aligned_seq}\n")
            f.write(f"#=GS {seq_id} AC {protein_id}\n")

        f.write("//\n")

    return len(aligned_seqs)


@click.command()
@click.option("--motifs", required=True, help="Input TSV of arrest motifs (gzipped)")
@click.option("--alignments-dir", required=True, help="Output directory for alignments")
@click.option("--summary", required=True, help="Output summary TSV")
@click.option(
    "--min-sequences", default=10, help="Minimum sequences per motif type for alignment"
)
@click.option(
    "--max-sequences", default=1000, help="Maximum sequences per motif type (sample if more)"
)
def main(motifs, alignments_dir, summary, min_sequences, max_sequences):
    """Cluster arrest motifs by type and create alignments."""

    click.echo(f"Clustering arrest motifs from {motifs}")

    # Load motifs
    df = pd.read_csv(motifs, sep="\t", compression="gzip")

    if df.empty:
        click.echo("No arrest motifs found!")
        Path(alignments_dir).mkdir(parents=True, exist_ok=True)
        pd.DataFrame().to_csv(summary, sep="\t", index=False)
        return

    click.echo(f"Loaded {len(df):,} arrest motifs")

    # Group by motif type
    motif_groups = defaultdict(list)
    for _, row in df.iterrows():
        motif_groups[row["motif_type"]].append({
            "protein_id": row["protein_id"],
            "context_sequence": row["context_sequence"],
            "motif_in_context": row["motif_in_context"],
            "motif_position": row["motif_position"],
            "c_terminal_distance": row["c_terminal_distance"],
            "is_at_cterm": row["is_at_cterm"],
        })

    # Create output directory
    alignments_path = Path(alignments_dir)
    alignments_path.mkdir(parents=True, exist_ok=True)

    # Create alignments for each motif type
    summary_data = []

    for motif_type, sequences in sorted(motif_groups.items()):
        n_total = len(sequences)

        if n_total < min_sequences:
            click.echo(f"  {motif_type}: {n_total} sequences (skipped, < {min_sequences})")
            summary_data.append({
                "motif_type": motif_type,
                "n_total": n_total,
                "n_aligned": 0,
                "alignment_created": False,
                "reason": f"too few sequences (< {min_sequences})",
            })
            continue

        # Sample if too many sequences
        if n_total > max_sequences:
            # Prefer sequences at C-terminus
            at_cterm = [s for s in sequences if s["is_at_cterm"]]
            near_cterm = [s for s in sequences if not s["is_at_cterm"]]

            # Take all C-terminal ones first, then sample from near-C-terminal
            if len(at_cterm) >= max_sequences:
                import random
                random.seed(42)
                sequences = random.sample(at_cterm, max_sequences)
            else:
                n_needed = max_sequences - len(at_cterm)
                import random
                random.seed(42)
                sequences = at_cterm + random.sample(
                    near_cterm, min(n_needed, len(near_cterm))
                )

        # Deduplicate by sequence
        seen_seqs = set()
        unique_sequences = []
        for s in sequences:
            if s["context_sequence"] not in seen_seqs:
                seen_seqs.add(s["context_sequence"])
                unique_sequences.append(s)

        if len(unique_sequences) < min_sequences:
            click.echo(
                f"  {motif_type}: {n_total} total, {len(unique_sequences)} unique "
                f"(skipped, < {min_sequences} unique)"
            )
            summary_data.append({
                "motif_type": motif_type,
                "n_total": n_total,
                "n_aligned": 0,
                "alignment_created": False,
                "reason": f"too few unique sequences (< {min_sequences})",
            })
            continue

        # Create alignment
        sto_path = alignments_path / f"{motif_type}.sto"
        n_aligned = create_stockholm_alignment(unique_sequences, motif_type, sto_path)

        click.echo(f"  {motif_type}: {n_total} total -> {n_aligned} aligned")

        summary_data.append({
            "motif_type": motif_type,
            "n_total": n_total,
            "n_aligned": n_aligned,
            "alignment_created": True,
            "reason": "",
        })

    # Save summary
    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv(summary, sep="\t", index=False)

    # Print summary
    click.echo("\n" + "=" * 60)
    click.echo("Arrest Motif Clustering Summary")
    click.echo("=" * 60)

    created = summary_df[summary_df["alignment_created"]]
    click.echo(f"Alignments created: {len(created)}")

    for _, row in created.iterrows():
        click.echo(f"  {row['motif_type']}: {row['n_aligned']} sequences")

    skipped = summary_df[~summary_df["alignment_created"]]
    if not skipped.empty:
        click.echo(f"\nSkipped ({len(skipped)}):")
        for _, row in skipped.iterrows():
            click.echo(f"  {row['motif_type']}: {row['reason']}")


if __name__ == "__main__":
    main()
