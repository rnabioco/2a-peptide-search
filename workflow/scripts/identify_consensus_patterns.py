#!/usr/bin/env python3
"""
Identify consensus patterns for major GP motif clusters.

For each cluster:
1. Extract consensus sequence from conservation analysis
2. Classify by motif type (RAGP, QAPP, RAPG, etc.)
3. Calculate quality metrics
4. Create Stockholm format alignments for HMM building
"""

import gzip
import re
import subprocess
import tempfile
from pathlib import Path

import click
import pandas as pd
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


MOTIF_PATTERNS = {
    "RAGP": r"RAGP",
    "RAPG": r"RAPG",
    "QAPP": r"QAPP",
    "QAGP": r"QAGP",
    "HGPP": r"HGPP",
    "HPGP": r"HPGP",
    "NPGP": r"NPGP",
    "RGP": r"[RK][AG]GP",  # R/K-A/G-GP
    "NGP": r"N[PG]GP",  # N-P/G-GP
    "GP": r"GP",  # Generic GP
}


def classify_motif_type(consensus_seq):
    """Classify consensus sequence by motif type."""
    # Try specific patterns first
    for motif_name, pattern in MOTIF_PATTERNS.items():
        if motif_name != "GP" and re.search(pattern, consensus_seq):
            return motif_name

    # Default to GP
    if "GP" in consensus_seq:
        return "GP"

    return "Unknown"


def calculate_pattern_quality(row, min_conservation=0.6):
    """Calculate quality score for a consensus pattern."""
    # Quality based on:
    # - Cluster size (larger = better)
    # - Mean conservation (higher = better)
    # - Presence of known motif patterns

    size_score = min(row["cluster_size"] / 100, 1.0)  # Normalize to 0-1
    conservation_score = row["mean_conservation"]

    # Bonus for known motif types
    motif_bonus = 0.2 if row["motif_type"] != "GP" else 0.0

    quality = size_score * 0.4 + conservation_score * 0.4 + motif_bonus * 0.2

    return quality


def align_cluster_sequences(cluster_id, sequences):
    """Align sequences for a cluster using MUSCLE."""
    # Write sequences to temp file
    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".fasta") as tmp_in:
        tmp_in_name = tmp_in.name
        for i, seq in enumerate(sequences):
            tmp_in.write(f">seq{i}|cluster{cluster_id}\n{seq}\n")

    # Create output file
    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".afa") as tmp_out:
        tmp_out_name = tmp_out.name

    # Run MUSCLE
    try:
        subprocess.run(
            ["muscle", "-align", tmp_in_name, "-output", tmp_out_name],
            check=True,
            capture_output=True,
        )
    except subprocess.CalledProcessError as e:
        click.echo(f"Error aligning cluster {cluster_id}: {e}", err=True)
        return None

    # Read alignment
    try:
        alignment = AlignIO.read(tmp_out_name, "fasta")
    except Exception as e:
        click.echo(f"Error reading alignment for cluster {cluster_id}: {e}", err=True)
        return None
    finally:
        # Cleanup
        import os

        os.unlink(tmp_in_name)
        os.unlink(tmp_out_name)

    return alignment


def save_stockholm_alignment(alignment, output_file, cluster_id):
    """Save alignment in Stockholm format."""
    # Convert to Stockholm format using esl-reformat
    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".afa") as tmp:
        tmp_name = tmp.name
        AlignIO.write(alignment, tmp, "fasta")

    # Convert to Stockholm
    try:
        with open(output_file, "w") as out:
            subprocess.run(
                ["esl-reformat", "stockholm", tmp_name], check=True, stdout=out
            )
    except subprocess.CalledProcessError as e:
        click.echo(f"Error converting to Stockholm: {e}", err=True)
        return False
    finally:
        import os

        os.unlink(tmp_name)

    return True


@click.command()
@click.option(
    "--conservation",
    "conservation_file",
    required=True,
    help="Input TSV of conservation data",
)
@click.option(
    "--clusters",
    "clusters_file",
    required=True,
    help="Input TSV of cluster assignments",
)
@click.option("--motifs", "motifs_file", required=True, help="Input TSV of GP motifs")
@click.option(
    "--consensus",
    "consensus_file",
    required=True,
    help="Output TSV of consensus patterns",
)
@click.option(
    "--alignments-dir", required=True, help="Output directory for Stockholm alignments"
)
@click.option("--min-conservation", default=0.6, help="Minimum conservation threshold")
@click.option("--min-cluster-size", default=5, help="Minimum cluster size")
@click.option("--top-n", default=20, help="Number of top clusters to process")
def main(
    conservation_file,
    clusters_file,
    motifs_file,
    consensus_file,
    alignments_dir,
    min_conservation,
    min_cluster_size,
    top_n,
):
    """Identify consensus patterns for major GP motif clusters."""

    click.echo("=" * 60)
    click.echo("Consensus Pattern Identification")
    click.echo("=" * 60)

    # Load conservation data
    click.echo(f"Loading conservation data from {conservation_file}...")
    conservation_df = pd.read_csv(conservation_file, sep="\t", compression="gzip")

    # Filter by conservation and size
    filtered = conservation_df[
        (conservation_df["mean_conservation"] >= min_conservation)
        & (conservation_df["cluster_size"] >= min_cluster_size)
    ].copy()

    click.echo(f"Clusters passing filters: {len(filtered)}/{len(conservation_df)}")

    # Classify motif types
    filtered["motif_type"] = filtered["consensus_sequence"].apply(classify_motif_type)

    # Calculate quality scores
    filtered["quality_score"] = filtered.apply(calculate_pattern_quality, axis=1)

    # Sort by quality
    filtered = filtered.sort_values("quality_score", ascending=False)

    # Take top N
    top_clusters = filtered.head(top_n)

    click.echo(f"\nProcessing top {len(top_clusters)} clusters:")
    click.echo(
        f"Quality range: {top_clusters['quality_score'].min():.3f} - {top_clusters['quality_score'].max():.3f}"
    )

    # Load cluster assignments and motifs
    click.echo(f"\nLoading clusters from {clusters_file}...")
    clusters_df = pd.read_csv(clusters_file, sep="\t", compression="gzip")

    click.echo(f"Loading motifs from {motifs_file}...")
    motifs_df = pd.read_csv(motifs_file, sep="\t", compression="gzip")

    # Merge to get sequences
    merged = pd.merge(
        clusters_df,
        motifs_df[["protein_id", "gp_index", "context_sequence"]],
        on=["protein_id", "gp_index"],
        how="left",
    )

    # Create output directory
    alignments_dir = Path(alignments_dir)
    alignments_dir.mkdir(parents=True, exist_ok=True)

    # Process each cluster
    consensus_patterns = []

    for idx, row in top_clusters.iterrows():
        cluster_id = row["cluster_id"]
        click.echo(
            f"\n[{idx + 1}/{len(top_clusters)}] Processing cluster {cluster_id}..."
        )
        click.echo(
            f"  Size: {row['cluster_size']}, Conservation: {row['mean_conservation']:.3f}"
        )
        click.echo(f"  Motif type: {row['motif_type']}")
        click.echo(f"  Consensus: {row['consensus_sequence']}")

        # Get sequences for this cluster
        cluster_seqs = merged[merged["cluster_id"] == cluster_id][
            "context_sequence"
        ].tolist()

        if not cluster_seqs:
            click.echo(f"  Warning: No sequences found for cluster {cluster_id}")
            continue

        # Align sequences
        alignment = align_cluster_sequences(cluster_id, cluster_seqs)

        if alignment is None:
            click.echo(f"  Warning: Failed to align cluster {cluster_id}")
            continue

        # Save Stockholm alignment
        output_file = alignments_dir / f"cluster_{cluster_id}.sto"
        success = save_stockholm_alignment(alignment, output_file, cluster_id)

        if success:
            click.echo(f"  Saved alignment: {output_file}")

            # Add to consensus patterns
            consensus_patterns.append(
                {
                    "cluster_id": cluster_id,
                    "cluster_size": row["cluster_size"],
                    "mean_conservation": row["mean_conservation"],
                    "consensus_sequence": row["consensus_sequence"],
                    "motif_type": row["motif_type"],
                    "quality_score": row["quality_score"],
                    "gp_position": row.get("gp_position", -1),
                    "alignment_file": str(output_file),
                }
            )
        else:
            click.echo(f"  Warning: Failed to save Stockholm format")

    # Save consensus patterns summary
    if consensus_patterns:
        consensus_df = pd.DataFrame(consensus_patterns)
        consensus_df.to_csv(consensus_file, sep="\t", index=False)

        click.echo("\n" + "=" * 60)
        click.echo("Consensus Patterns Summary")
        click.echo("=" * 60)
        click.echo(f"Total patterns identified: {len(consensus_patterns)}")

        # Group by motif type
        motif_counts = consensus_df["motif_type"].value_counts()
        click.echo("\nMotif type distribution:")
        for motif, count in motif_counts.items():
            click.echo(f"  {motif}: {count}")

        # Top patterns
        click.echo("\nTop 5 consensus patterns:")
        for i, row in consensus_df.head(5).iterrows():
            click.echo(
                f"  {row['cluster_id']:3d} ({row['motif_type']:6s}): {row['consensus_sequence'][:40]}"
            )
            click.echo(
                f"      Size: {row['cluster_size']:4d}, Conservation: {row['mean_conservation']:.3f}, Quality: {row['quality_score']:.3f}"
            )

        click.echo(f"\nOutput files:")
        click.echo(f"  Consensus patterns: {consensus_file}")
        click.echo(f"  Alignments: {alignments_dir}")
    else:
        click.echo("\nNo consensus patterns identified!")


if __name__ == "__main__":
    main()
