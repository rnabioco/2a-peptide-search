#!/usr/bin/env python3
"""
Analyze conservation patterns within GP motif clusters.

For each cluster:
1. Align sequences using MUSCLE
2. Calculate per-position conservation (Shannon entropy)
3. Identify conserved positions
4. Generate consensus sequences
5. Create sequence logos (placeholder for issue #12)
"""

import gzip
import subprocess
import tempfile
from collections import Counter
from pathlib import Path

import click
import pandas as pd
from Bio import AlignIO, SeqIO
from Bio.Align import AlignInfo
import numpy as np


def calculate_shannon_entropy(alignment_column):
    """Calculate Shannon entropy for an alignment column."""
    # Count amino acids
    counts = Counter(alignment_column.upper())

    # Remove gaps
    if "-" in counts:
        del counts["-"]

    if not counts:
        return 0.0

    total = sum(counts.values())

    # Calculate entropy
    entropy = 0.0
    for count in counts.values():
        p = count / total
        if p > 0:
            entropy -= p * np.log2(p)

    return entropy


def calculate_conservation_score(alignment_column):
    """
    Calculate conservation score (0-1, higher = more conserved).

    Uses normalized inverse Shannon entropy.
    """
    max_entropy = np.log2(20)  # Maximum entropy for 20 amino acids
    entropy = calculate_shannon_entropy(alignment_column)

    # Normalize: 1 = fully conserved, 0 = maximum diversity
    conservation = 1 - (entropy / max_entropy)

    return conservation


def align_sequences(sequences, method="muscle"):
    """Align sequences using MUSCLE."""
    # Write sequences to temp file
    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".fasta") as tmp_in:
        tmp_in_name = tmp_in.name
        for i, seq in enumerate(sequences):
            tmp_in.write(f">seq{i}\n{seq}\n")

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
        click.echo(f"Error running MUSCLE: {e}", err=True)
        return None

    # Read alignment
    try:
        alignment = AlignIO.read(tmp_out_name, "fasta")
    except Exception as e:
        click.echo(f"Error reading alignment: {e}", err=True)
        return None
    finally:
        # Cleanup
        import os

        os.unlink(tmp_in_name)
        os.unlink(tmp_out_name)

    return alignment


def analyze_cluster(cluster_id, sequences, min_size=5):
    """Analyze conservation for a single cluster."""
    if len(sequences) < min_size:
        return None

    # Align sequences
    alignment = align_sequences(sequences)

    if alignment is None or len(alignment) == 0:
        return None

    # Calculate per-position conservation
    alignment_length = alignment.get_alignment_length()
    conservation_scores = []
    entropy_scores = []
    consensus = []

    for i in range(alignment_length):
        column = alignment[:, i]

        # Calculate scores
        conservation = calculate_conservation_score(column)
        entropy = calculate_shannon_entropy(column)

        conservation_scores.append(conservation)
        entropy_scores.append(entropy)

        # Get consensus residue (most common, excluding gaps)
        counts = Counter(c.upper() for c in column if c != "-")
        if counts:
            consensus_aa = counts.most_common(1)[0][0]
            consensus.append(consensus_aa)
        else:
            consensus.append("-")

    # Find GP position in consensus
    consensus_seq = "".join(consensus)
    gp_pos = consensus_seq.find("GP")

    # Calculate mean conservation
    mean_conservation = np.mean(conservation_scores)

    # Identify highly conserved positions (>80% conservation)
    conserved_positions = [
        i for i, score in enumerate(conservation_scores) if score > 0.8
    ]

    return {
        "cluster_id": cluster_id,
        "cluster_size": len(sequences),
        "alignment_length": alignment_length,
        "mean_conservation": mean_conservation,
        "consensus_sequence": consensus_seq,
        "gp_position": gp_pos,
        "conserved_positions": ",".join(map(str, conserved_positions)),
        "conservation_scores": ",".join(f"{s:.3f}" for s in conservation_scores),
        "entropy_scores": ",".join(f"{s:.3f}" for s in entropy_scores),
    }


def create_placeholder_logo(cluster_id, output_dir):
    """Create placeholder PNG for sequence logo (to be implemented in issue #12)."""
    from PIL import Image, ImageDraw, ImageFont

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    output_file = output_dir / f"cluster_{cluster_id}.png"

    # Create simple placeholder image
    width, height = 800, 200
    img = Image.new("RGB", (width, height), color="white")
    draw = ImageDraw.Draw(img)

    # Add text
    text = f"Sequence Logo for Cluster {cluster_id}\n(To be generated in issue #12)"
    try:
        font = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 24)
    except:
        font = ImageFont.load_default()

    # Draw text in center
    bbox = draw.textbbox((0, 0), text, font=font)
    text_width = bbox[2] - bbox[0]
    text_height = bbox[3] - bbox[1]
    x = (width - text_width) / 2
    y = (height - text_height) / 2
    draw.text((x, y), text, fill="gray", font=font)

    # Save
    img.save(output_file)

    return str(output_file)


@click.command()
@click.option(
    "--clusters",
    "clusters_file",
    required=True,
    help="Input TSV of cluster assignments",
)
@click.option("--motifs", "motifs_file", required=True, help="Input TSV of GP motifs")
@click.option(
    "--conservation",
    "conservation_file",
    required=True,
    help="Output TSV of conservation data",
)
@click.option("--logos-dir", default=None, help="Output directory for sequence logos")
@click.option("--min-cluster-size", default=5, help="Minimum cluster size to analyze")
@click.option("--top-n", default=20, help="Number of top clusters to analyze")
def main(
    clusters_file, motifs_file, conservation_file, logos_dir, min_cluster_size, top_n
):
    """Analyze conservation patterns within GP motif clusters."""

    click.echo("=" * 60)
    click.echo("Cluster Conservation Analysis")
    click.echo("=" * 60)

    # Load data
    click.echo(f"Loading clusters from {clusters_file}...")
    clusters_df = pd.read_csv(clusters_file, sep="\t", compression="gzip")

    click.echo(f"Loading motifs from {motifs_file}...")
    motifs_df = pd.read_csv(motifs_file, sep="\t", compression="gzip")

    # Merge to get sequences for each cluster
    merged = pd.merge(
        clusters_df,
        motifs_df[["protein_id", "gp_index", "context_sequence"]],
        on=["protein_id", "gp_index"],
        how="left",
    )

    # Get cluster sizes
    cluster_sizes = merged.groupby("cluster_id").size().sort_values(ascending=False)
    click.echo(f"\nTotal clusters: {len(cluster_sizes)}")
    click.echo(f"Analyzing top {top_n} clusters (min size: {min_cluster_size})")

    # Analyze top N clusters
    conservation_results = []

    for i, (cluster_id, size) in enumerate(cluster_sizes.head(top_n).items(), 1):
        click.echo(f"\n[{i}/{top_n}] Analyzing cluster {cluster_id} (n={size})...")

        # Get sequences for this cluster
        cluster_seqs = merged[merged["cluster_id"] == cluster_id][
            "context_sequence"
        ].tolist()

        # Analyze conservation
        result = analyze_cluster(cluster_id, cluster_seqs, min_cluster_size)

        if result:
            conservation_results.append(result)
            click.echo(f"  Mean conservation: {result['mean_conservation']:.3f}")
            click.echo(f"  Consensus: {result['consensus_sequence']}")

            # Create placeholder logo
            if logos_dir:
                logo_file = create_placeholder_logo(cluster_id, logos_dir)
                click.echo(f"  Placeholder logo: {logo_file}")
        else:
            click.echo(f"  Skipped (insufficient data)")

    # Save conservation results
    if conservation_results:
        conservation_df = pd.DataFrame(conservation_results)
        conservation_df.to_csv(
            conservation_file, sep="\t", index=False, compression="gzip"
        )

        click.echo("\n" + "=" * 60)
        click.echo("Conservation Analysis Summary")
        click.echo("=" * 60)
        click.echo(f"Clusters analyzed: {len(conservation_results)}")
        click.echo(
            f"Mean conservation range: {conservation_df['mean_conservation'].min():.3f} - {conservation_df['mean_conservation'].max():.3f}"
        )

        # Show top conserved clusters
        click.echo("\nTop 5 most conserved clusters:")
        top_conserved = conservation_df.nlargest(5, "mean_conservation")
        for _, row in top_conserved.iterrows():
            click.echo(
                f"  Cluster {row['cluster_id']}: {row['mean_conservation']:.3f} (n={row['cluster_size']})"
            )
    else:
        click.echo("\nNo clusters analyzed!")


if __name__ == "__main__":
    main()
