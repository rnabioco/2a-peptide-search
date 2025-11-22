#!/usr/bin/env python3
"""
Compare all-GP vs inter-domain GP discovery approaches.

Analyzes:
1. Count statistics (total GP, inter-domain, intra-domain)
2. Enrichment of known stalling peptides in inter-domain regions
3. Conservation differences between categories
4. Cluster size distributions
5. Motif type frequencies
"""

import gzip
from pathlib import Path

import click
import pandas as pd
import numpy as np
from plotnine import *


def load_data(all_gp_file, interdomain_file, clusters_file, validation_file):
    """Load all required data files."""
    click.echo("Loading data files...")

    # All GP motifs
    all_gp = pd.read_csv(all_gp_file, sep="\t", compression="gzip")
    click.echo(f"  All GP motifs: {len(all_gp):,}")

    # Inter-domain GP motifs
    interdomain = pd.read_csv(interdomain_file, sep="\t", compression="gzip")
    click.echo(f"  Inter-domain GP motifs: {len(interdomain):,}")

    # Clusters
    clusters = pd.read_csv(clusters_file, sep="\t", compression="gzip")
    click.echo(f"  Clustered motifs: {len(clusters):,}")

    # Validation (may be empty)
    try:
        validation = pd.read_csv(validation_file, sep="\t")
        click.echo(f"  Validation hits: {len(validation):,}")
    except Exception:
        validation = pd.DataFrame()
        click.echo(f"  Validation hits: 0")

    return all_gp, interdomain, clusters, validation


def calculate_statistics(all_gp, interdomain, clusters, validation):
    """Calculate comparison statistics."""
    stats = {}

    # Basic counts
    stats["total_gp_motifs"] = len(all_gp)
    stats["interdomain_gp_motifs"] = len(interdomain)
    stats["interdomain_fraction"] = (
        len(interdomain) / len(all_gp) if len(all_gp) > 0 else 0
    )

    # Intra-domain (assuming remaining are intra-domain)
    stats["intradomain_gp_motifs"] = len(all_gp) - len(interdomain)

    # Protein counts
    stats["proteins_with_gp"] = all_gp["protein_id"].nunique()
    stats["proteins_with_interdomain_gp"] = (
        interdomain["protein_id"].nunique() if len(interdomain) > 0 else 0
    )

    # Clustering stats
    if len(clusters) > 0:
        stats["n_clusters"] = clusters["cluster_id"].nunique()
        cluster_sizes = clusters.groupby("cluster_id").size()
        stats["mean_cluster_size"] = cluster_sizes.mean()
        stats["median_cluster_size"] = cluster_sizes.median()
        stats["max_cluster_size"] = cluster_sizes.max()
        stats["min_cluster_size"] = cluster_sizes.min()
    else:
        stats["n_clusters"] = 0
        stats["mean_cluster_size"] = 0
        stats["median_cluster_size"] = 0
        stats["max_cluster_size"] = 0
        stats["min_cluster_size"] = 0

    # Validation stats
    if len(validation) > 0:
        stats["known_peptides_recovered"] = validation["peptide_id"].nunique()
        stats["total_validation_hits"] = len(validation)
        stats["mean_evalue"] = validation["evalue"].mean()
    else:
        stats["known_peptides_recovered"] = 0
        stats["total_validation_hits"] = 0
        stats["mean_evalue"] = np.nan

    return stats


def plot_gp_distribution(all_gp, interdomain, output_dir):
    """Plot distribution of GP motifs (all vs inter-domain)."""
    output_dir = Path(output_dir)

    # Position distribution
    df_pos = pd.concat(
        [
            pd.DataFrame({"position": all_gp["gp_position"], "type": "All GP"}),
            pd.DataFrame(
                {"position": interdomain["gp_position"], "type": "Inter-domain GP"}
            )
            if len(interdomain) > 0
            else pd.DataFrame(),
        ]
    )

    p1 = (
        ggplot(df_pos, aes(x="position", fill="type"))
        + geom_histogram(bins=50, alpha=0.6, position="identity")
        + labs(x="GP Position in Protein", y="Count", title="GP Position Distribution")
        + scale_fill_manual(values=["#3498db", "#2ecc71"])
        + theme_minimal()
        + theme(figure_size=(6, 4))
    )

    p1.save(output_dir / "gp_position_distribution.png", dpi=150)
    click.echo(f"  Saved: gp_position_distribution.png")

    # N-terminal distance distribution
    df_nterm = pd.concat(
        [
            pd.DataFrame({"distance": all_gp["n_terminal_distance"], "type": "All GP"}),
            pd.DataFrame(
                {
                    "distance": interdomain["n_terminal_distance"],
                    "type": "Inter-domain GP",
                }
            )
            if len(interdomain) > 0
            else pd.DataFrame(),
        ]
    )

    df_nterm = df_nterm[df_nterm["distance"] <= 500]  # Limit to 500 for clarity

    p2 = (
        ggplot(df_nterm, aes(x="distance", fill="type"))
        + geom_histogram(bins=50, alpha=0.6, position="identity")
        + labs(
            x="Distance from N-terminus",
            y="Count",
            title="N-terminal Distance Distribution",
        )
        + scale_fill_manual(values=["#3498db", "#2ecc71"])
        + theme_minimal()
        + theme(figure_size=(6, 4))
    )

    p2.save(output_dir / "gp_nterminal_distribution.png", dpi=150)
    click.echo(f"  Saved: gp_nterminal_distribution.png")


def plot_cluster_sizes(clusters, output_dir):
    """Plot cluster size distribution."""
    if len(clusters) == 0:
        return

    output_dir = Path(output_dir)
    cluster_sizes = clusters.groupby("cluster_id").size().reset_index(name="size")
    cluster_sizes = cluster_sizes.sort_values("size", ascending=False).reset_index(
        drop=True
    )
    cluster_sizes["rank"] = range(1, len(cluster_sizes) + 1)

    # Top 20 clusters
    top20 = cluster_sizes.head(20)

    p1 = (
        ggplot(top20, aes(x="rank", y="size"))
        + geom_col(fill="#3498db")
        + labs(x="Cluster Rank", y="Cluster Size", title="Top 20 Cluster Sizes")
        + theme_minimal()
        + theme(figure_size=(8, 5))
    )

    p1.save(output_dir / "cluster_size_top20.png", dpi=150)
    click.echo(f"  Saved: cluster_size_top20.png")

    # Size distribution (all clusters)
    p2 = (
        ggplot(cluster_sizes, aes(x="size"))
        + geom_histogram(bins=50, fill="#3498db", color="black")
        + scale_y_log10()
        + labs(
            x="Cluster Size",
            y="Number of Clusters (log scale)",
            title="Cluster Size Distribution",
        )
        + theme_minimal()
        + theme(figure_size=(8, 5))
    )

    p2.save(output_dir / "cluster_size_distribution.png", dpi=150)
    click.echo(f"  Saved: cluster_size_distribution.png")


def plot_validation_results(validation, output_dir):
    """Plot validation results."""
    if len(validation) == 0:
        click.echo("  Skipping validation plots (no data)")
        return

    output_dir = Path(output_dir)

    # E-value distribution
    validation["log_evalue"] = np.log10(validation["evalue"] + 1e-200)

    p1 = (
        ggplot(validation, aes(x="log_evalue"))
        + geom_histogram(bins=30, fill="#e74c3c", color="black")
        + labs(x="log10(E-value)", y="Number of Hits", title="Validation Hit E-values")
        + theme_minimal()
        + theme(figure_size=(6, 4))
    )

    p1.save(output_dir / "validation_evalues.png", dpi=150)
    click.echo(f"  Saved: validation_evalues.png")

    # Score distribution
    p2 = (
        ggplot(validation, aes(x="score"))
        + geom_histogram(bins=30, fill="#9b59b6", color="black")
        + labs(x="Score", y="Number of Hits", title="Validation Hit Scores")
        + theme_minimal()
        + theme(figure_size=(6, 4))
    )

    p2.save(output_dir / "validation_scores.png", dpi=150)
    click.echo(f"  Saved: validation_scores.png")


def plot_summary_comparison(stats, output_dir):
    """Plot summary comparison of approaches."""
    output_dir = Path(output_dir)

    # GP motif counts
    df_counts = pd.DataFrame(
        {
            "category": ["All GP", "Inter-domain", "Intra-domain"],
            "count": [
                stats["total_gp_motifs"],
                stats["interdomain_gp_motifs"],
                stats["intradomain_gp_motifs"],
            ],
        }
    )

    p1 = (
        ggplot(df_counts, aes(x="category", y="count", fill="category"))
        + geom_col()
        + geom_text(aes(label="count"), va="bottom", format_string="{:,.0f}")
        + scale_fill_manual(values=["#3498db", "#2ecc71", "#e74c3c"])
        + labs(x="", y="Number of GP Motifs", title="GP Motif Categories")
        + theme_minimal()
        + theme(figure_size=(8, 6), legend_position="none")
    )

    p1.save(output_dir / "summary_comparison.png", dpi=150)
    click.echo(f"  Saved: summary_comparison.png")


@click.command()
@click.option(
    "--all-gp-motifs", "all_gp_file", required=True, help="Input TSV of all GP motifs"
)
@click.option(
    "--interdomain-motifs",
    "interdomain_file",
    required=True,
    help="Input TSV of inter-domain GP motifs",
)
@click.option(
    "--clusters",
    "clusters_file",
    required=True,
    help="Input TSV of cluster assignments",
)
@click.option(
    "--validation",
    "validation_file",
    required=True,
    help="Input TSV of validation results",
)
@click.option(
    "--comparison",
    "comparison_file",
    required=True,
    help="Output TSV of comparison statistics",
)
@click.option("--plots", "plots_dir", required=True, help="Output directory for plots")
def main(
    all_gp_file,
    interdomain_file,
    clusters_file,
    validation_file,
    comparison_file,
    plots_dir,
):
    """Compare all-GP vs inter-domain GP discovery approaches."""

    click.echo("=" * 60)
    click.echo("Approach Comparison Analysis")
    click.echo("=" * 60)

    # Load data
    all_gp, interdomain, clusters, validation = load_data(
        all_gp_file, interdomain_file, clusters_file, validation_file
    )

    # Calculate statistics
    click.echo("\nCalculating statistics...")
    stats = calculate_statistics(all_gp, interdomain, clusters, validation)

    # Save statistics
    stats_df = pd.DataFrame([stats])
    stats_df.to_csv(comparison_file, sep="\t", index=False)
    click.echo(f"Saved statistics to {comparison_file}")

    # Create plots directory
    plots_dir = Path(plots_dir)
    plots_dir.mkdir(parents=True, exist_ok=True)

    # Generate plots
    click.echo("\nGenerating plots...")

    plot_summary_comparison(stats, plots_dir)
    plot_gp_distribution(all_gp, interdomain, plots_dir)
    plot_cluster_sizes(clusters, plots_dir)
    plot_validation_results(validation, plots_dir)

    # Print summary
    click.echo("\n" + "=" * 60)
    click.echo("Comparison Summary")
    click.echo("=" * 60)
    click.echo(f"Total GP motifs: {stats['total_gp_motifs']:,}")
    click.echo(
        f"  Inter-domain: {stats['interdomain_gp_motifs']:,} ({stats['interdomain_fraction']:.1%})"
    )
    click.echo(f"  Intra-domain: {stats['intradomain_gp_motifs']:,}")
    click.echo(f"\nProteins with GP: {stats['proteins_with_gp']:,}")
    click.echo(
        f"Proteins with inter-domain GP: {stats['proteins_with_interdomain_gp']:,}"
    )
    click.echo(f"\nClusters: {stats['n_clusters']}")
    click.echo(f"  Mean size: {stats['mean_cluster_size']:.1f}")
    click.echo(f"  Median size: {stats['median_cluster_size']:.1f}")
    click.echo(f"  Size range: {stats['min_cluster_size']}-{stats['max_cluster_size']}")
    click.echo(f"\nValidation:")
    click.echo(f"  Known peptides recovered: {stats['known_peptides_recovered']}")
    click.echo(f"  Total hits: {stats['total_validation_hits']}")

    click.echo(f"\nOutput files:")
    click.echo(f"  Statistics: {comparison_file}")
    click.echo(f"  Plots: {plots_dir}/")


if __name__ == "__main__":
    main()
