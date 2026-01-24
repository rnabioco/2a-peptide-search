#!/usr/bin/env python3
"""
Compare all-GP vs inter-domain GP vs arrest motif discovery approaches.

Analyzes:
1. Count statistics (total GP, inter-domain, intra-domain)
2. Enrichment of known stalling peptides in inter-domain regions
3. Conservation differences between categories
4. Cluster size distributions
5. Motif type frequencies
6. Arrest motif approach recovery and overlap
"""

import gzip
from pathlib import Path

import click
import pandas as pd
import numpy as np
from plotnine import *


def load_data(
    all_gp_file,
    interdomain_file,
    clusters_file,
    gp_validation_file,
    arrest_motifs_file,
    arrest_validation_file,
    arrest_summary_file,
):
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

    # GP validation (may be empty or not exist)
    try:
        gp_validation = pd.read_csv(gp_validation_file, sep="\t")
        click.echo(f"  GP validation hits: {len(gp_validation):,}")
    except Exception:
        gp_validation = pd.DataFrame()
        click.echo(f"  GP validation hits: 0")

    # Arrest motifs
    arrest_motifs = pd.read_csv(arrest_motifs_file, sep="\t", compression="gzip")
    click.echo(f"  Arrest motifs: {len(arrest_motifs):,}")

    # Arrest validation
    arrest_validation = pd.read_csv(arrest_validation_file, sep="\t")
    click.echo(f"  Arrest validation hits: {len(arrest_validation):,}")

    # Arrest summary
    arrest_summary = pd.read_csv(arrest_summary_file, sep="\t")
    click.echo(f"  Arrest motif types: {len(arrest_summary):,}")

    return (
        all_gp,
        interdomain,
        clusters,
        gp_validation,
        arrest_motifs,
        arrest_validation,
        arrest_summary,
    )


def calculate_statistics(
    all_gp, interdomain, clusters, gp_validation, arrest_validation, arrest_summary
):
    """Calculate comparison statistics for all approaches."""
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

    # GP validation stats
    if len(gp_validation) > 0:
        stats["gp_known_peptides_recovered"] = gp_validation["peptide_id"].nunique()
        stats["gp_total_validation_hits"] = len(gp_validation)
        stats["gp_mean_evalue"] = gp_validation["evalue"].mean()
    else:
        stats["gp_known_peptides_recovered"] = 0
        stats["gp_total_validation_hits"] = 0
        stats["gp_mean_evalue"] = np.nan

    # Arrest motif stats
    stats["arrest_motif_types"] = len(arrest_summary)
    stats["arrest_total_motifs"] = arrest_summary["n_total"].sum()
    stats["arrest_aligned_motifs"] = arrest_summary["n_aligned"].sum()

    # Arrest validation stats
    if len(arrest_validation) > 0:
        # Extract unique peptide IDs (gene names) from validation
        stats["arrest_known_peptides_recovered"] = arrest_validation["gene"].nunique()
        stats["arrest_total_validation_hits"] = len(arrest_validation)
        stats["arrest_correct_motif_matches"] = arrest_validation["motif_match"].sum()
        stats["arrest_mean_evalue"] = arrest_validation["evalue"].mean()
    else:
        stats["arrest_known_peptides_recovered"] = 0
        stats["arrest_total_validation_hits"] = 0
        stats["arrest_correct_motif_matches"] = 0
        stats["arrest_mean_evalue"] = np.nan

    return stats


def calculate_approach_overlap(gp_validation, arrest_validation):
    """Calculate overlap between GP approach and arrest motif approach."""
    overlap_stats = {}

    # Get unique gene names from each approach
    if len(gp_validation) > 0 and "peptide_id" in gp_validation.columns:
        # Extract gene names from peptide_id format
        gp_genes = set()
        for pid in gp_validation["peptide_id"].unique():
            # Peptide ID format varies - extract gene if possible
            parts = str(pid).split("|")
            if len(parts) >= 3:
                gp_genes.add(parts[2])  # gene is typically 3rd field
            else:
                gp_genes.add(str(pid))
    else:
        gp_genes = set()

    if len(arrest_validation) > 0 and "gene" in arrest_validation.columns:
        arrest_genes = set(arrest_validation["gene"].unique())
    else:
        arrest_genes = set()

    # Calculate overlap
    overlap_stats["gp_unique_genes"] = len(gp_genes)
    overlap_stats["arrest_unique_genes"] = len(arrest_genes)
    overlap_stats["both_approaches"] = len(gp_genes & arrest_genes)
    overlap_stats["gp_only"] = len(gp_genes - arrest_genes)
    overlap_stats["arrest_only"] = len(arrest_genes - gp_genes)
    overlap_stats["union_total"] = len(gp_genes | arrest_genes)

    # Lists for detailed reporting
    overlap_stats["genes_in_both"] = sorted(gp_genes & arrest_genes)
    overlap_stats["genes_gp_only"] = sorted(gp_genes - arrest_genes)
    overlap_stats["genes_arrest_only"] = sorted(arrest_genes - gp_genes)

    return overlap_stats


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


def plot_approach_recovery_comparison(stats, overlap_stats, output_dir):
    """Plot comparison of recovery rates between GP and arrest motif approaches."""
    output_dir = Path(output_dir)

    # Recovery rate comparison
    # Assume total known peptides is 47 (from arrest analysis summary)
    total_known = 47

    gp_recovered = stats.get("gp_known_peptides_recovered", 0)
    arrest_recovered = stats.get("arrest_known_peptides_recovered", 0)

    df_recovery = pd.DataFrame(
        {
            "approach": ["GP Approach", "Arrest Motif"],
            "recovered": [gp_recovered, arrest_recovered],
            "rate": [
                gp_recovered / total_known * 100 if total_known > 0 else 0,
                arrest_recovered / total_known * 100 if total_known > 0 else 0,
            ],
        }
    )

    p1 = (
        ggplot(df_recovery, aes(x="approach", y="recovered", fill="approach"))
        + geom_col()
        + geom_text(aes(label="recovered"), va="bottom", size=12)
        + scale_fill_manual(values=["#3498db", "#e74c3c"])
        + labs(
            x="",
            y="Known Peptides Recovered",
            title=f"Recovery of Known Stalling Peptides (n={total_known})",
        )
        + theme_minimal()
        + theme(figure_size=(6, 5), legend_position="none")
    )

    p1.save(output_dir / "approach_recovery_comparison.png", dpi=150)
    click.echo(f"  Saved: approach_recovery_comparison.png")

    # Overlap Venn-style bar chart
    df_overlap = pd.DataFrame(
        {
            "category": ["GP Only", "Both Approaches", "Arrest Only"],
            "count": [
                overlap_stats["gp_only"],
                overlap_stats["both_approaches"],
                overlap_stats["arrest_only"],
            ],
        }
    )
    # Preserve order
    df_overlap["category"] = pd.Categorical(
        df_overlap["category"],
        categories=["GP Only", "Both Approaches", "Arrest Only"],
        ordered=True,
    )

    p2 = (
        ggplot(df_overlap, aes(x="category", y="count", fill="category"))
        + geom_col()
        + geom_text(aes(label="count"), va="bottom", size=12)
        + scale_fill_manual(values=["#3498db", "#9b59b6", "#e74c3c"])
        + labs(
            x="",
            y="Number of Known Peptides",
            title="Overlap Between Discovery Approaches",
        )
        + theme_minimal()
        + theme(figure_size=(7, 5), legend_position="none")
    )

    p2.save(output_dir / "approach_overlap.png", dpi=150)
    click.echo(f"  Saved: approach_overlap.png")

    # Recovery rate bar chart (percentage)
    p3 = (
        ggplot(df_recovery, aes(x="approach", y="rate", fill="approach"))
        + geom_col()
        + geom_text(aes(label="rate"), va="bottom", format_string="{:.1f}%", size=11)
        + scale_fill_manual(values=["#3498db", "#e74c3c"])
        + labs(x="", y="Recovery Rate (%)", title="Recovery Rate Comparison")
        + theme_minimal()
        + theme(figure_size=(6, 5), legend_position="none")
    )

    p3.save(output_dir / "approach_recovery_rates.png", dpi=150)
    click.echo(f"  Saved: approach_recovery_rates.png")


def plot_arrest_motif_summary(arrest_summary, output_dir):
    """Plot arrest motif type distribution."""
    output_dir = Path(output_dir)

    if len(arrest_summary) == 0:
        click.echo("  Skipping arrest motif plots (no data)")
        return

    # Sort by total count
    df = arrest_summary.copy()
    df = df.sort_values("n_total", ascending=True)
    df["motif_type"] = pd.Categorical(
        df["motif_type"], categories=df["motif_type"].tolist(), ordered=True
    )

    p1 = (
        ggplot(df, aes(x="motif_type", y="n_total"))
        + geom_col(fill="#e74c3c")
        + geom_col(aes(y="n_aligned"), fill="#3498db", alpha=0.7)
        + coord_flip()
        + labs(
            x="Arrest Motif Type",
            y="Count",
            title="Arrest Motif Counts (red=total, blue=aligned)",
        )
        + theme_minimal()
        + theme(figure_size=(6, 5))
    )

    p1.save(output_dir / "arrest_motif_types.png", dpi=150)
    click.echo(f"  Saved: arrest_motif_types.png")


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
    "--gp-validation",
    "gp_validation_file",
    required=True,
    help="Input TSV of GP approach validation results",
)
@click.option(
    "--arrest-motifs",
    "arrest_motifs_file",
    required=True,
    help="Input TSV of arrest motifs",
)
@click.option(
    "--arrest-validation",
    "arrest_validation_file",
    required=True,
    help="Input TSV of arrest motif validation results",
)
@click.option(
    "--arrest-summary",
    "arrest_summary_file",
    required=True,
    help="Input TSV of arrest motif summary",
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
    gp_validation_file,
    arrest_motifs_file,
    arrest_validation_file,
    arrest_summary_file,
    comparison_file,
    plots_dir,
):
    """Compare all-GP vs inter-domain GP vs arrest motif discovery approaches."""

    click.echo("=" * 60)
    click.echo("Approach Comparison Analysis")
    click.echo("=" * 60)

    # Load data
    (
        all_gp,
        interdomain,
        clusters,
        gp_validation,
        arrest_motifs,
        arrest_validation,
        arrest_summary,
    ) = load_data(
        all_gp_file,
        interdomain_file,
        clusters_file,
        gp_validation_file,
        arrest_motifs_file,
        arrest_validation_file,
        arrest_summary_file,
    )

    # Calculate statistics
    click.echo("\nCalculating statistics...")
    stats = calculate_statistics(
        all_gp, interdomain, clusters, gp_validation, arrest_validation, arrest_summary
    )

    # Calculate approach overlap
    click.echo("Calculating approach overlap...")
    overlap_stats = calculate_approach_overlap(gp_validation, arrest_validation)

    # Save statistics (combine stats and overlap)
    combined_stats = {**stats, **{k: v for k, v in overlap_stats.items() if not isinstance(v, list)}}
    stats_df = pd.DataFrame([combined_stats])
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
    plot_validation_results(gp_validation, plots_dir)
    plot_approach_recovery_comparison(stats, overlap_stats, plots_dir)
    plot_arrest_motif_summary(arrest_summary, plots_dir)

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

    click.echo(f"\n" + "-" * 40)
    click.echo("Recovery of Known Stalling Peptides")
    click.echo("-" * 40)
    click.echo(f"GP Approach:")
    click.echo(f"  Known peptides recovered: {stats['gp_known_peptides_recovered']}")
    click.echo(f"  Total validation hits: {stats['gp_total_validation_hits']}")
    click.echo(f"\nArrest Motif Approach:")
    click.echo(f"  Known peptides recovered: {stats['arrest_known_peptides_recovered']}")
    click.echo(f"  Total validation hits: {stats['arrest_total_validation_hits']}")
    click.echo(f"  Correct motif matches: {stats['arrest_correct_motif_matches']}")
    click.echo(f"  Motif types analyzed: {stats['arrest_motif_types']}")

    click.echo(f"\n" + "-" * 40)
    click.echo("Overlap Between Approaches")
    click.echo("-" * 40)
    click.echo(f"GP approach unique genes: {overlap_stats['gp_unique_genes']}")
    click.echo(f"Arrest approach unique genes: {overlap_stats['arrest_unique_genes']}")
    click.echo(f"Found by both approaches: {overlap_stats['both_approaches']}")
    click.echo(f"GP only: {overlap_stats['gp_only']}")
    click.echo(f"Arrest only: {overlap_stats['arrest_only']}")
    click.echo(f"Union (total unique): {overlap_stats['union_total']}")

    if overlap_stats['genes_in_both']:
        click.echo(f"\nGenes found by both: {', '.join(overlap_stats['genes_in_both'])}")
    if overlap_stats['genes_gp_only']:
        click.echo(f"GP only genes: {', '.join(overlap_stats['genes_gp_only'])}")
    if overlap_stats['genes_arrest_only']:
        click.echo(f"Arrest only genes: {', '.join(overlap_stats['genes_arrest_only'])}")

    click.echo(f"\nOutput files:")
    click.echo(f"  Statistics: {comparison_file}")
    click.echo(f"  Plots: {plots_dir}/")


if __name__ == "__main__":
    main()
