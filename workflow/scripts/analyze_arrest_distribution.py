#!/usr/bin/env python3
"""
Analyze the distribution of arrest motifs across databases.

Generates summary statistics and comparisons showing:
- Motif counts per database/organism
- Enriched/depleted motif types
- Host-phage comparison
"""

from pathlib import Path

import click
import numpy as np
import pandas as pd


def classify_database(db_name):
    """Classify database as host (bacteria/archaea) or phage."""
    host_dbs = {"bacteria", "archaea"}
    phage_dbs = {
        "uniprot_viruses",
        "ncbi_viral_refseq",
        "inphared_proteins",
        "imgvr",
        "gpd",
        "millardlab_proteins",
    }

    if db_name in host_dbs:
        return "host"
    elif db_name in phage_dbs:
        return "phage"
    else:
        return "other"


@click.command()
@click.option(
    "--summaries",
    multiple=True,
    required=True,
    help="Input motif summary TSV files (one per database)",
)
@click.option(
    "--cross-summary",
    required=True,
    help="Cross-search summary TSV",
)
@click.option(
    "--matrix-out",
    required=True,
    help="Output distribution matrix TSV",
)
@click.option(
    "--host-phage-out",
    required=True,
    help="Output host-phage comparison TSV",
)
@click.option(
    "--summary-out",
    required=True,
    help="Output summary text file",
)
def main(summaries, cross_summary, matrix_out, host_phage_out, summary_out):
    """Analyze arrest motif distribution across databases."""

    click.echo("=" * 60)
    click.echo("Arrest Motif Distribution Analysis")
    click.echo("=" * 60)

    # Load all database summaries
    all_data = []

    for summary_file in summaries:
        # Extract database name from path
        # Expected path: .../arrest_analysis/{database}/motif_summary.tsv
        path = Path(summary_file)
        db_name = path.parent.name

        try:
            df = pd.read_csv(summary_file, sep="\t")
            df["database"] = db_name
            df["db_type"] = classify_database(db_name)
            all_data.append(df)
            click.echo(f"  Loaded {db_name}: {len(df)} motif types")
        except Exception as e:
            click.echo(f"  Warning: Could not load {summary_file}: {e}")

    if not all_data:
        click.echo("Error: No data loaded!")
        # Create empty output files
        pd.DataFrame().to_csv(matrix_out, sep="\t", index=False)
        pd.DataFrame().to_csv(host_phage_out, sep="\t", index=False)
        Path(summary_out).write_text("No data available\n")
        return

    # Combine all data
    combined = pd.concat(all_data, ignore_index=True)

    click.echo(f"\nTotal records: {len(combined)}")
    click.echo(f"Databases: {combined['database'].nunique()}")
    click.echo(f"Motif types: {combined['motif_type'].nunique()}")

    # =========================================================================
    # 1. Create distribution matrix (motif types x databases)
    # =========================================================================

    click.echo("\nCreating distribution matrix...")

    # Pivot to create matrix
    matrix = combined.pivot_table(
        index="motif_type",
        columns="database",
        values="n_total",
        aggfunc="sum",
        fill_value=0,
    )

    # Add row totals
    matrix["total"] = matrix.sum(axis=1)

    # Sort by total count
    matrix = matrix.sort_values("total", ascending=False)

    # Save matrix
    matrix.to_csv(matrix_out, sep="\t")
    click.echo(f"  Saved to {matrix_out}")

    # =========================================================================
    # 2. Host vs Phage comparison
    # =========================================================================

    click.echo("\nCreating host-phage comparison...")

    # Aggregate by db_type
    host_phage = combined.groupby(["motif_type", "db_type"]).agg({
        "n_total": "sum",
        "n_aligned": "sum",
        "alignment_created": "sum",
    }).reset_index()

    # Pivot for comparison
    comparison = host_phage.pivot_table(
        index="motif_type",
        columns="db_type",
        values="n_total",
        fill_value=0,
    ).reset_index()

    # Calculate ratios and enrichment
    if "host" in comparison.columns and "phage" in comparison.columns:
        comparison["total"] = comparison["host"] + comparison["phage"]
        comparison["phage_fraction"] = comparison["phage"] / comparison["total"]

        # Log2 fold change (with pseudocount)
        comparison["log2_phage_vs_host"] = np.log2(
            (comparison["phage"] + 1) / (comparison["host"] + 1)
        )
    else:
        # Handle missing columns
        if "host" not in comparison.columns:
            comparison["host"] = 0
        if "phage" not in comparison.columns:
            comparison["phage"] = 0
        comparison["total"] = comparison.get("host", 0) + comparison.get("phage", 0)
        comparison["phage_fraction"] = 0
        comparison["log2_phage_vs_host"] = 0

    # Sort by total
    comparison = comparison.sort_values("total", ascending=False)

    # Save comparison
    comparison.to_csv(host_phage_out, sep="\t", index=False)
    click.echo(f"  Saved to {host_phage_out}")

    # =========================================================================
    # 3. Load cross-search summary
    # =========================================================================

    click.echo("\nLoading cross-search summary...")

    try:
        cross_df = pd.read_csv(cross_summary, sep="\t")
        click.echo(f"  Cross-searches: {len(cross_df)} pairs")
    except Exception as e:
        click.echo(f"  Warning: Could not load cross-summary: {e}")
        cross_df = pd.DataFrame()

    # =========================================================================
    # 4. Generate text summary
    # =========================================================================

    click.echo("\nGenerating summary report...")

    with open(summary_out, "w") as f:
        f.write("=" * 70 + "\n")
        f.write("ARREST MOTIF DISTRIBUTION ANALYSIS\n")
        f.write("=" * 70 + "\n\n")

        # Overall statistics
        f.write("OVERALL STATISTICS\n")
        f.write("-" * 40 + "\n")
        f.write(f"Total databases analyzed: {combined['database'].nunique()}\n")
        f.write(f"Total motif types: {combined['motif_type'].nunique()}\n")
        f.write(f"Total motifs found: {combined['n_total'].sum():,}\n")
        f.write(f"Alignments created: {combined['n_aligned'].sum():,}\n\n")

        # Per-database summary
        f.write("PER-DATABASE SUMMARY\n")
        f.write("-" * 40 + "\n")

        db_summary = combined.groupby(["database", "db_type"]).agg({
            "n_total": "sum",
            "n_aligned": "sum",
            "motif_type": "count",
        }).rename(columns={"motif_type": "n_motif_types"})

        for (db, db_type), row in db_summary.iterrows():
            f.write(f"{db} ({db_type}):\n")
            f.write(f"  Total motifs: {row['n_total']:,}\n")
            f.write(f"  Aligned: {row['n_aligned']:,}\n")
            f.write(f"  Motif types: {row['n_motif_types']}\n")

        f.write("\n")

        # Motif type distribution
        f.write("MOTIF TYPE DISTRIBUTION\n")
        f.write("-" * 40 + "\n")

        motif_summary = combined.groupby("motif_type")["n_total"].sum().sort_values(ascending=False)
        for motif, count in motif_summary.items():
            f.write(f"  {motif}: {count:,}\n")

        f.write("\n")

        # Host vs Phage
        f.write("HOST VS PHAGE COMPARISON\n")
        f.write("-" * 40 + "\n")

        if "host" in comparison.columns and "phage" in comparison.columns:
            total_host = comparison["host"].sum()
            total_phage = comparison["phage"].sum()
            f.write(f"Total host (bacteria/archaea): {total_host:,}\n")
            f.write(f"Total phage/viral: {total_phage:,}\n")
            if total_host > 0:
                f.write(f"Phage/Host ratio: {total_phage/total_host:.2f}\n")

            f.write("\nMotifs enriched in phage (log2 FC > 1):\n")
            enriched = comparison[comparison["log2_phage_vs_host"] > 1].sort_values(
                "log2_phage_vs_host", ascending=False
            )
            for _, row in enriched.iterrows():
                f.write(
                    f"  {row['motif_type']}: log2FC={row['log2_phage_vs_host']:.2f} "
                    f"(host={row['host']:,}, phage={row['phage']:,})\n"
                )

            f.write("\nMotifs enriched in host (log2 FC < -1):\n")
            depleted = comparison[comparison["log2_phage_vs_host"] < -1].sort_values(
                "log2_phage_vs_host"
            )
            for _, row in depleted.iterrows():
                f.write(
                    f"  {row['motif_type']}: log2FC={row['log2_phage_vs_host']:.2f} "
                    f"(host={row['host']:,}, phage={row['phage']:,})\n"
                )

        f.write("\n")

        # Cross-search summary
        if not cross_df.empty:
            f.write("CROSS-DATABASE SEARCH SUMMARY\n")
            f.write("-" * 40 + "\n")

            for _, row in cross_df.iterrows():
                f.write(
                    f"  {row['source_db']} -> {row['target_db']}: "
                    f"{row['n_hits']:,} hits ({row['n_unique_sequences']:,} unique)\n"
                )

    click.echo(f"  Saved to {summary_out}")

    # Print final summary to stdout
    click.echo("\n" + "=" * 60)
    click.echo("Analysis Complete")
    click.echo("=" * 60)
    click.echo(f"Distribution matrix: {matrix_out}")
    click.echo(f"Host-phage comparison: {host_phage_out}")
    click.echo(f"Summary report: {summary_out}")


if __name__ == "__main__":
    main()
