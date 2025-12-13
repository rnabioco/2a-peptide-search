#!/usr/bin/env python3
"""
Validate discovered GP motif clusters against known stalling peptides.

Searches discovered cluster HMMs against a database of known stalling peptides
(SecM, TnaC, MifM, CydA, etc.) to assess recovery rate and motif classification.
"""

import gzip
import re
import subprocess
import tempfile
from pathlib import Path

import click
import pandas as pd
from Bio import SeqIO


def run_hmmsearch(hmm_file, fasta_file, evalue=10.0):
    """Run hmmsearch and parse results."""
    # Create temp file for tblout
    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".tblout") as tmp:
        tmp_tblout = tmp.name

    # Run hmmsearch
    cmd = [
        "hmmsearch",
        "--tblout",
        tmp_tblout,
        "-E",
        str(evalue),
        "--noali",
        hmm_file,
        fasta_file,
    ]

    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        click.echo(f"Error running hmmsearch: {e}", err=True)
        return []

    # Parse tblout
    hits = []
    with open(tmp_tblout, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue

            fields = line.split()
            if len(fields) < 9:
                continue

            target_name = fields[0]
            query_name = fields[2]
            full_evalue = float(fields[4])
            score = float(fields[5])

            hits.append(
                {
                    "target": target_name,
                    "query": query_name,
                    "evalue": full_evalue,
                    "score": score,
                }
            )

    # Cleanup
    import os

    os.unlink(tmp_tblout)

    return hits


def extract_peptide_info(fasta_file):
    """Extract information about known stalling peptides from FASTA."""
    peptides = {}

    for record in SeqIO.parse(fasta_file, "fasta"):
        # Try to extract motif type from description
        desc = record.description
        motif_type = "Unknown"

        # Look for common motif patterns in description
        if "RAGP" in desc or "SecM" in desc:
            motif_type = "RAGP"
        elif "TnaC" in desc:
            motif_type = "GP"
        elif "MifM" in desc:
            motif_type = "GP"
        elif "QAPP" in desc:
            motif_type = "QAPP"
        elif "RAPG" in desc:
            motif_type = "RAPG"

        # Extract organism if present
        organism = "Unknown"
        org_match = re.search(r"OS=([^=]+?)(?:\s+[A-Z]{2}=|$)", desc)
        if org_match:
            organism = org_match.group(1).strip()

        peptides[record.id] = {
            "peptide_id": record.id,
            "description": record.description,
            "sequence": str(record.seq),
            "length": len(record.seq),
            "motif_type": motif_type,
            "organism": organism,
        }

    return peptides


@click.command()
@click.option(
    "--hmms",
    "hmms_pattern",
    required=True,
    help="Pattern for HMM files (e.g., results/prokaryotic/models/initial/cluster_*.hmm)",
)
@click.option(
    "--known-peptides",
    "known_peptides_file",
    required=True,
    help="FASTA file of known stalling peptides",
)
@click.option(
    "--validation",
    "validation_file",
    required=True,
    help="Output TSV of validation hits",
)
@click.option("--summary", "summary_file", required=True, help="Output text summary")
@click.option(
    "--evalue",
    default=10.0,
    help="E-value threshold for hmmsearch (relaxed for validation)",
)
@click.option("--threads", default=4, help="Number of threads")
def main(
    hmms_pattern, known_peptides_file, validation_file, summary_file, evalue, threads
):
    """Validate discovered motifs against known stalling peptides."""

    click.echo("=" * 60)
    click.echo("Validation Against Known Stalling Peptides")
    click.echo("=" * 60)

    # Find all HMM files
    from glob import glob

    hmm_files = sorted(glob(hmms_pattern))

    click.echo(f"Found {len(hmm_files)} HMM files")

    if not hmm_files:
        click.echo(
            f"Error: No HMM files found matching pattern: {hmms_pattern}", err=True
        )
        return

    # Load known peptides info
    click.echo(f"Loading known peptides from {known_peptides_file}...")
    peptides_info = extract_peptide_info(known_peptides_file)
    click.echo(f"Known peptides: {len(peptides_info)}")

    # Search each HMM against known peptides
    all_hits = []

    for i, hmm_file in enumerate(hmm_files, 1):
        # Extract cluster ID from filename
        cluster_match = re.search(r"cluster_(\d+)", hmm_file)
        if cluster_match:
            cluster_id = int(cluster_match.group(1))
        else:
            cluster_id = i

        click.echo(f"\n[{i}/{len(hmm_files)}] Searching with cluster {cluster_id}...")

        # Run hmmsearch
        hits = run_hmmsearch(hmm_file, known_peptides_file, evalue)

        # Add cluster info to hits
        for hit in hits:
            hit["cluster_id"] = cluster_id
            hit["hmm_file"] = hmm_file

            # Add peptide info
            peptide_id = hit["target"]
            if peptide_id in peptides_info:
                hit.update(peptides_info[peptide_id])

            all_hits.append(hit)

        if hits:
            click.echo(f"  Found {len(hits)} hits")
        else:
            click.echo(f"  No hits")

    # Convert to DataFrame
    if all_hits:
        hits_df = pd.DataFrame(all_hits)

        # Sort by peptide, then by evalue
        hits_df = hits_df.sort_values(["peptide_id", "evalue"])

        # Save validation results
        hits_df.to_csv(validation_file, sep="\t", index=False)

        # Calculate statistics
        n_peptides_hit = hits_df["peptide_id"].nunique()
        n_total_peptides = len(peptides_info)
        recovery_rate = n_peptides_hit / n_total_peptides

        # Best hit for each peptide (lowest e-value)
        best_hits = hits_df.loc[hits_df.groupby("peptide_id")["evalue"].idxmin()]

        # Generate summary
        summary_lines = []
        summary_lines.append("=" * 60)
        summary_lines.append("Validation Summary")
        summary_lines.append("=" * 60)
        summary_lines.append(f"Total known stalling peptides: {n_total_peptides}")
        summary_lines.append(f"Peptides with hits: {n_peptides_hit}")
        summary_lines.append(f"Recovery rate: {recovery_rate:.1%}")
        summary_lines.append(f"Total hits: {len(hits_df)}")
        summary_lines.append(
            f"Average hits per peptide: {len(hits_df) / n_peptides_hit:.1f}"
        )
        summary_lines.append("")

        # Hits by motif type
        if "motif_type" in best_hits.columns:
            summary_lines.append("Recovery by motif type:")
            motif_recovery = best_hits["motif_type"].value_counts()
            for motif, count in motif_recovery.items():
                summary_lines.append(f"  {motif}: {count}")
            summary_lines.append("")

        # Top scoring hits
        summary_lines.append("Top 10 hits:")
        for _, row in best_hits.nsmallest(10, "evalue").iterrows():
            summary_lines.append(
                f"  {row['peptide_id']:20s} -> Cluster {row['cluster_id']:3d}  E={row['evalue']:.2e}  Score={row['score']:.1f}"
            )
        summary_lines.append("")

        # Peptides not hit
        peptides_not_hit = set(peptides_info.keys()) - set(
            hits_df["peptide_id"].unique()
        )
        if peptides_not_hit:
            summary_lines.append(f"Peptides not recovered ({len(peptides_not_hit)}):")
            for peptide_id in sorted(peptides_not_hit):
                info = peptides_info[peptide_id]
                summary_lines.append(
                    f"  {peptide_id:20s} ({info['motif_type']}) - {info['organism']}"
                )
        else:
            summary_lines.append("All known peptides recovered!")

        summary_lines.append("")
        summary_lines.append("=" * 60)

        # Write summary
        summary_text = "\n".join(summary_lines)
        with open(summary_file, "w") as f:
            f.write(summary_text)

        # Print summary
        click.echo("\n" + summary_text)

        click.echo(f"\nOutput files:")
        click.echo(f"  Validation hits: {validation_file}")
        click.echo(f"  Summary: {summary_file}")

    else:
        click.echo("\nNo hits found!")

        # Still write empty files
        pd.DataFrame().to_csv(validation_file, sep="\t", index=False)

        with open(summary_file, "w") as f:
            f.write("Validation Summary\n")
            f.write("=" * 60 + "\n")
            f.write(f"Total known stalling peptides: {len(peptides_info)}\n")
            f.write("Peptides with hits: 0\n")
            f.write("Recovery rate: 0.0%\n")
            f.write("\nNo hits found!\n")


if __name__ == "__main__":
    main()
