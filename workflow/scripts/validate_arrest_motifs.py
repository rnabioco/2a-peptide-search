#!/usr/bin/env python3
"""
Validate arrest motif HMMs against known bacterial stalling peptides.

Runs hmmscan on known SecM, TnaC, and other characterized stalling peptides
against HMMs built from extracted arrest motifs to measure recovery rate.
"""

import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

import click
import pandas as pd
from Bio import SeqIO


def run_hmmscan(hmm_db, fasta_file, evalue=10.0):
    """Run hmmscan and parse results."""
    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".tblout") as tmp:
        tmp_tblout = tmp.name

    cmd = [
        "hmmscan",
        "--tblout", tmp_tblout,
        "-E", str(evalue),
        "--noali",
        hmm_db,
        fasta_file,
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)

    hits = []
    if result.returncode == 0:
        with open(tmp_tblout) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.split()
                if len(fields) >= 9:
                    hits.append({
                        "hmm_name": fields[0],  # target (HMM name = motif type)
                        "peptide_id": fields[2],  # query (sequence name)
                        "evalue": float(fields[4]),
                        "score": float(fields[5]),
                        "bias": float(fields[6]),
                    })

    Path(tmp_tblout).unlink(missing_ok=True)
    return hits


@click.command()
@click.option("--hmm-db", required=True, help="HMM database file (arrest_motifs.hmm)")
@click.option("--known-peptides", required=True, help="FASTA of known stalling peptides")
@click.option("--output", required=True, help="Output validation report TSV")
@click.option("--summary", required=True, help="Output summary file")
@click.option("--evalue", default=10.0, help="E-value threshold for hmmscan")
def main(hmm_db, known_peptides, output, summary, evalue):
    """Validate arrest motif HMMs against known stalling peptides."""

    click.echo("=" * 60)
    click.echo("Validating Arrest Motif HMMs Against Known Stalling Peptides")
    click.echo("=" * 60)

    # Load known peptides info
    click.echo(f"\nLoading known peptides from {known_peptides}...")
    known = {}
    for record in SeqIO.parse(known_peptides, "fasta"):
        parts = record.id.split("|")
        motif = parts[0] if parts else "Unknown"
        organism = parts[1] if len(parts) > 1 else "Unknown"
        gene = parts[2] if len(parts) > 2 else "Unknown"
        known[record.id] = {
            "motif_type": motif,
            "organism": organism,
            "gene": gene,
            "sequence": str(record.seq),
        }
    click.echo(f"Loaded {len(known)} known stalling peptides")

    # Check HMM database exists
    if not Path(hmm_db).exists():
        click.echo(f"ERROR: HMM database not found: {hmm_db}")
        with open(summary, "w") as f:
            f.write("ERROR: HMM database not found\n")
        pd.DataFrame().to_csv(output, sep="\t", index=False)
        return

    # Count HMMs in database
    n_hmms = 0
    with open(hmm_db) as f:
        for line in f:
            if line.startswith("NAME"):
                n_hmms += 1
    click.echo(f"HMM database contains {n_hmms} models")

    # Run hmmscan
    click.echo(f"\nSearching known peptides against arrest motif HMMs...")
    hits = run_hmmscan(hmm_db, known_peptides, evalue)
    click.echo(f"Found {len(hits)} hits")

    # Process results
    summary_lines = []
    summary_lines.append("=" * 60)
    summary_lines.append("Arrest Motif Validation Summary")
    summary_lines.append("=" * 60)
    summary_lines.append("")

    if hits:
        hits_df = pd.DataFrame(hits)

        # Add peptide info
        hits_df["known_motif_type"] = hits_df["peptide_id"].map(
            lambda x: known.get(x, {}).get("motif_type", "Unknown")
        )
        hits_df["organism"] = hits_df["peptide_id"].map(
            lambda x: known.get(x, {}).get("organism", "Unknown")
        )
        hits_df["gene"] = hits_df["peptide_id"].map(
            lambda x: known.get(x, {}).get("gene", "Unknown")
        )

        # Check if HMM motif matches known motif (expected match)
        hits_df["motif_match"] = hits_df["hmm_name"].str.replace("_arrest", "") == hits_df["known_motif_type"]

        # Sort by peptide, then evalue
        hits_df = hits_df.sort_values(["peptide_id", "evalue"])

        # Save all hits
        hits_df.to_csv(output, sep="\t", index=False)

        # Best hit per peptide
        best_hits = hits_df.loc[hits_df.groupby("peptide_id")["evalue"].idxmin()]

        # Generate summary
        n_matched = best_hits["peptide_id"].nunique()
        n_total = len(known)
        n_correct_motif = best_hits["motif_match"].sum()

        summary_lines.append(f"Total known peptides: {n_total}")
        summary_lines.append(f"Matched by HMMs: {n_matched} ({100*n_matched/n_total:.1f}%)")
        summary_lines.append(f"Correct motif match: {n_correct_motif} ({100*n_correct_motif/n_matched:.1f}% of matched)")
        summary_lines.append(f"Total hits: {len(hits_df)}")
        summary_lines.append("")

        # Recovery by known motif type
        summary_lines.append("Recovery by known motif type:")
        known_motif_types = sorted(set(known[p]["motif_type"] for p in known))
        for motif in known_motif_types:
            total = sum(1 for p in known if known[p]["motif_type"] == motif)
            found = (best_hits["known_motif_type"] == motif).sum()
            correct = ((best_hits["known_motif_type"] == motif) & best_hits["motif_match"]).sum()
            summary_lines.append(f"  {motif:6s}: {found:2d}/{total:2d} recovered ({100*found/total:.0f}%), {correct} correct motif")
        summary_lines.append("")

        # HMM hits summary
        summary_lines.append("Hits by HMM model:")
        hmm_summary = best_hits.groupby("hmm_name").agg({
            "peptide_id": "count",
            "known_motif_type": lambda x: ", ".join(sorted(set(x))),
            "motif_match": "sum",
        }).reset_index()
        hmm_summary.columns = ["hmm_name", "n_hits", "known_types", "n_correct"]
        hmm_summary = hmm_summary.sort_values("n_hits", ascending=False)

        for _, row in hmm_summary.iterrows():
            summary_lines.append(
                f"  {row['hmm_name']:12s}: {row['n_hits']:2d} hits ({row['n_correct']} correct), "
                f"from known types: {row['known_types']}"
            )
        summary_lines.append("")

        # Top 10 best hits
        summary_lines.append("Top 10 best hits (by E-value):")
        for _, row in best_hits.nsmallest(10, "evalue").iterrows():
            match_str = "✓" if row["motif_match"] else "✗"
            summary_lines.append(
                f"  {row['gene']:12s} ({row['known_motif_type']}) -> {row['hmm_name']} "
                f"E={row['evalue']:.2e} {match_str}"
            )
        summary_lines.append("")

        # Not found
        found_peptides = set(best_hits["peptide_id"])
        not_found = [p for p in known if p not in found_peptides]
        if not_found:
            summary_lines.append(f"Peptides NOT found ({len(not_found)}):")
            for p in sorted(not_found)[:20]:
                summary_lines.append(f"  {known[p]['gene']} ({known[p]['motif_type']})")
            if len(not_found) > 20:
                summary_lines.append(f"  ... and {len(not_found) - 20} more")

    else:
        summary_lines.append("No hits found!")
        summary_lines.append("")
        summary_lines.append("This may indicate:")
        summary_lines.append("  - No arrest motifs were extracted from the database")
        summary_lines.append("  - The E-value threshold is too stringent")
        summary_lines.append("  - The HMM models don't capture the known peptide patterns")
        pd.DataFrame().to_csv(output, sep="\t", index=False)

    summary_lines.append("")
    summary_lines.append("=" * 60)

    summary_text = "\n".join(summary_lines)
    with open(summary, "w") as f:
        f.write(summary_text)
    click.echo("\n" + summary_text)


if __name__ == "__main__":
    main()
