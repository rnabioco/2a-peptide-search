#!/usr/bin/env python3
"""
Validate GP motif clusters against known bacterial stalling peptides.

Uses hmmscan to search known peptides against HMMs built from cluster alignments.
"""

import subprocess
import tempfile
from pathlib import Path

import click
import pandas as pd
from Bio import SeqIO


def build_cluster_hmms(alignments_dir, hmm_db_path):
    """Build HMMs from cluster alignments and concatenate into a database."""
    alignments_dir = Path(alignments_dir)
    sto_files = sorted(alignments_dir.glob("cluster_*.sto"))

    if not sto_files:
        raise ValueError(f"No .sto files found in {alignments_dir}")

    click.echo(f"Building HMMs from {len(sto_files)} cluster alignments...")

    # Build individual HMMs and concatenate
    hmm_contents = []
    for sto_file in sto_files:
        cluster_id = sto_file.stem.replace("cluster_", "")

        with tempfile.NamedTemporaryFile(suffix=".hmm", delete=False) as tmp_hmm:
            tmp_hmm_path = tmp_hmm.name

        result = subprocess.run(
            ["hmmbuild", "-n", f"cluster_{cluster_id}", tmp_hmm_path, str(sto_file)],
            capture_output=True,
            text=True,
        )

        if result.returncode == 0:
            with open(tmp_hmm_path) as f:
                hmm_contents.append(f.read())

        Path(tmp_hmm_path).unlink(missing_ok=True)

    # Write concatenated HMM database
    with open(hmm_db_path, "w") as f:
        f.write("".join(hmm_contents))

    # Press the HMM database
    subprocess.run(["hmmpress", "-f", hmm_db_path], capture_output=True)

    click.echo(f"Built HMM database with {len(hmm_contents)} models")
    return len(hmm_contents)


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
                        "cluster": fields[0],  # target (HMM name)
                        "peptide_id": fields[2],  # query (sequence name)
                        "evalue": float(fields[4]),
                        "score": float(fields[5]),
                        "bias": float(fields[6]),
                    })

    Path(tmp_tblout).unlink(missing_ok=True)
    return hits


@click.command()
@click.option("--known-peptides", required=True, help="FASTA of known stalling peptides")
@click.option("--alignments-dir", required=True, help="Directory with cluster .sto alignments")
@click.option("--output", required=True, help="Output validation report TSV")
@click.option("--summary", required=True, help="Output summary file")
@click.option("--evalue", default=10.0, help="E-value threshold for hmmscan")
def main(known_peptides, alignments_dir, output, summary, evalue):
    """Validate clusters against known stalling peptides using hmmscan."""

    click.echo("=" * 60)
    click.echo("Validating Clusters Against Known Stalling Peptides")
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

    # Build HMM database from cluster alignments
    with tempfile.NamedTemporaryFile(suffix=".hmm", delete=False) as tmp_db:
        hmm_db_path = tmp_db.name

    try:
        n_hmms = build_cluster_hmms(alignments_dir, hmm_db_path)

        # Run hmmscan
        click.echo(f"\nSearching known peptides against {n_hmms} cluster HMMs...")
        hits = run_hmmscan(hmm_db_path, known_peptides, evalue)
        click.echo(f"Found {len(hits)} hits")

        # Process results
        if hits:
            hits_df = pd.DataFrame(hits)

            # Add peptide info
            hits_df["motif_type"] = hits_df["peptide_id"].map(lambda x: known.get(x, {}).get("motif_type", "Unknown"))
            hits_df["organism"] = hits_df["peptide_id"].map(lambda x: known.get(x, {}).get("organism", "Unknown"))
            hits_df["gene"] = hits_df["peptide_id"].map(lambda x: known.get(x, {}).get("gene", "Unknown"))

            # Extract cluster ID from HMM name
            hits_df["cluster_id"] = hits_df["cluster"].str.extract(r"cluster_(\d+)").astype(int)

            # Sort by peptide, then evalue
            hits_df = hits_df.sort_values(["peptide_id", "evalue"])

            # Save results
            hits_df.to_csv(output, sep="\t", index=False)

            # Best hit per peptide
            best_hits = hits_df.loc[hits_df.groupby("peptide_id")["evalue"].idxmin()]

            # Generate summary
            n_matched = best_hits["peptide_id"].nunique()
            n_total = len(known)

            summary_lines = []
            summary_lines.append("=" * 60)
            summary_lines.append("Known Stalling Peptide Validation Summary (hmmscan)")
            summary_lines.append("=" * 60)
            summary_lines.append("")
            summary_lines.append(f"Total known peptides: {n_total}")
            summary_lines.append(f"Matched to clusters: {n_matched} ({100*n_matched/n_total:.1f}%)")
            summary_lines.append(f"Total hits: {len(hits_df)}")
            summary_lines.append("")

            summary_lines.append("Clusters containing known peptides:")
            cluster_summary = best_hits.groupby("cluster_id").agg({
                "peptide_id": "count",
                "motif_type": lambda x: ", ".join(sorted(set(x))),
            }).reset_index()
            cluster_summary.columns = ["cluster_id", "n_known", "motif_types"]
            cluster_summary = cluster_summary.sort_values("n_known", ascending=False)

            for _, row in cluster_summary.iterrows():
                summary_lines.append(
                    f"  Cluster {int(row['cluster_id']):5d}: {row['n_known']:2d} known ({row['motif_types']})"
                )
            summary_lines.append("")

            summary_lines.append("Recovery by motif type:")
            for motif in sorted(set(known[p]["motif_type"] for p in known)):
                total = sum(1 for p in known if known[p]["motif_type"] == motif)
                found = (best_hits["motif_type"] == motif).sum()
                summary_lines.append(f"  {motif:6s}: {found:2d}/{total:2d} ({100*found/total:.0f}%)")
            summary_lines.append("")

            summary_lines.append("Top 10 hits:")
            for _, row in best_hits.nsmallest(10, "evalue").iterrows():
                summary_lines.append(
                    f"  {row['gene']:10s} ({row['motif_type']}) -> Cluster {row['cluster_id']:5d}  E={row['evalue']:.2e}"
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

            summary_lines.append("")
            summary_lines.append("=" * 60)

        else:
            summary_lines = ["No hits found!"]
            pd.DataFrame().to_csv(output, sep="\t", index=False)

        summary_text = "\n".join(summary_lines)
        with open(summary, "w") as f:
            f.write(summary_text)
        click.echo("\n" + summary_text)

    finally:
        # Cleanup HMM database files
        for ext in ["", ".h3f", ".h3i", ".h3m", ".h3p"]:
            Path(hmm_db_path + ext).unlink(missing_ok=True)


if __name__ == "__main__":
    main()
