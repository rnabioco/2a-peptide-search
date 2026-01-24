#!/usr/bin/env python3
"""
Refine arrest motif HMMs using hits from multiple databases.

Collects hits from bacterial and viral databases, filters by quality,
deduplicates, and builds refined HMMs with broader sequence diversity.
This iterative refinement improves sensitivity for detecting divergent
arrest peptides.
"""

import gzip
import re
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

import click
import pandas as pd


def parse_tblout(tblout_file):
    """Parse HMMER tblout format and return hits as list of dicts."""
    hits = []

    open_func = gzip.open if str(tblout_file).endswith(".gz") else open
    mode = "rt" if str(tblout_file).endswith(".gz") else "r"

    with open_func(tblout_file, mode) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 9:
                continue

            hits.append({
                "target_name": parts[0],
                "target_accession": parts[1],
                "query_name": parts[2],
                "query_accession": parts[3],
                "evalue": float(parts[4]),
                "score": float(parts[5]),
                "bias": float(parts[6]),
            })

    return hits


def parse_stockholm(sto_file):
    """Parse Stockholm alignment file and return sequences.

    Returns:
        dict: mapping of sequence ID to aligned sequence
    """
    sequences = {}

    open_func = gzip.open if str(sto_file).endswith(".gz") else open
    mode = "rt" if str(sto_file).endswith(".gz") else "r"

    with open_func(sto_file, mode) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line == "//":
                continue

            # Sequence line: "seqid aligned_sequence"
            parts = line.split(None, 1)
            if len(parts) == 2:
                seq_id, aligned_seq = parts
                if seq_id not in sequences:
                    sequences[seq_id] = aligned_seq
                else:
                    # Handle multi-line sequences (append)
                    sequences[seq_id] += aligned_seq

    return sequences


def extract_motif_type_from_query(query_name):
    """Extract motif type from query name (e.g., 'RAGP' from 'RAGP_arrest')."""
    match = re.match(r"([A-Z]{4})(?:_arrest)?", query_name)
    return match.group(1) if match else query_name


def write_stockholm(sequences, motif_type, output_path):
    """Write sequences in Stockholm format."""
    with open(output_path, "w") as f:
        f.write("# STOCKHOLM 1.0\n")
        f.write(f"#=GF ID {motif_type}_arrest_refined\n")
        f.write(f"#=GF DE Refined arrest motif {motif_type}\n")
        f.write(f"#=GF SQ {len(sequences)}\n")
        f.write("\n")

        for seq_id, seq in sequences.items():
            # Remove gaps for ungapped sequence
            ungapped = seq.replace("-", "").replace(".", "")
            f.write(f"{seq_id} {seq}\n")
            # Add accession annotation
            base_id = seq_id.split("/")[0]
            f.write(f"#=GS {seq_id} AC {base_id}\n")

        f.write("//\n")

    return len(sequences)


def deduplicate_sequences(sequences, identity_threshold=0.95):
    """Remove highly similar sequences to reduce redundancy.

    Simple deduplication based on ungapped sequence identity.
    """
    unique = {}
    seen_seqs = {}

    for seq_id, aligned_seq in sequences.items():
        ungapped = aligned_seq.replace("-", "").replace(".", "").upper()

        # Check if we've seen a very similar sequence
        is_duplicate = False
        for seen_seq, seen_id in seen_seqs.items():
            # Simple identity check (exact match after some normalization)
            if ungapped == seen_seq:
                is_duplicate = True
                break

        if not is_duplicate:
            unique[seq_id] = aligned_seq
            seen_seqs[ungapped] = seq_id

    return unique


def build_hmm(alignment_file, output_hmm, name):
    """Build HMM from Stockholm alignment using hmmbuild."""
    cmd = ["hmmbuild", "-n", name, output_hmm, alignment_file]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        return True
    except subprocess.CalledProcessError as e:
        click.echo(f"Warning: hmmbuild failed for {name}: {e.stderr}")
        return False


@click.command()
@click.option(
    "--bacteria-alignments",
    required=True,
    help="Directory containing bacterial arrest motif alignments",
)
@click.option(
    "--cross-alignments",
    multiple=True,
    help="Cross-search alignment files (Stockholm, gzipped)",
)
@click.option(
    "--cross-tblouts",
    multiple=True,
    help="Cross-search tblout files",
)
@click.option(
    "--output-dir",
    required=True,
    help="Output directory for refined alignments",
)
@click.option(
    "--output-hmm",
    required=True,
    help="Output combined HMM database",
)
@click.option(
    "--summary",
    required=True,
    help="Output summary TSV",
)
@click.option(
    "--evalue-threshold",
    default=1e-3,
    help="E-value threshold for including cross-search hits",
)
@click.option(
    "--min-sequences",
    default=10,
    help="Minimum sequences per motif type for building HMM",
)
def main(
    bacteria_alignments,
    cross_alignments,
    cross_tblouts,
    output_dir,
    output_hmm,
    summary,
    evalue_threshold,
    min_sequences,
):
    """Refine arrest motif HMMs using hits from multiple databases."""

    click.echo("=" * 60)
    click.echo("Arrest Motif HMM Refinement")
    click.echo("=" * 60)

    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Collect sequences by motif type
    motif_sequences = defaultdict(dict)

    # 1. Load original bacterial alignments
    bacteria_path = Path(bacteria_alignments)
    click.echo(f"\nLoading bacterial alignments from {bacteria_path}")

    for sto_file in bacteria_path.glob("*.sto"):
        motif_type = sto_file.stem
        sequences = parse_stockholm(sto_file)
        click.echo(f"  {motif_type}: {len(sequences)} sequences")

        for seq_id, seq in sequences.items():
            # Add source prefix to avoid ID collisions
            prefixed_id = f"bact|{seq_id}"
            motif_sequences[motif_type][prefixed_id] = seq

    # 2. Parse cross-search tblout files to get hit quality information
    click.echo(f"\nParsing {len(cross_tblouts)} cross-search tblout files")
    hit_evalues = {}  # target_name -> best evalue

    for tblout_file in cross_tblouts:
        hits = parse_tblout(tblout_file)
        for hit in hits:
            target = hit["target_name"]
            evalue = hit["evalue"]
            if target not in hit_evalues or evalue < hit_evalues[target]:
                hit_evalues[target] = evalue

    click.echo(f"  Found {len(hit_evalues)} unique targets in cross-searches")

    # 3. Load cross-search alignments and filter by E-value
    click.echo(f"\nLoading {len(cross_alignments)} cross-search alignment files")

    for sto_file in cross_alignments:
        sequences = parse_stockholm(sto_file)
        click.echo(f"  {sto_file}: {len(sequences)} sequences")

        # Group sequences by motif type (from sequence coordinate comments or inference)
        for seq_id, seq in sequences.items():
            # Extract base ID (without coordinates)
            base_id = seq_id.split("/")[0]

            # Check if this hit passes E-value threshold
            if base_id in hit_evalues and hit_evalues[base_id] <= evalue_threshold:
                # Infer motif type from sequence - look for known motifs
                ungapped = seq.replace("-", "").replace(".", "").upper()

                # Find which motif type this sequence belongs to
                for motif in ["RAGP", "RAPG", "QAPP", "QGPP", "HAPP", "HGPP", "RAPP", "RPPP"]:
                    if motif in ungapped:
                        prefixed_id = f"cross|{seq_id}"
                        motif_sequences[motif][prefixed_id] = seq
                        break

    # 4. Deduplicate and create refined alignments
    click.echo("\nDeduplicating and building refined alignments")

    summary_data = []
    hmm_files = []

    for motif_type in sorted(motif_sequences.keys()):
        sequences = motif_sequences[motif_type]
        n_original = len(sequences)

        # Deduplicate
        unique_sequences = deduplicate_sequences(sequences)
        n_unique = len(unique_sequences)

        click.echo(f"  {motif_type}: {n_original} -> {n_unique} unique sequences")

        if n_unique < min_sequences:
            click.echo(f"    Skipping (< {min_sequences} sequences)")
            summary_data.append({
                "motif_type": motif_type,
                "n_original": n_original,
                "n_unique": n_unique,
                "n_bacteria": sum(1 for k in sequences if k.startswith("bact|")),
                "n_cross": sum(1 for k in sequences if k.startswith("cross|")),
                "hmm_built": False,
                "reason": f"too few sequences (< {min_sequences})",
            })
            continue

        # Write refined alignment
        sto_path = output_path / f"{motif_type}.refined.sto"
        write_stockholm(unique_sequences, motif_type, sto_path)

        # Build HMM
        hmm_path = output_path / f"{motif_type}.hmm"
        success = build_hmm(str(sto_path), str(hmm_path), f"{motif_type}_refined")

        if success:
            hmm_files.append(hmm_path)

        summary_data.append({
            "motif_type": motif_type,
            "n_original": n_original,
            "n_unique": n_unique,
            "n_bacteria": sum(1 for k in sequences if k.startswith("bact|")),
            "n_cross": sum(1 for k in sequences if k.startswith("cross|")),
            "hmm_built": success,
            "reason": "" if success else "hmmbuild failed",
        })

    # 5. Concatenate all HMMs into single database
    click.echo(f"\nCombining {len(hmm_files)} HMMs into database")

    if hmm_files:
        with open(output_hmm, "w") as out:
            for hmm_file in hmm_files:
                with open(hmm_file) as f:
                    out.write(f.read())

        # Press the database
        subprocess.run(["hmmpress", "-f", output_hmm], capture_output=True)
        click.echo(f"  Created {output_hmm}")
    else:
        # Create empty HMM file
        Path(output_hmm).touch()
        click.echo("  Warning: No HMMs built")

    # 6. Save summary
    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv(summary, sep="\t", index=False)

    # Print final summary
    click.echo("\n" + "=" * 60)
    click.echo("Refinement Summary")
    click.echo("=" * 60)

    successful = summary_df[summary_df["hmm_built"]]
    click.echo(f"HMMs built: {len(successful)}/{len(summary_df)}")

    if not successful.empty:
        click.echo("\nMotif types with refined HMMs:")
        for _, row in successful.iterrows():
            click.echo(
                f"  {row['motif_type']}: {row['n_unique']} sequences "
                f"({row['n_bacteria']} bacteria, {row['n_cross']} cross-database)"
            )

    failed = summary_df[~summary_df["hmm_built"]]
    if not failed.empty:
        click.echo(f"\nSkipped ({len(failed)}):")
        for _, row in failed.iterrows():
            click.echo(f"  {row['motif_type']}: {row['reason']}")


if __name__ == "__main__":
    main()
