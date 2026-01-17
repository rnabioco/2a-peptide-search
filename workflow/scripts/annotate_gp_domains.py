#!/usr/bin/env python3
"""
Annotate domains in GP-containing proteins using Pfam/HMMER.

Runs hmmscan to identify protein domains, enabling classification of GP motifs
as inter-domain, intra-domain, or in unstructured regions.
"""

import gzip
import re
import subprocess
import tempfile
from pathlib import Path

import click
import pandas as pd
from Bio import SeqIO


def run_hmmscan(fasta_file, pfam_db, evalue=1e-3, cpus=8):
    """Run hmmscan against Pfam database."""
    click.echo(f"Running hmmscan with {cpus} CPUs...")

    # Create temporary file for output
    with tempfile.NamedTemporaryFile(
        mode="w", delete=False, suffix=".domtblout"
    ) as tmp:
        tmp_output = tmp.name

    # Run hmmscan
    cmd = [
        "hmmscan",
        "--cpu",
        str(cpus),
        "--domtblout",
        tmp_output,
        "--cut_ga",  # Use Pfam gathering thresholds
        pfam_db,
        fasta_file,
    ]

    click.echo(f"Command: {' '.join(cmd)}")

    try:
        # Don't use text=True to avoid UTF-8 decode errors with binary output
        result = subprocess.run(cmd, check=True, capture_output=True)
        click.echo("hmmscan completed successfully")
    except subprocess.CalledProcessError as e:
        click.echo(f"Error running hmmscan: {e}", err=True)
        # Decode stderr with error handling
        stderr = e.stderr.decode('utf-8', errors='replace') if e.stderr else ''
        click.echo(f"stderr: {stderr}", err=True)
        raise

    return tmp_output


def parse_domtblout(domtblout_file):
    """Parse HMMER domtblout format."""
    annotations = []

    with open(domtblout_file, "r") as f:
        for line in f:
            # Skip comments
            if line.startswith("#"):
                continue

            # Parse domain table output
            # Format: target name, accession, tlen, query name, accession, qlen,
            #         E-value, score, bias, #, of, c-Evalue, i-Evalue, score, bias,
            #         from, to, from, to, from, to, acc, description
            fields = line.split()
            if len(fields) < 23:
                continue

            target_name = fields[0]  # Pfam domain name
            target_acc = fields[1]  # Pfam accession
            query_name = fields[3]  # Protein ID
            full_evalue = float(fields[6])
            domain_evalue = float(fields[12])
            ali_from = int(fields[17])  # Alignment start in query
            ali_to = int(fields[18])  # Alignment end in query

            # Get description (everything after position 22)
            description = " ".join(fields[22:])

            annotations.append(
                {
                    "protein_id": query_name,
                    "domain_name": target_name,
                    "domain_accession": target_acc,
                    "domain_description": description,
                    "ali_start": ali_from,
                    "ali_end": ali_to,
                    "full_evalue": full_evalue,
                    "domain_evalue": domain_evalue,
                }
            )

    return annotations


def download_pfam_if_needed(pfam_url, output_path):
    """Download Pfam database if not present."""
    output_path = Path(output_path)

    if output_path.exists():
        click.echo(f"Pfam database found at {output_path}")
        return str(output_path)

    click.echo(f"Downloading Pfam from {pfam_url}...")
    output_path.parent.mkdir(parents=True, exist_ok=True)

    import urllib.request

    urllib.request.urlretrieve(pfam_url, output_path)
    click.echo(f"Downloaded to {output_path}")

    # If gzipped, decompress
    if str(output_path).endswith(".gz"):
        import gzip
        import shutil

        click.echo("Decompressing...")
        with gzip.open(output_path, "rb") as f_in:
            with open(str(output_path).replace(".gz", ""), "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        output_path = Path(str(output_path).replace(".gz", ""))

    # Press the HMM database
    click.echo("Pressing HMM database...")
    subprocess.run(["hmmpress", str(output_path)], check=True)

    return str(output_path)


@click.command()
@click.option(
    "--sequences",
    "sequences_file",
    required=True,
    help="Input FASTA of GP-containing sequences (can be gzipped)",
)
@click.option(
    "--annotations",
    "annotations_file",
    required=True,
    help="Output TSV of domain annotations (gzipped)",
)
@click.option(
    "--pfam-db", default="data/pfam/Pfam-A.hmm", help="Path to Pfam HMM database"
)
@click.option(
    "--pfam-url",
    default="http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz",
    help="URL to download Pfam if not present",
)
@click.option("--evalue", default=1e-3, help="E-value threshold")
@click.option("--cpus", default=8, help="Number of CPUs for hmmscan")
def main(sequences_file, annotations_file, pfam_db, pfam_url, evalue, cpus):
    """Annotate domains in GP-containing proteins using Pfam."""

    click.echo("=" * 60)
    click.echo("Domain Annotation with Pfam")
    click.echo("=" * 60)

    # Download Pfam if needed
    pfam_db = download_pfam_if_needed(pfam_url, pfam_db)

    # If input is gzipped, decompress to temp file
    if sequences_file.endswith(".gz"):
        click.echo("Decompressing input sequences...")
        with tempfile.NamedTemporaryFile(
            mode="w", delete=False, suffix=".fasta"
        ) as tmp:
            tmp_fasta = tmp.name
            with gzip.open(sequences_file, "rt") as f_in:
                tmp.write(f_in.read())
        sequences_file = tmp_fasta

    # Count sequences
    n_sequences = sum(1 for _ in SeqIO.parse(sequences_file, "fasta"))
    click.echo(f"Input sequences: {n_sequences:,}")

    # Run hmmscan
    domtblout = run_hmmscan(sequences_file, pfam_db, evalue, cpus)

    # Parse results
    click.echo("Parsing domain annotations...")
    annotations = parse_domtblout(domtblout)

    # Convert to DataFrame
    df = pd.DataFrame(annotations)

    # Calculate domain length
    if not df.empty:
        df["domain_length"] = df["ali_end"] - df["ali_start"] + 1

    # Save output
    click.echo(f"Saving annotations to {annotations_file}...")
    df.to_csv(annotations_file, sep="\t", index=False, compression="gzip")

    # Summary statistics
    click.echo("\n" + "=" * 60)
    click.echo("Domain Annotation Summary")
    click.echo("=" * 60)
    click.echo(f"Total sequences: {n_sequences:,}")

    if not df.empty:
        n_annotated = df["protein_id"].nunique()
        n_domains = len(df)
        click.echo(
            f"Sequences with domains: {n_annotated:,} ({n_annotated / n_sequences * 100:.1f}%)"
        )
        click.echo(f"Total domain hits: {n_domains:,}")
        click.echo(f"Average domains per protein: {n_domains / n_annotated:.2f}")
        click.echo(f"\nTop 10 most common domains:")
        top_domains = df["domain_name"].value_counts().head(10)
        for domain, count in top_domains.items():
            click.echo(f"  {domain}: {count}")
    else:
        click.echo("No domains found!")

    # Cleanup
    import os

    if domtblout:
        os.unlink(domtblout)
    if sequences_file.endswith(".fasta") and "tmp" in sequences_file:
        os.unlink(sequences_file)


if __name__ == "__main__":
    main()
