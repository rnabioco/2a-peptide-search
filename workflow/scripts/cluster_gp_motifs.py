#!/usr/bin/env python3
"""
Cluster GP motifs by sequence context.

Uses sequence similarity of the context around GP to identify
related stalling peptide families.
"""

import gzip
import subprocess
import tempfile
from pathlib import Path

import click
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def write_fasta_for_clustering(motifs_df, fasta_path):
    """Write context sequences to FASTA for clustering."""
    records = []
    for idx, row in motifs_df.iterrows():
        seq_id = f"{row['protein_id']}_GP{row['gp_index']}_pos{row['gp_position']}"
        record = SeqRecord(
            Seq(row['context_sequence']),
            id=seq_id,
            description=f"{row['protein_description']}"
        )
        records.append(record)

    SeqIO.write(records, fasta_path, 'fasta')
    return len(records)


def run_mmseqs_clustering(fasta_path, output_prefix, identity=0.7, coverage=0.8, threads=8):
    """Run MMseqs2 clustering."""

    click.echo(f"Running MMseqs2 clustering (identity={identity}, coverage={coverage})...")

    # Create MMseqs2 database
    db_path = f"{output_prefix}_db"
    subprocess.run([
        'mmseqs', 'createdb',
        str(fasta_path),
        db_path
    ], check=True)

    # Cluster
    cluster_path = f"{output_prefix}_cluster"
    tmp_dir = f"{output_prefix}_tmp"
    Path(tmp_dir).mkdir(exist_ok=True)

    subprocess.run([
        'mmseqs', 'cluster',
        db_path,
        cluster_path,
        tmp_dir,
        '--min-seq-id', str(identity),
        '-c', str(coverage),
        '--threads', str(threads)
    ], check=True)

    # Create TSV output
    cluster_tsv = f"{output_prefix}_cluster.tsv"
    subprocess.run([
        'mmseqs', 'createtsv',
        db_path,
        db_path,
        cluster_path,
        cluster_tsv
    ], check=True)

    return cluster_tsv


@click.command()
@click.option('--interdomain', required=True, help='Inter-domain GP motifs TSV')
@click.option('--clusters-out', required=True, help='Output clusters TSV')
@click.option('--representatives-out', required=True, help='Output cluster representatives FASTA')
@click.option('--identity', default=0.7, help='Sequence identity threshold')
@click.option('--coverage', default=0.8, help='Coverage threshold')
@click.option('--threads', default=8, help='Number of threads')
def main(interdomain, clusters_out, representatives_out, identity, coverage, threads):
    """Cluster GP motifs by sequence context."""

    click.echo("Loading inter-domain GP motifs...")
    motifs_df = pd.read_csv(interdomain, sep='\t', compression='gzip')

    click.echo(f"Found {len(motifs_df):,} inter-domain GP motifs")

    # Create temporary directory for clustering
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        # Write sequences for clustering
        fasta_path = tmpdir / "gp_motifs.fasta"
        click.echo(f"Writing sequences to {fasta_path}...")
        n_seqs = write_fasta_for_clustering(motifs_df, fasta_path)
        click.echo(f"Wrote {n_seqs:,} sequences")

        # Run clustering
        output_prefix = tmpdir / "mmseqs_out"
        cluster_tsv = run_mmseqs_clustering(
            fasta_path, output_prefix,
            identity, coverage, threads
        )

        # Parse clustering results
        click.echo("Parsing clustering results...")
        cluster_df = pd.read_csv(
            cluster_tsv, sep='\t',
            names=['representative', 'member']
        )

        # Map back to original motifs
        seq_id_to_motif = {}
        for idx, row in motifs_df.iterrows():
            seq_id = f"{row['protein_id']}_GP{row['gp_index']}_pos{row['gp_position']}"
            seq_id_to_motif[seq_id] = idx

        # Assign cluster IDs
        clusters = {}
        cluster_id = 1

        for rep in cluster_df['representative'].unique():
            members = cluster_df[cluster_df['representative'] == rep]['member'].tolist()

            clusters[cluster_id] = {
                'representative': rep,
                'members': members,
                'size': len(members)
            }
            cluster_id += 1

        # Add cluster assignments to motifs
        motif_clusters = []
        for cid, cinfo in clusters.items():
            for member in cinfo['members']:
                if member in seq_id_to_motif:
                    motif_idx = seq_id_to_motif[member]
                    motif = motifs_df.iloc[motif_idx].to_dict()
                    motif['cluster_id'] = cid
                    motif['cluster_size'] = cinfo['size']
                    motif['is_representative'] = (member == cinfo['representative'])
                    motif_clusters.append(motif)

        # Save clusters
        click.echo(f"\nSaving {len(clusters)} clusters...")
        cluster_df = pd.DataFrame(motif_clusters)
        cluster_df.to_csv(clusters_out, sep='\t', index=False, compression='gzip')

        # Extract representative sequences
        click.echo("Extracting cluster representatives...")
        representatives = []
        for cid, cinfo in clusters.items():
            rep_id = cinfo['representative']
            if rep_id in seq_id_to_motif:
                motif_idx = seq_id_to_motif[rep_id]
                motif = motifs_df.iloc[motif_idx]

                record = SeqRecord(
                    Seq(motif['context_sequence']),
                    id=f"cluster_{cid}",
                    description=f"size={cinfo['size']} rep={rep_id}"
                )
                representatives.append(record)

        SeqIO.write(representatives, representatives_out, 'fasta')

    # Print summary
    click.echo("\n" + "="*60)
    click.echo("GP Motif Clustering Summary")
    click.echo("="*60)
    click.echo(f"Input motifs: {len(motifs_df):,}")
    click.echo(f"Number of clusters: {len(clusters):,}")
    click.echo(f"Clustered motifs: {len(motif_clusters):,}")

    # Cluster size distribution
    cluster_sizes = [c['size'] for c in clusters.values()]
    click.echo(f"\nCluster size distribution:")
    click.echo(f"  Mean: {sum(cluster_sizes)/len(cluster_sizes):.1f}")
    click.echo(f"  Median: {sorted(cluster_sizes)[len(cluster_sizes)//2]}")
    click.echo(f"  Max: {max(cluster_sizes)}")
    click.echo(f"  Singleton clusters: {sum(1 for s in cluster_sizes if s == 1)}")


if __name__ == '__main__':
    main()
