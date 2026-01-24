#!/usr/bin/env python3
"""
Parse MMseqs2 clustering results and combine with GP motif data.
"""

import click
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


@click.command()
@click.option('--cluster-tsv', required=True, help='MMseqs2 cluster TSV')
@click.option('--motifs', required=True, help='GP motifs TSV (all or inter-domain)')
@click.option('--clusters-out', required=True, help='Output clusters TSV')
@click.option('--representatives-out', required=True, help='Output representatives FASTA')
def main(cluster_tsv, motifs, clusters_out, representatives_out):
    """Parse MMseqs2 clustering and assign cluster IDs."""

    click.echo("Loading clustering results...")
    cluster_df = pd.read_csv(cluster_tsv, sep='\t', names=['representative', 'member'])

    click.echo("Loading GP motifs...")
    motifs_df = pd.read_csv(motifs, sep='\t', compression='gzip')

    click.echo("Creating sequence IDs...")
    # Create sequence ID mapping (replace | with _ to match FASTA IDs)
    motifs_df['seq_id'] = (motifs_df['protein_id'].str.replace('|', '_', regex=False) + '_GP' +
                            motifs_df['gp_index'].astype(str) + '_pos' +
                            motifs_df['gp_position'].astype(str))

    click.echo("Assigning cluster IDs...")
    # Assign numeric cluster IDs to representatives (much faster than loop)
    rep_to_cluster = pd.DataFrame({
        'representative': cluster_df['representative'].unique()
    })
    rep_to_cluster['cluster_id'] = range(1, len(rep_to_cluster) + 1)

    # Map representatives to cluster IDs, then map members via representatives
    cluster_df = cluster_df.merge(rep_to_cluster, on='representative')

    # Map by SEQUENCE CONTENT instead of seq_id to handle deduplication
    # After deduplication, only unique sequences exist in cluster file,
    # but we need to assign ALL motifs (including duplicates) to clusters

    # Get unique seq_id -> context_sequence mapping from motifs
    seq_id_to_sequence = motifs_df[['seq_id', 'context_sequence']].drop_duplicates('seq_id')

    # Add sequence content to cluster members
    cluster_with_seqs = cluster_df.merge(
        seq_id_to_sequence,
        left_on='member',
        right_on='seq_id',
        how='left'
    )

    # Create sequence -> cluster_id mapping
    # Each unique sequence maps to one cluster (via its representative)
    seq_to_cluster = cluster_with_seqs.set_index('context_sequence')['cluster_id'].to_dict()

    # Map ALL motifs by their sequence content (handles duplicates correctly)
    motifs_df['cluster_id'] = motifs_df['context_sequence'].map(seq_to_cluster)

    n_unique_clustered = len(cluster_df)
    click.echo(f"Clustered {n_unique_clustered} unique sequences")

    click.echo("Computing cluster sizes...")
    # Get cluster sizes using value_counts (faster than groupby.size)
    # dropna=False to include NaN counts for debugging
    cluster_sizes = motifs_df['cluster_id'].value_counts(dropna=True).to_dict()
    motifs_df['cluster_size'] = motifs_df['cluster_id'].map(cluster_sizes)

    n_clusters = len(rep_to_cluster)
    n_assigned = motifs_df['cluster_id'].notna().sum()
    n_unassigned = motifs_df['cluster_id'].isna().sum()
    click.echo(f"Found {len(seq_to_cluster)} unique sequences in {n_clusters} clusters")
    click.echo(f"Assigned: {n_assigned}, Unassigned: {n_unassigned}")
    if cluster_sizes:
        click.echo(f"Largest cluster: {max(cluster_sizes.values())} members")
    else:
        click.echo("Warning: No clusters found!")

    click.echo("Saving clusters...")
    # Save clusters
    motifs_df.to_csv(clusters_out, sep='\t', index=False, compression='gzip')

    click.echo("Extracting representative sequences...")
    # Extract representative sequences - use merge instead of loop
    # Get unique representatives with their cluster IDs
    rep_df = cluster_df.drop_duplicates('representative')[['representative', 'cluster_id']]
    rep_df = rep_df.rename(columns={'representative': 'seq_id'})

    # Merge with motifs to get sequences
    rep_with_seqs = rep_df.merge(
        motifs_df[['seq_id', 'context_sequence']].drop_duplicates('seq_id'),
        on='seq_id'
    )

    # Create SeqRecords
    rep_seqs = []
    for _, row in rep_with_seqs.iterrows():
        size = cluster_sizes.get(row['cluster_id'], 0)
        record = SeqRecord(
            Seq(row['context_sequence']),
            id=f"cluster_{row['cluster_id']}",
            description=f"representative={row['seq_id']} size={size}"
        )
        rep_seqs.append(record)

    SeqIO.write(rep_seqs, representatives_out, 'fasta')
    click.echo(f"Wrote {len(rep_seqs)} representative sequences")


if __name__ == '__main__':
    main()
