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
    
    # Create sequence ID mapping
    motifs_df['seq_id'] = (motifs_df['protein_id'] + '_GP' + 
                            motifs_df['gp_index'].astype(str) + '_pos' + 
                            motifs_df['gp_position'].astype(str))
    
    # Assign cluster IDs
    cluster_map = {}
    cluster_id = 0
    for rep in cluster_df['representative'].unique():
        cluster_id += 1
        members = cluster_df[cluster_df['representative'] == rep]['member'].tolist()
        for member in members:
            cluster_map[member] = cluster_id
    
    motifs_df['cluster_id'] = motifs_df['seq_id'].map(cluster_map)
    
    # Get cluster sizes
    cluster_sizes = motifs_df.groupby('cluster_id').size().to_dict()
    motifs_df['cluster_size'] = motifs_df['cluster_id'].map(cluster_sizes)
    
    click.echo(f"Found {len(cluster_map)} sequences in {cluster_id} clusters")
    click.echo(f"Largest cluster: {max(cluster_sizes.values())} members")
    
    # Save clusters
    motifs_df.to_csv(clusters_out, sep='\t', index=False, compression='gzip')
    
    # Extract representative sequences
    representatives = cluster_df.groupby('representative').first().reset_index()
    rep_seqs = []
    
    for _, row in representatives.iterrows():
        seq_id = row['representative']
        motif_row = motifs_df[motifs_df['seq_id'] == seq_id].iloc[0]
        
        record = SeqRecord(
            Seq(motif_row['context_sequence']),
            id=f"cluster_{motif_row['cluster_id']}",
            description=f"representative={seq_id} size={cluster_sizes[motif_row['cluster_id']]}"
        )
        rep_seqs.append(record)
    
    SeqIO.write(rep_seqs, representatives_out, 'fasta')
    click.echo(f"Wrote {len(rep_seqs)} representative sequences")


if __name__ == '__main__':
    main()
