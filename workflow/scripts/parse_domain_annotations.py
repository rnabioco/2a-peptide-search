#!/usr/bin/env python3
"""
Parse HMMER domtblout output and combine with GP motif data.
"""

import gzip
import click
import pandas as pd


def parse_domtblout(domtblout_file):
    """Parse HMMER domtblout format."""
    annotations = []
    
    with open(domtblout_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
                
            fields = line.split()
            if len(fields) < 23:
                continue
                
            annotations.append({
                'protein_id': fields[0],
                'protein_length': int(fields[2]),
                'domain_id': fields[1],
                'domain_name': fields[3],
                'domain_start': int(fields[17]),
                'domain_end': int(fields[18]),
                'evalue': float(fields[6]),
                'score': float(fields[7]),
            })
    
    return pd.DataFrame(annotations)


@click.command()
@click.option('--domtblout', required=True, help='HMMER domtblout file')
@click.option('--motifs', required=True, help='GP motifs TSV file')
@click.option('--annotations', required=True, help='Output annotations TSV')
def main(domtblout, motifs, annotations):
    """Parse domain annotations and combine with GP motif data."""
    
    click.echo("Parsing domtblout file...")
    domains_df = parse_domtblout(domtblout)
    
    click.echo(f"Found {len(domains_df)} domain annotations")
    click.echo(f"Covering {domains_df['protein_id'].nunique()} proteins")
    
    # Group domains by protein
    protein_domains = domains_df.groupby('protein_id').apply(
        lambda x: x.to_dict('records')
    ).to_dict()
    
    click.echo("Loading GP motifs...")
    motifs_df = pd.read_csv(motifs, sep='\t', compression='gzip')
    
    # Add domain information to each motif
    motifs_df['domains'] = motifs_df['protein_id'].map(
        lambda pid: protein_domains.get(pid, [])
    )
    motifs_df['n_domains'] = motifs_df['domains'].apply(len)
    
    # Save
    motifs_df.to_csv(annotations, sep='\t', index=False, compression='gzip')
    
    click.echo(f"Saved annotations to {annotations}")
    click.echo(f"Proteins with domains: {(motifs_df['n_domains'] > 0).sum()}")


if __name__ == '__main__':
    main()
