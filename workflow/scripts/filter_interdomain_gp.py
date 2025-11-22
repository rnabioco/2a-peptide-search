#!/usr/bin/env python3
"""
Filter GP motifs based on domain boundary annotations.

Identifies GP motifs that occur:
- Between domains (inter-domain)
- Within domains (intra-domain)
- In unstructured regions
"""

import gzip
import pandas as pd
import click


def classify_gp_position(gp_pos, protein_length, domains):
    """
    Classify GP position relative to domains.

    Returns:
        classification: 'inter-domain', 'intra-domain', 'unstructured'
        details: dict with additional info
    """
    if not domains:
        return 'unstructured', {'reason': 'no_domains'}

    # Sort domains by start position
    sorted_domains = sorted(domains, key=lambda x: x['start'])

    # Check if GP is within any domain
    for domain in sorted_domains:
        if domain['start'] <= gp_pos < domain['end']:
            return 'intra-domain', {
                'domain_id': domain['id'],
                'domain_name': domain['name'],
                'position_in_domain': gp_pos - domain['start'],
                'domain_length': domain['end'] - domain['start']
            }

    # Check if GP is between domains
    for i in range(len(sorted_domains) - 1):
        domain1 = sorted_domains[i]
        domain2 = sorted_domains[i + 1]

        if domain1['end'] <= gp_pos < domain2['start']:
            return 'inter-domain', {
                'upstream_domain': domain1['name'],
                'downstream_domain': domain2['name'],
                'distance_to_upstream': gp_pos - domain1['end'],
                'distance_to_downstream': domain2['start'] - gp_pos,
                'linker_length': domain2['start'] - domain1['end']
            }

    # GP is in unstructured region (before first or after last domain)
    if gp_pos < sorted_domains[0]['start']:
        return 'unstructured', {
            'reason': 'before_first_domain',
            'distance_to_domain': sorted_domains[0]['start'] - gp_pos
        }
    elif gp_pos >= sorted_domains[-1]['end']:
        return 'unstructured', {
            'reason': 'after_last_domain',
            'distance_to_domain': gp_pos - sorted_domains[-1]['end']
        }

    return 'unstructured', {'reason': 'unknown'}


@click.command()
@click.option('--motifs', required=True, help='GP motifs TSV')
@click.option('--annotations', required=True, help='Domain annotations TSV')
@click.option('--interdomain-out', required=True, help='Output inter-domain GP motifs')
@click.option('--intradomain-out', required=True, help='Output intra-domain GP motifs')
@click.option('--stats-out', required=True, help='Output statistics file')
def main(motifs, annotations, interdomain_out, intradomain_out, stats_out):
    """Filter GP motifs by domain boundaries."""

    click.echo("Loading GP motifs...")
    motifs_df = pd.read_csv(motifs, sep='\t', compression='gzip')

    click.echo("Loading domain annotations...")
    # Assuming annotations format: protein_id, domain_id, domain_name, start, end
    annot_df = pd.read_csv(annotations, sep='\t', compression='gzip')

    # Group domains by protein
    protein_domains = {}
    for _, row in annot_df.iterrows():
        protein_id = row['protein_id']
        if protein_id not in protein_domains:
            protein_domains[protein_id] = []
        protein_domains[protein_id].append({
            'id': row['domain_id'],
            'name': row['domain_name'],
            'start': row['start'],
            'end': row['end']
        })

    click.echo("Classifying GP motif positions...")

    classifications = []
    inter_domain_motifs = []
    intra_domain_motifs = []

    for _, motif in motifs_df.iterrows():
        protein_id = motif['protein_id']
        gp_pos = motif['gp_position']
        protein_len = motif['protein_length']

        domains = protein_domains.get(protein_id, [])

        classification, details = classify_gp_position(gp_pos, protein_len, domains)

        motif_with_class = motif.to_dict()
        motif_with_class['classification'] = classification
        motif_with_class.update({f'detail_{k}': v for k, v in details.items()})

        classifications.append(motif_with_class)

        if classification == 'inter-domain':
            inter_domain_motifs.append(motif_with_class)
        elif classification == 'intra-domain':
            intra_domain_motifs.append(motif_with_class)

    # Save classified motifs
    click.echo(f"\nSaving inter-domain motifs ({len(inter_domain_motifs)})...")
    pd.DataFrame(inter_domain_motifs).to_csv(
        interdomain_out, sep='\t', index=False, compression='gzip'
    )

    click.echo(f"Saving intra-domain motifs ({len(intra_domain_motifs)})...")
    pd.DataFrame(intra_domain_motifs).to_csv(
        intradomain_out, sep='\t', index=False, compression='gzip'
    )

    # Calculate statistics
    stats = {
        'total_gp_motifs': len(classifications),
        'inter_domain': sum(1 for c in classifications if c['classification'] == 'inter-domain'),
        'intra_domain': sum(1 for c in classifications if c['classification'] == 'intra-domain'),
        'unstructured': sum(1 for c in classifications if c['classification'] == 'unstructured'),
    }

    stats_df = pd.DataFrame([stats])
    stats_df['inter_domain_pct'] = stats_df['inter_domain'] / stats_df['total_gp_motifs'] * 100
    stats_df['intra_domain_pct'] = stats_df['intra_domain'] / stats_df['total_gp_motifs'] * 100
    stats_df['unstructured_pct'] = stats_df['unstructured'] / stats_df['total_gp_motifs'] * 100

    stats_df.to_csv(stats_out, sep='\t', index=False)

    # Print summary
    click.echo("\n" + "="*60)
    click.echo("GP Motif Classification Summary")
    click.echo("="*60)
    click.echo(f"Total GP motifs: {stats['total_gp_motifs']:,}")
    click.echo(f"Inter-domain: {stats['inter_domain']:,} ({stats_df['inter_domain_pct'].values[0]:.1f}%)")
    click.echo(f"Intra-domain: {stats['intra_domain']:,} ({stats_df['intra_domain_pct'].values[0]:.1f}%)")
    click.echo(f"Unstructured: {stats['unstructured']:,} ({stats_df['unstructured_pct'].values[0]:.1f}%)")


if __name__ == '__main__':
    main()
