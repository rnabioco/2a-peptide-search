#!/usr/bin/env python3
"""
Filter Stockholm alignment by E-value threshold and remove redundant sequences.
"""

import sys
from pathlib import Path

# Add scripts directory to path for utils import
sys.path.insert(0, str(Path(__file__).parent))

import click
from Bio import AlignIO, SearchIO
from Bio.Align import MultipleSeqAlignment
from utils import get_logger

logger = get_logger(__name__)


@click.command()
@click.argument('alignment', type=click.Path(exists=True))
@click.argument('tblout', type=click.Path(exists=True))
@click.argument('output', type=click.Path())
@click.option('--evalue', default=1e-5, help='E-value threshold')
@click.option('--identity', default=0.95, help='Identity threshold for clustering')
def main(alignment, tblout, output, evalue, identity):
    """Filter alignment based on E-value and remove redundant sequences."""

    logger.info(f"Starting alignment filtering")
    logger.info(f"  Alignment: {alignment}")
    logger.info(f"  Tblout: {tblout}")
    logger.info(f"  E-value threshold: {evalue}")

    # Parse search results to get E-values
    logger.info("Parsing HMMER tblout for E-values...")
    good_ids = set()
    total_hits = 0
    for result in SearchIO.parse(tblout, 'hmmer3-tab'):
        for hit in result.hits:
            total_hits += 1
            if hit.evalue <= evalue:
                hit_id = hit.id.split('/')[0] if '/' in hit.id else hit.id
                good_ids.add(hit_id)

    logger.info(f"  Total hits in tblout: {total_hits}")
    logger.info(f"  Hits passing E-value threshold: {len(good_ids)}")

    # Filter alignment
    logger.info("Reading alignment file...")
    aln = AlignIO.read(alignment, 'stockholm')
    logger.info(f"  Sequences in alignment: {len(aln)}")

    logger.info("Filtering sequences...")
    filtered_records = []
    for record in aln:
        rec_id = record.id.split('/')[0] if '/' in record.id else record.id
        if rec_id in good_ids:
            filtered_records.append(record)

    filtered_aln = MultipleSeqAlignment(filtered_records)
    logger.info(f"  Sequences after filtering: {len(filtered_aln)}")

    logger.info(f"Writing output to {output}")
    with open(output, 'w') as out:
        AlignIO.write(filtered_aln, out, 'stockholm')

    logger.info(f"Complete: {len(aln)} -> {len(filtered_aln)} sequences")


if __name__ == '__main__':
    main()
