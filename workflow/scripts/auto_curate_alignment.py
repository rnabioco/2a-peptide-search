#!/usr/bin/env python3
"""
Automatically curate Stockholm alignments for 2A peptide models.

Applies multiple quality filters to produce publication-ready alignments
without manual intervention:
- E-value/bit score filtering
- Sequence length filtering (remove truncated sequences)
- Gap percentage filtering
- C-terminal motif validation (PGP/NPG conservation)
- Outlier removal based on bit score distribution
"""

import re
import sys
from pathlib import Path

# Add scripts directory to path for utils import
sys.path.insert(0, str(Path(__file__).parent))

import click
import numpy as np
from Bio import AlignIO, SearchIO
from Bio.Align import MultipleSeqAlignment
from utils import get_logger

logger = get_logger(__name__)


def parse_hmmer_scores(tblout_path):
    """Parse HMMER tblout to get E-values and bit scores per hit."""
    scores = {}
    try:
        for result in SearchIO.parse(tblout_path, 'hmmer3-tab'):
            for hit in result.hits:
                # Extract base ID without coordinates
                hit_id = hit.id.split('/')[0] if '/' in hit.id
                scores[hit_id] = {
                    'evalue': hit.evalue,
                    'bitscore': hit.bitscore,
                }
                # Also store with full ID for coordinate-based matching
                scores[hit.id] = scores[hit_id]
    except Exception as e:
        click.echo(f"Warning: Could not parse tblout: {e}", err=True)
    return scores


def calculate_gap_percentage(seq_str):
    """Calculate percentage of gaps in a sequence."""
    gaps = seq_str.count('-') + seq_str.count('.')
    return gaps / len(seq_str) if len(seq_str) > 0 else 1.0


def get_ungapped_length(seq_str):
    """Get length of sequence without gaps."""
    return len(seq_str.replace('-', '').replace('.', ''))


def check_cterminal_motif(seq_str, motif_patterns=None):
    """
    Check for conserved C-terminal motif.

    2A peptides should end with PGP, NPGP, or similar motifs.
    Returns True if a valid C-terminal motif is found.
    """
    if motif_patterns is None:
        # Common 2A C-terminal patterns
        motif_patterns = [
            r'[NHP]PGP$',      # NPGP, PPGP, HPGP at very end
            r'[NHP]PGP[A-Z]$', # With one trailing residue
            r'PGP$',           # Minimal PGP
            r'PGP[A-Z]{1,3}$', # PGP with short tail
        ]

    # Remove gaps and get C-terminal region
    ungapped = seq_str.replace('-', '').replace('.', '').upper()
    if len(ungapped) < 4:
        return False

    c_term = ungapped[-10:]  # Check last 10 residues

    for pattern in motif_patterns:
        if re.search(pattern, c_term):
            return True
    return False


def filter_by_bitscore_percentile(records, scores, percentile=10):
    """Remove sequences with bit scores below the given percentile."""
    if not scores:
        return records

    # Get bit scores for records we have
    record_scores = []
    for rec in records:
        rec_id = rec.id.split('/')[0] if '/' in rec.id else rec.id
        if rec_id in scores:
            record_scores.append((rec, scores[rec_id]['bitscore']))
        elif rec.id in scores:
            record_scores.append((rec, scores[rec.id]['bitscore']))
        else:
            # Keep records without scores (conservative)
            record_scores.append((rec, float('inf')))

    if not record_scores:
        return records

    # Calculate percentile threshold
    bitscores = [s for _, s in record_scores if s != float('inf')]
    if not bitscores:
        return records

    threshold = np.percentile(bitscores, percentile)

    return [rec for rec, score in record_scores if score >= threshold]


@click.command()
@click.argument('alignment', type=click.Path(exists=True))
@click.argument('output', type=click.Path())
@click.option('--tblout', type=click.Path(exists=True), default=None,
              help='HMMER tblout file for score-based filtering')
@click.option('--evalue', default=1e-5, help='E-value threshold')
@click.option('--min-length', default=15, help='Minimum ungapped sequence length')
@click.option('--max-gap-pct', default=0.5, help='Maximum gap percentage (0-1)')
@click.option('--require-motif/--no-require-motif', default=True,
              help='Require C-terminal PGP motif')
@click.option('--bitscore-percentile', default=10,
              help='Remove sequences below this bit score percentile')
@click.option('--min-sequences', default=3,
              help='Minimum sequences required in output')
def main(alignment, output, tblout, evalue, min_length, max_gap_pct,
         require_motif, bitscore_percentile, min_sequences):
    """
    Auto-curate alignment for 2A peptide HMM building.

    Applies quality filters to produce a clean alignment suitable for
    building production HMM models without manual curation.
    """
    logger.info("Starting auto-curation")
    logger.info(f"  Alignment: {alignment}")
    logger.info(f"  Output: {output}")
    logger.info(f"  E-value threshold: {evalue}")
    logger.info(f"  Min length: {min_length}")
    logger.info(f"  Max gap %: {max_gap_pct*100:.0f}%")
    logger.info(f"  Require C-term motif: {require_motif}")

    # Parse scores if tblout provided
    scores = parse_hmmer_scores(tblout) if tblout else {}
    if tblout:
        logger.info(f"  Loaded {len(scores)} scores from tblout")

    # Read alignment
    try:
        aln = AlignIO.read(alignment, 'stockholm')
    except Exception as e:
        logger.error(f"Error reading alignment: {e}")
        sys.exit(1)

    initial_count = len(aln)
    logger.info(f"Input: {initial_count} sequences")

    filtered_records = []
    filter_stats = {
        'evalue': 0,
        'length': 0,
        'gaps': 0,
        'motif': 0,
    }

    logger.info("Applying filters...")
    for record in aln:
        seq_str = str(record.seq)
        rec_id = record.id.split('/')[0] if '/' in record.id else record.id

        # E-value filter (if scores available)
        if scores:
            rec_score = scores.get(rec_id) or scores.get(record.id)
            if rec_score and rec_score['evalue'] > evalue:
                filter_stats['evalue'] += 1
                continue

        # Length filter
        ungapped_len = get_ungapped_length(seq_str)
        if ungapped_len < min_length:
            filter_stats['length'] += 1
            continue

        # Gap percentage filter
        gap_pct = calculate_gap_percentage(seq_str)
        if gap_pct > max_gap_pct:
            filter_stats['gaps'] += 1
            continue

        # C-terminal motif filter
        if require_motif and not check_cterminal_motif(seq_str):
            filter_stats['motif'] += 1
            continue

        filtered_records.append(record)

    # Bit score percentile filter
    if scores and bitscore_percentile > 0:
        before_bitscore = len(filtered_records)
        filtered_records = filter_by_bitscore_percentile(
            filtered_records, scores, bitscore_percentile
        )
        filter_stats['bitscore'] = before_bitscore - len(filtered_records)

    # Report filtering stats
    logger.info("Filtered out:")
    logger.info(f"  - E-value > {evalue}: {filter_stats['evalue']}")
    logger.info(f"  - Length < {min_length}: {filter_stats['length']}")
    logger.info(f"  - Gap % > {max_gap_pct*100:.0f}%: {filter_stats['gaps']}")
    if require_motif:
        logger.info(f"  - Missing C-term motif: {filter_stats['motif']}")
    if 'bitscore' in filter_stats:
        logger.info(f"  - Low bit score: {filter_stats['bitscore']}")

    # Check minimum sequences
    if len(filtered_records) < min_sequences:
        logger.warning(
            f"Only {len(filtered_records)} sequences passed filters "
            f"(minimum: {min_sequences}). Relaxing filters..."
        )
        # Fall back to less strict filtering
        filtered_records = []
        for record in aln:
            seq_str = str(record.seq)
            # Only apply length filter
            if get_ungapped_length(seq_str) >= min_length:
                filtered_records.append(record)

        if len(filtered_records) < min_sequences:
            logger.error(
                f"Cannot produce alignment with >= {min_sequences} sequences"
            )
            # Output whatever we have
            if not filtered_records:
                filtered_records = list(aln)[:min_sequences]

    # Write output
    filtered_aln = MultipleSeqAlignment(filtered_records)

    logger.info(f"Writing output to {output}")
    with open(output, 'w') as out:
        AlignIO.write(filtered_aln, out, 'stockholm')

    logger.info(f"Complete: {initial_count} -> {len(filtered_aln)} sequences")

    return 0


if __name__ == '__main__':
    main()
