#!/usr/bin/env python3
"""
Sanitize FASTA files by removing invalid IUPAC characters.

Processes large gzipped FASTA files in streaming mode to handle
files with hundreds of millions of lines efficiently.

Valid protein IUPAC characters:
- Standard: A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y
- Ambiguous: B (D/N), Z (E/Q), X (any), J (I/L)
- Special: U (selenocysteine), O (pyrrolysine), * (stop)

Invalid characters (gaps, etc.) are removed from sequences.
"""

import gzip
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

import click
from utils import get_logger

logger = get_logger(__name__)

# Valid IUPAC protein characters (case-insensitive)
VALID_PROTEIN_CHARS = set('ACDEFGHIKLMNPQRSTVWYBZXJUOacdefghiklmnpqrstvwybzxjuo*')

# Regex to remove invalid characters (anything not in valid set)
INVALID_CHAR_PATTERN = re.compile(r'[^ACDEFGHIKLMNPQRSTVWYBZXJUOacdefghiklmnpqrstvwybzxjuo*]')


def open_file(path, mode='rt'):
    """Open file, handling gzip transparently."""
    if str(path).endswith('.gz'):
        return gzip.open(path, mode)
    return open(path, mode)


def sanitize_sequence(seq: str) -> str:
    """Remove invalid characters from sequence."""
    return INVALID_CHAR_PATTERN.sub('', seq)


@click.command()
@click.argument('input_fasta', type=click.Path(exists=True))
@click.argument('output_fasta', type=click.Path())
@click.option('--report-interval', default=10_000_000,
              help='Report progress every N lines')
@click.option('--dry-run', is_flag=True,
              help='Count invalid characters without writing output')
def main(input_fasta, output_fasta, report_interval, dry_run):
    """
    Sanitize FASTA file by removing invalid IUPAC characters.

    Handles gzipped files automatically based on .gz extension.
    Processes in streaming mode for memory efficiency.
    """
    logger.info(f"Sanitizing FASTA: {input_fasta}")
    logger.info(f"Output: {output_fasta}")
    if dry_run:
        logger.info("DRY RUN - no output will be written")

    lines_processed = 0
    sequences_processed = 0
    chars_removed = 0
    sequences_with_invalid = 0
    current_header = None

    # Determine compression for output
    out_gzip = str(output_fasta).endswith('.gz')

    with open_file(input_fasta, 'rt') as infile:
        outfile = None
        if not dry_run:
            outfile = gzip.open(output_fasta, 'wt') if out_gzip else open(output_fasta, 'w')

        try:
            for line in infile:
                lines_processed += 1

                if line.startswith('>'):
                    # Header line - pass through unchanged
                    current_header = line.strip()
                    sequences_processed += 1
                    if outfile:
                        outfile.write(line)
                else:
                    # Sequence line - sanitize
                    original = line.rstrip('\n\r')
                    sanitized = sanitize_sequence(original)

                    removed = len(original) - len(sanitized)
                    if removed > 0:
                        chars_removed += removed
                        if removed == len(original) - len(sanitized):
                            # First time we see invalid chars in this sequence
                            sequences_with_invalid += 1
                            if sequences_with_invalid <= 10:
                                logger.warning(
                                    f"Invalid chars in {current_header[:50]}...: "
                                    f"removed {removed} chars"
                                )

                    if outfile and sanitized:  # Don't write empty lines
                        outfile.write(sanitized + '\n')

                # Progress report
                if lines_processed % report_interval == 0:
                    logger.info(
                        f"Progress: {lines_processed:,} lines, "
                        f"{sequences_processed:,} sequences, "
                        f"{chars_removed:,} invalid chars removed"
                    )

        finally:
            if outfile:
                outfile.close()

    # Final summary
    logger.info("=" * 60)
    logger.info(f"Complete: {lines_processed:,} lines processed")
    logger.info(f"  Sequences: {sequences_processed:,}")
    logger.info(f"  Invalid characters removed: {chars_removed:,}")
    logger.info(f"  Sequences with invalid chars: {sequences_with_invalid:,}")

    if chars_removed > 0:
        logger.warning(
            f"Removed {chars_removed:,} invalid characters from "
            f"{sequences_with_invalid:,} sequences"
        )
    else:
        logger.info("No invalid characters found - file is clean")

    return 0


if __name__ == '__main__':
    main()
