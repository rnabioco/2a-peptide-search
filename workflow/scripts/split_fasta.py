#!/usr/bin/env python3
"""Split a large FASTA file into chunks for parallel processing.

This script is designed to be called by Snakemake to split phage genome
files into fixed number of chunks for parallel ORF prediction.

Snakemake provides:
- snakemake.input.genomes: Input gzipped FASTA file
- snakemake.output.chunks: List of output chunk files
- snakemake.params.seqs_per_chunk: Target sequences per chunk (approximate)
- snakemake.log: Log file path
"""

import gzip
import sys
from pathlib import Path

from Bio import SeqIO


def split_fasta(input_file, output_files, log_file):
    """Split input FASTA into the specified output chunk files.

    Distributes sequences evenly across chunks rather than using a fixed
    sequences-per-chunk, since we need exactly the specified number of
    output files.

    Parameters
    ----------
    input_file : str
        Path to input gzipped FASTA file
    output_files : list
        List of output chunk file paths
    log_file : str
        Path to log file
    """
    n_chunks = len(output_files)

    with open(log_file, "w") as log:
        log.write(f"Splitting {input_file} into {n_chunks} chunks\n")

        # First pass: count sequences
        log.write("Counting sequences...\n")
        n_seqs = 0
        with gzip.open(input_file, "rt") as handle:
            for _ in SeqIO.parse(handle, "fasta"):
                n_seqs += 1

        log.write(f"Total sequences: {n_seqs}\n")

        if n_seqs == 0:
            log.write("WARNING: No sequences found in input file\n")
            # Create empty chunk files
            for output_file in output_files:
                Path(output_file).parent.mkdir(parents=True, exist_ok=True)
                with open(output_file, "w") as f:
                    pass
            return

        # Calculate sequences per chunk (distribute evenly)
        base_per_chunk = n_seqs // n_chunks
        remainder = n_seqs % n_chunks

        log.write(f"Base sequences per chunk: {base_per_chunk}\n")
        log.write(f"Remainder: {remainder} (these chunks get 1 extra)\n")

        # Second pass: write sequences to chunks
        log.write("Writing sequences to chunks...\n")

        # Ensure output directories exist
        for output_file in output_files:
            Path(output_file).parent.mkdir(parents=True, exist_ok=True)

        # Open all output files
        handles = [open(f, "w") for f in output_files]

        try:
            with gzip.open(input_file, "rt") as in_handle:
                chunk_idx = 0
                seqs_in_chunk = 0
                # First 'remainder' chunks get base_per_chunk + 1
                target_for_chunk = base_per_chunk + (1 if chunk_idx < remainder else 0)

                for record in SeqIO.parse(in_handle, "fasta"):
                    # Write to current chunk
                    SeqIO.write(record, handles[chunk_idx], "fasta")
                    seqs_in_chunk += 1

                    # Check if we need to move to next chunk
                    if seqs_in_chunk >= target_for_chunk and chunk_idx < n_chunks - 1:
                        chunk_idx += 1
                        seqs_in_chunk = 0
                        target_for_chunk = (
                            base_per_chunk + (1 if chunk_idx < remainder else 0)
                        )

            log.write(f"Successfully created {n_chunks} chunk files\n")

            # Log chunk sizes
            for i, output_file in enumerate(output_files):
                size = Path(output_file).stat().st_size
                log.write(f"  Chunk {i + 1}: {size} bytes\n")

        finally:
            for handle in handles:
                handle.close()


if __name__ == "__main__":
    # When called from Snakemake
    if "snakemake" in dir():
        split_fasta(
            input_file=snakemake.input.genomes,
            output_files=snakemake.output.chunks,
            log_file=str(snakemake.log),
        )
    else:
        # CLI usage for testing
        import argparse

        parser = argparse.ArgumentParser(
            description="Split FASTA file into chunks"
        )
        parser.add_argument("input", help="Input gzipped FASTA file")
        parser.add_argument("output_dir", help="Output directory for chunks")
        parser.add_argument(
            "-n", "--num-chunks", type=int, default=100, help="Number of chunks"
        )
        parser.add_argument(
            "-l", "--log", default="/dev/stderr", help="Log file path"
        )
        args = parser.parse_args()

        output_files = [
            f"{args.output_dir}/chunk_{i}.fasta" for i in range(1, args.num_chunks + 1)
        ]
        split_fasta(args.input, output_files, args.log)
