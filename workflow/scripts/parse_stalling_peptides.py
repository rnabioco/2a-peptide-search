#!/usr/bin/env python3
"""
Parse stalling peptides from TSV and create FASTA file.
Extracts peptide regions ending at the stalling motif.
"""

import pandas as pd
import sys
from pathlib import Path


def extract_peptide_window(seq_aa, position, motif, upstream=30):
    """
    Extract peptide sequence ending at the motif position.

    Args:
        seq_aa: Full protein sequence
        position: Position string like "68..71"
        motif: Expected motif sequence (for validation)
        upstream: Number of residues upstream of motif start to include

    Returns:
        Extracted peptide sequence ending with motif
    """
    if pd.isna(seq_aa) or pd.isna(position):
        return None

    # Parse position (format: "68..71", 1-indexed)
    try:
        start_pos, end_pos = map(int, position.split(".."))
    except (ValueError, AttributeError):
        return None

    # Convert to 0-indexed Python coordinates
    motif_start = start_pos - 1
    motif_end = end_pos  # Python slicing is exclusive at end

    # Verify the motif matches what's expected
    extracted_motif = seq_aa[motif_start:motif_end]
    if extracted_motif.upper() != motif.upper():
        print(f"Warning: Motif mismatch. Expected {motif}, got {extracted_motif} at position {position}",
              file=sys.stderr)

    # Extract window: upstream context + motif
    window_start = max(0, motif_start - upstream)
    window_end = motif_end

    peptide = seq_aa[window_start:window_end]

    return peptide


def main():
    input_tsv = Path("resources/stalling-peptides/pmid-38565864.tsv")
    output_fasta = Path("resources/stalling-peptides/known_stalling_peptides.fasta")

    # Read TSV (first column is '#' as index column name)
    df = pd.read_csv(input_tsv, sep="\t")
    # Rename the first column if it's '#'
    if df.columns[0] == '#':
        df = df.rename(columns={'#': 'id'})

    # Filter to rows with sequence data
    df_with_seq = df[df["seq_aa"].notna() & df["position"].notna()].copy()

    print(f"Processing {len(df_with_seq)} sequences with position data...", file=sys.stderr)

    # Extract peptides
    peptides = []
    for idx, row in df_with_seq.iterrows():
        peptide = extract_peptide_window(
            row["seq_aa"],
            row["position"],
            row["query_motif"],
            upstream=30
        )

        if peptide:
            # Create FASTA header with metadata
            header = (
                f">{row['query_motif']}|{row['organism'].replace(' ', '_')}|"
                f"{row['gene']}|{row['locus_tag']}|pos_{row['position']}|"
                f"invitro_{row['in vitro analysis']}"
            )

            peptides.append((header, peptide))

    # Write FASTA
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    with open(output_fasta, "w") as f:
        for header, seq in peptides:
            f.write(f"{header}\n")
            # Write sequence in 60-character lines
            for i in range(0, len(seq), 60):
                f.write(f"{seq[i:i+60]}\n")

    print(f"Wrote {len(peptides)} peptide sequences to {output_fasta}", file=sys.stderr)
    print(f"Sequences range from {min(len(s) for _, s in peptides)} to "
          f"{max(len(s) for _, s in peptides)} residues", file=sys.stderr)


if __name__ == "__main__":
    main()
