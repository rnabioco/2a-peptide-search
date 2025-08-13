# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a bioinformatics research project focused on systematic discovery and analysis of 2A peptides using hidden Markov models (HMMs). 2A peptides are short (~15 residue) sequences that cause ribosomal skipping, primarily found in viral proteins.

The project uses computational methods to:
- Build and refine HMM models for two classes of 2A peptides
- Search large protein databases (UniProt, UniParc, MGnify, IMGVR) for new instances
- Analyze and curate sequence alignments
- Generate visualizations and reports

## Key Commands

### Building HMM Models
```bash
# From curated-models/ directory
./build-models.sh
```
This creates `2A-class-1.hmm` and `2A-class-2.hmm` profile models from Stockholm format alignments.

### Running Database Searches
```bash
# Basic hmmsearch against a protein database
./src/hmmsearch.sh

# Search MGnify database (requires cluster submission)
./src/hmmsearch.mgy.sh

# Search with specific parameters
hmmsearch --cpu 12 -o results.out -A results.sto model.hmm database.fa.gz
```

### Processing Results
```bash
# Extract sequences around hits from database
python src/expand-hits.py alignment.sto database.fa.gz

# Merge multiple Stockholm alignments
cd db-searches/
./merge-alignments.sh
```

### Python Dependencies
The Python scripts require:
- BioPython (`Bio.SeqIO`, `Bio.AlignIO`)
- pyfaidx (for FASTA indexing)
- click (for CLI interfaces)

## Architecture

### Directory Structure
- `curated-models/` - HMM models and seed alignments for class 1 and 2 2A peptides
- `db-searches/` - Results from searching major protein databases, including merged alignments
- `results/YYYY-MM-DD/` - Time-stamped analysis results with notes and intermediate files
- `src/` - Core scripts for database searches and sequence processing
- `data/ref/` - Reference databases (downloaded separately)
- `img/` - Sequence logos and visualizations

### Key Files
- `curated-models/2A-class-1.hmm` and `2A-class-2.hmm` - Main HMM models
- `db-searches/combined-class-*.sto.gz` - Merged alignments across all database searches
- `src/expand-hits.py` - Extracts extended sequences around HMM hits
- `src/hmmsearch.sh` - Template for database searches

### Analysis Workflow
1. **Model Building**: Start with curated seed alignments in Stockholm format
2. **Database Search**: Use `hmmsearch` to find homologs in protein databases
3. **Results Processing**: Extract hit sequences with flanking regions using `expand-hits.py`
4. **Alignment Merging**: Combine results from multiple databases using `merge-alignments.sh`
5. **Curation**: Manual review and refinement of alignments for next iteration

### File Formats
- **Stockholm (.sto)**: Multiple sequence alignments, often gzipped
- **HMM**: HMMER profile models
- **FASTA**: Protein sequence databases (typically gzipped)
- **TSV**: Tabular results from hmmsearch and analysis scripts

## R Project Integration
This is an R project (`.Rproj` file present) with Quarto analysis documents in `results/` subdirectories for statistical analysis and visualization.