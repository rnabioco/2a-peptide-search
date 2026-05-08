# Legacy Directory

This directory contains archived historical results, scripts, and models from earlier phases of the 2A peptide discovery project. These files are preserved for reference but are **not actively used** by the current Snakemake pipeline.

## Directory Structure

### `db-searches-historical/`
Historical HMM search results from early analysis phases:
- `*.sto.gz` - Stockholm alignment files from database searches
- `merge-alignments.sh` - Script used to combine alignments
- `README.md` - Original documentation for these results

These results predate the automated pipeline and used different filtering criteria.

### `curated-models-original/`
Original manually curated HMM models:
- `2A-class-1.hmm`, `2A-class-2.hmm` - Original HMM models
- `2A-class-1.sto.gz`, `2A-class-2.sto.gz` - Source alignments
- `build-models.sh` - Script used to build models

These models served as the foundation for the current seed alignments in `resources/seed-alignments/`.

### `old-scripts/src/`
Superseded scripts from before pipeline automation:
- `expand-hits.py` - Now replaced by `workflow/scripts/expand_hits.py`
- `hmmsearch.sh`, `hmmsearch.mgy.sh` - Manual search scripts (now automated in `workflow/rules/search.smk`)

## Cleanup Policy

These files are retained for:
1. **Reproducibility** - Understanding how original models were built
2. **Reference** - Comparing current results to historical baselines
3. **Recovery** - Source material if seed alignments need reconstruction

Files here should **not** be deleted without team discussion, as they represent the project's analytical history.

## Migration Notes

The current pipeline (`workflow/`) supersedes all legacy scripts. Key differences:
- Automated workflow via Snakemake (vs. manual shell scripts)
- Pixi for dependency management (vs. ad-hoc conda environments)
- Structured output organization (vs. flat file layout)
- Automated curation with manual override option (vs. fully manual)
