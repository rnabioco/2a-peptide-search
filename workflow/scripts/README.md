# Workflow Scripts

Python scripts called by Snakemake rules for data processing and analysis.

## Script Index

| Script | Purpose | Used By |
|--------|---------|---------|
| `analyze_gp_conservation.py` | Analyze conservation patterns within GP motif clusters | `prokaryotic_gp.smk:analyze_cluster_conservation` |
| `annotate_gp_domains.py` | Annotate GP-containing sequences with Pfam domains | (utility) |
| `auto_curate_alignment.py` | Apply automated quality filters to alignments | `refine.smk:auto_curate_alignment` |
| `cluster_gp_motifs.py` | Cluster GP motifs by sequence similarity | (utility) |
| `collect_statistics.py` | Collect hit statistics from tblout files for reporting | `report.smk:collect_statistics` |
| `compare_gp_approaches.py` | Compare seed-based vs GP-discovery results | `prokaryotic_analysis.smk:compare_approaches` |
| `expand_hits.py` | Expand partial hits to full-length sequences | `refine.smk:expand_hits` |
| `extract_gp_motifs.py` | Extract GP-containing sequences with context | `prokaryotic_gp.smk:extract_gp_motifs` |
| `filter_alignment.py` | Filter alignments by quality criteria | `refine.smk:filter_alignment` |
| `filter_interdomain_gp.py` | Filter GP motifs occurring between domains | `prokaryotic_gp.smk:filter_interdomain_gp` |
| `generate_skylign_logos.py` | Generate sequence logos via Skylign API | `logos.smk:generate_skylign_logos` |
| `identify_consensus_patterns.py` | Identify consensus patterns from clusters | `prokaryotic_gp.smk:identify_consensus_patterns` |
| `merge_alignments.py` | Merge multiple Stockholm alignments | `refine.smk:merge_alignments` |
| `parse_domain_annotations.py` | Parse hmmscan domtblout output | `prokaryotic_gp.smk:parse_domain_annotations` |
| `parse_gp_clusters.py` | Parse MMseqs2 clustering results | `prokaryotic_gp.smk:parse_gp_clusters` |
| `parse_stalling_peptides.py` | Parse known stalling peptide sequences | (utility) |
| `sanitize_fasta.py` | Remove invalid IUPAC characters from FASTA | `prokaryotic_download.smk:sanitize_imgvr` |
| `split_peptides_by_motif.py` | Split peptides into files by motif type | `prokaryotic_seeds.smk:split_known_peptides_by_motif` |
| `validate_stalling_peptides.py` | Validate HMMs against known stalling peptides | `prokaryotic_analysis.smk:validate_against_known_peptides` |

## Quarto Reports

| File | Purpose |
|------|---------|
| `prokaryotic_discovery_report.qmd` | Quarto document for prokaryotic discovery report |

## Utilities

The `utils/` subdirectory contains shared modules:
- `__init__.py` - Package initialization
- `logging_config.py` - Standardized logging configuration

## Development Guidelines

### CLI Interface
All scripts use [Click](https://click.palletsprojects.com/) for command-line interfaces:

```python
import click

@click.command()
@click.option('--input', '-i', required=True, help='Input file')
@click.option('--output', '-o', required=True, help='Output file')
def main(input, output):
    """Script description."""
    pass

if __name__ == '__main__':
    main()
```

### Snakemake Integration
Scripts can be called via Snakemake's `script:` directive, which provides access to `snakemake.input`, `snakemake.output`, `snakemake.params`, and `snakemake.log`.

### Plotting
Use **plotnine** (ggplot2 for Python) for all data visualization:

```python
from plotnine import ggplot, aes, geom_point
```

### Dependencies
All dependencies are managed in `pixi.toml` at the project root. Do not create separate conda environment files.
