---
name: quarto
description: Create and edit Quarto documents (.qmd files) following project standards. Use when asked to create a new Quarto document, analysis notebook, or report, or when editing existing .qmd files.
allowed-tools: Read, Write, Edit, Glob, Grep
---

# Quarto Document Creation

## Overview
This skill creates Quarto documents (.qmd files) following project-specific standards for bioinformatics analysis. This project primarily uses Python for data analysis.

## Chunk Options Format
Always use hash-pipe style for ALL chunk options, including labels:

```python
#| label: load-data
#| echo: false
#| warning: false
```

## Chunk Label Conventions
- Always include informative labels describing what the chunk does
- Use the `#| label:` hash-pipe option
- Use kebab-case: `load-results`, `process-alignments`
- Prefix conventions:
  - Figures: `fig-` (e.g., `fig-evalue-dist`, `fig-conservation`)
  - Tables: `tbl-` (e.g., `tbl-summary-stats`, `tbl-top-hits`)
  - Setup chunks: `setup-` or descriptive name (e.g., `load-packages`)

## Standard Document Structure

```yaml
---
title: "Analysis Title"
subtitle: "Brief description"
date: today
execute:
  echo: false
  warning: false
  message: false
---
```

## Python Package Loading
Standard imports for this project:

```python
#| label: setup
#| echo: false

import pandas as pd
import numpy as np
from pathlib import Path
import plotnine as p9
from Bio import SeqIO, AlignIO
```

## Directory Configuration
Use environment variables for cluster/local flexibility:

```python
#| label: setup-directories

import os
results_dir = Path(os.environ.get("RESULTS_DIR", "../scratch/results"))
data_dir = Path(os.environ.get("RESOURCES_DIR", "../scratch/data"))
```

## Plotting with plotnine
This project uses **plotnine** (ggplot2 for Python), NOT matplotlib or seaborn:

```python
#| label: fig-example
#| fig-cap: "Example plot"

from plotnine import *

tbl = pd.DataFrame({'x': [1,2,3], 'y': [4,5,6]})

(ggplot(tbl, aes('x', 'y'))
 + geom_point()
 + theme_minimal()
 + labs(x='X axis', y='Y axis'))
```

## Common Chunk Options Reference

| Option | Values | Description |
|--------|--------|-------------|
| `echo` | true/false | Show code in output |
| `eval` | true/false | Execute the code |
| `warning` | true/false | Show warnings |
| `fig-width` | number | Figure width in inches |
| `fig-height` | number | Figure height in inches |
| `fig-cap` | string | Figure caption |
| `tbl-cap` | string | Table caption |

## Figure Dimensions
Use the golden ratio (≈1.618) for figure dimensions. Default to 8×5 inches:

```python
#| fig-width: 8
#| fig-height: 5
```

Common sizes:
- Standard: 8×5 (default)
- Wide: 10×6
- Square: 6×6 (for heatmaps)

## Figure Chunks Example

```python
#| label: fig-evalue-distribution
#| fig-cap: "Distribution of E-values from HMM search"
#| fig-width: 8
#| fig-height: 5

tbl = results_df.query('evalue < 1e-3')

(ggplot(tbl, aes(x='evalue'))
 + geom_histogram(bins=50)
 + scale_x_log10()
 + theme_minimal()
 + labs(x='E-value', y='Count', title='HMM Search Results'))
```

## Tables with pandas/gt
For simple tables, use pandas with styling. For publication tables, consider gt:

```python
#| label: tbl-top-hits
#| tbl-cap: "Top 10 HMM hits by E-value"

top_hits = results_df.nsmallest(10, 'evalue')[['accession', 'description', 'evalue', 'score']]
top_hits
```

## Mermaid Diagrams
Use mermaid for workflow diagrams:

````markdown
```{mermaid}
%%| label: fig-workflow
%%| fig-cap: "Analysis workflow"
flowchart LR
    A[Input] --> B[Process] --> C[Output]
```
````

## Summary Callouts
Wrap key findings in callouts:

```markdown
::: {.callout-note}
## Summary

- Key finding 1
- Key finding 2
:::
```

Callout types:
- `.callout-note` - General summaries
- `.callout-tip` - Recommendations
- `.callout-warning` - Caveats
- `.callout-important` - Critical information

## Session Info
Include session info in collapsible callout:

````markdown
::: {.callout-note collapse="true"}
## Session Info

```{python}
#| label: session-info

import session_info
session_info.show()
```
:::
````

## Project-Specific Notes
- Analysis documents go in `analysis/` directory
- Use `plotnine` for plotting (NOT matplotlib/seaborn)
- Use `Bio` (Biopython) for sequence handling
- Results go to `scratch/results/`
- Render with: `cd analysis && quarto render`
