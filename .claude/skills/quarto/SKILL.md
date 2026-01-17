---
name: quarto
description: Create and edit Quarto documents (.qmd files) following project standards. Use when asked to create a new Quarto document, analysis notebook, or report, or when editing existing .qmd files.
allowed-tools: Read, Write, Edit, Glob, Grep
---

# Quarto Document Creation

## Overview
This skill creates Quarto documents (.qmd files) following project-specific standards for bioinformatics analysis. This project primarily uses R with tidyverse for data analysis.

## Chunk Options Format
Always use hash-pipe style for ALL chunk options, including labels:

```r
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

## R Package Loading
Standard imports for this project:

```r
#| label: setup
#| echo: false
#| message: false

library(tidyverse)
library(here)
library(fs)
library(ggseqlogo)  # For sequence logos
```

## Directory Configuration
Use environment variables for cluster/local flexibility:

```r
#| label: setup-directories

results_dir <- Sys.getenv("RESULTS_DIR", "../scratch/results")
data_dir <- Sys.getenv("RESOURCES_DIR", "../scratch/data")
```

## Plotting with ggplot2
This project uses **ggplot2** from the tidyverse for all visualizations:

```r
#| label: fig-example
#| fig-cap: "Example plot"

tbl <- tibble(x = 1:3, y = 4:6)

ggplot(tbl, aes(x, y)) +
    geom_point() +
    theme_minimal() +
    labs(x = "X axis", y = "Y axis")
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

```r
#| label: fig-evalue-distribution
#| fig-cap: "Distribution of E-values from HMM search"
#| fig-width: 8
#| fig-height: 5

tbl <- results_df |>
    filter(evalue < 1e-3)

ggplot(tbl, aes(x = evalue)) +
    geom_histogram(bins = 50) +
    scale_x_log10() +
    theme_minimal() +
    labs(x = "E-value", y = "Count")
```

## Tables with gt
For simple tables, use tibbles directly. For publication tables, use gt:

```r
#| label: tbl-top-hits
#| tbl-cap: "Top 10 HMM hits by E-value"

top_hits <- results_df |>
    slice_min(evalue, n = 10) |>
    select(accession, description, evalue, score)

top_hits
```

For polished tables, use the gt package:

```r
#| label: tbl-formatted-hits

library(gt)

top_hits |>
    gt() |>
    fmt_scientific(columns = evalue) |>
    fmt_number(columns = score, decimals = 1)
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

```{r}
#| label: session-info

sessionInfo()
```
:::
````

## Project-Specific Notes
- Analysis documents go in `analysis/` directory
- Use `ggplot2` for plotting (tidyverse)
- Use `ggseqlogo` for sequence logo visualization
- Results go to `scratch/results/`
- Render with: `cd analysis && quarto render`
