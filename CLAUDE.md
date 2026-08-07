# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

This is the source repository for *Introduction to Computational Biology*, a course taught by Michael Love. Materials consist of R Markdown (`.Rmd`) and Quarto (`.qmd`) lecture notes, R scripts, and Stan models, organized by topic into module directories.

Course website: http://biodatascience.github.io/compbio

## Rendering Documents

```r
# Render a single .Rmd file
rmarkdown::render("module/file.Rmd")

# Render a .qmd file
quarto::quarto_render("module/file.qmd")
```

Rendered HTML files and cache directories (`*_cache/`, `*_files/`) are committed alongside source files.

## Module Structure

Each module directory follows a consistent pattern:
- Lecture notes: `topic.Rmd` or `topic.qmd` → rendered to `topic.html`
- Homework: `topic_HW.Rmd` → rendered to `topic_HW.html`
- Supporting R scripts for standalone computations (not knitr chunks)

| Directory | Topic |
|-----------|-------|
| `bioc/`   | Bioconductor objects, GenomicRanges, annotation, string manipulation |
| `eda/`    | Exploratory data analysis with dplyr and ggplot2 |
| `model/`  | EM algorithm, mixture models, motif finding |
| `hier/`   | Hierarchical models for variance (limma-style shrinkage, Stan) |
| `dist/`   | Distance metrics, batch effects, variance stabilization |
| `multiple/` | Multiple testing: IDR, local FDR, multtest |
| `markov/` | Hidden Markov Models |
| `net/`    | Network analysis |
| `signal/` | Signal processing |
| `github/` | Git/GitHub workflow assignments |

## Key R Packages

- **Bioconductor**: `GenomicRanges`, `SummarizedExperiment`, `BSgenome`, `AnnotationHub`, `DESeq2`, `limma`
- **Tidyverse**: `dplyr`, `ggplot2`, `readr`
- **Stan**: `rstan` (see `hier/simple.stan` for a hierarchical normal model)
- **Genomics utilities**: `rtracklayer`, `Biostrings`

## Rmd → Qmd Migration

As of 2026, migration to Quarto (`.qmd`) is in progress. Current status (11 `.Rmd`, 17 `.qmd`):

- **Already migrated to `.qmd`**: `bioc/` (4 files), `dist/` (5 files), `model/` (4 files), `hier/` (2 files), `eda/brain_RNA.qmd`, `eda/EDA.qmd`
- **Still `.Rmd`**: `markov/`, `multiple/`, `net/`, `eda/NAs.Rmd`, and all `_HW.Rmd` homework files

When migrating a file, rename `.Rmd` → `.qmd` and update the YAML header (replace `output: html_document` with `format: html`).

## Depositing HTML to gh-pages

Rendered HTML files are published via a parallel sibling repo at `../compbio`, which is checked out on the `gh-pages` branch. The directory structure mirrors this repo (e.g. `hier/` → `../compbio/hier/`).

To deposit HTML after rendering:

```bash
cp module/file.html ../compbio/module/file.html
```

For example, after rendering all files in `model/` and `hier/`:

```bash
cp model/*.html ../compbio/model/
cp hier/*.html ../compbio/hier/
```

Then commit and push in `../compbio` to publish.

## Architecture Notes

- The repo is course materials, not a package — there is no `DESCRIPTION`, `NAMESPACE`, or test suite.
- `bioc/my_scale_counts.R` is a reusable utility for normalizing RNA-seq counts by mapped read depth; it is sourced by other documents.
- `hier/simple.stan` implements a hierarchical normal shrinkage prior; companion R scripts (`simple_stan.R`, `simple_hierarchical.R`, `hierarchical_iterate.R`) demonstrate manual EM-style iteration vs. Stan.
- Data files (`.bed`, `.tsv`, `.csv.gz`) live alongside the `.Rmd` that uses them; no separate data directory.
