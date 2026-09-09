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

HTML files are gitignored in this repo (see `.gitignore`). Cache directories (`*_cache/`, `*_files/`) are also gitignored. Only source `.qmd`/`.Rmd` files are committed here.

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

As of 2026, all lecture files have been migrated to Quarto (`.qmd`). Current status (5 `.Rmd`, 23 `.qmd`):

- **All lecture files**: migrated to `.qmd` across all modules
- **Still `.Rmd`**: only `_HW.Rmd` homework files (intentionally kept as `.Rmd` with `execute: eval: false` to prevent standalone execution)

When migrating a file, rename `.Rmd` → `.qmd` and update the YAML header:
- Lecture files: replace `output: html_document` with `format:\n  html:\n    embed-resources: true`
- HW files: also add `execute:\n  eval: false` (they depend on cross-document context)

Known issues when rendering:
- After loading `plotgardener`, `TxDb.*`, or `org.Hs.eg.db`, `keepStandardChromosomes()` may become unavailable. Use `GenomeInfoDb::keepStandardChromosomes()` and `GenomeInfoDb::seqlevelsStyle()` instead.
- `multiple/multtest.qmd` uses `recount3` to load ERP020977 (HipSci macrophage RNA-seq, 317 naive samples); requires internet access at render time.
- `library(Seqinfo)` in `bioc/objects.qmd` is correct: `Seqinfo` is a standalone Bioconductor package (separate from `GenomeInfoDb`) that exposes `Seqinfo(genome=)` for fetching chromosome lengths by genome build.

## Depositing HTML to gh-pages

Rendered HTML files are published via a parallel sibling repo at `../compbio` (this repo is `compbio_src`; the gh-pages repo is the sibling `compbio`), which is checked out on the `gh-pages` branch. The directory structure mirrors this repo (e.g. `hier/` → `../compbio/hier/`).

To deposit HTML after rendering, use `\cp` to bypass the `cp -i` alias:

```bash
\cp module/file.html ../compbio/module/file.html
```

For example, after rendering all files in `model/` and `hier/`:

```bash
\cp model/*.html ../compbio/model/
\cp hier/*.html ../compbio/hier/
```

**Do not deposit `_HW.html` files** — homework files are not published to the course website.

The user handles all commits and pushes in `../compbio` — do not do this automatically.

## plyranges / tidyomics Style

When working with GRanges objects, prefer `plyranges` idioms over base GenomicRanges:

- Use `keepStandardChromosomes(pruning.mode="coarse")` rather than `filter(seqnames %in% std_chroms)`.
- Call `unstrand()` before `reduce_ranges()` when the goal is strand-agnostic region definitions (e.g. "does a variant fall in an exon"). Without it, `reduce_ranges()` reduces per strand and overlapping + / - ranges are never merged.
- Use `gaps(..., ignore.strand=TRUE)` to compute intergenic complements; without it, `gaps()` returns spurious full-chromosome ranges on `"+"` and `"-"` because those strands have no coverage in an unstranded input.
- Use `mutate()` to modify range components (`start`, `width`) as well as metadata columns, and to factorize columns in the pipe rather than via `$<-` after the fact.
- Use `slice_sample(n=n, weight_by=col)` in place of `sample.int(..., prob=...)` + subsetting.
- Use `bind_ranges(..., .id="col")` instead of `c()` when combining labeled groups; the argument names become the level values of the id column automatically.
- Use `select(-col)` to drop temporary columns (e.g. intermediate weight columns).
- Use `unique(col)` for factor levels derived from the data rather than hardcoding a character vector.

## Writing Style

- Do not use bold (`**text**`) in prose
- Do not use em dashes (`—`); use a comma, parentheses, or rewrite the sentence instead

## Architecture Notes

- The repo is course materials, not a package — there is no `DESCRIPTION`, `NAMESPACE`, or test suite.
- `bioc/my_scale_counts.R` is a reusable utility for normalizing RNA-seq counts by mapped read depth; it is sourced by other documents.
- `hier/simple.stan` implements a hierarchical normal shrinkage prior; companion R scripts (`simple_stan.R`, `simple_hierarchical.R`, `hierarchical_iterate.R`) demonstrate manual EM-style iteration vs. Stan.
- Data files (`.bed`, `.tsv`, `.csv.gz`) live alongside the `.Rmd` that uses them; no separate data directory.
