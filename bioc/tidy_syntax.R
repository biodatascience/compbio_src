# bioc/tidy_syntax.R
# Tidy / plyranges examples using the three prepared data objects:
#   Task 1: simulated variant locations (GRanges)
#   Task 2: asthmatic and non-asthmatic airway epithelial RNA-seq
#           (SummarizedExperiment, VST-normalized and counts)
#   Task 3: paired isoform exons comparison (GRangesList)

suppressPackageStartupMessages({
  library(here)
  library(SummarizedExperiment)
  library(EnsDb.Hsapiens.v86)
  library(tibble)
  library(dplyr)
  library(plyranges)
  library(plyxp)
})

# a short-hand for this longer function name
keepChroms <- function(x) {
  GenomeInfoDb::keepStandardChromosomes(x, pruning.mode = "coarse")
}

# ─────────────────────────────────────────────────────────────────────────────
# Task 1: locate these variants with respect to genomic annotations
# ─────────────────────────────────────────────────────────────────────────────

variants <- readRDS(here("bioc", "variants.rds"))
table(variants$type)

edb <- EnsDb.Hsapiens.v86

exon_regions <- exons(edb) |>
  keepChroms() |>
  unstrand() |>
  reduce_ranges()

gene_bodies <- genes(edb) |>
  keepChroms() |>
  unstrand() |>
  reduce_ranges()

intron_regions <- gene_bodies |>
  setdiff(exon_regions)

prom_regions <- genes(edb) |>
  keepChroms() |>
  promoters(upstream = 2000, downstream = 0) |>
  trim() |>
  unstrand() |>
  reduce_ranges() |>
  filter_by_non_overlaps(gene_bodies)

intergenic_regions <- c(gene_bodies, prom_regions) |>
  reduce_ranges() |>
  gaps(ignore.strand = TRUE)

findOverlaps(variants, exon_regions)
table(variants[variants %over% exon_regions]$type)

# tidy syntax with plyranges

pred_var <- variants %>%
  mutate(predicted_type = case_when(
    count_overlaps(., exon_regions)       > 0 ~ "exonic",
    count_overlaps(., intron_regions)     > 0 ~ "intronic",
    count_overlaps(., prom_regions)       > 0 ~ "promoter",
    TRUE ~ "intergenic"
  ))

pred_var |>
  select(type, predicted_type, .drop_ranges=TRUE) |>
  as_tibble() |>
  count(type, predicted_type)

# ─────────────────────────────────────────────────────────────────────────────
# Task 2: compute CPM and TPM from the counts matrix
# ─────────────────────────────────────────────────────────────────────────────

vsd <- readRDS(here("bioc", "asthma_vst_data_counts.rds"))

plotPCA(vsd, "treatment")
rowRanges(vsd)$gene_kb <- width(rowRanges(vsd)) / 1e3

xp <- vsd |> new_plyxp()

counts <- assay(vsd, "counts")

counts_per_kb <- counts / rowData(vsd)$gene_kb
assay(vsd, "counts_per_kb") <- counts_per_kb

cpm <- sweep(counts, 2, colSums(counts) / 1e6, FUN = "/")
assay(vsd, "cpm") <- cpm

tpm <- sweep(counts_per_kb, 2, colSums(counts_per_kb) / 1e6, FUN = "/")
assay(vsd, "tpm") <- tpm

keep <- rowData(vsd)$LRTPvalue < 1e-3 & rowData(vsd)$treatment_HRV16_vs_Vehicle > 0
upreg_expr_profile <- rowMeans(assay(vsd[keep, ], "vst"))

# tidy syntax with plyxp
# plyxp keeps all computation inside the object: no temporary matrices floating
# in the workspace, column summaries via cols() are scoped to the pipe and
# never stored separately, and results land directly in assays(xp) without an
# explicit assay<- assignment.

xp <- xp |>
  mutate(
    cols(col_sum = colSums(counts)),
    cpm = counts / .cols$col_sum * 1e6
  )

assay(vsd, "cpm")[1:5,1:5]
assay(xp, "cpm")[1:5,1:5]

xp <- xp |>
  mutate(
    counts_per_kb = counts / .rows$gene_kb,
    cols(cpk_sum = colSums(counts_per_kb)),
    tpm = counts_per_kb / .cols$cpk_sum * 1e6
  )

assay(vsd, "tpm")[1:5,1:5]
assay(xp, "tpm")[1:5,1:5]



# ─────────────────────────────────────────────────────────────────────────────
# Task 3: paired isoform exons
# ─────────────────────────────────────────────────────────────────────────────

ebt_sub <- readRDS(here("bioc", "two_isoform_exons.rds"))
