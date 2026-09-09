# bioc/tidy_syntax_make_data.R
# Task 1: simulated variant locations drawn from promoter / intergenic /
#         intronic / exonic regions  (plyranges / tidyomics style)
# Task 2: normalized RNA-seq SummarizedExperiment (airway epithelial cells,
#         SRP046226 from recount2, all processing from bioc/objects.qmd)
# Task 3: GRangesList of paired isoforms with skipped-exon events, drawn
#         from two-isoform genes in EnsDb.Hsapiens.v86

suppressPackageStartupMessages({
  library(EnsDb.Hsapiens.v86)
  library(plyranges)
  library(GenomeInfoDb)
  library(here)
  library(SummarizedExperiment)
  library(DESeq2)
  library(stringr)
  library(rtracklayer)
  library(Seqinfo)
  library(splicelogic)
})

# ─────────────────────────────────────────────────────────────────────────────
# Task 1: variant locations  (plyranges / tidyomics style)
# ─────────────────────────────────────────────────────────────────────────────

edb <- EnsDb.Hsapiens.v86

# Stranded gene ranges, standard chromosomes only
g <- genes(edb) |>
  GenomeInfoDb::keepStandardChromosomes(pruning.mode = "coarse")

# Promoters: 2 kb upstream of the TSS (strand-aware), clipped to chr bounds
prom_regions <- g |>
  promoters(upstream = 2000, downstream = 0) |>
  trim() |>
  unstrand() |>
  reduce_ranges()

# Exonic regions (merged, unstranded)
# unstrand() before reduce_ranges() so overlapping + and - exons are merged
exon_regions <- exons(edb) |>
  GenomeInfoDb::keepStandardChromosomes(pruning.mode = "coarse") |>
  unstrand() |>
  reduce_ranges()

# Gene bodies (merged, unstranded); introns = gene body minus exons
gene_bodies <- g |>
  unstrand() |>
  reduce_ranges()

intron_regions <- gene_bodies |> 
  setdiff(exon_regions)

# Intergenic: complement of (gene bodies union promoters) within chr bounds
genic_union <- c(gene_bodies, prom_regions) |>
  reduce_ranges()

# gaps() computes the complement for each distinct (seqname, strand) pair.
# ignore.strand=TRUE avoids spurious full-chromosome gaps on "+" and "-"
# that would otherwise appear because genic_union has no ranges on those strands.
intergenic_regions <- gaps(genic_union, ignore.strand = TRUE)

# Sample n point variants uniformly at random from within a set of ranges
sample_from_regions <- function(regions, n) {
  stopifnot(n <= length(regions))
  regions |>
    mutate(prop = width/mean(width)) |>
    slice_sample(n = n, weight_by = prop) |>
    mutate(start = start + as.integer(runif(n, 0, width)),
           width = 1L) |>
    unstrand() |>
    select(-prop)
}

set.seed(5)
n_per_class <- 500

variants <- bind_ranges(
  promoter   = sample_from_regions(prom_regions,       n_per_class),
  intergenic = sample_from_regions(intergenic_regions, n_per_class),
  intronic   = sample_from_regions(intron_regions,     n_per_class),
  exonic     = sample_from_regions(exon_regions,       n_per_class),
  .id = "type"
) |>
  mutate(type = factor(type, levels = unique(type))) |>
  sort()

saveRDS(variants, here("bioc", "variants.rds"))
message("Task 1 data: ", length(variants), " variants -> bioc/variants.rds")

# Version without true labels (for use as a prediction task)
variants |>
  select(-type) |>
  saveRDS(here("bioc", "variants_no_labels.rds"))
message("  unlabeled copy -> bioc/variants_no_labels.rds")

# ─────────────────────────────────────────────────────────────────────────────
# Task 2: normalize the airway epithelial SummarizedExperiment
# ─────────────────────────────────────────────────────────────────────────────

url  <- "http://duffel.rail.bio/recount/SRP046226/rse_gene.Rdata"
file <- here("bioc", "asthma.rda")
if (!file.exists(file)) download.file(url, file)
load(file)

source(here("bioc", "my_scale_counts.R"))
rse <- my_scale_counts(rse_gene)

# Parse sample characteristics into condition and treatment columns
rse$condition <- sapply(rse$characteristics, `[`, 3)
rse$treatment  <- sapply(rse$characteristics, `[`, 4)

rse$condition <- rse$condition |>
  str_remove("disease state: ") |>
  str_replace("-", ".") |>
  factor()

rse$treatment <- rse$treatment |>
  str_remove("treatment: ") |>
  factor(levels=c("Vehicle","HRV16"))

# Attach hg38 chromosome lengths to seqinfo
si_hg38 <- Seqinfo::Seqinfo(genome = "hg38")
seqinfo(rse) <- si_hg38[seqlevels(rse)]

# Add gene biotype from GENCODE v25 GTF (same annotation used by recount2)
gtf_url  <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_25/gencode.v25.annotation.gtf.gz"
gtf_file <- here("bioc", "gencode.v25.gtf.gz")
if (!file.exists(gtf_file)) download.file(gtf_url, gtf_file)
genes_gtf <- rtracklayer::import(gtf_file, feature.type = "gene")

m <- match(rownames(rse), genes_gtf$gene_id)
mcols(rse)$gene_type <- genes_gtf$gene_type[m]

# Variance-stabilizing transformation
keep <- rowSums(assay(rse) >= 10) >= 12
dds <- DESeqDataSet(rse[keep,], ~ condition + treatment)
dds <- DESeq(dds, test="LRT", reduced = ~ condition)
vsd <- vst(dds, blind=FALSE)
rowData(vsd) <- rowData(vsd)[,c("gene_id","symbol","gene_type","treatment_HRV16_vs_Vehicle","LRTPvalue")]

saveRDS(vsd, here("bioc", "airway_vst_data.rds"))
message("Task 2 data: VST-normalized SE -> bioc/airway_vst_data.rds")

# ─────────────────────────────────────────────────────────────────────────────
# Task 3: compare transcript sets
# ─────────────────────────────────────────────────────────────────────────────

# Find all genes with exactly two isoforms
txdf <- AnnotationDbi::select(edb,
               keys    = keys(edb, "GENEID"),
               columns = c("GENEID", "TXID"),
               keytype = "GENEID")

two_iso_genes <- txdf |>
  dplyr::count(GENEID) |>
  dplyr::filter(n == 2) |>
  dplyr::pull(GENEID)

# Pick genes at random
set.seed(5)
selected_genes <- sample(two_iso_genes, 100)

# GRangesList of exons per transcript (up to 200 transcripts before filtering)
ebt <- exonsBy(edb, by = "tx") |>
  GenomeInfoDb::keepStandardChromosomes(pruning.mode = "coarse")

mcols(ebt)$tx_id   <- names(ebt)
mcols(ebt)$gene_id <- mapIds(edb, names(ebt), "GENEID", "TXID")

ebt_sub <- ebt[mcols(ebt)$gene_id %in% selected_genes]
ebt_sub <- ebt_sub[order(mcols(ebt_sub)$gene_id)]

# Assign direction: +1 to the isoform with fewer exons, -1 to the other.
# Genes where both isoforms have equal exon counts are dropped.
n_exons   <- lengths(ebt_sub)
direction <- setNames(rep(NA_integer_, length(ebt_sub)), mcols(ebt_sub)$tx_id)

for (gid in unique(mcols(ebt_sub)$gene_id)) {
  txs  <- mcols(ebt_sub)$tx_id[mcols(ebt_sub)$gene_id == gid]
  cnts <- n_exons[txs]
  if (cnts[1] == cnts[2]) next
  direction[txs[which.min(cnts)]] <-  1L
  direction[txs[which.max(cnts)]] <- -1L
}

ebt_sub <- ebt_sub[!is.na(direction[mcols(ebt_sub)$tx_id])]

# Propagate tx_id, gene_id, and direction into each element's mcols
exons_gr <- unlist(ebt_sub, use.names = FALSE) |>
  mutate(
    tx_id     = rep(mcols(ebt_sub)$tx_id,   lengths(ebt_sub)),
    gene_id   = rep(mcols(ebt_sub)$gene_id, lengths(ebt_sub)),
    direction = rep(direction[mcols(ebt_sub)$tx_id], lengths(ebt_sub))
  )

se_exons <- exons_gr |>
  preprocess(coef_col="direction") |>
  find_skipped_exons()

se_isoforms <- c(se_exons$event_tx_id, se_exons$tx_id)
ebt_sub <- relist(exons_gr, ebt_sub)
ebt_sub <- ebt_sub[names(ebt_sub) %in% se_isoforms]

saveRDS(ebt_sub, here("bioc", "two_isoform_exons.rds"))
message("Task 3 data: GRangesList of transcript pairs -> bioc/two_isoform_exons.rds")
