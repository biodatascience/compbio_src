library(EnsDb.Hsapiens.v86)
library(ensembldb)
txp <- transcripts(EnsDb.Hsapiens.v86)
seqlevelsStyle(txp) <- "UCSC"

library(plyranges)
query <- data.frame(seqnames="chr1", start=100e6, width=2e6) |>
  as_granges()

txp |>
  filter_by_overlaps(query) |>
  write_bed("chr1_100Mb_v86_txps.bed")
