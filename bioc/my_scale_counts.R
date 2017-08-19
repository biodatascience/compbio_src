my_scale_counts <- function(rse_gene, round=TRUE) {
  cts <- assays(rse_gene)$counts
  # mapped_read_count includes the count for both reads of a pair
  paired <- ifelse(colData(rse_gene)$paired_end, 2, 1)
  x <- (colData(rse_gene)$mapped_read_count / paired) / colSums(cts)
  assays(rse_gene)$counts <- t(t(assays(rse_gene)$counts) * x)
  if (round) {
    assays(rse_gene)$counts <- round(assays(rse_gene)$counts)
  }
  rse_gene
}
