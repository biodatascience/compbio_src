library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
suppressMessages({
  g <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
  g <- keepStandardChromosomes(g, pruning.mode="coarse")
})
pro <- promoters(g, upstream=1500, downstream=1500)
seq <- getSeq(Hsapiens, pro)
gc <- as.vector( letterFrequency(seq, letters="GC") / 3000 )
cpg <- as.vector( vcountPDict(PDict("CG"), seq) / 3000 )

library(plyranges)
g <- g |>
  mutate(norm_cpg = cpg / (gc/2)^2)

library(org.Hs.eg.db)
g$refseq <- mapIds(org.Hs.eg.db, g$gene_id, "REFSEQ", "ENTREZID")
load("~/Downloads/Housekeeping_GenesHuman.RData")
head(Housekeeping_Genes)
g <- g |>
  mutate(housekeeping = refseq %in% Housekeeping_Genes$Refseq)
table(g$housekeeping)

library(ggplot2)
library(tibble)
tab <- g |>
  plyranges::select(housekeeping, norm_cpg, .drop_ranges=TRUE) |>
  as_tibble() |>
  filter(norm_cpg < 1)

f <- prop.table(table(g$housekeeping))["TRUE"]
tab |>
  ggplot(aes(norm_cpg)) + 
  geom_histogram(bins = 100, aes(y=..density..)) + 
  geom_density(data = tab |> filter(housekeeping), 
               aes(y = ..density.. * f),
               col="magenta")
