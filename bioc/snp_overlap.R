# read in the SNPs

# GWAS Catalogue 
# https://www.ebi.ac.uk/gwas/studies/GCST90624624
# gzcat GCST90624624.tsv.gz | awk -F'\t' 'NR==1 || $8 < 5e-8' > GCST90624624-filtered.tsv
# gzip GCST90624624-filtered.tsv
library(readr)
library(here)
snps0 <- read_delim(here("bioc","GCST90624624-filtered.tsv.gz"))
library(dplyr)
snps <- snps0 |>
  select(seqnames = chromosome, start=base_pair_location, beta, p_value) |>
  mutate(width = 1, strand = "*", end = start)

library(GenomicRanges)
snps <- as(snps, "GRanges")
# WARNING! DIDN'T CHECK GENOME BUILD
# it is GRCh37...
library(yaml)
gwas_metadata <- read_yaml(here("bioc","GCST90624624.tsv.gz-meta.yaml"))
gwas_metadata[["genome_assembly"]]

# this code deals with changes in genome build...
GenomeInfoDb::seqlevelsStyle(snps) <- "UCSC"
genome(snps) <- "hg19"
library(easylift)
# downloaded from:
# https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
snps <- snps |> easylift::easylift(
  to="hg38", 
  chain=here("bioc","hg19ToHg38.over.chain.gz")
)
GenomeInfoDb::seqlevelsStyle(snps) <- "NCBI"
genome(snps) <- "GRCh38"
# all done...

library(plyranges)
chr17 <- as(data.frame(seqnames="17",start=1,end=83257441), "GRanges")
snps17 <- snps |> filter_by_overlaps(chr17)

if (FALSE) {
  # load genes
  library(ensembldb)
  library(EnsDb.Hsapiens.v86) # this pkg is about 75 Mb
  library(GenomeInfoDb)
  edb <- EnsDb.Hsapiens.v86
  g <- genes(edb) |> keepStandardChromosomes(pruning.mode="coarse")
  ebg <- exonsBy(edb, by="gene") |> keepStandardChromosomes(pruning.mode="coarse")
  save(g, ebg, file=here("bioc","genes.rda"))
}
load(here("bioc","genes.rda"))

ovlps <- snps17 |>
  join_overlap_inner(g)

table(ovlps$symbol)
# PRKCA is the correct answer (according to papers)

res <- snps %>%
  mutate(ovlp = count_overlaps(., g)) |>
  mcols()

table(res$ovlp > 0)

tiles <- tileGenome(
  seqlengths(snps), 
  tilewidth = 1e6, 
  cut.last.tile.in.chrom = TRUE
)

tiles <- tiles |> 
  filter_by_overlaps(snps) %>%
  mutate(tile_id = seq_along(.))

pruned_snps <- snps |>
  join_overlap_inner(tiles) |>
  group_by(tile_id) |>
  slice_min(p_value) |>
  ungroup() |>
  filter(!duplicated(tile_id))

res <- pruned_snps %>%
  mutate(ovlp = count_overlaps(., g)) |>
  mcols()
table(res$ovlp > 0)

pruned_snps |>
  join_nearest(g) |>
  as_tibble() |>
  pull(symbol)

# ebg is a GRangesList so we can use basic GRanges functions
# plyranges only works with GRanges...
fo <- findOverlaps(pruned_snps, ebg)
fo
ebg[subjectHits(fo)]

# can always unlist...
plyranges::find_overlaps(pruned_snps, unlist(ebg))
