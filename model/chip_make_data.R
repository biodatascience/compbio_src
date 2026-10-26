# Make the ChIP-seq peak sequences used in the motif finding class activity.
# Peaks are ENCODE (mostly Snyder lab, Yale) TFBS narrowPeak files from AnnotationHub (hg19).
# For each TF we keep the top peaks by signal, take 200 bp centered on the peak
# summit, and lift over to hg38. Control sequences are the same windows shifted
# 2 kb to the right. Output: chip_seqs.rds, a list with one element per TF,
# each a list of two DNAStringSets (peaks, control).
library(AnnotationHub)
library(plyranges)
library(rtracklayer)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg38)
library(here)

ah <- AnnotationHub()
q <- query(ah, c("Tfbs", "StdPk", "hg19"))

tfs <- c(
  cJun  = "wgEncodeSydhTfbsK562CjunStdPk.narrowPeak.gz",
  NFYA  = "wgEncodeSydhTfbsK562NfyaStdPk.narrowPeak.gz",
  cMyc  = "wgEncodeSydhTfbsK562CmycStdPk.narrowPeak.gz",
  YY1   = "wgEncodeSydhTfbsGm12878Yy1StdPk.narrowPeak.gz",
  STAT1 = "wgEncodeSydhTfbsK562Stat1Ifng30StdPk.narrowPeak.gz",
  CEBPB = "wgEncodeSydhTfbsHepg2CebpbForsklnStdPk.narrowPeak.gz",
  Pol2  = "wgEncodeSydhTfbsK562Pol2StdPk.narrowPeak.gz"
)
ntop <- 500
win_width <- 200
chain_file <- tempfile(fileext=".chain")
R.utils::gunzip(here("bioc", "hg19ToHg38.over.chain.gz"), destname=chain_file, remove=FALSE)
chain <- import.chain(chain_file)

get_chip_and_ctl_seqs <- function(title) {
  peaks <- q[q$title == title][[1]]
  # narrow to the summit, then take a window centered on it
  peaks <- peaks |>
    keepStandardChromosomes(pruning.mode="coarse") |>
    mutate(start=start + peak, width=1) |>
    anchor_center() |>
    mutate(width=win_width)
  # lift over all peaks and drop those that did not map to a single window,
  # and only then take the top peaks, so that we end up with ntop peaks
  peaks <- unlist(liftOver(peaks, chain))
  peaks <- peaks |>
    filter(width == win_width) |>
    keepStandardChromosomes(pruning.mode="coarse") |>
    slice_max(signalValue, n=ntop, with_ties=FALSE) |>
    arrange(desc(signalValue))
  control <- peaks |> shift_right(2000)
  list(peaks=getSeq(Hsapiens, peaks), control=getSeq(Hsapiens, control))
}

chip_seqs <- lapply(tfs, get_chip_and_ctl_seqs)
sapply(chip_seqs, function(x) length(x$peaks))
saveRDS(chip_seqs, file=here("model", "chip_seqs.rds"))
