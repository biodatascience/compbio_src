library(here)
library(DESeq2)
library(matrixStats)
library(ggplot2)
load(here("bioc","geuvadis.rda"))

meanSdPlot <- function(x) {
  df <- data.frame(mean = rowMeans2(x), sd = rowSds(x))
  ggplot(df, aes(mean, sd)) +
    geom_hex() +
    geom_smooth()
}

ntd <- normTransform(dds)
meanSdPlot(assay(ntd))

x <- assay(ntd)[,1]
y <- assay(ntd)[,2]
plot(.5*(x + y), y - x,
     cex=.5, col=rgb(0,0,0,.1), pch=20)
abline(h=0, col="red", lwd=3)

vsd <- vst(dds)
