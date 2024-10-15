library(curatedBladderData)
library(affy)
data(GSE13507_eset)
e <- GSE13507_eset
e <- e[,which(e$summarystage == "superficial")]
e <- e[,-which(colnames(e)=="GSM340606")]
library(matrixStats)
rm <- rowMeans(exprs(e))
rv <- rowVars(exprs(e))
library(limma)

evaluate <- function(n) {
  sample.idx <- sample(ncol(e), n)
  e.sub <- e[,sample.idx]
  rv.sub <- rowVars(exprs(e.sub))
  design <- model.matrix(~1, data.frame(row.names=1:n))
  fit <- lmFit(exprs(e.sub), design)
  fit <- eBayes(fit)
  true <- sqrt(rv) # the SD over the whole dataset
  samp.sd <- sqrt(rv.sub)
  ebayes.sd <- sqrt(fit$s2.post)
  rmse <- function(x) sqrt(mean((x - true)^2))
  mad <- function(x) median(abs(x - true))
  c(sd = rmse(samp.sd),
    eb = rmse(ebayes.sd),
    pr = rmse(fit$s2.prior),
    sd_mad = mad(samp.sd),
    eb_mad = mad(ebayes.sd),
    pr_mad = mad(fit$s2.prior)
    )
}

boxplot(t(replicate(50, evaluate(n=3))))
boxplot(t(replicate(50, evaluate(n=10))))
boxplot(t(replicate(50, evaluate(n=50))))
