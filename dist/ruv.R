library(limma)
library(ruv)

# sample size, number of genes, number of factors
n <- 100
m <- 10e3
k <- 3

# W: n x k
W <- matrix(rnorm(n * k), ncol=k)

# spike and slab model: prob slab = non-zero effects
p_slab <- .25
# alpha: k x m
alpha <- t(replicate(k, {
  ifelse(
    rbinom(m, size=1, prob=p_slab) == 0,
    rep(0,m),
    rnorm(m)
  )
}))

# make the alphas decreasing in effect
alpha <- alpha * 5:1

p_de <- .05 # percent DE genes
# X: n x 2
X <- cbind(rep(1,n), rep(0:1,each=n/2)) # second column is "treatment"
cor(W, X[,2])

# beta: 2 x m
beta <- rbind(rnorm(m),
              c(rep(0, (1 - p_de) * m),
                rnorm(p_de * m, 0, 1)))

# data model
# Y: n x m (samples x genes) 
# NOTE! this is not our normal genomics orientation
Y <- X %*% beta + W %*% alpha + rnorm(n * m, 0, .5)

# standard limma fit
fit <- lmFit(t(Y), X)
efit <- eBayes(fit)
tt <- topTable(efit, coef = 2, number = m, sort.by = "none")

plot(beta[2,], tt$t, cex=.1)
plot(beta[2,], tt$P.Value, cex=.1, log="y")

# RUV fit using known control genes
known_controls <- sample(which(beta[2,] == 0), 1000)
rfit <- RUV2(Y, X[,2], ctl = known_controls, k = 2 * k)

image(cor(rfit$W, W), zlim=c(-1,1), axes = FALSE, 
      xlab="estimated W", ylab="true W",
      col=colorRampPalette(c("red","white","blue"))(99))
axis(1, at = seq(0, 1, length.out = 2*k), labels = 1:(2*k))
axis(2, at = seq(0, 1, length.out = k), labels = 1:k, las=1)

# batch corrected
fit_bc <- lmFit(t(Y), cbind(X, rfit$W[,1:k]))
efit_bc <- eBayes(fit_bc)
tt_bc <- topTable(efit_bc, coef = 2, number = m, sort.by = "none")

panel.cor <- function(x, y, digits = 2,
                      prefix = "", ...) {
  usr <- par("usr"); on.exit(par(usr = usr))
  par(usr = c(0, 1, 0, 1))
  r <- round(cor(x, y), digits)
  txt <- paste0(prefix, r)
  text(0.5, 0.5, txt, cex = 1.5)
}

cols <- ifelse(colSums(abs(alpha)) > 0, rgb(1,0,0,.5), rgb(0,0,0,.5))
pairs(cbind(beta[2,], tt$logFC, tt_bc$logFC), 
      labels = c("beta","limma orig","limma + ruv"),
      col = cols, cex=.5,
      lower.panel = panel.cor)

fdr_alpha <- 0.1
table(orig = tt$adj.P.Val < fdr_alpha,
      with_ruv = tt_bc$adj.P.Val < fdr_alpha,
      actually_de = beta[2,] != 0)

# RUV versions
# ruv_logFC <- rfit$betahat[1,]
# ruv_padj <- p.adjust(variance_adjust(rfit)$p.ebayes[1,], method = "BH")

# pairs(cbind(beta[2,], tt$logFC, ruv_logFC),
#       labels = c("beta","limma orig","ruv2 (correct SE)"),
#       col = cols, cex=.5,
#       lower.panel = panel.cor)

# table(orig = tt$adj.P.Val < fdr_alpha,
#       with_ruv = ruv_padj < fdr_alpha,
#       actually_de = beta[2,] != 0)