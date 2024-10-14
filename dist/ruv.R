library(limma)
library(ruv)

# sample size, number of genes, number of factors
n <- 20
m <- 400
k <- 5

W <- matrix(rnorm(n * k), ncol=k)

# spike and slab model: prob slab = non-zero effects
p_slab <- .5
alpha <- t(replicate(k, {
  ifelse(
    rbinom(m, size=1, prob=p_slab) == 1,
    rep(0,m),
    rnorm(m)
  )
}))

# make the alphas decreasing in effect
alpha <- alpha * (5:1 / 2.5)

p_de <- .25 # percent DE genes
X <- cbind(rep(1,n), W[,1] + rnorm(n)) # partial confound
#X <- cbind(rep(1,n), rep(0:1,each=n/2)) # orthogonal to W
cor(W[,1], X[,2])
beta <- rbind(rnorm(m), 
              c(rep(0, (1 - p_de) * m), 
                rnorm(p_de * m, 0, .5)))

# data model
Y <- X %*% beta + W %*% alpha + rnorm(n * m, 0, .5)

# standard limma fit
fit <- lmFit(t(Y), X)
efit <- eBayes(fit)
tt <- topTable(efit, coef = 2, number = m, sort.by = "none")

plot(beta[2,], tt$t)

# RUV fit using genes with high p-value, more factors than needed
rfit <- RUV2(Y, X[,2], ctl=which(tt$P.Value > .5), k=2*k) # from pval
#rfit <- RUV2(Y, X[,2], ctl=which(beta[2,] == 0), k=2*k) # oracle

image(cor(rfit$W, W), zlim=c(-1,1), axes = FALSE, 
      xlab="estimated W", ylab="true W",
      col=colorRampPalette(c("red","white","blue"))(99))
axis(1, at = seq(0, 1, length.out = 2*k), labels = 1:(2*k))
axis(2, at = seq(0, 1, length.out = k), labels = 1:k, las=1)

# batch corrected
fit_bc <- lmFit(t(Y), cbind(X, rfit$W[,1:k]))
efit_bc <- eBayes(fit_bc)
tt_bc <- topTable(efit_bc, coef = 2, number = m, sort.by = "none")

pairs(cbind(beta[2,], tt$logFC, tt_bc$logFC), 
      labels = c("beta","limma orig","limma + ruv"),
      lower.panel = panel.cor) # defined below

table(orig = tt$adj.P.Val < .1, 
      with_ruv = tt_bc$adj.P.Val < .1,
      actually_de = beta[2,] != 0)

######

panel.cor <- function(x, y, digits = 2, 
                      prefix = "", cex.cor, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- round(cor(x, y), digits) 
  txt <- paste0(prefix, r) 
  if (missing(cex.cor)) cex.cor <- 0.8 / strwidth(txt) 
  text(0.5, 0.5, txt, cex = cex.cor * abs(r)) 
}
