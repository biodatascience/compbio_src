# Stan example 3: genes nested in groups, slopes over time
#
#   beta_k ~ N(mu_beta, sigma_group)      k = 1, ..., K groups
#   b_i ~ N(beta_{group(i)}, sigma_gene)  i = 1, ..., G genes
#   x_i | b_i ~ N(b_i, x_sd_i)
#
# where x_i and x_sd_i are the estimated slope and its standard error
# from a per-gene linear regression of expression on time.
#
# The change from the previous examples: instead of one level of
# shrinkage, there are two. Each gene's slope is pulled toward its
# group's slope, and each group's slope is pulled toward the overall
# mean. Genes in the same group (e.g. a pathway, or a co-regulated
# module) have similar slopes, so they can borrow strength from each
# other.
#
# This is a "two-stage" analysis: first summarize each gene by its
# slope estimate and SE, then fit the hierarchical model to those
# summaries, as in examples 1 and 2. We could instead give Stan every
# observation and fit the regressions inside the model; the results
# are very similar, because the slope estimate and its SE carry the
# information about the slope. One difference: each SE here is
# estimated from only 16 points per gene, so it is itself noisy, and
# the model treats it as known (compare limma, which shrinks these
# variance estimates).
#
# In the Stan file, `beta[group]` expands the K group slopes to a
# vector of length G, one per gene. The model is written in the
# simple "centered" form, b ~ normal(beta[group], sigma_gene). Stan
# will give some warnings for this model (low E-BFMI, low ESS for
# sigma_gene): the SEs are large relative to sigma_gene, so b and
# sigma_gene are strongly dependent in the posterior (a "funnel"
# shape). See the comment in the Stan file for the "non-centered"
# form, which is the same model but samples better here. In examples
# 1 and 2 many observations have fairly small SEs, and the centered
# form samples well.

set.seed(1)
K <- 6       # groups
per <- 10    # genes per group
G <- K * per
n_time <- 16
time <- seq(0, 1, length=n_time)

mu_beta <- 1
sigma_group <- 1.5
sigma_gene <- 0.3
sigma_y <- 1

group <- rep(1:K, each=per)
beta <- rnorm(K, mu_beta, sigma_group)
b <- rnorm(G, beta[group], sigma_gene)
a <- rnorm(G, 5, 1)

# one row per observation
dat_long <- data.frame(gene = rep(1:G, each=n_time),
                       time = rep(time, times=G))
dat_long$y <- rnorm(nrow(dat_long),
                    a[dat_long$gene] + b[dat_long$gene] * dat_long$time,
                    sigma_y)

# trajectories, one panel per group, one line per gene
cols <- palette.colors(K + 1)[-1]
y_mat <- matrix(dat_long$y, nrow=n_time)  # time x gene
par(mfrow=c(2,3), mar=c(4.5,4.5,2,1))
for (k in 1:K) {
  matplot(time, y_mat[,group == k], type="l", lty=1, col=cols[k],
          ylim=range(y_mat), xlab="time", ylab="expression",
          main=paste("group", k))
}
par(mfrow=c(1,1))

# per-gene least squares slope and its standard error
ols <- t(sapply(1:G, function(i) {
  fit_i <- lm(y ~ time, data=dat_long, subset = gene == i)
  summary(fit_i)$coefficients["time", c("Estimate","Std. Error")]
}))
b_ols <- ols[,1]
b_se <- ols[,2]

dat <- list(G = G, K = K, x = b_ols, x_sd = b_se, group = group)

library(rstan)
library(here)
fit <- stan(file = here("hier","stan_3_groups.stan"),
            data = dat)
check_hmc_diagnostics(fit)

print(fit, pars=c("mu_beta","sigma_group","sigma_gene"))
c(mu_beta=mu_beta, sigma_group=sigma_group, sigma_gene=sigma_gene)

post <- extract(fit)
