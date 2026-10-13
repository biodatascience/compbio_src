# Stan example 3: genes nested in groups, slopes over time
#
#   beta_k ~ N(mu_beta, sigma_group)            k = 1, ..., K groups
#   b_i ~ N(beta_{group(i)}, sigma_gene)        i = 1, ..., G genes
#   a_i ~ N(mu_a, sigma_a)
#   y_ij ~ N(a_i + b_i * time_j, sigma_y)       j = 1, ..., T times
#
# The change from the previous examples: instead of one level of
# shrinkage, there are two. Each gene's slope is pulled toward its
# group's slope, and each group's slope is pulled toward the overall
# mean. Genes in the same group (e.g. a pathway, or a co-regulated
# module) have similar slopes, so they can borrow strength from each
# other.
#
# The data are given to Stan in "long" format, one row per
# observation, with integer indices saying which gene each row belongs
# to and which group each gene belongs to. In Stan, `b[gene]` expands
# the G gene slopes to a vector of length N, and `beta[group]` expands
# the K group slopes to a vector of length G.
#
# In the Stan file the slopes are written in "non-centered" form,
# e.g. b = beta[group] + sigma_gene * z with z ~ N(0,1), instead of
# b ~ normal(beta[group], sigma_gene). It is the same model, but easier
# for the sampler when the data per gene are weak and sigma_gene is
# small: then b and sigma_gene are strongly dependent in the posterior
# (a "funnel" shape), while z and sigma_gene are not. In examples 1 and
# 2 there are many observations with fairly small SEs, and the
# ordinary "centered" form samples better.

set.seed(1)
K <- 6       # groups
per <- 10    # genes per group
G <- K * per
n_time <- 8
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

dat <- list(N = nrow(dat_long), G = G, K = K,
            y = dat_long$y, time = dat_long$time,
            gene = dat_long$gene, group = group)

library(rstan)
library(here)
fit <- stan(file = here("hier","stan_3_groups.stan"),
            data = dat, iter = 4000,
            control = list(adapt_delta = 0.95))
check_hmc_diagnostics(fit)

print(fit, pars=c("mu_beta","sigma_group","sigma_gene","sigma_y"))
c(mu_beta=mu_beta, sigma_group=sigma_group,
  sigma_gene=sigma_gene, sigma_y=sigma_y)

post <- extract(fit)
