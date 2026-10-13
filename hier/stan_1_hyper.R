# Stan example 1: learn the prior from the data
#
#   theta_i ~ N(mu_theta, sigma_theta)  i = 1, ..., n
#   x_i | theta_i ~ N(theta_i, x_sd_i)
#
# Two changes from stan_0_simple.R:
# 1) mu_theta and sigma_theta are parameters, so the amount of
#    shrinkage is learned from the data rather than fixed in advance
# 2) each observation has its own known standard error x_sd_i
#    (e.g. genes with different numbers of reads)

set.seed(1)
n <- 500
mu_theta <- 1
sigma_theta <- 0.7
theta <- rnorm(n, mu_theta, sigma_theta)
x_sd <- rgamma(n, 2, 2)
x <- rnorm(n, theta, x_sd)
dat <- list(n=n, x=x, x_sd=x_sd)

library(rstan)
library(here)
fit <- stan(file = here("hier","stan_1_hyper.stan"),
            data = dat)
check_hmc_diagnostics(fit)

# the hyperparameters, compared to the true values
print(fit, pars=c("mu_theta","sigma_theta"))
c(mu_theta=mu_theta, sigma_theta=sigma_theta)

post <- extract(fit)
