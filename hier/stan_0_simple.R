# Stan example 0: a simple model with a fixed prior
#
#   theta_i ~ N(0, 1)            i = 1, ..., n
#   x_i | theta_i ~ N(theta_i, x_sd)
#
# The prior is fixed (mean 0, SD 1) and x_sd is known, so Stan only
# needs to estimate the theta_i. Nothing is learned about the
# distribution of the thetas from the data; see stan_1_hyper.R for that.

set.seed(1)
n <- 50
theta <- rnorm(n, mean=0, sd=1)
x_sd <- .75
x <- rnorm(n, theta, x_sd)
dat <- list(n=n, x_sd=x_sd, x=x)

library(rstan)
library(here)
fit <- stan(file = here("hier","stan_0_simple.stan"), 
            data = dat)
stan_plot(fit)

tab <- summary(fit)$summary

plot(theta, x)
points(theta, tab[1:n,1], col="red")

# data
sqrt(mean((x - theta)^2))
# prior mean
sqrt(mean(theta^2))
# posterior means
sqrt(mean((tab[1:n,1] - theta)^2))

# optimal shrinkage estimator from hier_simulation.R, where here each
# theta has a single observation (one replicate, unlike the n
# replicates there), mu_theta = 0, sigma2_theta = 1 and
# sigma2_x = x_sd^2:
#   theta_hat = (1 - B) * x + B * mu_theta
#   B = sigma2_x / (sigma2_x + sigma2_theta)
B <- x_sd^2 / (x_sd^2 + 1)
B
theta_B <- (1 - B) * x + B * 0
sqrt(mean((theta_B - theta)^2))

# for this model the posterior mean is exactly the optimal shrinkage
# estimator, so Stan should agree up to Monte Carlo error
plot(theta_B, tab[1:n,1], xlab="optimal B estimator",
     ylab="Stan posterior mean")
abline(0, 1, col="red")

