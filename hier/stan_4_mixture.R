# Stan example 4: mixture of two normals
#
#   z_i ~ Categorical(pi)              i = 1, ..., n
#   x_i | z_i = c ~ N(mu_c, sigma2_c)  c = 1, 2
#
# The same model we fit with EM in model/EM.qmd, now with priors
# and a posterior. The new piece is the discrete latent variable
# z_i: Stan's sampler (HMC) needs continuous parameters, so the Stan
# file sums z out of the likelihood.

set.seed(5)
n <- 400
pi <- c(0.7, 0.3)
mu <- c(0, 5)
sigma2 <- c(1, 4)
z <- sample(1:2, n, replace=TRUE, prob=pi)
x <- rnorm(n, mu[z], sqrt(sigma2[z]))
hist(x, breaks=40, col="grey")
dat <- list(n=n, x=x)

# some quantiles
quantile(x, c(.1, .25, .5, .75, .9))

# Like EM, the sampler can get stuck in a local mode, e.g. one wide
# component covering all the data and one narrow component in the middle.
# Start each chain from rough values based on the data.
init_fn <- function() {
  list(pi = c(.5, .5),
       mu = quantile(x, c(.25, .75)) + rnorm(2, 0, .25),
       sigma2 = c(2, 2))
}

library(rstan)
library(here)
fit <- stan(file = here("hier","stan_4_mixture.stan"),
            data = dat, init = init_fn)
check_hmc_diagnostics(fit)

print(fit)
rbind(pi, mu, sigma2)

post <- extract(fit)

# Gibbs sampler for the same model
#
# Gibbs sampling draws each unknown from its distribution
# conditional on all the others, cycling through:
#
#   z | pi, mu, sigma2, x
#   pi | z
#   mu | z, sigma2, x
#   sigma2 | z, mu, x
#
# Unlike HMC, Gibbs has no trouble with the discrete z: we just
# sample it. Each conditional has a closed form because the priors
# are conjugate:
#
#   pi ~ Dirichlet(1, 1)        (the flat simplex prior in Stan)
#   mu_c ~ N(0, 10^2)
#   sigma2_c ~ Inv-Gamma(1, 1)

# z_i given parameters: the E step from EM, but we draw z_i
# instead of keeping the probability
draw_z <- function(x, pi, mu, sigma2) {
  lik1 <- pi[1] * dnorm(x, mu[1], sqrt(sigma2[1]))
  lik2 <- pi[2] * dnorm(x, mu[2], sqrt(sigma2[2]))
  p_z2 <- lik2 / (lik1 + lik2)
  1 + rbinom(length(x), 1, p_z2)
}

# pi given z: Beta(1 + n_1, 1 + n_2), from counting the assignments
draw_pi <- function(z) {
  n1 <- sum(z == 1)
  n2 <- sum(z == 2)
  p1 <- rbeta(1, 1 + n1, 1 + n2)
  c(p1, 1 - p1)
}

# mu_c given z and sigma2_c: normal prior times normal likelihood.
# Precisions add, and the mean is a precision-weighted average of
# the prior mean (0) and the data in component c
draw_mu <- function(x, z, sigma2, prior_sd=10) {
  mu <- numeric(2)
  for (c in 1:2) {
    x_c <- x[z == c]
    prec <- 1 / prior_sd^2 + length(x_c) / sigma2[c]
    mean <- (sum(x_c) / sigma2[c]) / prec
    mu[c] <- rnorm(1, mean, sqrt(1 / prec))
  }
  mu
}

# sigma2_c given z and mu_c: Inv-Gamma(a + n_c/2, b + SS_c/2),
# where SS_c is the sum of squares around mu_c.
# If 1/s ~ Gamma, then s ~ Inv-Gamma
draw_sigma2 <- function(x, z, mu, a=1, b=1) {
  sigma2 <- numeric(2)
  for (c in 1:2) {
    x_c <- x[z == c]
    ss <- sum((x_c - mu[c])^2)
    sigma2[c] <- 1 / rgamma(1, a + length(x_c) / 2, b + ss / 2)
  }
  sigma2
}

# Cycle through the conditionals, storing the parameter draws.
# Like Stan, we start from the data-based values in init_fn. From
# there the labels stay put, because these two components are well
# separated
gibbs <- function(x, n_iter=2000) {
  init <- init_fn()
  pi <- init$pi
  mu <- init$mu
  sigma2 <- init$sigma2
  draws <- list(pi = matrix(NA, n_iter, 2),
                mu = matrix(NA, n_iter, 2),
                sigma2 = matrix(NA, n_iter, 2))
  for (t in 1:n_iter) {
    z <- draw_z(x, pi, mu, sigma2)
    pi <- draw_pi(z)
    mu <- draw_mu(x, z, sigma2)
    sigma2 <- draw_sigma2(x, z, mu)
    draws$pi[t,] <- pi
    draws$mu[t,] <- mu
    draws$sigma2[t,] <- sigma2
  }
  draws
}

g <- gibbs(x)

# trace plots: check that the chain has settled down
par(mfrow=c(3,1))
matplot(g$pi, type="l", lty=1, ylab="pi")
matplot(g$mu, type="l", lty=1, ylab="mu")
matplot(g$sigma2, type="l", lty=1, ylab="sigma2")
par(mfrow=c(1,1))

# drop the first half as warmup, and compare to Stan
keep <- 1001:2000
sapply(g, function(m) colMeans(m[keep,]))
sapply(post[c("pi","mu","sigma2")], colMeans)

# compare the posterior for mu_2
plot(density(post$mu[,2]), main="mu_2", lwd=2)
lines(density(g$mu[keep,2]), col="red", lwd=2)
legend("topright", c("Stan","Gibbs"), col=c("black","red"), lwd=2)

# compare the posterior for sigma2_2
plot(density(post$sigma2[,2]), main="sigma2_2", lwd=2)
lines(density(g$sigma2[keep,2]), col="red", lwd=2)
legend("topright", c("Stan","Gibbs"), col=c("black","red"), lwd=2)

# Optimization: find the posterior mode with optim()
#
# Instead of drawing from the posterior, we can find its peak, the
# maximum a posteriori (MAP) estimate. As in the Stan file, z is
# summed out of the likelihood. optim() works best without
# constraints, so it searches over transformed parameters:
#   theta = (logit pi_1, mu_1, mu_2, log sigma2_1, log sigma2_2)
neg_log_post <- function(theta, x) {
  pi1 <- plogis(theta[1])
  pi <- c(pi1, 1 - pi1)
  mu <- theta[2:3]
  sigma2 <- exp(theta[4:5])
  lik <- pi[1] * dnorm(x, mu[1], sqrt(sigma2[1])) +
         pi[2] * dnorm(x, mu[2], sqrt(sigma2[2]))
  log_prior <- sum(dnorm(mu, 0, 10, log=TRUE)) +
    sum(dgamma(1 / sigma2, 1, 1, log=TRUE) - 2 * log(sigma2))  # Inv-Gamma(1,1)
  # log posterior = log likelihood + log prior (up to a constant).
  # optim() minimizes by default, so return the negative
  -(sum(log(lik)) + log_prior)
}

init <- init_fn()
theta0 <- unname(c(qlogis(init$pi[1]), init$mu, log(init$sigma2)))
opt <- optim(theta0, neg_log_post, x=x, method="BFGS", hessian=TRUE)
opt$convergence  # 0 means it converged

# back-transform to the original scale. The mode is not the same as
# the posterior mean from Stan or Gibbs, because the posterior is skewed
c(pi1 = plogis(opt$par[1]),
  mu = opt$par[2:3],
  sigma2 = exp(opt$par[4:5]))

rbind(pi, mu, sigma2)

# The Hessian is the matrix of second derivatives of neg_log_post
# at the mode, i.e. the curvature of the peak. Its inverse
# approximates the posterior covariance of theta (a normal
# approximation to the posterior), so a sharper peak means smaller
# standard errors
opt$hessian
se <- sqrt(diag(solve(opt$hessian)))

# mu is not transformed, so its SEs compare directly to Stan
se[2:3]
apply(post$mu, 2, sd)

