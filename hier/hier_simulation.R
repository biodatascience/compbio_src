# In-class exercise: shrinkage estimation in a simple hierarchical model
#
# The model:
#
#   theta_i ~ N(mu_theta, sigma2_theta)    i = 1, ..., m   (e.g. genes)
#   x_ij | theta_i ~ N(theta_i, sigma2_x)  j = 1, ..., n   (replicates)
#
# Our goal is to estimate each theta_i. The obvious estimator is the
# per-gene average, xbar_i. A shrinkage estimator pulls xbar_i toward
# the center of the distribution of the thetas:
#
#   theta_hat_i = (1 - B) * xbar_i + B * mu_theta
#
# where B in [0,1] is the shrinkage factor. B = 0 is no shrinkage
# (just use xbar_i), B = 1 is complete shrinkage (every gene gets
# the same estimate, mu_theta).

library(matrixStats)

set.seed(5)
m <- 1000
n <- 5
mu_theta <- 2
sigma2_theta <- 1.5
sigma2_x <- 4

theta <- rnorm(m, mu_theta, sqrt(sigma2_theta))
# row i of x has the n replicates for theta_i
x <- matrix(rnorm(m * n, theta, sqrt(sigma2_x)), nrow=m, ncol=n)
xbar <- rowMeans(x)

plot(theta, xbar, col=rgb(0,0,0,.3))
abline(0, 1, col="red")
abline(h=mu_theta, v=mu_theta, lty=2)

# Is the spread of xbar larger or smaller than the spread of
# theta? Why? What should var(xbar) be, in terms of the parameters?
var(theta)
var(xbar)

# The spread of xbar is larger. Each xbar_i is theta_i plus
# independent noise with variance sigma2_x / n, so the two sources of
# variance add:
#
#   var(xbar_i) = sigma2_theta + sigma2_x / n = 1.5 + 4/5 = 2.3
#
# This is why the most extreme xbar_i tend to overshoot their theta_i,
# and why pulling them back toward the middle can help.
sigma2_theta + sigma2_x / n

shrink <- function(xbar, B, mu) (1 - B) * xbar + B * mu

## Part 1: optimal shrinkage factor

# Here we pretend mu_theta, sigma2_theta and sigma2_x are known.
# Try a grid of shrinkage factors and compute the mean squared error
# of the estimates, averaging over all m genes.

B_grid <- seq(from=0, to=1, by=0.01)
mse <- sapply(B_grid, function(B) {
  mean((shrink(xbar, B, mu_theta) - theta)^2)
})
plot(B_grid, mse, type="l", xlab="B", ylab="MSE")
B_grid[which.min(mse)]

# Work out the optimal B analytically.
# See hier_shrink_derivation.typ
# (render to PDF with `typst compile`)

B_formula <- (sigma2_x / n) / (sigma2_x / n + sigma2_theta)
abline(v=B_formula, col="red")
B_formula
# theoretical MSE(B) = (1 - B)^2 * sigma2_x / n + B^2 * sigma2_theta
# (here curve() uses x for B), compare to the simulated MSE
curve((1 - x)^2 * sigma2_x / n + x^2 * sigma2_theta,
      add=TRUE, col="red", lty=2)

# How does the optimal B change when you increase n? When you
# increase sigma2_x? When you increase sigma2_theta? Make a guess,
# then change the parameters at the top and re-run to check.

## Part 2: bias and variance of the shrinkage estimator

# Hold the thetas fixed and redraw the data `x` many times.
# This lets us compute, for each gene, the bias and the variance of
# theta_hat_i over repeated experiments

nrep <- 500
xbar_rep <- replicate(nrep, {
  x_new <- matrix(rnorm(m * n, theta, sqrt(sigma2_x)), nrow=m, ncol=n)
  rowMeans(x_new)
})
dim(xbar_rep) # m genes x nrep experiments

bias_and_var <- function(B) {
  est <- shrink(xbar_rep, B, mu_theta)
  data.frame(theta = theta,
             bias = rowMeans(est) - theta,
             var = rowVars(est))
}

# Look at a few fixed shrinkage factors
B_values <- c(0, 0.25, 0.5, 0.75)
par(mfrow=c(2,4), mar=c(4.5,4.5,2,1))
for (B in B_values) {
  bv <- bias_and_var(B)
  plot(bv$theta, bv$bias, col=rgb(0,0,0,.3), ylim=c(-3.5,3.5),
       xlab="theta", ylab="bias", main=paste("B =", B))
  abline(h=0, v=mu_theta, lty=2)
}
for (B in B_values) {
  bv <- bias_and_var(B)
  plot(bv$theta, bv$var, col=rgb(0,0,0,.3), ylim=c(0,1),
       xlab="theta", ylab="variance", main=paste("B =", B))
}
par(mfrow=c(1,1))

# Average squared bias and average variance across genes, over the
# grid of B values
tradeoff <- t(sapply(B_grid, function(B) {
  bv <- bias_and_var(B)
  c(bias2 = mean(bv$bias^2), var = mean(bv$var))
}))
tradeoff <- data.frame(B = B_grid, tradeoff)
tradeoff$mse <- tradeoff$bias2 + tradeoff$var

plot(tradeoff$B, tradeoff$mse, type="l", lwd=2,
     ylim=c(0, max(tradeoff$mse)), xlab="B", ylab="")
lines(tradeoff$B, tradeoff$bias2, col="red", lwd=2)
lines(tradeoff$B, tradeoff$var, col="blue", lwd=2)
legend("top", c("MSE","bias^2","variance"),
       col=c("black","red","blue"), lwd=2, inset=.05)
abline(v=B_formula, lty=2)

# theoretical curves (dashed): bias^2 = B^2 * sigma2_theta (here using
# the observed spread of the thetas), variance = (1 - B)^2 * sigma2_x / n,
# and their sum, the MSE(B) from Part 1
curve(x^2 * mean((theta - mu_theta)^2), add=TRUE, col="red", lty=2)
curve((1 - x)^2 * sigma2_x / n, add=TRUE, col="blue", lty=2)
curve(x^2 * sigma2_theta + (1 - x)^2 * sigma2_x / n, add=TRUE, lty=2)
