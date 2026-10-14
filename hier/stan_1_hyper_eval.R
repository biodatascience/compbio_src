# Evaluation for Stan example 1: run after stan_1_hyper.R,
# which leaves the simulated data and the fit in the workspace.

hist(post$sigma_theta, breaks=40, col="grey",
     main="posterior of sigma_theta")
abline(v=sigma_theta, col="red", lwd=2)

theta_hat <- colMeans(post$theta)

# noisier observations (light colors) are pulled further toward
# mu_theta
cols <- hcl.colors(10, "viridis")[cut(log(x_sd), 10)]
plot(x, theta_hat, col=cols, pch=20,
     xlab="x (data)", ylab="posterior mean of theta")
abline(0, 1, lty=2)
abline(h=mu_theta, lty=3)

# the shrinkage factor B_i for each observation, from the posterior
# mean: theta_hat = (1 - B) * x + B * mu_theta. Compare to the
# formula from hier_shrink_derivation,
#   B = x_sd^2 / (x_sd^2 + sigma_theta^2)
# (only for x not too close to mu_theta, where this ratio is
# unstable)
mu_theta_hat <- mean(post$mu_theta)
ok <- abs(x - mu_theta_hat) > 0.5
B_hat <- (x - theta_hat) / (x - mu_theta_hat)
plot(x_sd[ok], B_hat[ok], ylim=c(0,1), xlab="x_sd", ylab="shrinkage factor B")
curve(x^2 / (x^2 + sigma_theta^2), add=TRUE, col="red", lwd=2)

# "oracle" posterior mean using the true mu_theta and sigma_theta
B <- x_sd^2 / (x_sd^2 + sigma_theta^2)
theta_oracle <- (1 - B) * x + B * mu_theta

rmse <- function(est) sqrt(mean((est - theta)^2))
c(data = rmse(x),
  posterior = rmse(theta_hat),
  oracle = rmse(theta_oracle))
