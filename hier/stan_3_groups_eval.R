# Evaluation for Stan example 3: run after stan_3_groups.R,
# which leaves the simulated data and the fit in the workspace.

b_hat <- colMeans(post$b)
beta_hat <- colMeans(post$beta)

# OLS vs posterior slopes, against the truth. The true group slopes
# are the horizontal dotted lines. The posterior pulls each gene
# toward its group's estimated slope
par(mfrow=c(1,2))
lims <- range(c(b, b_ols, b_hat))
plot(b, b_ols, col=cols[group], pch=20, xlim=lims, ylim=lims,
     xlab="true slope", ylab="estimate", main="OLS")
abline(0, 1, lty=2)
abline(h=beta, col=cols, lty=3)
plot(b, b_hat, col=cols[group], pch=20, xlim=lims, ylim=lims,
     xlab="true slope", ylab="estimate", main="posterior mean")
abline(0, 1, lty=2)
abline(h=beta, col=cols, lty=3)
par(mfrow=c(1,1))

rmse <- function(est) sqrt(mean((est - b)^2))
c(ols = rmse(b_ols),
  posterior = rmse(b_hat))

# group slopes: posterior intervals vs true values
ci <- apply(post$beta, 2, quantile, c(.025, .975))
plot(1:K, beta_hat, ylim=range(c(ci, beta)), pch=19, col=cols,
     xlab="group", ylab="group slope")
arrows(1:K, ci[1,], 1:K, ci[2,], angle=90, code=3, length=.05, col=cols)
points(1:K, beta, pch=4, cex=1.5)
