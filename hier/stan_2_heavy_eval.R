# Evaluation for Stan example 2: run after stan_2_heavy.R,
# which leaves the simulated data and the fit in the workspace.

# the true nu is 3
hist(post_t$nu, breaks=40, col="grey", main="posterior of nu")
abline(v=3, col="red", lwd=2)

theta_norm <- colMeans(post_norm$theta)
theta_t <- colMeans(post_t$theta)

# compare the shrinkage: the normal prior is a straight line for a
# given x_sd, while the t prior shrinks small x and leaves large x
par(mfrow=c(1,2))
plot(x, theta_norm, col=rgb(0,0,0,.3), main="normal prior",
     xlab="x (data)", ylab="posterior mean")
abline(0, 1, col="red")
plot(x, theta_t, col=rgb(0,0,0,.3), main="t prior",
     xlab="x (data)", ylab="posterior mean")
abline(0, 1, col="red")
par(mfrow=c(1,1))

rmse <- function(est, idx=seq_len(n)) sqrt(mean((est[idx] - theta[idx])^2))

# overall
c(data = rmse(x),
  normal = rmse(theta_norm),
  t = rmse(theta_t))

# among the 50 largest true effects: any shrinkage pulls these toward
# 0, so the data x is hard to beat here, but the t prior does better
# than the normal prior
top <- order(-abs(theta))[1:50]
c(data = rmse(x, top),
  normal = rmse(theta_norm, top),
  t = rmse(theta_t, top))

# among the rest
c(data = rmse(x, -top),
  normal = rmse(theta_norm, -top),
  t = rmse(theta_t, -top))

# optional: compare to adaptive shrinkage (ash), which estimates a
# flexible unimodal prior instead of assuming a parametric family
if (requireNamespace("ashr", quietly=TRUE)) {
  fit_ash <- ashr::ash(x, x_sd, method="shrink")
  theta_ash <- fit_ash$result$PosteriorMean
  print(c(ash = rmse(theta_ash), ash_top = rmse(theta_ash, top)))
}
