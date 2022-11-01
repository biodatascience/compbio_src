# simple stan model
n <- 1000
theta <- rnorm(n)
x_sd <- 1
x <- rnorm(n, theta, x_sd)
dat <- list(n=n, x_sd=x_sd, x=x)

library(rstan)
fit <- stan(file = "simple.stan", data = dat)
stan_plot(fit)

tab <- summary(fit)$summary

plot(theta, x)
points(theta, tab[1:n,1], col="red")

# data
sqrt(mean((theta - x)^2))
# prior mean
sqrt(mean(theta^2))
# posterior means
sqrt(mean((theta - tab[1:n,1])^2))
