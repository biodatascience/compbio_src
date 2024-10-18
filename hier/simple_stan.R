# simple stan model
n <- 1000
theta <- rnorm(n)
x_sd <- 1
x <- rnorm(n, theta, x_sd)
dat <- list(n=n, x_sd=x_sd, x=x)

library(rstan)
library(here)
fit <- stan(file = here("hier","simple.stan"), 
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
