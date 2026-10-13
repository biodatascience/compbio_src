# Stan example 2: a heavy-tailed prior
#
#   theta_i ~ t_nu(mu_theta, sigma_theta)  i = 1, ..., n
#   x_i | theta_i ~ N(theta_i, x_sd_i)
#
# The change from stan_1_hyper.R: the thetas follow a Student t
# distribution, and the degrees of freedom nu is also learned.
# This matches the simulation in simple_hierarchical.R, where most
# thetas are near 0 but a few are very large (e.g. most genes do not
# change, a few change a lot).
#
# A normal prior has to choose a single sigma_theta to describe both
# the many small effects and the few large ones. It ends up shrinking
# the large effects too much. The t prior can shrink small values
# strongly while leaving large values mostly alone.

set.seed(1)
n <- 1000
theta <- rt(n, df=3)
x_sd <- rgamma(n, 2, 2)
x <- rnorm(n, theta, x_sd)
dat <- list(n=n, x=x, x_sd=x_sd)

library(rstan)
library(here)

# fit the normal prior model from example 1, and the t prior model
fit_norm <- stan(file = here("hier","stan_1_hyper.stan"),
                 data = dat)
fit_t <- stan(file = here("hier","stan_2_heavy.stan"),
              data = dat)
check_hmc_diagnostics(fit_norm)
check_hmc_diagnostics(fit_t)

print(fit_norm, pars=c("mu_theta","sigma_theta"))
print(fit_t, pars=c("mu_theta","sigma_theta","nu"))

post_norm <- extract(fit_norm)
post_t <- extract(fit_t)
