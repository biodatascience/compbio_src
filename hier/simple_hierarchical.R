# simple hierarchical model
n <- 1000
mu <- rnorm(n)
y_sd <- 1
y <- rnorm(n, mu, y_sd)
var(y)

# add two points of complexity: 1) heavy tails 2) est of SD
set.seed(1)
n <- 5000
mu <- rt(n, df=3)
y_sd <- rgamma(n, 2, 1)
y <- rnorm(n, mu, y_sd)
df <- 10
y_sd_est <- y_sd/df * rchisq(n, df=df)

plot(y_sd, y_sd_est)

library(ggplot2)
suppressPackageStartupMessages(library(dplyr))
library(tibble)

dat <- tibble(mu, y, y_sd, y_sd_est)

# show how large ests can come from mu or SD
dat %>% 
  slice(1:50) %>%
  arrange(y_sd) %>%
  mutate(idx = seq_along(y)) %>%
  ggplot(aes(idx, y)) + 
  geom_point(shape=4,size=3,stroke=2) + 
  geom_pointrange(aes(idx, mu, 
                      ymin=mu-2*y_sd, 
                      ymax=mu+2*y_sd),
                  col="blue")

# add ranking by data and true value 'mu'
dat <- dat %>% 
  arrange(-abs(y)) %>%
  mutate(rank_y = seq_along(y),
         rank = rank(-abs(mu)))

# top ranked by data
dat %>%
  slice(1:500) %>%
  ggplot(aes(rank_y, y)) + 
  geom_point(size=.5) + 
  geom_point(aes(rank_y, mu, color=log(y_sd)))

# add z-score, and show rank by z-score
dat %>%
  mutate(z = y/y_sd_est,
         rank_z = rank(-abs(z))) %>%
  arrange(rank_z) %>%
  slice(1:1000) %>%
  ggplot(aes(rank_z, z)) + 
  geom_point(size=.5) + 
  geom_point(aes(rank_z, mu, color=log(y_sd/y_sd_est))) +
  coord_cartesian(ylim=c(-15,15))

# fit a hierarchical model: adaptive shrinkage, "ash"
library(ashr)
fit <- with(
  dat, ash(y, y_sd_est, method="shrink"))

# add the posterior mean to the dataset
dat <- dat %>%
  mutate(posterior_mean_mu = fit$result$PosteriorMean,
         rank_ash = rank(-abs(posterior_mean_mu)))

# rank by ash posterior mean
dat %>%
  arrange(rank_ash) %>%
  slice(1:500) %>%
  ggplot(aes(rank_ash, posterior_mean_mu)) + 
  geom_point(size=.5) + 
  geom_point(aes(rank_ash, mu, col=log(y_sd)))

# make a CAT plot: concordance at the top
cuts <- seq(from=25,to=1000,by=25)
cat <- matrix(nrow=length(cuts), ncol=2)
var <- c("rank_y","rank_ash")
colnames(cat) <- var
for (i in seq_along(cuts)) {
  for (j in 1:2) {
    # here .data[[x]] is used to avoid copy-pasting code
    cat[i,j] <- dat %>% 
      filter(.data[[var[j]]] <= cuts[i]) %>%
      summarize(cat = mean(rank <= cuts[i])) %>%
      pull(cat)
  }
}

# CAT plot
plot(cuts, cat[,"rank_y"], type="o", col=2, ylab="", ylim=c(0,.8))
points(cuts, cat[,"rank_ash"], type="o", col=4)
legend("bottomright", var, pch=1, lty=1, col=c(2,4), inset = .05)

