# simple hierarchical model
n <- 1000
theta <- rnorm(n)
x_sd <- 1
x <- rnorm(n, theta, x_sd)
var(x)

# add two points of complexity: 1) heavy tails 2) est of SD
set.seed(1)
n <- 5000
theta <- rt(n, df=3)
x_sd <- rgamma(n, 2, 1)
x <- rnorm(n, theta, x_sd)
df <- 10
x_sd_est <- x_sd/df * rchisq(n, df=df)

plot(x_sd, x_sd_est)

library(ggplot2)
suppressPackageStartupMessages(library(dplyr))
library(tibble)

dat <- tibble(theta, x, x_sd, x_sd_est)

# show how large ests can come from theta or SD
# black 'x' - the data
# blue dot - the true value
# blue line - the SD of the data
dat %>% 
  slice(1:50) %>%
  arrange(x_sd) %>%
  mutate(idx = seq_along(x)) %>%
  ggplot(aes(idx, x)) + 
  geom_point(shape=4,size=3,stroke=2) + 
  geom_pointrange(aes(idx, theta, 
                      ymin=theta-2*x_sd, 
                      ymax=theta+2*x_sd),
                  col="blue")

# add ranking by data and true value 'theta'
dat <- dat %>% 
  arrange(-abs(x)) %>%
  mutate(rank_x = seq_along(x),
         rank = rank(-abs(theta)))

# top ranked by data
# in blue is the theta
# -- a mix of high SD ones and true large 'x'
dat %>%
  slice(1:500) %>%
  ggplot(aes(rank_x, x)) + 
  geom_point(size=.5) + 
  geom_point(aes(rank_x, theta, color=log(x_sd)))

# fit a hierarchical model: adaptive shrinkage, "ash"
library(ashr)
fit <- with(
  dat, ash(x, x_sd_est, method="shrink"))

# add the posterior mean to the dataset
dat <- dat %>%
  mutate(posterior_mean_theta = fit$result$PosteriorMean,
         rank_ash = rank(-abs(posterior_mean_theta)))

# rank by ash posterior mean
# -- less blue points on the left side
dat %>%
  arrange(rank_ash) %>%
  slice(1:500) %>%
  ggplot(aes(rank_ash, posterior_mean_theta)) + 
  geom_point(size=.5) + 
  geom_point(aes(rank_ash, theta, col=log(x_sd)))

# make a CAT plot: concordance at the top
cuts <- seq(from=25,to=1000,by=25)
cat <- matrix(nrow=length(cuts), ncol=2)
var <- c("rank_x","rank_ash")
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
plot(cuts, cat[,"rank_x"], type="o", col=2, ylab="", ylim=c(0,.8))
points(cuts, cat[,"rank_ash"], type="o", col=4)
legend("bottomright", var, pch=1, lty=1, col=c(2,4), inset = .05)

