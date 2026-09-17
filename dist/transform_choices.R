library(SummarizedExperiment)
library(matrixStats)
library(ggplot2)
library(patchwork)

# we will use this function later to assess trends 
# in the mean and variance across features
myMeanSdPlot <- function(x) {
  df <- data.frame(mean = rowMeans2(x), sd = rowSds(x))
  ggplot(df, aes(mean, sd)) +
    geom_hex() +
    geom_smooth()
}

# computes a PCA on the given assay, writing the loadings
# for PC1..PCk into rowData and the sample scores for PC1..PCk into colData
myPcaFunction <- function(se, assay_name = 1, num_pcs = 4) {
  pca <- prcomp(t(assay(se, assay_name)), center = TRUE, scale. = FALSE)
  pc_cols <- paste0("PC", seq_len(num_pcs))
  rowData(se)[, pc_cols] <- pca$rotation[, seq_len(num_pcs)]
  colData(se)[, pc_cols] <- pca$x[, seq_len(num_pcs)]
  se
}

# simulates data for 'se' on some original scale from one of a few
# distributions, spiking a continuous biological signal (colData(se)$biol_signal)
# into pct_spike% of rows (rowData(se)$spiked), each with an effect_size
# drawn from effect_scale * {-1,-.1,.1,1} (effect_scale scales how large
# the log2 fold change from biol_signal is). rowData(se)$feat_signal in
# [0,1] positions each row's baseline mean along mu_range (for binomial,
# mu_range is on the probability scale and gets multiplied by params$n).
# spiked rows are drawn only from the bottom low_var_frac fraction of
# features ranked by their (untransformed) variance h(mu_base), so the
# signal is hidden in low-count features (Poisson/NB) or
# extreme-probability features (binomial) -- exactly where an untransformed,
# Euclidean-distance PCA contributes least. A separate technical factor
# (colData(se)$tech_signal) is spiked into the top pct_tech_spike% of rows
# by h(mu_base) instead -- the highest-variance features, which dominate an
# untransformed PCA -- to mimic a technical artifact riding on the most
# highly expressed/highest-count features. assay 1 ("counts") holds the
# noisy draws, assay 2 ("truth") the noiseless mean mu_ij.
mySimulator <- function(se, dist = c("poisson", "binomial", "nb"),
                         params = list(), pct_spike = 20, mu_range,
                         effect_scale = 0.25, low_var_frac = 0.5,
                         pct_tech_spike = 2) {
  dist <- match.arg(dist)
  n_features <- nrow(se)
  n_samples <- ncol(se)
  feat_signal <- rowData(se)$feat_signal
  biol_signal <- colData(se)$biol_signal
  tech_signal <- colData(se)$tech_signal

  if (dist == "binomial") {
    p_base <- mu_range[1] + feat_signal * (mu_range[2] - mu_range[1])
    mu_base <- p_base * params$n
  } else {
    mu_base <- mu_range[1] + feat_signal * (mu_range[2] - mu_range[1])
  }

  h_of_mu <- switch(dist,
    poisson  = mu_base,
    binomial = mu_base * (1 - mu_base / params$n),
    nb       = mu_base + mu_base^2 / params$size
  )

  n_spike <- round(pct_spike / 100 * n_features)
  pool_size <- round(low_var_frac * n_features)
  if (n_spike > pool_size) {
    stop("pct_spike is too large for low_var_frac: need pct_spike <= low_var_frac * 100")
  }
  eligible <- order(h_of_mu)[seq_len(pool_size)]
  spike_idx <- sample(eligible, n_spike)
  spiked <- rep(FALSE, n_features)
  spiked[spike_idx] <- TRUE
  effect_size <- rep(0, n_features)
  effect_size[spike_idx] <- effect_scale * sample(c(-1, -0.1, 0.1, 1), n_spike, replace = TRUE)

  n_tech_spike <- round(pct_tech_spike / 100 * n_features)
  tech_idx <- order(h_of_mu, decreasing = TRUE)[seq_len(n_tech_spike)]
  tech_spiked <- rep(FALSE, n_features)
  tech_spiked[tech_idx] <- TRUE
  tech_effect_size <- rep(0, n_features)
  tech_effect_size[tech_idx] <- effect_scale * sample(c(-1, -0.1, 0.1, 1), n_tech_spike, replace = TRUE)

  # log2(mu_ij) = log2(mu_base_i) + effect_size_i * biol_signal_j + tech_effect_size_i * tech_signal_j
  log2_mu <- matrix(log2(mu_base), nrow = n_features, ncol = n_samples) +
    outer(effect_size, biol_signal) +
    outer(tech_effect_size, tech_signal)
  mu_true <- 2^log2_mu

  if (dist == "binomial") {
    prob_true <- pmin(pmax(mu_true / params$n, 1e-4), 1 - 1e-4)
    mu_true <- prob_true * params$n
  }

  sim <- switch(dist,
    poisson = matrix(rpois(length(mu_true), lambda = mu_true), nrow = n_features),
    binomial = matrix(rbinom(length(mu_true), size = params$n, prob = mu_true / params$n), nrow = n_features),
    nb = matrix(rnbinom(length(mu_true), size = params$size, mu = mu_true), nrow = n_features)
  )

  assays(se) <- list(counts = sim, truth = mu_true)
  rowData(se)$spiked <- spiked
  rowData(se)$effect_size <- effect_size
  rowData(se)$tech_spiked <- tech_spiked
  rowData(se)$tech_effect_size <- tech_effect_size
  se
}

# For a random variable X with E(X) = mu and Var(X) = h(mu), the delta
# method (see dist/vst_math.qmd) says a transform g with
# g(mu) = integral of C / sqrt(h(mu)) dmu gives Var(g(X)) approximately
# constant. Below we derive and check g for a few distributions.

n_features <- 2000
n_samples <- 20

## ---- Poisson: h(mu) = mu -> g(mu) = 2*sqrt(mu), i.e. sqrt transform ----
poisson_vst <- function(x) sqrt(x)

se_pois <- SummarizedExperiment(
  assays = list(counts = matrix(NA_real_, n_features, n_samples)),
  rowData = DataFrame(feat_signal = runif(n_features)),
  colData = DataFrame(biol_signal = runif(n_samples), tech_signal = runif(n_samples))
)
se_pois <- mySimulator(se_pois, dist = "poisson", mu_range = c(1, 500))

myMeanSdPlot(assay(se_pois, "counts"))                  # raw: variance grows with mean
myMeanSdPlot(poisson_vst(assay(se_pois, "counts")))      # after VST: flat
assay(se_pois, "vst") <- poisson_vst(assay(se_pois, "counts"))

## ---- Binomial(n=100): h(mu) = mu*(1-mu/n), a hump (not monotonic!) ----
## with p = mu/n: g(p) = 2*asin(sqrt(p)), the arcsine-sqrt transform
binomial_vst <- function(x, n) asin(sqrt(x / n))

n_binom <- 100
se_binom <- SummarizedExperiment(
  assays = list(counts = matrix(NA_real_, n_features, n_samples)),
  rowData = DataFrame(feat_signal = runif(n_features)),
  colData = DataFrame(biol_signal = runif(n_samples), tech_signal = runif(n_samples))
)
se_binom <- mySimulator(se_binom, dist = "binomial",
                         params = list(n = n_binom),
                         mu_range = c(0.01, 0.99))

myMeanSdPlot(assay(se_binom, "counts"))                       # raw: humped variance
myMeanSdPlot(binomial_vst(assay(se_binom, "counts"), n_binom)) # after VST: flat
assay(se_binom, "vst") <- binomial_vst(assay(se_binom, "counts"), n_binom)

# TODO: myFourPlots(se_binom) to check whether PCA recovers biol_signal vs
# tech_signal (see Poisson example above)

## ---- Negative Binomial(size=100): h(mu) = mu + mu^2/r ----
## g(mu) = 2*sqrt(r)*asinh(sqrt(mu/r)), which interpolates sqrt (mu << r,
## Poisson-like) and log (mu >> r) -- the same idea as DESeq2's vst() 
## but here with a fixed dispersion 1/r instead of one that
## depends on the mean
nb_vst <- function(x, r) 2 * sqrt(r) * asinh(sqrt(x / r))

r_nb <- 100
se_nb <- SummarizedExperiment(
  assays = list(counts = matrix(NA_real_, n_features, n_samples)),
  rowData = DataFrame(feat_signal = runif(n_features)),
  colData = DataFrame(biol_signal = runif(n_samples), tech_signal = runif(n_samples))
)
se_nb <- mySimulator(se_nb, dist = "nb", 
                     params = list(size = r_nb),
                     mu_range = c(1, 500))

myMeanSdPlot(assay(se_nb, "counts"))
myMeanSdPlot(nb_vst(assay(se_nb, "counts"), r_nb))
assay(se_nb, "vst") <- nb_vst(assay(se_nb, "counts"), r_nb)

# TODO: myFourPlots(se_nb) to check whether PCA recovers biol_signal vs
# tech_signal (see Poisson example above)
