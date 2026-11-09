# A simple functional interface for fitting Hidden Markov Models
# (HMMs) with the Baum-Welch (EM) algorithm, following the notation of
# the lecture notes in hmm.qmd: pi (starting distribution), A
# (transition matrix), B (emission parameters), and alpha, beta, gamma,
# xi for the forward/backward quantities.
#
# Main function:
#
#   hmm_fit(y, n_states, family, init_pi, init_A, init_B, size=NULL,
#           update_B=TRUE, emission_weight=1)
#
# where family is one of "normal", "poisson", "binomial",
# "multinomial", and init_pi, init_A, init_B are functions that
# produce starting values:
#
#   init_pi(n)    -> vector of length n summing to 1
#   init_A(n)     -> n x n matrix with rows summing to 1
#   init_B(y, n)  -> list of emission parameters, per family:
#     normal:      list(mu=, sigma=)  each length n
#     poisson:     list(lambda=)      length n
#     binomial:    list(p=)           length n (also give 'size' to hmm_fit)
#     multinomial: list(prob=)        n x K matrix, rows summing to 1,
#                                     with y taking values 1..K (or a factor)
#
# Zeros in the initial A stay zero during EM (e.g. forbidding the
# 1 -> 3 transition in the copy number example).
#
# With update_B=FALSE, the emission parameters are held fixed at
# init_B, and EM only updates pi and A.
#
# emission_weight (w) multiplies the log emission probabilities,
# log B_{i,y_t}, before running forward-backward and Viterbi, i.e. B is
# raised to the power w. With w=1 (default) this is the usual HMM.
# Why would we want w < 1? The HMM assumes that the Y_t are independent
# given the hidden states. Sometimes we know they are not: e.g. if Y_t
# is the 3-mer of bases t-1, t, t+1, then Y_t and Y_{t+1} share two
# bases, and each base appears in three of the Y's. The HMM then counts
# the evidence from each base about three times over, and so the
# emissions are overconfident relative to the transitions A. The path
# will switch states too often, whenever a few Y's in a row happen to
# look a bit more like another state. Down-weighting the log emissions
# (e.g. w = 1/2 or 1/3) puts the emissions and transitions back on a
# more even footing. (In statistics this is related to "composite" or
# "tempered" likelihoods.) Note that with w != 1, the log likelihood
# reported is no longer the log likelihood of the data under the model.

# Emission families --------------------------------------------------
#
# Each family has two functions:
#   logdens(y, B, size) returns an n x T matrix of log P(Y_t=y_t | X_t=i)
#   mstep(y, gamma, size, B) returns updated B, using gamma (n x T) as
#     weights, i.e. a weighted maximum likelihood estimate per state
#   rand(x, B, size) draws y given a vector of hidden states x

hmm_families <- list(
  normal = list(
    logdens = function(y, B, size) {
      n <- length(B$mu)
      matrix(dnorm(rep(y, each=n), B$mu, B$sigma, log=TRUE), nrow=n)
    },
    mstep = function(y, gamma, size, B) {
      w <- rowSums(gamma)
      mu <- as.vector(gamma %*% y) / w
      sigma <- sqrt(rowSums(gamma * outer(mu, y, "-")^2) / w)
      list(mu=mu, sigma=sigma)
    },
    rand = function(x, B, size) {
      rnorm(length(x), B$mu[x], B$sigma[x])
    }
  ),
  poisson = list(
    logdens = function(y, B, size) {
      n <- length(B$lambda)
      matrix(dpois(rep(y, each=n), B$lambda, log=TRUE), nrow=n)
    },
    mstep = function(y, gamma, size, B) {
      list(lambda = as.vector(gamma %*% y) / rowSums(gamma))
    },
    rand = function(x, B, size) {
      rpois(length(x), B$lambda[x])
    }
  ),
  binomial = list(
    logdens = function(y, B, size) {
      n <- length(B$p)
      matrix(dbinom(rep(y, each=n), rep(size, each=n), B$p, log=TRUE),
             nrow=n)
    },
    mstep = function(y, gamma, size, B) {
      list(p = as.vector(gamma %*% y) / as.vector(gamma %*% size))
    },
    rand = function(x, B, size) {
      rbinom(length(x), size, B$p[x])
    }
  ),
  multinomial = list(
    logdens = function(y, B, size) {
      log(B$prob[, y, drop=FALSE])
    },
    mstep = function(y, gamma, size, B) {
      K <- ncol(B$prob)
      # n x K weighted counts: sum of gamma over the t where y_t = k
      counts <- vapply(seq_len(K), function(k) rowSums(gamma[, y == k, drop=FALSE]),
                       numeric(nrow(gamma)))
      list(prob = counts / rowSums(counts))
    },
    rand = function(x, B, size) {
      # inverse CDF: count how many cumulative probabilities u exceeds
      cum_prob <- t(apply(B$prob, 1, cumsum))
      u <- runif(length(x))
      1L + as.integer(rowSums(u > cum_prob[x, , drop=FALSE]))
    }
  )
)

# Forward-backward ---------------------------------------------------
#
# Uses the scaling of Rabiner and Juang: each column of alpha is
# scaled to sum to 1, with scaling constants c_t, and beta is scaled
# by the same constants. Then gamma = alpha * beta exactly, and the
# log-likelihood is sum(log(c_t)).
#
# The emission matrix b (n x T) is also scaled per column (subtracting
# the max log density at each t) to avoid underflow; this does not
# change gamma or xi, and we add the offset back to the log-likelihood.

hmm_forward_backward <- function(logb, pi, A) {
  n <- nrow(logb)
  T <- ncol(logb)
  m <- apply(logb, 2, max)
  b <- exp(logb - rep(m, each=n))

  alpha <- matrix(0, nrow=n, ncol=T)
  c <- numeric(T)
  alpha[,1] <- pi * b[,1]
  c[1] <- sum(alpha[,1])
  alpha[,1] <- alpha[,1] / c[1]
  for (t in 2:T) {
    alpha[,t] <- as.vector(alpha[,t-1] %*% A) * b[,t]
    c[t] <- sum(alpha[,t])
    alpha[,t] <- alpha[,t] / c[t]
  }

  beta <- matrix(0, nrow=n, ncol=T)
  beta[,T] <- 1
  for (t in (T-1):1) {
    beta[,t] <- as.vector(A %*% (b[,t+1] * beta[,t+1])) / c[t+1]
  }

  gamma <- alpha * beta

  # xi[i,j,t] = alpha[i,t] A[i,j] b[j,t+1] beta[j,t+1] / c[t+1]
  # we only need the sum over t for the M-step
  xi_sum <- A * (alpha[,-T,drop=FALSE] %*%
                 t(b[,-1,drop=FALSE] * beta[,-1,drop=FALSE] /
                   rep(c[-1], each=n)))

  list(gamma=gamma, xi_sum=xi_sum,
       loglik=sum(log(c)) + sum(m))
}

# Viterbi ------------------------------------------------------------
#
# Computed in log space: V[i,t] is the log probability of the most
# probable path ending in state i at time t, and Ptr points back along
# that path.

hmm_viterbi <- function(logb, pi, A) {
  n <- nrow(logb)
  T <- ncol(logb)
  logA <- log(A)
  V <- matrix(0, nrow=n, ncol=T)
  Ptr <- matrix(0L, nrow=n, ncol=T)
  V[,1] <- log(pi) + logb[,1]
  for (t in 2:T) {
    scores <- V[,t-1] + logA # scores[i,j]: from state i to state j
    Ptr[,t] <- apply(scores, 2, which.max)
    V[,t] <- apply(scores, 2, max) + logb[,t]
  }
  x_hat <- integer(T)
  x_hat[T] <- which.max(V[,T])
  for (t in (T-1):1) {
    x_hat[t] <- Ptr[x_hat[t+1], t+1]
  }
  x_hat
}

# Simulation ---------------------------------------------------------
#
# Draw a hidden chain x of length T from pi and A, then draw the
# observations y given x using the emission family and parameters B
# (same format as returned by init_B). Returns list(x=, y=).

hmm_simulate <- function(T, pi, A, B,
                         family=c("normal","poisson","binomial","multinomial"),
                         size=NULL) {
  family <- match.arg(family)
  n <- length(pi)
  if (family == "binomial") {
    if (is.null(size)) stop("binomial family requires 'size'")
    size <- rep_len(size, T)
  }
  x <- integer(T)
  x[1] <- sample(n, 1, prob=pi)
  for (t in seq_len(T)[-1]) {
    x[t] <- sample(n, 1, prob=A[x[t-1],])
  }
  y <- hmm_families[[family]]$rand(x, B, size)
  list(x=x, y=y)
}

# Main fitting function ----------------------------------------------

hmm_fit <- function(y, n_states,
                    family=c("normal","poisson","binomial","multinomial"),
                    init_pi, init_A, init_B,
                    size=NULL, update_B=TRUE, emission_weight=1,
                    max_iter=200, tol=1e-6, verbose=FALSE) {
  family <- match.arg(family)
  fam <- hmm_families[[family]]
  n <- n_states
  T <- length(y)
  stopifnot(T >= 2)

  if (family == "binomial") {
    if (is.null(size)) stop("binomial family requires 'size'")
    size <- rep_len(size, T)
  }
  if (family == "multinomial") {
    y <- as.integer(y) # factors become their level index 1..K
  }

  pi <- init_pi(n)
  A <- init_A(n)
  B <- init_B(y, n)
  stopifnot(length(pi) == n, abs(sum(pi) - 1) < 1e-8)
  stopifnot(nrow(A) == n, ncol(A) == n, all(abs(rowSums(A) - 1) < 1e-8))

  loglik <- numeric(0)
  converged <- FALSE
  for (iter in seq_len(max_iter)) {
    # E-step
    logb <- emission_weight * fam$logdens(y, B, size)
    fb <- hmm_forward_backward(logb, pi, A)
    loglik[iter] <- fb$loglik
    if (verbose) message("iter ", iter, ": loglik = ", round(fb$loglik, 4))
    if (iter > 1 && abs(loglik[iter] - loglik[iter-1]) < tol) {
      converged <- TRUE
      break
    }
    # M-step
    pi <- fb$gamma[,1]
    A <- fb$xi_sum / rowSums(fb$xi_sum)
    if (update_B) B <- fam$mstep(y, fb$gamma, size, B)
  }

  if (!converged) {
    # make gamma and Viterbi consistent with the final parameters
    logb <- emission_weight * fam$logdens(y, B, size)
    fb <- hmm_forward_backward(logb, pi, A)
    warning("EM did not converge in ", max_iter, " iterations")
  }

  list(pi=pi, A=A, B=B, family=family,
       loglik=loglik, gamma=fb$gamma,
       viterbi=hmm_viterbi(logb, pi, A),
       iterations=length(loglik), converged=converged)
}

# Demo: two state HMM with binomial outcome --------------------------

if (FALSE) {

  # simulate the hidden chain and observed successes out of 'size'
  # trials, where p depends on the state
  A_true <- matrix(c(.99, .01,
                     .01, .99), nrow=2, byrow=TRUE)
  p_true <- c(.1, .7)
  size <- 10
  set.seed(123)
  sim <- hmm_simulate(T=1000, pi=c(1, 0), A=A_true, B=list(p=p_true),
                      family="binomial", size=size)
  x <- sim$x
  y <- sim$y
  plot(y, col=x, ylim=c(-1.5, size), yaxt="n")
  axis(2, at=seq(0, size, by=2))
  # blue triangles at the bottom mark the time points where the state changes
  switch_t <- which(diff(x) != 0) + 1
  points(switch_t, rep(-1, length(switch_t)), pch=17, col="blue")

  # initialization functions
  init_pi <- function(n) rep(1/n, n)
  init_A <- function(n) {
    A <- matrix(.1/(n-1), nrow=n, ncol=n)
    diag(A) <- .9
    A
  }
  init_B <- function(y, n) list(p=seq(.3, .5, length.out=n))

  fit <- hmm_fit(y, n_states=2, family="binomial",
                 init_pi=init_pi, init_A=init_A, init_B=init_B,
                 size=size)

  round(fit$B$p, 4)
  round(fit$A, 4)
  fit$iterations
  plot(fit$loglik, type="b", xlab="iteration", ylab="log likelihood")

  # Viterbi path vs truth, and posterior probability of state 2
  table(truth=x, viterbi=fit$viterbi)
  old_par <- par(mfrow=c(2,1), mar=c(3,4,2,1))
  plot(y, col=fit$viterbi, main="data colored by Viterbi state")
  # black solid line: posterior probability of state 2 at each t, gamma[2,t]
  plot(fit$gamma[2,], type="l", main="P(X_t = 2 | Y, theta)",
       xlab="", ylab=expression(gamma[list(2,t)]))
  # red dashed line: the true hidden state, shifted to 0 (state 1) or 1 (state 2)
  lines(x - 1, col="red", lty=2)
  par(old_par) # restore the previous mfrow and mar
}
