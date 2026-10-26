# EM for motif finding on a DNAStringSet, following the model in motif.qmd:
# every window of width w is either a motif (theta) or background (theta_bg).
# Differences from motif.qmd, so that it can be used on real peaks:
#  * takes a DNAStringSet, uses the central `trim` bp of each sequence
#  * searches both strands (reverse complements are added as extra sequences)
#  * a small pseudocount keeps entries of theta from becoming exactly 0
#  * optionally starts from a k-mer (`seed`) instead of random values
# Returns theta (4 x w, rows A,C,G,T), lambda, loglik, and the number of
# motif windows (lambda * number of windows).
library(Biostrings)

seed_theta <- function(seed, w, strength=.7) {
  k <- nchar(seed)
  stopifnot(k <= w)
  theta <- matrix(.25, nrow=4, ncol=w, dimnames=list(c("A","C","G","T"), NULL))
  off <- floor((w - k) / 2)
  for (j in seq_len(k)) {
    theta[, off + j] <- (1 - strength) / 3
    theta[substr(seed, j, j), off + j] <- strength
  }
  theta
}

em_motif <- function(seqs, w, n_starts=10, niter=30, trim=100, seed=NULL,
                     pseudo=.01) {
  mid <- width(seqs)[1] / 2
  seqs <- subseq(seqs, mid - trim/2 + 1, width=trim)
  seqs <- c(seqs, reverseComplement(seqs))
  m <- as.matrix(seqs)
  x <- matrix(match(m, c("A","C","G","T")), nrow=nrow(m))
  n_win <- ncol(x) - w + 1
  win <- matrix(nrow=nrow(x) * n_win, ncol=w)
  for (j in seq_len(w)) {
    win[,j] <- as.vector(t(x[, j - 1 + seq_len(n_win)]))
  }
  theta_bg <- tabulate(x, 4) / length(x)
  lik_bg_win <- apply(matrix(theta_bg[win], ncol=w), 1, prod)
  is_nt <- lapply(1:4, function(k) win == k)
  run_em <- function(theta, lambda) {
    for (t in seq_len(niter)) {
      lik <- rep(1, nrow(win))
      for (j in seq_len(w)) lik <- lik * theta[cbind(win[,j], j)]
      f_motif <- lambda * lik
      f_bg <- (1 - lambda) * lik_bg_win
      loglik <- sum(log(f_motif + f_bg))
      post_motif <- f_motif / (f_motif + f_bg)
      lambda <- mean(post_motif)
      counts <- t(sapply(1:4, function(k) colSums(post_motif * is_nt[[k]])))
      theta <- sweep(counts + pseudo, 2, colSums(counts + pseudo), "/")
    }
    list(theta=theta, lambda=lambda, loglik=loglik)
  }
  if (is.null(seed)) {
    fits <- lapply(seq_len(n_starts), function(s) {
      theta <- matrix(runif(4 * w), nrow=4)
      run_em(sweep(theta, 2, colSums(theta), "/"), lambda=1/n_win)
    })
    fit <- fits[[which.max(sapply(fits, `[[`, "loglik"))]]
  } else {
    fit <- run_em(seed_theta(seed, w), lambda=1/n_win)
  }
  rownames(fit$theta) <- c("A","C","G","T")
  fit$n_motif <- fit$lambda * nrow(win)
  fit
}
