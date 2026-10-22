// mixture of two normals
//
// Stan cannot sample discrete parameters, so the component
// assignment z is summed out of the likelihood:
//   p(x_i) = pi_1 N(x_i | mu_1, sigma2_1) + pi_2 N(x_i | mu_2, sigma2_2)
// log_mix(p, lp1, lp2) computes log(p exp(lp1) + (1-p) exp(lp2)).
data {
  int n;
  vector[n] x;
}
parameters {
  simplex[2] pi;                  // mixing proportions, sum to 1
  vector[2] mu;
  vector<lower=0>[2] sigma2;      // variances
}
model {
  mu ~ normal(0, 10);
  sigma2 ~ inv_gamma(1, 1);       // conjugate, to match the Gibbs sampler
  for (i in 1:n) {
    target += log_mix(pi[1],
                      normal_lpdf(x[i] | mu[1], sqrt(sigma2[1])),
                      normal_lpdf(x[i] | mu[2], sqrt(sigma2[2])));
  }
}
