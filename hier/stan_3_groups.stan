// two-level hierarchical model for per-gene slope estimates:
// genes within groups, groups within the population
//
// This is the "centered" form, which is easiest to read. For this
// model the "non-centered" form samples better: the SEs are large
// relative to sigma_gene, so b and sigma_gene are strongly dependent
// in the posterior (a "funnel"), and the centered form gives
// warnings (low E-BFMI, low ESS for sigma_gene, high Rhat). The
// non-centered form is the same model, written with helper variables
// w and z:
//   parameters: vector[K] w; vector[G] z;
//   transformed parameters:
//     vector[K] beta = mu_beta + sigma_group * w;
//     vector[G] b = beta[group] + sigma_gene * z;
//   model: w ~ normal(0, 1); z ~ normal(0, 1);
data {
  int G;                                 // number of genes
  int K;                                 // number of groups
  vector[G] x;                           // estimated slope of each gene
  vector<lower=0>[G] x_sd;               // standard error of each slope
  array[G] int<lower=1, upper=K> group;  // group of each gene
}
parameters {
  real mu_beta;
  real<lower=0> sigma_group;
  real<lower=0> sigma_gene;
  vector[K] beta;                        // group slopes
  vector[G] b;                           // gene slopes
}
model {
  mu_beta ~ normal(0, 5);
  sigma_group ~ normal(0, 2);
  sigma_gene ~ normal(0, 1);
  beta ~ normal(mu_beta, sigma_group);
  b ~ normal(beta[group], sigma_gene);
  x ~ normal(b, x_sd);
}
