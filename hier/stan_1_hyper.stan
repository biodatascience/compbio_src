// hierarchical normal model with learned hyperparameters
// and a known standard error for each observation
data {
  int n;
  vector[n] x;
  vector<lower=0>[n] x_sd;
}
parameters {
  real mu_theta;
  real<lower=0> sigma_theta;
  vector[n] theta;
}
model {
  mu_theta ~ normal(0, 5);
  sigma_theta ~ normal(0, 2);
  theta ~ normal(mu_theta, sigma_theta);
  x ~ normal(theta, x_sd);
}
