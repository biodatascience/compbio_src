// hierarchical model with a heavy-tailed (Student t) prior
// and a known standard error for each observation
data {
  int n;
  vector[n] x;
  vector<lower=0>[n] x_sd;
}
parameters {
  real mu_theta;
  real<lower=0> sigma_theta;
  real<lower=1> nu;
  vector[n] theta;
}
model {
  mu_theta ~ normal(0, 5);
  sigma_theta ~ normal(0, 2);
  nu ~ gamma(2, 0.1);
  theta ~ student_t(nu, mu_theta, sigma_theta);
  x ~ normal(theta, x_sd);
}
