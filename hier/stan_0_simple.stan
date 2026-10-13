data {
  int n;
  real x_sd;
  vector[n] x;
}
parameters {
  vector[n] theta;
}
model {
  x ~ normal(theta, x_sd);
  theta ~ normal(0, 1);
}
