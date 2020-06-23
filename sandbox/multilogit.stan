// code shameless taken from Stan documentation; comments added so those may be wrong
data {
  int K; // outcomes
  int N; // data points
  int D; // number of predictors
  int y[N]; // observed outcomes
  matrix[N, D] x; // predictor values for each data point
}
parameters {
  matrix[D, K] beta; // coefficients of each predictor for each outcome
}

model {
  matrix[N, K] x_beta = x * beta; // "probability" (softmax scale) of each outcome (columns) for each data point (rows)

  to_vector(beta) ~ normal(0, 2); // priors for the effect of each predictor

  for (n in 1:N)
    y[n] ~ categorical_logit(x_beta[n]');
}
