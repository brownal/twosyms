//
data {
  int P; // data points
  vector[P] y; // observed times to rescue
  vector[P] L; // light values for each data point
  vector[P] N; // nitrogen values for each data point
  vector[P] X;
  vector[P] kCO2;
  vector[P] kNPQ;
  vector[P] astar;
  vector[P] jSGm;
  vector[P] jCPm;
  vector[P] kROS;
  vector[P] b;
}
parameters {
  real alpha; // intercept (mean time to rescue)

  real beta_L;
  real beta_N;
  real beta_X;
  real beta_kCO2;
  real beta_kNPQ;
  real beta_astar;
  real beta_jSGm;
  real beta_jCPm;
  real beta_kROS;
  real beta_b;

  real gamma_LxL;
  real gamma_LxN;
  real gamma_LxX;
  real gamma_LxkCO2;
  real gamma_LxkNPQ;
  real gamma_Lxastar;
  real gamma_LxjSGm;
  real gamma_LxjCPm;
  real gamma_LxkROS;
  real gamma_Lxb;

  real gamma_NxN;
  real gamma_NxX;
  real gamma_NxkCO2;
  real gamma_NxkNPQ;
  real gamma_Nxastar;
  real gamma_NxjSGm;
  real gamma_NxjCPm;
  real gamma_NxkROS;
  real gamma_Nxb;

  real gamma_XxX;
  real gamma_XxkCO2;
  real gamma_XxkNPQ;
  real gamma_Xxastar;
  real gamma_XxjSGm;
  real gamma_XxjCPm;
  real gamma_XxkROS;
  real gamma_Xxb;

  real gamma_kCO2xkCO2;
  real gamma_kCO2xkNPQ;
  real gamma_kCO2xastar;
  real gamma_kCO2xjSGm;
  real gamma_kCO2xjCPm;
  real gamma_kCO2xkROS;
  real gamma_kCO2xb;

  real gamma_kNPQxkNPQ;
  real gamma_kNPQxastar;
  real gamma_kNPQxjSGm;
  real gamma_kNPQxjCPm;
  real gamma_kNPQxkROS;
  real gamma_kNPQxb;

  real gamma_astarxastar;
  real gamma_astarxjSGm;
  real gamma_astarxjCPm;
  real gamma_astarxkROS;
  real gamma_astarxb;

  real gamma_jSGmxjSGm;
  real gamma_jSGmxjCPm;
  real gamma_jSGmxkROS;
  real gamma_jSGmxb;

  real gamma_jCPmxjCPm;
  real gamma_jCPmxkROS;
  real gamma_jCPmxb;

  real gamma_kROSxkROS;
  real gamma_kROSxb;

  real gamma_bxb;

  real<lower=0> sigma; // standard deviation
}

model {
  alpha ~ normal(0, 5);

  beta_L ~ normal(0, 5);
  beta_N ~ normal(0, 5);
  beta_X ~ normal(0, 5);
  beta_kCO2 ~ normal(0, 5);
  beta_kNPQ ~ normal(0, 5);
  beta_astar ~ normal(0, 5);
  beta_jSGm ~ normal(0, 5);
  beta_jCPm ~ normal(0, 5);
  beta_kROS ~ normal(0, 5);
  beta_b ~ normal(0, 5);

  gamma_LxL ~ normal(0, 5);
  gamma_LxN ~ normal(0, 5);
  gamma_LxX ~ normal(0, 5);
  gamma_LxkCO2 ~ normal(0, 5);
  gamma_LxkNPQ ~ normal(0, 5);
  gamma_Lxastar ~ normal(0, 5);
  gamma_LxjSGm ~ normal(0, 5);
  gamma_LxjCPm ~ normal(0, 5);
  gamma_LxkROS ~ normal(0, 5);
  gamma_Lxb ~ normal(0, 5);

  gamma_NxN ~ normal(0, 5);
  gamma_NxX ~ normal(0, 5);
  gamma_NxkCO2 ~ normal(0, 5);
  gamma_NxkNPQ ~ normal(0, 5);
  gamma_Nxastar ~ normal(0, 5);
  gamma_NxjSGm ~ normal(0, 5);
  gamma_NxjCPm ~ normal(0, 5);
  gamma_NxkROS ~ normal(0, 5);
  gamma_Nxb ~ normal(0, 5);

  gamma_XxX ~ normal(0, 5);
  gamma_XxkCO2 ~ normal(0, 5);
  gamma_XxkNPQ ~ normal(0, 5);
  gamma_Xxastar ~ normal(0, 5);
  gamma_XxjSGm ~ normal(0, 5);
  gamma_XxjCPm ~ normal(0, 5);
  gamma_XxkROS ~ normal(0, 5);
  gamma_Xxb ~ normal(0, 5);

  gamma_kCO2xkCO2 ~ normal(0, 5);
  gamma_kCO2xkNPQ ~ normal(0, 5);
  gamma_kCO2xastar ~ normal(0, 5);
  gamma_kCO2xjSGm ~ normal(0, 5);
  gamma_kCO2xjCPm ~ normal(0, 5);
  gamma_kCO2xkROS ~ normal(0, 5);
  gamma_kCO2xb ~ normal(0, 5);

  gamma_kNPQxkNPQ ~ normal(0, 5);
  gamma_kNPQxastar ~ normal(0, 5);
  gamma_kNPQxjSGm ~ normal(0, 5);
  gamma_kNPQxjCPm ~ normal(0, 5);
  gamma_kNPQxkROS ~ normal(0, 5);
  gamma_kNPQxb ~ normal(0, 5);

  gamma_astarxastar ~ normal(0, 5);
  gamma_astarxjSGm ~ normal(0, 5);
  gamma_astarxjCPm ~ normal(0, 5);
  gamma_astarxkROS ~ normal(0, 5);
  gamma_astarxb ~ normal(0, 5);

  gamma_jSGmxjSGm ~ normal(0, 5);
  gamma_jSGmxjCPm ~ normal(0, 5);
  gamma_jSGmxkROS ~ normal(0, 5);
  gamma_jSGmxb ~ normal(0, 5);

  gamma_jCPmxjCPm ~ normal(0, 5);
  gamma_jCPmxkROS ~ normal(0, 5);
  gamma_jCPmxb ~ normal(0, 5);

  gamma_kROSxkROS ~ normal(0, 5);
  gamma_kROSxb ~ normal(0, 5);

  gamma_bxb ~ normal(0, 5);

  sigma ~ exponential(1);

  y ~ normal(alpha
   + beta_L * L
   + beta_N * N
   + beta_X * X
   + beta_kCO2 * kCO2
   + beta_kNPQ * kNPQ
   + beta_astar * astar
   + beta_jSGm * jSGm
   + beta_jCPm * jCPm
   + beta_kROS * kROS
   + beta_b * b
   + gamma_LxL * L .* L
   + gamma_LxN * L .* N
   + gamma_LxX * L .* X
   + gamma_LxkCO2 * L .* kCO2
   + gamma_LxkNPQ * L .* kNPQ
   + gamma_Lxastar * L .* astar
   + gamma_LxjSGm * L .* jSGm
   + gamma_LxjCPm * L .* jCPm
   + gamma_LxkROS * L .* kROS
   + gamma_Lxb * L .* b
   + gamma_NxN * N .* N
   + gamma_NxX * N .* X
   + gamma_NxkCO2 * N .* kCO2
   + gamma_NxkNPQ * N .* kNPQ
   + gamma_Nxastar * N .* astar
   + gamma_NxjSGm * N .* jSGm
   + gamma_NxjCPm * N .* jCPm
   + gamma_NxkROS * N .* kROS
   + gamma_Nxb * N .* b
   + gamma_XxX * X .* X
   + gamma_XxkCO2 * X .* kCO2
   + gamma_XxkNPQ * X .* kNPQ
   + gamma_Xxastar * X .* astar
   + gamma_XxjSGm * X .* jSGm
   + gamma_XxjCPm * X .* jCPm
   + gamma_XxkROS * X .* kROS
   + gamma_Xxb * X .* b
   + gamma_kCO2xkCO2 * kCO2 .* kCO2
   + gamma_kCO2xkNPQ * kCO2 .* kNPQ
   + gamma_kCO2xastar * kCO2 .* astar
   + gamma_kCO2xjSGm * kCO2 .* jSGm
   + gamma_kCO2xjCPm * kCO2 .* jCPm
   + gamma_kCO2xkROS * kCO2 .* kROS
   + gamma_kCO2xb * kCO2 .* b
   + gamma_kNPQxkNPQ * kNPQ .* kNPQ
   + gamma_kNPQxastar * kNPQ .* astar
   + gamma_kNPQxjSGm * kNPQ .* jSGm
   + gamma_kNPQxjCPm * kNPQ .* jCPm
   + gamma_kNPQxkROS * kNPQ .* kROS
   + gamma_kNPQxb * kNPQ .* b
   + gamma_astarxastar * astar .* astar
   + gamma_astarxjSGm * astar .* jSGm
   + gamma_astarxjCPm * astar .* jCPm
   + gamma_astarxkROS * astar .* kROS
   + gamma_astarxb * astar .* b
   + gamma_jSGmxjSGm * jSGm .* jSGm
   + gamma_jSGmxjCPm * jSGm .* jCPm
   + gamma_jSGmxkROS * jSGm .* kROS
   + gamma_jSGmxb * jSGm .* b
   + gamma_jCPmxjCPm * jCPm .* jCPm
   + gamma_jCPmxkROS * jCPm .* kROS
   + gamma_jCPmxb * jCPm .* b
   + gamma_kROSxkROS * kROS .* kROS
   + gamma_kROSxb * kROS .* b
   + gamma_bxb * b .* b, sigma);

}
