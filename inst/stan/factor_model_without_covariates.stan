// Copyright 2021 Google LLC
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     https://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
// Bayesian Factor Model without observed covariates
// Follows from https://arxiv.org/abs/2103.16244 (Pinkney 2021)
// Counterfactual Model: y = F'beta + X'gamma + Delta + kappa
//   - y Observed outcome
//   - F Latent factor components
//   - beta Latent factor unit traits
//   - Delta Time specific component
//   - kappa Unit specific component

// Load factor functions
#include factor_functions.stan

data{
  int<lower=2> L;                      // number of factors
  int<lower=1> N;                      // number of pre-intervention periods
  row_vector[N] y_treated_pre;         // pre intervention outcome for the
                                       //   treated unit
  int<lower=0> J;                      // number of donors
  row_vector[N] y_donors_pre[J];       // matrix of order [J,N] with
                                       //   pre-intervention outcome for donors
  int<lower=1> N_pred;                 // number of post-intervention periods
  row_vector[N_pred] y_donors_post[J]; // post-intervention outcome for donors
}

transformed data{
  int T = N + N_pred;
  int j_plus_1 = J+1;
  int<lower=1> M = L * (T - L) + L * (L - 1) / 2;
  row_vector[j_plus_1] j_ones = rep_row_vector(1, j_plus_1);
  vector[T] t_ones = rep_vector(1.0, T);

  row_vector[T] y_donors[J];

  row_vector[N] y_donors_pre_std[J];
  row_vector[N_pred] y_donors_post_std[J];
  vector[J] mean_y_donors_pre;
  vector[J] sd_y_donors_pre;
  real mean_y = mean(y_treated_pre);
  real sd_y = sd(y_treated_pre);
  row_vector[N] y_std = (y_treated_pre - mean_y) / sd_y;

  for (j in 1:J) {
    mean_y_donors_pre[j] = mean(y_donors_pre[j]);
    sd_y_donors_pre[j] = sd(y_donors_pre[j]);
    y_donors_pre_std[j] =
      (y_donors_pre[j] - mean_y_donors_pre[j]) / sd_y_donors_pre[j];
    y_donors_post_std[j] =
      (y_donors_post[j] - mean_y_donors_pre[j]) / sd_y_donors_pre[j];
    y_donors[j] = append_col(y_donors_pre_std[j], y_donors_post_std[j]);
  }
}

parameters{
  vector[T] raw_b;
  real<lower=0> sigma_b;
  row_vector [j_plus_1] raw_c;
  real<lower=0> sigma_c;

  matrix[j_plus_1, L] beta_off;
  vector<lower=0, upper=1>[L] lambda;
  real<lower=0, upper=1> eta;
  vector<lower=0, upper=1>[j_plus_1] tau;

  row_vector[N_pred] y_missing;

  real<lower=0> sigma;

  vector<lower=0>[L] F_diag;
  vector[M] F_lower;
}

transformed parameters{
  matrix[L, j_plus_1] beta = make_beta(j_plus_1,
                                       beta_off,
                                       lambda,
                                       eta,
                                       tau);
  vector[T] b = raw_b * sigma_b;             // time random effects
  row_vector [j_plus_1] c = raw_c * sigma_c; // unit random effects
}

model{
  to_vector(beta_off) ~ std_normal();
  F_diag ~ std_normal();
  F_lower ~ normal(0, 2);
  raw_b ~ std_normal();
  sigma_b ~ std_normal();
  raw_c ~ std_normal();
  sigma_c ~ std_normal();
  sigma ~ std_normal();
  {
    matrix[T, L] F = make_F(T, F_diag, F_lower);
    row_vector[T] Y_target[1];
    row_vector[T] Y_temp[j_plus_1];
    Y_target[1] = append_col(y_std,
                             y_missing);

    Y_temp = append_array(Y_target,
                          y_donors);

    for (j in 1:j_plus_1)
        Y_temp[j]' ~ normal_id_glm(F,
                                   b + c[j],
                                   beta[ , j],
                                   sigma);
  }
}

generated quantities{
  vector[T] synth_out;
  {
    matrix[T, L] F_ = make_F(T, F_diag, F_lower);
    matrix[T, j_plus_1] Synth_ = F_ * beta +
                                 b * j_ones +
                                 t_ones * c ;

    for (t in 1:T)
      synth_out[t] = normal_rng(Synth_[t,1], sigma) * sd_y + mean_y ;
  }
}
