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
// Simple Bayesian Synthetic control model with multiple treated units and
// a Gaussian Process.
// Model follows this paper: ...
// TODO(jvives): Specify the paper.
// This implementation does not include staggered adoption or covariates.
data {
   int<lower=1> N;                    // number of pre-intervention periods
   int<lower=1> I;                    // number of treated units
   array[I] vector[N] y;                    // outcome
   int<lower=0> K;                    // number of donors
   matrix[N,K] X;                     // pre-intervention outcome for donors
   int<lower=1> N_pred;               // number of post-intervention periods
   matrix[N_pred,K] X_pred;           // post-intervention outcome for donors
}

transformed data{
   matrix[N, K] X_std;
   matrix[N_pred, K] X_pred_std;
   array[I] vector[N] y_std;
   array[I] real mean_y;
	array[I] real sd_y;
   array[N + N_pred] real time ;
   int<lower=1> sumN = N + N_pred;

   for(t in 1:sumN){
     time[t] = t;
   }

   {
     vector[K] mean_X;
  	  vector[K] sd_X;

  	 for (k in 1:K) {
      		mean_X[k] = mean(X[,k]);
      		sd_X[k] = sd(X[,k]);
      		X_std[,k] = (X[,k] - mean_X[k]) / sd_X[k];
      		X_pred_std[,k] = (X_pred[,k] - mean_X[k]) / sd_X[k];
      	}

     for (i in 1:I){
        mean_y[i] = mean(y[i]);
	      sd_y[i] = sd(y[i]);
        y_std[i] = (y[i] - mean_y[i]) / sd_y[i];
     }
   }
}

parameters {
   array[I] real<lower=0> sigma;
   array[I] simplex[K] beta;
   array[I] real<lower=0> rho;
   array[I] real<lower=0> alpha;
   array[I] vector[sumN] eta;
}

transformed parameters {
  array[I] vector[sumN] f;
  for(i in 1:I){
  matrix[sumN,sumN] K_matrix = gp_exp_quad_cov(time,alpha[i],rho[i]) +
                          diag_matrix(rep_vector(1e-9, sumN));
  matrix[sumN,sumN] L_K = cholesky_decompose(K_matrix);
  f[i] = L_K * eta[i];
  }
}

model {
   for(i in 1:I){
      sigma[i] ~ normal(0,1);
      rho[i] ~ normal(0,3);
      alpha[i] ~ normal(0,1);
      eta[i] ~ normal(0,1);
      target += normal_lpdf(y_std[i] | X_std*beta[i] +
                            f[i][1:N], sigma[i]);
   }
}

generated quantities {
   array[I] vector[N] y_sim;
   array[I] vector[N_pred] y_pred;
   for(i in 1:I){
      for (n in 1:N) {
      y_sim[i][n] = normal_rng(X_std[n,]*beta[i] +
                               f[i][n], sigma[i]) * sd_y[i] + mean_y[i]; //
     }
     for (j in 1:N_pred) {
      y_pred[i][j] = normal_rng(X_pred_std[j,]*beta[i] +
                                f[i][N + j], sigma[i]) * sd_y[i] +
                                mean_y[i]; //
     }
   }
}
