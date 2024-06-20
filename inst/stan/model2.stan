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
// Simple Bayesian Synthetic control model with a Gaussian Process.
// Model follows this paper: ...
// TODO(jvives): Specify the paper.
// staggered adoption, but allows for predictors/covariates.
data {
   int<lower=1> N;                    // number of pre-intervention periods
   vector[N] y;                       // outcome
   int<lower=0> K;                    // number of donors
   matrix[N,K] X;                     //  pre-intervention outcome for donors
   int<lower=1> N_pred;               // number of post-intervention periods
   matrix[N_pred,K] X_pred;           //  post-intervention outcome for donors
}

transformed data{
   matrix[N, K] X_std;
   matrix[N_pred, K] X_pred_std;
	 vector[K] mean_X;
	 vector[K] sd_X;
	 real mean_y = mean(y);
	 real sd_y = sd(y);
	 real time[N + N_pred] ;
	 vector[N] y_std = (y - mean_y) / sd_y;
	 int<lower=1> sumN = N + N_pred;

	 for(t in 1:sumN){
	    time[t] = t;
	 }

	 for (k in 1:K) {
		mean_X[k] = mean(X[,k]);
		sd_X[k] = sd(X[,k]);
		X_std[,k] = (X[,k] - mean_X[k]) / sd_X[k];
		X_pred_std[,k] = (X_pred[,k] - mean_X[k]) / sd_X[k];
	}
}

parameters {
   real<lower=0> sigma;
   simplex[K] beta;
   real<lower=0> rho;
   real<lower=0> alpha;
   vector[sumN] eta;
}

transformed parameters {
  vector[sumN] f;
  {
  matrix[sumN,sumN] K_matrix = cov_exp_quad(time,alpha,rho) +
                          diag_matrix(rep_vector(1e-9, sumN));
  matrix[sumN,sumN] L_K = cholesky_decompose(K_matrix);
  f = L_K * eta;
  }
}

model {
   // Priors.
   rho ~ normal(0,3);
   alpha ~ normal(0,1);
   sigma ~ normal(0,1);
   eta ~ normal(0,1);
     target += normal_lpdf(y_std | X_std*beta + f[1:N],
                           sigma);
}

generated quantities {
   vector[N] y_sim;
   vector[N_pred] y_pred;
   for (i in 1:N) {
      y_sim[i] = normal_rng(X_std[i,]*beta + f[i], sigma) * sd_y + mean_y; //
   }
   for (j in 1:N_pred) {
      y_pred[j] = normal_rng(X_pred_std[j,]*beta + f[N + j], sigma) * sd_y + mean_y; //
   }
}
