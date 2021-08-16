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
// Simple Bayesian Synthetic control model.
// Model follows this paper: ...
// TODO(jvives): Document the paper.
// This implementation does not include multiple treated units or
// staggered adoption or additional covariates.
data {
   int<lower=1> N;                    // number of pre-intervention periods
   vector[N] y;                       // outcome
   int<lower=0> K;                    // number of donors
   matrix[N,K] X;                     // pre-intervention outcome for donors
   int<lower=1> N_pred;               // number of post-intervention periods
   matrix[N_pred,K] X_pred;           // post-intervention outcome for donors
}

transformed data{ // normalize using pre-treatment values
   matrix[N, K] X_std;
   matrix[N_pred, K] X_pred_std;
	 vector[K] mean_X;
	 vector[K] sd_X;
	 real mean_y = mean(y);
	 real sd_y = sd(y);
	 vector[N] y_std = (y - mean_y) / sd_y;

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
}

model {
   // Priors.
   sigma ~ normal(0,1);
     target += normal_lpdf(y_std | X_std*beta,
                           sigma);
}

generated quantities {
   vector[N] y_sim;
   vector[N_pred] y_pred;
   for (i in 1:N) {
      y_sim[i] = normal_rng(X_std[i,]*beta, sigma) * sd_y + mean_y; //
   }
   for (j in 1:N_pred) {
      y_pred[j] = normal_rng(X_pred_std[j,]*beta, sigma) * sd_y + mean_y; //
   }
}
