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
// Simple Bayesian Synthetic control model with multiple treated units.
// Model follows this paper: ...
// TODO(jvives): Specify the paper.
// This implementation does not include staggered adoption or covariates.
data {
   int<lower=1> N;                    // number of pre-intervention periods
   int<lower=1> I;                    // number of treated units
   vector[N] y[I];                    // outcome
   int<lower=0> K;                    // number of donors
   matrix[N,K] X;                     // pre-intervention outcome for donors
   int<lower=1> N_pred;               // number of post-intervention periods
   matrix[N_pred,K] X_pred;           // post-intervention outcome for donors

}

transformed data{
   matrix[N, K] X_std;
   matrix[N_pred, K] X_pred_std;
   vector[N] y_std[I];
   real mean_y[I];
	 real sd_y[I];
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
   real<lower=0> sigma[I];
   simplex[K] beta[I];
}

model {
   for(i in 1:I){
      sigma[i] ~ normal(0,1);
      target += normal_lpdf(y_std[i] | X_std*beta[i],
                            sigma[i]);
   }
}

generated quantities {
   vector[N] y_sim[I];
   vector[N_pred] y_pred[I];
   for(i in 1:I){
      for (n in 1:N) {
      y_sim[i][n] = normal_rng(X_std[n,]*beta[i], sigma[i]) * sd_y[i] +
                      mean_y[i]; //
     }
     for (j in 1:N_pred) {
      y_pred[i][j] = normal_rng(X_pred_std[j,]*beta[i], sigma[i]) * sd_y[i] +
                    mean_y[i]; //
     }
   }
}
