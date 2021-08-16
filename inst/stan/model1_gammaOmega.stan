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
// Simple Bayesian Synthetic control model with predictor match.
// Model follows this paper: ...
// TODO(jvives): Specify the paper.
// This implementation does not include multiple treated units or
// staggered adoption, but allows for predictors/covariates.
data {
   // Data requirements
   int<lower=1> K;                     // number of predictors
   vector[K] X1;                       // treated unit outcome predictors
   int<lower=0> J;                     // number of donors
   matrix[K,J] X0;                     // donor units outcome predictors
   int<lower=1> T_post;                // number of post-intervention periods
   matrix[T_post,J] X0_pred;           // post-intervention outcome for donors
   vector<lower=0>[K] vs;              // predictor weights
}

transformed data{
   // Normalize design matrices using pre-treatment values.
   matrix[K, J] X0_std;
   matrix[T_post, J] X0_pred_std;
	 vector[J] mean_X0;
	 vector[J] sd_X0;
	 real mean_X1 = mean(X1);
	 real sd_X1 = sd(X1);
	 vector[K] X1_std = (X1 - mean_X1) / sd_X1;     // standarized X1
   vector[K] vs_std;                   // standarize the V weights
   for (k in 1:K) {
    vs_std[k] = pow(sd(X0[k,]), -1);
   }
   for (j in 1:J) {
		mean_X0[j] = mean(X0[,j]);
		sd_X0[j] = sd(X0[,j]);
		X0_std[,j] = (X0[,j] - mean_X0[j]) / sd_X0[j];           // standarized X0
		X0_pred_std[,j] = (X0_pred[,j] - mean_X0[j]) / sd_X0[j]; // standarized X0pred
	}
}

parameters {
   real<lower=0> sigma;              // std. dev predictors
   simplex[J] w;                     // synthetic control weights, Dir(1) prior
   simplex[K] gamma;                 // scaling vector for predictor weights
}

transformed parameters{
  // reparametrization
  // Transform variance parameters
  vector[K] Omega;
  vector[K] Gamma;
  for (k in 1:K) {
    Gamma[k] = pow(gamma[k], -1);
  }
  Omega = sigma*Gamma;
}

model {
   // Priors.
   sigma ~ normal(0,1);
   gamma ~ dirichlet(vs_std);
   target += normal_lpdf(X1_std | X0_std*w, Omega);
}

generated quantities {
   vector[K] X1_sim;                    // conditional posterior for X1
   vector[T_post] X1_pred;              // prediction for post-treatment period
   for (i in 1:K) {
      X1_sim[i] = normal_rng(X0_std[i,]*w, sigma) * sd_X1 + mean_X1; //
   }
   // TODO(jvives): Decide whether to use vs also in generated quantities.
   for (j in 1:T_post) {
      X1_pred[j] = normal_rng(X0_pred_std[j,]*w, sigma) * sd_X1 + mean_X1; //
   }
}
