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
// Auxiliary functions for stan factor model.
// Functions:
//  make_F Creates the factor matrix F of size T x L.
//  make_beta Creates the unit loadings of size L x 1 for each J+1 unit.
functions{
  // @description Given the factor loadings creates the T x L factor matrix.
  // @param T Number of time periods.
  // @param diagonal_loadings Loadings for each factor.
  // @param lower_tri_loadings Cross-factor loadings.
  matrix make_F(int T,
                vector diagonal_loadings,
                vector lower_tri_loadings) {
  int L = num_elements(diagonal_loadings);
  int M = num_elements(lower_tri_loadings);
  matrix[T, L] F;

  int idx = 0; // Index for the lower diagonal

  for (j in 1:L) {
    F[j, j] = diagonal_loadings[j];
    for (i in (j + 1):T) {
      idx += 1;
      F[i, j] = lower_tri_loadings[idx];
    }
  }
  for (j in 1:(L - 1)) {
    for (i in (j + 1):L) F[j, i] = 0;
  }

  return F;
  }

  // @description Given prior parameters generates the J+1xL matrix of unit
  // weights for each factor from the prior model.
  // @param J Number of donor units.
  // @param off Initial values of beta.
  // @param off Beta matrix to be updated.
  // @param eta Mixing parameter of half-cauchy distribution for each unit.
  // @param lambda Local shrinkage parameter common for each unit.
  // @param tau Global shrinkage parameter.
  matrix make_beta (int J,
                    matrix off,
                    vector lambda,
                    real eta,
                    vector tau) {
    int L = cols(off);
    vector[L] cache = (tan(0.5 * pi() * lambda) *
                       tan(0.5 * pi() * eta));

    vector[J] tau_ = tan(0.5 * pi() * tau);
    matrix[J, L] out;

    for (j in 1:J)
      out[j] = off[j] * tau_[j];

    return diag_pre_multiply(cache, out');
  }
}
