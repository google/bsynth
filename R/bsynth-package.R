## Copyright 2021 Google LLC
##
## Licensed under the Apache License, Version 2.0 (the "License");
## you may not use this file except in compliance with the License.
## You may obtain a copy of the License at
##
##     https://www.apache.org/licenses/LICENSE-2.0
##
## Unless required by applicable law or agreed to in writing, software
## distributed under the License is distributed on an "AS IS" BASIS,
## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
## See the License for the specific language governing permissions and
## limitations under the License.
##
#' The 'bsynth' package.
#'
#' @description Provides causal inference with a Bayesian synthetic
#'     control method.
#'

#' @name bsynth-package
#' @aliases bsynth
#' @useDynLib bsynth, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling
#' @importFrom rlang :=
#'
#' @references
#' Stan Development Team (2020). RStan: the R interface to Stan. R package version 2.21.2. https://mc-stan.org
#'
"_PACKAGE"
globalVariables(c(".", "i_idx", "t_idx", "draw", "idx", "UB", "LB", "tau",
                  "tau_LB", "tau_UB", "y_hat", "diff_draw", "y_synth",
                  "y_sim", "y_pred", "X1_sim", "X1_pred", "status", ".tmp_id"))
