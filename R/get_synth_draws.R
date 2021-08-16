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
#' @description
#' Helper function (not exported) to get the draws as a tidy data frame when
#'   you only have a single treated unit.
#' @param fit Stan object with the fitted model.
#' @param pre_data Data.frame with data before the intervention.
#' @param post_data Data.frame with data after the intervention.
#' @param time Name of the time period variable.
#' @param oucome Name of the outcome variable.
.get_synth_draws <- function(fit, pre_data, post_data, time, outcome) {
  y_sim_draws <- .get_par_long(fit = fit, par = y_sim)
  dateXwalk <- pre_data %>%
    dplyr::mutate(idx = 1:dplyr::n()) %>%
    dplyr::select(idx, !!time)
  y_hat <- dplyr::inner_join(y_sim_draws, dateXwalk, by = "idx") %>%
    dplyr::rename(y_synth = y_sim)

  y_pred_draws <- .get_par_long(fit = fit, par = y_pred)

  dateXwalk <- post_data %>%
    dplyr::mutate(idx = 1:dplyr::n()) %>%
    dplyr::select(idx, !!time)
  y_pred_hat <- dplyr::inner_join(y_pred_draws, dateXwalk, by = "idx") %>%
    dplyr::rename(y_synth = y_pred)
  y_synth <- dplyr::bind_rows(y_hat, y_pred_hat) %>%
    dplyr::select(-idx)

  pre_outcome <- pre_data %>%
    dplyr::select(!!outcome, !!time)
  post_outcome <- post_data %>%
    dplyr::select(!!outcome, !!time)

  y <- dplyr::bind_rows(pre_outcome, post_outcome)
  y_synth <- dplyr::full_join(y_synth, y, by = rlang::as_name(time))
  return(y_synth)
}

# TODO(jvives): Unify get_synth_draws functions into one.
#' @description
#' Helper function (not exported) to get the draws as a tidy data frame when
#'   you only have a single treated unit.
#' This function uses the variable definitions of the predictor_match
#  stan model. Otherwise it is the same as get_synth_draws.
#' @details params are the same as in get_synth_draws().
.get_synth_draws_predictor_match <- function(fit, pre_data,
                                             post_data, time, outcome) {
  X1_sim_draws <- .get_par_long(fit = fit, par = X1_sim)
  dateXwalk <- pre_data %>%
    dplyr::mutate(idx = 1:dplyr::n()) %>%
    dplyr::select(idx, !!time)
  y_hat <- dplyr::inner_join(X1_sim_draws, dateXwalk, by = "idx") %>%
    dplyr::rename(y_synth = X1_sim)

  X1_pred_draws <- .get_par_long(fit = fit, par = X1_pred)
  dateXwalk <- post_data %>%
    dplyr::mutate(idx = 1:dplyr::n()) %>%
    dplyr::select(idx, !!time)
  X1_pred_hat <- dplyr::inner_join(X1_pred_draws, dateXwalk, by = "idx") %>%
    dplyr::rename(y_synth = X1_pred)

  y_synth <- dplyr::bind_rows(y_hat, X1_pred_hat) %>%
    dplyr::select(-idx)

  pre_outcome <- pre_data %>%
    dplyr::select(!!outcome, !!time)
  post_outcome <- post_data %>%
    dplyr::select(!!outcome, !!time)

  y <- dplyr::bind_rows(pre_outcome, post_outcome)
  y_synth <- dplyr::full_join(y_synth, y, by = rlang::as_name(time))
  return(y_synth)
}

#' @description
#' Helper function (not exported) to get the draws as a tidy data frame when
#'   you have multiple treated units.
#' @param data Data.frame with the input data.
#' @param intervention Name of the variable that identifies the intervention
#' time.
#' @param treated_ids Name of the variable that identifies the treated units.
#' @details other params are the same as in get_synth_draws().
.get_synth_draws3d <- function(fit, data, id, treated_ids, time, outcome,
                               intervention) {
  y_sim_draws <-
    .get_draws3d(
      fit = fit,
      data = data,
      id = id,
      treated_ids = treated_ids,
      time = time,
      outcome = outcome,
      intervention = intervention,
      period = "pre"
    )

  y_pred_draws <-
    .get_draws3d(
      fit = fit,
      data = data,
      id = id,
      treated_ids = treated_ids,
      time = time,
      outcome = outcome,
      intervention = intervention,
      period = "post"
    )

  y_draws <- dplyr::bind_rows(y_sim_draws, y_pred_draws)

  return(y_draws)
}

.get_draws3d <- function(fit, data, id, treated_ids, time, outcome,
                         intervention, period = c("pre", "post")) {
  period <- match.arg(period)
  if (period == "pre") {
    y_sim_draws <- rstan::extract(fit,
      pars = "y_sim"
    )[[1]]
    data <- data %>%
      dplyr::filter(!!time < intervention)
  } else {
    y_sim_draws <- rstan::extract(fit,
      pars = "y_pred"
    )[[1]]
    data <- data %>%
      dplyr::filter(!!time >= intervention)
  }

  dimnames(y_sim_draws) <- list(
    "draw" = seq(1, dim(y_sim_draws)[[1]]),
    "i_idx" = seq(1, dim(y_sim_draws)[[2]]),
    "t_idx" = seq(1, dim(y_sim_draws)[[3]])
  )

  y_sim_draws <- y_sim_draws %>%
    cubelyr::as.tbl_cube() %>%
    tibble::as_tibble() %>%
    dplyr::rename(y_hat = 4)

  wide_df_treated <- data %>%
    dplyr::filter(!!id %in% treated_ids) %>%
    dplyr::select(
      !!id,
      !!time,
      !!outcome
    ) %>%
    tidyr::pivot_wider(names_from = !!id, values_from = !!outcome) %>%
    dplyr::arrange(time)


  iXwalk <- wide_df_treated %>%
    dplyr::select(-!!time) %>%
    colnames() %>%
    dplyr::tibble(id = .) %>%
    dplyr::mutate(i_idx = 1:dplyr::n())

  tXwalk <- wide_df_treated %>%
    dplyr::select(!!time) %>%
    dplyr::mutate(t_idx = 1:dplyr::n())

  y_sim_draws <- y_sim_draws %>%
    dplyr::inner_join(iXwalk, by = "i_idx") %>%
    dplyr::inner_join(tXwalk, by = "t_idx") %>%
    dplyr::select(-i_idx, -t_idx)

  return(y_sim_draws)
}
