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
#'
#' Returns Data Frame Ready for Plotting with Confidence Intervals
#'
#' This function processes data frames containing synthetic and observed outcomes,
#' calculates confidence intervals for the synthetic outcomes, and returns a
#' combined data frame suitable for plotting the results.
#'
#' @param y_synth_draws A data frame containing draws from the Stan fit object.
#' @param pre_data A data frame with data before the intervention.
#' @param post_data A data frame with data after the intervention.
#' @param time The name of the time period variable (as a string).
#' @param outcome The name of the outcome variable (as a string).
#' @param ci The width of the credible confidence interval (default: 0.75).
#'
#' @return A data frame containing:
#'   * `time`: The time period.
#'   * `outcome`: The observed outcome.
#'   * `y_synth`: The mean synthetic outcome.
#'   * `LB`: The lower bound of the confidence interval for the synthetic outcome.
#'   * `UB`: The upper bound of the confidence interval for the synthetic outcome.
#'   * `tau`: The difference between the observed and synthetic outcomes.
#'   * `tau_LB`: The lower bound of the confidence interval for `tau`.
#'   * `tau_UB`: The upper bound of the confidence interval for `tau`.
#'
.get_plot_df <- function(y_synth_draws, pre_data,
                         post_data, time, outcome, ci = 0.75) {
  y_synth <- y_synth_draws %>%
    dplyr::group_by(!!time) %>%
    dplyr::summarise(
      LB = stats::quantile(y_synth, (1 - ci) / 2),
      UB = stats::quantile(y_synth, 1 - (1 - ci) / 2),
      y_synth = mean(y_synth)
    )

  all_data <- dplyr::bind_rows(pre_data, post_data)

  df_plot_all <- dplyr::inner_join(y_synth, all_data,
    by = rlang::as_name(time)
  ) %>%
    dplyr::mutate(
      tau = !!outcome - y_synth,
      tau_LB = !!outcome - UB,
      tau_UB = !!outcome - LB
    ) %>%
    dplyr::select(!!time, !!outcome, y_synth, LB, UB, tau, tau_LB, tau_UB)
  return(df_plot_all)
}

#' Prepare Data Frame for Plotting with Multiple Treated Units
#'
#' This function processes data for multiple treated units, calculating synthetic outcomes,
#' confidence intervals, and treatment effects. It combines this information into a data
#' frame suitable for plotting the results.
#'
#' @param y_synth_draws A data frame containing synthetic outcome draws for each treated unit and time period.
#' @param data A data frame with the original data, including outcomes for treated units.
#' @param treated_ids A vector of identifiers for the treated units.
#' @param id The name of the variable in `data` that identifies units (as a string).
#' @param time The name of the time period variable (as a string).
#' @param outcome The name of the outcome variable (as a string).
#' @param ci The width of the credible confidence interval (default: 0.75).
#'
#' @return A data frame containing:
#'   * `time`: The time period.
#'   * `id`: The unit identifier (including "Average" for the average treatment effect).
#'   * `outcome`: The observed outcome (for treated units).
#'   * `y_synth`: The mean synthetic outcome (for treated units and the average).
#'   * `LB`: The lower bound of the confidence interval for the synthetic outcome.
#'   * `UB`: The upper bound of the confidence interval for the synthetic outcome.
#'   * `tau`: The treatment effect (difference between observed and synthetic outcomes).
#'   * `tau_LB`: The lower bound of the confidence interval for the treatment effect.
#'   * `tau_UB`: The upper bound of the confidence interval for the treatment effect.
#'
.get_plot_df2 <- function(y_synth_draws, data, treated_ids,
                          id, time, outcome, ci = 0.75) {
  data_treated <- data %>%
    dplyr::filter(!!id %in% treated_ids) %>%
    dplyr::select(!!id, !!time, !!outcome) %>%
    dplyr::mutate(!!rlang::as_label(id) := as.character(!!id))

  ate <- data_treated %>%
    dplyr::inner_join(y_synth_draws, by = c(
      rlang::as_name(id),
      rlang::as_name(time)
    )) %>%
    dplyr::mutate(diff = !!outcome - y_hat) %>%
    dplyr::group_by(!!time, draw) %>%
    dplyr::summarise(diff_draw = mean(diff)) %>%
    dplyr::group_by(!!time) %>%
    dplyr::summarise(
      tau = mean(diff_draw),
      tau_LB = stats::quantile(diff_draw, (1 - ci) / 2),
      tau_UB = stats::quantile(diff_draw, 1 - (1 - ci) / 2)
    ) %>%
    dplyr::mutate(id = "Average")

  y_synth_i <- y_synth_draws %>%
    dplyr::group_by(!!time, !!id) %>%
    dplyr::summarise(
      LB = stats::quantile(y_hat, (1 - ci) / 2),
      UB = stats::quantile(y_hat, 1 - (1 - ci) / 2),
      y_synth = mean(y_hat)
    )

  df_plot_i <- y_synth_i %>%
    dplyr::inner_join(data_treated, by = c(
      rlang::as_name(id),
      rlang::as_name(time)
    )) %>%
    dplyr::mutate(
      tau = !!outcome - y_synth,
      tau_LB = !!outcome - UB,
      tau_UB = !!outcome - LB
    )
  df_plot <- dplyr::bind_rows(df_plot_i, ate)
  return(df_plot)
}
