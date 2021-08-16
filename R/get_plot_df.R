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
#' @description
#' Returns data.frame ready for plotting with confidence intervals.
#' @param y_synth_draws Data.frame with each draw from the stan fit object.
#' @param pre_data Data.frame with data before the intervention.
#' @param post_data Data.frame with data after the intervention.
#' @param time Name of the time period variable.
#' @param oucome Name of the outcome variable.
#' @param ci Width of the credible confidence interval.
.get_plot_df <- function(y_synth_draws, pre_data,
                         post_data, time, outcome, ci = 0.75) {
  y_synth <- y_synth_draws %>%
    dplyr::group_by(!!time) %>%
    dplyr::summarise(
      LB = quantile(y_synth, (1 - ci) / 2),
      UB = quantile(y_synth, 1 - (1 - ci) / 2),
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

#' @description
#' Helper function to transform the data into a data.frame that can be
#' plotted for the case with more than one treated unit.
#' @param treated_ids Identifiers for the treated units.
#' @param data Data.frame with the input data.
#' @details other params same as get_plot_df() function.
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
      tau_LB = quantile(diff_draw, (1 - ci) / 2),
      tau_UB = quantile(diff_draw, 1 - (1 - ci) / 2)
    ) %>%
    dplyr::mutate(id = "Average")

  y_synth_i <- y_synth_draws %>%
    dplyr::group_by(!!time, !!id) %>%
    dplyr::summarise(
      LB = quantile(y_hat, (1 - ci) / 2),
      UB = quantile(y_hat, 1 - (1 - ci) / 2),
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
