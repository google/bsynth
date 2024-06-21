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
#' Convert Data to Wide Format
#'
#' This internal helper function transforms data from a long format, where each
#' row represents an observation for a specific unit and time, to a wide format,
#' where each row represents a time period and each column represents a unit's outcome.
#' It specifically focuses on separating treated and untreated units.
#'
#' @param data A data frame containing the input data.
#' @param id The name of the variable in `data` that identifies units (as a string).
#' @param time The name of the time period variable (as a string).
#' @param outcome The name of the outcome variable (as a string).
#' @param treatment The name of the variable in `data` that indicates treatment status (as a string).
#'
#' @return A data frame in wide format, where each row corresponds to a time period,
#'   and columns include the time variable, the treatment indicator, and the outcome
#'   values for each treated unit and all untreated units.
#'
.makeWide <- function(data, id, time, outcome, treatment) {
  data <- data %>%
    dplyr::mutate(.tmp_id = as.integer(as.factor(!!id)))

  treated <- data %>%
    dplyr::filter(status == "Treated") %>%
    dplyr::select(.tmp_id) %>%
    dplyr::distinct() %>%
    dplyr::pull(.tmp_id)

  wide_df_treated <- data %>%
    dplyr::filter(.tmp_id %in% treated) %>%
    dplyr::select(!!time, !!treatment, !!outcome)

  wide_df_untreated <- data %>%
    dplyr::filter(!(.tmp_id %in% treated)) %>%
    dplyr::select(!!time, !!outcome, !!id) %>%
    tidyr::pivot_wider(
      names_from = !!id,
      values_from = !!outcome
    )

  dplyr::inner_join(wide_df_treated, wide_df_untreated,
    by = rlang::as_name(time)
  ) %>% return()
}
