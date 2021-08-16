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
#' Converts data to wide format.
#' @param data Data.frame with the input data.
#' @param time Variable name for the time period.
#' @param id Variable name for the units.
#' @param oucome Name of the outcome variable.
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
