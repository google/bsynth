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
#' Get Parameter Estimates in Long Format
#'
#' Helper function to get the long dataset of draws given a stan fit object.
#'
#' @param fit Stan object with the fitted model.
#' @param par Variable to do the long table for.expand_more
#'
#' @return A tibble containing the parameter estimates in long format.
#'
.get_par_long <- function(fit, par) {
  par <- rlang::enquo(par)
  long_tlb <- fit %>%
    as.data.frame(pars = rlang::quo_text(par)) %>%
    dplyr::mutate(draw = 1:dplyr::n()) %>%
    tidyr::pivot_longer(
      names_to = "idx",
      values_to = rlang::quo_text(par), -draw
    ) %>%
    dplyr::mutate(idx = as.numeric(gsub(
      glue::glue("{rlang::quo_text(par)}\\[(\\d+)\\]"),
      "\\1",
      idx
    ))) %>%
    dplyr::as_tibble()
  return(long_tlb)
}
