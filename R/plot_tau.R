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
#' Plots the treatment effect estimate given the model.
#' @param data Data.frame with the input data.
#' @param x Name of the x axis variable (e.g. time period).
#' @param y Name of the y axis (e.g. treatement effect).
#' @param ymin, ymax Range values of the y variable.
#' @param xintercept Value of the time of the intervention to plot dashed line.
#' @param facet Variable to split the plots by.
#' @param id Variable name of the units.
#' @param subset Set of units to use for plot, all if NULL.
.plot_tau <- function(data, x, y, ymin, ymax, xintercept,
                      facet, id, subset = NULL) {
  if (!is.null(subset)) {
    data <- data %>%
      dplyr::filter(!!id %in% subset)
  }
  tau_plot <- ggplot2::ggplot(data = data, ggplot2::aes(x = !!x)) +
    ggplot2::geom_line(ggplot2::aes(y = {{ y }})) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = {{ ymin }}, ymax = {{ ymax }}),
      color = "gray",
      alpha = 0.2
    ) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(
      legend.position = "none",
      panel.border = ggplot2::element_blank(),
      axis.line = ggplot2::element_line()
    ) +
    ggplot2::geom_vline(xintercept = xintercept, linetype = "dashed")

  if (!missing(facet)) {
    tau_plot <- tau_plot + ggplot2::facet_grid(cols = dplyr::vars(!!facet))
  }
  return(tau_plot)
}
