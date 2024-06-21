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
#' Plot Treatment Effect Estimate
#'
#' This internal helper function creates a plot to visualize the estimated
#' treatment effect over time. It allows for faceting by a specified variable
#' and optional subsetting of units to include in the plot.
#'
#' @param data A data frame containing the data to be plotted.
#' @param x The name of the x-axis variable (typically the time period) (as a string).
#' @param y The name of the y-axis variable (typically the treatment effect) (as a string).
#' @param ymin The name of the variable containing the lower bound of the confidence interval (as a string).
#' @param ymax The name of the variable containing the upper bound of the confidence interval (as a string).
#' @param xintercept The time point of the intervention to be marked with a vertical dashed line.
#' @param facet (Optional) The name of the variable to facet the plot by (as a string).
#' @param id The name of the variable identifying the units (as a string).
#' @param subset (Optional) A vector specifying a subset of units to include in the plot. If NULL, all units are included.
#'
#' @return A ggplot object displaying the treatment effect plot.
#'
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
