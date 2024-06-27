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
#' Time Tiles Plot of Intervention Impact
#'
#' This function creates a time tiles plot visualizing when and which units are
#' affected by an intervention. Each tile represents a unit at a specific time point,
#' with the color indicating the treatment status.
#'
#' @param data A data frame containing the input data.
#' @param time The name of the time period variable (as a string).
#' @param id The name of the unit identifier variable (as a string).
#' @param status The name of the variable that identifies the treatment status (as a string).
#'
#' @return A ggplot object displaying the time tiles plot.
#'
#' @export
time_tiles <- function(data, time, id, status) {
  tiles_plot <- ggplot2::ggplot(data, ggplot2::aes(
    x = {{ time }},
    y = {{ id }},
    fill = {{ status }}
  )) +
    ggplot2::geom_tile(color = "white", linewidth = 1) +
    ggplot2::scale_y_discrete(limits = rev) +
    ggplot2::scale_fill_manual(values = c("#4285F4", "#F4B400", "#DB4437")) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(
      legend.position = "bottom",
      panel.border = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      legend.title = ggplot2::element_blank()
    )
  return(tiles_plot)
}
