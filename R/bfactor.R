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
#' Create a Bayesian Synthetic Control Object Using Panel Data
#'
#' @description
#'
#' A Bayesian Factor Model has raw data and draws from the posterior
#' distribution. This is represented by an R6 Class.
#'
#' Code and theory based on
#' [Pinkney 2021](https://arxiv.org/abs/2103.16244).
#'
#' public methods:
#'
#' * `initialize()` initializes the variables and model parameters
#' * `fit()`  fits the stan model and returns a fit object
#' * `updateWidth` updates the width of the credible interval
#' * `placeboPlot` generates a counterfactual placebo plot
#' * `effectPlot` returns a plot of the treatment effect over time
#' * `summarizeLift`returns descriptive statistics of the lift estimate
#' * `biasDraws` returns a plot of the relative bias in a LFM
#' * `liftDraws` returns a plot of the posterior lift distribution
#' * `liftBias` returns a plot of the relative bias given a lift offset
#'
#' @param time Name of the variable in the data frame that
# `   identifies the time period (e.g. year, month, week etc).
#' @param id Name of the variable in the data frame that
#'   identifies the units (e.g. country, region etc).
#' @param treated Name of the variable in the data frame that contains the
#'   treatment assignment of the intervention.
#' @param data Long data.frame object with fields outcome, time, id,
#'   and treatment indicator.
#' @param outcome Name of the outcome variable.
#' @param ci_width Credible interval's width.  This number is in the
#'   (0,1) interval.
#' @param covariates Dataframe with a column for id and the other columns
#'   Defaults to NULL if no covariates should be included in the model.
#' @export
bayesianFactor <- R6::R6Class(
  classname = "bayesianFactor",
  private = list(
    time = NULL,
    id = NULL,
    treated = NULL,
    data = NULL,
    treated_ids = NULL,
    time_tiles = NULL,
    intervention = NULL,
    outcome = NULL,
    ci_width = NULL,
    stan_model = NULL,
    covariates = NULL,
    fitted = NULL,
    plot_data = NULL,
    y_synth_draws = NULL
  ),
  active = list(
    #' @field timeTiles ggplot2 object that shows when
    #'   the intervention happened.
    timeTiles = function() {
      return(private$time_tiles)
    },

    #' @field plotData tibble with the observed outcome and the
    #'   counterfactual data.
    plotData = function() {
      return(private$plot_data)
    },

    #' @field interventionTime returns the intervention time period.
    interventionTime = function() {
      return(private$intervention)
    },

    #' @field synthetic ggplot2 object that shows the
    #'   observed and counterfactual outcomes over time.
    synthetic = function() {
      df_plot <- private$plot_data %>%
        dplyr::rename(
          Observed = !!private$outcome,
          Synthetic = y_synth
        ) %>%
        dplyr::select(!!private$time, Observed, Synthetic, LB, UB) %>%
        tidyr::pivot_longer(cols = c(Observed, Synthetic))

      synthetic_plot <- ggplot2::ggplot(
        data = df_plot,
        ggplot2::aes(x = !!private$time)
        ) +
        ggplot2::geom_line(ggplot2::aes(y = value, linetype = name)) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = LB, ymax = UB),
          color = "gray",
          alpha = 0.2
        ) +
        ggplot2::theme_bw(base_size = 14) +
        ggplot2::theme(
          legend.title = ggplot2::element_blank(),
          legend.position = c(0.9, .1),
          legend.background = ggplot2::element_rect(
            fill =
              ggplot2::alpha("white", 0)
          ),
          panel.border = ggplot2::element_blank(),
          axis.line = ggplot2::element_line()
        ) +
        ggplot2::geom_vline(
          xintercept = private$intervention,
          linetype = "dashed"
        )

      return(synthetic_plot)
    }
  ),
  public = list(
    #' @description
    #' Create a new bayesianFactor object.
    #' @details params described in the data structure section of the
    #'    documentation of the R6 class at the top of the file.
    #' @return A new `bayesianFactor` object.
    initialize = function(data,
                          time,
                          id,
                          treated,
                          outcome,
                          ci_width = 0.75,
                          covariates) {
      stopifnot((ci_width > 0 & ci_width < 1))
      private$time <- rlang::enquo(time)
      private$id <- rlang::enquo(id)
      private$treated <- rlang::enquo(treated)
      if (!missing(covariates)) {
        # TODO(jvives): ADD SOME CHECKS (no missing data, sorted, etc)
        # Check that covariates contain no missing values.
        if (covariates %>% dplyr::select(-{{id}}) %>% any(is.na())) {
          stop("Covariates have missing data.")
        }

        # Check that covariates are ordered the same way as input data.
        # if (covariates %>% dplyr::select({{id}})
        # stopif(covariates %>%
        #  dplyr::select(-{{id}}) %>%
        #  any(is.na()))

        private$stan_model <- stanmodels$factor_model_with_covariates
      } else {
        private$stan_model <- stanmodels$factor_model_without_covariates
      }

      private$data <- data %>%
        dplyr::mutate(
          status = dplyr::case_when(
            {{ treated }} == 1 ~ "Treated",
            {{ treated }} == 0 ~ "Untreated",
            is.na({{ treated }}) ~ "N/A"
          ),
          status = factor(status, levels = c(
            "Treated",
            "Untreated",
            "N/A"
          ))
        )
      if (!(data %>% dplyr::pull({{ id }}) %>% is.factor())) {
        private$data <- private$data %>%
          dplyr::mutate(!!rlang::quo_name(private$id) := factor({{ id }}))
      }

      private$treated_ids <- private$data %>%
        dplyr::filter(!!private$treated == 1) %>%
        dplyr::select({{ id }}) %>%
        dplyr::distinct() %>%
        dplyr::pull({{ id }})

      private$time_tiles <- time_tiles(
        data = private$data,
        time = !!private$time,
        id = !!private$id,
        status = status
      )

      private$intervention <- private$data %>%
        dplyr::filter(status == "Treated") %>%
        dplyr::summarise(!!rlang::as_label(private$time) := min(!!private$time)) %>%
        dplyr::pull(!!private$time)

      private$outcome <- rlang::enquo(outcome)
      private$updateWidth(ci_width)

      # TODO(jvives): Add support to handle more than one treated unit.
      if (length(private$treated_ids) != 1) {
        stop("Support for more than one treated unit is not yet avaialble.")
      }
    },

    #' @description
    #' Fit Stan model.
    #' @param L Number of factors.
    #' @param ... other arguments passed to [rstan::sampling()].
    fit = function(L = 8, ...) {
      y_treated_pre <- private$data %>%
        dplyr::filter(
          !!private$time < private$intervention,
          !!private$id %in% private$treated_ids
        ) %>%
        dplyr::pull(!!private$outcome)

      y_donors_pre <- private$data %>%
        dplyr::filter(
          !!private$time < private$intervention,
          !(!!private$id %in% private$treated_ids)
        ) %>%
        dplyr::select(!!private$id, !!private$outcome, !!private$time) %>%
        tidyr::pivot_wider(
          names_from = !!private$id,
          values_from = !!private$outcome
        ) %>%
        dplyr::arrange(!!private$time) %>%
        dplyr::select(-!!private$time) %>%
        dplyr::select(sort(tidyselect::peek_vars())) %>%
        as.matrix() %>%
        t()

      y_donors_post <- private$data %>%
        dplyr::filter(
          !!private$time >= private$intervention,
          !(!!private$id %in% private$treated_ids)
        ) %>%
        dplyr::select(!!private$id, !!private$outcome, !!private$time) %>%
        tidyr::pivot_wider(
          names_from = !!private$id,
          values_from = !!private$outcome
        ) %>%
        dplyr::arrange(!!private$time) %>%
        dplyr::select(-!!private$time) %>%
        dplyr::select(sort(tidyselect::peek_vars())) %>%
        as.matrix() %>%
        t()

      if (!is.null(private$covariates)) {
        stan_data <- list(
          L = L,
          N = length(y_treated_pre),
          y_treated_pre = y_treated_pre,
          J = nrow(y_donors_pre),
          y_donors_pre = y_donors_pre,
          N_pred = ncol(y_donors_post),
          y_donors_post = y_donors_post,
          P = nrow(private$covariates),
          X = private$covariates
        )
      } else {
        stan_data <- list(
          L = L,
          N = length(y_treated_pre),
          y_treated_pre = y_treated_pre,
          J = nrow(y_donors_pre),
          y_donors_pre = y_donors_pre,
          N_pred = ncol(y_donors_post),
          y_donors_post = y_donors_post
        )
      }

      private$fitted <-
        rstan::sampling(private$stan_model,
          data = stan_data,
          ...
        )

      timeXwalk <- private$data %>%
        dplyr::select(!!private$time) %>%
        dplyr::distinct() %>%
        dplyr::arrange(!!private$time) %>%
        dplyr::mutate(idx = seq_len(dplyr::n()))

      private$y_synth_draws <-
        .get_par_long(fit = private$fitted, par = synth_out) %>%
        dplyr::inner_join(timeXwalk, by = "idx") %>%
        dplyr::rename(y_synth = synth_out)

      # TODO(jvives): Rename synth_out to y_sim.
      # private$fitted <- private$fitted %>%
      #  dplyr::rename_all(
      #    funs(stringr::str_replace_all(., "synth_out", "y_sim"))
      #  )

      wide_df <- .makeWide(
        data = private$data,
        id = private$id,
        time = private$time,
        outcome = private$outcome,
        treatment = private$treated
      )

      pre_data <- wide_df %>%
        dplyr::filter(!!private$time < private$intervention)

      post_data <- wide_df %>%
        dplyr::filter(!!private$time >= private$intervention)

      # Add outcome variable to y_synth_draws.
      # TODO(jvives): Integrate this step in get_synth_draws function.
      pre_outcome <- pre_data %>%
        dplyr::select(!!private$outcome, !!private$time)
      post_outcome <- post_data %>%
        dplyr::select(!!private$outcome, !!private$time)

      y <- dplyr::bind_rows(pre_outcome, post_outcome)
      private$y_synth_draws <- dplyr::full_join(
        private$y_synth_draws, y,
        by = rlang::as_name(private$time))

      private$plot_data <-
        .get_plot_df(
          y_synth_draws = private$y_synth_draws,
          pre_data = pre_data,
          post_data = post_data,
          ci = private$ci_width,
          time = private$time,
          outcome = private$outcome
        )
    },

    #' @description
    #' Update the width of the credible interval.
    #' @param ci_width New width for the credible interval. This number
    #'     should be in the (0,1) interval.
    updateWidth = function(ci_width = 0.75) {
      stopifnot((ci_width > 0 & ci_width < 1))
      private$ci_width <- ci_width

      wide_df <- .makeWide(
        data = private$data,
        id = private$id,
        time = private$time,
        outcome = private$outcome,
        treatment = private$treated
      )

      pre_data <- wide_df %>%
        dplyr::filter(!!private$time < private$intervention)

      post_data <- wide_df %>%
        dplyr::filter(!!private$time >= private$intervention)

      private$plot_data <-
        .get_plot_df(
          y_synth_draws = private$y_synth_draws,
          pre_data = pre_data,
          post_data = post_data,
          ci = private$ci_width,
          time = private$time,
          outcome = private$outcome
        )
    },

    #' @description summarizeLift returns descriptive statistics of
    #' the lift estimate.
    summarizeLift = function() {
      if (is.null(private$lift_draws)) {
        stop("You first need to run the `liftDraws()` method.")
      }

      return({
        c(
          point = mean(private$lift_draws$lift),
          lower_bound = quantile(
            private$lift_draws$lift,
            (1 - private$ci_width) / 2
          ),
          upper_bound = quantile(
            private$lift_draws$lift,
            1 - (1 - private$ci_width) / 2
          )
        )
      })
    },

    #' @description effectPlot returns ggplot2 object that shows the
    #'   effect of the intervention over time.
    effectPlot = function() {
      tau_plot <- .plot_tau(
        data = private$plot_data,
        x = private$time,
        y = tau,
        ymin = tau_LB,
        ymax = tau_UB,
        xintercept = private$intervention
      )

      return(tau_plot)
    },

    #' @description
    #' Plots lift.
    #' @param from First period to consider when calculating lift. If infinite,
    #' set to the time of the intervention.
    #' @param to Last period to consider when calculating lift. If infinite, set
    #' to the last period.
    #' @param ... other arguments passed to vizdraws::vizdraws().
    #' @return vizdraws object with the posterior distribution of the lift.
    liftDraws = function(from, to, ...) {
      if (inherits(from, "Date")) {
        if (is.infinite(from)) {
          from <- private$intervention
        }
        if (is.infinite(to)) {
          to <- private$data %>%
            pull(!!private$time) %>%
            max()
        }
      }
      # TODO(jvives): Add check for small values of sum_y0.
      lift <- private$y_synth_draws %>%
        dplyr::filter(!!private$time >= from, !!private$time <= to) %>%
        dplyr::group_by(draw) %>%
        dplyr::summarise(
          sum_y0 = sum(y_synth),
          sum_y1 = sum(!!private$outcome),
          lift = (sum_y1 - sum_y0) / sum_y0
        ) %>% pull(lift)

      vizdraws::vizdraws(
        posterior = lift,
        percentage = TRUE,
        quantity = TRUE,
        tense = "past",
        display_mode_name = TRUE,
        title = "Contrafactual Lift",
        xlab = "Effect of the Intervention",
        units = "the contrafacutal lift",
        ...
      ) %>% return()
    },

    #' @description
    #' Plot bias magnitude in terms of lift for period (firstT, lastT)
    #' @param firstT Start of the time period to compute relative bias over.
    #' Must be after the intervention.
    #' @param lastT End of the time period to compute relative bias over.
    #' Must be after the intervention.
    #'     over. They must be after the intervention.
    #' @param offset Target lift %.
    #' @param ... other arguments passed to vizdraws::vizdraws().
    #' @returns vizdraws object with the relative bias with offset.
    liftBias = function(firstT, lastT, offset, ...) {
      # Calculate the gaps.
      # TODO(jvives): Add capability for this to work with covariates.
      mad_pre_gaps <- private$y_synth_draws %>%
        dplyr::filter(!!private$time < private$intervention) %>%
        dplyr::mutate(gap = y_synth - !!private$outcome) %>%
        dplyr::group_by(draw) %>%
        dplyr::summarise(mad = stats::mad(gap))

      # Get all time periods in TL.
      Ts <- private$y_synth_draws %>%
        dplyr::filter(!!private$time <= lastT, !!private$time >= firstT) %>%
        dplyr::select(!!private$time) %>%
        unique() %>%
        dplyr::pull()

      # Get time periods between intervention and first_T.
      Ts_from_T0 <- private$y_synth_draws %>%
        dplyr::filter(!!private$time <= firstT,
                      !!private$time >= private$intervention) %>%
        dplyr::select(!!private$time) %>%
        dplyr::n_distinct() %>%
        as.numeric()
      T_from_T0 <- Ts_from_T0 + 1

      # Calculate the estimated treatment effects (taus).
      # TODO(jvives): Consolidate functions into one.
      taus_TL <- private$y_synth_draws %>%
        dplyr::filter(!!private$time >= private$intervention) %>%
        dplyr::filter(!!private$time == Ts) %>%
        dplyr::group_by(draw) %>%
        dplyr::mutate(gapTauTL = abs(sum(y_synth - !!private$outcome))) %>%
        dplyr::mutate(sumTauTL = abs(sum(!!private$outcome)))

      # Divide by first period tau. For now just first period.
      TL <- length(Ts)
      if (TL >= 2) {
        rel_bias <- dplyr::inner_join(mad_pre_gaps, taus_TL, by = "draw") %>%
          dplyr::mutate(
            bias_ratio = abs(T_from_T0 * 0.5 * TL * (TL - 1) * mad - sumTauTL *
                             offset) / gapTauTL) %>%
          dplyr::filter(bias_ratio < 2.)
      } else {
        rel_bias <- dplyr::inner_join(mad_pre_gaps, taus_TL, by = "draw") %>%
          dplyr::mutate(
            bias_ratio = abs(T_from_T0 * mad - sumTauTL * offset) / gapTauTL) %>%
          dplyr::filter(bias_ratio < 2.)
      }

      # Return the relative bias vizdraw plot.
      biasDrawsViz <- vizdraws::vizdraws(
        posterior = rel_bias$bias_ratio,
        percentage = TRUE,
        quantity = TRUE,
        tense = "past",
        display_mode_name = TRUE,
        title = sprintf(
          "Upper bias changing a %#.2f %% Lift",
          offset * 100
        ),
        xlab = "Bias Relative to Lift",
        units = "Bias Relative to Lift",
        breaks = c(0.3, 1.),
        break_names = c(
          "Unlikely",
          "Close to change",
          "Changes bucket"
        ),
        colors = c("#4daf4a", "#377eb8", "#e41a1c"),
        ...
      )
      return(biasDrawsViz)
    },

    #' @description
    #' Plots relative upper bias / tau for a time period (firstT, lastT).
    #' @param small_bias Threshold value for considering the bias "small".
    #' @param firstT,lastT Time periods to compute relative bias over, they must
    #'     after the intervention.
    #' @return vizdraw object with the posterior distribution of relative bias.
    #' Bias is scaled by the time periods.
    biasDraws = function(small_bias = 0.3, firstT, lastT) {
      # Calculate the gaps.
      # TODO(jvives): Add capability for it to work with covariates.
      mad_pre_gaps <- private$y_synth_draws %>%
        dplyr::filter(!!private$time < private$intervention) %>%
        dplyr::mutate(gap = y_synth - !!private$outcome) %>%
        dplyr::group_by(draw) %>%
        dplyr::summarise(mad = stats::mad(gap))

      # Get all time periods in TL.
      Ts <- private$y_synth_draws %>%
        dplyr::filter(!!private$time <= lastT, !!private$time >= firstT) %>%
        dplyr::select(!!private$time) %>%
        unique() %>%
        dplyr::pull()

      # Get time periods between intervention and first_T.
      Ts_from_T0 <- private$y_synth_draws %>%
        dplyr::filter(
          !!private$time <= firstT,
          !!private$time >= private$intervention
        ) %>%
        dplyr::select(!!private$time) %>%
        dplyr::n_distinct() %>%
        as.numeric()
      T_from_T0 <- Ts_from_T0 + 1

      # Calculate the estimated treatment effects (taus).
      taus_TL <- private$y_synth_draws %>%
        dplyr::filter(!!private$time >= private$intervention) %>%
        dplyr::filter(!!private$time == Ts) %>%
        dplyr::group_by(draw) %>%
        dplyr::mutate(gapTauTL = abs(sum(y_synth - !!private$outcome)))

      # Divide by first period tau. For now just first period.
      TL <- length(Ts)
      if (TL >= 2) {
        rel_bias <- dplyr::inner_join(mad_pre_gaps, taus_TL, by = "draw") %>%
          dplyr::mutate(bias_ratio = T_from_T0 * 0.5 * TL * (TL - 1) * mad / gapTauTL) %>%
          dplyr::filter(bias_ratio < 5.)
      } else {
        rel_bias <- dplyr::inner_join(mad_pre_gaps, taus_TL, by = "draw") %>%
          dplyr::mutate(bias_ratio = T_from_T0 * mad / gapTauTL) %>%
          dplyr::filter(bias_ratio < 5.)
      }

      # Return the relative bias vizdraw plot.
      biasDrawsViz <- vizdraws::vizdraws(
        posterior = rel_bias$bias_ratio,
        percentage = TRUE,
        quantity = TRUE,
        tense = "past",
        display_mode_name = TRUE,
        title = "Counterfactual Upper Bias",
        xlab = "Relative Bias",
        units = "Relative Bias",
        breaks = c(small_bias, 1.),
        break_names = c(
          "Small Bias",
          "Some Bias",
          "Sign Change"
        ),
        colors = c("#4daf4a", "#377eb8", "#e41a1c")
      )
      return(biasDrawsViz)
    }
  ) # end public methods
) # end class
