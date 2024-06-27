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
#' A Bayesian Synthetic Control has raw data and draws from the posterior
#' distribution. This is represented by an R6 Class.
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
#' Data structure:
#' @param data Long data.frame object with fields outcome, time, id,
#'   and treatment indicator.
#' @param covariates Data.frame with time dependent covariates for
#'   for each unit and time field.
#'   Defaults to NULL if no covariates should be included in the model.
#' @param predictor_match_covariates0 data.frame with time independent
#'   covariates on each row and column indicating the control unit names
#'   (dim k x J+1).
#' @param predictor_match_covariates1 Vector with time independent
#'   covariates for the treated unit (dim k x 1).
#' @param predictor_match Logical that indicates whether or not to run
#'   the matching version of the Bayesian Synthetic Control. This option can
#'   not be used with gp, covariates or multiple treated units.
#' @param vs Vector of weights for the importance of the predictors used
#'   in creating the synthetic control.
#'   Defaults to equal weight for all predictors.
#' @param time Name of the variable in the data frame that
#'   identifies the time period (e.g. year, month, week etc).
#' @param id Name of the variable in the data frame that
#'   identifies the units (e.g. country, region etc).
#' @param treated Name of the variable in the data frame that contains
#'   the treatment assignment of the intervention.
#' @param outcome Name of the outcome variable.
#' @param ci_width Credible interval's width.  This number is in the
#' (0,1) interval.
#' @param intervention Intervention time period (e.g., year)
#'    in which the treatment occurred.
#' @param plot_data Tibble with the observed outcome and the
#'   counterfactual data.
#' @param time_Tiles ggplot2 object with a timeline for the intervention.
#' @export
bayesianSynth <- R6::R6Class(
  classname = "bayesianSynth",
  private = list(
    data = NULL,
    covariates = NULL,
    predictor_match_covariates0 = NULL,
    predictor_match_covariates1 = NULL,
    predictor_match = NULL,
    vs = NULL,
    time = NULL,
    id = NULL,
    treated = NULL,
    outcome = NULL,
    ci_width = NULL,
    intervention = NULL,
    fitted = NULL,
    plot_data = NULL,
    time_tiles = NULL,
    stan_model = NULL,
    y_synth_draws = NULL,
    lift_draws = NULL,
    treated_ids = NULL,
    mcmc_checks = NULL # a mcmc_checks object
  ),
  active = list(

    #' @field timeTiles ggplot2 object that shows when
    #'   the intervention happened.
    timeTiles = function() {
      return(private$time_tiles)
    },

    #' @field plotData returns tibble with the observed outcome and the
    #'   counterfactual data.
    plotData = function() {
      return(private$plot_data)
    },

    #' @field interventionTime returns intervention time period (e.g., year)
    #'   in which the treatment occurred.
    interventionTime = function() {
      return(private$intervention)
    },

    #' @field synthetic returns ggplot2 object that shows the
    #'   observed and counterfactual outcomes over time.
    synthetic = function() {
      df_plot <- private$plot_data %>%
        dplyr::rename(
          Observed = !!private$outcome,
          Synthetic = y_synth
        )

      if (length(private$treated_ids) == 1) {
        df_plot <- df_plot %>%
          dplyr::select(!!private$time, Observed, Synthetic, LB, UB)
      } else {
        df_plot <- df_plot %>%
          dplyr::select(
            !!private$id, !!private$time, Observed,
            Synthetic, LB, UB
          ) %>%
          dplyr::filter(!!private$id != "Average")
      }
      df_plot <- df_plot %>%
        tidyr::pivot_longer(cols = c(Observed, Synthetic))

      synthetic_plot <- ggplot2::ggplot(
        data = df_plot,
        ggplot2::aes(x = !!private$time)
      ) +
        ggplot2::geom_line(
          ggplot2::aes(y = value, linetype = name)
        ) +
        ggplot2::geom_ribbon(
          ggplot2::aes(ymin = LB, ymax = UB),
          color = "gray",
          alpha = 0.2
        ) +
        ggplot2::theme_bw(base_size = 14) +
        ggplot2::theme(
          legend.title = ggplot2::element_blank(),
          legend.position = c(0.9, .1),
          legend.background = ggplot2::element_rect(
            fill = ggplot2::alpha("white", 0)
          ),
          panel.border = ggplot2::element_blank(),
          axis.line = ggplot2::element_line()
        ) +
        ggplot2::geom_vline(
          xintercept = private$intervention,
          linetype = "dashed"
        )

      if (length(private$treated_ids) > 1) {
        synthetic_plot <- synthetic_plot +
          ggplot2::facet_wrap(dplyr::vars(!!private$id))
      }

      return(synthetic_plot)
    },

    #' @field checks returns MCMC checks.
    checks = function() {
      private$mcmc_checks
    },

    #' @field lift draws from the posterior distribution of the lift.
    lift = function() {
      private$lift_draws
    }
  ),
  public = list(
    #' @description
    #' Create a new bayesianSynth object.
    #' @param gp Logical that indicates whether or not to include a
    #'   Gaussian Process as part of the model.
    #' @return A new `bayesianSynth` object.
    initialize = function(data, time, id, treated, outcome, ci_width = 0.75,
                          gp = FALSE, covariates = NULL,
                          predictor_match = FALSE,
                          predictor_match_covariates0 = NULL,
                          predictor_match_covariates1 = NULL, vs = NULL) {
      stopifnot((ci_width > 0 & ci_width < 1))
      private$time <- rlang::enquo(time)
      private$id <- rlang::enquo(id)
      private$treated <- rlang::enquo(treated)
      private$predictor_match <- predictor_match
      private$predictor_match_covariates0 <- predictor_match_covariates0
      private$predictor_match_covariates1 <- predictor_match_covariates1
      private$vs <- vs # Predictor weights

      message("Transforming data")

      # Check treated has the right input.
      if (!setequal(
        data %>%
          dplyr::select(!!private$treated) %>% dplyr::distinct() %>%
          dplyr::pull({{ treated }}),
        c(1, 0)
      )) {
        stop("Treated identifier not binary 1 - 0.")
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
      private$ci_width <- ci_width

      if (length(private$treated_ids) == 1) {
        if (is.null(covariates)) {
          if (gp) {
            private$stan_model <- stanmodels$model2
          } else {
            if (predictor_match) {
              private$stan_model <- stanmodels$model1_gammaOmega
            } else {
              private$stan_model <- stanmodels$model1
            }
          }
        } else {
          if (gp) {
            private$stan_model <- stanmodels$model4
          } else {
            private$stan_model <- stanmodels$model3
          }
          private$covariates <- covariates %>%
            dplyr::arrange(!!private$time)
        }
      } else {
        if (is.null(covariates)) {
          if (gp) {
            private$stan_model <- stanmodels$model8
          } else {
            private$stan_model <- stanmodels$model5
          }
        } else {
          if (gp) {
            private$stan_model <- stanmodels$model7
          } else {
            private$stan_model <- stanmodels$model6
          }
          private$covariates <- covariates %>%
            dplyr::arrange(!!private$time)
        }
      }
    },

    #' @description
    #' Fit Stan model.
    #' @param ... other arguments passed to [rstan::sampling()].
    fit = function(...) {
      if (length(private$treated_ids) == 1) { ## Only one treated unit
        # TODO(jvives): Create utility function to make pre/post data frames.
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

        X <- pre_data %>%
          dplyr::select(-!!private$time, -!!private$treated, -!!private$outcome)
        X_pred <- post_data %>%
          dplyr::select(-!!private$time, -!!private$treated, -!!private$outcome)

        if (!is.null(private$covariates)) {
          M <- private$covariates %>%
            dplyr::filter(!!private$time < private$intervention)

          M_pred <- private$covariates %>%
            dplyr::filter(!!private$time >= private$intervention)

          stan_data <- list(
            N = nrow(pre_data),
            y = pre_data %>% dplyr::pull(!!private$outcome),
            K = ncol(X),
            X = X,
            M_K = ncol(M),
            M = M,
            N_pred = nrow(X_pred),
            X_pred = X_pred,
            M_pred = M_pred
          )
        } else {
          if (private$predictor_match) {
            X1 <- pre_data %>% dplyr::pull(!!private$outcome)
            if (!is.null(private$predictor_match_covariates0)) {
              # covariates are a k x J matrix
              # TODO(jvives): Make it general and add check for input order.
              X <- rbind(X, setNames(
                private$predictor_match_covariates0,
                names(X)
              ))
              X1 <- c(
                pre_data %>% dplyr::pull(!!private$outcome),
                private$predictor_match_covariates1
              )
            }
            # TODO(jvives): Move to initialize.
            # Set predictor weights if not given.
            if (is.null(private$vs)) {
              private$vs <- rep(1, nrow(X))
            }

            stan_data <- list(
              K = nrow(X),
              X1 = X1,
              J = ncol(X),
              X0 = X,
              T_post = nrow(X_pred),
              X0_pred = X_pred,
              vs = private$vs
            )
          } else {
            stan_data <- list(
              N = nrow(pre_data),
              y = pre_data %>% dplyr::pull(!!private$outcome),
              K = ncol(X),
              X = X,
              N_pred = nrow(X_pred),
              X_pred = X_pred
            )
          }
        }
      } else { ## More than 1 treated unit

        Y <- private$data %>%
          dplyr::filter(!!private$id %in% private$treated_ids) %>%
          dplyr::select(
            !!private$id,
            !!private$time,
            !!private$outcome
          ) %>%
          dplyr::filter(!!private$time < private$intervention) %>%
          tidyr::pivot_wider(
            names_from = !!private$id,
            values_from = !!private$outcome
          ) %>%
          dplyr::arrange(!!private$time) %>%
          dplyr::select(-!!private$time) %>%
          as.matrix() %>%
          t()

        donors <- private$data %>%
          dplyr::filter(!(!!private$id %in% private$treated_ids)) %>%
          dplyr::select(time, Y, id)

        X <- donors %>%
          dplyr::filter(!!private$time < !!private$intervention) %>%
          dplyr::arrange(!!private$time) %>%
          tidyr::pivot_wider(
            names_from = !!private$id,
            values_from = !!private$outcome
          ) %>%
          dplyr::select(-!!private$time)

        X_pred <- donors %>%
          dplyr::filter(!!private$time >= !!private$intervention) %>%
          dplyr::arrange(!!private$time) %>%
          tidyr::pivot_wider(
            names_from = !!private$id,
            values_from = !!private$outcome
          ) %>%
          dplyr::select(-!!private$time)

        if (!is.null(private$covariates)) {
          M <- private$covariates %>%
            dplyr::filter(
              !!private$id %in% private$treated_ids,
              !!private$time < private$intervention
            ) %>%
            dplyr::arrange(!!private$id, !!private$time) %>%
            dplyr::select(-!!private$time) %>%
            dplyr::group_by(!!private$id) %>%
            tidyr::nest() %>%
            dplyr::pull(data) %>%
            purrr::map(as.matrix)

          M_pred <- private$covariates %>%
            dplyr::filter(
              !!private$id %in% private$treated_ids,
              !!private$time >= private$intervention
            ) %>%
            dplyr::arrange(!!private$id, !!private$time) %>%
            dplyr::select(-!!private$time) %>%
            dplyr::group_by(!!private$id) %>%
            tidyr::nest() %>%
            dplyr::pull(data) %>%
            purrr::map(as.matrix)

          stan_data <- list(
            N = ncol(Y),
            I = nrow(Y),
            y = Y,
            K = ncol(X),
            X = X,
            M_K = ncol(M[[1]]),
            M = M,
            N_pred = nrow(X_pred),
            X_pred = X_pred,
            M_pred = M_pred
          )
        } else {
          stan_data <- list(
            N = ncol(Y),
            I = nrow(Y),
            y = Y,
            K = ncol(X),
            X = X,
            N_pred = nrow(X_pred),
            X_pred = X_pred
          )
        }
      }

      private$fitted <-
        rstan::sampling(
          private$stan_model,
          data = stan_data,
          ...
        )

      if (length(private$treated_ids) == 1) {
        if (private$predictor_match) {
          private$y_synth_draws <- .get_synth_draws_predictor_match(
            fit = private$fitted,
            pre_data = pre_data,
            post_data = post_data,
            time = private$time,
            outcome = private$outcome
          )
        } else {
          private$y_synth_draws <- .get_synth_draws(
            fit = private$fitted,
            pre_data = pre_data,
            post_data = post_data,
            time = private$time,
            outcome = private$outcome
          )
        }

        private$plot_data <-
          .get_plot_df(
            y_synth_draws = private$y_synth_draws,
            pre_data = pre_data,
            post_data = post_data,
            ci = private$ci_width,
            time = private$time,
            outcome = private$outcome
          )
      } else {
        private$y_synth_draws <- .get_synth_draws3d(
          fit = private$fitted,
          data = private$data,
          id = private$id,
          treated_ids = private$treated_ids,
          time = private$time,
          outcome = private$outcome,
          intervention = private$intervention
        )

        # add plot data
        private$plot_data <- .get_plot_df2(
          y_synth_draws = private$y_synth_draws,
          data = private$data,
          treated_ids = private$treated_ids,
          id = private$id,
          time = private$time,
          outcome = private$outcome,
          ci = private$ci_width
        )
      }
    },

    #' @description
    #' Update the width of the credible interval.
    #' @param ci_width New width for the credible interval. This number should
    #' be in the (0,1) interval.
    updateWidth = function(ci_width = 0.75) {
      stopifnot(exprs = {
        ci_width > 0
        ci_width < 1
      })
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
    #' @description returns descriptive statistics of the lift estimate.
    #'
    summarizeLift = function() {
      if (is.null(private$lift_draws)) {
        stop("You first need to run the `liftDraws()` method.")
      }

      return({
        c(
          point = mean(private$lift_draws$lift),
          lower_bound = stats::quantile(
            private$lift_draws$lift,
            (1 - private$ci_width) / 2
          ),
          upper_bound = stats::quantile(
            private$lift_draws$lift,
            1 - (1 - private$ci_width) / 2
          )
        )
      })
    },
    #' @description effect ggplot2 object that shows the
    #'   effect of the intervention over time.
    #' @param facet Boolean that is TRUE if we want to divide the plot for each
    #' unit.
    #' @param subset Set of units to use in the effect plot.
    effectPlot = function(facet = TRUE, subset = NULL) {
      if (length(private$treated_ids) == 1) {
        tau_plot <- .plot_tau(
          data = private$plot_data,
          x = private$time,
          y = tau,
          ymin = tau_LB,
          ymax = tau_UB,
          xintercept = private$intervention
        )
      } else {
        if (is.null(subset)) {
          tau_plot <- .plot_tau(
            data = private$plot_data,
            x = private$time,
            y = tau,
            ymin = tau_LB,
            ymax = tau_UB,
            xintercept = private$intervention,
            facet = private$id
          )
        } else {
          if (facet) {
            tau_plot <- .plot_tau(
              data = private$plot_data,
              x = private$time,
              y = tau,
              ymin = tau_LB,
              ymax = tau_UB,
              xintercept = private$intervention,
              facet = private$id,
              id = private$id,
              subset = subset
            )
          } else {
            tau_plot <- .plot_tau(
              data = private$plot_data,
              x = private$time,
              y = tau,
              ymin = tau_LB,
              ymax = tau_UB,
              xintercept = private$intervention,
              id = private$id,
              subset = subset
            )
          }
        }
      }

      return(tau_plot)
    },

    #' @description
    #' Plot placebo intervention.
    #' @param periods Positive number of periods for the placebo intervention.
    #' @param ... other arguments passed to [rstan::sampling()].
    #' @return ggplot2 object for placebo treatment effect.
    placeboPlot = function(periods, ...) {
      stopifnot(periods > 0)

      treated <- private$data %>%
        dplyr::filter(status == "Treated") %>%
        dplyr::select(!!private$id) %>%
        dplyr::distinct() %>%
        dplyr::pull(!!private$id)

      keep <- private$data %>%
        dplyr::filter(!!private$time < private$intervention) %>%
        dplyr::select(!!private$time) %>%
        dplyr::distinct() %>%
        dplyr::slice_max(n = periods, order_by = !!private$time) %>%
        dplyr::pull(!!private$time)

      wide_df <- .makeWide(
        data = private$data,
        id = private$id,
        time = private$time,
        outcome = private$outcome,
        treatment = private$treated
      ) %>%
        dplyr::filter(!!private$time < private$intervention) %>%
        dplyr::mutate(
          !!rlang::as_label(private$treated) :=
            dplyr::case_when(!!private$time %in% keep ~ 1, TRUE ~ 0)
        )

      pre_data <- wide_df %>%
        dplyr::filter(!!private$treated == 0)
      post_data <- wide_df %>%
        dplyr::filter(!!private$treated == 1)

      X <- pre_data %>%
        dplyr::select(-!!private$time, -!!private$treated, -!!private$outcome)
      X_pred <- post_data %>%
        dplyr::select(-!!private$time, -!!private$treated, -!!private$outcome)

      if (!is.null(private$covariates)) {
        max_t_pre <- pre_data %>%
          dplyr::pull(!!private$time) %>%
          max()
        max_t_post <- post_data %>%
          dplyr::pull(!!private$time) %>%
          max()

        M <- private$covariates %>%
          dplyr::filter(!!private$time <= max_t_pre)
        M_pred <- private$covariates %>%
          dplyr::filter(
            !!private$time > max_t_pre,
            !!private$time <= max_t_post
          )

        stan_data <- list(
          N = nrow(pre_data),
          y = pre_data %>% dplyr::pull(!!private$outcome),
          K = ncol(X),
          X = X,
          M_K = ncol(M),
          M = M,
          N_pred = nrow(X_pred),
          X_pred = X_pred,
          M_pred = M_pred
        )
      } else {
        stan_data <- list(
          N = nrow(pre_data),
          y = pre_data %>% dplyr::pull(!!private$outcome),
          K = ncol(X),
          X = X,
          N_pred = nrow(X_pred),
          X_pred = X_pred
        )
      }

      placebo_fitted <-
        rstan::sampling(
          private$stan_model,
          data = stan_data,
          ...
        )

      placebo_y_synth_draws <- .get_synth_draws(
        fit = placebo_fitted,
        pre_data = pre_data,
        post_data = post_data,
        time = private$time,
        outcome = private$outcome
      )

      plot_placebo_data <-
        .get_plot_df(
          y_synth_draws = placebo_y_synth_draws,
          pre_data = pre_data,
          post_data = post_data,
          ci = private$ci_width,
          time = private$time,
          outcome = private$outcome
        )

      tau_plot <- .plot_tau(
        data = plot_placebo_data,
        x = private$time,
        y = tau,
        ymin = tau_LB,
        ymax = tau_UB,
        xintercept = min(keep)
      )
      return(tau_plot)
    },
    #' @description
    #' Plots relative upper bias / tau for a time period (firstT, lastT).
    #'
    #' @return vizdraw object with the posterior distribution of relative bias.
    #' Bias is scaled by the time periods.
    #'
    #' @param small_bias Threshold value for considering the bias "small".
    #' @param firstT Start of the time period to compute relative bias over.
    #' Must be after the intervention.
    #' @param lastT End of the time period to compute relative bias over.
    #' Must be after the intervention.
    biasDraws = function(small_bias = 0.3, firstT, lastT) {
      # Calculate the gaps. This won't work with covariates.
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

      # Calculate the estimated Treatment Effects (taus).
      taus_TL <- private$y_synth_draws %>%
        dplyr::filter(!!private$time >= private$intervention) %>%
        dplyr::filter(!!private$time == Ts) %>%
        dplyr::group_by(draw) %>%
        dplyr::mutate(gapTauTL = abs(sum(y_synth - !!private$outcome)))

      # Adjust the relative bias by the number of time periods that passed
      # since the intervention
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
      vizdraws::vizdraws(
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
      ) %>% return()
    },
    #' @description
    #' Plots lift.
    #' @param from First period to consider when calculating lift. If infinite,
    #' set to the time of the intervention.
    #' @param to Last period to consider when calculating lift. If infinite, set
    #' to the last period.
    #' @return vizdraws object with the posterior distribution of the lift.
    #' @param ... other arguments passed to vizdraws::vizdraws().
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
      private$lift_draws <- private$y_synth_draws %>%
        dplyr::filter(!!private$time >= from, !!private$time <= to) %>%
        dplyr::group_by(draw) %>%
        dplyr::summarise(
          sum_y0 = sum(y_synth),
          sum_y1 = sum(!!private$outcome),
          lift = ifelse(abs(sum_y0) > 1e-9, (sum_y1 - sum_y0) / sum_y0, 0.)
        )

      vizdraws::vizdraws(
        posterior = private$lift_draws$lift,
        percentage = TRUE,
        quantity = TRUE,
        tense = "past",
        display_mode_name = TRUE,
        title = "Counterfactual Lift",
        xlab = "Effect of the Intervention",
        units = "the counterfactual lift",
        ...
      ) %>% return()
    },
    # TODO(jvives): Add utility functions to make liftBias and Draws
    # more compact.
    #' @description
    #' Plot Bias magnitude in terms of lift for period (firstT, lastT)
    #' pre_MADs / y0 relative to lift thresholds.
    #' @param firstT start of the time period to compute relative bias
    #'     over. They must be after the intervention.
    #' @param lastT end of the Time period to compute relative bias
    #'     over. They must be after the intervention.
    #' @param offset Target lift %.
    #' @param ... other arguments passed to vizdraws::vizdraws().
    #' @returns vizdraws object with the relative bias with offset.
    liftBias = function(firstT, lastT, offset, ...) {
      # Calculate the gaps. // This wont work with covariates
      mad_pre_gaps <- private$y_synth_draws %>%
        dplyr::filter(!!private$time < private$intervention) %>%
        dplyr::mutate(gap = y_synth - !!private$outcome) %>%
        dplyr::group_by(draw) %>%
        dplyr::summarise(mad = stats::mad(gap))

      # Get all time periods in TL
      Ts <- private$y_synth_draws %>%
        dplyr::filter(!!private$time <= lastT, !!private$time >= firstT) %>%
        dplyr::select(!!private$time) %>%
        unique() %>%
        dplyr::pull()

      # Get time periods between intervention and first_T
      Ts_from_T0 <- private$y_synth_draws %>%
        dplyr::filter(
          !!private$time <= firstT,
          !!private$time >= private$intervention
        ) %>%
        dplyr::select(!!private$time) %>%
        dplyr::n_distinct() %>%
        as.numeric()
      T_from_T0 <- Ts_from_T0 + 1

      # Calculate the estimated Treatment Effects (taus)
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
          dplyr::mutate(bias_ratio = abs(T_from_T0 * 0.5 * TL * (TL - 1) * mad - sumTauTL * offset) / gapTauTL) %>%
          dplyr::filter(bias_ratio < 2.)
      } else {
        rel_bias <- dplyr::inner_join(mad_pre_gaps, taus_TL, by = "draw") %>%
          dplyr::mutate(bias_ratio = abs(T_from_T0 * mad - sumTauTL * offset) / gapTauTL) %>%
          dplyr::filter(bias_ratio < 2.)
      }

      # Return the relative bias vizdraw plot
      vizdraws::vizdraws(
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
      ) %>% return()
    },
    # TODO(jvives): finish implicit weight plots function to work generally
    #' @description
    #' Plot implicit weight distribution across draws.
    #' @return ggplot object with weight distribution per unit.
    weightDraws = function(){
      betas <- private$fitted %>%
        as.data.frame() %>%
        dplyr::select(contains('beta'))

      beta_names <- private$data %>%
        dplyr::filter(!!private$treated == 0) %>%
        dplyr::select(!!private$id) %>%
        unique() %>% pull()

      treated_name <- private$data %>%
        dplyr::filter(!!private$treated == 1) %>%
        dplyr::select(!!private$id) %>%
        unique() %>% pull()

      donor_names <- beta_names[beta_names %in% treated_name == FALSE]

      names(betas) <- donor_names
      melt_betas <- tidyr::gather(betas, ID, weight)

      melt_betas %>% ggplot2::ggplot(ggplot2::aes(x=weight, y=ID, fill=ID)) +
        ggridges::geom_density_ridges() +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = "none") %>% return()
    },
    ## TODO(jvives): finish correlation plot function to work generally
    #' @description
    #' Plots correlations between weights across draws.
    #' @return ggplot heatmap object with correlations.
    weightCorr = function(){
      betas <- private$fitted %>%
        as.data.frame() %>%
        dplyr::select(contains('beta'))

      beta_names <- private$data %>%
        dplyr::filter(!!private$treated == 0) %>%
        dplyr::select(!!private$id) %>%
        unique() %>% pull()

      treated_name <- private$data %>%
        dplyr::filter(!!private$treated == 1) %>%
        dplyr::select(!!private$id) %>%
        unique() %>% pull()

      donor_names <- beta_names[beta_names %in% treated_name == FALSE]

      names(betas) <- donor_names

      cormat <- round(cor(betas),3)
      diag(cormat) <- NA

      melted_cormat <- melt(cormat)
      ggplot2::ggplot(data = melted_cormat,
                      ggplot2::aes(x=X1,
                                   y=X2,
                                   fill=value)) +
        ggplot2::xlab('Units') +
        ggplot2::ylab('Units') +
        ggplot2::labs(fill = 'Corr') +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_gradient2(low='red',
                                      mid = 'white',
                                      high='green') +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
                                                           vjust = 0.5,
                                                           hjust=1)) %>%
        return()
    }
  )
)
