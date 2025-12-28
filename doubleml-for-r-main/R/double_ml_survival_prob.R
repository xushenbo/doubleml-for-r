#' @title Double machine learning for survival probability
#' @name DoubleMLSurvivalProb
#'
#' @description
#' Double machine learning for survival probability estimation with censoring.
#' This class implements the efficient influence function (UEIF) based
#' estimator for survival probability using censoring martingale theory.
#'
#' @format [R6::R6Class] object inheriting from [DoubleML].
#'
#' @family DoubleML
#' @details
#' The survival probability model estimates the causal effect of treatment on
#' survival probability at specific time points. The model handles right-censored
#' survival data and uses double machine learning to estimate nuisance functions.
#'
#' The key nuisance functions are:
#' - Treatment propensity: P(A=1|X)
#' - Censoring survival functions: G(t|X,A) = P(C > t|X,A)
#' - Event survival functions: S(t|X,A) = P(T > t|X,A)
#'
#' @usage NULL
#'
#' @examples
#' \donttest{
#' # Example with SuperLearner framework
#' library(DoubleML)
#' library(SuperLearner)
#' library(survSuperLearner)
#'
#' # Create survival data
#' set.seed(123)
#' n <- 500
#' x1 <- rnorm(n)
#' x2 <- rnorm(n)
#' a <- rbinom(n, 1, plogis(0.5 * x1 + 0.3 * x2))
#'
#' # Generate survival times
#' lambda <- exp(0.5 * a + 0.2 * x1 + 0.1 * x2)
#' time <- rexp(n, lambda)
#' cens_time <- runif(n, 0, 10)
#' observed_time <- pmin(time, cens_time)
#' event <- as.numeric(time <= cens_time)
#'
#' # Create data frame
#' data <- data.frame(
#'   time = observed_time,
#'   event = event,
#'   a = a,
#'   x1 = x1,
#'   x2 = x2
#' )
#'
#' # Create DoubleMLData
#' obj_dml_data <- DoubleMLSurvivalData$new(
#'   data,
#'   time_col = "time",
#'   event_col = "event",
#'   d_cols = "a",
#'   x_cols = c("x1", "x2")
#' )
#'
#' # Define learners
#' ml_g_surv <- "survSL.km"  # For survival functions
#' ml_g_cens <- "survSL.km"  # For censoring functions
#' ml_m <- "SL.glm"          # For propensity score
#'
#' # Create and fit model
#' dml_surv_prob <- DoubleMLSurvivalProb$new(
#'   obj_dml_data,
#'   ml_g_surv,
#'   ml_g_cens,
#'   ml_m,
#'   admin_cens = 5
#' )
#' dml_surv_prob$fit()
#' dml_surv_prob$summary()
#' }

#' Survival Probability Censoring Martingale Estimation Function
#'
#' Core statistical estimation function for survival probability using efficient
#' influence function (UEIF) based approach with censoring martingale theory.
#' This function computes the doubly robust efficient influence function for
#' survival probability estimation with right-censored data.
#'
#' @param id Vector of individual identifiers
#' @param a Treatment assignment vector (0/1)
#' @param time Observed time (min of event time and censoring time)
#' @param event Event indicator (0=censored, 1+=event)
#' @param bw Balancing weights (inverse propensity weights)
#' @param tilt Tilting weights for target population
#' @param G.a0 Censoring survival function matrix for A=0 (n × n_times)
#' @param G.a1 Censoring survival function matrix for A=1 (n × n_times)
#' @param S.a0 Event survival function matrix for A=0 (n × n_times)
#' @param S.a1 Event survival function matrix for A=1 (n × n_times)
#' @param freq.time Time grid frequency. If NULL, uses unique event times
#' @param admin.cens Administrative censoring time
#'
#' @return List with two components:
#' \describe{
#'   \item{ueif.a1}{Influence function matrix for A=1 (n × n_times)}
#'   \item{ueif.a0}{Influence function matrix for A=0 (n × n_times)}
#' }
#'
#' @details
#' The function implements a 3-term efficient influence function decomposition:
#' \itemize{
#'   \item Term 1: Observed survival indicator weighted by balancing weights
#'   \item Term 2: Predicted survival function weighted by tilting weights
#'   \item Term 3: Censoring martingale correction for doubly robust property
#' }
#'
#' The UEIF ensures the estimator is doubly robust: consistent if either the
#' outcome model (S,G) or the propensity model is correctly specified.
#'
#' @keywords internal
double_ml_survival_prob <- function(id, a, time, event, bw, tilt, G.a0, G.a1, S.a0, S.a1, freq.time=NULL, admin.cens)
{
  n <- length(id)
  causes <- sort(unique(event[event!=0]))
  ncauses <- length(causes)
  if(is.null(freq.time)) {s <- sort(unique(time))} else{
    s <- seq(freq.time, admin.cens, freq.time)}
  ns <- length(s)
  ds <- diff(c(0, s))

  S.a0 <- t(na.locf(t(ifelse(S.a0 < 1e-3, 1e-3, S.a0))))
  S.a1 <- t(na.locf(t(ifelse(S.a1 < 1e-3, 1e-3, S.a1))))
  G.a0 <- t(na.locf(t(ifelse(G.a0 < 1e-3, 1e-3, G.a0))))
  G.a1 <- t(na.locf(t(ifelse(G.a1 < 1e-3, 1e-3, G.a1))))
  G.dHazard.a0 <- t(apply(cbind(0, -log(G.a0)), 1, diff))
  G.dHazard.a1 <- t(apply(cbind(0, -log(G.a1)), 1, diff))

  Yt <- do.call(cbind, lapply(1:ns, function(u)  ifelse(time >= s[u], 1, 0)))
  dNct <- do.call(cbind, lapply(1:ns, function(u)  (dplyr::near(s[u], time))*(event == 0)))

  cens.martingale.integral.a0 <- t(apply((dNct - Yt*G.dHazard.a0)/(S.a0*G.a0), 1, cumsum))
  term1.a0 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*do.call(cbind, lapply(1:ns, function(u)  ifelse(time > s[u], 1, 0)))/G.a0 # cbind(1, G.a0[, 1:(ns-1)])
  term2.a0 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*S.a0
  term3.a0 <- matrix(bw*(a==0),ncol=ns,nrow=n,byrow=FALSE)*S.a0*(cens.martingale.integral.a0-1)
  ueif.a0 <- 1/mean(tilt)*(term1.a0+term2.a0+term3.a0)

  cens.martingale.integral.a1 <- t(apply((dNct - Yt*G.dHazard.a1)/(S.a1*G.a1), 1, cumsum))
  term1.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*do.call(cbind, lapply(1:ns, function(u)  ifelse(time > s[u], 1, 0)))/G.a1 # cbind(1, G.a1[, 1:(ns-1)])
  term2.a1 <- matrix(tilt,ncol=ns,nrow=n,byrow=FALSE)*S.a1
  term3.a1 <- matrix(bw*(a==1),ncol=ns,nrow=n,byrow=FALSE)*S.a1*(cens.martingale.integral.a1-1)
  ueif.a1 <- 1/mean(tilt)*(term1.a1+term2.a1+term3.a1)
  return(list(ueif.a1=ueif.a1, ueif.a0=ueif.a0))
}

#' @export
DoubleMLSurvivalProb = R6Class("DoubleMLSurvivalProb",
  inherit = DoubleMLSurvival, 
  public = list(
    #' @description
    #' Creates a new instance of this R6 class.
    #'
    #' @param data (`DoubleMLSurvivalData`) \cr
    #' The `DoubleMLSurvivalData` object providing the survival data and 
    #' specifying the variables of the survival model.
    #'
    #' @param ml_g_surv (character vector) \cr
    #' survSuperLearner library for event survival functions S(t|X,A).
    #' Use survSuperLearner library names like c("survSL.km", "survSL.coxph").
    #'
    #' @param ml_g_cens (character vector) \cr
    #' survSuperLearner library for censoring survival functions G(t|X,A).
    #' Use survSuperLearner library names like c("survSL.km", "survSL.coxph").
    #'
    #' @param ml_m (character vector) \cr
    #' SuperLearner library for treatment propensity P(A=1|X).
    #' Use SuperLearner library names like c("SL.glm", "SL.ranger").
    #'
    #' @param admin_cens (`numeric(1)`) \cr
    #' Administrative censoring time (restriction time for analysis).
    #'
    #' @param time_col (`character(1)`) \cr
    #' Name of the time column in the data.
    #'
    #' @param event_col (`character(1)`) \cr
    #' Name of the event column in the data (0=censored, 1+=event).
    #'
    #' @param tau (`numeric(1)` or `numeric(vector)`) \cr
    #' Evaluation time point(s) for survival probability estimation. Default is admin_cens.
    #'
    #' @param freq_time (`numeric(1)`, optional) \cr
    #' Time grid frequency. If NULL, uses unique event times.
    #'
    #' @param ml_framework (`character(1)`) \cr
    #' Machine learning framework to use. Only "SuperLearner" is supported.
    #' Default is "SuperLearner".
    #'
    #' @param n_folds (`integer(1)`) \cr
    #' Number of folds for cross-fitting. Default is 5.
    #'
    #' @param n_rep (`integer(1)`) \cr
    #' Number of repetitions for sample splitting. Default is 1.
    #'
    #' @param score (`character(1)`) \cr
    #' Score function specification. Default is "ueif".
    #'
    #' @param dml_procedure (`character(1)`) \cr
    #' DML procedure. Default is "dml2".
    #'
    #' @param draw_sample_splitting (`logical(1)`) \cr
    #' Whether to draw sample splitting. Default is TRUE.
    #'
    #' @param apply_cross_fitting (`logical(1)`) \cr
    #' Whether to apply cross-fitting. Default is TRUE.
    initialize = function(data,
                          ml_g_surv,
                          ml_g_cens,
                          ml_m,
                          admin_cens,
                          time_col,
                          event_col,
                          tau = NULL,
                          freq_time = NULL,
                          ml_framework = "SuperLearner",
                          n_folds = 5,
                          n_rep = 1,
                          score = "ueif",
                          dml_procedure = "dml2",
                          draw_sample_splitting = TRUE,
                          apply_cross_fitting = TRUE) {

      # Store survival-specific parameters
      private$admin_cens_ <- admin_cens
      private$time_col_ <- time_col
      private$event_col_ <- event_col
      private$freq_time_ <- freq_time
      private$ml_framework_ <- ml_framework
      private$tau_ <- if (is.null(tau)) admin_cens else tau

      # Validate tau values
      private$validate_tau()

      # Store ML learners
      private$ml_g_surv_ <- ml_g_surv
      private$ml_g_cens_ <- ml_g_cens
      private$ml_m_ <- ml_m

      # Check data structure
      private$check_data(data)

      # Check score function
      private$check_score(score)

      # Validate ML framework
      if (ml_framework != "SuperLearner") {
        stop("Only SuperLearner framework is supported")
      }
      
      # Initialize parent DoubleML class
      super$initialize_double_ml(
        data,
        n_folds,
        n_rep,
        score,
        dml_procedure,
        draw_sample_splitting,
        apply_cross_fitting
      )
    },
    
    #' @description
    #' Get full survival probability inference curve with confidence bands.
    #'
    #' @param level (`numeric(1)`) \cr
    #' Confidence level. Default is 0.95.
    #'
    #' @param n_bootstrap (`integer(1)`) \cr
    #' Number of bootstrap samples for confidence bands. Default is 500.
    #'
    #' @return Data frame with time points, estimates, and confidence intervals.
    get_inference_curve = function(level = 0.95, n_bootstrap = 500) {
      if (is.null(private$all_ueif_)) {
        stop("Model must be fitted first. Call fit() method.")
      }

      # Use survival_prob_inference to get full curve with confidence bands
      survival_prob_inference(
        ueif.a1 = private$all_ueif_$ueif.a1,
        ueif.a0 = private$all_ueif_$ueif.a0,
        time = self$get_time_grid(),
        level = level,
        n_bootstrap = n_bootstrap
      )
    }
  ),

  active = list(
    #' @field coef Treatment effect estimate
    coef = function() private$coef_,

    #' @field se Standard error of treatment effect estimate
    se = function() private$se_,

    #' @field tau Evaluation time point(s)
    tau = function(value) {
      if (missing(value)) {
        return(private$tau_)
      } else {
        private$tau_ <- value
        private$validate_tau()
      }
    },

    #' @field survival_prob_estimates All survival probability estimates and standard errors
    survival_prob_estimates = function() private$survival_prob_estimates_,

    #' @field admin_cens Administrative censoring time
    admin_cens = function() private$admin_cens_
  ),
  
  private = list(
    admin_cens_ = NULL,
    freq_time_ = NULL,
    ml_framework_ = NULL,
    time_col_ = NULL,
    event_col_ = NULL,
    tau_ = NULL,
    coef_ = NULL,
    se_ = NULL,
    survival_prob_estimates_ = NULL,
    all_ueif_ = NULL,
    ml_g_surv_ = NULL,
    ml_g_cens_ = NULL,
    ml_m_ = NULL,
    
    check_data = function(obj_dml_data) {
      # Check that required columns exist in the data
      if (!private$time_col_ %in% names(obj_dml_data$data)) {
        stop(paste("Time column", private$time_col_, "not found in data"))
      }
      if (!private$event_col_ %in% names(obj_dml_data$data)) {
        stop(paste("Event column", private$event_col_, "not found in data"))
      }
      invisible(NULL)
    },
    
    check_score = function(score) {
      valid_scores = "ueif"
      if (!score %in% valid_scores) {
        stop(paste("Invalid score. Valid scores:", paste(valid_scores, collapse = ", ")))
      }
      return(invisible(NULL))
    },

    validate_tau = function() {
      if (any(private$tau_ <= 0)) {
        stop("tau must be positive")
      }
      if (any(private$tau_ > private$admin_cens_)) {
        stop("tau must be <= admin_cens")
      }
      invisible(NULL)
    },

    # Score elements extraction
    score_elements = function(ueif_a1, ueif_a0, tau, time_grid) {
      # For single tau: extract influence function at specific tau
      if (length(tau) == 1) {
        # Use survival_prob_extract_average to get point estimates at tau
        results <- survival_prob_extract_average(
          train.ueif.a1 = ueif_a1,
          train.ueif.a0 = ueif_a0,
          estimation.ueif.a1 = ueif_a1,
          estimation.ueif.a0 = ueif_a0,
          train.time = time_grid,
          estimation.time = time_grid,
          tau = tau
        )

        # Store detailed estimates
        private$survival_prob_estimates_ <- list(
          tau = tau,
          estimate.a0 = results$estimate.a0,
          estimate.a1 = results$estimate.a1,
          estimate.diff = results$estimate.diff,
          se.a0 = results$estimate.a0.se,
          se.a1 = results$estimate.a1.se,
          se.diff = results$estimate.diff.se
        )

        # Extract individual-level influence functions at tau
        # Find time index closest to tau
        tau_idx <- which.min(abs(time_grid - tau))

        # Individual influence functions at tau
        ueif_a1_at_tau <- ueif_a1[, tau_idx]
        ueif_a0_at_tau <- ueif_a0[, tau_idx]

        # Score for treatment effect: ψ(W) = ψ_1(W) - ψ_0(W)
        psi_b <- ueif_a1_at_tau - ueif_a0_at_tau

        # Jacobian component (for variance): ∂θ/∂ψ = -1 for linear functional
        psi_a <- rep(-1, nrow(ueif_a1))

      } else {
        # Handle multiple tau values
        estimates <- lapply(tau, function(tau_i) {
          survival_prob_extract_average(
            train.ueif.a1 = ueif_a1,
            train.ueif.a0 = ueif_a0,
            estimation.ueif.a1 = ueif_a1,
            estimation.ueif.a0 = ueif_a0,
            train.time = time_grid,
            estimation.time = time_grid,
            tau = tau_i
          )
        })

        # Store detailed estimates for all tau
        private$survival_prob_estimates_ <- list(
          tau = tau,
          estimate.a0 = sapply(estimates, function(x) x$estimate.a0),
          estimate.a1 = sapply(estimates, function(x) x$estimate.a1),
          estimate.diff = sapply(estimates, function(x) x$estimate.diff),
          se.a0 = sapply(estimates, function(x) x$estimate.a0.se),
          se.a1 = sapply(estimates, function(x) x$estimate.a1.se),
          se.diff = sapply(estimates, function(x) x$estimate.diff.se)
        )

        # Extract individual-level influence functions at first tau
        # (base class expects single vectors for psi_a/psi_b)
        tau_idx <- which.min(abs(time_grid - tau[1]))
        ueif_a1_at_tau <- ueif_a1[, tau_idx]
        ueif_a0_at_tau <- ueif_a0[, tau_idx]

        # Use first tau for psi_a, psi_b (base class expects vectors)
        # For multiple tau, user should access survival_prob_estimates directly
        psi_b <- ueif_a1_at_tau - ueif_a0_at_tau
        psi_a <- rep(-1, nrow(ueif_a1))
      }

      return(list(
        psi_a = psi_a,
        psi_b = psi_b
      ))
    },

    nuisance_est = function(smpls, ...) {
      # Extract full data
      full_data <- self$data$data
      n_obs <- self$data$n_obs

      # Get time grid
      if (is.null(private$freq_time_)) {
        time_grid <- sort(unique(full_data[[private$time_col_]]))
      } else {
        time_grid <- seq(private$freq_time_, private$admin_cens_, private$freq_time_)
      }
      n_times <- length(time_grid)

      # Initialize prediction matrices for all observations
      S_a0_pred <- matrix(NA, nrow = n_obs, ncol = n_times)
      S_a1_pred <- matrix(NA, nrow = n_obs, ncol = n_times)
      G_a0_pred <- matrix(NA, nrow = n_obs, ncol = n_times)
      G_a1_pred <- matrix(NA, nrow = n_obs, ncol = n_times)
      propensity_pred <- rep(NA, n_obs)

      # Loop over folds for cross-fitting
      for (i_fold in seq_along(smpls$train_ids)) {
        train_ids <- smpls$train_ids[[i_fold]]
        test_ids <- smpls$test_ids[[i_fold]]

        # Get train/test data
        train_data <- full_data[train_ids, ]
        test_data <- full_data[test_ids, ]

        # Fit nuisance functions on train, predict on test (pass consistent time_grid)
        ml_results <- private$ml_nuisance_est_superlearner(train_data, test_data, time_grid)

        # Store predictions for test fold
        S_a0_pred[test_ids, ] <- ml_results$S_a0
        S_a1_pred[test_ids, ] <- ml_results$S_a1
        G_a0_pred[test_ids, ] <- ml_results$G_a0
        G_a1_pred[test_ids, ] <- ml_results$G_a1
        propensity_pred[test_ids] <- ml_results$propensity
      }

      # Calculate balancing weights (inverse propensity weights)
      bw <- ifelse(full_data[[self$data$d_cols]] == 1,
                   1 / propensity_pred,
                   1 / (1 - propensity_pred))

      # Calculate tilting weights (for ATE, uniform weights)
      tilt <- rep(1, n_obs)

      # Calculate UEIF using double_ml_survival_prob function
      ueif_results <- double_ml_survival_prob(
        id = 1:n_obs,
        a = full_data[[self$data$d_cols]],
        time = full_data[[private$time_col_]],
        event = full_data[[private$event_col_]],
        bw = bw,
        tilt = tilt,
        G.a0 = G_a0_pred,
        G.a1 = G_a1_pred,
        S.a0 = S_a0_pred,
        S.a1 = S_a1_pred,
        freq.time = private$freq_time_,
        admin.cens = private$admin_cens_
      )

      # Store UEIF for inference curves
      private$all_ueif_ <- ueif_results

      # Extract score components using score_elements
      score_comp <- private$score_elements(
        ueif_a1 = ueif_results$ueif.a1,
        ueif_a0 = ueif_results$ueif.a0,
        tau = private$tau_,
        time_grid = time_grid
      )

      # Return required format for base class
      return(list(
        psi_a = score_comp$psi_a,
        psi_b = score_comp$psi_b,
        preds = list(
          S_a0 = S_a0_pred,
          S_a1 = S_a1_pred,
          G_a0 = G_a0_pred,
          G_a1 = G_a1_pred,
          propensity = propensity_pred
        ),
        models = NULL  # Can store fitted models if needed
      ))
    },

    ml_nuisance_est_superlearner = function(train_data, test_data, time_grid) {

      # Estimate propensity scores P(A=1|X) using SuperLearner
      propensity_raw <- fit_propensity_score(
        data = train_data,
        treatment_col = self$data$d_cols,
        covariate_cols = self$data$x_cols,
        sl_library = private$ml_m_,
        new_data = test_data
      )

      # Apply winsorization to propensity scores
      propensity <- winsorize_values(propensity_raw, lower = 1e-3,
                                    upper = 1 - 1e-3)

      # Estimate both survival S(t|X,A) and censoring G(t|X,A) functions using survSuperLearner
      # Single call returns both event.SL.predict and cens.SL.predict
      surv_results <- fit_survival_functions(
        data = train_data,
        time_col = private$time_col_,
        event_col = private$event_col_,
        treatment_col = self$data$d_cols,
        covariate_cols = self$data$x_cols,
        event_sl_lib = private$ml_g_surv_,
        cens_sl_lib = private$ml_g_cens_,
        new_data = test_data,
        new_times = time_grid
      )

      # Apply winsorization to survival functions
      s_a0 <- winsorize_values(surv_results$S_a0, lower = 1e-3,
                              upper = 1 - 1e-3)
      s_a1 <- winsorize_values(surv_results$S_a1, lower = 1e-3,
                              upper = 1 - 1e-3)

      # Apply winsorization to censoring functions
      g_a0 <- winsorize_values(surv_results$G_a0, lower = 1e-3,
                              upper = 1 - 1e-3)
      g_a1 <- winsorize_values(surv_results$G_a1, lower = 1e-3,
                              upper = 1 - 1e-3)

      list(
        G_a0 = g_a0,
        G_a1 = g_a1,
        S_a0 = s_a0,
        S_a1 = s_a1,
        propensity = propensity
      )
    }
  )
)
