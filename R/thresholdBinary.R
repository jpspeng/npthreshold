#' Non-parametric estimation of the covariate-adjusted threshold-response function
#'
#' This function estimates the covariate-adjusted probability of a binary outcome above
#' a range of thresholds for a continuous marker.
#'
#' @param data A \code{data.frame} containing the dataset to be analyzed.
#' @param covariates A \code{character} vector specifying the names of columns to be used as adjustment covariates.
#' @param outcome A \code{character} string specifying the name of the binary outcome column.
#' @param marker A \code{character} string specifying the name of the column representing the continuous marker.
#' @param Delta A \code{character} string specifying the name of the column indicating whether the outcome is observed.
#' If \code{NULL}, all outcome data is assumed to be available.
#' @param weights A \code{character} string specifying the name of the column containing the observation weights.
#' If \code{NULL}, all data are weighted equally.
#' @param threshold_list A \code{numeric} vector of threshold values for estimation.
#' If \code{NULL}, estimation is performed using the minimum value of the marker.
#' @param sl_library A \code{character} vector specifying the Super Learner librarIes to be used for estimation.
#' The default is \code{c("SL.mean", "SL.glm")}.
#' @param method A \code{character} string specifying the method used by Super Learner. The default is \code{"method.CC_nloglik"}.
#' @param cvControl A list specifying the cross-validation control parameters for Super Learner. The default is \code{list(V = 5)}.
#'
#' @return A \code{data.frame} containing the following columns for each threshold in \code{threshold_list}:
#' \itemize{
#'   \item \code{threshold}: The threshold value.
#'   \item \code{estimate}: The estimated value at the threshold.
#'   \item \code{se}: The standard error of the estimate.
#'   \item \code{ci_lo}: The lower bound of the confidence interval.
#'   \item \code{ci_hi}: The upper bound of the confidence interval.
#'   \item \code{estimate_monotone}: The monotone-adjusted estimate.
#'   \item \code{ci_lo_monotone}: The lower bound of the monotone-adjusted confidence interval.
#'   \item \code{ci_hi_monotone}: The upper bound of the monotone-adjusted confidence interval.
#'   \item \code{n_in_bin}: The number of observations above the threshold.
#'   \item \code{n_events_in_bin}: The number of events above the threshold.
#' }
#'
#' @import data.table isotone SuperLearner
#' @export
thresholdBinary <- function(data,
                            covariates,
                            outcome,
                            marker,
                            Delta = NULL,
                            weights = NULL,
                            threshold_list = NULL,
                            sl_library = c("SL.mean", "SL.glm"),
                            method = "method.CC_nloglik",
                            cvControl = list(V = 5)){

  data <- data.frame(data)

  if (is.null(Delta)){
    data$Delta <- 1
    Delta <- "Delta"
  }

  if (is.null(weights)){
    data$weights <- 1
    weights <- "weights"
  }

  data$marker <- ifelse(data[,weights] == 0, 0, data[,marker])

  if (any(is.na(data[[marker]]))) {
    stop("NAs in marker for threshold CR")
  }

  # if threshold list is null, then estimate on whole dataset
  if (is.null(threshold_list)){
    threshold_list <- min(data_select[[marker]], na.rm = T)
  }

  results <- data.frame()

  n <- nrow(data)

  for (threshold in threshold_list) {

    covariates_matrix <- model.matrix(~ . - 1, data = data[, covariates])

    res_temp <- tmleThreshold.auto(threshold = threshold,
                                   W = covariates_matrix,
                                   A = data[, marker],
                                   Y = data[, outcome],
                                   Delta = data[, Delta],
                                   weights = data[, weights],
                                   sl_library = sl_library)

    results <- rbind(results, res_temp)

  }

  n_in_bin <- sapply(threshold_list, function(thresh) {
    sum( data[[marker]] >= thresh)
  })

  n_events_in_bin <- sapply(threshold_list, function(thresh) {
    sum(data[data[[marker]] >= thresh, outcome] == 1)
  })

  results$n_in_bin <- n_in_bin
  results$n_events_in_bin <- n_events_in_bin

  n_event_cutoff <- 10

  weights_for_iso <- sqrt(n_in_bin)
  weights_for_iso[n_events_in_bin < n_event_cutoff] <- 0.01
  weights_for_iso <- weights_for_iso / sum(weights_for_iso)

  results$estimate_monotone <- -isotone::gpava(results$threshold,
                                               -results$estimate,
                                               weights = weights_for_iso)$x

  results$ci_lo_monotone <- results$estimate_monotone - (1.96 * results$se)
  results$ci_hi_monotone <- results$estimate_monotone + (1.96 * results$se)

  results

}


tmleThreshold.auto <- function(threshold, W, A, Y, Delta = rep(1, length(A)),
                               weights = rep(1, length(A)),
                               sl_library = c("SL.mean", "SL.glm"),
                               method = "method.CC_nloglik",
                               cvControl = list(V = 5),
                               run_mono = FALSE) {
  if (is.null(sl_library)) {
    # Define the SuperLearner library
    sl_library <- c("SL.mean", "SL.glm", "SL.gam", "SL.ranger", "SL.earth")
  }

  n <- length(A)
  W <- as.matrix(W)
  colnames(W) <- paste0("W", 1:ncol(W))

  Av <- as.numeric(A >= threshold)

  # Task setup for SuperLearner with cvControl
  lrnr_Qv <- SuperLearner::SuperLearner(Y = Y[Delta == 1 & Av == 1], X = as.data.frame(W[Delta == 1 & Av == 1, ]),
                          newX = as.data.frame(W), family = binomial(), SL.library = sl_library,
                          method = method, cvControl = cvControl)

  if (all(Av == 1)) {
    lrnr_gv <- SuperLearner::SuperLearner(Y = Av, X = as.data.frame(W), family = binomial(), SL.library = "SL.mean",
                            method = method, cvControl = cvControl[1])
  } else {
    lrnr_gv <- SuperLearner::SuperLearner(Y = Av, X = as.data.frame(W), family = binomial(), SL.library = sl_library,
                            method = method, cvControl = cvControl[1])
  }

  if (all(Delta == 1)) {
    lrnr_Gv <- SuperLearner::SuperLearner(Y = Delta, X = as.data.frame(cbind(W, Av)), family = binomial(), SL.library = "SL.mean",
                            method = method, cvControl = cvControl[1])
  } else {
    lrnr_Gv <- SuperLearner::SuperLearner(Y = Delta, X = as.data.frame(cbind(W, Av)), family = binomial(), SL.library = sl_library,
                            method = method, cvControl = cvControl[1])
  }

  Qv <- lrnr_Qv$SL.predict
  gv <- lrnr_gv$SL.predict
  Gv <- lrnr_Gv$SL.predict

  out_v <- suppressWarnings(tmle.inefficient(threshold, A, Delta, Y,
                                             gv, Qv, Gv,
                                             weights))
  psi <- out_v$psi
  IF <- out_v$EIF
  se <- sd(IF) / sqrt(n)
  CI_left <- psi - 1.96 * se
  CI_right <- psi + 1.96 * se

  list(threshold = threshold,
       estimate = psi,
       se = se,
       ci_lo = CI_left,
       ci_hi = CI_right)
}


#The binary-treatment-based TMLE (binTMLE):
tmle.inefficient <- function(threshold, A, Delta, Y,
                             gv, Qv, Gv,
                             weights = rep(1, length(A))) {
  n <- length(A)
  gv <- truncate_pscore_adaptive(A, gv)
  Gv <- truncate_pscore_adaptive(Delta, Gv)
  Av <- as.numeric(A >= threshold)
  Qv <- pmax(Qv, 1e-8)
  Qv <- pmin(Qv, 1-1e-8)
  eps <- coef(glm.fit(as.matrix(rep(1,n)), Y,
                      offset = qlogis(Qv), weights = weights * Delta*Av/(gv*Gv),
                      family = binomial(), start = 0))
  Qv_star <- plogis(qlogis(Qv) + eps )
  psi <- weighted.mean(Qv_star, weights)
  EIF <- Delta*Av/(gv*Gv) * (Y - Qv_star) + Qv_star - psi
  EIF <- weights * EIF
  return(list(psi = psi, EIF = EIF))
}


truncate_pscore_adaptive <- function(A, pi, min_trunc_level = 1e-8) {
  risk_function <- function(cutoff, level) {
    pi <- pmax(pi, cutoff)

    alpha <- A/pi #Riesz-representor
    alpha1 <- 1/pi
    mean(alpha^2 - 2*(alpha1))
  }
  cutoff <- optim(1e-3, fn = risk_function, method = "Brent", lower = min_trunc_level, upper = 0.5, level = 1)$par

  pi <- pmax(pi, cutoff)
  return(pi)
}

