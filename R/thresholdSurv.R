#' Non-parametric estimation of the covariate-adjusted threshold-response function
#'
#' This method estimates \eqn{E[E[Y|A \ge v, W]]} for a range of thresholds \eqn{v}, where \eqn{Y}
#' is a binary outcome of interest, \eqn{A} is a continuous biomarker of interest, and
#' \eqn{W} are baseline variables. To account for missingness in the marker \eqn{A}, the user can
#' provide inverse probability observation weights in the function call.
#'
#' @param data A \code{data.frame} containing the dataset to be analyzed.
#' @param covariates A \code{character} vector specifying the names of columns to be used as adjustment covariates.
#' @param failure_time A \code{character} string indicating the name of the column representing the failure time variable.
#' @param event_type A \code{character} string specifying the name of the column representing the event type variable.
#' @param marker A \code{character} string for the name of the column representing the marker variable.
#' @param weights A \code{character} string specifying the name of the column containing the weights for estimation.
#' @param threshold_list A \code{numeric} vector of threshold values for estimation. If \code{NULL}, estimation is performed on the entire dataset.
#' @param tf A \code{numeric} value indicating the reference time point for the analysis.
#' @param nbins_time A \code{numeric} value specifying the number of time bins to be used.
#' @param nbins_threshold (NOT USED HERE - WHAT IS THIS FOR?)
#' This learner is used as \code{binomial_learner} for \code{\link[sl3]{Lrnr_pooled_hazards}}.
#' Default uses auto ensemble learning by calling `causalutils::get_autoML()`.
#' @param learner.censoring_time A \code{binomial} \code{\link[sl3]{Lrnr_base}} learner object used for fitting conditional hazard model for censoring.
#' This learner is used as \code{binomial_learner} for \code{\link[sl3]{Lrnr_pooled_hazards}}.
#' Default uses auto ensemble learning by calling `causalutils::get_autoML()`.
#' @param learner.event_type A \code{binomial} \code{\link[sl3][Lrnr_base]} learner object used for fitting the conditional probability distribution for failure time of the failure event type.
#' For nonbinary categorical event type, this is used as \code{binomial_learner} in \code{\link[sl3]{Lrnr_independent_binomial}}.
#' Default uses auto ensemble learning by calling `causalutils::get_autoML()`.
#' @param learner.treatment A \code{binomial} \code{\link[sl3]{Lrnr_base}} learner object used for fitting the propensity score model for the treatment mechanism.
#' For nonbinary categorical treatment, this is used as \code{binomial_learner} in \code{\link[sl3]{Lrnr_independent_binomial}}.
#' Default uses auto ensemble learning by calling `causalutils::get_autoML()`.
#' @param verbose A logical value indicating whether to display progress and diagnostic messages during the computation. Default is FALSE.
#'
#' @return A data.frame containing the following columns for each threshold in `threshold_list`:
#' - `estimate`: The estimated value.
#' - `se`: The standard error of the estimate.
#' - `ci_lo`: The lower bound of the confidence interval.
#' - `ci_hi`: The upper bound of the confidence interval.
#' - `estimate_monotone`: The monotone-adjusted estimate.
#' - `ci_lo_monotone`: The lower bound of the monotone-adjusted confidence interval.
#' - `ci_hi_monotone`: The upper bound of the monotone-adjusted confidence interval.
#' - `n_in_bin`: The number of observations above the threshold.
#' - `n_events_in_bin`: The number of events above the threshold.
#'
#' @import sl3 data.table origami
#' @examples
#' @export
thresholdSurv <- function(data,
                          covariates,
                          failure_time,
                          event_type,
                          marker,
                          tf,
                          weights = NULL,
                          threshold_list = NULL,
                          nbins_time = 20,
                          nbins_threshold = 20,
                          verbose = FALSE,
                          learner.treatment = NULL,
                          learner.event_type = NULL,
                          learner.failure_time = NULL,
                          learner.censoring_time = NULL
) {

  data <- data.table::as.data.table(data)
  tf <- as.numeric(tf)

  # set all weights to equal 1 if weights are NULL
  if (is.null(weights)){
    data$weights <- 1
    weights <- "weights"
  }

  # subset to relevant variables
  data_select <- data[, c(covariates, failure_time,
                          event_type, marker, weights), with = FALSE]

  tmp_marker <- ifelse(data_select[[weights]] == 0, 0, data_select[[marker]])
  data.table::set(data_select, , marker, tmp_marker)

  if (any(is.na(data_select[[marker]]))) {
    stop("NAs in marker for threshold CR")
  }

  # discretize
  time_grid <- unique(quantile(data_select[[failure_time]],
                               seq(0, 1,length = nbins_time+1),
                               type = 1))

  # add time of interest to grid
  time_grid <- sort(union(time_grid, tf))
  failure_time_discrete <- findInterval(data_select[[failure_time]],
                                        time_grid, all.inside = TRUE)
  tf_discrete <- findInterval(tf, time_grid, all.inside = FALSE)
  data_select[[failure_time]] <- failure_time_discrete

  # if threshold list is null, then estimate on whole dataset
  if (is.null(threshold_list)){
    threshold_list <- min(data_select[[marker]], na.rm = T)
  }

  out_list <- list()

  for(threshold in threshold_list ) {
    if (verbose) print(paste0("THRESHOLD: ", threshold))

    try({
      data_select[["trt_temp"]] <- 1*(data_select[[marker]] >= threshold)
      if(sum(data_select[["trt_temp"]]) == 0){
        stop(sprintf("There were zero observations above the threshold %s.",
                     threshold))
      }

      survout <- survtmle3_discrete(data_select[[failure_time]],
                                    data_select[[event_type]],
                                    data_select[["trt_temp"]],
                                    data_select[, covariates, with = FALSE],
                                    weights = data_select[[weights]],
                                    learner.treatment =  learner.treatment,
                                    learner.failure_time =  learner.failure_time,
                                    learner.censoring_time = learner.censoring_time,
                                    learner.event_type = learner.event_type,
                                    target_failure_time = tf_discrete,
                                    target_treatment = c(1),
                                    target_event_type = 1,
                                    failure_time.stratify_by_time = FALSE,
                                    censoring_time.stratify_by_time = FALSE,
                                    cross_fit = FALSE,
                                    cross_validate = FALSE,
                                    calibrate = FALSE,
                                    verbose = verbose,
                                    max_iter = 100,
                                    tol = 1e-7
      )
    })

    if (is.null(threshold_list)){
      survout$threshold <- NA
    }
    else{
      survout$threshold <- threshold
    }
    out_list[[paste0(threshold)]] <- survout
  }

  res <- rbindlist(out_list)

  # make results table
  res$estimate <- unlist(res$estimates )
  res$se <- unlist(res$se )
  res$times <- tf

  n_in_bin <- sapply(threshold_list, function(thresh) {
    sum( data_select[[marker]] >= thresh)
  })

  n_events_in_bin <- sapply(threshold_list, function(thresh) {
    sum(data_select[data_select[[marker]] >= thresh , event_type, with = FALSE] == 1)
  })

  res$n_in_bin <- n_in_bin
  res$n_events_in_bin <- n_events_in_bin

  n_event_cutoff <- 10

  weights_for_iso <- sqrt(n_in_bin)
  weights_for_iso[n_events_in_bin < n_event_cutoff] <- 0.01
  weights_for_iso <- weights_for_iso / sum(weights_for_iso)

  res$estimate_monotone <- -isotone::gpava(res$threshold,
                                               -res$estimate,
                                               weights = weights_for_iso)$x

  res$ci_lo <- unlist(res$CI)[ c(TRUE,FALSE) ]
  res$ci_hi <- unlist(res$CI)[ c(FALSE,TRUE) ]

  res$ci_lo_monotone <- res$estimate_monotone - (1.96 * res$se)
  res$ci_hi_monotone <- res$estimate_monotone + (1.96 * res$se)

  res[,c("threshold", "n_in_bin", "n_events_in_bin", "estimate", "se",
         "ci_lo", "ci_hi", "estimate_monotone", "ci_lo_monotone", "ci_hi_monotone")]
}

