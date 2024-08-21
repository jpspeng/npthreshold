#' Non-parametric estimation of the covariate-adjusted threshold-response function
#'
#' https://arxiv.org/pdf/2107.11459
#'
#' @param data Dataframe
#' @param covariates Vector of strings, representing column names for the adjustment covariates
#' @param failure_time A string for the column name of the failure time variable
#' @param event_type A string for the column name fo the event type variable
#' @param marker A string for the column name of the marker variable
#' @param weights A string for the column name of the weights used in estimation
#' @param threshold_list A vector of numeric variables of thresholds to estimate.
#' If NULL, then estimate on the entire dataset
#' @param tf A numeric for reference time point
#' @param nbins_time A numeric for time bins
#' @param nbins_threshold Threshold of time bins ?? (NOT USED HERE)
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
#' @return A dataframe with the estimates, estimates assuming monotonicity,
#' standard errors, confidence intervals, number above the threshold, number of events
#' above the threshold for each threshold in threshold_list.
#'
#' @examples
thresholdSurv <- function(data,
                          covariates,
                          failure_time,
                          event_type,
                          marker,
                          weights,
                          tf,
                          threshold_list = NULL,
                          nbins_time = 20,
                          nbins_threshold = 20,
                          verbose = FALSE,
                          learner.treatment = NULL,
                          learner.event_type = NULL,
                          learner.failure_time = NULL,
                          learner.censoring_time = NULL
) {

  data <- as.data.table(data)
  tf <- as.numeric(tf)

  # subset to relevant variables
  data_select <- data[, c(covariates, failure_time,
                          event_type, marker, weights), with = FALSE]

  tmp_marker <- ifelse(data_select[[weights]] == 0, 0, data_select[[marker]])
  set(data_select, , marker, tmp_marker)

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

    survout$threshold <- threshold
    out_list[[paste0(threshold)]] <- survout
  }

  output <- rbindlist(out_list)

  # make results table
  output$estimate <- unlist(output$estimates )
  output$se <- unlist(output$se )
  output$times <- tf

  n_in_bin <- sapply(threshold_list, function(thresh) {
    sum( data_select[[marker]] >= thresh)
  })

  n_events_in_bin <- sapply(threshold_list, function(thresh) {
    sum(data_select[data_select[[marker]] >= thresh , event_type, with = FALSE] == 1)
  })

  output$n_in_bin <- n_in_bin
  output$n_events_in_bin <- n_events_in_bin

  n_event_cutoff <- 10

  weights_for_iso <- sqrt(n_in_bin)
  weights_for_iso[n_events_in_bin < n_event_cutoff] <- 0.01
  weights_for_iso <- weights_for_iso / sum(weights_for_iso)

  output$estimate_monotone <- -isotone::gpava(output$threshold,
                                               -output$estimate,
                                               weights = weights_for_iso)$x

  output$ci_lo <- unlist(output$CI)[ c(TRUE,FALSE) ]
  output$ci_hi <- unlist(output$CI)[ c(FALSE,TRUE) ]

  output[,c("threshold", "n_in_bin", "n_events_in_bin", "estimate",
            "estimate_monotone", "se", "ci_lo", "ci_hi")]
}

