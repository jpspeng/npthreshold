#' Plot threshold-response estimates
#'
#' This function creates a plot of estimates and confidence intervals for \eqn{E[E[Y|A \ge v, W]]} across various thresholds.
#'
#' @param res A \code{data.frame} obtained from the \code{thresholdSurv()} function, containing the estimates and related statistics.
#' @param monotone A \code{logical} value. Set to \code{TRUE} to display results under the monotonicity assumption.
#'
#' @return A \code{ggplot2} object displaying the estimates and confidence intervals for \eqn{E[E[Y|A \ge v, W]]} across the specified thresholds.
#' @export
graphthresh <- function(res,
                        monotone = F){

  if (monotone){
    y_var <- "estimate_monotone"
  }
  else{
    y_var <- "estimate"
  }
  ggplot2::ggplot(res, aes(x = factor(threshold), y = !!sym(y_var))) +
    geom_point() +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.2) +
    labs(x = "Thresholds", y = "Estimates (CI)") +
    theme_minimal()
}
