#' Function to plot estimates and confidence intervals across different thresholds
#'
#' @param res Dataframe obtained from running thresholdSurv() function
#' @param monotone TRUE to show the results with the monotonicity assumption
#'
#' @return A ggplot2 graph
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
