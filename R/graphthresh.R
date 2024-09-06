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
                        monotone = F,
                        ylabel = "Estimate",
                        xlabel = "Thresholds"){

  if (monotone){
    y_var <- "estimate_monotone"
  }
  else{
    y_var <- "estimate"
  }
  ggplot2::ggplot(res, ggplot2::aes(x = threshold, y = !!sym(y_var))) +
    geom_point() +
    geom_line() +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.3, color = NA) +
    labs(x = "Thresholds", y = "Estimates (CI)") +
    theme_minimal() +
    scale_y_continuous(n.breaks = 10) +
    scale_x_continuous(
      limits = c(min(res$threshold), max(res$threshold))) +
    xlab(xlabel) +
    ylab(ylabel)

}

add_rcdf <- function(ggthresh,
                     data,
                     weights,
                     marker){

  scale_coef <- max(ggthresh$data$ci_hi, na.rm = T)

  RCDF <- function(a) {
    sum(data[[weights]] * (data[[marker]] >= a)) /
      sum(data[[weights]]) * scale_coef
  }

  RCDF <- Vectorize(RCDF)

  col <- c(col2rgb("olivedrab3")) # orange, darkgoldenrod2
  col <- rgb(col[1], col[2], col[3], alpha = 255 * 0.4, maxColorValue = 255)

  ggthresh +
    stat_function(fun = RCDF, color = col,
                  geom = "area", fill = col, alpha = 0.2) +
    scale_y_continuous(
      sec.axis = sec_axis(~ . / scale_coef,
                          name = "Reverse CDF"), n.breaks = 10)  +
    theme(plot.title = element_text(size = 25),
          axis.text.x = element_text(angle = 0, hjust = 1, size = 18),
          axis.text.y = element_text(angle = 0, hjust = 1, size = 18)) +
    scale_x_continuous(limits = c(min(data[[marker]]), max(data[[marker]])))

}
