#' Plot threshold-response estimates
#'
#' This function creates a plot of estimates and confidence intervals for
#' \eqn{E[E[Y|A \ge v, W]]} across various thresholds.
#'
#' @param res A \code{data.frame} obtained from the \code{thresholdSurv()}
#' function, containing the estimates and related statistics.
#' @param monotone A \code{logical} value. Set to \code{TRUE} to display results
#' under the monotonicity assumption.
#'
#' @return A \code{ggplot2} object displaying the estimates and confidence
#' intervals for \eqn{E[E[Y|A \ge v, W]]} across the specified thresholds.
#' @export
graphthresh <- function(res,
                        type = "raw",
                        ylabel = "Estimate",
                        xlabel = "Thresholds",
                        cutoffs = NA,
                        # exp10 = F,
                        # lod_limit = NA,
                        # lod_label = "LLOQ/2",
                        # include_lod = T,
                        plot_density = F,
                        data = NA,
                        weights = NA,
                        marker = NA,
                        annotate = "",
                        event = NA,
                        time_var = NA,
                        tf = NA){

  if (type == "monotone"){
    y_var <- "estimate_monotone"
    ci_lo_var <- "ci_lo_monotone"
    ci_hi_var <- "ci_hi_monotone"
    scale_coef <- max(res[[ci_hi_var]], na.rm = T)
    yright <- min(res[[y_var]])
  }
  else if (type == "raw"){
    y_var <- "estimate"
    ci_lo_var <- "ci_lo"
    ci_hi_var <- "ci_hi"
    scale_coef <- max(res[[ci_hi_var]], na.rm = T)
    yright <- min(res[[y_var]])
  }
  else if (type == "ve"){
    y_var <- "ve_monotone"
    ci_lo_var <- "ve_monotone_ci_lo"
    ci_hi_var <- "ve_monotone_ci_hi"
    scale_coef <- 1
    yright <- max(res[[y_var]])
  }

  if (!is.na(cutoffs)){
    # implement cutoffs
  }

  ggthresh <- ggplot2::ggplot(res, ggplot2::aes(x = threshold, y = !!sym(y_var))) +
    geom_point(size = 0.2) +
    geom_line() +
    geom_ribbon(aes(ymin = !!sym(ci_lo_var), ymax = !!sym(ci_hi_var)), alpha = 0.3, color = NA) +
    labs(x = "Thresholds", y = "Estimates (CI)") +
    theme_minimal() +
    scale_y_continuous(n.breaks = 10) +
    scale_x_continuous(
      limits = c(min(res$threshold), max(res$threshold))) +
    xlab(xlabel) +
    ylab(ylabel)  +
    theme(plot.title = element_text(hjust = 0.5))

  if (plot_density){

    RCDF <- function(a) {
      sum(data[[weights]] * (data[[marker]] >= a)) /
        sum(data[[weights]]) * scale_coef
    }

    data_event <- data %>%
      filter(.data[[event]] == 1 & .data[[time_var]] <= tf)

    data_event$y_inter <- approx(res$threshold, res[[y_var]],
                                 xout = data_event[,marker],
                                 yright = yright)$y

    RCDF <- Vectorize(RCDF)

    col <- c(col2rgb("olivedrab3")) # orange, darkgoldenrod2
    col <- rgb(col[1], col[2], col[3], alpha = 255 * 0.4, maxColorValue = 255)

    ggthresh <- ggthresh +
      geom_point(aes_string(x = marker, y = "y_inter"), data = data_event, color = "blue") +
      stat_function(fun = RCDF, color = col,
                    geom = "area", fill = col, alpha = 0.2) +
      scale_y_continuous(
        sec.axis = sec_axis(~ . / scale_coef,
                            name = "Reverse CDF"), n.breaks = 10)  +
      theme(plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 0, hjust = 1),
            axis.text.y = element_text(angle = 0, hjust = 1)) +
      geom_vline(xintercept = max(data_event[,marker]),
                 linetype = "dotted",
                 color = "red")

  }

  # if (exp10){
  #   min_value <- min(data[, marker], na.rm = TRUE)
  #
  #
  #   if (min_value < ceiling(min(data[, marker], na.rm = TRUE)) - 0.1){
  #     breaks <- c(min_value,
  #                 ceiling(min_value):floor((max(data[, marker], na.rm = TRUE))))
  #   }
  #   else{
  #     breaks <- c(min_value,
  #                 (ceiling(min_value) + 1):floor((max(data[, marker], na.rm = TRUE))))
  #   }
  #
  #   plt <- plt +
  #     scale_x_continuous(
  #       limits = c(min_value,
  #                  max(data[, marker], na.rm = TRUE)),
  #       breaks = breaks,
  #       labels = label_parsed(custom_labels(breaks, lod_limit, lod_label))
  #     )
  #
  # }
  # else{
  #   min_value <- min(data[, marker], na.rm = TRUE)
  #
  #   breaks <- (ceiling(min_value)):floor((max(data[, marker], na.rm = TRUE)))
  #
  #   plt <- plt +
  #     scale_x_continuous(limits = c(min(data[,marker]), max(data[,marker])),
  #                        breaks = breaks,
  #                        labels = label_parsed(custom_labels(breaks, lod_limit,
  #                                                            lod_label))
  #     )
  # }

  x_annotate_loc <- max(data[,marker]) * 0.85

  ggthresh  +
    annotate("text",
             x = x_annotate_loc,
             y = max(res[,ci_hi_var]) * 0.95,
             label = annotate)

}


# custom_lod_label <- function(breaks, lod_limit, lod_label) {
#   labels <- ifelse(breaks == lod_limit, lod_label, sprintf("10^%s  ", round(breaks,1)))
#   return(labels)
# }
