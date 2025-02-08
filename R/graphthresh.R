#' Plot threshold-response estimates
#'
#' This function creates a plot of estimates and confidence intervals for
#' \eqn{E[E[Y|A \ge v, W]]} across various thresholds.
#'
#' @param res A \code{data.frame} obtained from the \code{thresholdSurv()}
#' @param type Either "raw", "monotone", or "ve"
#' @param ylabel Label for x-axis of graph
#' @param xlabel Label for y-axis of graph
#' @param plot_density T to include cumulative density function on graph
#' @param plot_endpoints T to include points for disease endpoints on graph
#' @param data Original dataset (needed to graph cdf)
#' @param weights String for weight variable (needed to graph cdf)
#' @param marker String for marker variable (needed to graph cdf)
#' @param annotate String for text to annotate in the corner of the graph
#' @param event String for event variable (needed to graph individual disease
#' endpoints)
#' @param time_var String for time to event variable (needed to graph individual
#' disease endpoints)
#' @param tf Reference timepoint for the analysis
#' @return A \code{ggplot2} object displaying the estimates and confidence
#' intervals for \eqn{E[E[Y|A \ge v, W]]} across the specified thresholds.
#' @import ggplot2
#' @export
graphthresh <- function(res,
                        type = "raw",
                        ylabel = "Estimate",
                        xlabel = "Thresholds",
                        cutoffs = NULL,
                        exp10 = F,
                        # lod_limit = NA,
                        # lod_label = "LLOQ/2",
                        # include_lod = T,
                        plot_density = F,
                        plot_endpoints = F,
                        data = NA,
                        weights = NA,
                        marker = NA,
                        annotate = NA,
                        event = NA,
                        time_var = NA,
                        tf = NA,
                        ylim = NULL){

  res <- data.frame(res)

  if (type == "monotone"){
    y_var <- "estimate_monotone"
    ci_lo_var <- "ci_lo_monotone"
    ci_hi_var <- "ci_hi_monotone"
  }
  else if (type == "raw"){
    y_var <- "estimate"
    ci_lo_var <- "ci_lo"
    ci_hi_var <- "ci_hi"
  }
  else if (type == "ve"){
    y_var <- "ve_monotone"
    ci_lo_var <- "ve_monotone_ci_lo"
    ci_hi_var <- "ve_monotone_ci_hi"
  }

  if (!is.null(cutoffs)){
    res[,ci_lo_var] <- pmax(res[, ci_lo_var], cutoffs[1])
    res[,ci_hi_var] <- pmin(res[, ci_hi_var], cutoffs[2])
  }

  if (type == "ve"){
    scale_coef <- 1
    yright <- max(res[[y_var]])
  }
  else{
    if (is.null(ylim)){
      scale_coef <- max(res[[ci_hi_var]], na.rm = T)
    }
    else{
      scale_coef <- ylim[2]
    }
    yright <- min(res[[y_var]])
  }


  ggthresh <- ggplot(res, aes(x = threshold, y = !!rlang::sym(y_var))) +
    geom_point(size = 0.2) +
    geom_line() +
    geom_ribbon(aes(ymin = !!rlang::sym(ci_lo_var), ymax = !!rlang::sym(ci_hi_var)), alpha = 0.3, color = NA) +
    labs(x = "Thresholds", y = "Estimates (CI)") +
    theme_minimal() +
    scale_y_continuous(n.breaks = 10) +
    xlab(xlabel) +
    ylab(ylabel)  +
    theme(plot.title = element_text(hjust = 0.5))

  if (plot_density){

    RCDF <- function(a) {
      sum(data[[weights]] * (data[[marker]] >= a)) /
        sum(data[[weights]]) * scale_coef
    }

    RCDF <- Vectorize(RCDF)

    col <- c(col2rgb("olivedrab3")) # orange, darkgoldenrod2
    col <- rgb(col[1], col[2], col[3], alpha = 255, maxColorValue = 255)

    ggthresh <- ggthresh +
      stat_function(fun = RCDF, color = col, geom = "line") +
      scale_y_continuous(
        sec.axis = sec_axis(~ . / scale_coef,
                            name = "Reverse CDF"), n.breaks = 10)  +
      theme(plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 0, hjust = 1),
            axis.text.y = element_text(angle = 0, hjust = 1))
  }

  if (plot_endpoints){
    data_event <- data %>%
      filter(.data[[event]] == 1 & .data[[time_var]] <= tf)

    data_event$y_inter <- approx(res$threshold, res[[y_var]],
                                 xout = data_event[,marker],
                                 yright = yright)$y
    ggthresh <- ggthresh +
      geom_point(aes_string(x = marker, y = "y_inter"), data = data_event, color = "blue") +
      geom_vline(xintercept = max(data_event[,marker]),
                 linetype = "dotted",
                 color = "red")

  }

  if (exp10){
    min_value <- min(data[, marker], na.rm = TRUE)
    max_value <- max(data[, marker], na.rm = TRUE)

    if (ceiling(min_value) == floor(max_value) |
        max_value - min_value < 1){
      breaks <- seq(ceiling(min_value * 10) / 10,
                    floor(max_value * 10) / 10, by = 0.1)

      ggthresh <- ggthresh +
        scale_x_continuous(
          limits = c(min_value, max_value),
          breaks = breaks,
          labels = parse(text = paste0("10^", breaks))
        )
    }
    else{
      breaks <- seq(ceiling(min_value), floor(max_value))

      ggthresh <- ggthresh +
        scale_x_continuous(
          limits = c(min_value, max_value),
          breaks = breaks,
          labels = 10^breaks
        )
    }
  }

  if (!is.null(ylim)){
    ggthresh <- ggthresh + ylim(ylim)
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

  if (!is.na(annotate)){
    x_annotate_loc <- max(data[,marker], na.rm = T) * 0.85

    ggthresh <- ggthresh  +
      annotate("text",
               x = x_annotate_loc,
               y = max(res[,ci_hi_var]) * 0.95,
               label = annotate)
  }

  ggthresh
}


# custom_lod_label <- function(breaks, lod_limit, lod_label) {
#   labels <- ifelse(breaks == lod_limit, lod_label, sprintf("10^%s  ", round(breaks,1)))
#   return(labels)
# }
