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