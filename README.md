# Non-parametric estimation of the covariate-adjusted threshold-response function

This (in-development) package implements the non-parametric estimation of the covariate-adjusted threshold-response function. More information about these methods is provided in the paper: https://arxiv.org/pdf/2107.11459.

## Install

```{r}
devtools::install_github("jpspeng/npthreshold")
```
## Using the package 

```{r}
library(npthreshold)

res <- thresholdSurv(data = thresh_sample,
                     covariates = c("W1", "W2"), 
                     failure_time = "time",
                     event_type = "event", 
                     marker = "A", 
                     weights = NULL, 
                     threshold_list = c(40, 50, 60, 70, 80),
                     tf = 15, 
                     learner.treatment = Lrnr_glm$new(),
                     learner.event_type = Lrnr_glm$new(),
                     learner.failure_time = Lrnr_glm$new(),
                     learner.censoring_time = Lrnr_glm$new(), 
                     verbose = F)

graphthresh(res)
```
