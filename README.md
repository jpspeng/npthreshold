# Non-parametric estimation of the covariate-adjusted threshold-response function

This (in-development) package implements the non-parametric estimation of the covariate-adjusted threshold-response function. More information about these methods is provided in the paper: https://arxiv.org/pdf/2107.11459.

## Install

```{r}
devtools::install_github("jpspeng/npthreshold")
```
## Using the package 

### With time-to-event outcomes (account for right-censoring)  

```{r}
library(npthreshold)
library(sl3)

# estimates the threshold-response function across specified thresholds
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

# creates of graph of this estimated function with confidence intervals
graphthresh(res)
```
### With binary outcomes 

```{r}
library(npthreshold)
library(SuperLearner)

# estimates the threshold-response function across specified thresholds
res <- thresholdBinary(data = thresh_sample,
                       covariates = c("W1", "W2"), 
                       outcome = "event", 
                       marker = "A", 
                       threshold_list = c(40, 50, 60, 70, 80),
                       Delta = NULL, # no missing outcomes
                       weights = NULL, # all equally weighted
                       sl_library = c("SL.mean", "SL.glm"))

# creates of graph of this estimated function with confidence intervals
graphthresh(res)
```

