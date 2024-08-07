#-----------------------------------------------
# obligatory to append to the top of each script
renv::activate(project = here::here(".."))

#-----------------------------------------------
# load parameters

source(here::here("code", "params.R"))


# specify learners
learner.treatment <-  Lrnr_glm$new() #Lrnr_glm$new()
learner.event_type <- Lrnr_hal9001$new(max_degree = 2, num_knots = 20, family = "binomial", smoothness_orders =0,
                                       #lambda = 0, fit_control = list(cv_select = FALSE),
                                       formula = ~   h(risk_score, k = 10) +   h(treatment, pf = 0) + h(failure_time, pf = 0) + h(treatment,failure_time, pf = 0))
stack.failure <-  Lrnr_hal9001$new(max_degree = 2, num_knots = 20, smoothness_orders =0, formula = ~  h(risk_score, k = 10) + h(treatment, pf = 0) + h(t, pf = 0) + h(treatment,t, pf = 0))
stack.censoring <- stack.failure


learner.treatment <-  Lrnr_glm$new()
learner.censoring_time <-  Lrnr_glm$new()
learner.failure_time <-  Lrnr_glm$new()
learner.event_type <-  Lrnr_glm$new()


thresholdSurv <- function(data,
                          covariates,
                          failure_time, # T
                          event_type, #  J (binary) J = 0 (censoring), J =1 , J=2, would be competing events
                          marker,
                          weights,
                          Trt = NULL,
                          threshold_list,
                          t_reference,
                          nbins_time = 20,
                          nbins_threshold = 20,
                          learner.treatment,
                          learner.event_type,
                          learner.failure_time,
                          learner.censoring_time
) {


}

run_competing_risk_analysis = function(data,
                                       covariates,
                                       failure_time, # T
                                       event_type, #  J (binary)
                                       marker,
                                       weights,
                                       Trt,
                                       TwophasesampIndD29,
                                       viral_load,
                                       threshold_list,
                                       nbins_time = 20
){







    ########################
    #### Run competing risk analysis
    ########################

    # Reference time for CR
    # this should be argument
    tf <- tf_by_variant[paste0(variant)]
    # Variant/event type for CR
    event_type_target <- event_type


    # make datasets of placebo and treated analysis
    subset_treated <- data[[Trt]] == 1
    data_placebo <- data[!subset_treated]
    data_treated <- data[subset_treated]
    subset <- which(data_treated[[TwophasesampIndD29]] == 1);  data_treated <- data_treated[subset]
    # subset to relevant variables
    data_treated <- data_treated[, c(covariates, failure_time, event_type_target, marker, weights), with = FALSE]
    data_placebo <-  data_placebo[, c(covariates, failure_time, event_type_target), with = FALSE]
    tmp_marker <- ifelse(data_treated[[weights]] == 0, 0, data_treated[[marker]])
    set(data_treated, , marker, tmp_marker)

    if(any(is.na(data_treated[[marker]]))) {
      print(marker)
      print(table(is.na(data_treated[[marker]])))
      stop("NAs in marker for threshold CR")
    }


    data_placebo[[weights]] <- 1
    print(colnames(data_placebo))
    output_control = run_survtmle3_control(tf = tf,
                          data_surv = data_placebo,
                          covariates = covariates,
                          failure_time = failure_time,
                          event_type = event_type_target,
                          weights = weights,
                          nbins_time = nbins_time)
    #print(c(fit_discrete$est[1,1], est_placebo, unlist(output_control$est)))
    est_placebo <- unlist(output_control$est)
    se_placebo <- unlist(output_control$se)
    print(est_placebo)









    output_treated <- run_survtmle3(tf = tf,
                                    data_surv = data_treated,
                                    covariates = covariates,
                                    failure_time = failure_time,
                                    event_type = event_type_target,
                                    weights = weights,
                                    marker = marker,
                                    nbins_time = nbins_time,
                                    threshold_list = threshold_list)

    # make survival dataset
    #output_treated <- output_treated[-ncol(output_treated)]
    output_treated$estimates <- unlist( output_treated$estimates )
    output_treated$se <- unlist( output_treated$se )
    output_treated$times <- tf



    n_in_bin <- sapply(threshold_list, function(thresh) {
      sum( data_treated[[marker]] >= thresh)
    })

    n_events_in_bin <- sapply(threshold_list, function(thresh) {
      sum(data_treated[data_treated[[marker]] >= thresh , event_type_target, with = FALSE] == 1)
    })

    output_treated$n_in_bin <- n_in_bin
    output_treated$n_events_in_bin <- n_events_in_bin

    n_event_cutoff <- 10


    weights_for_iso <- sqrt(n_in_bin)
    weights_for_iso[n_events_in_bin < n_event_cutoff] <- 0.01
    weights_for_iso <- weights_for_iso / sum(weights_for_iso)


    output_treated$estimates_monotone <- -isotone::gpava(output_treated$threshold, -output_treated$estimates, weights = weights_for_iso)$x

    #output_treated$estimates_monotone <- -as.stepfun(isoreg(output_treated$threshold, -output_treated$estimates))(output_treated$threshold)









    # delta method log(est) ~ sd(IF)/est
    #output_treated$estimates_placebo <- est_placebo
    #output_treated$se_placebo <- se_placebo
    #output_treated$estimates_log_RR <-  log(output_treated$estimates) - log(est_placebo)
    #output_treated$se_log_RR <- sqrt((output_treated$se/output_treated$estimates)^2 + (se_placebo/est_placebo)^2)

    #saveRDS(output_treated, file = here::here(paste0("output/vaccine_", marker, "_", failure_time, "_", event_type_target, ".RDS")))
    #fwrite(data_treated,  here::here(paste0("data_clean/data_treated_", marker, "_", failure_time, "_", event_type_target, ".csv")))
    #fwrite(data_placebo,   here::here(paste0("data_clean/data_placebo_", marker, "_", failure_time, "_", event_type_target, ".csv")))

}






run_survtmle3 <- function(tf, data_survival, covariates, failure_time, event_type, weights , marker = NULL, nbins_time = 30, threshold_list) {

  tf <- as.numeric(tf)

  data_survival <- as.data.table(data_survival)

  # effectivelly removed observations with no weights
  data_survival <- (data_survival[, c(covariates, failure_time, event_type, marker, weights), with = FALSE ])


  # discretize
  time_grid <- unique(quantile(data_survival[[failure_time]], seq(0,1, length = nbins_time+1), type = 1))
  # add time of interest to grid
  time_grid <- sort(union(time_grid, tf))
  failure_time_discrete <- findInterval(data_survival[[failure_time]], time_grid, all.inside = TRUE)
  tf_discrete <- findInterval(tf, time_grid, all.inside = FALSE)
  data_survival[[failure_time]] <-failure_time_discrete



  treatment <- "treatment"
  out_list <- list()
  for(threshold in threshold_list ) {
    print(paste0("THRESHOLD: ", threshold))
    try({
      data_survival[[treatment]] <- 1*(data_survival[[marker]] >= threshold)
      if(sum(data_survival[[treatment]]) == 0){
        print(dim(data_survival))
        print(quantile(data_survival[[marker]]))
        print(marker)
        print(event_type)
        stop("There were zero observation above the threshold. ")
      }
      #print(mean(1*(data_survival[[marker]] >= threshold)))

      survout <- survtmle3_discrete(data_survival[[failure_time]], data_survival[[event_type]],
                                    data_survival[[treatment]], data_survival[, covariates, with = FALSE],
                                    weights = data_survival[[weights]],
                                    learner.treatment =  learner.treatment,
                                    learner.failure_time =  stack.failure,
                                    learner.censoring_time = stack.censoring,
                                    learner.event_type = learner.event_type,
                                    target_failure_time = tf_discrete,
                                    target_treatment = c(1),
                                    target_event_type = 1,
                                    failure_time.stratify_by_time = FALSE,
                                    censoring_time.stratify_by_time = FALSE,
                                    cross_fit = FALSE,
                                    cross_validate = FALSE,
                                    calibrate = FALSE,
                                    verbose = TRUE, max_iter = 100,
                                    tol = 1e-7
      )
    })

    survout$threshold <- threshold
    out_list[[paste0(threshold)]] <- survout
    # })
  }



  output <- rbindlist(out_list)
  return(output)

}


run_survtmle3_control <- function(tf, data_survival, covariates, failure_time, event_type , weights,  nbins_time = 30) {

  tf <- as.numeric(tf)

  data_survival <- as.data.table(data_survival)

  # effectivelly removed observations with no weights
  #data_survival <- (data_survival[, c(covariates, failure_time, event_type, marker, weights), with = FALSE ])


  # discretize
  time_grid <- unique(quantile(data_survival[[failure_time]], seq(0,1, length = nbins_time+1), type = 1))
  # add time of interest to grid
  time_grid <- sort(union(time_grid, tf))
  failure_time_discrete <- findInterval(data_survival[[failure_time]], time_grid, all.inside = TRUE)
  tf_discrete <- findInterval(tf, time_grid, all.inside = FALSE)
  data_survival[[failure_time]] <-failure_time_discrete


  treatment <- "treatment"


  out_list <- list()



    data_survival[[treatment]] <- 1
    #print(mean(1*(data_survival[[marker]] >= threshold)))

    survout <- survtmle3_discrete(data_survival[[failure_time]], data_survival[[event_type]],
                                  data_survival[[treatment]], data_survival[, covariates, with = FALSE],
                                  weights = data_survival[[weights]],
                                  learner.treatment =  learner.treatment,
                                  learner.failure_time =  stack.failure,
                                  learner.censoring_time = stack.censoring,
                                  learner.event_type = learner.event_type,
                                  target_failure_time = tf_discrete,
                                  target_treatment = c(1),
                                  target_event_type = 1,
                                  failure_time.stratify_by_time = FALSE,
                                  censoring_time.stratify_by_time = FALSE,
                                  cross_fit = FALSE,
                                  cross_validate = FALSE,
                                  calibrate = FALSE,
                                  verbose = TRUE, max_iter = 100,
                                  tol = 1e-7
    )


  return(list(est = unlist(survout$estimates, use.names= FALSE), se = unlist(survout$se, use.names = FALSE)))


}




