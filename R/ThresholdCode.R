
sl3_Learner <- Lrnr_gam$new()




tmleThreshold.auto <- function(thresholds, W , A , Y,  Delta = rep(1, length(A)), weights = rep(1, length(A)), sl3_Learner = NULL) {
  require(data.table)
  require(sl3)
  require(earth)
  require(mgcv)
  require(ranger)
  lrnr_glm <-  Lrnr_glmnet$new()
  if(ncol(as.matrix(W) == 1)) {
    lrnr_glm <-  Lrnr_glm$new()
  }
  if(is.null(sl3_Learner)) {
    lrnr_stack <- Stack$new(list(
      Lrnr_mean$new(),
      Lrnr_gam$new(),
      Lrnr_earth$new(degree = 1),
      Lrnr_earth$new(degree = 2),
      lrnr_glm,
      Lrnr_ranger$new()
    ))
    sl3_Learner <- make_learner(Pipeline, Lrnr_cv$new(lrnr_stack), Lrnr_cv_selector$new(loss_squared_error))
  }
  n <- length(A)
  W <- as.matrix(W)
  colnames(W) <- paste0("W", 1:ncol(W))
  results <- data.table()

  for(threshold in thresholds) {
    Av <- as.numeric(A >= threshold)
    data <- data.table(W, A = A, Av = Av, Delta = Delta, Y = Y, weights = weights)

    task_Qv <- sl3_Task$new(data, covariates = colnames(W), outcome = "Y", weights = "weights")
    task_gv <- sl3_Task$new(data, covariates = colnames(W), outcome = "Av", weights = "weights")
    task_Gv <-  sl3_Task$new(data, covariates =  c(colnames(W), "Av"), outcome = "Delta", weights = "weights")
    lrnr_Qv <- sl3_Learner$train(task_Qv[Delta == 1 & Av ==1])
    lrnr_gv <- sl3_Learner$train(task_gv)
    if(all(Delta==1)) {
      lrnr_Gv <- Lrnr_mean$new()$train(task_Gv)
    } else {
      lrnr_Gv <- sl3_Learner$train(task_Gv)
    }
    Qv <- lrnr_Qv$predict(task_Qv)
    gv <- lrnr_gv$predict(task_gv)
    Gv <- lrnr_Gv$predict(task_Gv)

    out_v <- suppressWarnings(tmle.inefficient(threshold, A, Delta, Y,
                                               gv, Qv, Gv,
                                               weights))
    psi <- out_v$psi
    IF <- out_v$EIF
    se <- sd(IF)/sqrt(n)
    CI_left <- psi - 1.96 * se
    CI_right <- psi + 1.96 * se
    results <-rbind(results,
                    data.table(threshold = threshold, psi = psi,
                               se = se,  CI_left = CI_left, CI_right = CI_right))
  }


  psi <- results$psi
  psi_mono <- as.stepfun(isoreg(thresholds, psi))(results$threshold)
  results$psi_mono <- psi_mono
  results$CI_left_mono <- psi_mono - 1.96 * results$se
  results$CI_right_mono <- psi_mono + 1.96 * results$se

  return(results)



}


# Sequential regression TMLE
tmle.efficient <- function(threshold, A, Delta, Y,
                           gv, Q, Qv, G,
                           weights = rep(1, length(A)), bound = 0.005) {
  n <- length(A)
  gv <- pmax(gv, bound)
  G <- pmax(G, bound)
  Av <- as.numeric(A >= threshold)
  # Step 1
  Q <- pmax(Q, 1e-8)
  Q <- pmin(Q, 1-1e-8)
  eps <- coef(glm.fit(Av, Y, offset = qlogis(Q),
                      weights = weights * Delta/(gv*G),
                      family = binomial(), start = 0))
  Q_star <- plogis(qlogis(Q) + eps * Av )
  # Step 2 (sequential regression)
  Qv <- pmax(Qv, 1e-8)
  Qv <- pmin(Qv, 1-1e-8)
  eps <- coef(glm.fit(as.matrix(rep(1,n)), Q_star,
                      offset = qlogis(Qv), weights = weights * Av/(gv),
                      family = binomial(), start = 0))
  Qv_star  <- plogis(qlogis(Qv) + eps )

  psi <- weighted.mean(Qv_star, weights)
  EIF <- Delta*Av/(gv*G) * (Y - Q_star) + Av/gv * (Q_star - Qv_star) + Qv_star - psi
  EIF <- weights * EIF
  return(list(psi = psi, EIF = EIF))
}


#The binary-treatment-based TMLE (binTMLE):

tmle.inefficient <- function(threshold, A, Delta, Y,
                             gv, Qv, Gv,
                             weights = rep(1, length(A))) {
  n <- length(A)
  gv <- truncate_pscore_adaptive(A, gv)
  Gv <- truncate_pscore_adaptive(Delta, Gv)
  Av <- as.numeric(A >= threshold)
  Qv <- pmax(Qv, 1e-8)
  Qv <- pmin(Qv, 1-1e-8)
  eps <- coef(glm.fit(as.matrix(rep(1,n)), Y,
                      offset = qlogis(Qv), weights = weights * Delta*Av/(gv*Gv),
                      family = binomial(), start = 0))
  Qv_star <- plogis(qlogis(Qv) + eps )
  psi <- weighted.mean(Qv_star, weights)
  EIF <- Delta*Av/(gv*Gv) * (Y - Qv_star) + Qv_star - psi
  EIF <- weights * EIF
  return(list(psi = psi, EIF = EIF))
}



truncate_pscore_adaptive <- function(A, pi, min_trunc_level = 1e-8) {
  risk_function <- function(cutoff, level) {
    pi <- pmax(pi, cutoff)

    alpha <- A/pi #Riesz-representor
    alpha1 <- 1/pi
    mean(alpha^2 - 2*(alpha1))
  }
  cutoff <- optim(1e-3, fn = risk_function, method = "Brent", lower = min_trunc_level, upper = 0.5, level = 1)$par

  pi <- pmax(pi, cutoff)
  return(pi)
}

