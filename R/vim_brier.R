vim_brier <- function(time,
                      event,
                      X,
                      X_reduced,
                      approx_times,
                      landmark_times){
  event.SL.library <- cens.SL.library <- c("survSL.km",
                                           "survSL.coxph",
                                           "survSL.expreg",
                                           "survSL.weibreg",
                                           "survSL.loglogreg",
                                           "survSL.gam",
                                           "survSL.rfsrc")

  # time <- train$y
  # event <- train$delta
  # X <- train[,1:dimension]
  # X_reduced <- train[,2:dimension]
  # approx_times <- approx_times
  # landmark_times <- landmark_times

  # NOTE: Ted's function breaks with just a single covariate
  # also there are weird namespace issues with predict method
  # also why is new.times a mandatory argument
  # currently we assume approx_times include landmark times
  fit <- survSuperLearner::survSuperLearner(time = time,
                                            event = event,
                                            X = X,
                                            newX = X,
                                            new.times = c(1),
                                            event.SL.library = event.SL.library,
                                            cens.SL.library = cens.SL.library,
                                            verbose = FALSE,
                                            obsWeights = NULL,
                                            control = list(initWeightAlg = "survSL.rfsrc"))
  # oracle prediction, as well as nuisance predictions
  f_hat <- 1-survSuperLearner:::predict.survSuperLearner(fit,
                                                       newdata = X,
                                                       new.times = landmark_times)$event.SL.predict
  nuis_preds <- survSuperLearner:::predict.survSuperLearner(fit,
                                                            newdata = X,
                                                            new.times = approx_times)
  S_hat <- nuis_preds$event.SL.predict
  G_hat <- nuis_preds$cens.SL.predict
  fit_reduced <- survSuperLearner::survSuperLearner(time = time,
                                            event = event,
                                            X = X_reduced,
                                            newX = X_reduced,
                                            new.times = c(1),
                                            event.SL.library = event.SL.library,
                                            cens.SL.library = cens.SL.library,
                                            verbose = FALSE,
                                            obsWeights = NULL,
                                            control = list(initWeightAlg = "survSL.rfsrc"))
  # reduced oracle predictions
  fs_hat <- 1-survSuperLearner:::predict.survSuperLearner(fit_reduced,
                                                        newdata = X_reduced,
                                                        new.times = landmark_times)$event.SL.predict

  n <- length(time)
  # calculate integral for KM influence function
  # check left vs right continuous stuff here
  int.vals <- t(sapply(1:n, function(i) {
    vals <- diff(1/S_hat[i,])* 1/ G_hat[i,-ncol(G_hat)]
    if(any(approx_times[-1] > time[i])){
      vals[approx_times[-1] > time[i]] <- 0
    }
    c(0,cumsum(vals))
  }))
  # again, check left v right continuous here, especially with right = TRUE/FALSE
  S_hat_Y <- sapply(1:n, function(i) stepfun(approx_times, c(1,S_hat[i,]), right = FALSE)(time[i]))
  G_hat_Y <- sapply(1:n, function(i) stepfun(approx_times, c(1,G_hat[i,]), right = TRUE)(time[i]))
  # places to hold the value of the influence function, as well as the actual estimate
  IF.vals <- matrix(NA, nrow=n, ncol=length(landmark_times))
  brier <- rep(NA, length(landmark_times))

  for(i in 1:length(landmark_times)) {
    t0 <- landmark_times[i]
    k <- min(which(approx_times >= t0))
    F_hat_k <- 1-S_hat[,k]
    S_hat_k <- S_hat[,k]
    f_hat_k <- f_hat[,i]
    fs_hat_k <- fs_hat[,i]
    inner.func.1 <- ifelse(time <= t0 & event == 1, 1/(S_hat_Y * G_hat_Y), 0 )
    inner.func.2 <- int.vals[,k]
    KM.if <- -S_hat_k * ( inner.func.1 - inner.func.2)

    phi0 <- 2 * KM.if * (f_hat_k - F_hat_k)
    phi_tilde_0 <- -(f_hat_k - F_hat_k)^2 - mean(-(f_hat_k - F_hat_k)^2)

    phi0s <- 2 * KM.if * (fs_hat_k - F_hat_k)
    phi_tilde_0s <- -(fs_hat_k - F_hat_k)^2 - mean(-(fs_hat_k - F_hat_k)^2)

    if.func <- phi0 + phi_tilde_0 - phi0s - phi_tilde_0s

    brier[i] <- mean(-(f_hat_k - F_hat_k)^2) - mean(-(fs_hat_k - F_hat_k)^2) + mean(if.func)
    IF.vals[,i] <- if.func
  }

  return(data.frame(t = landmark_times,
                    brier = brier))#,
                    #IF.vals = IF.vals))
}
