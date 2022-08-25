vim_brier <- function(time,
                      event,
                      X,
                      X_reduced,
                      approx_times,
                      landmark_times,
                      nuisance,
                      holdout){
  event.SL.library <- cens.SL.library <- c("survSL.km",
                                           "survSL.coxph",
                                           "survSL.expreg",
                                           "survSL.weibreg",
                                           "survSL.loglogreg",
                                           "survSL.gam",
                                           "survSL.rfsrc")

  time <- train$y
  event <- train$delta
  X <- train[,1:dimension]
  X_reduced <- train[,2:dimension]
  approx_times <- approx_times
  landmark_times <- landmark_times

  # NOTE: Ted's function breaks with just a single covariate
  # also there are weird namespace issues with predict method
  # also why is new.times a mandatory argument
  # currently we assume approx_times include landmark times
  dimension <- 4
  time_holdout <- holdout$y
  event_holdout <- holdout$delta
  X_holdout <-holdout[,1:dimension]
  X_reduced_holdout <- holdout[,2:dimension]

  if (nuisance == "survSL"){
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
    # reverse_fit <- survSuperLearner::survSuperLearner(time = time,
    #                                           event = 1-event,
    #                                           X = X,
    #                                           newX = X,
    #                                           new.times = c(1),
    #                                           event.SL.library = event.SL.library,
    #                                           cens.SL.library = cens.SL.library,
    #                                           verbose = FALSE,
    #                                           obsWeights = NULL,
    #                                           control = list(initWeightAlg = "survSL.rfsrc"))
    # oracle prediction, as well as nuisance predictions
    f_hat <- survSuperLearner:::predict.survSuperLearner(fit,
                                                         newdata = X_holdout,
                                                         new.times = landmark_times)$event.SL.predict
    nuis_preds <- survSuperLearner:::predict.survSuperLearner(fit,
                                                              newdata = X_holdout,
                                                              new.times = approx_times)
    # nuis_preds_reverse <- survSuperLearner:::predict.survSuperLearner(reverse_fit,
    #                                                           newdata = X,
    #                                                           new.times = approx_times)
    S_hat <- nuis_preds$event.SL.predict
    G_hat <- nuis_preds$cens.SL.predict
    # S_hat_reverse <- nuis_preds_reverse$cens.SL.predict
    # G_hat_reverse <- nuis_preds_reverse$event.SL.predict
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
    fs_hat <- survSuperLearner:::predict.survSuperLearner(fit_reduced,
                                                          newdata = X_reduced_holdout,
                                                          new.times = landmark_times)$event.SL.predict
  } else if (nuisance == "param_AFT"){
    df_full <- data.frame(time = time,
                          event = event,
                          X)
    df_reduced <- data.frame(time = time,
                             event = event,
                             X_reduced)
    event_fit <- survival::survreg(survival::Surv(time, event) ~ ., data = df_full,
                                   dist = "lognormal")
    cens_fit <- survival::survreg(survival::Surv(time, 1-event) ~ ., data = df_full,
                                  dist = "lognormal")
    event_fit_reduced <- survival::survreg(survival::Surv(time, event) ~ ., data = df_reduced,
                                           dist = "lognormal")
    event_q_pred <- predict(event_fit,
                            newdata = X_holdout,
                            type = 'quantile', p = seq(0, .999, by=.001))
    cens_q_pred <- predict(cens_fit,
                           newdata = X_holdout,
                           type = 'quantile', p = seq(0, .999, by=.001))
    event_q_pred_reduced <- predict(event_fit_reduced,
                                    newdata = X_reduced_holdout,
                                    type = 'quantile', p = seq(0, .999, by=.001))
    # this is here to handle exact 0 times - see survSuperLearner code
    pos.pred <- rep(1, nrow(X))
    f_hat <- try(t(sapply(1:nrow(event_q_pred), function(j) {
      pos.pred[j] * (1-stats::approx(event_q_pred[j,], seq(0, .999, by=.001),
                                     xout = landmark_times,
                                     method = 'linear',
                                     rule = 2)$y)
    })), silent = TRUE)
    fs_hat <- try(t(sapply(1:nrow(event_q_pred_reduced), function(j) {
      pos.pred[j] * (1-stats::approx(event_q_pred_reduced[j,], seq(0, .999, by=.001),
                                     xout = landmark_times,
                                     method = 'linear',
                                     rule = 2)$y)
    })), silent = TRUE)
    S_hat <- try(t(sapply(1:nrow(event_q_pred), function(j) {
      pos.pred[j] * (1-stats::approx(event_q_pred[j,], seq(0, .999, by=.001),
                                     xout = approx_times,
                                     method = 'linear',
                                     rule = 2)$y)
    })), silent = TRUE)
    G_hat <- try(t(sapply(1:nrow(cens_q_pred), function(j) {
      pos.pred[j] * (1-stats::approx(cens_q_pred[j,], seq(0, .999, by=.001),
                                     xout = approx_times,
                                     method = 'linear',
                                     rule = 2)$y)
    })), silent = TRUE)
  }

  time <- time_holdout
  event <- event_holdout
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
  #IF.vals <- matrix(NA, nrow=n, ncol=length(landmark_times))
  brier <- rep(NA, length(landmark_times))
  brier_plug <- rep(NA, length(landmark_times))
  S_t <- rep(NA, length(landmark_times))
  G_t <- rep(NA, length(landmark_times))
  var_est <- rep(NA, length(landmark_times))
  for(i in 1:length(landmark_times)) {
    t0 <- landmark_times[i]
    k <- min(which(approx_times >= t0))
    #F_hat_k <- 1-S_hat[,k]
    S_hat_k <- S_hat[,k]
    G_hat_k <- G_hat[,k]
    f_hat_k <- f_hat[,i]
    fs_hat_k <- fs_hat[,i]
    inner.func.1 <- ifelse(time <= t0 & event == 1, 1/(S_hat_Y * G_hat_Y), 0 )
    inner.func.2 <- int.vals[,k]
    KM.if <- -S_hat_k * ( inner.func.1 - inner.func.2)

    # should phi_0 and phi_0s be negative or positive...
    phi0 <- 2 * KM.if * (f_hat_k - S_hat_k)
    phi_tilde_0 <- -(f_hat_k - S_hat_k)^2 - mean(-(f_hat_k - S_hat_k)^2)

    phi0s <- 2 * KM.if * (fs_hat_k - S_hat_k)
    phi_tilde_0s <- -(fs_hat_k - S_hat_k)^2 - mean(-(fs_hat_k - S_hat_k)^2)

    # eta0 <- -mean(S_hat_k*(1 - S_hat_k))
    # eta0_if_tilde <- -S_hat_k*(1 - S_hat_k) - mean(-S_hat_k*(1 - S_hat_k))
    # eta0_if <- -KM.if + 2*KM.if*S_hat_k

    if.func <- phi0 + phi_tilde_0 - phi0s - phi_tilde_0s
    # if.func_eta <- eta0_if_tilde + eta0_if
    # mu <- mean(if.func)
    # sigma <- sd(if.func)
    # if.func.z <- (if.func-mu)/sigma
    # if.func <- if.func[if.func.z > -4 & if.func.z < 4]
    # print(sum(if.func.z < -5 | if.func.z > 5)/length(if.func.z)*100)

    brier[i] <- mean(-(f_hat_k - S_hat_k)^2) - mean(-(fs_hat_k - S_hat_k)^2) + mean(if.func)
    brier_plug[i] <- mean(-(f_hat_k - S_hat_k)^2) - mean(-(fs_hat_k - S_hat_k)^2)
    #IF.vals[,i] <- if.func
    S_t[i] <- mean(S_hat_k)
    G_t[i] <- mean(G_hat_k)
    var_est[i] <- mean(if.func^2)#var(if.func)
  }

  return(data.frame(t = landmark_times,
                    brier = brier,
                    brier_plug = brier_plug,
                    S_t = S_t,
                    G_t = G_t,
                    var_est = var_est))
}
