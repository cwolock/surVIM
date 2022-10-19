estimate_AUC_newxfit <- function(time,
                         event,
                         approx_times,
                         t0,
                         preds,
                         S_hat,
                         G_hat,
                         preds_holdout,
                         S_hat_holdout){

  # time <- holdout$y
  # event <- holdout$delta
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

  k <- min(which(approx_times >= t0))
  S_hat_k <- S_hat[,k]
  S_hat_k_holdout <- S_hat_holdout[,k]
  G_hat_k <- G_hat[,k]
  f_hat_k <- preds # oracle is CDF, not survival
  f_hat_k_holdout <- preds_holdout
  inner.func.1 <- ifelse(time <= t0 & event == 1, 1/(S_hat_Y * G_hat_Y), 0 )
  inner.func.2 <- int.vals[,k]
  KM.if <- -S_hat_k * ( inner.func.1 - inner.func.2)

  calc_phi_01 <- function(j){
    fx <- f_hat_k[j]
    varphi_x <- KM.if[j]
    int <- mean(ifelse(f_hat_k_holdout > fx, 1, 0) * (1 - S_hat_k_holdout) -
                  ifelse(f_hat_k_holdout <= fx, 1, 0) * S_hat_k_holdout)
    return(varphi_x*int)
  }

  calc_phi_tilde_01 <- function(j){
    fx <- f_hat_k[j]
    Sx <- S_hat_k[j]
    int <- mean(ifelse(f_hat_k_holdout > fx, 1, 0) * (1 - S_hat_k_holdout) * Sx +
                  ifelse(f_hat_k_holdout <= fx, 1, 0) * S_hat_k_holdout * (1 - Sx))
    return(int)
  }

  phi_01 <- unlist(lapply(1:n, FUN = calc_phi_01))

  phi_tilde_01_uncentered <- unlist(lapply(1:n, FUN = calc_phi_tilde_01))

  phi_tilde_01 <- phi_tilde_01_uncentered - mean(phi_tilde_01_uncentered)

  calc_phi_02 <- function(j){
    varphi_x <- KM.if[j]
    int <- mean((1 - S_hat_k_holdout) - S_hat_k_holdout)
    return(varphi_x * int)
  }

  calc_phi_tilde_02 <- function(j){
    Sx <- S_hat_k[j]
    int <- mean((1 - S_hat_k_holdout) * Sx + S_hat_k_holdout * (1 - Sx))
    return(int)
  }

  phi_02 <- unlist(lapply(1:n, FUN = calc_phi_02))

  phi_tilde_02_uncentered <- unlist(lapply(1:n, FUN = calc_phi_tilde_02))

  phi_tilde_02 <- phi_tilde_02_uncentered - mean(phi_tilde_02_uncentered)

  V_1 <- mean(phi_tilde_01_uncentered)/2
  V_2 <- mean(phi_tilde_02_uncentered)/2

  if_func <- (phi_01 + phi_tilde_01)/V_2 - V_1/(V_2^2)*(phi_02 + phi_tilde_02)
  plug_in <- V_1/V_2
  one_step <- V_1/V_2 + mean(if_func)

  return(list(one_step = one_step,
              plug_in = plug_in,
              if_func = if_func))
}
