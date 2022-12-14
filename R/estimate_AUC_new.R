estimate_AUC_new<- function(time,
                            event,
                            approx_times,
                            t0,
                            preds,
                            S_hat,
                            G_hat,
                            preds_holdout,
                            S_hat_holdout,
                            time_holdout,
                            event_holdout,
                            G_hat_holdout){

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

  # n_extra <- length(time_holdout)
  # # calculate integral for KM influence function
  # # check left vs right continuous stuff here
  # int.vals_holdout <- t(sapply(1:n_extra, function(i) {
  #   vals <- diff(1/S_hat_holdout[i,])* 1/ G_hat_holdout[i,-ncol(G_hat_holdout)]
  #   if(any(approx_times[-1] > time_holdout[i])){
  #     vals[approx_times[-1] > time_holdout[i]] <- 0
  #   }
  #   c(0,cumsum(vals))
  # }))
  # # again, check left v right continuous here, especially with right = TRUE/FALSE
  # S_hat_Y_holdout <- sapply(1:n_extra, function(i) stepfun(approx_times, c(1,S_hat_holdout[i,]), right = FALSE)(time_holdout[i]))
  # G_hat_Y_holdout <- sapply(1:n_extra, function(i) stepfun(approx_times, c(1,G_hat_holdout[i,]), right = TRUE)(time_holdout[i]))

  k <- min(which(approx_times >= t0))
  S_hat_k <- S_hat[,k]
  S_hat_k_holdout <- S_hat_holdout[,k]
  G_hat_k <- G_hat[,k]
  f_hat_k <- preds # oracle is CDF, not survival
  f_hat_k_holdout <- preds_holdout


  inner.func.1 <- ifelse(time <= t0 & event == 1, 1/(S_hat_Y * G_hat_Y), 0 )
  inner.func.2 <- int.vals[,k]
  KM.if <- -S_hat_k * ( inner.func.1 - inner.func.2)
  #
  #
  #   inner.func.1 <- ifelse(time_holdout <= t0 & event_holdout == 1, 1/(S_hat_Y_holdout * G_hat_Y_holdout), 0 )
  #   inner.func.2 <- int.vals_holdout[,k]
  #   KM.if_holdout <- -S_hat_k_holdout * ( inner.func.1 - inner.func.2)

  #print(KM.if_holdout)

  # calc_phi_01 <- function(j){
  #   fx <- f_hat_k[j]
  #   varphi_x <- KM.if[j]
  #   int <- mean(ifelse(f_hat_k_holdout > fx, 1, 0) * (1 - S_hat_k_holdout) -
  #                 ifelse(f_hat_k_holdout <= fx, 1, 0) * S_hat_k_holdout)
  #   return(varphi_x*int)
  # }
  #
  # calc_phi_01_extra <- function(j){
  #   fx <- f_hat_k[j]
  #   varphi_x <- KM.if[j]
  #   int <- mean(ifelse(f_hat_k_holdout > fx, 1, 0) * (1 - KM.if_holdout) -
  #                 ifelse(f_hat_k_holdout <= fx, 1, 0) * KM.if_holdout)
  #   return(varphi_x*int/2)
  # }
  #
  # calc_phi_tilde_01 <- function(j){
  #   fx <- f_hat_k[j]
  #   Sx <- S_hat_k[j]
  #   int <- mean(ifelse(f_hat_k_holdout > fx, 1, 0) * (1 - S_hat_k_holdout) * Sx +
  #                 ifelse(f_hat_k_holdout <= fx, 1, 0) * S_hat_k_holdout * (1 - Sx))
  #   return(int)
  # }

  calc_phi_01 <- function(j){
    fx <- f_hat_k[j]
    varphi_x <- KM.if[j]
    int <- mean(ifelse(f_hat_k > fx, 1, 0) * (1 - S_hat_k) -
                  ifelse(f_hat_k <= fx, 1, 0) * S_hat_k)
    return(varphi_x*int)
  }

  calc_phi_01_extra <- function(j){
    fx <- f_hat_k[j]
    varphi_x <- KM.if[j]
    int <- mean(ifelse(f_hat_k > fx, 1, 0) * (- KM.if) -
                  ifelse(f_hat_k <= fx, 1, 0) * KM.if)
    return(varphi_x*int/2)
  }

  # calc_phi_bar_01 <- function(j){
  #   fx <- f_hat_k[j]
  #   Sx <- S_hat_k[j]
  #   varphi_x <- KM.if[j]
  #   int <- mean(ifelse(f_hat_k > fx, 1, 0) * (1 - S_hat_k - KM.if) * (Sx + varphi_x) +
  #                 ifelse(f_hat_k <= fx, 1, 0) * (S_hat_k + KM.if) * (1 - Sx - varphi_x))
  #   return(int)
  # }
  #
  # calc_phi_01 <- function(j){
  #   fx <- f_hat_k[j]
  #   Sx <- S_hat_k[j]
  #   varphi_x <- KM.if[j]
  #   int <- mean(ifelse(f_hat_k > fx, 1, 0) * (1 - S_hat_k) * (Sx + varphi_x) +
  #                 ifelse(f_hat_k <= fx, 1, 0) * (S_hat_k) * (1 - Sx - varphi_x))
  #   return(int)
  # }

  calc_phi_tilde_01 <- function(j){
    fx <- f_hat_k[j]
    Sx <- S_hat_k[j]
    varphi_x <- KM.if[j]
    int <- mean(ifelse(f_hat_k > fx, 1, 0) * (1 - S_hat_k) * (Sx) +
                  ifelse(f_hat_k <= fx, 1, 0) * (S_hat_k) * (1 - Sx))
    return(int)
  }
  #
  phi_01 <- unlist(lapply(1:n, FUN = calc_phi_01))
  #phi_bar_01 <- unlist(lapply(1:n, FUN = calc_phi_bar_01))
  #phi_tilde_01 <- unlist(lapply(1:n, FUN = calc_phi_tilde_01))
  phi_01_extra <- unlist(lapply(1:n, FUN = calc_phi_01_extra))

  phi_tilde_01_uncentered <- unlist(lapply(1:n, FUN = calc_phi_tilde_01))

  phi_tilde_01 <- phi_tilde_01_uncentered - mean(phi_tilde_01_uncentered)

  # calc_phi_02 <- function(j){
  #   varphi_x <- KM.if[j]
  #   int <- mean((1 - S_hat_k_holdout) - S_hat_k_holdout)
  #   return(varphi_x * int)
  # }
  #
  # calc_phi_02_extra <- function(j){
  #   varphi_x <- KM.if[j]
  #   int <- mean((1 - KM.if_holdout) - KM.if_holdout)
  #   return(varphi_x * int/2)
  # }
  #
  # calc_phi_tilde_02 <- function(j){
  #   Sx <- S_hat_k[j]
  #   int <- mean((1 - S_hat_k_holdout) * Sx + S_hat_k_holdout * (1 - Sx))
  #   return(int)
  # }

  calc_phi_02 <- function(j){
    varphi_x <- KM.if[j]
    int <- mean((1 - S_hat_k) - S_hat_k)
    return(varphi_x * int)
  }

  calc_phi_02_extra <- function(j){
    varphi_x <- KM.if[j]
    int <- mean((-KM.if) - KM.if)
    return(varphi_x * int/2)
  }

  # calc_phi_bar_02 <- function(j){
  #   Sx <- S_hat_k[j]
  #   varphi_x <- KM.if[j]
  #   int <- mean((1 - S_hat_k - KM.if) * (Sx+varphi_x) + (S_hat_k+KM.if) * (1 - Sx-varphi_x))
  #   return(int)
  # }
  #
  # calc_phi_02 <- function(j){
  #   Sx <- S_hat_k[j]
  #   varphi_x <- KM.if[j]
  #   int <- mean((1 - S_hat_k) * (Sx+varphi_x) + (S_hat_k) * (1 - Sx-varphi_x))
  #   return(int)
  # }

  calc_phi_tilde_02 <- function(j){
    Sx <- S_hat_k[j]
    int <- mean((1 - S_hat_k_holdout) * Sx + S_hat_k_holdout * (1 - Sx))
    return(int)
  }

  phi_02 <- unlist(lapply(1:n, FUN = calc_phi_02))
  phi_02_extra <- unlist(lapply(1:n, FUN = calc_phi_02_extra))
  phi_tilde_02_uncentered <- unlist(lapply(1:n, FUN = calc_phi_tilde_02))

  phi_tilde_02 <- phi_tilde_02_uncentered - mean(phi_tilde_02_uncentered)

  V_1 <- mean(phi_tilde_01_uncentered)/2
  V_2 <- mean(phi_tilde_02_uncentered)/2
  #
  if_func_1 <- phi_01 + phi_tilde_01 + phi_01_extra
  if_func_2 <- phi_02 + phi_tilde_02 + phi_02_extra
  #
  V_1_os <- V_1 + mean(if_func_1)
  V_2_os <- V_2 + mean(if_func_2)

  one_step <- V_1_os/V_2_os

  if_func <- (phi_01 + phi_tilde_01)/V_2 - V_1/(V_2^2)*(phi_02 + phi_tilde_02)
  plug_in <- V_1/V_2
  #one_step <- V_1/V_2 + mean(if_func)

  return(list(one_step = one_step,
              plug_in = plug_in,
              if_func = if_func,
              numerator = V_1_os,
              denominator = V_2_os))
}
