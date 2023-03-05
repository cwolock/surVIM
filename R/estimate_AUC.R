estimate_AUC<- function(time,
                        event,
                        approx_times,
                        tau,
                        preds,
                        S_hat,
                        G_hat){

  n <- length(time)

  KM_IFs <- calc_KM_IF(time = time,
                       event = event,
                       S_hat = S_hat,
                       G_hat = G_hat,
                       approx_times = approx_times)

  k <- min(which(approx_times >= tau))
  S_hat_k <- S_hat[,k]
  G_hat_k <- G_hat[,k]
  KM_IFs <- KM_IFs[,k]

  calc_phi_01 <- function(j){
    fx <- preds[j]
    varphi_x <- KM_IFs[j]
    int <- mean(ifelse(preds > fx, 1, 0) * (1 - S_hat_k) -
                  ifelse(preds <= fx, 1, 0) * S_hat_k)
    return(varphi_x*int)
  }

  calc_phi_01_extra <- function(j){
    fx <- preds[j]
    varphi_x <- KM_IFs[j]
    int <- mean(ifelse(preds > fx, 1, 0) * (- KM_IFs) -
                  ifelse(preds <= fx, 1, 0) * KM_IFs)
    return(varphi_x*int/2)
  }

  calc_phi_tilde_01 <- function(j){
    fx <- preds[j]
    Sx <- S_hat_k[j]
    varphi_x <- KM_IFs[j]
    int <- mean(ifelse(preds > fx, 1, 0) * (1 - S_hat_k) * (Sx) +
                  ifelse(preds <= fx, 1, 0) * (S_hat_k) * (1 - Sx))
    return(int)
  }



  calc_phi_02 <- function(j){
    varphi_x <- KM_IFs[j]
    int <- mean((1 - S_hat_k) - S_hat_k)
    return(varphi_x * int)
  }

  calc_phi_02_extra <- function(j){
    varphi_x <- KM_IFs[j]
    int <- mean((-KM_IFs) - KM_IFs)
    return(varphi_x * int/2)
  }

  calc_phi_tilde_02 <- function(j){
    Sx <- S_hat_k[j]
    int <- mean((1 - S_hat_k) * Sx + S_hat_k * (1 - Sx))
    return(int)
  }

  phi_01 <- unlist(lapply(1:n, FUN = calc_phi_01))
  phi_01_extra <- unlist(lapply(1:n, FUN = calc_phi_01_extra))
  phi_tilde_01_uncentered <- unlist(lapply(1:n, FUN = calc_phi_tilde_01))
  phi_02 <- unlist(lapply(1:n, FUN = calc_phi_02))
  phi_02_extra <- unlist(lapply(1:n, FUN = calc_phi_02_extra))
  phi_tilde_02_uncentered <- unlist(lapply(1:n, FUN = calc_phi_tilde_02))

  phi_tilde_01 <- phi_tilde_01_uncentered - mean(phi_tilde_01_uncentered)
  phi_tilde_02 <- phi_tilde_02_uncentered - mean(phi_tilde_02_uncentered)

  V_1 <- mean(phi_tilde_01_uncentered)/2
  V_2 <- mean(phi_tilde_02_uncentered)/2

  if_func_1 <- phi_01 + phi_tilde_01 + phi_01_extra
  if_func_2 <- phi_02 + phi_tilde_02 + phi_02_extra

  V_1_os <- V_1 + mean(if_func_1)
  V_2_os <- V_2 + mean(if_func_2)

  one_step <- V_1_os/V_2_os

  EIF <- (phi_01 + phi_tilde_01)/V_2 - V_1/(V_2^2)*(phi_02 + phi_tilde_02)
  plug_in <- V_1/V_2

  return(list(one_step = one_step,
              plug_in = plug_in,
              EIF = EIF,
              numerator = V_1_os,
              denominator = V_2_os))
}
