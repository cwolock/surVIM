estimate_cindex <- function(time,
                            event,
                            approx_times,
                            preds,
                            S_hat,
                            G_hat,
                            tau){

  n <- length(time)

  KM_IFs <- calc_KM_IF(time = time,
                       event = event,
                       S_hat = S_hat,
                       G_hat = G_hat,
                       approx_times = approx_times)

  S_hat_k <- S_hat[,approx_times <= tau]
  KM_IFs_k <- KM_IFs[,approx_times <= tau]

  k <- length(approx_times)

  calc_phi_01 <- function(j){
    fx <- preds[j]
    varphi_x <- KM_IFs_k[j,]
    exceed_probs1 <- -rowSums(sweep(S_hat_k[,-k], MARGIN=2, diff(varphi_x), `*`))
    exceed_probs2 <- -rowSums(sweep(t(diff(t(S_hat_k))), MARGIN=2, varphi_x[-k], `*`))
    int <- mean(ifelse(fx > preds, 1, 0)* exceed_probs1 + ifelse(preds > fx, 1, 0)* exceed_probs2)
    return(int)
  }

  calc_phi_01_extra <- function(j){
    fx <- preds[j]
    varphi_x <- KM_IFs_k[j,]
    exceed_probs1 <- -rowSums(sweep(KM_IFs_k[,-k], MARGIN=2, diff(varphi_x), `*`))
    exceed_probs2 <- -rowSums(sweep(t(diff(t(KM_IFs_k))), MARGIN=2, varphi_x[-k], `*`))
    int <- mean(ifelse(fx > preds, 1, 0)* exceed_probs1 + ifelse(preds > fx, 1, 0)* exceed_probs2)
    return(int)
  }

  calc_phi_tilde_01 <- function(j){
    fx <- preds[j]
    Sx <- S_hat_k[j,]
    exceed_probs1 <- -rowSums(sweep(S_hat_k[,-k], MARGIN=2, diff(Sx), `*`))
    exceed_probs2 <- -rowSums(sweep(t(diff(t(S_hat_k))), MARGIN=2, Sx[-k], `*`))
    int <- mean(ifelse(fx > preds, 1, 0)* exceed_probs1 + ifelse(preds > fx, 1, 0)* exceed_probs2)
    return(int)
  }

  calc_phi_02 <- function(j){
    varphi_x <- KM_IFs_k[j,]
    exceed_probs1 <- -rowSums(sweep(S_hat_k[,-k], MARGIN=2, diff(varphi_x), `*`))
    exceed_probs2 <- -rowSums(sweep(t(diff(t(S_hat_k))), MARGIN=2, varphi_x[-k], `*`))
    int <- mean(exceed_probs1 + exceed_probs2)
    return(int)
  }

  calc_phi_02_extra <- function(j){
    varphi_x <- KM_IFs_k[j,]
    exceed_probs1 <- -rowSums(sweep(KM_IFs_k[,-k], MARGIN=2, diff(varphi_x), `*`))
    exceed_probs2 <- -rowSums(sweep(t(diff(t(KM_IFs_k))), MARGIN=2, varphi_x[-k], `*`))
    int <- mean(exceed_probs1 + exceed_probs2)
    return(int)
  }

  calc_phi_tilde_02 <- function(j){
    Sx <- S_hat_k[j,]
    exceed_probs1 <- -rowSums(sweep(S_hat_k[,-k], MARGIN=2, diff(Sx), `*`))
    exceed_probs2 <- -rowSums(sweep(t(diff(t(S_hat_k))), MARGIN=2, Sx[-k], `*`))
    int <- mean(exceed_probs1 + exceed_probs2)
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
              EIF = EIF))
}
