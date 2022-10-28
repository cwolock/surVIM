estimate_cindex_newxfit_newapprox <- function(time,
                                    event,
                                    approx_times,
                                    preds,
                                    S_hat,
                                    G_hat,
                                    tau,
                                    preds_holdout,
                                    S_hat_holdout){
  sigma <- 0.1

  approxLoss <- function(x, sigma) { ## sigmoid function for loss
    1 / (1 + exp(x / sigma))
  }

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

  calc_KM_if <- function(i){
    t0 <- approx_times[i]
    S_hat_i <- S_hat[,i]
    inner.func.1 <- ifelse(time <= t0 & event == 1, 1/(S_hat_Y * G_hat_Y), 0 )
    inner.func.2 <- int.vals[,i]
    KM.if <- -S_hat_i * ( inner.func.1 - inner.func.2)
    return(KM.if)
  }

  # need to calculate KM IF at every point in approx_times so we can integrate
  KM_ifs <- matrix(unlist(lapply(1:length(approx_times), FUN = calc_KM_if)),
                   nrow = n)

  # approx_times_k <- approx_times[approx_times <= tau]
  # approx_times_k <- approx_times_k[-1]
  S_hat_k <- S_hat[,approx_times <= tau]
  S_hat_k_holdout <- S_hat_holdout[,approx_times <= tau]
  KM_ifs_k <- KM_ifs[,approx_times <= tau]

  k <- length(approx_times)

  # calc_phi_01 <- function(j){
  #   fx <- preds[j]
  #   varphi_x <- KM_ifs_k[j,]
  #   exceed_probs1 <- -rowSums(sweep(S_hat_k_holdout[,-k], MARGIN=2, diff(varphi_x), `*`))
  #   exceed_probs2 <- -rowSums(sweep(t(diff(t(S_hat_k_holdout))), MARGIN=2, varphi_x[-k], `*`))
  #   int <- mean(ifelse(fx > preds_holdout, 1, 0) * exceed_probs1 + ifelse(preds_holdout >= fx, 1, 0) * exceed_probs2)
  #   return(int)
  # }

  calc_phi_01 <- function(j){
    fx <- preds[j]
    varphi_x <- KM_ifs_k[j,]
    exceed_probs1 <- -rowSums(sweep(S_hat_k_holdout[,-k], MARGIN=2, diff(varphi_x), `*`))
    exceed_probs2 <- -rowSums(sweep(t(diff(t(S_hat_k_holdout))), MARGIN=2, varphi_x[-k], `*`))
    int <- mean((1 / (1 + exp(preds_holdout - fx / sigma)))* exceed_probs1 +
                  (1 / (1 + exp(fx - preds_holdout / sigma)))* exceed_probs2)
    return(int)
  }

  phi_01 <- unlist(lapply(1:n, FUN = calc_phi_01))

  # calc_phi_tilde_01 <- function(j){
  #   fx <- preds[j]
  #   Sx <- S_hat_k[j,]
  #   exceed_probs1 <- -rowSums(sweep(S_hat_k_holdout[,-k], MARGIN=2, diff(Sx), `*`))
  #   exceed_probs2 <- -rowSums(sweep(t(diff(t(S_hat_k_holdout))), MARGIN=2, Sx[-k], `*`))
  #   int <- mean(ifelse(fx > preds_holdout, 1, 0) * exceed_probs1 +
  #                 ifelse(preds_holdout >= fx, 1, 0) * exceed_probs2)
  #   return(int)
  # }

  calc_phi_tilde_01 <- function(j){
    fx <- preds[j]
    Sx <- S_hat_k[j,]
    exceed_probs1 <- -rowSums(sweep(S_hat_k_holdout[,-k], MARGIN=2, diff(Sx), `*`))
    exceed_probs2 <- -rowSums(sweep(t(diff(t(S_hat_k_holdout))), MARGIN=2, Sx[-k], `*`))
    int <- mean((1 / (1 + exp(preds_holdout - fx / sigma)))* exceed_probs1 +
                  (1 / (1 + exp(fx - preds_holdout / sigma))) * exceed_probs2)
    return(int)
  }

  phi_tilde_01_uncentered <- unlist(lapply(1:n, FUN = calc_phi_tilde_01))

  calc_phi_02 <- function(j){
    varphi_x <- KM_ifs_k[j,]
    exceed_probs1 <- -rowSums(sweep(S_hat_k_holdout[,-k], MARGIN=2, diff(varphi_x), `*`))
    exceed_probs2 <- -rowSums(sweep(t(diff(t(S_hat_k_holdout))), MARGIN=2, varphi_x[-k], `*`))
    int <- mean(exceed_probs1 + exceed_probs2)
    return(int)
  }

  phi_02 <- unlist(lapply(1:n, FUN = calc_phi_02))

  calc_phi_tilde_02 <- function(j){
    Sx <- S_hat_k[j,]
    exceed_probs1 <- -rowSums(sweep(S_hat_k_holdout[,-k], MARGIN=2, diff(Sx), `*`))
    exceed_probs2 <- -rowSums(sweep(t(diff(t(S_hat_k_holdout))), MARGIN=2, Sx[-k], `*`))
    int <- mean(exceed_probs1 + exceed_probs2)
    return(int)
  }

  phi_tilde_02_uncentered <- unlist(lapply(1:n, FUN = calc_phi_tilde_02))

  # calc_phi_02 <- function(j){
  #   varphi_x <- KM_ifs_k[j,length(approx_times)]
  #   return(varphi_x)
  # }
  # phi_02_new_new <- unlist(lapply(1:n, FUN = calc_phi_02))
  #
  # calc_phi_tilde_02 <- function(j){
  #   Sx <- S_hat_k[j,length(approx_times)]
  #   return(Sx)
  # }
  #
  # #phi_tilde_02_uncentered_new_new <- 0.5 - 0.5*(mean(S_hat_k[,length(approx_times)]))^2
  # phi_tilde_02_uncentered_new_new <- unlist(lapply(1:n, FUN = calc_phi_tilde_02))

  phi_tilde_01 <- phi_tilde_01_uncentered - mean(phi_tilde_01_uncentered)
  phi_tilde_02 <- phi_tilde_02_uncentered - mean(phi_tilde_02_uncentered)
  # phi_tilde_01_new <- phi_tilde_01_uncentered_new - mean(phi_tilde_01_uncentered_new)
  # phi_tilde_02_new <- phi_tilde_02_uncentered_new - mean(phi_tilde_02_uncentered_new)
  # phi_tilde_02_new_new <- phi_tilde_02_uncentered_new_new - mean(phi_tilde_02_uncentered_new_new)
  # note that including phi_tilde in the correction is a bit silly because we've
  # centered that term - it has no contribution to the correction (nor should it)
  # if we're underestimating, (both one-step and plug-in), then that means
  # the ratio V_1/V_2 is too small. So this gets amplfied in the one-step (we subtract
  # off something larger, making the estimate smaller)
  V_1 <- mean(phi_tilde_01_uncentered)/2
  V_2 <- mean(phi_tilde_02_uncentered)/2
  # V_1_new <- mean(phi_tilde_01_uncentered_new)/2
  # V_2_new <- mean(phi_tilde_02_uncentered_new)/2
  # V_2_new_new <- mean(phi_tilde_02_uncentered_new_new)
  # if_func <- (phi_01 + phi_tilde_01)/V_2 - V_1/(V_2^2)*(phi_02 + phi_tilde_02)
  # plug_in <- V_1/V_2
  # one_step <- V_1/V_2 + mean(if_func)

  if_func <- (phi_01 + phi_tilde_01)/V_2 - V_1/(V_2^2)*(phi_02 + phi_tilde_02)
  plug_in <- V_1/V_2
  one_step <- V_1/V_2 + mean(if_func)
#
#   if_func_new_new <- 2*(phi_01_new + phi_tilde_01_new)/(1-V_2_new_new^2) +
#     4*V_2_new_new*V_1_new/(1-V_2_new_new^2)*(phi_02_new_new + phi_tilde_02_new_new)
#   plug_in_new_new <- 2*V_1_new/(1 - V_2_new_new^2)
#   one_step_new_new <- plug_in_new_new + mean(if_func_new_new)

  return(list(one_step = one_step,
              plug_in = plug_in,
              if_func = if_func))
}
