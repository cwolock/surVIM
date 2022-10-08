estimate_cindex <- function(time,
                            event,
                            approx_times,
                            preds,
                            S_hat,
                            G_hat,
                            tau){

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

  approx_times_k <- approx_times[approx_times <= tau]
  approx_times_k <- approx_times_k[-1]
  S_hat_k <- S_hat[,approx_times <= tau]
  #G_hat_k <- G_hat[,approx_times <= tau]
  KM_ifs_k <- KM_ifs[,approx_times <= tau]

  calc_phi_01 <- function(j){
    fx <- preds[j]
    varphi_x <- KM_ifs_k[j,]
    exceed_probs1 <- -rowSums(sweep(S_hat_k[,-1], MARGIN=2, diff(varphi_x), `*`))
    exceed_probs2 <- -rowSums(sweep(t(diff(t(S_hat_k))), MARGIN=2, varphi_x[-1], `*`))
    # int <- -mean(ifelse(fx > preds, 1, 0) * rowSums(S_hat_k[,-1] * diff(varphi_x)) +
    #                ifelse(preds > fx, 1, 0) *  rowSums(varphi_x[-1] * t(diff(t(S_hat_k)))))
    #int <- mean(ifelse(fx > preds, 1, 0) * exceed_probs - ifelse(preds > fx, 1, 0) * exceed_probs) # old version
    int <- mean(ifelse(fx > preds, 1, 0) * exceed_probs1 + ifelse(preds > fx, 1, 0) * exceed_probs2)
    return(int)
  }

  phi_01 <- unlist(lapply(1:n, FUN = calc_phi_01))

  calc_phi_tilde_01 <- function(j){
    #tx <- time[j]
    #deltax <- event[j]
    #l <- min(which(approx_times >= tx))
    #wx <- G_hat_k[j,l]
    fx <- preds[j]
    Sx <- S_hat_k[j,]
    # int1 <- mean(ifelse(fx > preds & tx < time, 1, 0) * deltax / (wx * G_hat_k[,l]))
    # int2 <- mean(ifelse(tx < time, 1, 0) * deltax / (wx * G_hat_k[,l]))# * rowSums(S_hat_k[,-1] * diff(Sx))) #+
    #                #ifelse(preds > fx, 1, 0) *  rowSums(Sx[-1] * t(diff(t(S_hat_k)))))
    exceed_probs1 <- -rowSums(sweep(S_hat_k[,-1], MARGIN=2, diff(Sx), `*`)) # old one, but doesn't properly account for tau
    exceed_probs2 <- -rowSums(sweep(t(diff(t(S_hat_k))), MARGIN=2, Sx[-1], `*`))
    #exceed_probs <-  -rowSums(S_hat_k[,-1] * diff(Sx))
    int <- mean(ifelse(fx > preds, 1, 0) * exceed_probs1 +
                  ifelse(preds > fx, 1, 0) * exceed_probs2)#-rowSums(Sx[-1] * t(diff(t(S_hat_k)))))
    # -rowSums(S_hat_k[,-1] * diff(Sx))
    # -rowSums(sweep(t(diff(t(S_hat_k))), MARGIN=2, Sx[-1], `*`))
    # -rowSums(Sx[-1] * t(diff(t(S_hat_k))))
    return(int)
  }

  phi_tilde_01_uncentered <- unlist(lapply(1:n, FUN = calc_phi_tilde_01))

  calc_phi_02 <- function(j){
    varphi_x <- KM_ifs_k[j,]
    exceed_probs1 <- -rowSums(sweep(S_hat_k[,-1], MARGIN=2, diff(varphi_x), `*`))
    exceed_probs2 <- -rowSums(sweep(t(diff(t(S_hat_k))), MARGIN=2, varphi_x[-1], `*`))
    int <- mean(exceed_probs1 + exceed_probs2)
    return(int)
  }

  phi_02 <- unlist(lapply(1:n, FUN = calc_phi_02))

  calc_phi_tilde_02 <- function(j){
    #tx <- time[j]
    #deltax <- event[j]
    #l <- min(which(approx_times >= tx))
    #wx <- G_hat_k[j,l]
    Sx <- S_hat_k[j,]
    exceed_probs1 <- -rowSums(sweep(S_hat_k[,-1], MARGIN=2, diff(Sx), `*`)) # old one, but doesn't properly account for tau
    exceed_probs2 <- -rowSums(sweep(t(diff(t(S_hat_k))), MARGIN=2, Sx[-1], `*`))
    #exceed_probs <-  -rowSums(S_hat_k[,-1] * diff(Sx))
    int <- mean(exceed_probs1 + exceed_probs2)#-rowSums(Sx[-1] * t(diff(t(S_hat_k)))))
    return(int)
  }

  phi_tilde_02_uncentered <- unlist(lapply(1:n, FUN = calc_phi_tilde_02))


  phi_tilde_01 <- phi_tilde_01_uncentered - mean(phi_tilde_01_uncentered)
  phi_tilde_02 <- phi_tilde_02_uncentered - mean(phi_tilde_02_uncentered)
  V_1 <- mean(phi_tilde_01_uncentered)/2
  V_2 <- mean(phi_tilde_02_uncentered)/2
  if_func <- (phi_01 + phi_tilde_01)/V_2 - V_1/(V_2^2)*(phi_02 + phi_tilde_02)
  plug_in <- V_1/V_2
  one_step <- V_1/V_2 + mean(if_func)

  # plug_in <- mean(phi_tilde_01)#/mean(phi_tilde_01_b)
  # # UnoC(Surv(time, event), Surv(dtest$time, dtest$event), lpnew = preds)
  #
  # phi_tilde_01 <- phi_tilde_01 - plug_in
  #
  # if_func <- phi_01 + phi_tilde_01
  #
  # one_step <- plug_in + mean(if_func)

  return(list(one_step = one_step,
              plug_in = plug_in,
              if_func = if_func))
}
