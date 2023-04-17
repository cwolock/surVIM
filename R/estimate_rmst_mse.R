estimate_rmst_mse <- function(time,
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
  pi_k <- S_hat_k + KM_IFs_k

  approx_times_k <- approx_times[approx_times <= tau]
  approx_times_k <- approx_times_k[-1]

  calc_phi_01_combined <- function(i){
    -sum((preds[i] - approx_times_k)^2 * diff(pi_k[i,])) +
      (preds[i] - tau)^2*pi_k[i,which(approx_times == tau)]
  }

  phi_01_combined <- unlist(lapply(1:n, FUN = calc_phi_01_combined))

  # k <- length(approx_times)
  # k <- sum(approx_times <= tau)



  # n <- length(time)
  # calculate integral for KM influence function
  # check left vs right continuous stuff here
  # int.vals <- t(sapply(1:n, function(i) {
  #   vals <- diff(1/S_hat[i,])* 1/ G_hat[i,-ncol(G_hat)]
  #   if(any(approx_times[-1] > time[i])){
  #     vals[approx_times[-1] > time[i]] <- 0
  #   }
  #   c(0,cumsum(vals))
  # }))
  # again, check left v right continuous here, especially with right = TRUE/FALSE
  # S_hat_Y <- sapply(1:n, function(i) stepfun(approx_times, c(1,S_hat[i,]), right = FALSE)(time[i]))
  # G_hat_Y <- sapply(1:n, function(i) stepfun(approx_times, c(1,G_hat[i,]), right = TRUE)(time[i]))
  # places to hold the value of the influence function, as well as the actual estimate

  # calc_KM_if <- function(i){
  #   t0 <- approx_times[i]
  #   S_hat_i <- S_hat[,i]
  #   inner.func.1 <- ifelse(time <= t0 & event == 1, 1/(S_hat_Y * G_hat_Y), 0 )
  #   inner.func.2 <- int.vals[,i]
  #   KM.if <- -S_hat_i * ( inner.func.1 - inner.func.2)
  #   return(KM.if)
  # }

  # need to calculate KM IF at every point in approx_times so we can integrate
  # KM_ifs <- matrix(unlist(lapply(1:length(approx_times), FUN = calc_KM_if)),
                   # nrow = n)

  # calc_phi_01 <- function(i){
  #   -sum(-(rep(preds[i], length(approx_times)-1)-approx_times[-1])^2*diff(KM_ifs[i,])*ifelse(approx_times[-1] <= tau, 1, 0)) +
  #     (preds[i] - tau)^2*KM_ifs[i,which(approx_times == tau)]
  # }

  # t_wedge <- ifelse(approx_times <= tau, approx_times, tau)
  # t_wedge <- t_wedge[-1]

  # calc_phi_01 <- function(i){
  #   sum(-(rep(f_k[i], length(approx_times)-1) - t_wedge)^2 * diff(KM_ifs[i,]))
  # }

  # calc_phi_tilde_01 <- function(i){
  #   sum(-(rep(f_k[i], length(approx_times)-1) - t_wedge)^2 * diff(S_hat[i,]))
  # }
#
#   approx_times_k <- approx_times[approx_times <= tau]
#   approx_times_k <- approx_times_k[-1]
#   S_hat_k <- S_hat[,approx_times <= tau]
#   KM_ifs_k <- KM_ifs[,approx_times <= tau]

  calc_phi_01 <- function(i){
    -sum((preds[i] - approx_times_k)^2 * diff(KM_IFs_k[i,])) + (preds[i] - tau)^2*KM_IFs[i,which(approx_times == tau)]
  }

  phi_01 <- unlist(lapply(1:n, FUN = calc_phi_01))

  calc_phi_tilde_01 <- function(i){
    -sum((preds[i] - approx_times_k)^2 * diff(S_hat_k[i,])) + (preds[i] - tau)^2*S_hat[i,which(approx_times == tau)]
  }
  phi_tilde_01 <- unlist(lapply(1:n, FUN = calc_phi_tilde_01))

  plug_in <- mean(phi_tilde_01)

  phi_tilde_01 <- phi_tilde_01 - plug_in

  if_func <- phi_01 + phi_tilde_01

  one_step <- mean(phi_01_combined)#plug_in + mean(if_func)
  return(list(one_step = one_step,
              plug_in = plug_in,
              EIF = if_func))
}
