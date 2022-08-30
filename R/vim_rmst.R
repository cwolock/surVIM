vim_rmst <- function(time,
                      event,
                      X,
                      X_reduced,
                      approx_times,
                      f_hat,
                      fs_hat,
                      S_hat,
                      G_hat,
                      holdout,
                      taus){

  dimension <- 4
  time_holdout <- holdout$y
  event_holdout <- holdout$delta
  X_holdout <-holdout[,1:dimension]
  X_reduced_holdout <- holdout[,2:dimension]


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

  calc_KM_if <- function(i){
    t0 <- approx_times[i]
    S_hat_i <- S_hat[,i]
    inner.func.1 <- ifelse(time <= t0 & event == 1, 1/(S_hat_Y * G_hat_Y), 0 )
    inner.func.2 <- int.vals[,i]
    KM.if <- -S_hat_i * ( inner.func.1 - inner.func.2)
    return(KM.if)
  }

  KM_ifs <- matrix(unlist(lapply(1:length(approx_times), FUN = calc_KM_if)),
                   nrow = n)

  plug_ins <- rep(NA, length(taus))
  one_steps <- rep(NA, length(taus))
  var_ests <- rep(NA, length(taus))
  for (tau in taus){

    f_k <- f_hat[,which(taus == tau)]

    t_wedge <- ifelse(approx_times <= tau, approx_times, tau)
    t_wedge <- t_wedge[-1]

    calc_phi_01 <- function(i){
      (rep(f_k[i], length(approx_times)-1) - t_wedge)^2 * diff(KM_ifs[i,])
    }

    phi_01 <- matrix(unlist(lapply(1:n, FUN = calc_phi_01)),
                     nrow = n)

    phi_01 <- rowSums(phi_01)

    calc_phi_tilde_01 <- function(i){
      (rep(f_k[i], length(approx_times)-1) - t_wedge)^2 * diff(S_hat[i,])
    }

    phi_tilde_01 <- matrix(unlist(lapply(1:n, FUN = calc_phi_tilde_01)),
                           nrow = n)
    phi_tilde_01 <- rowSums(phi_tilde_01)

    plug_in <- mean(phi_tilde_01)

    phi_tilde_01 <- phi_tilde_01 - plug_in

    if.func <- phi_01 + phi_tilde_01

    one_step <- plug_in + mean(if.func)

    var_est <- mean(if.func^2)

    plug_ins[which(taus == tau)] <- plug_in
    one_steps[which(taus == tau)] <- one_step
    var_ests[which(taus == tau)] <- var_est
  }

  return(data.frame(tau = taus,
                    plug_in = plug_ins,
                    one_step = one_steps,
                    var_est = var_ests))
}
