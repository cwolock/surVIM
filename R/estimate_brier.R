estimate_brier <- function(time,
                           event,
                           approx_times,
                           t0,
                           preds,
                           S_hat,
                           G_hat){

  n <- length(time)
  # calculate integral for KM influence function
  # check left vs right continuous stuff here
  G_hat[G_hat < 0.01] <- 0.01 # truncate!
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

  k <- min(which(approx_times >= t0))
  S_hat_k <- S_hat[,k]
  G_hat_k <- G_hat[,k]
  inner.func.1 <- ifelse(time <= t0 & event == 1, 1/(S_hat_Y * G_hat_Y), 0 )
  inner.func.2 <- int.vals[,k]
  KM.if <- -S_hat_k * ( inner.func.1 - inner.func.2)
  f_hat_k <- preds

  # should phi_0 and phi_0s be negative or positive...
  # phi0_old <- 2 * KM.if * (f_hat_k - S_hat_k)
  # phi_tilde_0_old <- -(f_hat_k - S_hat_k)^2 - mean(-(f_hat_k - S_hat_k)^2)

  phi0 <- 2*f_hat_k*KM.if - KM.if
  phi_tilde_0 <- 2*f_hat_k*S_hat_k - f_hat_k^2 - S_hat_k - mean(2*f_hat_k*S_hat_k - f_hat_k^2 - S_hat_k)

  if_func <- phi0 + phi_tilde_0

  one_step <- mean(2*f_hat_k*S_hat_k - f_hat_k^2 - S_hat_k) + mean(if_func)
  plug_in <- mean(2*f_hat_k*S_hat_k - f_hat_k^2 - S_hat_k)

  return(list(one_step = one_step,
              plug_in = plug_in,
              if_func = if_func))
}
