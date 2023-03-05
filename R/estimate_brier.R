estimate_brier <- function(time,
                           event,
                           approx_times,
                           tau,
                           preds,
                           S_hat,
                           G_hat){

  KM_IFs <- calc_KM_IF(time = time,
                       event = event,
                       S_hat = S_hat,
                       G_hat = G_hat,
                       approx_times = approx_times)

  k <- min(which(approx_times >= tau))
  S_hat_k <- S_hat[,k]
  G_hat_k <- G_hat[,k]
  KM_IFs <- KM_IFs[,k]

  # calculate two components of EIF
  phi0 <- 2*preds*KM_IFs - KM_IFs
  phi_tilde_0 <- 2*preds*S_hat_k - preds^2 - S_hat_k - mean(2*preds*S_hat_k - preds^2 - S_hat_k)

  EIF <- phi0 + phi_tilde_0

  one_step <- mean(2*preds*S_hat_k - preds^2 - S_hat_k) + mean(EIF)
  plug_in <- mean(2*preds*S_hat_k - preds^2 - S_hat_k)

  return(list(one_step = one_step,
              plug_in = plug_in,
              EIF = EIF))
}
