estimate_rsquared <- function(time,
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
  phi_01 <- 2*preds*KM_IFs - KM_IFs
  phi_tilde_01_uncentered <- 2*preds*S_hat_k - preds^2 - S_hat_k
  phi_tilde_01 <- phi_tilde_01_uncentered - mean(phi_tilde_01_uncentered)

  phi_01_combined <- 2*preds*(S_hat_k + KM_IFs) - preds^2 - (S_hat_k + KM_IFs)
  phi_02_combined <- (S_hat_k + KM_IFs)*(1 - S_hat_k - KM_IFs)

  phi_02 <- KM_IFs - 2*KM_IFs
  phi_tilde_02_uncentered <- S_hat_k - S_hat_k^2
  phi_tilde_02 <- phi_tilde_02_uncentered - mean(phi_tilde_02_uncentered)

  V_1 <- mean(phi_tilde_01_uncentered)
  V_2 <- mean(phi_tilde_02_uncentered)

  EIF <- (phi_01 + phi_tilde_01)/V_2 - V_1/(V_2^2)*(phi_02 + phi_tilde_02)

  V_1_alternative <- mean(phi_01_combined)
  V_2_alternative <- mean(phi_02_combined)

  one_step <- V_1_alternative/V_2_alternative
  plug_in <- V_1/V_2



  return(list(one_step = one_step,
              plug_in = plug_in,
              EIF = EIF,
              numerator = V_1_alternative,
              denominator = V_2_alternative))
}
