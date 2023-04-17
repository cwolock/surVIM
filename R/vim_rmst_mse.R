vim_rmst_mse <- function(time,
                         event,
                         approx_times,
                         tau,
                         f_hat,
                         fs_hat,
                         S_hat,
                         G_hat,
                         folds,
                         sample_split,
                         ss_folds){

  V <- length(unique(folds))

  CV_full_plug_ins <- rep(NA, V)
  CV_reduced_plug_ins <- rep(NA, V)
  CV_full_one_steps <- rep(NA, V)
  CV_reduced_one_steps <- rep(NA, V)
  CV_one_steps <- rep(NA, V)
  CV_plug_ins <- rep(NA, V)
  CV_var_ests <- rep(NA, V)
  split_one_step_fulls <- rep(NA, V)
  split_plug_in_fulls <- rep(NA, V)
  split_one_step_reduceds <- rep(NA, V)
  split_plug_in_reduceds <- rep(NA, V)
  split_var_est_fulls <- rep(NA, V)
  split_var_est_reduceds <- rep(NA, V)
  for (j in 1:V){

    time_holdout <- time[folds == j]
    event_holdout <- event[folds == j]

    V_0 <- surVIM:::estimate_rmst_mse(time = time_holdout,
                                      event = event_holdout,
                                      approx_times = approx_times,
                                      tau = tau,
                                      preds = f_hat[[j]],
                                      S_hat = S_hat[[j]],
                                      G_hat = G_hat[[j]])
    V_0s <- surVIM:::estimate_rmst_mse(time = time_holdout,
                                       event = event_holdout,
                                       approx_times = approx_times,
                                       tau = tau,
                                       preds = fs_hat[[j]],
                                       S_hat = S_hat[[j]],
                                       G_hat = G_hat[[j]])

    CV_full_one_steps[j] <- V_0$one_step
    CV_full_plug_ins[j] <- V_0$plug_in
    CV_reduced_one_steps[j] <- V_0s$one_step
    CV_reduced_plug_ins[j] <- V_0s$plug_in
    CV_one_steps[j] <-  V_0$one_step -V_0s$one_step
    CV_plug_ins[j] <-  V_0$plug_in -V_0s$plug_in
    split_one_step_fulls[j] <- V_0$one_step
    split_plug_in_fulls[j] <- V_0$plug_in
    split_one_step_reduceds[j] <- V_0s$one_step
    split_plug_in_reduceds[j] <- V_0s$plug_in
    split_var_est_fulls[j] <- mean(V_0$EIF^2)
    split_var_est_reduceds[j] <- mean(V_0s$EIF^2)
    EIF <- V_0$EIF - V_0s$EIF
    CV_var_ests[j] <- mean(EIF^2)
  }

  if (sample_split){
    folds_0 <- sort(unique(folds[ss_folds == 0]))
    folds_1 <- sort(unique(folds[ss_folds == 1]))
    one_step <- mean(split_one_step_fulls[folds_0]) -
      mean(split_one_step_reduceds[folds_1])
    full_one_step <- mean(split_one_step_fulls[folds_0])
    full_plug_in <- mean(split_plug_in_fulls[folds_0])
    reduced_one_step <- mean(split_one_step_reduceds[folds_1])
    reduced_plug_in <- mean(split_plug_in_reduceds[folds_1])
    var_est <- mean(split_var_est_fulls[folds_0]) +
      mean(split_var_est_reduceds[folds_1])
  } else{
    one_step <- mean(CV_one_steps)
    var_est <- mean(CV_var_ests)
    full_one_step <- mean(CV_full_one_steps)
    reduced_one_step <- mean(CV_reduced_one_steps)
    full_plug_in <- mean(CV_full_plug_ins)
    reduced_plug_in <- mean(CV_reduced_plug_ins)
  }

  return(data.frame(tau = tau,
                    full_one_step = full_one_step,
                    reduced_one_step = reduced_one_step,
                    one_step = one_step,
                    full_plug_in = full_plug_in,
                    reduced_plug_in = reduced_plug_in,
                    var_est = var_est))
}
