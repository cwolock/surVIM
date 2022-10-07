CV_brier <- function(time,
                     event,
                     approx_times = approx_times,
                     t0 = t0,
                     f_hat,
                     fs_hat,
                     S_hat ,
                     G_hat,
                     folds,
                     landmark_times){
  i <- which(landmark_times == t0)
  V <- length(unique(folds))
  folds <- as.numeric(factor(rank(folds)))
  CV_fulls <- rep(NA, length(unique(folds)))
  CV_reduceds <- rep(NA, length(unique(folds)))
  CV_one_steps <- rep(NA, length(unique(folds)))
  CV_plug_ins <- rep(NA, length(unique(folds)))
  CV_var_ests <- rep(NA, length(unique(folds)))
  split_one_step_fulls <- rep(NA, length(unique(folds)))
  split_plug_in_fulls <- rep(NA, length(unique(folds)))
  split_one_step_reduceds <- rep(NA, length(unique(folds)))
  split_plug_in_reduceds <- rep(NA, length(unique(folds)))
  split_var_est_fulls <- rep(NA, length(unique(folds)))
  split_var_est_reduceds <- rep(NA, length(unique(folds)))
  for (j in 1:V){
    time_holdout <- time[folds == j]
    event_holdout <- event[folds == j]
    V_0 <- estimate_brier(time = time_holdout,
                          event = event_holdout,
                          approx_times = approx_times,
                          t0 = t0,
                          preds = f_hat[[j]][,i],
                          S_hat = S_hat[[j]],
                          G_hat = G_hat[[j]])
    V_0s <- estimate_brier(time = time_holdout,
                           event = event_holdout,
                           approx_times = approx_times,
                           t0 = t0,
                           preds = fs_hat[[j]][,i],
                           S_hat = S_hat[[j]],
                           G_hat = G_hat[[j]])
    CV_fulls[j] <- V_0$one_step
    CV_reduceds[j] <- V_0s$one_step
    CV_one_steps[j] <-  V_0$one_step -V_0s$one_step
    CV_plug_ins[j] <-  V_0$plug_in -V_0s$plug_in
    split_one_step_fulls[j] <- V_0$one_step
    split_plug_in_fulls[j] <- V_0$plug_in
    split_one_step_reduceds[j] <- V_0s$one_step
    split_plug_in_reduceds[j] <- V_0s$plug_in
    split_var_est_fulls[j] <- mean(V_0$if_func^2)
    split_var_est_reduceds[j] <- mean(V_0s$if_func^2)
    if_func <- V_0$if_func - V_0s$if_func
    CV_var_ests[j] <- mean(if_func^2)
  }

  return(list(CV_fulls = CV_fulls,
              CV_reduceds = CV_reduceds,
              CV_one_steps = CV_one_steps,
              CV_plug_ins = CV_plug_ins,
              split_one_step_fulls =split_one_step_fulls,
              split_plug_in_fulls =split_plug_in_fulls,
              split_one_step_reduceds =split_one_step_reduceds,
              split_plug_in_reduceds =split_plug_in_reduceds,
              split_var_est_fulls =split_var_est_fulls,
              split_var_est_reduceds =split_var_est_reduceds,
              CV_var_ests = CV_var_ests))
}
