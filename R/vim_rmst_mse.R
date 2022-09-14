vim_rmst_mse <- function(time,
                     event,
                     approx_times,
                     taus,
                     f_hat,
                     fs_hat,
                     S_hat,
                     G_hat,
                     folds){

  one_step <- rep(NA, length(taus))
  plug_in <- rep(NA, length(taus))
  var_est <- rep(NA, length(taus))
  full <- rep(NA, length(taus))
  reduced <- rep(NA, length(taus))
  for(i in 1:length(taus)) {
    tau <- taus[i]
    CV_fulls <- rep(NA, length(unique(folds)))
    CV_reduceds <- rep(NA, length(unique(folds)))
    CV_one_steps <- rep(NA, length(unique(folds)))
    CV_plug_ins <- rep(NA, length(unique(folds)))
    CV_var_ests <- rep(NA, length(unique(folds)))
    for (j in 1:length(unique(folds))){
      time_holdout <- time[folds == j]
      event_holdout <- event[folds == j]
      V_0 <- estimate_rmst_mse(time = time_holdout,
                               event = event_holdout,
                               approx_times = approx_times,
                               tau = tau,
                               preds = f_hat[[j]][,i],
                               S_hat = S_hat[[j]],
                               G_hat = G_hat[[j]])
      V_0s <- estimate_rmst_mse(time = time_holdout,
                                event = event_holdout,
                                approx_times = approx_times,
                                tau = tau,
                                preds = fs_hat[[j]][,i],
                                S_hat = S_hat[[j]],
                                G_hat = G_hat[[j]])
      CV_fulls[j] <- V_0$plug_in
      CV_reduceds[j] <- V_0s$plug_in
      CV_one_steps[j] <-  V_0$one_step -V_0s$one_step
      CV_plug_ins[j] <-  V_0$plug_in -V_0s$plug_in
      if_func <- V_0$if_func - V_0s$if_func
      CV_var_ests[j] <- mean(if_func^2)
    }


    full[i] <- mean(CV_fulls)#V_0$one_step
    reduced[i] <- mean(CV_reduceds)#V_0s$one_step
    one_step[i] <- mean(CV_one_steps)#V_0$one_step - V_0s$one_step
    plug_in[i] <- mean(CV_plug_ins)#V_0$plug_in - V_0s$plug_in
    #if_func <- V_0$if_func - V_0s$if_func
    var_est[i] <- mean(CV_var_ests)#mean(if_func^2)
  }

  return(data.frame(t = taus,
                    full = full,
                    reduced = reduced,
                    one_step = one_step,
                    plug_in = plug_in,
                    var_est = var_est))
}
