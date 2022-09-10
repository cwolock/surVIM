vim_brier <- function(time,
                      event,
                      approx_times,
                      landmark_times,
                      f_hat,
                      fs_hat,
                      S_hat,
                      G_hat,
                      folds){

  one_step <- rep(NA, length(landmark_times))
  plug_in <- rep(NA, length(landmark_times))
  var_est <- rep(NA, length(landmark_times))
  full <- rep(NA, length(landmark_times))
  reduced <- rep(NA, length(landmark_times))
  for(i in 1:length(landmark_times)) {
    t0 <- landmark_times[i]
    CV_fulls <- rep(NA, length(unique(folds)))
    CV_reduceds <- rep(NA, length(unique(folds)))
    CV_one_steps <- rep(NA, length(unique(folds)))
    CV_plug_ins <- rep(NA, length(unique(folds)))
    CV_var_ests <- rep(NA, length(unique(folds)))
    for (j in 1:length(unique(folds))){
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

  return(data.frame(t = landmark_times,
                    full = full,
                    reduced = reduced,
                    one_step = one_step,
                    plug_in = plug_in,
                    var_est = var_est))
}
