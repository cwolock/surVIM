vim_cindex <- function(time,
                       event,
                       approx_times,
                       tau,
                       f_hat,
                       fs_hat,
                       S_hat,
                       G_hat,
                       folds,
                       sample_split){

  # CV_fulls <- rep(NA, length(unique(folds)))
  # CV_reduceds <- rep(NA, length(unique(folds)))
  # CV_one_steps <- rep(NA, length(unique(folds)))
  # CV_plug_ins <- rep(NA, length(unique(folds)))
  # CV_var_ests <- rep(NA, length(unique(folds)))
  # split_one_step_fulls <- rep(NA, length(unique(folds)))
  # split_plug_in_fulls <- rep(NA, length(unique(folds)))
  # split_one_step_reduceds <- rep(NA, length(unique(folds)))
  # split_plug_in_reduceds <- rep(NA, length(unique(folds)))
  # split_var_est_fulls <- rep(NA, length(unique(folds)))
  # split_var_est_reduceds <- rep(NA, length(unique(folds)))
  # for (j in 1:length(unique(folds))){
  #   time_holdout <- time[folds == j]
  #   event_holdout <- event[folds == j]
  #   if (length(unique(folds)) > 1){
  #     subfolds <- ss_folds[folds != j]
  #     if (sample_split){
  #       curr_folds <- which(subfolds == unique(ss_folds[folds == j]))
  #       preds_holdout <- unlist(lapply(f_hat, function(x) x[,i])[-j])[curr_folds]
  #       S_hat_holdout = do.call(rbind, S_hat[-j])[curr_folds,]
  #       preds_holdout_reduced = unlist(lapply(fs_hat, function(x) x[,i])[-j])[curr_folds]
  #     } else{
  #       preds_holdout <- unlist(lapply(f_hat, function(x) x[,i])[-j])
  #       S_hat_holdout = do.call(rbind, S_hat[-j])
  #       preds_holdout_reduced = unlist(lapply(fs_hat, function(x) x[,i])[-j])
  #     }
  #
  #   } else{
  #     preds_holdout <- f_hat[[j]][,i]
  #     S_hat_holdout = S_hat[[j]]
  #     preds_holdout_reduced = fs_hat[[j]][,i]
  #   }
  #   V_0 <- estimate_cindex_newxfit_newapprox(time = time_holdout,
  #                          event = event_holdout,
  #                          approx_times = approx_times,
  #                          tau = tau,
  #                          preds = f_hat[[j]][,i],
  #                          preds_holdout = preds_holdout,
  #                          S_hat = S_hat[[j]],
  #                          S_hat_holdout = S_hat_holdout,
  #                          G_hat = G_hat[[j]])
  #   V_0s <- estimate_cindex_newxfit_newapprox(time = time_holdout,
  #                           event = event_holdout,
  #                           approx_times = approx_times,
  #                           tau = tau,
  #                           preds = fs_hat[[j]][,i],
  #                           preds_holdout = preds_holdout_reduced,
  #                           S_hat = S_hat[[j]],
  #                           S_hat_holdout = S_hat_holdout,
  #                           G_hat = G_hat[[j]])
  #   CV_fulls[j] <- V_0$one_step
  #   CV_reduceds[j] <- V_0s$one_step
  #   CV_one_steps[j] <-  V_0$one_step -V_0s$one_step
  #   CV_plug_ins[j] <-  V_0$plug_in -V_0s$plug_in
  #   # split_one_step_fulls[j] <- V_0$one_step
  #   # split_plug_in_fulls[j] <- V_0$plug_in
  #   # split_one_step_reduceds[j] <- V_0s$one_step
  #   # split_plug_in_reduceds[j] <- V_0s$plug_in
  #   # split_var_est_fulls[j] <- mean(V_0$if_func^2)
  #   # split_var_est_reduceds[j] <- mean(V_0s$if_func^2)
  #   if_func <- V_0$if_func - V_0s$if_func
  #   CV_var_ests[j] <- mean(if_func^2)
  # }
  #
  # if (sample_split){
  #   ss_folds <- sample(rep(1:2, length = length(unique(folds))))
  #   one_step <- mean(split_one_step_fulls[ss_folds == 1]) - mean(split_one_step_reduceds[ss_folds == 2])
  #   plug_in <- mean(split_plug_in_fulls[ss_folds == 1]) - mean(split_plug_in_reduceds[ss_folds == 2])
  #   full <- mean(split_one_step_fulls[ss_folds == 1])
  #   reduced <- mean(split_one_step_reduceds[ss_folds == 2])
  #   var_est <- mean(split_var_est_fulls[ss_folds == 1]) + mean(split_var_est_reduceds[ss_folds == 2])
  # } else{
  #   one_step <- mean(CV_one_steps)
  #   plug_in <- mean(CV_plug_ins)
  #   var_est <- mean(CV_var_ests)
  #   full <- mean(CV_fulls)#V_0$one_step
  #   reduced <- mean(CV_reduceds)#V_0s$one_step
  # }
  #
  #
  # return(data.frame(tau = tau,
  #                   full = full,
  #                   reduced = reduced,
  #                   one_step = one_step,
  #                   plug_in = plug_in,
  #                   var_est = var_est))

  # one_step <- rep(NA, length(landmark_times))
  # plug_in <- rep(NA, length(landmark_times))
  # var_est <- rep(NA, length(landmark_times))
  # full <- rep(NA, length(landmark_times))
  # reduced <- rep(NA, length(landmark_times))
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
  for (j in 1:length(unique(folds))){
    time_holdout <- time[folds == j]
    event_holdout <- event[folds == j]
    if (length(unique(folds)) > 1){
      subfolds <- ss_folds[folds != j]
      if (sample_split){
        curr_folds <- which(subfolds == unique(ss_folds[folds == j]))
        preds_holdout <- unlist(f_hat[-j])[curr_folds]
        S_hat_holdout = do.call(rbind, S_hat[-j])[curr_folds,]
        preds_holdout_reduced = unlist(fs_hat[-j])[curr_folds]
      } else{
        preds_holdout <- unlist(f_hat[-j])
        S_hat_holdout = do.call(rbind, S_hat[-j])
        preds_holdout_reduced = unlist(fs_hat[-j])
      }

    } else{
      preds_holdout <- f_hat[[j]]
      S_hat_holdout = S_hat[[j]]
      preds_holdout_reduced = fs_hat[[j]]
    }
    # fix bug here in how the holdout data is picked with sample splitting
    V_0 <- estimate_cindex(time = time_holdout,
                           event = event_holdout,
                           approx_times = approx_times,
                           tau = tau,
                           preds = f_hat[[j]][,i],
                           preds_holdout = preds_holdout,
                           S_hat = S_hat[[j]],
                           S_hat_holdout = S_hat_holdout,
                           G_hat = G_hat[[j]])
    V_0s <- estimate_cindex(time = time_holdout,
                            event = event_holdout,
                            approx_times = approx_times,
                            tau = tau,
                            preds = fs_hat[[j]][,i],
                            preds_holdout = preds_holdout_reduced,
                            S_hat = S_hat[[j]],
                            S_hat_holdout = S_hat_holdout,
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


  if (sample_split){
    one_step <- mean(split_one_step_fulls[sort(unique(folds[ss_folds == 0]))]) -
      mean(split_one_step_reduceds[sort(unique(folds[ss_folds == 1]))])
    plug_in <- mean(split_plug_in_fulls[sort(unique(folds[ss_folds == 0]))]) -
      mean(split_plug_in_reduceds[sort(unique(folds[ss_folds == 1]))])
    full <- mean(split_one_step_fulls[sort(unique(folds[ss_folds == 0]))])
    reduced <- mean(split_one_step_reduceds[sort(unique(folds[ss_folds == 1]))])
    var_est <- mean(split_var_est_fulls[sort(unique(folds[ss_folds == 0]))]) +
      mean(split_var_est_reduceds[sort(unique(folds[ss_folds == 1]))])
  } else{
    one_step <- mean(CV_one_steps)
    plug_in <- mean(CV_plug_ins)
    var_est <- mean(CV_var_ests)
    full <- mean(CV_fulls)#V_0$one_step
    reduced <- mean(CV_reduceds)#V_0s$one_step
  }


  return(data.frame(tau = tau,
                    full = full,
                    reduced = reduced,
                    one_step = one_step,
                    plug_in = plug_in,
                    var_est = var_est))
}
