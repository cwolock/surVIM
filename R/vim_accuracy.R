vim_accuracy <- function(time,
                         event,
                         approx_times,
                         landmark_times,
                         f_hat,
                         fs_hat,
                         S_hat,
                         G_hat,
                         folds,
                         sample_split,
                         ss_folds,
                         bootstrap = FALSE,
                         B = 1000){
  n <- length(time)
  one_step <- rep(NA, length(landmark_times))
  plug_in <- rep(NA, length(landmark_times))
  var_est <- rep(NA, length(landmark_times))
  full <- rep(NA, length(landmark_times))
  reduced <- rep(NA, length(landmark_times))
  boot_var_est <- rep(NA, length(landmark_times))
  for(i in 1:length(landmark_times)) {
    t0 <- landmark_times[i]
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
    CV_boot_var_ests<- rep(NA, V)
    for (j in 1:length(unique(folds))){
      time_holdout <- time[folds == j]
      event_holdout <- event[folds == j]
      V_0 <- estimate_accuracy(time = time_holdout,
                               event = event_holdout,
                               approx_times = approx_times,
                               t0 = t0,
                               preds = f_hat[[j]][,i],
                               S_hat = S_hat[[j]],
                               G_hat = G_hat[[j]])
      V_0s <- estimate_accuracy(time = time_holdout,
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

      if (bootstrap){
        boot_CV_ests <- rep(NA, B)
        for (k in 1:B){
          boot_inds <- sample(1:n, size = n, replace = TRUE)
          V_0 <- estimate_brier(time = time_holdout[boot_inds],
                                event = event_holdout[boot_inds],
                                approx_times = approx_times,
                                t0 = t0,
                                preds = f_hat[[j]][boot_inds,i],
                                S_hat = S_hat[[j]][boot_inds,],
                                G_hat = G_hat[[j]][boot_inds,])
          V_0s <- estimate_brier(time = time_holdout[boot_inds],
                                 event = event_holdout[boot_inds],
                                 approx_times = approx_times,
                                 t0 = t0,
                                 preds = fs_hat[[j]][boot_inds,i],
                                 S_hat = S_hat[[j]][boot_inds,],
                                 G_hat = G_hat[[j]][boot_inds,])
          boot_CV_ests[k] <- V_0$plug_in -V_0s$plug_in
        }
        CV_boot_var_ests[j] <- var(boot_CV_ests)
      }

    }

    if (sample_split){
      one_step[i] <- mean(split_one_step_fulls[sort(unique(folds[ss_folds == 0]))]) -
        mean(split_one_step_reduceds[sort(unique(folds[ss_folds == 1]))])
      plug_in[i] <- mean(split_plug_in_fulls[sort(unique(folds[ss_folds == 0]))]) -
        mean(split_plug_in_reduceds[sort(unique(folds[ss_folds == 1]))])
      full[i] <- mean(split_one_step_fulls[sort(unique(folds[ss_folds == 0]))])
      reduced[i] <- mean(split_one_step_reduceds[sort(unique(folds[ss_folds == 1]))])
      var_est[i] <- mean(split_var_est_fulls[sort(unique(folds[ss_folds == 0]))]) +
        mean(split_var_est_reduceds[sort(unique(folds[ss_folds == 1]))])
    } else{
      one_step[i] <- mean(CV_one_steps)
      plug_in[i] <- mean(CV_plug_ins)
      var_est[i] <- mean(CV_var_ests)
      full[i] <- mean(CV_fulls)#V_0$one_step
      reduced[i] <- mean(CV_reduceds)#V_0s$one_step
    }


    if (bootstrap){
      boot_var_est[i] <- mean(CV_boot_var_ests)
    }
  }

  return(data.frame(t = landmark_times,
                    full = full,
                    reduced = reduced,
                    one_step = one_step,
                    plug_in = plug_in,
                    var_est = var_est,
                    boot_var_est = boot_var_est))
}
