vim_AUC_new <- function(time,
                    event,
                    approx_times,
                    landmark_times,
                    f_hat,
                    fs_hat,
                    S_hat,
                    G_hat,
                    folds,
                    sample_split,
                    ss_folds){
  n <- length(time)
  one_step <- rep(NA, length(landmark_times))
  plug_in <- rep(NA, length(landmark_times))
  var_est <- rep(NA, length(landmark_times))
  full_one_step <- rep(NA, length(landmark_times))
  reduced_one_step <- rep(NA, length(landmark_times))
  full_plug_in <- rep(NA, length(landmark_times))
  reduced_plug_in <- rep(NA, length(landmark_times))
  for(i in 1:length(landmark_times)) {
    t0 <- landmark_times[i]
    CV_full_plug_ins <- rep(NA, length(unique(folds)))
    CV_reduced_plug_ins <- rep(NA, length(unique(folds)))
    CV_full_one_steps <- rep(NA, length(unique(folds)))
    CV_reduced_one_steps <- rep(NA, length(unique(folds)))
    CV_one_steps <- rep(NA, length(unique(folds)))
    CV_plug_ins <- rep(NA, length(unique(folds)))
    CV_var_ests <- rep(NA, length(unique(folds)))
    CV_numerators <- rep(NA, length(unique(folds)))
    CV_denominators <- rep(NA, length(unique(folds)))
    split_numerator_fulls <- rep(NA, length(unique(folds)))
    split_denominator_fulls <- rep(NA, length(unique(folds)))
    split_numerator_reduceds <- rep(NA, length(unique(folds)))
    split_denominator_reduceds <- rep(NA, length(unique(folds)))
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
        if (sample_split & length(unique(folds)) > 2){ # xfit and sample split
          other_js <- unique(folds[ss_folds == unique(ss_folds[folds == j]) & folds != j])
          # print(j)
          # print(other_js)
          #curr_folds <- which(subfolds == unique(ss_folds[folds == j]))
          # print(ss_folds)
          # print(subfolds)
          # print(curr_folds)
          preds_holdout <- unlist(lapply(f_hat, function(x) x[,i])[other_js])
          S_hat_holdout = do.call(rbind, S_hat[other_js])
          G_hat_holdout = do.call(rbind, G_hat[other_js])
          preds_holdout_reduced = unlist(lapply(fs_hat, function(x) x[,i])[other_js])
          # preds_holdout2 <- unlist(lapply(f_hat, function(x) x[,i])[-j])[curr_folds]
          # S_hat_holdout2 = do.call(rbind, S_hat[-j])[curr_folds,]
          # preds_holdout_reduced2 = unlist(lapply(fs_hat, function(x) x[,i])[-j])[curr_folds]
        } else if (sample_split & length(unique(folds)) <= 2){ # sample split only, no xfit
          preds_holdout <- f_hat[[j]][,i]
          S_hat_holdout = S_hat[[j]]
          G_hat_holdout = G_hat[[j]]
          preds_holdout_reduced = fs_hat[[j]][,i]
        } else{ # xfit, no sample split
          preds_holdout <- unlist(lapply(f_hat, function(x) x[,i])[-j])
          S_hat_holdout = do.call(rbind, S_hat[-j])
          G_hat_holdout = do.call(rbind, G_hat[-j])
          preds_holdout_reduced = unlist(lapply(fs_hat, function(x) x[,i])[-j])
          time_notholdout = time[folds != j]
          event_notholdout = event[folds != j]
        }
      } else{ # no xfit, no sample split
        preds_holdout <- f_hat[[j]][,i]
        S_hat_holdout = S_hat[[j]]
        G_hat_holdout = G_hat[[j]]
        preds_holdout_reduced = fs_hat[[j]][,i]
        time_notholdout = time[folds == j]
        event_notholdout = event[folds == j]
      }
      # print(length(preds_holdout))
      # print(dim(S_hat_holdout))
      # print(length(preds_holdout_reduced))
      # fix bug here in how the holdout data is picked with sample splitting
      V_0 <- estimate_AUC_new(time = time_holdout,
                          event = event_holdout,
                          approx_times = approx_times,
                          t0 = t0,
                          preds = f_hat[[j]][,i],
                          preds_holdout = preds_holdout,
                          S_hat = S_hat[[j]],
                          S_hat_holdout = S_hat_holdout,
                          G_hat = G_hat[[j]],
                          G_hat_holdout = G_hat_holdout,
                          time_holdout = time_notholdout,
                          event_holdout = event_notholdout)
      V_0s <- estimate_AUC_new(time = time_holdout,
                           event = event_holdout,
                           approx_times = approx_times,
                           t0 = t0,
                           preds = fs_hat[[j]][,i],
                           preds_holdout = preds_holdout_reduced,
                           S_hat = S_hat[[j]],
                           S_hat_holdout = S_hat_holdout,
                           G_hat = G_hat[[j]],
                           G_hat_holdout = G_hat_holdout,
                           time_holdout = time_notholdout,
                           event_holdout = event_notholdout)
      CV_full_one_steps[j] <- V_0$one_step
      CV_full_plug_ins[j] <- V_0$plug_in
      CV_reduced_one_steps[j] <- V_0s$one_step
      CV_reduced_plug_ins[j] <- V_0s$plug_in
      CV_one_steps[j] <-  V_0$one_step -V_0s$one_step
      CV_plug_ins[j] <-  V_0$plug_in -V_0s$plug_in
      CV_numerators[j] <- V_0$numerator - V_0s$numerator
      CV_denominators[j] <- V_0$denominator
      split_numerator_fulls[j] <- V_0$numerator
      split_denominator_fulls[j] <- V_0$denominator
      split_numerator_reduceds[j] <- V_0s$numerator
      split_denominator_reduceds[j] <- V_0s$denominator
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
      one_step[i] <- mean(split_numerator_fulls[sort(unique(folds[ss_folds == 0]))])/mean(split_denominator_fulls[sort(unique(folds[ss_folds == 0]))]) -
        mean(split_numerator_reduceds[sort(unique(folds[ss_folds == 1]))])/mean(split_denominator_reduceds[sort(unique(folds[ss_folds == 1]))])
      # one_step[i] <- mean(split_one_step_fulls[sort(unique(folds[ss_folds == 0]))]) -
      #   mean(split_one_step_reduceds[sort(unique(folds[ss_folds == 1]))])
      plug_in[i] <- mean(split_plug_in_fulls[sort(unique(folds[ss_folds == 0]))]) -
        mean(split_plug_in_reduceds[sort(unique(folds[ss_folds == 1]))])
      full_one_step[i] <- mean(split_one_step_fulls[sort(unique(folds[ss_folds == 0]))])
      full_plug_in[i] <- mean(split_plug_in_fulls[sort(unique(folds[ss_folds == 0]))])
      reduced_one_step[i] <- mean(split_one_step_reduceds[sort(unique(folds[ss_folds == 1]))])
      reduced_plug_in[i] <- mean(split_plug_in_reduceds[sort(unique(folds[ss_folds == 1]))])
      var_est[i] <- mean(split_var_est_fulls[sort(unique(folds[ss_folds == 0]))]) +
        mean(split_var_est_reduceds[sort(unique(folds[ss_folds == 1]))])
    } else{
      one_step[i] <- mean(CV_numerators)/mean(CV_denominators)#mean(CV_one_steps)
      plug_in[i] <- mean(CV_plug_ins)
      var_est[i] <- mean(CV_var_ests)
      full_one_step[i] <- mean(CV_full_one_steps)#V_0$one_step
      reduced_one_step[i] <- mean(CV_reduced_one_steps)#V_0s$one_step
      full_plug_in[i] <- mean(CV_full_plug_ins)#V_0$one_step
      reduced_plug_in[i] <- mean(CV_reduced_plug_ins)#V_0s$one_step
    }

  }

  return(data.frame(t = landmark_times,
                    full_one_step = full_one_step,
                    reduced_one_step = reduced_one_step,
                    one_step = one_step,
                    plug_in = plug_in,
                    full_plug_in = full_plug_in,
                    reduced_plug_in = reduced_plug_in,
                    var_est = var_est))
}
