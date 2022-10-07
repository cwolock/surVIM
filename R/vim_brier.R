vim_brier <- function(time,
                      event,
                      approx_times,
                      landmark_times,
                      f_hat,
                      fs_hat,
                      S_hat,
                      G_hat,
                      folds,
                      ss_folds,
                      sample_split){

  V <- length(unique(folds))
  # if (sample_split){
  #   V <- length(unique(folds))/2
  # } else{
  #   V <- length(unique(folds))
  # }

  # if (sample_split){
  #   full_folds <- sort(unique(folds[ss_folds == 0]))
  #   reduced_folds <- sort(unique(folds[ss_folds == 1]))
  # } else{
  #   full_folds <- NULL
  #   reduced_folds <- NULL
  # }
  # print(full_folds)
  # print(reduced_folds)

  one_step <- rep(NA, length(landmark_times))
  plug_in <- rep(NA, length(landmark_times))
  var_est <- rep(NA, length(landmark_times))
  full <- rep(NA, length(landmark_times))
  reduced <- rep(NA, length(landmark_times))
  for(i in 1:length(landmark_times)) {

    t0 <- landmark_times[i]

    CV_fulls <- rep(NA, V)
    CV_reduceds <- rep(NA, V)
    CV_one_steps <- rep(NA, V)
    CV_plug_ins <- rep(NA, V)
    CV_var_ests <- rep(NA, V)
    split_one_step_fulls <- rep(NA, V)
    split_plug_in_fulls <- rep(NA, V)
    split_one_step_reduceds <- rep(NA, V)
    split_plug_in_reduceds <- rep(NA, V)
    split_var_est_fulls <- rep(NA, V)
    split_var_est_reduceds <- rep(NA, V)

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



    # for (j in 1:V){
    #   print(full_folds[j])
    #   print(reduced_folds[j])
    #   if (sample_split){
    #     V_0 <- estimate_brier(time = time[folds == full_folds[j]],
    #                           event = time[folds == full_folds[j]],
    #                           approx_times = approx_times,
    #                           t0 = t0,
    #                           preds = f_hat[[full_folds[j]]][,i],
    #                           S_hat = S_hat[[full_folds[j]]],
    #                           G_hat = G_hat[[full_folds[j]]])
    #     V_0s <- estimate_brier(time = time[folds == reduced_folds[j]],
    #                            event = event[folds == reduced_folds[j]],
    #                            approx_times = approx_times,
    #                            t0 = t0,
    #                            preds = fs_hat[[reduced_folds[j]]][,i],
    #                            S_hat = S_hat[[reduced_folds[j]]],
    #                            G_hat = G_hat[[reduced_folds[j]]])
    #   } else{
    #     time_holdout <- time[folds == j]
    #     event_holdout <- event[folds == j]
    #     V_0 <- estimate_brier(time = time_holdout,
    #                           event = event_holdout,
    #                           approx_times = approx_times,
    #                           t0 = t0,
    #                           preds = f_hat[[j]][,i],
    #                           S_hat = S_hat[[j]],
    #                           G_hat = G_hat[[j]])
    #     V_0s <- estimate_brier(time = time_holdout,
    #                            event = event_holdout,
    #                            approx_times = approx_times,
    #                            t0 = t0,
    #                            preds = fs_hat[[j]][,i],
    #                            S_hat = S_hat[[j]],
    #                            G_hat = G_hat[[j]])
    #   }
    #
    #   CV_fulls[j] <- V_0$one_step
    #   CV_reduceds[j] <- V_0s$one_step
    #   CV_one_steps[j] <-  V_0$one_step -V_0s$one_step
    #   CV_plug_ins[j] <-  V_0$plug_in -V_0s$plug_in
    #   split_one_step_fulls[j] <- V_0$one_step
    #   split_plug_in_fulls[j] <- V_0$plug_in
    #   split_one_step_reduceds[j] <- V_0s$one_step
    #   split_plug_in_reduceds[j] <- V_0s$plug_in
    #   split_var_est_fulls[j] <- mean(V_0$if_func^2)
    #   split_var_est_reduceds[j] <- mean(V_0s$if_func^2)
    #   if_func <- V_0$if_func - V_0s$if_func
    #   CV_var_ests[j] <- mean(if_func^2)
    # }
    #
    # if (sample_split){
    #   one_step[i] <- mean(split_one_step_fulls) - mean(split_one_step_reduceds)
    #   plug_in[i] <- mean(split_plug_in_fulls) - mean(split_plug_in_reduceds)
    #   full[i] <- mean(split_one_step_fulls)
    #   reduced[i] <- mean(split_one_step_reduceds)
    #   var_est[i] <- mean(split_var_est_fulls) + mean(split_var_est_reduceds)
    # } else{
    #   one_step[i] <- mean(CV_one_steps)
    #   plug_in[i] <- mean(CV_plug_ins)
    #   var_est[i] <- mean(CV_var_ests)
    #   full[i] <- mean(CV_fulls)#V_0$one_step
    #   reduced[i] <- mean(CV_reduceds)#V_0s$one_step
    # }

    # if (sample_split){
    #   one_step[i] <- mean(CV_brier_full$split_one_step_fulls) - mean(CV_brier_reduced$split_one_step_reduceds)
    #   plug_in[i] <-mean(CV_brier_full$split_plug_in_fulls) - mean(CV_brier_reduced$split_plug_in_reduceds)
    #   full[i] <- mean(CV_brier_full$split_one_step_fulls)
    #   reduced[i] <- mean(CV_brier_reduced$split_one_step_reduceds)
    #   var_est[i] <- mean(CV_brier_full$split_var_est_fulls) + mean(CV_brier_reduced$split_var_est_reduceds)
    # } else{
    #   one_step[i] <- mean(CV_brier_full$CV_one_steps)
    #   plug_in[i] <- mean(CV_brier_full$CV_plug_ins)
    #   var_est[i] <- mean(CV_brier_full$CV_var_ests)
    #   full[i] <- mean(CV_brier_full$CV_fulls)#V_0$one_step
    #   reduced[i] <- mean(CV_brier_full$CV_reduceds)#V_0s$one_step
    # }




    # CV_brier_full <- CV_brier(time = time[ss_folds == 0],
    #                           event= event[ss_folds == 0],
    #                           approx_times = approx_times,
    #                           t0 = t0,
    #                           f_hat= f_hat[ss_folds == 0],
    #                           fs_hat= fs_hat[ss_folds == 0],
    #                           S_hat = S_hat[folds%in% which(ss_folds == 1)],
    #                           G_hat= G_hat[folds%in% which(ss_folds == 1)],
    #                           folds = folds[folds%in% which(ss_folds == 1)],
    #                           landmark_times = landmark_times)
    #
    # if (sample_split){
    #   CV_brier_reduced <- CV_brier(time = time[folds%in% which(ss_folds == 2)],
    #                                event= event[folds%in% which(ss_folds == 2)],
    #                                approx_times = approx_times,
    #                                t0 = t0,
    #                                f_hat= f_hat[folds%in% which(ss_folds == 2)],
    #                                fs_hat= fs_hat[folds%in% which(ss_folds == 2)],
    #                                S_hat = S_hat[folds%in% which(ss_folds == 2)],
    #                                G_hat= G_hat[folds%in% which(ss_folds == 2)],
    #                                folds = folds[folds%in% which(ss_folds == 2)],
    #                                landmark_times = landmark_times)
    #   one_step[i] <- mean(CV_brier_full$split_one_step_fulls) - mean(CV_brier_reduced$split_one_step_reduceds)
    #   plug_in[i] <-mean(CV_brier_full$split_plug_in_fulls) - mean(CV_brier_reduced$split_plug_in_reduceds)
    #   full[i] <- mean(CV_brier_full$split_one_step_fulls)
    #   reduced[i] <- mean(CV_brier_reduced$split_one_step_reduceds)
    #   var_est[i] <- mean(CV_brier_full$split_var_est_fulls) + mean(CV_brier_reduced$split_var_est_reduceds)
    # } else{
    #   one_step[i] <- mean(CV_brier_full$CV_one_steps)
    #   plug_in[i] <- mean(CV_brier_full$CV_plug_ins)
    #   var_est[i] <- mean(CV_brier_full$CV_var_ests)
    #   full[i] <- mean(CV_brier_full$CV_fulls)#V_0$one_step
    #   reduced[i] <- mean(CV_brier_full$CV_reduceds)#V_0s$one_step
    # }
  }

  return(data.frame(t = landmark_times,
                    full = full,
                    reduced = reduced,
                    one_step = one_step,
                    plug_in = plug_in,
                    var_est = var_est))
}
