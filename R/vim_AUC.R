vim_AUC <- function(time,
                    event,
                    approx_times,
                    landmark_times,
                    f_hat,
                    fs_hat,
                    S_hat,
                    G_hat,
                    folds){

  # time <- holdout$y
  # event <- holdout$delta

  one_step <- rep(NA, length(landmark_times))
  plug_in <- rep(NA, length(landmark_times))
  var_est <- rep(NA, length(landmark_times))
  full <- rep(NA, length(landmark_times))
  reduced <- rep(NA, length(landmark_times))
  for(i in 1:length(landmark_times)) {
    t0 <- landmark_times[i]
    V_0 <- estimate_AUC(time = time,
                        event = event,
                        approx_times = approx_times,
                        t0 = t0,
                        preds = f_hat[,i],
                        S_hat = S_hat,
                        G_hat = G_hat)
    V_0s <- estimate_AUC(time = time,
                         event = event,
                         approx_times = approx_times,
                         t0 = t0,
                         preds = fs_hat[,i],
                         S_hat = S_hat,
                         G_hat = G_hat)
    full[i] <- V_0$one_step
    reduced[i] <- V_0s$one_step
    one_step[i] <- V_0$one_step - V_0s$one_step
    plug_in[i] <- V_0$plug_in - V_0s$plug_in
    if_func <- V_0$if_func - V_0s$if_func
    var_est[i] <- mean(if_func^2)
  }

  return(data.frame(t = landmark_times,
                    full = full,
                    reduced = reduced,
                    one_step = one_step,
                    plug_in = plug_in,
                    var_est = var_est))
}
