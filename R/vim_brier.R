vim_brier <- function(time,
                      event,
                      X,
                      X_reduced,
                      approx_times,
                      landmark_times,
                      f_hat,
                      fs_hat,
                      S_hat,
                      G_hat,
                      holdout){

  # NOTE: Ted's function breaks with just a single covariate
  # also there are weird namespace issues with predict method
  # also why is new.times a mandatory argument
  dimension <- 4
  time_holdout <- holdout$y
  event_holdout <- holdout$delta
  X_holdout <-holdout[,1:dimension]
  X_reduced_holdout <- holdout[,2:dimension]


  time <- time_holdout
  event <- event_holdout
  n <- length(time)
  # calculate integral for KM influence function
  # check left vs right continuous stuff here
  int.vals <- t(sapply(1:n, function(i) {
    vals <- diff(1/S_hat[i,])* 1/ G_hat[i,-ncol(G_hat)]
    if(any(approx_times[-1] > time[i])){
      vals[approx_times[-1] > time[i]] <- 0
    }
    c(0,cumsum(vals))
  }))
  # again, check left v right continuous here, especially with right = TRUE/FALSE
  S_hat_Y <- sapply(1:n, function(i) stepfun(approx_times, c(1,S_hat[i,]), right = FALSE)(time[i]))
  G_hat_Y <- sapply(1:n, function(i) stepfun(approx_times, c(1,G_hat[i,]), right = TRUE)(time[i]))
  # places to hold the value of the influence function, as well as the actual estimate
  one_step <- rep(NA, length(landmark_times))
  plug_in <- rep(NA, length(landmark_times))
  var_est <- rep(NA, length(landmark_times))

  for(i in 1:length(landmark_times)) {
    t0 <- landmark_times[i]
    k <- min(which(approx_times >= t0))
    S_hat_k <- S_hat[,k]
    G_hat_k <- G_hat[,k]
    f_hat_k <- f_hat[,i]
    #fs_hat_k <- fs_hat[,i]
    inner.func.1 <- ifelse(time <= t0 & event == 1, 1/(S_hat_Y * G_hat_Y), 0 )
    inner.func.2 <- int.vals[,k]
    KM.if <- -S_hat_k * ( inner.func.1 - inner.func.2)

    # should phi_0 and phi_0s be negative or positive...
    # phi0_old <- 2 * KM.if * (f_hat_k - S_hat_k)
    # phi_tilde_0_old <- -(f_hat_k - S_hat_k)^2 - mean(-(f_hat_k - S_hat_k)^2)
    #
    # phi0s_old <- 2 * KM.if * (fs_hat_k - S_hat_k)
    # phi_tilde_0s_old <- -(fs_hat_k - S_hat_k)^2 - mean(-(fs_hat_k - S_hat_k)^2)

    phi0 <- 2*f_hat_k*KM.if - KM.if
    phi_tilde_0 <- 2*f_hat_k*S_hat_k - f_hat_k^2 - S_hat_k - mean(2*f_hat_k*S_hat_k - f_hat_k^2 - S_hat_k)
    # phi0s <- 2*fs_hat_k*KM.if - KM.if
    # phi_tilde_0s <- 2*fs_hat_k*S_hat_k - fs_hat_k^2 - S_hat_k

    # eta0 <- -mean(S_hat_k*(1 - S_hat_k))
    # eta0_if_tilde <- -S_hat_k*(1 - S_hat_k) - mean(-S_hat_k*(1 - S_hat_k))
    # eta0_if <- -KM.if + 2*KM.if*S_hat_k

    if.func <- phi0 + phi_tilde_0# - phi0s - phi_tilde_0s
    # if.func_old <- phi0_old + phi_tilde_0_old - phi0s_old - phi_tilde_0s_old
    # if.func_eta <- eta0_if_tilde + eta0_if
    # mu <- mean(if.func)
    # sigma <- sd(if.func)
    # if.func.z <- (if.func-mu)/sigma
    # if.func <- if.func[if.func.z > -4 & if.func.z < 4]
    # print(sum(if.func.z < -5 | if.func.z > 5)/length(if.func.z)*100)

    one_step[i] <- mean(2*f_hat_k*S_hat_k - f_hat_k^2 - S_hat_k) + mean(if.func)
    plug_in[i] <- mean(2*f_hat_k*S_hat_k - f_hat_k^2 - S_hat_k)
    plug_in2 <- mean(2*f_hat_k^2 - f_hat_k^2 - f_hat_k)
    print(plug_in)
    print(plug_in2)
    # brier_old[i] <- mean(-(f_hat_k - S_hat_k)^2) - mean(-(fs_hat_k - S_hat_k)^2) + mean(if.func_old)
    # brier_plug_old[i] <- mean(-(f_hat_k - S_hat_k)^2) - mean(-(fs_hat_k - S_hat_k)^2)
    #IF.vals[,i] <- if.func
    var_est[i] <- mean(if.func^2)#var(if.func)
  }

  return(data.frame(t = landmark_times,
                    one_step = one_step,
                    plug_in = plug_in,
                    var_est = var_est))
}
