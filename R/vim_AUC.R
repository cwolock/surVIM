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
  #IF.vals <- matrix(NA, nrow=n, ncol=length(landmark_times))
  AUC <- rep(NA, length(landmark_times))
  AUC_plug <- rep(NA, length(landmark_times))
  S_t <- rep(NA, length(landmark_times))
  G_t <- rep(NA, length(landmark_times))
  var_est <- rep(NA, length(landmark_times))
  for(i in 1:length(landmark_times)) {
    t0 <- landmark_times[i]
    k <- min(which(approx_times >= t0))
    S_hat_k <- S_hat[,k]
    G_hat_k <- G_hat[,k]
    f_hat_k <- 1- f_hat[,i] # oracle is CDF, not survival
    #fs_hat_k <- fs_hat[,i]
    inner.func.1 <- ifelse(time <= t0 & event == 1, 1/(S_hat_Y * G_hat_Y), 0 )
    inner.func.2 <- int.vals[,k]
    KM.if <- -S_hat_k * ( inner.func.1 - inner.func.2)

    calc_phi_01 <- function(j){
      fx <- f_hat_k[j]
      varphi_x <- KM.if[j]
      int <- mean(ifelse(f_hat_k >= fx, 1, 0) * (1 - S_hat_k) - ifelse(f_hat_k < fx, 1, 0) * S_hat_k)
      return(varphi_x * int)
    }

    calc_phi_tilde_01 <- function(j){
      fx <- f_hat_k[j]
      Sx <- S_hat_k[j]
      int <- mean(ifelse(f_hat_k >= fx, 1, 0) * (1 - S_hat_k) * S_x - ifelse(f_hat_k < fx, 1, 0) * S_hat_k * (1 - Sx))
      return(int)
    }

    phi_01 <- unlist(lapply(1:nrow(X_holdout),
                            FUN = calc_phi_01))

    phi_tilde_01_uncentered <- unlist(lapply(1:nrow(X_holdout),
                            FUN = calc_phi_tilde_01))

    phi_tilde_01 <- phi_tilde_01_uncentered - 2*mean(phi_tilde_01)

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

    brier[i] <- mean(2*f_hat_k*S_hat_k - f_hat_k^2 - S_hat_k) + mean(if.func) #-
      #mean(2*fs_hat_k*S_hat_k - fs_hat_k^2 - S_hat_k)
    brier_plug[i] <- mean(2*f_hat_k*S_hat_k - f_hat_k^2 - S_hat_k)# -
      #mean(2*fs_hat_k*S_hat_k - fs_hat_k^2 - S_hat_k)
    # brier_old[i] <- mean(-(f_hat_k - S_hat_k)^2) - mean(-(fs_hat_k - S_hat_k)^2) + mean(if.func_old)
    # brier_plug_old[i] <- mean(-(f_hat_k - S_hat_k)^2) - mean(-(fs_hat_k - S_hat_k)^2)
    #IF.vals[,i] <- if.func
    S_t[i] <- mean(S_hat_k)
    G_t[i] <- mean(G_hat_k)
    var_est[i] <- mean(if.func^2)#var(if.func)
  }

  return(data.frame(t = landmark_times,
                    brier = brier,
                    brier_plug = brier_plug,
                    S_t = S_t,
                    G_t = G_t,
                    var_est = var_est))
}
