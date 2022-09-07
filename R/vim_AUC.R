vim_AUC <- function(time,
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

  dimension <- 5
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
  top_one_step <- rep(NA, length(landmark_times))
  top_plug_in <- rep(NA, length(landmark_times))
  top_var_est <- rep(NA, length(landmark_times))
  bottom_one_step <- rep(NA, length(landmark_times))
  bottom_plug_in <- rep(NA, length(landmark_times))
  bottom_var_est <- rep(NA, length(landmark_times))
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
      int_vals <- ifelse(f_hat_k > fx, 1, 0) * (1 - S_hat_k) - ifelse(f_hat_k < fx, 1, 0) * S_hat_k
      int <- mean(ifelse(f_hat_k > fx, 1, 0) * (1 - S_hat_k) - ifelse(f_hat_k < fx, 1, 0) * S_hat_k)
      return(varphi_x*int)
    }

    calc_phi_tilde_01 <- function(j){
      fx <- f_hat_k[j]
      Sx <- S_hat_k[j]
      int_vals <- ifelse(f_hat_k > fx, 1, 0) * (1 - S_hat_k) * Sx + ifelse(f_hat_k < fx, 1, 0) * S_hat_k * (1 - Sx)
      int <- mean(ifelse(f_hat_k > fx, 1, 0) * (1 - S_hat_k) * Sx + ifelse(f_hat_k < fx, 1, 0) * S_hat_k * (1 - Sx))
      return(int)
    }

    phi_01 <- unlist(lapply(1:nrow(X_holdout),
                            FUN = calc_phi_01))

    phi_tilde_01_uncentered <- unlist(lapply(1:nrow(X_holdout),
                            FUN = calc_phi_tilde_01))

    phi_tilde_01 <- phi_tilde_01_uncentered - mean(phi_tilde_01_uncentered)

    calc_phi_02 <- function(j){
      varphi_x <- KM.if[j]
      int <- mean((1 - S_hat_k) - S_hat_k)
      return(varphi_x * int)
    }

    calc_phi_tilde_02 <- function(j){
      Sx <- S_hat_k[j]
      int <- mean((1 - S_hat_k) * Sx + S_hat_k * (1 - Sx))
      return(int)
    }

    phi_02 <- unlist(lapply(1:nrow(X_holdout),
                            FUN = calc_phi_02))

    phi_tilde_02_uncentered <- unlist(lapply(1:nrow(X_holdout),
                                             FUN = calc_phi_tilde_02))

    theta_0 <- mean(S_hat_k)
    V_2_new <- theta_0 * (1 - theta_0)
    phi_tilde_02_new <- S_hat_k - theta_0 - 2*theta_0*S_hat_k + 2*theta_0^2
    phi_02_new <- KM.if * (1 - 2*theta_0)


    phi_tilde_02 <- phi_tilde_02_uncentered - mean(phi_tilde_02_uncentered)

    V_1 <- mean(phi_tilde_01_uncentered)
    V_2 <- mean(phi_tilde_02_uncentered)

    # if.func <- (phi_01 + phi_tilde_01)
    # plug_in[i] <- mean(phi_tilde_01_uncentered)
    # one_step[i] <- mean(phi_tilde_01_uncentered) + mean(if.func)
    top_plug_in[i] <- V_1
    top_one_step[i] <- V_1 + mean(phi_01 + phi_tilde_01)
    top_var_est[i] <- mean((phi_01 + phi_tilde_01)^2)
    bottom_plug_in[i] <- V_2_new#V_2
    bottom_one_step[i] <- V_2_new + mean(phi_tilde_02_new + phi_02_new)#V_2 + mean(phi_02 + phi_tilde_02)
    bottom_var_est[i] <- mean((phi_02_new + phi_tilde_02_new)^2)
    # if.func <- (phi_01 + phi_tilde_01)/V_2 - V_1/(V_2^2)*(phi_02 + phi_tilde_02)
    # plug_in[i] <- V_1/V_2
    # one_step[i] <- V_1/V_2 + mean(if.func)
    #
    # var_est[i] <- mean(if.func^2)
  }

  return(data.frame(t = landmark_times,
                    top_one_step = top_one_step,
                    top_plug_in = top_plug_in,
                    top_var_est = top_var_est,
                    bottom_one_step = bottom_one_step,
                    bottom_plug_in = bottom_plug_in,
                    bottom_var_est = bottom_var_est))
}
