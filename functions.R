# library(plyr)
library(ggplot2)
Rcpp::sourceCpp('cluster_py.cpp')

rdir <- function(R, param, log_scale=FALSE) {
  H   <- length(param)

  sim <- matrix(rgamma(R*H,param), ncol=H, byrow=TRUE)
  sum <- sim %*% rep(1,H)
  sim/as.vector(sum)
}

logEPPF_PY <- function(theta, alpha, frequencies){
  
  # Sample size
  n   <- sum(frequencies)
  # Number of frequencies
  K   <- length(frequencies)
  # Loglikelihood
  loglik <- sum(log(theta + 1:(K-1)*alpha)) - lgamma(theta + n) + lgamma(theta + 1) + sum(lgamma(frequencies - alpha)) - K*lgamma(1 - alpha)
  
  loglik
}

logEPPF_DP <- function(theta,  frequencies){
  
  # Sample size
  n   <- sum(frequencies)
  # Number of frequencies
  K   <- length(frequencies)
  
  # Loglikelihood
  loglik <- K*log(theta) - lgamma(theta + n) + lgamma(theta) 
  
  loglik 
}

max_EPPF_PY <- function(frequencies){
  start <- c(1,0.5) # Initialization of the maximization algorithm
  out <- nlminb(start= start, 
                function(param) -logEPPF_PY(theta=param[1], alpha=param[2], frequencies = frequencies), 
                lower=c(1e-10,1e-10), upper=c(Inf,1-1e-10))  
  return(out)
}


max_EPPF_DP <- function(frequencies){
  start <- 1 # Initialization of the maximization algorithm
  out <- nlminb(start= start, 
                function(param) - logEPPF_DP(theta=param, frequencies = frequencies), 
                lower=1e-10, upper=Inf)  
  return(out)
}


tau1_dp <- function(m1, n, theta, N){
  m1*(n + theta -1)/(N + theta -1)
}

tau1_py <- function(m1, n, theta, alpha, N){
  #exp(log(m1) + lgamma(theta - 1 + alpha + N) - lgamma(theta + n - 1 + alpha) - lgamma(theta + N) + lgamma(theta + n))
  exp(log(m1) + lgamma(theta + alpha + N) - lgamma(theta + n - 1 + alpha) - lgamma(theta + N + 1) + lgamma(theta + n)
      + log(theta+N) - log(theta+N+alpha-1))
}

tau1_py_sim2 <- function(frequencies, theta, alpha, N, R = 100, verbose=TRUE){
  
  # Main quantities
  m1  <- sum(frequencies == 1)
  K_n <- length(frequencies)
  
  # First sample the probabilities
  cat("Pre-computing the sampling probabilities. ")
  
  post_probs <- rdir(R, c(frequencies - alpha, theta + K_n*alpha))
  
  cat("Done. \n")
  
  tau1_py_internal <-  function(prob, frequencies, N, m1){
    
    freq_temp    <- c(frequencies,0) # Add the extra value 
    index_obs    <- freq_temp == 1
    
    # This vector has length m1 + 1. The probability each singleton has
    prob_m1 <- c(prob[index_obs], sum(prob[!index_obs]))
    
    # Most of the times you do not even sample a singleton
    X_pred <- sample(1:length(prob_m1), size = N - sum(frequencies), prob = prob_m1, replace=TRUE)
    
    m1 - length(unique(X_pred[X_pred != m1 + 1]))
    
  }
  
  # Apply the above function to each row
  sim <- plyr::aaply(post_probs, 1, function(x) tau1_py_internal(x, frequencies, N, m1), .progress = "text")
  
  cat("Estimated mean: ", round(mean(sim),2),". Monte Carlo se: ", round(sd(sim)/sqrt(R),2),".\n",sep="")
  
  sim
}



tau1_py_sim <- function(frequencies, theta, alpha, N, R = 100, verbose=TRUE){
  
  # Main quantities
  n   <- sum(frequencies)
  m1  <- sum(frequencies == 1)
  K_n <- length(frequencies)
  
  K_sim <- plyr::aaply(matrix(0,R,2), 1, function(x) cluster_py_C(N - n, 1 - alpha, theta+n), .progress = "text")
  #K_sim <- replicate(R,cluster_py_C(N - n, 1 - alpha, theta+n))
  
  sim <- rhyper(R, (theta + n)/(1-alpha) - 1, K_sim, m1)
  sim
  
  cat("Estimated mean: ", round(mean(sim),2),". Monte Carlo se: ", round(sd(sim)/sqrt(R),2),".\n",sep="")
  sim
}


expected_cl_py <- function(n, alpha, theta){
  n <- as.integer(n)
  if(alpha==0) {
      out <- theta * sum(1/(theta - 1 + 1:n))
  } else {
      out <- 1/alpha*exp(lgamma(theta + alpha + n) - lgamma(theta + alpha) - lgamma(theta + n) + lgamma(theta + 1)) - theta/alpha
  }

  return(out)
}

expected_m1 <- function(n, alpha, theta){
  n/(n + theta - 1)*(theta + alpha*expected_cl_py(n-1,alpha,theta))
}

model_checking_PY <- function(frequencies, percentage=0.75, step = 100){
  
  # Total number of observations
  N   <- sum(frequencies)
  K_N <- length(frequencies)
  
  points_full   <- factor(rep(1:length(frequencies), frequencies)) # it is important they are factor
  points_full   <- points_full[sample(N)] # Randomly permute the data
  freq_full     <- as.numeric(table(points_full))
  
  freq_observed  <- as.numeric(table(points_full[1:round(percentage*N)]))
  n              <- sum(freq_observed)
  m1             <- sum(freq_observed==1)
  
  nn <- seq.int(from = n, to = N, length.out=step)
  
  fit_PY    <- max_EPPF_PY(freq_observed[freq_observed >0])

  est_tau1  <- numeric(step)
  true_tau1 <- numeric(step)
  
  for(i in 1:step){
    
    points_full1  <- points_full[1:nn[i]]
    
    # Observed frequencies
    freq_full1 <- as.numeric(table(points_full1))
    
    # True tau1
    est_tau1[i]   <- tau1_py(m1 = m1, n = n, theta = fit_PY$par[1], alpha = fit_PY$par[2], N = nn[i])
    true_tau1[i]  <- sum((freq_full1 == 1) & (freq_observed == 1))
  } 
  
  data_plot <- data.frame(nn = nn, est_tau1 = est_tau1, true_tau1 = true_tau1)
  p <- ggplot(data_plot, aes(x = nn, y = est_tau1))  + theme_bw() + xlab("# of observations") + ylab(expression(tau[1])) + geom_point(aes(y = true_tau1))+ geom_line(col="gold")
  p
}

model_checking_DP <- function(frequencies, percentage=0.75, step = 100){
  
  # Total number of observations
  N   <- sum(frequencies)
  K_N <- length(frequencies)
  
  points_full   <- factor(rep(1:length(frequencies), frequencies)) # it is important they are factor
  points_full   <- points_full[sample(N)] # Randomly permute the data
  freq_full     <- as.numeric(table(points_full))
  
  freq_observed  <- as.numeric(table(points_full[1:round(percentage*N)]))
  n              <- sum(freq_observed)
  m1             <- sum(freq_observed==1)
  
  nn <- seq.int(from = n, to = N, length.out=step)
  
  fit_DP    <- max_EPPF_DP(freq_observed[freq_observed >0])
  
  est_tau1  <- numeric(step)
  true_tau1 <- numeric(step)
  lower_tau1 <- numeric(step)
  upper_tau1 <- numeric(step)
  
  for(i in 1:step){
    
    points_full1  <- points_full[1:nn[i]]
    
    # Observed frequencies
    freq_full1 <- as.numeric(table(points_full1))
    
    # True tau1
    est_tau1[i]   <- tau1_dp(m1 = m1, n = n, theta = fit_DP$par[1], N = nn[i])
    true_tau1[i]  <- sum((freq_full1 == 1) & (freq_observed == 1))
    lower_tau1[i] <- qhyper(0.025, fit_DP$par[1] + n - 1, nn[i] -n, m1)
    upper_tau1[i] <- qhyper(0.975, fit_DP$par[1] + n - 1, nn[i] - n, m1)
    
  } 
  
  data_plot <- data.frame(nn = nn, est_tau1 = est_tau1, true_tau1 = true_tau1, lower_tau1 = lower_tau1, upper_tau1 = upper_tau1)
  p <- ggplot(data_plot, aes(x = nn, y = est_tau1, ymin = lower_tau1, ymax = upper_tau1))  + theme_bw() + xlab("# of observations") + ylab(expression(tau[1])) + geom_point(aes(y = true_tau1))+ geom_line(col="blue") + geom_ribbon(alpha=0.1)
  p
}

frequency_check_PY <- function(frequencies){
  
  n      <- sum(frequencies)
  fit_PY <- max_EPPF_PY(frequencies)
  sigma  <- fit_PY$par[2]
  
  
  M_l <- as.numeric(table(factor(frequencies, levels=1:n)))
  P_l <- M_l / sum(M_l)
  
  idx <- which(P_l > 0)
  
  data_plot <- data.frame(Size = idx, P_l = P_l[idx], Theoretical = exp(log(sigma) + lgamma(idx - sigma) - lgamma(1-sigma) - lfactorial(idx)))
  p <- ggplot(data=data_plot,aes(x = Size, y = P_l))  + geom_point() +  geom_line(aes(y = Theoretical),color="blue", linetype="dashed")+ scale_y_log10() + scale_x_log10() + theme_bw() + xlab(expression(M[ln])) + ylab(expression(M[ln]/K[n]))
  p                             
}


