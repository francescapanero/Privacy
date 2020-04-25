# library(plyr)

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

dataset_creation <- function(n, alpha, H, N){
  
  # Expanded representation
  probs =  1/(1:H)^(alpha) #
  #points_obs  <- sample(1:H, n, prob = probs, replace=TRUE)
  #points_test <- sample(1:H, N - n, prob = probs, replace=TRUE)
  #points_full <- factor(c(points_obs,points_test),levels=1:H)
  #points_obs  <- factor(points_obs, levels=1:H)
  
  points_full <- factor(sample(1:H, N, prob = probs, replace=TRUE),levels=1:H)
  points_obs  <- sample(points_full, n, replace=FALSE)
  
  # Observed frequencies
  freq_full     <- as.numeric(table(points_full))
  freq_observed <- as.numeric(table(points_obs))
  
  # True tau1
  true_tau1 <- sum((freq_full == 1) & (freq_observed == 1))
  freq_observed <- freq_observed[freq_observed >0]
  
  list(frequencies = freq_observed, n = sum(freq_observed), H = H, N = sum(freq_full), 
       K_n = sum(freq_observed>0), K_N = sum(freq_full>0), 
       m1 = sum(freq_observed==1), percentage = n/N, alpha = alpha, true_tau1 = true_tau1)
  
}

tau1_dp <- function(m1, n, theta, N){
  m1*(n + theta -1)/(N + theta -1)
}

tau1_py <- function(m1, n, theta, alpha, N){
  #exp(log(m1) + lgamma(theta - 1 + alpha + N) - lgamma(theta + n - 1 + alpha) - lgamma(theta + N) + lgamma(theta + n))
  exp(log(m1) + lgamma(theta + alpha + N) - lgamma(theta + n - 1 + alpha) - lgamma(theta + N + 1) + lgamma(theta + n)
      + log(theta+N) - log(theta+N+alpha-1))
}

tau1_py_sim <- function(frequencies, theta, alpha, N, R = 100, verbose=TRUE){
  
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
  sim <- aaply(post_probs, 1, function(x) tau1_py_internal(x, frequencies, N, m1), .progress = "text")
  
  cat("Estimated mean: ", round(mean(sim),2),". Monte Carlo se: ", round(sd(sim)/sqrt(R),2),".\n",sep="")
  
  sim
}



