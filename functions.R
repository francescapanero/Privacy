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
  sim <- plyr::aaply(post_probs, 1, function(x) tau1_py_internal(x, frequencies, N, m1), .progress = "text")
  
  cat("Estimated mean: ", round(mean(sim),2),". Monte Carlo se: ", round(sd(sim)/sqrt(R),2),".\n",sep="")
  
  sim
}


# The following function estimates tau1 with PY, DP and Binomial approximation, and corresponding CI
# For the PY estimation, we have two cases:
# 'max_EPPF' uses as point estimate the max of the logEPPF
# 'posterior_samples' samples values from posterior (computed with post_sim_param) and returns 
# tau1_PY estimates and CIs for the last 10 sampled values
# INPUT:
# data: population dataset
# id_uniq_pop: id of the uniques in the population
# percentage: the desired size of the sample
# iter: number of samples you want to take (each of size percentage*N)
# nsamples: number of samples for posterior sampling (input in post_sim_param)
# burnin: not actually a burnin, but how many samples you want to keep
# case: 'max_EPPF' or 'posterior_samples'
# ...: the variables according to which we cross-classify
# OUPUT: 
# table with estimates and CI
loop_estimate <- function(data, id_uniq_pop, percentage, iter, nsamples, burnin=5, case, ...){
  
  N = dim(data)[1]
  output = vector(mode='list', length=iter)
  n = as.integer(percentage*N)
  
  for(i in 1:iter){
    # Create a sample
    ind <- sample(1:N, n, replace = FALSE)
    sample <- data[ind,]
    # Find uniques in sample
    id_uniq_sample <- sample %>% add_count(...) %>% filter(n==1) %>% select(id)
    id_uniq_sample <- id_uniq_sample$id
    # Find frequencies counts in the sample
    table_samp <- sample %>% group_by(...) %>% count()
    freq_samp <- table_samp$n
    # Find uniques in the sample
    m1 <- sum(freq_samp==1)
    # Find tau1_true
    id_uniq_sample_pop <- intersect(id_uniq_sample, id_uniq_pop)
    tau1_true <- length(id_uniq_sample_pop)
    
    if(case == 'max_EPPF'){
      # PY
      out_PY <- max_EPPF_PY(freq_samp)
      tau1_PY <- tau1_py(m1, n, out_PY$par[1], out_PY$par[2], N)
      PY_sim <- tau1_py_sim(freq_samp, out_PY$par[1], out_PY$par[2], N)
      PY_lower <- quantile(PY_sim, 0.005)
      PY_upper <- quantile(PY_sim, 0.995)
      
      # # DP
      # out_DP <- max_EPPF_DP(freq_samp)
      # tau1_DP <- tau1_dp(m1, n, out_DP$par, N)
      # DP_sim  <- tau1_py_sim(freq_samp, out_DP$par, 0, N)
      # DP_lower <- quantile(DP_sim, 0.005)
      # DP_upper <- quantile(DP_sim, 0.995)
      # # Binomial Approx
      # tau1_py_binom <- m1 * (n / N)^(1 - out_PY$par[2])
      # lower_py_binom <- qbinom(0.005, m1, (n / N)^(1 - out_PY$par[2]))
      # upper_py_binom <- qbinom(0.995, m1, (n / N)^(1 - out_PY$par[2]))
      
      output[[i]] <- list(tau1_true=tau1_true, tau1_PY = tau1_PY, #tau1_DP = tau1_DP,  tau1_binom = tau1_py_binom,
                     CI_PY = c(PY_lower, PY_upper), #CI_DP = c(DP_lower, DP_upper), CI_binom = c(lower_py_binom, upper_py_binom),
                     theta=out_PY$par[1], alpha=out_PY$par[2])
    }
    
    if(case == 'posterior_samples'){
      # PY
      out_PY_post <- post_sim_param(freq_samp, nsamples)
      out_PY <- out_PY_post[(nsamples-burnin):nsamples,]
      tau1_PY <- apply(out_PY, 1 , function(x) tau1_py(m1, n, x[1], x[2], N))
      PY_sim <- apply(out_PY, 1, function(x) tau1_py_sim(freq_samp, x[1], x[2], N))
      PY_lower <- apply(PY_sim, 2, function(x) quantile(x, 0.005))
      PY_upper <- apply(PY_sim, 2, function(x) quantile(x, 0.995))
      
      # # DP
      # out_DP <- max_EPPF_DP(freq_samp)
      # tau1_DP <- tau1_dp(m1, n, out_DP$par, N)
      # DP_sim  <- tau1_py_sim(freq_samp, out_DP$par, 0, N)
      # DP_lower <- quantile(DP_sim, 0.005)
      # DP_upper <- quantile(DP_sim, 0.995)
      # # Binomial Approx
      # tau1_py_binom <- m1 * (n / N)^(1 - out_PY$par[2])
      # lower_py_binom <- qbinom(0.005, m1, (n / N)^(1 - out_PY$par[2]))
      # upper_py_binom <- qbinom(0.995, m1, (n / N)^(1 - out_PY$par[2]))
      
      output[[i]] <- list(tau1_true=tau1_true, tau1_PY = mean(tau1_PY), #tau1_DP = tau1_DP,  tau1_binom = tau1_py_binom,
                          CI_PY = c(mean(PY_lower), mean(PY_upper)),
                          #CI_DP = c(DP_lower, DP_upper), CI_binom = c(lower_py_binom, upper_py_binom),
                          theta=mean(out_PY), theta_sample = out_PY_post[,1],
                          alpha=mean(out_PY), alpha_sample = out_PY_post[,2])
    }
  }
  return(output)
}


# This function samples from the posterior of the parameters theta and alpha of the PY
# at the moment the priors are fixed:
# theta \sim LogNormal(0,1)
# alpha \sim Beta(2,2)
# INPUT: 
# freq: frequencies of your sample
# nsamples: desired number of samples
# OUPUT:
# matrix (nsamples X 2) with estimates
post_sim_param <- function(freq, nsamples){
  
  param_estimate <- arms(c(500, 0.9), 
       function(param, a=2, b=2, mu=0, sigma=1, frequencies=freq)  -logEPPF_PY(param[1], param[2], frequencies) +
                                                ((a-1)*log(param[2]) + (b-1)*log(1-param[2]) - lbeta(a,b) -
                                               (log(param[1])-mu)^2/(2*sigma^2)-log(param[1]*sigma*sqrt(2*pi))) ,
       function(param) (param[1]>=1e-10)*(param[1]<=1e+16)*(param[2]>=1e-10)*(param[2]<=1-1e-10),
       nsamples)

  return(param_estimate)
}





