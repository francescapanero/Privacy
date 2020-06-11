rm(list=ls())
library(knitr)

source('functions.R')

dataset_creation_zipf <- function(n, zipf_param, H, N){
  
  # Expanded representation
  probs =  1/(1:H)^(zipf_param) #
  
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
       m1 = sum(freq_observed==1), percentage = n/N, zipf_param = zipf_param, true_tau1 = true_tau1)
  
}

# -------------------------------------------
# Scenario 1 - zipf_param = 2
# -------------------------------------------
   
H          <- 10^7
N          <- 10^6
n          <- 100000
zipf_param <- 2

set.seed(123)
dataset <- dataset_creation_zipf(n = n, zipf_param = zipf_param, H = H, N = N)

kable(data.frame(n = dataset$n, N = dataset$N, percentage = dataset$percentage, 
                 K_n =dataset$K_n,  K_N = dataset$K_N, H = dataset$H, m1 =dataset$m1, 
                 true_tau1 =dataset$true_tau1, 
                 zipf_param = dataset$zipf_param))

model_checking_PY(dataset$frequencies)
model_checking_DP(dataset$frequencies)

# 1 ] DP estimation
max_EPPF_DP(dataset$frequencies)

# Theta estimation via MLE
theta_hat <- max_EPPF_DP(dataset$frequencies)$par

# Estimation of tau1_dp
dataset$tau1_dp <- tau1_dp(dataset$m1, dataset$n, theta_hat, dataset$N)

dataset$true_tau1
dataset$tau1_dp

# EXACT uncertainty quantification using MonteCarlo
sim_dp  <- tau1_py_sim(dataset$frequencies, theta_hat, 0, dataset$N, R = 1000)
dataset$lower_dp <- quantile(sim_dp, 0.025)
dataset$upper_dp <- quantile(sim_dp, 0.975)

# # Binomial approximation
# dataset$tau1_dp_binom <- dataset$m1 * dataset$n / dataset$N
# dataset$lower_dp_binom <- qbinom(0.025,dataset$m1, dataset$n / dataset$N)
# dataset$upper_dp_binom <- qbinom(0.975,dataset$m1, dataset$n / dataset$N)

# 2 ] PY estimation

max_EPPF_PY(dataset$frequencies)

param_hat <- max_EPPF_PY(dataset$frequencies)
theta_hat <- param_hat$par[1]
alpha_hat <- param_hat$par[2]

dataset$tau1_py <- tau1_py(dataset$m1, dataset$n, theta_hat, alpha_hat, dataset$N)

dataset$true_tau1
dataset$tau1_py 

# Uncertainty quantification
sim_py <- tau1_py_sim(dataset$frequencies, theta_hat, alpha_hat, dataset$N, R = 1000)
dataset$lower_py <- quantile(sim_py, 0.025)
dataset$upper_py <- quantile(sim_py, 0.975)


# # Binomial approximation
# dataset$tau1_py_binom <- dataset$m1 * (dataset$n / dataset$N)^(1 - alpha_hat)
# dataset$lower_py_binom <- qbinom(0.025,dataset$m1, (dataset$n / dataset$N)^(1 - alpha_hat))
# dataset$upper_py_binom <- qbinom(0.975,dataset$m1, (dataset$n / dataset$N)^(1 - alpha_hat))

# Results -----------------------------------

kable(data.frame(true_tau1 = dataset$true_tau1, 
                 tau_dp       = dataset$tau1_dp,  CI_dp = paste("[", dataset$lower_dp,", ", dataset$upper_dp,"]",sep=""),
                 tau_py       = dataset$tau1_py, CI_py = paste("[", dataset$lower_py,", ", dataset$upper_py,"]",sep="")
))

# -------------------------------------------
# Scenario 1 - zipf_param = 1.5
# -------------------------------------------

H          <- 10^7
N          <- 10^6
n          <- 100000
zipf_param <- 1.5

set.seed(123)
dataset <- dataset_creation_zipf(n = n, zipf_param = zipf_param, H = H, N = N)

kable(data.frame(n = dataset$n, N = dataset$N, percentage = dataset$percentage, 
                 K_n =dataset$K_n,  K_N = dataset$K_N, H = dataset$H, m1 =dataset$m1, 
                 true_tau1 =dataset$true_tau1, 
                 zipf_param = dataset$zipf_param))

model_checking_PY(dataset$frequencies, percentage = 0.8)
model_checking_DP(dataset$frequencies)

# 1 ] DP estimation

max_EPPF_DP(dataset$frequencies)

# Theta estimation via MLE
theta_hat <- max_EPPF_DP(dataset$frequencies)$par

# Estimation of tau1_dp
dataset$tau1_dp <- tau1_dp(dataset$m1, dataset$n, theta_hat, dataset$N)

dataset$true_tau1
dataset$tau1_dp

# EXACT uncertainty quantification using MonteCarlo
sim_dp  <- tau1_py_sim(dataset$frequencies, theta_hat, 0, dataset$N, R = 1000)
dataset$lower_dp <- quantile(sim_dp, 0.025)
dataset$upper_dp <- quantile(sim_dp, 0.975)

# # Binomial approximation
# dataset$tau1_dp_binom <- dataset$m1 * dataset$n / dataset$N
# dataset$lower_dp_binom <- qbinom(0.025,dataset$m1, dataset$n / dataset$N)
# dataset$upper_dp_binom <- qbinom(0.975,dataset$m1, dataset$n / dataset$N)

# 2 ] PY estimation

param_hat <- max_EPPF_PY(dataset$frequencies)
theta_hat <- param_hat$par[1]
alpha_hat <- param_hat$par[2]

dataset$tau1_py <- tau1_py(dataset$m1, dataset$n, theta_hat, alpha_hat, dataset$N)

dataset$true_tau1
dataset$tau1_py 

# Uncertainty quantification
sim_py <- tau1_py_sim(dataset$frequencies, theta_hat, alpha_hat, dataset$N, R = 1000)
dataset$lower_py <- quantile(sim_py, 0.025)
dataset$upper_py <- quantile(sim_py, 0.975)

# # Binomial approximation
# dataset$tau1_py_binom <- dataset$m1 * (dataset$n / dataset$N)^(1 - alpha_hat)
# dataset$lower_py_binom <- qbinom(0.025,dataset$m1, (dataset$n / dataset$N)^(1 - alpha_hat))
# dataset$upper_py_binom <- qbinom(0.975,dataset$m1, (dataset$n / dataset$N)^(1 - alpha_hat))

# Results -----------------------------------

kable(data.frame(true_tau1 = dataset$true_tau1, 
                 tau_dp       = dataset$tau1_dp,  CI_dp = paste("[", dataset$lower_dp,", ", dataset$upper_dp,"]",sep=""),
                 tau_py       = dataset$tau1_py, CI_py = paste("[", dataset$lower_py,", ", dataset$upper_py,"]",sep="")
))

# -------------------------------------------
# Scenario 1 - zipf_param = 1.3
# -------------------------------------------

H          <- 10^7
N          <- 10^6
n          <- 20000
zipf_param <- 1.3

set.seed(123)
dataset <- dataset_creation_zipf(n = n, zipf_param = zipf_param, H = H, N = N)

kable(data.frame(n = dataset$n, N = dataset$N, percentage = dataset$percentage, 
                 K_n =dataset$K_n,  K_N = dataset$K_N, H = dataset$H, m1 =dataset$m1, 
                 true_tau1 =dataset$true_tau1, 
                 zipf_param = dataset$zipf_param))

# 1 ] DP estimation

max_EPPF_DP(dataset$frequencies)

# Theta estimation via MLE
theta_hat <- max_EPPF_DP(dataset$frequencies)$par

# Estimation of tau1_dp
dataset$tau1_dp <- tau1_dp(dataset$m1, dataset$n, theta_hat, dataset$N)

dataset$true_tau1
dataset$tau1_dp

# EXACT uncertainty quantification using MonteCarlo
sim_dp  <- tau1_py_sim(dataset$frequencies, theta_hat, 0, dataset$N, R = 1000)
dataset$lower_dp <- quantile(sim_dp, 0.025)
dataset$upper_dp <- quantile(sim_dp, 0.975)

# # Binomial approximation
# dataset$tau1_dp_binom <- dataset$m1 * dataset$n / dataset$N
# dataset$lower_dp_binom <- qbinom(0.025,dataset$m1, dataset$n / dataset$N)
# dataset$upper_dp_binom <- qbinom(0.975,dataset$m1, dataset$n / dataset$N)

# 2 ] PY estimation

param_hat <- max_EPPF_PY(dataset$frequencies)
theta_hat <- param_hat$par[1]
alpha_hat <- param_hat$par[2]

dataset$tau1_py <- tau1_py(dataset$m1, dataset$n, theta_hat, alpha_hat, dataset$N)

dataset$true_tau1
dataset$tau1_py 

# Uncertainty quantification
sim_py <- tau1_py_sim(dataset$frequencies, theta_hat, alpha_hat, dataset$N, R = 1000)
dataset$lower_py <- quantile(sim_py, 0.025)
dataset$upper_py <- quantile(sim_py, 0.975)


# # Binomial approximation
# dataset$tau1_py_binom <- dataset$m1 * (dataset$n / dataset$N)^(1 - alpha_hat)
# dataset$lower_py_binom <- qbinom(0.025,dataset$m1, (dataset$n / dataset$N)^(1 - alpha_hat))
# dataset$upper_py_binom <- qbinom(0.975,dataset$m1, (dataset$n / dataset$N)^(1 - alpha_hat))


# Results -----------------------------------
kable(data.frame(true_tau1 = dataset$true_tau1, 
                 tau_dp       = dataset$tau1_dp,  CI_dp = paste("[", dataset$lower_dp,", ", dataset$upper_dp,"]",sep=""),
                 tau_py       = dataset$tau1_py, CI_py = paste("[", dataset$lower_py,", ", dataset$upper_py,"]",sep="")
))

# -------------------------------------------
# Scenario 1 - zipf_param = 1.1
# -------------------------------------------

H          <- 10^7
N          <- 10^6
n          <- 100000
zipf_param <- 1.1

set.seed(123)
dataset <- dataset_creation_zipf(n = n, zipf_param = zipf_param, H = H, N = N)

kable(data.frame(n = dataset$n, N = dataset$N, percentage = dataset$percentage, 
                 K_n =dataset$K_n,  K_N = dataset$K_N, H = dataset$H, m1 =dataset$m1, 
                 true_tau1 =dataset$true_tau1, 
                 zipf_param = dataset$zipf_param))

# 1 ] DP estimation

max_EPPF_DP(dataset$frequencies)

# Theta estimation via MLE
theta_hat <- max_EPPF_DP(dataset$frequencies)$par

# Estimation of tau1_dp
dataset$tau1_dp <- tau1_dp(dataset$m1, dataset$n, theta_hat, dataset$N)

dataset$true_tau1
dataset$tau1_dp

# EXACT uncertainty quantification using MonteCarlo
sim_dp  <- tau1_py_sim(dataset$frequencies, theta_hat, 0, dataset$N, R = 1000)
dataset$lower_dp <- quantile(sim_dp, 0.025)
dataset$upper_dp <- quantile(sim_dp, 0.975)

# # Binomial approximation
# dataset$tau1_dp_binom <- dataset$m1 * dataset$n / dataset$N
# dataset$lower_dp_binom <- qbinom(0.025,dataset$m1, dataset$n / dataset$N)
# dataset$upper_dp_binom <- qbinom(0.975,dataset$m1, dataset$n / dataset$N)

# 2 ] PY estimation

param_hat <- max_EPPF_PY(dataset$frequencies)
theta_hat <- param_hat$par[1]
alpha_hat <- param_hat$par[2]

dataset$tau1_py <- tau1_py(dataset$m1, dataset$n, theta_hat, alpha_hat, dataset$N)

dataset$true_tau1
dataset$tau1_py 

# Uncertainty quantification
sim_py <- tau1_py_sim(dataset$frequencies, theta_hat, alpha_hat, dataset$N, R = 1000)
dataset$lower_py <- quantile(sim_py, 0.025)
dataset$upper_py <- quantile(sim_py, 0.975)


# # Binomial approximation
# dataset$tau1_py_binom <- dataset$m1 * (dataset$n / dataset$N)^(1 - alpha_hat)
# dataset$lower_py_binom <- qbinom(0.025,dataset$m1, (dataset$n / dataset$N)^(1 - alpha_hat))
# dataset$upper_py_binom <- qbinom(0.975,dataset$m1, (dataset$n / dataset$N)^(1 - alpha_hat))


# Results -----------------------------------

kable(data.frame(true_tau1 = dataset$true_tau1, 
                 tau_dp       = dataset$tau1_dp,  CI_dp = paste("[", dataset$lower_dp,", ", dataset$upper_dp,"]",sep=""),
                 tau_py       = dataset$tau1_py, CI_py = paste("[", dataset$lower_py,", ", dataset$upper_py,"]",sep="")
))

