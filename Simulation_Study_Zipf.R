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

# Parameter estimation
out_PY  <- max_EPPF_PY(dataset$frequencies)
tau1_PY <- tau1_py(dataset$m1, dataset$n, out_PY$par[1], out_PY$par[2], dataset$N)
PY_sim  <- tau1_py_sim(dataset$frequencies, out_PY$par[1], out_PY$par[2],  dataset$N)
#PY_sim2  <- tau1_py_sim2(dataset$frequencies, out_PY$par[1], out_PY$par[2],  dataset$N)
PY_lower <- quantile(PY_sim, 0.01/2)
PY_upper <- quantile(PY_sim, 1 - 0.01/2)

# Dirichlet process estimation
out_DP    <- max_EPPF_DP(dataset$frequencies)
tau1_DP   <- tau1_dp(dataset$m1, dataset$n, out_DP$par[1], dataset$N)
DP_lower <- qhyper(0.01/2, out_DP$par[1] + dataset$n - 1, dataset$N - dataset$n, dataset$m1)
DP_upper <- qhyper(1 - 0.01/2, out_DP$par[1] + dataset$n - 1, dataset$N - dataset$n, dataset$m1)

# Summary
kable(data.frame(n = dataset$n, N = dataset$N, percentage = round(dataset$percentage,2),
                 m1 = dataset$m1,
                 K_n = dataset$K_n,
                 true_tau1 = dataset$true_tau1,
                 tau1_py = tau1_PY, CI_PY = paste("[", PY_lower,", ", PY_upper,"]",sep=""),
                 tau1_dp = tau1_DP, CI_DP = paste("[", DP_lower,", ", DP_upper,"]",sep="")))


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

set.seed(1234)
model_checking_PY(dataset$frequencies)
model_checking_DP(dataset$frequencies)

# Parameter estimation
out_PY  <- max_EPPF_PY(dataset$frequencies)
tau1_PY <- tau1_py(dataset$m1, dataset$n, out_PY$par[1], out_PY$par[2], dataset$N)
PY_sim  <- tau1_py_sim(dataset$frequencies, out_PY$par[1], out_PY$par[2],  dataset$N, R = 1000)
# PY_sim2  <- tau1_py_sim2(dataset$frequencies, out_PY$par[1], out_PY$par[2],  dataset$N)

#hist(PY_sim)
#hist(PY_sim2)

PY_lower <- quantile(PY_sim, 0.01/2)
PY_upper <- quantile(PY_sim, 1 - 0.01/2)

# Dirichlet process estimation
out_DP    <- max_EPPF_DP(dataset$frequencies)
tau1_DP   <- tau1_dp(dataset$m1, dataset$n, out_DP$par[1], dataset$N)
DP_lower <- qhyper(0.01/2, out_DP$par[1] + dataset$n - 1, dataset$N - dataset$n, dataset$m1)
DP_upper <- qhyper(1 - 0.01/2, out_DP$par[1] + dataset$n - 1, dataset$N - dataset$n, dataset$m1)

# Summary
kable(data.frame(n = dataset$n, N = dataset$N, percentage = round(dataset$percentage,2),
                 m1 = dataset$m1,
                 K_n = dataset$K_n,
                 true_tau1 = dataset$true_tau1,
                 tau1_py = tau1_PY, CI_PY = paste("[", PY_lower,", ", PY_upper,"]",sep=""),
                 tau1_dp = tau1_DP, CI_DP = paste("[", DP_lower,", ", DP_upper,"]",sep="")))

# -------------------------------------------
# Scenario 1 - zipf_param = 1.3
# -------------------------------------------

H          <- 10^7
N          <- 10^6
n          <- 100000
zipf_param <- 1.3

set.seed(123)
dataset <- dataset_creation_zipf(n = n, zipf_param = zipf_param, H = H, N = N)

kable(data.frame(n = dataset$n, N = dataset$N, percentage = dataset$percentage, 
                 K_n =dataset$K_n,  K_N = dataset$K_N, H = dataset$H, m1 =dataset$m1, 
                 true_tau1 =dataset$true_tau1, 
                 zipf_param = dataset$zipf_param))

set.seed(1234)
model_checking_PY(dataset$frequencies)
model_checking_DP(dataset$frequencies)

# Parameter estimation
out_PY  <- max_EPPF_PY(dataset$frequencies)
tau1_PY <- tau1_py(dataset$m1, dataset$n, out_PY$par[1], out_PY$par[2], dataset$N)
PY_sim  <- tau1_py_sim(dataset$frequencies, out_PY$par[1], out_PY$par[2],  dataset$N, R = 1000)
# PY_sim2  <- tau1_py_sim2(dataset$frequencies, out_PY$par[1], out_PY$par[2],  dataset$N)

#hist(PY_sim)
#hist(PY_sim2)

PY_lower <- quantile(PY_sim, 0.01/2)
PY_upper <- quantile(PY_sim, 1 - 0.01/2)

# Dirichlet process estimation
out_DP    <- max_EPPF_DP(dataset$frequencies)
tau1_DP   <- tau1_dp(dataset$m1, dataset$n, out_DP$par[1], dataset$N)
DP_lower <- qhyper(0.01/2, out_DP$par[1] + dataset$n - 1, dataset$N - dataset$n, dataset$m1)
DP_upper <- qhyper(1 - 0.01/2, out_DP$par[1] + dataset$n - 1, dataset$N - dataset$n, dataset$m1)

# Summary
kable(data.frame(n = dataset$n, N = dataset$N, percentage = round(dataset$percentage,2),
                 m1 = dataset$m1,
                 K_n = dataset$K_n,
                 true_tau1 = dataset$true_tau1,
                 tau1_py = tau1_PY, CI_PY = paste("[", PY_lower,", ", PY_upper,"]",sep=""),
                 tau1_dp = tau1_DP, CI_DP = paste("[", DP_lower,", ", DP_upper,"]",sep="")))

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
set.seed(1234)
model_checking_PY(dataset$frequencies)
model_checking_DP(dataset$frequencies)

# Parameter estimation
out_PY  <- max_EPPF_PY(dataset$frequencies)
tau1_PY <- tau1_py(dataset$m1, dataset$n, out_PY$par[1], out_PY$par[2], dataset$N)
PY_sim  <- tau1_py_sim(dataset$frequencies, out_PY$par[1], out_PY$par[2],  dataset$N, R = 1000)
# PY_sim2  <- tau1_py_sim2(dataset$frequencies, out_PY$par[1], out_PY$par[2],  dataset$N)

#hist(PY_sim)
#hist(PY_sim2)

PY_lower <- quantile(PY_sim, 0.01/2)
PY_upper <- quantile(PY_sim, 1 - 0.01/2)

# Dirichlet process estimation
out_DP    <- max_EPPF_DP(dataset$frequencies)
tau1_DP   <- tau1_dp(dataset$m1, dataset$n, out_DP$par[1], dataset$N)
DP_lower <- qhyper(0.01/2, out_DP$par[1] + dataset$n - 1, dataset$N - dataset$n, dataset$m1)
DP_upper <- qhyper(1 - 0.01/2, out_DP$par[1] + dataset$n - 1, dataset$N - dataset$n, dataset$m1)

# Summary
kable(data.frame(n = dataset$n, N = dataset$N, percentage = round(dataset$percentage,2),
                 m1 = dataset$m1,
                 K_n = dataset$K_n,
                 true_tau1 = dataset$true_tau1,
                 tau1_py = tau1_PY, CI_PY = paste("[", PY_lower,", ", PY_upper,"]",sep=""),
                 tau1_dp = tau1_DP, CI_DP = paste("[", DP_lower,", ", DP_upper,"]",sep="")))

