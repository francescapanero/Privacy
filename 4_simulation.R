rm(list = ls())

library(knitr)
library(reticulate) # Library that interface R to python
use_condaenv() # Use python anaconda
random <- import("numpy.random") # Import python libraries

source("2_functions.R")

dataset_creation_zipf <- function(n, zipf_param, N) {

  # Old implementation using truncated zipf laws
  # probs <- 1 / (1:H)^(zipf_param) 
  # points_full <- factor(sample(1:H, N, prob = probs, replace = TRUE), levels = 1:H)
  
  # Using python zipf sampler
  points_full <- factor(random$zipf(zipf_param, N))
  points_obs <- sample(points_full, n, replace = FALSE)

  # Observed frequencies
  freq_full <- as.numeric(table(points_full))
  freq_observed <- as.numeric(table(points_obs))

  # True tau1
  true_tau1 <- sum((freq_full == 1) & (freq_observed == 1))
  freq_observed <- freq_observed[freq_observed > 0]

  list(
    frequencies = freq_observed, n = sum(freq_observed), N = sum(freq_full),
    K_n = sum(freq_observed > 0), K_N = sum(freq_full > 0),
    m1 = sum(freq_observed == 1), percentage = n / N, zipf_param = zipf_param, true_tau1 = true_tau1
  )
}

dataset_creation_geom <- function(n, N, p) {
  
  points_full <- factor(rgeom(N, p = p))
  points_obs <- sample(points_full, n, replace = FALSE)
  
  # Observed frequencies
  freq_full <- as.numeric(table(points_full))
  freq_observed <- as.numeric(table(points_obs))
  
  # True tau1
  true_tau1 <- sum((freq_full == 1) & (freq_observed == 1))
  freq_observed <- freq_observed[freq_observed > 0]
  
  list(
    frequencies = freq_observed, n = sum(freq_observed), N = sum(freq_full),
    K_n = sum(freq_observed > 0), K_N = sum(freq_full > 0),
    m1 = sum(freq_observed == 1), percentage = n / N, true_tau1 = true_tau1
  )
}

dataset_creation_zipfH <- function(n, N, zipf_param, H) {
  
  # Old implementation using truncated zipf laws
  points_full <- factor(sample(1:H, N, prob=1/(1:H)^zipf_param, replace=TRUE))
  points_obs <- sample(points_full, n, replace = FALSE)
  
  # Observed frequencies
  freq_full <- as.numeric(table(points_full))
  freq_observed <- as.numeric(table(points_obs))
  
  # True tau1
  true_tau1 <- sum((freq_full == 1) & (freq_observed == 1))
  freq_observed <- freq_observed[freq_observed > 0]
  
  list(
    frequencies = freq_observed, n = sum(freq_observed), N = sum(freq_full),
    K_n = sum(freq_observed > 0), K_N = sum(freq_full > 0),
    m1 = sum(freq_observed == 1), percentage = n / N, true_tau1 = true_tau1
  )
}


# -------------------------------------------
# Scenario 1 - Zipf
# -------------------------------------------

N <- 1000000L # Important to us L, otherwise is not recognized as integer
n <- 100000L

# List of potential parameters
zipf_param_list <- c(1.0526, 1.1765, 1.3333, 1.5385, 1.8182, 2.2222, 2.8571, 4, 6.6667, 20)
zipf_param <- zipf_param_list[4]

set.seed(123)
dataset <- dataset_creation_zipf(n = n, zipf_param = zipf_param,  N = N)

out_PY <- max_EPPF_PY(dataset$frequencies)

kable(data.frame(
  n = dataset$n, N = dataset$N, percentage = dataset$percentage,
  K_n = dataset$K_n, K_N = dataset$K_N, m1 = dataset$m1,
  true_tau1 = dataset$true_tau1,
  zipf_param = dataset$zipf_param,
  K_n_hat = expected_cl_py(dataset$n, out_PY$par[2], out_PY$par[1]),
  m1_hat = expected_m_py(1, dataset$n, out_PY$par[2], out_PY$par[1])
))

# Comparison between M_l and the expected values
M_l <- as.numeric(table(factor(dataset$frequencies, levels = 1:dataset$n)))

# PY comparison
tab <- rbind(PY = expected_m_py(1:15, dataset$n, out_PY$par[2], out_PY$par[1]),
             Data = M_l[1:15])
colnames(tab) <- 1:15
kable(tab, digits=0)

frequency_check_PY(dataset$frequencies)

# PY estimation
tau1_PY <- tau1_py(dataset$m1, dataset$n, out_PY$par[1], out_PY$par[2], dataset$N)
PY_sim <- tau1_py_sim(dataset$frequencies, out_PY$par[1], out_PY$par[2], dataset$N)
PY_lower <- quantile(PY_sim, 0.01 / 2)
PY_upper <- quantile(PY_sim, 1 - 0.01 / 2)

# Dirichlet process estimation
out_DP <- max_EPPF_DP(dataset$frequencies)
tau1_DP <- tau1_dp(dataset$m1, dataset$n, out_DP$par[1], dataset$N)
DP_lower <- qhyper(0.01 / 2, out_DP$par[1] + dataset$n - 1, dataset$N - dataset$n, dataset$m1)
DP_upper <- qhyper(1 - 0.01 / 2, out_DP$par[1] + dataset$n - 1, dataset$N - dataset$n, dataset$m1)

# Summary
kable(data.frame(
  n = dataset$n, N = dataset$N, percentage = round(dataset$percentage, 2),
  m1 = dataset$m1,
  K_n = dataset$K_n,
  true_tau1 = dataset$true_tau1,
  tau1_py = tau1_PY, CI_PY = paste("[", PY_lower, ", ", PY_upper, "]", sep = ""),
  tau1_dp = tau1_DP, CI_DP = paste("[", DP_lower, ", ", DP_upper, "]", sep = "")
))


# -------------------------------------------
# Scenario 2 - Geometric distribution
# -------------------------------------------

N <- 1000000L # Important to us L, otherwise is not recognized as integer
n <- 100000L

set.seed(123)
dataset <- dataset_creation_geom(n = n, N = N, p = 0.01)

out_PY <- max_EPPF_PY(dataset$frequencies)
out_DP <- max_EPPF_DP(dataset$frequencies)

kable(data.frame(
  n = dataset$n, N = dataset$N, percentage = dataset$percentage,
  K_n = dataset$K_n, K_N = dataset$K_N, m1 = dataset$m1,
  true_tau1 = dataset$true_tau1,
  K_n_hat = expected_cl_py(dataset$n, 0, out_DP$par[1]),
  m1_hat = expected_m_py(1, dataset$n, out_PY$par[2], out_PY$par[1])
))

# Comparison between M_l and the expected values
M_l <- as.numeric(table(factor(dataset$frequencies, levels = 1:dataset$n)))

# PY comparison
tab <- rbind(PY = expected_m_py(1:15, dataset$n, out_PY$par[2], out_PY$par[1]),
             Data = M_l[1:15])
colnames(tab) <- 1:15
kable(tab, digits=0)
frequency_check_PY(dataset$frequencies)

# PY estimation
tau1_PY <- tau1_py(dataset$m1, dataset$n, out_PY$par[1], out_PY$par[2], dataset$N)
PY_sim <- tau1_py_sim(dataset$frequencies, out_PY$par[1], out_PY$par[2], dataset$N)
PY_lower <- quantile(PY_sim, 0.01 / 2)
PY_upper <- quantile(PY_sim, 1 - 0.01 / 2)

# Dirichlet process estimation
out_DP <- max_EPPF_DP(dataset$frequencies)
tau1_DP <- tau1_dp(dataset$m1, dataset$n, out_DP$par[1], dataset$N)
DP_lower <- qhyper(0.01 / 2, out_DP$par[1] + dataset$n - 1, dataset$N - dataset$n, dataset$m1)
DP_upper <- qhyper(1 - 0.01 / 2, out_DP$par[1] + dataset$n - 1, dataset$N - dataset$n, dataset$m1)

# Summary
kable(data.frame(
  n = dataset$n, N = dataset$N, percentage = round(dataset$percentage, 2),
  m1 = dataset$m1,
  K_n = dataset$K_n,
  true_tau1 = dataset$true_tau1,
  tau1_py = tau1_PY, CI_PY = paste("[", PY_lower, ", ", PY_upper, "]", sep = ""),
  tau1_dp = tau1_DP, CI_DP = paste("[", DP_lower, ", ", DP_upper, "]", sep = "")
))

# -------------------------------------------
# Scenario 3 - Truncated Zipf
# -------------------------------------------

N <- 1000000L # Important to us L, otherwise is not recognized as integer
n <- 100000L
H <- 10000000

# List of potential parameters
zipf_param_list <- c(1.0526, 1.1765, 1.3333, 1.5385, 1.8182, 2.2222, 2.8571, 4, 6.6667, 20)
zipf_param <- zipf_param_list[2]

set.seed(123)

dataset <- dataset_creation_zipfH(n = n, N = N, zipf_param = zipf_param, H = H)

out_PY <- max_EPPF_PY(dataset$frequencies)

kable(data.frame(
  n = dataset$n, N = dataset$N, percentage = dataset$percentage,
  K_n = dataset$K_n, K_N = dataset$K_N, m1 = dataset$m1,
  true_tau1 = dataset$true_tau1,
  K_n_hat = expected_cl_py(dataset$n, out_PY$par[2], out_PY$par[1]),
  m1_hat = expected_m_py(1, dataset$n, out_PY$par[2], out_PY$par[1])
))

# Comparison between M_l and the expected values
M_l <- as.numeric(table(factor(dataset$frequencies, levels = 1:dataset$n)))

# PY comparison
tab <- rbind(PY = expected_m_py(1:15, dataset$n, out_PY$par[2], out_PY$par[1]),
             Data = M_l[1:15])
colnames(tab) <- 1:15
kable(tab, digits=0)

# Graphical representation of the above table
frequency_check_PY(dataset$frequencies)

# PY estimation
tau1_PY <- tau1_py(dataset$m1, dataset$n, out_PY$par[1], out_PY$par[2], dataset$N)
PY_sim <- tau1_py_sim(dataset$frequencies, out_PY$par[1], out_PY$par[2], dataset$N)
PY_lower <- quantile(PY_sim, 0.01 / 2)
PY_upper <- quantile(PY_sim, 1 - 0.01 / 2)

# Dirichlet process estimation
out_DP <- max_EPPF_DP(dataset$frequencies)
tau1_DP <- tau1_dp(dataset$m1, dataset$n, out_DP$par[1], dataset$N)
DP_lower <- qhyper(0.01 / 2, out_DP$par[1] + dataset$n - 1, dataset$N - dataset$n, dataset$m1)
DP_upper <- qhyper(1 - 0.01 / 2, out_DP$par[1] + dataset$n - 1, dataset$N - dataset$n, dataset$m1)

# Summary
kable(data.frame(
  n = dataset$n, N = dataset$N, percentage = round(dataset$percentage, 2),
  m1 = dataset$m1,
  K_n = dataset$K_n,
  true_tau1 = dataset$true_tau1,
  tau1_py = tau1_PY, CI_PY = paste("[", PY_lower, ", ", PY_upper, "]", sep = ""),
  tau1_dp = tau1_DP, CI_DP = paste("[", DP_lower, ", ", DP_upper, "]", sep = "")
))
