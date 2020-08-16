rm(list = ls())

source("2_functions.R")
library(knitr)
library(ggplot2)
library(reticulate) # Library that interface R to python
use_condaenv() # Use python anaconda

random <- import("numpy.random") # Import python libraries



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

dataset_creation_probs <- function(n, N, probs) {
  
  # Old implementation using truncated zipf laws
  H <- length(probs)
  points_full <- factor(sample(1:H, N, prob=probs, replace=TRUE))
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

# Bethlehem and Skinner estimators
estim <- tau1_bs(dataset$frequencies, dataset$N)
tau1_bet <- estim[1]
tau1_skin <- estim[2]


# Summary
kable(data.frame(
  n = dataset$n, N = dataset$N, percentage = round(dataset$percentage, 2),
  m1 = dataset$m1,
  K_n = dataset$K_n,
  true_tau1 = dataset$true_tau1,
  tau1_py = tau1_PY, CI_PY = paste("[", PY_lower, ", ", PY_upper, "]", sep = ""),
  tau1_dp = tau1_DP, CI_DP = paste("[", DP_lower, ", ", DP_upper, "]", sep = ""),
  tau1_bet = tau1_bet, tau1_skin = tau1_skin
))

# Summary rounded
kable(data.frame(
  N = dataset$N, percentage = round(dataset$percentage, 2),
  m1 = dataset$m1,
  "# classes" = dataset$K_n,
  "true tau1" = dataset$true_tau1,
  "PY" = round(tau1_PY,0), "CI PY" = paste("[", round(PY_lower,0), ", ", round(PY_upper,0), "]", sep = ""),
  "DP" = round(tau1_DP,0), "CI DP" = paste("[", round(DP_lower,0), ", ", round(DP_upper,0), "]", sep = ""),
  "Bethlehem" = round(tau1_bet,0), "Skinner" = round(tau1_skin,0)
))

# Plot estimates and confidence intervals
type = c('PY', 'DP', 'B', 'S')
estimates = c(tau1_PY, tau1_DP, tau1_bet, tau1_skin)
df <- data.frame(type = factor(type, levels = type[order(estimates)]), estim = estimates,
                 lower_CI = c(PY_lower, DP_lower, NA, NA), upper_CI = c(PY_upper, DP_upper, NA, NA))

p <- ggplot(df, aes(type, estim, color=type))
p + geom_pointrange(aes(ymin = lower_CI, ymax = upper_CI)) + 
  theme(legend.position = "none") + xlab('') + ylab('estimate') + ggtitle(paste0('Zipf ', zipf_param)) +
  theme(plot.title = element_text(hjust = 0.5)) + geom_hline(yintercept=dataset$true_tau1)

# -------------------------------------------
# Scenario 2 - Geometric distribution
# -------------------------------------------

N <- 1000000L # Important to us L, otherwise is not recognized as integer
n <- 100000L

set.seed(123)
prob = 0.005
dataset <- dataset_creation_geom(n = n, N = N, p = prob)

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

# Bethlehem and Skinner estimators
estim <- tau1_bs(dataset$frequencies, dataset$N)
tau1_bet <- estim[1]
tau1_skin <- estim[2]

# Summary
kable(data.frame(
  n = dataset$n, N = dataset$N, percentage = round(dataset$percentage, 2),
  m1 = dataset$m1,
  K_n = dataset$K_n,
  true_tau1 = dataset$true_tau1,
  tau1_py = tau1_PY, CI_PY = paste("[", PY_lower, ", ", PY_upper, "]", sep = ""),
  tau1_dp = tau1_DP, CI_DP = paste("[", DP_lower, ", ", DP_upper, "]", sep = ""),
  tau1_bet = tau1_bet, tau1_skin = tau1_skin
))

# Plot estimates and confidence intervals
type = c('PY', 'DP', 'B', 'S')
estimates = c(tau1_PY, tau1_DP, tau1_bet, tau1_skin)
df <- data.frame(type = factor(type, levels = type[order(estimates)]), estim = estimates,
                 lower_CI = c(PY_lower, DP_lower, NA, NA), upper_CI = c(PY_upper, DP_upper, NA, NA))

p <- ggplot(df, aes(type, estim, color=type))
p + geom_pointrange(aes(ymin = lower_CI, ymax = upper_CI)) + 
  theme(legend.position = "none") + xlab('') + ylab('estimate') + ggtitle(paste0('Geometric ', prob)) +
  theme(plot.title = element_text(hjust = 0.5)) + geom_hline(yintercept=dataset$true_tau1)

# -------------------------------------------
# Scenario 3 - Custom probabilities
# -------------------------------------------

N <- 1000000L # Important to us L, otherwise is not recognized as integer
n <- 100000L
H <- 100000

set.seed(123)

# Define the probabilities of the cells directly
probs <- c(rep(50, 1000), rep(1, 10^4), rep(0.1,10^4)); probs <- probs/sum(probs)

# Generating the dataset
dataset <- dataset_creation_probs(n = n, N = N, probs = probs)

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

# Bethlehem and Skinner estimators
estim <- tau1_bs(dataset$frequencies, dataset$N)
tau1_bet <- estim[1]
tau1_skin <- estim[2]

# Summary
kable(data.frame(
  n = dataset$n, N = dataset$N, percentage = round(dataset$percentage, 2),
  m1 = dataset$m1,
  K_n = dataset$K_n,
  true_tau1 = dataset$true_tau1,
  tau1_py = tau1_PY, CI_PY = paste("[", PY_lower, ", ", PY_upper, "]", sep = ""),
  tau1_dp = tau1_DP, CI_DP = paste("[", DP_lower, ", ", DP_upper, "]", sep = ""),
  tau1_bet = tau1_bet, tau1_skin = tau1_skin
))

# Plot estimates and confidence intervals
type = c('PY', 'DP', 'B', 'S')
estimates = c(tau1_PY, tau1_DP, tau1_bet, tau1_skin)
df <- data.frame(type = factor(type, levels = type[order(estimates)]), estim = estimates,
                 lower_CI = c(PY_lower, DP_lower, NA, NA), upper_CI = c(PY_upper, DP_upper, NA, NA))

p <- ggplot(df, aes(type, estim, color=type))
p + geom_pointrange(aes(ymin = lower_CI, ymax = upper_CI)) + geom_hline(yintercept = dataset$true_tau1)
  theme(legend.position = "none") + xlab('') + ylab('estimate') + ggtitle(paste0('Custom probabilities')) +
  theme(plot.title = element_text(hjust = 0.5))

