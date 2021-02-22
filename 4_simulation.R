rm(list = ls())

source("2_functions.R")
library(xtable)
library(gridExtra)
library(knitr)
library(ggplot2)
library(dplyr)
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
  points_full <- factor(sample(1:H, N, prob = probs, replace = TRUE))
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

N <- 1000000L # Use L, otherwise is not recognized as integer
n <- 100000L

# List of potential parameters
zipf_param_list <- c(1.0526, 1.1765, 1.3333, 1.5385, 1.8182, 2.2222, 2.8571, 4)
zipf_param <- zipf_param_list[1]

set.seed(123)
dataset <- dataset_creation_zipf(n = n, zipf_param = zipf_param, N = N)

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
tab <- rbind(
  PY = expected_m_py(1:15, dataset$n, out_PY$par[2], out_PY$par[1]),
  Data = M_l[1:15]
)
colnames(tab) <- 1:15
kable(tab, digits = 0)

frequency_check_PY(dataset$frequencies)

# PY estimation
tau1_PY <- tau1_py(dataset$m1, dataset$n, out_PY$par[1], out_PY$par[2], dataset$N)
PY_sim <- tau1_py_sim(dataset$frequencies, out_PY$par[1], out_PY$par[2], dataset$N)
PY_lower <- quantile(PY_sim, 0.01 / 2)
PY_upper <- quantile(PY_sim, 1 - 0.01 / 2)

# Dirichlet process estimation
out_DP <- max_EPPF_DP(dataset$frequencies)
tau1_DP <- tau1_dp(dataset$m1, dataset$n, out_DP$par[1], dataset$N)
DP_lower <- qhyper2(0.01 / 2, out_DP$par[1] + dataset$n - 1, dataset$N - dataset$n, dataset$m1)
DP_upper <- qhyper2(1 - 0.01 / 2, out_DP$par[1] + dataset$n - 1, dataset$N - dataset$n, dataset$m1)

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
  "PY" = round(tau1_PY, 0),
  "CI PY" = paste("[", round(PY_lower, 0), ", ", round(PY_upper, 0), "]", sep = ""),
  "DP" = round(tau1_DP, 0),
  "CI DP" = paste("[", round(DP_lower, 0), ", ", round(DP_upper, 0), "]", sep = ""),
  "Bethlehem" = round(tau1_bet, 0),
  "Skinner" = round(tau1_skin, 0)
))

# Plot estimates and confidence intervals
type <- c("PY", "DP", "B", "S")
estimates <- c(tau1_PY, tau1_DP, tau1_bet, tau1_skin)
df <- data.frame(
  type = factor(type, levels = type[order(estimates)]), estim = estimates,
  lower_CI = c(PY_lower, DP_lower, NA, NA), upper_CI = c(PY_upper, DP_upper, NA, NA)
)

p <- ggplot(df, aes(type, estim, color = type))
p + geom_pointrange(aes(ymin = lower_CI, ymax = upper_CI)) +
  theme(legend.position = "none") + xlab("") + ylab("estimate") + ggtitle(paste0("Zipf ", zipf_param)) + 
  theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + geom_hline(yintercept = dataset$true_tau1) + 
  scale_y_log10()


# Compare different Zipf parameters

N <- 1000000L # Important to us L, otherwise is not recognized as integer
n <- 100000L
p <- list()
check <- list()
true_tau1 = m1 <- c()
tau1_bet = tau1_skin = tau1_PY = tau1_DP <- c()
PY_lower = PY_upper = DP_lower = DP_upper <- c()
theta_MLE_PY = alpha_MLE_PY = theta_MLE_DP <- c()
zipf_param_list <- c(1.25, 1.5, 1.75, 2)
for (i in 1:length(zipf_param_list)) {
  zipf_param <- zipf_param_list[i]
  set.seed(123)
  dataset <- dataset_creation_zipf(n = n, zipf_param = zipf_param, N = N)
  K_n <- dataset$K_n
  true_tau1[i] <- dataset$true_tau1
  m1[i] <- dataset$m1

  # PY estimation
  out_PY <- max_EPPF_PY(dataset$frequencies)
  tau1_PY[i] <- tau1_py(dataset$m1, dataset$n, out_PY$par[1], out_PY$par[2], dataset$N)
  PY_sim <- tau1_py_sim(dataset$frequencies, out_PY$par[1], out_PY$par[2], dataset$N)
  PY_lower[i] <- quantile(PY_sim, 0.01 / 2)
  PY_upper[i] <- quantile(PY_sim, 1 - 0.01 / 2)
  theta_MLE_PY[i] <- out_PY$par[1]
  alpha_MLE_PY[i] <- out_PY$par[2]

  # Dirichlet process estimation
  out_DP <- max_EPPF_DP(dataset$frequencies)
  tau1_DP[i] <- tau1_dp(dataset$m1, dataset$n, out_DP$par[1], dataset$N)
  DP_lower[i] <- qhyper2(0.01 / 2, out_DP$par[1] + dataset$n - 1, dataset$N - dataset$n, dataset$m1)
  DP_upper[i] <- qhyper2(1 - 0.01 / 2, out_DP$par[1] + dataset$n - 1, dataset$N - dataset$n, dataset$m1)
  theta_MLE_DP[i] <- out_DP$par
  
  # Comparison between M_l and the expected values
  M_l <- as.numeric(table(factor(dataset$frequencies, levels = 1:dataset$n)))
  
  # PY comparison
  tab <- rbind(
    PY = expected_m_py(1:15, dataset$n, out_PY$par[2], out_PY$par[1]),
    Data = M_l[1:15]
  )
  colnames(tab) <- 1:15
  kable(tab, digits = 0)
  check[[i]] <- frequency_check_PY(dataset$frequencies)
  
  # Bethlehem and Skinner estimators
  estim <- tau1_bs(dataset$frequencies, dataset$N)
  tau1_bet[i] <- estim[1]
  tau1_skin[i] <- estim[2]

  # Plot estimates and confidence intervals
  type <- c("PY", "DP", "B", "S")
  estimates <- c(tau1_PY[i], tau1_DP[i], tau1_bet[i], tau1_skin[i])
  df <- data.frame(
    type = factor(type), estim = estimates,
    lower_CI = c(PY_lower[i], DP_lower[i], NA, NA),
    upper_CI = c(PY_upper[i], DP_upper[i], NA, NA)
  )
  p[[i]] <- ggplot(df, aes(type, estim))
  p[[i]] <- p[[i]] + geom_pointrange(aes(ymin = lower_CI, ymax = upper_CI)) + theme_bw() +
    theme(legend.position = "none") + xlab("") + ylab(expression(tau[1])) + 
    ggtitle(paste0("Zipf ", round(zipf_param, 2))) + 
    theme(plot.title = element_text(hjust = 0.5, size=30)) + 
    geom_hline(yintercept = dataset$true_tau1, linetype = "dotted")
}

# plot estimates
a <- do.call(grid.arrange, c(p, ncol = 4))
ggsave(a, file = "zipf.eps", height = 5, width = 7 * 2.5, device = "eps")

# Summary tables

df_zipf <- data.frame(
  D = round(zipf_param_list, 2),
  C = m1, E = true_tau1, H = paste(as.integer(tau1_PY), " in [", round(PY_lower, 0), ", ", round(PY_upper, 0), "]", sep = ""),
  L = paste(as.integer(tau1_DP), " in [", round(DP_lower, 0), ", ", round(DP_upper, 0), "]", sep = ""),
  M = as.integer(tau1_bet), N = as.integer(tau1_skin)
)
colnames(df_zipf) <- c(
  "Zipf parameter", "m_1", "tau_1", "Pitman-Yor",
  "Dirichlet Process", "Bethlehem", "Skinner"
)
knitr::kable(df_zipf)
xtable(df_zipf)

df_zipf_param <- data.frame(
  A = round(zipf_param_list, 2)[1:2], O = round(theta_MLE_PY[1:2], 2), P = round(alpha_MLE_PY[1:2], 2), Q = round(theta_MLE_DP[1:2], 2),
  A = round(zipf_param_list, 2)[3:4], O = round(theta_MLE_PY[3:4], 2), P = round(alpha_MLE_PY[3:4], 2), Q = round(theta_MLE_DP[3:4], 2)
)
colnames(df_zipf_param) <- c(
  "Zipf param", "theta PY", "alpha PY", "theta DP",
  "Zipf param", "theta PY", "alpha PY", "theta DP"
)
knitr::kable(df_zipf_param)
xtable(df_zipf_param)

# -------------------------------------------
# Scenario 2 - Geometric distribution
# -------------------------------------------

N <- 1000000L #  L is important, otherwise is not recognized as integer
n <- 100000L

set.seed(123)
prob <- 0.0005
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
tab <- rbind(
  PY = expected_m_py(1:15, dataset$n, out_PY$par[2], out_PY$par[1]),
  Data = M_l[1:15]
)
colnames(tab) <- 1:15
kable(tab, digits = 0)
frequency_check_PY(dataset$frequencies)

# PY estimation
tau1_PY <- tau1_py(dataset$m1, dataset$n, out_PY$par[1], out_PY$par[2], dataset$N)
PY_sim <- tau1_py_sim(dataset$frequencies, out_PY$par[1], out_PY$par[2], dataset$N)
PY_lower <- quantile(PY_sim, 0.01 / 2)
PY_upper <- quantile(PY_sim, 1 - 0.01 / 2)

# Dirichlet process estimation
out_DP <- max_EPPF_DP(dataset$frequencies)
tau1_DP <- tau1_dp(dataset$m1, dataset$n, out_DP$par[1], dataset$N)
DP_lower <- qhyper2(0.01 / 2, out_DP$par[1] + dataset$n - 1, dataset$N - dataset$n, dataset$m1)
DP_upper <- qhyper2(1 - 0.01 / 2, out_DP$par[1] + dataset$n - 1, dataset$N - dataset$n, dataset$m1)

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
type <- c("PY", "DP", "B", "S")
estimates <- c(tau1_PY, tau1_DP, tau1_bet, tau1_skin)
df <- data.frame(
  type = factor(type, levels = type[order(estimates)]), estim = estimates,
  lower_CI = c(PY_lower, DP_lower, NA, NA), upper_CI = c(PY_upper, DP_upper, NA, NA)
)

p <- ggplot(df, aes(type, estim, color = type))
p + geom_pointrange(aes(ymin = lower_CI, ymax = upper_CI)) + theme_bw() +
  theme(legend.position = "none") + xlab("") + ylab(expression(tau[1])) + ggtitle(paste0("Geometric parameter: ", round(prob, 6))) + theme(plot.title = element_text(hjust = 0.5)) + geom_hline(yintercept = dataset$true_tau1, linetype = "dotted")

# plot together
N <- 1000000L # Important to us L, otherwise is not recognized as integer
n <- 100000L

check <- list()
p_g <- list()
true_tau1_g = m1_g <- c()
tau1_bet_g = tau1_skin_g = tau1_PY_g = tau1_DP_g <- c()
PY_lower_g = PY_upper_g = DP_lower_g = DP_upper_g <- c()
theta_MLE_PY_g = alpha_MLE_PY_g = theta_MLE_DP_g <- c()
p_g <- list()
geom_param_list <- c(0.0001, 0.001)
for (i in 1:length(geom_param_list)) {
  geom_param <- geom_param_list[i]
  set.seed(123)
  dataset <- dataset_creation_geom(n = n, N = N, p = geom_param)
  K_n <- dataset$K_n
  true_tau1_g[i] <- dataset$true_tau1
  m1_g[i] <- dataset$m1

  # Comparison between M_l and the expected values
  M_l <- as.numeric(table(factor(dataset$frequencies, levels = 1:dataset$n)))

  # PY comparison
  tab <- rbind(
    PY = expected_m_py(1:15, dataset$n, out_PY$par[2], out_PY$par[1]),
    Data = M_l[1:15]
  )
  colnames(tab) <- 1:15
  kable(tab, digits = 0)
  check[[i]] <- frequency_check_PY(dataset$frequencies)

  # PY estimation
  out_PY <- max_EPPF_PY(dataset$frequencies)
  tau1_PY_g[i] <- tau1_py(dataset$m1, dataset$n, out_PY$par[1], out_PY$par[2], dataset$N)
  PY_sim <- tau1_py_sim(dataset$frequencies, out_PY$par[1], out_PY$par[2], dataset$N)
  PY_lower_g[i] <- quantile(PY_sim, 0.01 / 2)
  PY_upper_g[i] <- quantile(PY_sim, 1 - 0.01 / 2)
  theta_MLE_PY_g[i] <- out_PY$par[1]
  alpha_MLE_PY_g[i] <- out_PY$par[2]
  
  # Dirichlet process estimation
  out_DP <- max_EPPF_DP(dataset$frequencies)
  tau1_DP_g[i] <- tau1_dp(dataset$m1, dataset$n, out_DP$par[1], dataset$N)
  DP_lower_g[i] <- qhyper2(0.01 / 2, out_DP$par[1] + dataset$n - 1, dataset$N - dataset$n, dataset$m1)
  DP_upper_g[i] <- qhyper2(1 - 0.01 / 2, out_DP$par[1] + dataset$n - 1, dataset$N - dataset$n, dataset$m1)
  theta_MLE_DP_g[i] <- out_DP$par
  
  # Bethlehem and Skinner estimators
  estim <- tau1_bs(dataset$frequencies, dataset$N)
  tau1_bet_g[i] <- estim[1]
  tau1_skin_g[i] <- estim[2]
  
  # Plot estimates and confidence intervals
  type <- c("PY", "DP", "B", "S")
  estimates <- c(tau1_PY_g[i], tau1_DP_g[i], tau1_bet_g[i], tau1_skin_g[i])
  df <- data.frame(
    type = factor(type), estim = estimates,
    lower_CI = c(PY_lower_g[i], DP_lower_g[i], NA, NA),
    upper_CI = c(PY_upper_g[i], DP_upper_g[i], NA, NA)
  )
  p_g[[i]] <- ggplot(df, aes(type, estim))
  p_g[[i]] <- p_g[[i]] + geom_pointrange(aes(ymin = lower_CI, ymax = upper_CI)) + theme_bw() +
    theme(legend.position = "none") + xlab("") + ylab(expression(tau[1])) + 
    ggtitle(paste0("Geometric ", geom_param)) + theme(plot.title = element_text(hjust = 0.5)) + 
    geom_hline(yintercept = dataset$true_tau1, linetype = "dotted")
}

# Plot estimates
b <- do.call(grid.arrange, c(p_g, ncol = 2))
ggsave(b, file = "geom.eps", height = 5, width = 7 * 2.5, device = "eps")

# Summary tables
df_geom <- data.frame(
  D = geom_param_list,
  C = m1_g, E = true_tau1_g,
  H = paste(as.integer(tau1_PY_g), " in [", round(PY_lower_g, 0), ", ", round(PY_upper_g, 0), "]", sep = ""),
  L = paste( as.integer(tau1_DP_g), " in [", round(DP_lower_g, 0), ", ", round(DP_upper_g, 0), "]", sep = ""),
  M = as.integer(tau1_bet_g), N = as.integer(tau1_skin_g)
)
colnames(df_geom) <- c(
  "Geometric parameter", "m_1", "tau_1", "Pitman-Yor",
  "Dirichlet Process", "Bethlehem", "Skinner"
)
knitr::kable(df_geom)
xtable(df_geom)

df_geom_param <- data.frame(
  A = geom_param_list, O = round(theta_MLE_PY_g, 2), P = round(alpha_MLE_PY_g, 2), Q = round(theta_MLE_DP_g, 2)
)
colnames(df_geom_param) <- c(
  "Geometric parameter", "theta PY", "alpha PY", "theta DP"
)
knitr::kable(df_geom_param)
xtable(df_geom_param)


# -------------------------------------------
# Scenario 3 - Custom probabilities
# -------------------------------------------

N <- 1000000L # Important to us L, otherwise is not recognized as integer
n <- 100000L
H <- 100000

set.seed(123)

# Define the probabilities of the cells directly
probs <- c(rep(50, 1000), rep(1, 10^4), rep(0.1, 10^4))
probs <- probs / sum(probs)

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
tab <- rbind(
  PY = expected_m_py(1:15, dataset$n, out_PY$par[2], out_PY$par[1]),
  Data = M_l[1:15]
)
colnames(tab) <- 1:15
kable(tab, digits = 0)

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
DP_lower <- qhyper2(0.01 / 2, out_DP$par[1] + dataset$n - 1, dataset$N - dataset$n, dataset$m1)
DP_upper <- qhyper2(1 - 0.01 / 2, out_DP$par[1] + dataset$n - 1, dataset$N - dataset$n, dataset$m1)

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
m1 <- dataset$m1
true_tau1 <- dataset$true_tau1
CI_PY <- paste("[", PY_lower, ", ", PY_upper, "]", sep = "")
CI_DP <- paste("[", DP_lower, ", ", DP_upper, "]", sep = "")

df_custom <- data.frame(A = m1, B = true_tau1, C = tau1_PY, D = CI_PY, E = tau1_DP, f = CI_DP, G = tau1_bet, H = tau1_skin)
colnames(df_custom) <- c("m1", "tau_1", "tau_1^PY", "95 CI PY", "tau_1^DP", "95 CI DP", "tau_1^B", "tau_1^S")
xtable(df_custom)

# Plot estimates and confidence intervals
type <- c("PY", "DP", "B", "S")
estimates <- c(tau1_PY, tau1_DP, tau1_bet, tau1_skin)
df <- data.frame(
  type = factor(type, levels = type[order(estimates)]), estim = estimates,
  lower_CI = c(PY_lower, DP_lower, NA, NA), upper_CI = c(PY_upper, DP_upper, NA, NA)
)

p <- ggplot(df, aes(type, color = type, estim))
p <- p + geom_pointrange(aes(ymin = lower_CI, ymax = upper_CI), size = 1.2) + theme_bw(base_size = 18) +
  theme(legend.position = "none") + xlab("") + ylab(expression(tau[1])) + ggtitle("Custom probabilities") + 
  theme(plot.title = element_text(hjust = 0.5)) + geom_hline(yintercept = dataset$true_tau1, linetype = "dotted")
p
ggsave(p, height = 5, width = 7 * 2.5, file = "custom.eps", device = "eps")
