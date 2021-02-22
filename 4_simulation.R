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
  theme(legend.position = "none") + xlab("") + ylab("estimate") + ggtitle(paste0("Zipf ", zipf_param)) + theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + geom_hline(yintercept = dataset$true_tau1) + scale_y_log10()


# plot together
N <- 1000000L # Important to us L, otherwise is not recognized as integer
n <- 100000L

p <- list()
check <- list()
K_n <- c()
tau1_bet <- c()
tau1_skin <- c()
tau1_PY <- c()
tau1_DP <- c()
PY_lower <- c()
PY_upper <- c()
DP_lower <- c()
DP_upper <- c()
true_tau1 <- c()
m1 <- c()
out_PY <- list()
out_DP <- list()
dataset <- list()
zipf_param_list <- c(1.0526, 1.1765, 1.3333, 1.5385, 1.8182, 2.2222, 2.8571, 4)
for (i in 1:length(zipf_param_list)) {
  zipf_param <- zipf_param_list[i]
  set.seed(123)
  dataset[[i]] <- dataset_creation_zipf(n = n, zipf_param = zipf_param, N = N)
  K_n[i] <- dataset[[i]]$K_n
  true_tau1[i] <- dataset[[i]]$true_tau1
  m1[i] <- dataset[[i]]$m1

  # # Comparison between M_l and the expected values
  # M_l <- as.numeric(table(factor(dataset[[i]]$frequencies, levels = 1:dataset[[i]]$n)))

  # # PY comparison
  # tab <- rbind(
  #   PY = expected_m_py(1:15, dataset[[i]]$n, out_PY$par[2], out_PY$par[1]),
  #   Data = M_l[1:15]
  # )
  # colnames(tab) <- 1:15
  # kable(tab, digits = 0)
  # check[[i]] <- frequency_check_PY(dataset[[i]]$frequencies)

  # PY estimation
  out_PY[[i]] <- max_EPPF_PY(dataset[[i]]$frequencies)
  tau1_PY[i] <- tau1_py(dataset[[i]]$m1, dataset[[i]]$n, out_PY[[i]]$par[1], out_PY[[i]]$par[2], dataset[[i]]$N)
  PY_sim <- tau1_py_sim(dataset[[i]]$frequencies, out_PY[[i]]$par[1], out_PY[[i]]$par[2], dataset[[i]]$N)
  PY_lower[i] <- quantile(PY_sim, 0.01 / 2)
  PY_upper[i] <- quantile(PY_sim, 1 - 0.01 / 2)

  # Dirichlet process estimation
  out_DP[[i]] <- max_EPPF_DP(dataset[[i]]$frequencies)
  tau1_DP[i] <- tau1_dp(dataset[[i]]$m1, dataset[[i]]$n, out_DP[[i]]$par[1], dataset[[i]]$N)
  DP_lower[i] <- qhyper(0.01 / 2, out_DP[[i]]$par[1] + dataset[[i]]$n - 1, dataset[[i]]$N - dataset[[i]]$n, dataset[[i]]$m1)
  DP_upper[i] <- qhyper(1 - 0.01 / 2, out_DP[[i]]$par[1] + dataset[[i]]$n - 1, dataset[[i]]$N - dataset[[i]]$n, dataset[[i]]$m1)

  # Bethlehem and Skinner estimators
  estim <- tau1_bs(dataset[[i]]$frequencies, dataset[[i]]$N)
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
  p[[i]] <- ggplot(df, aes(type, estim, color = type))
  p[[i]] <- p[[i]] + geom_pointrange(aes(ymin = lower_CI, ymax = upper_CI)) + theme_bw() +
    theme(legend.position = "none") + xlab("") + ylab(expression(tau[1])) + ggtitle(paste0("Zipf parameter: ", round(zipf_param, 2))) + theme(plot.title = element_text(hjust = 0.5)) + geom_hline(yintercept = dataset$true_tau1, linetype = "dotted")
}

# plot estimates
a <- do.call(grid.arrange, c(p, ncol = 4))
ggsave(a, file = "zipf.eps", device = "eps")
do.call(grid.arrange, c(check, ncol = 4))

# Summary tables
theta_MLE_PY <- c()
alpha_MLE_PY <- c()
theta_MLE_DP <- c()
for (i in 1:length(zipf_param_list)) {
  theta_MLE_PY[i] <- out_PY[[i]]$par[1]
  alpha_MLE_PY[i] <- out_PY[[i]]$par[2]
  theta_MLE_DP[i] <- out_DP[[i]]$par
}
df_zipf <- data.frame(
  D = round(zipf_param_list, 2),
  C = m1, E = true_tau1, G = as.integer(tau1_PY), H = paste("[", round(PY_lower, 0), ", ", round(PY_upper, 0), "]", sep = ""),
  I = as.integer(tau1_DP), L = paste("[", round(DP_lower, 0), ", ", round(DP_upper, 0), "]", sep = ""),
  M = as.integer(tau1_bet), N = as.integer(tau1_skin)
)
colnames(df_zipf) <- c(
  "Zipf parameter", "$m_1$", "$\\tau_1$", "$\\tau_1^{PY}$",
  "CI PY", "$\\tau_1^{DP}$", "CI DP", "$\\tau_1^B$", "$\\tau_1^S$"
)
knitr::kable(df_zipf, col.names = c(
  "Zipf parameter", "m_1", "tau_1", "tau_1 PY",
  "CI PY", "tau_1 DP", "CI DP", "tau_1 B", "tau_1 S"
))

df_zipf_param <- data.frame(
  A = round(zipf_param_list, 2)[1:4], O = theta_MLE_PY[1:4], P = alpha_MLE_PY[1:4], Q = theta_MLE_DP[1:4],
  A = round(zipf_param_list, 2)[5:8], O = theta_MLE_PY[5:8], P = alpha_MLE_PY[5:8], Q = theta_MLE_DP[5:8]
)
colnames(df_zipf_param) <- c(
  "zipf param", "theta PY param", "alpha PY param", "theta DP param",
  "zipf param", "theta PY param", "alpha PY param", "theta DP param"
)

xtable(df_zipf)
xtable(df_zipf_param)


# # ----------
# # Other implementation for minimax estimation, and comparison with MLE
# # --------
#
# alpha_hat <- c()
# theta_hat <- c()
# for(i in 1:length(zipf_param_list)){
#   freq = dataset[[i]]$frequencies
#   m <- rep(0, dataset[[i]]$n)
#   a = data.frame(freq) %>% group_by(freq) %>% summarise(count=n())
#   for(j in 1:dataset[[i]]$n){
#     if(j %in% a$freq) m[j] = a$count[match(j ,a$freq)]
#   }
#   ind_max = tail(which(m!=0),1)
#   cumsum_m = cumsum(rev(m[1:ind_max]))
#   alpha_hat[i] = PY_alpha(freq, ind_max, cumsum_m, dataset[[i]]$K_n)
#   theta_hat[i] = max_EPPF_PY_theta(dataset[[i]]$frequencies, alpha_hat[i])$par
#   print(i)
# }
#
# alpha_MLE <- c()
# for(i in 1:length(zipf_param_list)) alpha_MLE[i] = out_PY[[i]]$par[2]
# theta_MLE <- c()
# for(i in 1:length(zipf_param_list)) theta_MLE[i] = out_PY[[i]]$par[1]
# logEPPF_PY_minimax <- c()
# for(i in 1:length(zipf_param_list)) logEPPF_PY_minimax[i] = logEPPF_PY(theta_hat[i], alpha_hat[i], dataset[[i]]$frequencies)
# logEPPF_PY_MLE <- c()
# for(i in 1:length(zipf_param_list)) logEPPF_PY_MLE[i] = logEPPF_PY(theta_MLE[i], alpha_MLE[i], dataset[[i]]$frequencies)
# df1 <- data.frame(a=alpha_hat, b=alpha_MLE, c=theta_hat, d=theta_MLE, e=logEPPF_PY_minimax, f=logEPPF_PY_MLE)
# knitr::kable(df1, col.names =c("alpha minimax", "alpha MLE", "theta minimax",
#                                 "theta MLE", "log EPPF minimax", "log EPPF MLE"))


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
p <- list()
K_n_g <- c()
tau1_bet_g <- c()
tau1_skin_g <- c()
tau1_PY_g <- c()
tau1_DP_g <- c()
PY_lower_g <- c()
PY_upper_g <- c()
DP_lower_g <- c()
DP_upper_g <- c()
true_tau1_g <- c()
m1_g <- c()
out_PY_g <- list()
out_DP_g <- list()
dataset_g <- list()
geom_param_list <- c(0.0001, 0.001, 0.005, 0.05)
for (i in 1:length(geom_param_list)) {
  geom_param <- geom_param_list[i]
  set.seed(123)
  dataset_g[[i]] <- dataset_creation_geom(n = n, N = N, p = geom_param)
  K_n_g[i] <- dataset_g[[i]]$K_n
  true_tau1_g[i] <- dataset_g[[i]]$true_tau1
  m1_g[i] <- dataset_g[[i]]$m1

  # Comparison between M_l and the expected values
  # M_l <- as.numeric(table(factor(dataset$frequencies, levels = 1:dataset$n)))

  # # PY comparison
  # tab <- rbind(
  #   PY = expected_m_py(1:15, dataset$n, out_PY$par[2], out_PY$par[1]),
  #   Data = M_l[1:15]
  # )
  # colnames(tab) <- 1:15
  # kable(tab, digits = 0)
  # check[[i]] <- frequency_check_PY(dataset$frequencies)

  # PY estimation
  out_PY_g[[i]] <- max_EPPF_PY(dataset_g[[i]]$frequencies)
  tau1_PY_g[i] <- tau1_py(dataset_g[[i]]$m1, dataset_g[[i]]$n, out_PY_g[[i]]$par[1], out_PY_g[[i]]$par[2], dataset_g[[i]]$N)
  PY_sim <- tau1_py_sim(dataset_g[[i]]$frequencies, out_PY_g[[i]]$par[1], out_PY_g[[i]]$par[2], dataset_g[[i]]$N)
  PY_lower_g[i] <- quantile(PY_sim, 0.01 / 2)
  PY_upper_g[i] <- quantile(PY_sim, 1 - 0.01 / 2)

  # Dirichlet process estimation
  out_DP_g[[i]] <- max_EPPF_DP(dataset_g[[i]]$frequencies)
  tau1_DP_g[i] <- tau1_dp(dataset_g[[i]]$m1, dataset_g[[i]]$n, out_DP_g[[i]]$par[1], dataset_g[[i]]$N)
  DP_lower_g[i] <- qhyper(0.01 / 2, out_DP_g[[i]]$par[1] + dataset_g[[i]]$n - 1, dataset_g[[i]]$N - dataset_g[[i]]$n, dataset_g[[i]]$m1)
  DP_upper_g[i] <- qhyper(1 - 0.01 / 2, out_DP_g[[i]]$par[1] + dataset_g[[i]]$n - 1, dataset_g[[i]]$N - dataset_g[[i]]$n, dataset_g[[i]]$m1)

  # Bethlehem and Skinner estimators
  estim <- tau1_bs(dataset_g[[i]]$frequencies, dataset_g[[i]]$N)
  tau1_bet_g[i] <- estim[1]
  tau1_skin_g[i] <- estim[2]

  # Plot estimates and confidence intervals
  type <- c("PY", "DP", "B", "S")
  estimates <- c(tau1_PY, tau1_DP, tau1_bet, tau1_skin)
  df <- data.frame(
    type = factor(type), estim = estimates,
    lower_CI = c(PY_lower_g[i], DP_lower_g[i], NA, NA), upper_CI = c(PY_upper_g[i], DP_upper_g[i], NA, NA)
  )

  p_g[[i]] <- ggplot(df, aes(type, estim, color = type))
  p_g[[i]] <- p_g[[i]] + geom_pointrange(aes(ymin = lower_CI, ymax = upper_CI)) + theme_bw() + theme(legend.position = "none") + xlab("") + ylab(expression(tau[1])) + ggtitle(paste0("Geometric parameter: ", round(geom_param, 7))) + theme(plot.title = element_text(hjust = 0.5)) + geom_hline(yintercept = dataset$true_tau1, linetype = "dotted")
}

# Plot estimates
b <- do.call(grid.arrange, c(p_g, ncol = 2))
ggsave(b, file = "geom.eps", device = "eps")
do.call(grid.arrange, c(check, ncol = 2))

# Summary tables
theta_MLE_PY_g <- c()
alpha_MLE_PY_g <- c()
theta_MLE_DP_g <- c()
for (i in 1:length(geom_param_list)) {
  theta_MLE_PY_g[i] <- out_PY_g[[i]]$par[1]
  alpha_MLE_PY_g[i] <- out_PY_g[[i]]$par[2]
  theta_MLE_DP_g[i] <- out_DP_g[[i]]$par
}
df_geom <- data.frame(
  D = geom_param_list,
  C = m1_g, E = true_tau1_g, G = as.integer(tau1_PY_g),
  H = paste("[", round(PY_lower_g, 0), ", ", round(PY_upper_g, 0), "]", sep = ""),
  I = as.integer(tau1_DP_g),
  L = paste("[", round(DP_lower_g, 0), ", ", round(DP_upper_g, 0), "]", sep = ""),
  M = as.integer(tau1_bet_g), N = as.integer(tau1_skin_g)
)
colnames(df_geom) <- c(
  "pi", "$m_1$", "$\\tau_1$", "$\\tau_1^{PY}$",
  "C.I. PY", "$\\tau_1^{DP}$", "C.I. DP", "$\\tau_1^B$", "$\\tau_1^S$"
)
knitr::kable(df_geom, col.names = c(
  "pi", "m_1", "tau_1", "tau_1 PY",
  "C.I. PY", "tau_1 DP", "C.I. DP", "tau_1 B", "tau_1 S"
))

df_geom_param <- data.frame(
  A = round(geom_param_list, 2)[1:2], O = theta_MLE_PY_g[1:2], P = alpha_MLE_PY_g[1:2], Q = theta_MLE_DP_g[1:2],
  A = round(geom_param_list, 2)[3:4], O = theta_MLE_PY_g[3:4], P = alpha_MLE_PY_g[3:4], Q = theta_MLE_DP_g[3:4]
)
colnames(df_geom_param) <- c(
  "geom param", "theta PY param", "alpha PY param", "theta DP param",
  "geom param", "theta PY param", "alpha PY param", "theta DP param"
)

xtable(df_geom)
xtable(df_geom_param)
print(xtable(df_zipf), file = "table.tex", sanitize.colnames.function = identity)


# # ------
# # compare minimax e MLE estimation of parameters
# # ------
#
# alpha_hat_geom <- c()
# theta_hat_geom <- c()
# for(i in 1:length(geom_param_list)){
#   freq = dataset_g[[i]]$frequencies
#   m_g <- rep(0, dataset_g[[i]]$n)
#   a = data.frame(freq) %>% group_by(freq) %>% summarise(count=n())
#   for(j in 1:dataset_g[[i]]$n){
#     if(j %in% a$freq) m_g[j] = a$count[match(j ,a$freq)]
#   }
#   ind_max = tail(which(m_g!=0),1)
#   cumsum_m = cumsum(rev(m_g[1:ind_max]))
#   alpha_hat_geom[i] = PY_alpha(freq, ind_max, cumsum_m, dataset_g[[i]]$K_n)
#   theta_hat_geom[i] = max_EPPF_PY_theta(dataset_g[[i]]$frequencies, alpha_hat_geom[i])$par
#   print(i)
# }
#
# alpha_MLE_geom <- c()
# for(i in 1:length(geom_param_list)) alpha_MLE_geom[i] = out_PY_g[[i]]$par[2]
# theta_MLE_geom <- c()
# for(i in 1:length(geom_param_list)) theta_MLE_geom[i] = out_PY_g[[i]]$par[1]
# logEPPF_PY_minimax_geom <- c()
# for(i in 1:length(geom_param_list)) logEPPF_PY_minimax_geom[i] = logEPPF_PY(theta_hat_geom[i], alpha_hat_geom[i], dataset_g[[i]]$frequencies)
# logEPPF_PY_MLE_geom <- c()
# for(i in 1:length(geom_param_list)) logEPPF_PY_MLE_geom[i] = logEPPF_PY(theta_MLE_geom[i], alpha_MLE_geom[i], dataset_g[[i]]$frequencies)
# df1 <- data.frame(a=alpha_hat_geom, b=alpha_MLE_geom, c=theta_hat_geom, d=theta_MLE_geom, e=logEPPF_PY_minimax_geom, f=logEPPF_PY_MLE_geom)
# knitr::kable(df1, col.names =c("alpha minimax", "alpha MLE", "theta minimax",
#                                "theta MLE", "log EPPF minimax", "log EPPF MLE"))


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
  theme(legend.position = "none") + xlab("") + ylab(expression(tau[1])) + ggtitle("Custom probabilities") + theme(plot.title = element_text(hjust = 0.5)) + geom_hline(yintercept = dataset$true_tau1, linetype = "dotted")
p
ggsave(p, file = "custom.eps", device = "eps")


# ------------------------
# Scenario 4: Graph for introduction: Zipf 3 parameters
# ------------------------

# plot together
N <- 1000000L # Important to us L, otherwise is not recognized as integer
n <- 100000L
p_mix <- list()
K_n_mix <- c()
tau1_bet_mix <- c()
tau1_skin_mix <- c()
tau1_PY_mix <- c()
tau1_DP_mix <- c()
PY_lower_mix <- c()
PY_upper_mix <- c()
DP_lower_mix <- c()
DP_upper_mix <- c()
true_tau1_mix <- c()
m1_mix <- c()
out_PY_mix <- list()
dataset_mix <- list()
param_list <- c(2, 1.75, 1.5)
type_data <- c("Zipf", "Zipf", "Zipf")
for (i in 1:length(param_list)) {
  param <- param_list[i]
  set.seed(123)
  if (type_data[i] == "Geom") {
    dataset_mix[[i]] <- dataset_creation_geom(n = n, N = N, p = param)
  } else {
    dataset_mix[[i]] <- dataset_creation_zipf(n = n, zipf_param = param, N = N)
  }

  # PY estimation
  out_PY_mix[[i]] <- max_EPPF_PY(dataset_mix[[i]]$frequencies)
  tau1_PY_mix[i] <- tau1_py(dataset_mix[[i]]$m1, dataset_mix[[i]]$n, out_PY_mix[[i]]$par[1], out_PY_mix[[i]]$par[2], dataset_mix[[i]]$N)
  PY_sim <- tau1_py_sim(dataset_mix[[i]]$frequencies, out_PY_mix[[i]]$par[1], out_PY_mix[[i]]$par[2], dataset_mix[[i]]$N)
  PY_lower_mix[i] <- quantile(PY_sim, 0.01 / 2)
  PY_upper_mix[i] <- quantile(PY_sim, 1 - 0.01 / 2)

  # Dirichlet process estimation
  out_DP_mix <- max_EPPF_DP(dataset_mix[[i]]$frequencies)
  tau1_DP_mix[i] <- tau1_dp(dataset_mix[[i]]$m1, dataset_mix[[i]]$n, out_DP_mix$par[1], dataset_mix[[i]]$N)
  DP_lower_mix[i] <- qhyper(0.01 / 2, out_DP_mix$par[1] + dataset_mix[[i]]$n - 1, dataset_mix[[i]]$N - dataset_mix[[i]]$n, dataset_mix[[i]]$m1)
  DP_upper_mix[i] <- qhyper(1 - 0.01 / 2, out_DP_mix$par[1] + dataset_mix[[i]]$n - 1, dataset_mix[[i]]$N - dataset_mix[[i]]$n, dataset_mix[[i]]$m1)

  # Bethlehem and Skinner estimators
  estim <- tau1_bs(dataset_mix[[i]]$frequencies, dataset_mix[[i]]$N)
  tau1_bet_mix[i] <- estim[1]
  tau1_skin_mix[i] <- estim[2]

  # Plot estimates and confidence intervals
  type <- c("PY", "DP", "B", "S")
  estimates <- c(tau1_PY_mix[i], tau1_DP_mix[i], tau1_bet_mix[i], tau1_skin_mix[i])
  df <- data.frame(
    type = factor(type), estim = estimates,
    lower_CI = c(PY_lower_mix[i], DP_lower_mix[i], NA, NA), upper_CI = c(PY_upper_mix[i], DP_upper_mix[i], NA, NA)
  )

  p_mix[[i]] <- ggplot(df, aes(type, estim))
  p_mix[[i]] <- p_mix[[i]] + geom_pointrange(aes(ymin = lower_CI, ymax = upper_CI), size = 1.2) + theme_bw() + theme(legend.position = "none") + xlab("") + ylab(expression(tau[1])) + ggtitle(paste0(type_data[i], " parameter: ", round(param, 7))) + theme(plot.title = element_text(hjust = 0.5)) + geom_hline(yintercept = dataset_mix[[i]]$true_tau1, linetype = "dotted")
}

plot_mix <- do.call(grid.arrange, c(p_mix, ncol = length(param_list)))

ggsave(plot_mix, height = 5, width = 7 * 2.5, file = "mix.png", device = "png")
ggsave(plot_mix, height = 5, width = 7 * 2.5, file = "mix.eps", device = "eps")
