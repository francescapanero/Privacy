# setwd("~/Documents/Privacy_git/Privacy")

library(tidyverse)
library(dplyr)
library(knitr)

rm(list = ls())

load("data/IPMUS.RData")
source("2_functions.R")

# -------------------------------
# 5% dataset
# -------------------------------

dataset <- data_5perc

# MODEL CHECKING ---------------

out_PY <- max_EPPF_PY(dataset$frequencies)
out_DP <- max_EPPF_DP(dataset$frequencies)

# Comparison between M_l and the asymptotic formula
# frequency_check_PY(dataset$frequencies)

# Comparison between M_l and the expected values
M_l <- as.numeric(table(factor(dataset$frequencies, levels = 1:dataset$n)))

# PY comparison
tab <- rbind(PY = expected_m_py(1:15, dataset$n, out_PY$par[2], out_PY$par[1]),
             DP = expected_m_dp(1:15, dataset$n, out_DP$par[1]),
             Data = M_l[1:15])

colnames(tab) <- 1:15
kable(tab, digits=0)
frequency_check_PY(dataset$frequencies)

# ---------------------------
# PY estimation
# ---------------------------

alpha <- 0.01 # Credible intervals percentage

out_PY <- max_EPPF_PY(dataset$frequencies)
tau1_PY <- tau1_py(dataset$m1, dataset$n, out_PY$par[1], out_PY$par[2], dataset$N)
PY_sim <- tau1_py_sim(dataset$frequencies, out_PY$par[1], out_PY$par[2], dataset$N)

PY_lower <- quantile(PY_sim, alpha / 2)
PY_upper <- quantile(PY_sim, 1 - alpha / 2)

# Dirichlet process estimation
out_DP <- max_EPPF_DP(dataset$frequencies)
tau1_DP <- tau1_dp(dataset$m1, dataset$n, out_DP$par[1], dataset$N)
DP_lower <- qhyper(alpha / 2, out_DP$par[1] + dataset$n - 1, dataset$N - dataset$n, dataset$m1)
DP_upper <- qhyper(1 - alpha / 2, out_DP$par[1] + dataset$n - 1, dataset$N - dataset$n, dataset$m1)

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

# -------------------------------
# 10% dataset
# -------------------------------

dataset <- data_10perc

# MODEL CHECKING ---------------

out_PY <- max_EPPF_PY(dataset$frequencies)
out_DP <- max_EPPF_DP(dataset$frequencies)

# Comparison between M_l and the asymptotic formula
# frequency_check_PY(dataset$frequencies)

# Comparison between M_l and the expected values
M_l <- as.numeric(table(factor(dataset$frequencies, levels = 1:dataset$n)))

# PY comparison
tab <- rbind(PY = expected_m_py(1:15, dataset$n, out_PY$par[2], out_PY$par[1]),
             DP = expected_m_dp(1:15, dataset$n, out_DP$par[1]),
             Data = M_l[1:15])

colnames(tab) <- 1:15
kable(tab, digits=0)

# ---------------------------
# PY estimation
# ---------------------------

alpha <- 0.01 # Credible intervals percentage

out_PY <- max_EPPF_PY(dataset$frequencies)
tau1_PY <- tau1_py(dataset$m1, dataset$n, out_PY$par[1], out_PY$par[2], dataset$N)
PY_sim <- tau1_py_sim(dataset$frequencies, out_PY$par[1], out_PY$par[2], dataset$N)

PY_lower <- quantile(PY_sim, alpha / 2)
PY_upper <- quantile(PY_sim, 1 - alpha / 2)

# Dirichlet process estimation
out_DP <- max_EPPF_DP(dataset$frequencies)
tau1_DP <- tau1_dp(dataset$m1, dataset$n, out_DP$par[1], dataset$N)
DP_lower <- qhyper(alpha / 2, out_DP$par[1] + dataset$n - 1, dataset$N - dataset$n, dataset$m1)
DP_upper <- qhyper(1 - alpha / 2, out_DP$par[1] + dataset$n - 1, dataset$N - dataset$n, dataset$m1)

# Summary
kable(data.frame(
  n = dataset$n, N = dataset$N, percentage = round(dataset$percentage, 2),
  m1 = dataset$m1,
  K_n = dataset$K_n,
  true_tau1 = dataset$true_tau1,
  tau1_py = tau1_PY, CI_PY = paste("[", PY_lower, ", ", PY_upper, "]", sep = ""),
  tau1_dp = tau1_DP, CI_DP = paste("[", DP_lower, ", ", DP_upper, "]", sep = "")
))
