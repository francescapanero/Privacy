# setwd("~/Documents/Privacy_git/Privacy")

library(tidyverse)
library(dplyr)
library(knitr)

rm(list=ls())

load("data/IPMUS.RData")
source("functions.R")

# -------------------------------
# 5% dataset
# -------------------------------

dataset <- data_5perc
alpha <- 0.01

#table(dataset$frequencies)
model_checking_PY(dataset$frequencies, percentage = 0.80)
model_checking_DP(dataset$frequencies, percentage = 0.80)

# Parameter estimation
out_PY  <- max_EPPF_PY(dataset$frequencies)
tau1_PY <- tau1_py(dataset$m1, dataset$n, out_PY$par[1], out_PY$par[2], dataset$N)
PY_sim  <- tau1_py_sim(dataset$frequencies, out_PY$par[1], out_PY$par[2],  dataset$N)
PY_sim2  <- tau1_py_sim2(dataset$frequencies, out_PY$par[1], out_PY$par[2],  dataset$N)

hist(PY_sim)
hist(PY_sim2)

PY_lower <- quantile(PY_sim, alpha/2)
PY_upper <- quantile(PY_sim, 1 - alpha/2)

# Dirichlet process estimation
out_DP    <- max_EPPF_DP(dataset$frequencies)
tau1_DP   <- tau1_dp(dataset$m1, dataset$n, out_DP$par[1], dataset$N)
DP_lower <- qhyper(alpha/2, out_DP$par[1] + dataset$n - 1, dataset$N - dataset$n, dataset$m1)
DP_upper <- qhyper(1 - alpha/2, out_DP$par[1] + dataset$n - 1, dataset$N - dataset$n, dataset$m1)

# Summary
kable(data.frame(n = dataset$n, N = dataset$N, percentage = round(dataset$percentage,2),
                 m1 = dataset$m1,
                 K_n = dataset$K_n,
                 true_tau1 = dataset$true_tau1,
                 tau1_py = tau1_PY, CI_PY = paste("[", PY_lower,", ", PY_upper,"]",sep=""),
                 tau1_dp = tau1_DP, CI_DP = paste("[", DP_lower,", ", DP_upper,"]",sep="")))

# -------------------------------
# 10% dataset
# -------------------------------

dataset <- data_10perc
alpha <- 0.01

#table(dataset$frequencies)
model_checking_PY(dataset$frequencies, percentage = 0.80)
model_checking_DP(dataset$frequencies, percentage = 0.80)


# Parameter estimation
out_PY  <- max_EPPF_PY(dataset$frequencies)
tau1_PY <- tau1_py(dataset$m1, dataset$n, out_PY$par[1], out_PY$par[2], dataset$N)

PY_sim  <- tau1_py_sim(dataset$frequencies, out_PY$par[1], out_PY$par[2],  dataset$N, R = 1000)
# PY_sim2  <- tau1_py_sim2(dataset$frequencies, out_PY$par[1], out_PY$par[2],  dataset$N, R = 1000)

#hist(PY_sim)
#hist(PY_sim2)

PY_lower <- quantile(PY_sim, alpha/2)
PY_upper <- quantile(PY_sim, 1 - alpha/2)

# Dirichlet process estimation
out_DP    <- max_EPPF_DP(dataset$frequencies)
tau1_DP   <- tau1_dp(dataset$m1, dataset$n, out_DP$par[1], dataset$N)
DP_lower <- qhyper(alpha/2, out_DP$par[1] + dataset$n - 1, dataset$N - dataset$n, dataset$m1)
DP_upper <- qhyper(1 - alpha/2, out_DP$par[1] + dataset$n - 1, dataset$N - dataset$n, dataset$m1)

# Summary
kable(data.frame(n = dataset$n, N = dataset$N, percentage = round(dataset$percentage,2),
                 m1 = dataset$m1,
                 K_n = dataset$K_n,
                 true_tau1 = dataset$true_tau1,
                 tau1_py = tau1_PY, CI_PY = paste("[", PY_lower,", ", PY_upper,"]",sep=""),
                 tau1_dp = tau1_DP, CI_DP = paste("[", DP_lower,", ", DP_upper,"]",sep="")))