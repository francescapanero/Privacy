# setwd("~/Documents/Privacy_git/Privacy")

library(tidyverse)
library(dplyr)
library(knitr)

rm(list=ls())

load("data/IPMUS.RData")
source("functions.R")

# -------------------------------
# 1% dataset
# -------------------------------

dataset <- data_10perc
alpha <- 0.01

model_checking(data_10perc$frequencies)

# Parameter estimation
out_PY  <- max_EPPF_PY(dataset$frequencies)
tau1_PY <- tau1_py(dataset$m1, dataset$n, out_PY$par[1], out_PY$par[2], dataset$N)

# Dirichlet process estimation
out_DP    <- max_EPPF_DP(dataset$frequencies)
tau1_DP   <- tau1_dp(dataset$m1, dataset$n, out_DP$par[1], dataset$N)
DP1_sim   <- tau1_py_sim(dataset$frequencies, out_DP$par[1], 0, dataset$N)
DP1_lower <- quantile(DP1_sim, alpha/2)
DP1_upper <- quantile(DP1_sim, 1 - alpha/2)

# # Binomial approximation
# tau1_py_binom1  <- m1_state * (n0 / N)^(1 - out_PY$par[2])
# lower_py_binom1 <- qbinom(alpha/2, m1_state, (n0 / N)^(1 - out_PY$par[2]))
# upper_py_binom1 <- qbinom(1 - alpha/2, m1_state, (n0 / N)^(1 - out_PY$par[2]))

# Summary
kable(data.frame(n = dataset$n, N = dataset$N, percentage = round(dataset$percentage,2),
                 m1 = dataset$m1,
                 K_n = dataset$K_n,
                 true_tau1 = dataset$true_tau1,
                 tau1_py = tau1_PY, CI_PY = paste("[", PY1_lower,", ", PY1_upper,"]",sep=""),
                 #tau1_py_binapprox = tau1_py_binom1,
                 #CI_PY_binapprox = paste("[", lower_py_binom1,", ", upper_py_binom1,"]",sep=""),
                 tau1_dp = tau1_DP, CI_DP = paste("[", DP1_lower,", ", DP1_upper,"]",sep="")))

