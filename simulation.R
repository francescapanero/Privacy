rm(list=ls())
library(knitr)

source('functions.R')

# -------------------------------------------
# Scenario 1 - Low sigma
# -------------------------------------------
   
H     <- 3*10^6
N     <- 10^6
n     <- 10^5/5
alpha <- 0.6

dataset <- dataset_creation(n = n, alpha = alpha, H = H, N = N)

kable(data.frame(n = dataset$n, N = dataset$N, percentage = dataset$percentage, 
                 K_n =dataset$K_n,  K_N = dataset$K_N, H = dataset$H, m1 =dataset$m1, 
                 true_tau1 =dataset$true_tau1, 
                 alpha = dataset$alpha))


# CA application, Alternative dataset ----------------------

dataset  <- read.csv("california/conteggi_1.csv",header=FALSE)
dataset <- list( frequencies = dataset$V1[dataset$V1 > 0])
dataset$n <- sum(dataset$frequencies)
dataset$percentage <- 0.1
dataset$N <- 1150934
dataset$K_n <- length(dataset$frequencies)
dataset$K_N <- NA
dataset$H <- 3600000
dataset$true_tau1 <- read.csv("vero_tau/verotau_1.csv",header=FALSE)[1,1]
dataset$alpha <- NA
dataset$m1 <- sum(dataset$frequencies == 1)

kable(data.frame(n = dataset$n, N = dataset$N, percentage = dataset$percentage,
                 K_n =dataset$K_n,  K_N = dataset$K_N, H = dataset$H, m1 =dataset$m1,
                 true_tau1 =dataset$true_tau1,
                 alpha = dataset$alpha))

# -----------------------
# 1 ] DP estimation
# -----------------------

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

# plot(table(sim_dp))

# Binomial approximation
dataset$tau1_dp_binom <- dataset$m1 * dataset$n / dataset$N
dataset$lower_dp_binom <- qbinom(0.025,dataset$m1, dataset$n / dataset$N)
dataset$upper_dp_binom <- qbinom(0.975,dataset$m1, dataset$n / dataset$N)


# -----------------------
# 2 ] PY estimation
# -----------------------

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
# plot(table(sim_py))


# Binomial approximation
dataset$tau1_py_binom <- dataset$m1 * (dataset$n / dataset$N)^(1 - alpha_hat)
dataset$lower_py_binom <- qbinom(0.025,dataset$m1, (dataset$n / dataset$N)^(1 - alpha_hat))
dataset$upper_py_binom <- qbinom(0.975,dataset$m1, (dataset$n / dataset$N)^(1 - alpha_hat))


# Results -----------------------------------

#rm(list=setdiff(ls(), "dataset"))
kable(data.frame(n = dataset$n, N = dataset$N, percentage = dataset$percentage, K_n = dataset$K_n, K_N = dataset$K_N, H = dataset$H, m1 =dataset$m1, true_tau1 = dataset$true_tau1,alpha = dataset$alpha))

kable(data.frame(true_tau1 = dataset$true_tau1, 
                 tau_dp       = dataset$tau1_dp,  CI_dp = paste("[", dataset$lower_dp,", ", dataset$upper_dp,"]",sep=""),
                 tau_dp_binom = dataset$tau1_dp_binom,  CI_dp_binom = paste("[", dataset$lower_dp_binom,", ", dataset$upper_dp_binom,"]",sep=""),
                 tau_py       = dataset$tau1_py, CI_py = paste("[", dataset$lower_py,", ", dataset$upper_py,"]",sep=""),
                 tau_py_binom = dataset$tau1_py_binom,  CI_py_binom = paste("[", dataset$lower_py_binom,", ", dataset$upper_py_binom,"]",sep="")
))



save.image("scenario1.RData")
