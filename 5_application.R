library(tidyverse)
library(ggpubr)

rm(list = ls())

load("data/IPMUS.RData")
source("2_functions.R")

data_list <- list(data_5perc, data_10perc)

# MODEL CHECKING

## 5% case -------------------------------

dataset <- data_list[[1]]
out_PY <- max_EPPF_PY(dataset$frequencies)
out_DP <- max_EPPF_DP(dataset$frequencies)

# Comparison between M_l and the expected values
M_l <- as.numeric(table(factor(dataset$frequencies, levels = 1:dataset$n)))

# Table creation
tab <- rbind(
  PY = expected_m_py(1:10, dataset$n, out_PY$par[2], out_PY$par[1]),
  DP = expected_m_dp(1:10, dataset$n, out_DP$par[1]),
  Data = M_l[1:10]
)

colnames(tab) <- 1:10
kable(tab, digits = 0)
xtable(tab, digits=0)

p_check <- frequency_check_PY(dataset$frequencies)
p_check <- p_check + ggtitle(paste("Dataset percentage: 5%")) + theme(plot.title = element_text(hjust = 0.5, size = 10))
p_check

## 10% case -----------------------------------

dataset <- data_list[[2]]
out_PY <- max_EPPF_PY(dataset$frequencies)
out_DP <- max_EPPF_DP(dataset$frequencies)

# Comparison between M_l and the expected values
M_l <- as.numeric(table(factor(dataset$frequencies, levels = 1:dataset$n)))

# Table creation
tab <- rbind(
  PY = expected_m_py(1:10, dataset$n, out_PY$par[2], out_PY$par[1]),
  DP = expected_m_dp(1:10, dataset$n, out_DP$par[1]),
  Data = M_l[1:10]
)

colnames(tab) <- 1:10
kable(tab, digits = 0)
xtable(tab, digits=0)

p1_check <- frequency_check_PY(dataset$frequencies)
p1_check <- p1_check + ggtitle(paste("Dataset percentage: 10%")) + theme(plot.title = element_text(hjust = 0.5, size = 10))
p1_check


# Figure 3 in the paper
p_check_tog <- list()
p_check_tog[[1]] <- p_check
p_check_tog[[2]] <- p1_check
a <- do.call(grid.arrange, c(p_check_tog, ncol = 2))

ggsave(a, height = 3, width = 8, file = "img/check.eps", device="eps")

# Model estimation ----------------------------------------

# Initialization of the relevant quantities
Method <- c("Sam.", "Beth.", "Sk.", "Cam.", "Pitman-Yor")
N_param <- 2
N_method <- length(Method)

# Results
results <- NULL

# Performing the simulation
for (i in 1:N_param) {

  # Simulating the dataset from a Zipf law
  dataset <- data_list[[i]]
  n <- dataset$n
  N <- dataset$N
  K_n <- dataset$K_n
  true_tau1 <- dataset$true_tau1
  m1 <- dataset$m1
  
  set.seed(123)  # Select the appropriate seed
  # PY estimation
  out_PY <- max_EPPF_PY(dataset$frequencies)
  tau1_PY <- tau1_py(dataset$m1, dataset$n, out_PY$par[1], out_PY$par[2], dataset$N)
  PY_sim <- tau1_py_sim(dataset$frequencies, out_PY$par[1], out_PY$par[2], dataset$N)
  PY_lower <- quantile(PY_sim, 0.01 / 2)
  PY_upper <- quantile(PY_sim, 1 - 0.01 / 2)
  theta_MLE_PY <- out_PY$par[1]
  alpha_MLE_PY <- out_PY$par[2]
  
  # Dirichlet process estimation
  out_DP <- max_EPPF_DP(dataset$frequencies)
  tau1_DP <- tau1_dp(dataset$m1, dataset$n, out_DP$par[1], dataset$N)
  DP_lower <- qhyper2(0.01 / 2, out_DP$par[1] + dataset$n - 1, dataset$N - dataset$n, dataset$m1)
  DP_upper <- qhyper2(1 - 0.01 / 2, out_DP$par[1] + dataset$n - 1, dataset$N - dataset$n, dataset$m1)
  theta_MLE_DP <- out_DP$par
  
  # Bethlehem and Skinner and Camerlenghi estimators
  estim <- tau1_bs(dataset$frequencies, dataset$N)
  tau1_bet <- estim[1]
  tau1_skin <- estim[2]
  tau1_cam <- tau1_np_pois(dataset$N, dataset$n, dataset$frequencies)
  
  df <- data.frame(
    n = n,
    m1 = m1,
    tau1 = true_tau1,
    PY = paste(round(tau1_PY), " [", round(PY_lower), ", ", round(PY_upper), "]", sep = ""),
    Sam = paste(round(tau1_DP), " [", round(DP_lower), ", ", round(DP_upper), "]", sep = ""),
    Beth = round(tau1_bet), 
    Sk = round(tau1_skin), 
    Cam = round(tau1_cam),
    alpha_hat = round(alpha_MLE_PY,2),
    theta_hat = round(theta_MLE_PY,2)
  )
  results <- rbind(results, df)
}

knitr::kable(results)
xtable(results[,1:8])
