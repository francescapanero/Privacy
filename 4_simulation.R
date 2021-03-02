rm(list = ls())

source("2_functions.R")

# Figure 1 in the paper: comparison of B, S, C and DP on Zipf data ----------------------

N <- 1000000L # The "L" is crucial, otherwise it is not recognized as an integer
n <- 100000L

# Initialization of the relevant quantities
zipf_params <- c(rep(1.5, 100), rep(1.75,100), rep(2,100))
Method <- c("NB", "PB-1", "PB-2", "NEB")
N_param <- length(zipf_params)
N_method <- length(Method)
  
# Results
results <- NULL


random$seed(0L)

set.seed(123)  # Select the appropriate seed

# Simulation
for (i in 1:N_param) {
  
  # Setting the seed of the Python code

  # Simulating the dataset from a Zipf law
  dataset <- dataset_creation_zipf(n = n, zipf_param = zipf_params[i], N = N)
  K_n <- dataset$K_n
  true_tau1 <- dataset$true_tau1
  m1 <- dataset$m1
  
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
    Zipf = zipf_params[i],
    tau1 = rep(true_tau1, N_method),
    Method = Method, 
    estimate = c(tau1_DP, tau1_bet, tau1_skin, tau1_cam),
    lower_CI = c(DP_lower, NA, NA, NA),
    upper_CI = c(DP_upper, NA, NA, NA)
  )
  results <- rbind(results, df)
}


results$Zipf <- paste("Power-law c =", results$Zipf)
results$Zipf <- factor(results$Zipf, levels = c("Power-law c = 2", "Power-law c = 1.75", "Power-law c = 1.5"))
results$error <- results$estimate - results$tau1

p <- ggplot(data = results, aes(x = Method, y = error)) + geom_boxplot() + theme_bw() + theme(legend.position = "none") + xlab("") + ylab(expression(tau[1] - hat(tau)[1])) + facet_wrap(.~Zipf, ncol = 3, scales = "free_y") #+ geom_hline(linetype="dotted", aes(yintercept=tau1))
p

# Save the plot, if needed
ggsave(p, file = "img/zipf_fig1.eps", height = 0.75*3, width = 0.75*10, device = "eps")

# Upper part of Table 1 and Table 2: -----------------------------------
# Compare estimation of PY, DP, B, S, C on different Zipf datasets

N <- 1000000L # The "L" is crucial, otherwise it is not recognized as an integer
n <- 100000L

# Initialization of the relevant quantities
zipf_params <- c(1.25, 1.5, 1.75, 2)
Method <- c("Sam.", "Beth.", "Sk.", "Cam.", "Pitman-Yor")
N_param <- length(zipf_params)
N_method <- length(Method)

# Results
results <- NULL

# Performing the simulation
for (i in 1:N_param) {
  
  # Setting the seed of the Python code
  random$seed(0L)
  
  # Simulating the dataset from a Zipf law
  dataset <- dataset_creation_zipf(n = n, zipf_param = zipf_params[i], N = N)
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
    Zipf = zipf_params[i],
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

results[,9:10]

# Upper part of Table 1
# xtable(results[,1:8])

# Lower part of Table 1: ----------------------------------------------
# Compare estimation of PY, DP, B, S, C on different Geometric datasets
N <- 1000000L # The "L" is crucial, otherwise it is not recognized as an integer
n <- 100000L

# Initialization of the relevant quantities
geom_params <- c(0.0001, 0.001)
Method <- c("Sam.", "Beth.", "Sk.", "Cam.", "Pitman-Yor")
N_param <- length(geom_params)
N_method <- length(Method)

# Results
results <- NULL

# Performing the simulation
for (i in 1:N_param) {
  
  # Setting the seed of the R code
  set.seed(123)
  
  # Simulating the dataset from a Zipf law
  dataset <- dataset_creation_geom(n = n, N = N, p = geom_params[i])  
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
    Geom = geom_params[i],
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

# Lower part of Table 1
# xtable(results[,1:8])

results[,9:10]

