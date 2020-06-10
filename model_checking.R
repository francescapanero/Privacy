rm(list=ls())

library(knitr)

source('functions.R')


dataset_creation <- function(n, alpha, H, N){
  
  # Expanded representation
  probs =  1/(1:H)^(alpha) #
  
  points_full <- factor(sample(1:H, N, prob = probs, replace=TRUE),levels=1:H)
  points_obs  <- sample(points_full, n, replace=FALSE)
  
  # Observed frequencies
  freq_full     <- as.numeric(table(points_full))
  freq_observed <- as.numeric(table(points_obs))
  
  # True tau1
  true_tau1 <- sum((freq_full == 1) & (freq_observed == 1))
  freq_observed <- freq_observed[freq_observed >0]
  
  list(frequencies = freq_observed, n = sum(freq_observed), H = H, N = sum(freq_full), 
       K_n = sum(freq_observed>0), K_N = sum(freq_full>0), 
       m1 = sum(freq_observed==1), percentage = n/N, alpha = alpha, true_tau1 = true_tau1)
  
}

H     <- 3*10^6
N     <- 10^6 + 2382634
n     <- 10^5/5 + 2322
alpha <- 1.6

dataset <- dataset_creation(n = n, alpha = alpha, H = H, N = N)

kable(data.frame(n = dataset$n, N = dataset$N, percentage = dataset$percentage, 
                 K_n =dataset$K_n,  K_N = dataset$K_N, H = dataset$H, m1 =dataset$m1, 
                 true_tau1 =dataset$true_tau1, 
                 alpha = dataset$alpha))

model_checking(dataset$frequencies)
