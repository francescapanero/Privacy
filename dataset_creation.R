# --------------------------------
# IPMUS DATA
# --------------------------------

library(dplyr)
library(readr)

rm(list=ls())
# Using readr for faster loading
data_full   <- read_csv("data/usa_00002.csv")

source('~/Google Drive/University/Lavori/Disclosure risk/Privacy/functions.R', echo=TRUE)

# head(dataset)

# Subsetting the dataset so that AGE > 20. The dataset is overwritten
data_full <- data_full %>% filter (AGE > 20) %>% mutate(id = 1:n())


# Selection of the relevant variables and groups
# colnames(data_full)
#  [1] "YEAR"     "REGION"   "STATEICP" "CITY"     "HHINCOME" "VALUEH"   "FAMSIZE" 
#[8] "NCHILD"   "SEX"      "AGE"      "MARST"    "RACE"     "RACED"    "HCOVANY" 
#[15] "EDUC"     "EDUCD"    "SCHLTYPE" "EMPSTAT"  "EMPSTATD" "OCC"      "INCTOT"  
#[22] "FTOTINC"  "VETSTAT"  "VETSTATD" "id"  

dataset <- data_full %>% group_by(NCHILD, RACED, SEX, EDUCD, INCTOT) %>% count()

# Frequencies of the full dataset
model_checking(sort(dataset$n),subset=50)

freq_full     <- dataset$n
points_full   <- factor(rep(1:length(freq_full), freq_full)) # it is important they are factor

N <- sum(freq_full)
K_N <- length(freq_full)

# ---------------------------------------------
# 1% dataset
# ---------------------------------------------

set.seed(123)
percentage <- 0.01
n          <- round(N*percentage)

# Full dataset vs   
points_obs  <- sample(points_full, n, replace=FALSE)

# Observed frequencies
freq_full     <- as.numeric(table(points_full))
freq_observed <- as.numeric(table(points_obs))

# True tau1
true_tau1 <- sum((freq_full == 1) & (freq_observed == 1))
freq_observed <- freq_observed[freq_observed >0]

data_1perc <- list(frequencies = freq_observed, n = sum(freq_observed),  N = sum(freq_full), K_n = sum(freq_observed>0), K_N = sum(freq_full>0), 
     m1 = sum(freq_observed==1), percentage = n/N, true_tau1 = true_tau1)

# ---------------------------------------------
# 5% dataset
# ---------------------------------------------

set.seed(123)
percentage <- 0.05
n          <- round(N*percentage)

# Full dataset vs   
points_obs  <- sample(points_full, n, replace=FALSE)

# Observed frequencies
freq_full     <- as.numeric(table(points_full))
freq_observed <- as.numeric(table(points_obs))

# True tau1
true_tau1     <- sum((freq_full == 1) & (freq_observed == 1))
freq_observed <- freq_observed[freq_observed >0]

data_5perc <- list(frequencies = freq_observed, n = sum(freq_observed),  N = sum(freq_full), K_n = sum(freq_observed>0), K_N = sum(freq_full>0), 
                   m1 = sum(freq_observed==1), percentage = n/N, true_tau1 = true_tau1)

# ---------------------------------------------
# 10% dataset
# ---------------------------------------------

set.seed(123)
percentage <- 0.1
n          <- round(N*percentage)

# Full dataset vs   
points_obs  <- sample(points_full, n, replace=FALSE)

# Observed frequencies
freq_full     <- as.numeric(table(points_full))
freq_observed <- as.numeric(table(points_obs))

# True tau1
true_tau1 <- sum((freq_full == 1) & (freq_observed == 1))
freq_observed <- freq_observed[freq_observed >0]

data_10perc <- list(frequencies = freq_observed, n = sum(freq_observed),  N = sum(freq_full), K_n = sum(freq_observed>0), K_N = sum(freq_full>0), 
                   m1 = sum(freq_observed==1), percentage = n/N, true_tau1 = true_tau1)

# ---------------------------------------------
# 20% dataset
# ---------------------------------------------

set.seed(123)
percentage <- 0.2
n          <- round(N*percentage)

# Full dataset vs   
points_obs  <- sample(points_full, n, replace=FALSE)

# Observed frequencies
freq_full     <- as.numeric(table(points_full))
freq_observed <- as.numeric(table(points_obs))

# True tau1
true_tau1 <- sum((freq_full == 1) & (freq_observed == 1))
freq_observed <- freq_observed[freq_observed >0]

data_20perc <- list(frequencies = freq_observed, n = sum(freq_observed),  N = sum(freq_full), K_n = sum(freq_observed>0), K_N = sum(freq_full>0), 
                   m1 = sum(freq_observed==1), percentage = n/N, true_tau1 = true_tau1)

# -------------------------------------------

rm(n,N,K_N,freq_full,freq_observed,points_full,points_obs,dataset,true_tau1,percentage,data_full)
save.image("data/IPMUS.RData")


library(knitr)
kable(data.frame(n = data_1perc$n, N = data_1perc$N, percentage = data_1perc$percentage, 
                 K_n =data_1perc$K_n,  K_N = data_1perc$K_N,  m1 =data_1perc$m1, 
                 true_tau1 =data_1perc$true_tau1))

kable(data.frame(n = data_5perc$n, N = data_5perc$N, percentage = data_5perc$percentage, 
                 K_n =data_5perc$K_n,  K_N = data_5perc$K_N,  m1 =data_5perc$m1, 
                 true_tau1 =data_5perc$true_tau1))

kable(data.frame(n = data_10perc$n, N = data_10perc$N, percentage = data_10perc$percentage, 
                 K_n =data_10perc$K_n,  K_N = data_10perc$K_N,  m1 =data_10perc$m1, 
                 true_tau1 =data_10perc$true_tau1))

kable(data.frame(n = data_20perc$n, N = data_20perc$N, percentage = data_20perc$percentage, 
                 K_n =data_20perc$K_n,  K_N = data_20perc$K_N,  m1 =data_20perc$m1, 
                 true_tau1 =data_20perc$true_tau1))
