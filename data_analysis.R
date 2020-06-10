# setwd("~/Documents/Privacy_git/Privacy")

library(tidyverse)
library(dplyr)
library(knitr)
library(poweRlaw)
library(HI)

rm(list=ls())
source("functions.R")

# --------------------------------
# IPMUS DATA
# --------------------------------

dataset   <- read.csv("data/usa_00002.csv", header = T)

# head(dataset)

# Subsetting the dataset so that AGE > 20. The dataset is overwritten
dataset <-dataset[dataset$AGE > 20,]
dataset <- dataset %>% dplyr::mutate(id = 1:n())

# Selection of the relevant variables and groups
tablestate <- dataset %>% group_by(REGION, HHINCOME, VALUEH, FAMSIZE, NCHILD, SEX, AGE, MARST,
                                     HCOVANY, EDUC, SCHLTYPE, EMPSTAT, OCC,
                                     INCTOT, FTOTINC, VETSTAT) %>% count()
# Frequencies
freq_state <- sort(tablestate$n)

# Graphical representation of 
m_pl = displ$new(freq_state)
est = estimate_xmin(m_pl)
m_pl$setXmin(est)
plot(m_pl)
id_state <- dataset %>% add_count(REGION, HHINCOME, VALUEH, FAMSIZE, NCHILD, SEX, AGE, MARST,
                                    HCOVANY, EDUC, SCHLTYPE, EMPSTAT, OCC,
                                    INCTOT, FTOTINC, VETSTAT) %>% filter(n==1) %>% select(id)
id_state <- id_state$id

# ---------------
# Now let's take some SAMPLES
# ---------------

source('functions.R')

estimates <- loop_estimate(data=dataset, id_uniq_pop=id_state, percentage=0.01, iter=2, 
                           nsamples=10000, burnin=2, case='posterior_samples', REGION, HHINCOME, 
                            VALUEH, FAMSIZE, NCHILD, SEX, AGE, MARST,
                            HCOVANY, EDUC, SCHLTYPE, EMPSTAT, OCC,
                            INCTOT, FTOTINC, VETSTAT)

# check convergence
plot(estimates$theta_sample, type='l')
plot(estimates$alpha_sample, type='l')

iter = 2
for(i in 1:iter){
  print(c(estimates[[i]]$tau1_true,estimates[[i]]$CI_PY, estimates[[i]]$theta, estimates[[i]]$alpha))
}




