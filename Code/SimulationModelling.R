# ---------------------------------------------------------------------------- #
# - FILE NAME:   SimulationModelling.R         
# - DATE:        10/10/2021
# - DESCRIPTION: Code to re-analyse the data with the simulated matrices
# - AUTHORS:     Pol Capdevila Lanzaco
# ---------------------------------------------------------------------------- #

rm(list = ls(all=T))

library(plyr)
library(tidyr)
library(brms)
library(bayestestR)
library(tidybayes)
library(dplyr)
library(rstan)
library(Rcompadre)

#Working directories

path = gsub("/Code", "", dirname(rstudioapi::getActiveDocumentContext()$path))

DataPath = paste0(path,"/Data")
ResultPath = paste0(path, "/Results") 

#Load  and prepare data----

setwd(DataPath)

# Load simulated matrices

load("SimData.RData")

#These options help Stan run faster:

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()-1)

#Set modelling parameters

iter = 8000
thin = 0.0005 * iter
warmup = 0.1 * iter

# Define priors 

priors <- c(prior(normal(0, 10), class = Intercept),
            prior(normal(0, 1), class = b))

# Random with progression, stasis and shrinkage

random_m <- random_m %>% drop_na(GenT, Fec, rlwr, rupr, xt)

model_random <- brm(mvbind(scale(xt), scale(rupr), scale(rlwr)) ~ scale(GenT) + scale(Fec) + scale(GenT):scale(Fec)+ scale(Dimension),
                    iter = iter, thin = thin, warmup = warmup,
                    prior= priors, 
                    control = list(adapt_delta = .975, max_treedepth = 20),
                    data = random_m,
                    family = gaussian(), 
                    cores = 18)
# Random with progression

random_mp <- random_mp %>% drop_na(GenT, Fec, rlwr, rupr, xt)

model_random_p <- brm(mvbind(scale(xt), scale(rupr), scale(rlwr)) ~ scale(GenT) + scale(Fec) + scale(GenT):scale(Fec)+ scale(Dimension),
                      iter = iter, thin = thin, warmup = warmup,
                      prior= priors, 
                      control = list(adapt_delta = .975, max_treedepth = 20),
                      data = random_m,
                      family = gaussian(), 
                      cores = 18)

# Save the models

setwd(ResultPath)
save(model_random_p,model_random,
     file = "SimulationModels.RData")
