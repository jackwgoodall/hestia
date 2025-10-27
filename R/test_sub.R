source("R/hestia_functions.R")
source("R/simulation.R")

library(rstan)
library(tidyverse)

# Basic SIR 
if(FALSE) {# Simulate data from a basic SIR model for 100 households
  # Everyone starts susceptible
  # Two observations, first informative about infection (e.g. culture or PCR, 
  # second informative about recovery (e.g. IgG))
  # dat_sim <- sim_sir(eh_prob = 0.01, ih_prob = 0.05, n_hh = 100,
  #                    hh_size = 1:5, tmax = 100, gamma = 1/5,
  #                    covs_eh = c(0, 0), covs_ih = c(0, 0),
  #                    obs_prob = list(c(0.05, 0.95, 0.05),
  #                                    c(0.01, 0.01, 0.8)),
  #                    start_prob = c(1, 0, 0),
  #                    complete_enroll = TRUE)
  # 
  # saveRDS(dat_sim, "sim_SIR_nocov_9.18.25.rds")
  
  dat_sim <- readRDS("sim_SIR_nocov_9.18.25.rds")
  dat <- dat_sim$obs %>%
    rename(pcr = y1,
           igg = y2)
  
  # Basic SIR
  # For now let's set the value of gamma
  inf_process <- make_infection_model(transmit("S", "I"),
                                      transit("I", "R", gamma = NA))
  
  obs_process <- make_observation_model(pcr = c("S" = 0.05, "I" = 0.95, "R" = 0.05),
                                        igg = c("S" = 0.01, "I" = 0.01, "R" = 0.8))
  
  epsilon <- 1e-10
  init_probs <- c(1-2*epsilon, epsilon, epsilon)
  dat_stan <- make_stan_data(inf_model = inf_process, obs_model = obs_process, data = dat,
                             init_probs = init_probs, epsilon = epsilon)
  
  stan_fit <- stan(file = "stan/hmm.stan",
                   data = dat_stan,
                   iter = 1000,
                   chains = 4,
                   cores = 4,
                   init = rep(list(list(beta_ih = logit(0.03),
                                        beta_eh = logit(0.03),
                                        logit_params = array(logit(0.5)))), 4))
  
  saveRDS(stan_fit, "sim_SIR_nocov_9.18.25_init_multtest_stanfit.rds")
  }

# SIIR 
if(TRUE) {
  # dat_sim <- sim_siir(eh_prob = 0.01, ih_prob = 0.05, n_hh = 100,
  #                     hh_size = 1:5, tmax = 100, gamma = c(1/5, 1/30), split = c(0.7, 0.3),
  #                     covs_eh = c(0, 0), covs_ih = c(0, 0),
  #                     obs_prob = list(c(0.05, 0.95, 0.95, 0.05),
  #                                     c(0.01, 0.01, 0.01, 0.8)),
  #                     start_prob = c(1, 0, 0, 0),
  #                     complete_enroll = TRUE)
  # saveRDS(dat_sim, "sim_SIIR_nocov_10.27.25.rds")
  
  # dat_sim <- sim_siir(eh_prob = 0.01, ih_prob = 0.05, n_hh = 500,
  #                     hh_size = 1:5, tmax = 100, gamma = c(1/5, 1/30), split = c(0.7, 0.3),
  #                     covs_eh = c(0, 0), covs_ih = c(0, 0),
  #                     obs_prob = list(c(0.05, 0.95, 0.95, 0.05),
  #                                     c(0.01, 0.01, 0.01, 0.8)),
  #                     start_prob = c(1, 0, 0, 0),
  #                     complete_enroll = TRUE)
  # saveRDS(dat_sim, "sim_SIIR_nocov_500hh_10.27.25.rds")
  
  # dat_sim <- readRDS("sim_SIIR_nocov_10.27.25.rds")
  dat_sim <- readRDS("sim_SIIR_nocov_500hh_10.27.25.rds")
  dat <- dat_sim$obs %>%
    rename(pcr = y1,
           igg = y2)

  inf_process <- make_infection_model(transmit("S", c("Ia", "Ic"), source = c("Ia", "Ic"), split = c(0.7, 0.3)),
                                      transit("Ia", "R", gamma_a = 1/5),
                                      transit("Ic", "R", gamma_c = 1/30))
  
  # inf_process <- make_infection_model(transmit("S", c("Ia", "Ic"), source = c("Ia", "Ic"), split = c(0.7, 0.3)),
  #                                     transit("Ia", "R", gamma_a = NA),
  #                                     transit("Ic", "R", gamma_c = 1/30))
  
  obs_process <- make_observation_model(pcr = c("S" = 0.05, "Ia" = 0.95, "Ic" = 0.95, "R" = 0.05),
                                        igg = c("S" = 0.01, "Ia" = 0.01, "Ic" = 0.01, "R" = 0.8))
  
  epsilon <- 1e-10
  init_probs <- c(1-3*epsilon, epsilon, epsilon, epsilon)
  dat_stan <- make_stan_data(inf_model = inf_process, obs_model = obs_process, data = dat,
                             init_probs = init_probs, epsilon = epsilon)
  
  stan_fit <- stan(file = "stan/hmm.stan",
                   data = dat_stan,
                   iter = 1000,
                   chains = 4,
                   cores = 4,
                   init = rep(list(list(beta_ih = logit(0.03),
                                        beta_eh = logit(0.03),
                                        logit_params = array(rep(logit(0.5), dat_stan$n_params)),
                                        logit_mult_params = array(rep(logit(0.5), dat_stan$n_mult_params)))), 4))
  
  saveRDS(stan_fit, "sim_SIIR_nocov_500hh_10.27.25_stanfit.rds")
}



if(FALSE) {
  
  res_stan <- readRDS("sim_SIIR_nocov_10.27.25_stanfit.rds")
  ch <- rstan::extract(res_stan)
  plot(ch$ih_prob)
  plot(ch$eh_prob)
  plot(ch$params)
  
}