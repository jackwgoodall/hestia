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

# SEIR
if(FALSE) {
  
  dat_sim <- sim_seir(eh_prob = 0.01, ih_prob = 0.05, n_hh = 500,
                      hh_size = 1:5, tmax = 100, gamma = 1/3, sigma = 1/2,
                      covs_eh = c(0, 0), covs_ih = c(0, 0),
                      obs_prob = list(c(0.05, 0.5, 0.9, 0.05),
                                      c(0.03, 0.03, 0.30, 0.03),
                                      c(0.02, 0.02, 0.02, 0.95)),
                      start_prob = c(1, 0, 0, 0),
                      complete_enroll = TRUE)
  saveRDS(dat_sim, "sim_SEIR_nocov_500hh_11.17.25.rds")
  
  dat_sim <- readRDS("sim_SEIR_nocov_500hh_11.17.25.rds")
  dat <- dat_sim$obs %>%
    rename(pcr = y1,
           sym = y2,
           igg = y3)
  
  # inf_process <- make_infection_model(transmit(from = "S", to = "E", source = "I"),
  #                                     transit("E", "I", sigma = 1/2),
  #                                     transit("I", "R", gamma = 1/3))
  
  # inf_process <- make_infection_model(transmit(from = "S", to = "E", source = "I"),
  #                                     transit("E", "I", sigma = 1/2),
  #                                     transit("I", "R", gamma = NA))
  
  # Define infection process model
  inf_process <- make_infection_model(transmit(from = "S", to = "E", source = "I"),
                                      progress("E", "I", sigma = NA),
                                      progress("I", "R", gamma = NA))
  
  # Define observation process model
  obs_process <- make_observation_model(pcr = c("S" = 0.05, "E" = 0.50, "I" = 0.90, "R" = 0.05),
                                        sym = c("S" = 0.03, "E" = 0.03, "I" = 0.30, "R" = 0.03),
                                        igg = c("S" = 0.02, "E" = 0.02, "I" = 0.02, "R" = 0.95))
  
  # Run model
  res_stan <- run_model(inf_model = inf_process, obs_model = obs_process,
                        data = dat, init_probs = c(1-3*1e-10, 1e-10, 1e-10, 1e-10),
                        iter = 2000, cores = 4, save_chains = TRUE,
                        file = "stan/hmm.stan")

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
  
  saveRDS(stan_fit, "sim_SEIR_nocov_500hh_11.17.25_fitgs_stanfit.rds")
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
  
  # dat_sim <- sim_siir(eh_prob = 0.01, ih_prob = 0.05, n_hh = 500,
  #                     hh_size = 1:5, tmax = 100, gamma = c(1/5, 1/30), split = c(0.7, 0.3),
  #                     covs_eh = c(0, 0), covs_ih = c(0, 0),
  #                     obs_prob = list(c(0.05, 0.95, 0.95, 0.05),
  #                                     c(0.01, 0.01, 0.01, 0.8),
  #                                     c(0.03, 0.9, 0.03, 0.03)),
  #                     start_prob = c(1, 0, 0, 0),
  #                     complete_enroll = TRUE)
  # saveRDS(dat_sim, "sim_SIIR_nocov_500hh_infind_10.27.25.rds")
  
  # dat_sim <- readRDS("sim_SIIR_nocov_10.27.25.rds")
  # dat_sim <- readRDS("sim_SIIR_nocov_500hh_10.27.25.rds")
  dat_sim <- readRDS("sim_SIIR_nocov_500hh_infind_10.27.25.rds")
  dat <- dat_sim$obs %>%
    rename(pcr = y1,
           igg = y2,
           inf_ind = y3)

  # inf_process <- make_infection_model(transmit("S", c("Ia", "Ic"), source = c("Ia", "Ic"), split = c(0.7, 0.3)),
  #                                     transit("Ia", "R", gamma_a = 1/5),
  #                                     transit("Ic", "R", gamma_c = 1/30))
  
  # inf_process <- make_infection_model(transmit("S", c("Ia", "Ic"), source = c("Ia", "Ic"), split = "phi"),
  #                                     transit("Ia", "R", gamma_a = 1/5),
  #                                     transit("Ic", "R", gamma_c = 1/30))
  
  # inf_process <- make_infection_model(transmit("S", c("Ia", "Ic"), source = c("Ia", "Ic"), split = c(0.7, 0.3)),
  #                                     transit("Ia", "R", gamma_a = NA),
  #                                     transit("Ic", "R", gamma_c = 1/30))
  
  # Define infection process model
  inf_process <- make_infection_model(transmit("S", c("Ia", "Ic"), source = c("Ia", "Ic"), split = "phi"),
                                      progress("Ia", "R", gamma_a = NA),
                                      progress("Ic", "R", gamma_c = NA))
  
  # Define observation process model
  obs_process <- make_observation_model(pcr = c("S" = 0.05, "Ia" = 0.95, "Ic" = 0.95, "R" = 0.05),
                                        igg = c("S" = 0.01, "Ia" = 0.01, "Ic" = 0.01, "R" = 0.8),
                                        inf_ind = c("S" = 0.03, "Ia" = 0.9, "Ic" = 0.03, "R" = 0.03))
  
  # Run model
  # res_stan <- run_model(inf_model = inf_process, obs_model = obs_process,
  #                       data = dat, init_probs = c(1-3*1e-10, 1e-10, 1e-10, 1e-10),
  #                       iter = 2000, cores = 4, save_chains = TRUE,
  #                       file = "stan/hmm.stan")
  
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
  
  saveRDS(stan_fit, "sim_SIIR_nocov_fitsplit_fit2gamma_500hh_infind_10.27.25_stanfit.rds")
}

if(FALSE) {
  
  res_stan <- readRDS("sim_SIIR_nocov_fitsplit_500hh_infind_10.27.25_stanfit.rds")
  res_stan <- readRDS("sim_SIIR_nocov_fitsplit_fit2gamma_500hh_infind_10.27.25_stanfit.rds")
  res_stan <- readRDS("sim_SEIR_nocov_500hh_11.17.25_stanfit.rds")
  res_stan <- readRDS("sim_SEIR_nocov_500hh_11.17.25_fitg_stanfit.rds")
  res_stan <- readRDS("sim_SEIR_nocov_500hh_11.17.25_fitgs_stanfit.rds")
  ch <- rstan::extract(res_stan)
  plot(ch$ih_prob)
  plot(ch$eh_prob)
  plot(ch$params)
  
  dat_state <- dat %>%
    group_by(t, state) %>%
    summarize(tot = n()) %>%
    mutate(state = ifelse(state == 1, "S", ifelse(state == 2, "Is", ifelse(state == 3, "Ia", "R"))),
           state = factor(state, levels = c("S", "Is", "Ia", "R")))
  
  ggplot(dat_state, aes(x = t, y = tot, color = state)) +
    geom_line() +
    theme_bw() +
    labs(x = "Day", y = "Count", color = "State") +
    theme(legend.position = "bottom") +
    scale_color_brewer(palette =  "Set1")
  
  dat_obs <- dat %>%
    pivot_longer(c(pcr, igg, symp)) %>%
    group_by(t, name) %>%
    summarize(tot = sum(value))
  
  
  ggplot(dat_obs, aes(x = t, y = tot, color = name)) +
    geom_line() +
    theme_bw() +
    labs(x = "Day", y = "Count", color = "Observation") +
    theme(legend.position = "bottom") 
  
  res_SEIR <- data.frame(param = c("Intra-household", "Extra-household", "Latent period", "Recovery period"),
                             est = c(median(ch$ih_prob),
                                     median(ch$eh_prob),
                                     median(1/ch$params[,1]),
                                     median(1/ch$params[,2])),
                             ci_low = c(quantile(ch$ih_prob, 0.025),
                                        quantile(ch$eh_prob, 0.025),
                                        quantile(1/ch$params[,1], 0.025),
                                        quantile(1/ch$params[,2], 0.025)),
                             ci_high = c(quantile(ch$ih_prob, 0.975),
                                         quantile(ch$eh_prob, 0.975),
                                         quantile(1/ch$params[,1], 0.975),
                                         quantile(1/ch$params[,2], 0.975)))
  
  truth <- data.frame(param = c("Intra-household", "Extra-household", "Latent period", "Recovery period"),
                      yint = c(0.05, 0.01, 2, 3))
  
  ggplot(res_SEIR %>% mutate(sim = 1), aes(x = sim, y = est, ymin = ci_low, ymax = ci_high)) +
    geom_point() +
    geom_errorbar() +
    geom_hline(data = truth, aes(yintercept = yint), lty = "dashed") +
    facet_wrap(~param, scale = "free_y", nrow = 1) +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
    
  
}