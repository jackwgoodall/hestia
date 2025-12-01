
library(tidyverse)
library(rstan)

source("R/hestia_functions.R")
source("R/simulation.R")

gen_data <- FALSE
cluster_run <- TRUE
run_type <- "SIIR"

############### First SIR

# Everyone starts susceptible
# Two observations, first informative about infection (e.g. culture or PCR, 
# second informative about recovery (e.g. IgG))

if(gen_data) {
  if(run_type == "SIR") {
    for(i in 1:10) {
      dat_sim <- sim_sir(eh_prob = 0.01, ih_prob = 0.05, n_hh = 300,
                         hh_size = 1:5, tmax = 100, gamma = 1/5,
                         covs_eh = c(0, 0), covs_ih = c(0, 0),
                         obs_prob = list(c(0.05, 0.95, 0.05),
                                         c(0.01, 0.01, 0.8)),
                         start_prob = c(1, 0, 0),
                         complete_enroll = TRUE)
      
      saveRDS(dat_sim, paste0("tests/sims/data/sim_SIR_nocov_", i, ".rds"))
    }
  }
  
  if(run_type == "SEIR") {
    for(i in 1:10) {
      dat_sim <- sim_seir(eh_prob = 0.01, ih_prob = 0.05, n_hh = 500,
                          hh_size = 1:5, tmax = 100, gamma = 1/3, sigma = 1/2,
                          covs_eh = c(0, 0), covs_ih = c(0, 0),
                          obs_prob = list(c(0.05, 0.5, 0.9, 0.05),
                                          c(0.03, 0.03, 0.30, 0.03),
                                          c(0.02, 0.02, 0.02, 0.95)),
                          start_prob = c(1, 0, 0, 0),
                          complete_enroll = TRUE)
      
      saveRDS(dat_sim, paste0("tests/sims/data/sim_SEIR_nocov_", i, ".rds"))
    }
  }
  
  if(run_type == "SIIR") {
    for(i in 1:10) {
      print(i)
      dat_sim <- sim_siir(eh_prob = 0.01, ih_prob = 0.05, n_hh = 500,
                          hh_size = 1:5, tmax = 100, gamma = c(1/5, 1/3), split = c(0.7, 0.3),
                          covs_eh = c(0, 0), covs_ih = c(0, 0),
                          obs_prob = list(c(0.05, 0.95, 0.95, 0.05),
                                          c(0.01, 0.01, 0.01, 0.8),
                                          c(0.03, 1-1e-10, 0.03, 0.03)),
                          start_prob = c(1, 0, 0, 0),
                          complete_enroll = TRUE)
      
      saveRDS(dat_sim, paste0("tests/sims/data/sim_SIIR2_nocov_", i, ".rds"))
    }
  }
  
}

if(cluster_run) {
  .args  <- commandArgs(trailingOnly = T)
  snum <- as.numeric(.args)
  
  if(run_type == "SIR") {
    dat <- readRDS(paste0("tests/sims/data/sim_SIR_nocov_", snum, ".rds"))$obs %>%
      rename(pcr = y1,
             igg = y2)
    
    
    
    inf_process <- make_infection_model(transmit("S", "I"),
                                        transit("I", "R", gamma = NA))
    obs_process <- make_observation_model(pcr = c("S" = 0.05, "I" = 0.95, "R" = 0.05),
                                          igg = c("S" = 0.01, "I" = 0.01, "R" = 0.8))
    epsilon <- 1e-10
    init_probs <- c(1-2*epsilon, epsilon, epsilon)
    
    if(snum == 1) {
      res_stan <- run_model(inf_model = inf_process, obs_model = obs_process,
                            data = dat, init_probs = c(1-2*1e-10, 1e-10, 1e-10),
                            iter = 2000, cores = 4, save_chains = TRUE,
                            file = "stan/hmm.stan")  
    } else  {
      res_stan <- run_model(inf_model = inf_process, obs_model = obs_process,
                            data = dat, init_probs = c(1-2*1e-10, 1e-10, 1e-10),
                            iter = 2000, cores = 4, save_chains = TRUE, save_states = FALSE,
                            file = "stan/hmm.stan") 
    }
    
    saveRDS(res_stan$chains, paste0("tests/sims/model_results/res_SIR_nocov_", snum, ".rds"))
    
  }
  
  if(run_type == "SEIR") {
    dat <- readRDS(paste0("tests/sims/data/sim_SEIR_nocov_", snum, ".rds"))$obs %>%
      rename(pcr = y1,
             sym = y2,
             igg = y3)
    
    # Define infection process model
    inf_process <- make_infection_model(transmit(from = "S", to = "E", source = "I"),
                                        transit("E", "I", sigma = NA),
                                        transit("I", "R", gamma = NA))
    
    # Define observation process model
    obs_process <- make_observation_model(pcr = c("S" = 0.05, "E" = 0.50, "I" = 0.90, "R" = 0.05),
                                          sym = c("S" = 0.03, "E" = 0.03, "I" = 0.30, "R" = 0.03),
                                          igg = c("S" = 0.02, "E" = 0.02, "I" = 0.02, "R" = 0.95))
    epsilon <- 1e-10
    init_probs <- c(1-3*epsilon, epsilon, epsilon, epsilon)
    
    if(snum == 1) {
      res_stan <- run_model(inf_model = inf_process, obs_model = obs_process,
                            data = dat, init_probs = init_probs,
                            iter = 2000, cores = 4, save_chains = TRUE,
                            file = "stan/hmm.stan")  
    } else  {
      res_stan <- run_model(inf_model = inf_process, obs_model = obs_process,
                            data = dat, init_probs = init_probs,
                            iter = 2000, cores = 4, save_chains = TRUE, save_states = FALSE,
                            file = "stan/hmm.stan") 
    }
    
    saveRDS(res_stan$chains, paste0("tests/sims/model_results/res_SEIR_nocov_", snum, ".rds"))
    
  }
  
  if(run_type == "SIIR") {
    dat <- readRDS(paste0("tests/sims/data/sim_SIIR2_nocov_", snum, ".rds"))$obs %>%
      rename(pcr = y1,
             igg = y2,
             symp = y3)
    
    inf_process <- make_infection_model(transmit("S", c("Is", "Ia"), source = c("Is", "Ia"), split = "phi"),
                                        transit("Is", "R", gamma_s = NA),
                                        transit("Ia", "R", gamma_a = NA))
    
    obs_process <- make_observation_model(pcr = c("S" = 0.05, "Is" = 0.95, "Ia" = 0.95, "R" = 0.05),
                                          igg = c("S" = 0.01, "Is" = 0.01, "Ia" = 0.01, "R" = 0.8),
                                          symp = c("S" = 0.03, "Is" = 1-1e-10, "Ia" = 0.03,  "R" = 0.03))
    epsilon <- 1e-10
    init_probs <- c(1-3*epsilon, epsilon, epsilon, epsilon)
    
    
    res_stan <- run_model(inf_model = inf_process, obs_model = obs_process,
                          data = dat, init_probs = init_probs,
                          iter = 2000, cores = 4, save_chains = TRUE, save_states = FALSE,
                          file = "stan/hmm.stan") 
  
    
    saveRDS(res_stan$chains, paste0("tests/sims/model_results/res_SIIR2_nocov_", snum, ".rds"))
    
  }
  
}

if(FALSE) {
  
  #### SEIR model
  
  ch <- list()
  res_SEIR <- list()
  
  # Load in chains from simulation runs
  for(snum in 1:10) {
    print(snum)
    ch[[snum]] <- readRDS(paste0("tests/sims/model_results/res_SEIR_nocov_", snum, ".rds"))[c("ih_prob", "eh_prob", "params")]
  }
  
  # Get estiamtes from chains
  for(i in 1:10) {
    res_SEIR[[i]] <- data.frame(param = c("Intra-household", "Extra-household", "Latent period", "Recovery period"),
                               est = c(median(ch[[i]]$ih_prob),
                                       median(ch[[i]]$eh_prob),
                                       median(1/ch[[i]]$params[,1]),
                                       median(1/ch[[i]]$params[,2])),
                               ci_low = c(quantile(ch[[i]]$ih_prob, 0.025),
                                          quantile(ch[[i]]$eh_prob, 0.025),
                                          quantile(1/ch[[i]]$params[,1], 0.025),
                                          quantile(1/ch[[i]]$params[,2], 0.025)),
                               ci_high = c(quantile(ch[[i]]$ih_prob, 0.975),
                                           quantile(ch[[i]]$eh_prob, 0.975),
                                           quantile(1/ch[[i]]$params[,1], 0.975),
                                           quantile(1/ch[[i]]$params[,2], 0.975)),
                               sim = i)
  }
  res_SEIR <- bind_rows(res_SEIR)
  
  truth <- data.frame(param = c("Intra-household", "Extra-household", "Latent period", "Recovery period"),
                      yint = c(0.05, 0.01, 2, 3))
  
  custom_limits <- bind_rows(data.frame(param = c("Intra-household", "Extra-household", "Latent period", "Recovery period"),
                              y = c(0.05 - 0.02, 0.01 - 0.002, 2 - .3, 3 - .4),
                              x = 1),
                             data.frame(param = c("Intra-household", "Extra-household", "Latent period", "Recovery period"),
                                        y = c(0.05 + 0.02, 0.01 + 0.002, 2 + .3, 3 + .4),
                                        x = 1))
  
  # Plot it - absolute
  ggplot() +
    geom_point(data = res_SEIR, aes(x = sim, y = est)) +
    geom_errorbar(data = res_SEIR, aes(x = sim, ymin = ci_low, ymax = ci_high)) +
    geom_blank(data = custom_limits, aes(x = x, y = y)) +
    geom_hline(data = truth, aes(yintercept = yint), lty = "dashed") +
    facet_wrap(~param, scale = "free_y", nrow = 1) +
    theme_bw() +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    labs(y = "",
         x = "")
  
  # Plot it - relative
  ggplot(res_SEIR %>% left_join(truth)) +
    geom_point(aes(x = sim, y = est/yint)) +
    geom_errorbar(aes(x = sim, ymin = ci_low/yint, ymax = ci_high/yint)) +
    # geom_blank(data = custom_limits, aes(x = x, y = y)) +
    geom_hline(yintercept = 1, lty = "dashed") +
    facet_wrap(~param, nrow = 1) +
    theme_bw() +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    labs(y = "",
         x = "")
  
  
  #### SIIR model
  
  ch <- list()
  
  # Load in chains from simulation runs
  for(snum in 1:10) {
    print(snum)
    ch[[snum]] <- readRDS(paste0("tests/sims/model_results/res_SIIR2_nocov_", snum, ".rds"))[c("ih_prob", "eh_prob", "params", "mult_params")]
  }
  
  # Get estiamtes from chains
  res_SIIR <- list()
  for(i in 1:10) {
    res_SIIR[[i]] <- data.frame(param = c("Intra-household", "Extra-household", "Recovery period 1", "Recovery period 2", "Symptomatic proportion"),
                                est = c(median(ch[[i]]$ih_prob),
                                        median(ch[[i]]$eh_prob),
                                        median(1/ch[[i]]$params[,1]),
                                        median(1/ch[[i]]$params[,2]),
                                        median(ch[[i]]$mult_params)),
                                ci_low = c(quantile(ch[[i]]$ih_prob, 0.025),
                                           quantile(ch[[i]]$eh_prob, 0.025),
                                           quantile(1/ch[[i]]$params[,1], 0.025),
                                           quantile(1/ch[[i]]$params[,2], 0.025),
                                           quantile(ch[[i]]$mult_params, 0.025)),
                                ci_high = c(quantile(ch[[i]]$ih_prob, 0.975),
                                            quantile(ch[[i]]$eh_prob, 0.975),
                                            quantile(1/ch[[i]]$params[,1], 0.975),
                                            quantile(1/ch[[i]]$params[,2], 0.975),
                                            quantile(ch[[i]]$mult_params, 0.975)),
                                sim = i)
  }
  res_SIIR <- bind_rows(res_SIIR)
  
  truth <- data.frame(param = c("Intra-household", "Extra-household", "Recovery period 1", "Recovery period 2", "Symptomatic proportion"),
                      yint = c(0.05, 0.01, 5, 3, 0.7))
  
  custom_limits <- bind_rows(data.frame(param = c("Intra-household", "Extra-household", "Recovery period 1", "Recovery period 2", "Symptomatic proportion"),
                                        y = c(0.05 - 0.01, 0.01 - 0.0015, 5 - .8, 3 - .4, 0.7 + .06),
                                        x = 1),
                             data.frame(param = c("Intra-household", "Extra-household", "Recovery period 1", "Recovery period 2", "Symptomatic proportion"),
                                        y = c(0.05 + 0.01, 0.01 + 0.0015, 5 + .8, 3 + .4, 0.7 - 0.06),
                                        x = 1))
  
  # Plot it
  ggplot() +
    geom_point(data = res_SIIR, aes(x = sim, y = est)) +
    geom_errorbar(data = res_SIIR, aes(x = sim, ymin = ci_low, ymax = ci_high)) +
    geom_blank(data = custom_limits, aes(x = x, y = y)) +
    geom_hline(data = truth, aes(yintercept = yint), lty = "dashed") +
    facet_wrap(~param, scale = "free_y", nrow = 1) +
    theme_bw() +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    labs(y = "", x = "")
  
  # Plot it - relative
  ggplot(res_SIIR %>% left_join(truth)) +
    geom_point(aes(x = sim, y = est/yint)) +
    geom_errorbar(aes(x = sim, ymin = ci_low/yint, ymax = ci_high/yint)) +
    # geom_blank(data = custom_limits, aes(x = x, y = y)) +
    geom_hline(yintercept = 1, lty = "dashed") +
    facet_wrap(~param, nrow = 1) +
    theme_bw() +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    labs(y = "",
         x = "")
  
}

  
  
  






