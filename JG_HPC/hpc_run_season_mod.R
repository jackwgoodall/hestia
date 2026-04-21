library(cmdstanr)
library(splines)
library(tidyverse)

source("../R/hestia_functions.R")

load("data/bac_viral_sim.Rdata")

inf_process <- make_infection_model(transmit(from = "S", to = "I"),
                                    progress(from = "I", to = "S", gamma = NA))

obs_process <- make_observation_model(pcr = c("S" = 0.05, "I" = 0.95))

age_season_mod <- run_model(inf_model = inf_process, 
                            obs_model = obs_process, 
                            data = bac_viral_sim, 
                            file = file.path("inst", "stan", "hmm_tv_cov_reduce_sum.stan"),
                            init_probs = c(0.7, 
                                           0.3),
                            ih_cov = x_ih_all,   # 3D array [T, N, k]
                            eh_cov = x_eh_splines_all,   # 3D array [T, N, k]
                            save_chains = FALSE,
                            save_states = FALSE,
                            time_varying = TRUE,
                            iter         = 500,
                            chains            = 4,
                            parallel_chains   = 4,
                            threads_per_chain = 13,    
                            adapt_delta       = 0.9,
                            max_treedepth     = 12,
                            backend = "cmdstanr"
)

save(age_season_mod, file = "data/age_season_mod.Rdata")