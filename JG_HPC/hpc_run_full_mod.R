library(cmdstanr)
library(splines)
library(tidyverse)

source("../R/joint_model_functions.R")

load("data/bac_viral_sim.Rdata")

inf_process <- make_infection_model(transmit(from = "S", to = "I"),
                                    progress(from = "I", to = "S", gamma = NA))

obs_joint <- make_joint_obs_model(
   viral_pcr = list(
       Sv_Sb = c(1, 99),   # FPR ~ 0.01
       Sv_Ib = c(1, 99),
       Iv_Sb = c(95, 5),   # TPR ~ 0.95
       Iv_Ib = c(95, 5),
       Rv_Sb = c(1, 99),
       Rv_Ib = c(1, 99)
     ),
     bac_pcr = list(
       Sv_Sb = c(1, 99),   # FPR ~ 0.01
       Sv_Ib = c(95, 5),   # TPR ~ 0.95
       Iv_Sb = c(1, 99),
       Iv_Ib = c(95, 5),   
       Rv_Sb = c(1, 99),
       Rv_Ib = c(95, 5)    
     )
   )

age_season_mod2 <- run_model(inf_model = inf_process, 
                             obs_model = obs_process_beta, 
                             data = bac_viral_sim, 
                             file = file.path("..", "inst", "stan", "hmm_tv_cov_reduce_sum_obs_prior.stan"),
                             init_probs = c(0.7, 
                                            0.3),
                             ih_cov = x_ih_all,   # 3D array [T, N, k]
                             eh_cov = x_eh_splines_all,   # 3D array [T, N, k]
                             save_chains = FALSE,
                             save_states = FALSE,
                             time_varying = TRUE,
                             iter         = 1000,
                             chains            = 4,
                             parallel_chains   = 4,
                             threads_per_chain = 13,    
                             adapt_delta       = 0.9,
                             max_treedepth     = 12,
                             backend = "cmdstanr"
)

age_season_mod2$save_object(file = "data/age_season_mod2.RDS")

age_season_mod2$save_object(file = "data/age_season_mod2_comp.RDS", compress = "xz")
