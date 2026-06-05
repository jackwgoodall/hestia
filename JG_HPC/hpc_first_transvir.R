library(cmdstanr)
library(splines)
library(tidyverse)

source("../R/joint_model_functions.R")

load("STPN_hestia.Rdata")

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

obs_joint <- make_joint_obs_model(
  viral_pcr = list(
    
  
first_transvir_mod <- run_joint_model(obs_model = obs_joint,
                                  data = STPN_hestia,
                                  obs_cols = c("viral_pcr", "bacterial_pcr"),
                                  init_probs = c(90, 5, 1, 1, 1, 2),
                                  ih_cov = x_ih_all,
                                  eh_cov = x_eh_splines_all,
                                  file = file.path("..", "inst", "stan", "hmm_tv_cov_reduce_sum_joint.stan"),
                                  iter         = 1000,
                                  chains            = 4,
                                  parallel_chains   = 4,
                                  threads_per_chain = 13,    
                                  adapt_delta       = 0.9,
                                  max_treedepth     = 12
)

Sys.time()

first_transvir_mod$save_object(file = "data/first_transvir_mod.RDS")

first_transvir_mod$save_object(file = "data/first_transvir_mod.RDS", compress = "xz")

Sys.time()
