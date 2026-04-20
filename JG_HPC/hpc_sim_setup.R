library(splines)
library(mgcv)
library(dplyr)
source("../R/simulation.R")

# HPC setup 
## Setup viral
age_spec_viral <- list(
  name = "age",
  levels = 1:4,
  probs = c(0.154, 0.407, 0.365, 0.0742),   # sampling distribution
  baseline = 1,                           #  reference level
  coefs_eh = c(-0.17741, -0.38282, -0.41593),   # levels 2..5
  coefs_ih = c(-0.17741, -0.38282, -0.41593))


viral_base <- sim_sirs(eh_prob = 0.02, 
                       ih_prob = 0.03, 
                       n_hh = 52,
                       hh_size = 5:12, 
                       tmax = 480, 
                       rho = 1/10,
                       gamma = 1/6,
                       covs_eh = NULL, 
                       covs_ih = NULL,
                       cat_covs = list(age_spec_viral),
                       obs_prob = list(c(0.005, 0.90, 0.005), # S&+ , I&+
                                       c(0.01, 0.05, 0.8)), # S&+ , I&+
                       start_prob = c(0.74, 0.06, 0.2),
                       complete_enroll = TRUE,
                       sigma_pid = 0,                   # standard deviation by individual on logit scale
                       season_type = "spline",
                       season_k = NULL,                    # number of Fourier harmonics
                       season_period = 365,             # period in days
                       season_df = 6,                   # spline basis dimension (mgcv::s k)
                       season_knots = list(day = c(0.5, 365.5)),             # optional knots list for mgcv::s (e.g., list(day = c(0.5, 365.5)))
                       season_coefs_eh = c(0.101632014,  
                                           -0.003037513, 
                                           0.390994253, 
                                           -0.230518142),          # length 2*season_k (sin1, cos1, sin2, cos2, ...)
                       season_coefs_ih = NULL)


age_spec <- list(
  name = "age",
  levels = 1:4,
  probs = c(0.154, 0.407, 0.365, 0.0742),   # sampling distribution
  baseline = 1,                           #  reference level
  coefs_eh = c(-0.3988, -0.7536, -0.9236),   # levels 2..5
  coefs_ih = c(-0.3988, -0.7536, -0.9236))


bacterial_base <- sim_sis_from_existing(
  base_complete_obs = viral_base$complete_obs,
  base_x = viral_base$x,
  a_complete_obs = viral_base$complete_obs,
  eh_prob = 0.04, 
  ih_prob = 0.01, 
  gamma = 1/12,
  covs_ih = c(age_2 = -0.3988,
              age_3 = -0.7536,
              age_4 = -0.9236),
  covs_eh = c(age_2 = -0.3988,
              age_3 = -0.7536,
              age_4 = -0.9236),
  covs_gamma = c(age_2 = 0.3988,
                 age_3 = 0.7536,
                 age_4 = 0.9236),
  cross_eh_coef = 0.2,
  cross_ih_susc_coef = 0.2,
  cross_ih_trans_coef = 0.3,
  season_type = "spline",
  season_k = NULL,                
  season_period = 365,             # period in days
  season_df = 6,                   # spline basis dimension (mgcv::s k)
  season_knots = list(day = c(0.5, 365.5)),             # optional knots list for mgcv::s (e.g., list(day = c(0.5, 365.5)))
  season_coefs_eh = c(0.1578415,  0.1509227, -1.1588201, -0.7021107),          # length 2*season_k (sin1, cos1, sin2, cos2, ...)
  season_coefs_ih = NULL
)

## Make array

bac_viral_sim <- bacterial_base$obs %>%
  mutate(state = if_else(t %in% c(seq(7,700,7)), state, NA)) %>%
  rename(pcr = y1)

bacterial_base_covs <- bacterial_base$x[,2:4]

spline_basis <- splines::ns(1:max(bac_viral_sim$t), df = 5)
spline_basis <- as.matrix(spline_basis)

x_eh_splines <- array(dim = c(max(bac_viral_sim$t), nrow(bacterial_base_covs), ncol(spline_basis)))


for(n in 1:max(bac_viral_sim$t)) {
  x_eh_splines[n, , ] <- matrix(rep(spline_basis[n, ], each = nrow(bacterial_base_covs)),
                                nrow = nrow(bacterial_base_covs), 
                                ncol = ncol(spline_basis))
}

x_eh_splines_all <- array(dim = c(max(bac_viral_sim$t), nrow(bacterial_base_covs), ncol(spline_basis) + ncol(bacterial_base_covs)))

for(n in 1:max(bac_viral_sim$t)) {
  x_eh_splines_all[n, , ] <- cbind(x_eh_splines[n, , ], as.matrix(bacterial_base_covs))
}


x_ih_all <- array(dim = c(max(bac_viral_sim$t), nrow(bacterial_base_covs), ncol(bacterial_base_covs)))

for(n in 1:max(bac_viral_sim$t)) {
  x_ih_all[n, , ] <- as.matrix(bacterial_base_covs)
}


save(bac_viral_sim, bacterial_base_covs, x_eh_splines_all, x_ih_all, file = "data/bac_viral_sim.Rdata")
