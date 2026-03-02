
make_outcome <- function(df, col, p1) {
  mutate(df, )
}

sim_sir <- function(eh_prob = 0.01, ih_prob = 0.05, n_hh = 100,
                    hh_size = 1:5, tmax = 100, gamma = 1/5,
                    covs_eh = c(0, 0), covs_ih = c(0, 0),
                    obs_prob = list(c(0.05, 0.95, 0.05),
                                    c(0.01, 0.01, 0.8)),
                    start_prob = c(1, 0, 0),
                    complete_enroll = TRUE) {

  epsilon <- 1e-10

  hh_size <- sample(hh_size, n_hh, replace = TRUE) # household sizes
  enroll_per_hh <- numeric(n_hh) # number of participants enrolled per HH

  x <- matrix(nrow = sum(hh_size), ncol = length(covs_ih))

  for(i in 1:length(covs_ih)) {
    x[,i] <- rbinom(sum(hh_size), 1, 0.4)
  }

  # Create participant IDs
  part_ids <- list()
  for(i in 1:n_hh) {
    part_ids[[i]] <- 1:(hh_size[i])
    if(hh_size[i] == 1) {
      enroll_per_hh[i] <- 1
    } else {
      enroll_per_hh[i] <- sample(1:hh_size[i], 1)
    }
  }

  epsilon <- 1e-10
  enroll_ids <- list()

  # Infection state for all household members on all time steps
  complete_obs <- data.frame(t = numeric(),
                             part_id = numeric(),
                             enroll = numeric(),
                             state = numeric(),
                             hh_size = numeric(),
                             hh_id = numeric())

  last_x <- 0

  for(i in 1:length(hh_size)) {

    # Move HH members through SIR states
    for(d in 1:tmax) {
      wk <- base::ceiling(d/7)
      if(d == 1) {
        new_obs <- bind_rows(data.frame(t = rep(d, hh_size[i]),
                                        part_id = part_ids[[i]],
                                        enroll = c(rep(1, enroll_per_hh[i]), rep(0, hh_size[i]-enroll_per_hh[i])),
                                        state = sample(1:3, hh_size[i], replace = T, prob = start_prob),
                                        hh_size = rep(hh_size[i], hh_size[i]),
                                        hh_id = rep(i, hh_size[i])))

        complete_obs <- complete_obs %>%
          bind_rows(new_obs)

        prior <- new_obs$state

      } else {
        prior_inf <- sum(prior == 2)
        new_states <- rep(0, hh_size[i])
        for(part in 1:hh_size[i]) {
          if(prior[part] == 1) {
            eh_prob_x <- inv_logit(logit(eh_prob)+sum(x[last_x+part,]*covs_eh))
            ih_prob_x <- inv_logit(logit(ih_prob)+sum(x[last_x+part,]*covs_ih))
            no_inf_prob <- (1-eh_prob_x)*(1-ih_prob_x)^prior_inf
            new_states[part] = sample(x = c(1, 2, 3),
                                      size = 1,
                                      prob = c(no_inf_prob,
                                               (1-no_inf_prob),
                                               0))
          } else if(prior[part] == 2) {
            new_states[part] = sample(x = c(1, 2, 3),
                                      size = 1,
                                      prob = c(0,
                                               1-gamma,
                                               gamma))
          } else {
            new_states[part] = sample(x = c(1, 2, 3),
                                      size = 1,
                                      prob = c(0,
                                               0,
                                               1))
          }
        }

        new_obs <- bind_rows(data.frame(t = rep(d, hh_size[i]),
                                        part_id = part_ids[[i]],
                                        enroll = c(rep(1, enroll_per_hh[i]), rep(0, hh_size[i]-enroll_per_hh[i])),
                                        state = new_states,
                                        hh_size = rep(hh_size[i], hh_size[i]),
                                        hh_id = rep(i, hh_size[i])))

        complete_obs <- complete_obs %>%
          bind_rows(new_obs)

        prior <- new_states
      }
    }
    last_x <- last_x + hh_size[i]
  }

  complete_obs <- complete_obs %>%
    arrange(hh_id, t, part_id)

  outcome <- matrix(nrow = nrow(complete_obs), ncol = length(obs_prob))
  outcome_names <- paste0("y", 1:length(obs_prob))
  colnames(outcome) <- outcome_names
  for(i in 1:nrow(complete_obs)) {
    for(j in 1:length(obs_prob)) {
      p1 <- obs_prob[[j]][complete_obs$state[i]]
      outcome[i,j] <- sample(c(0,1), 1, prob = c(1-p1, p1))
    }
  }

  complete_obs <- complete_obs %>%
    bind_cols(as.data.frame(outcome))

  if(!all(c(covs_eh, covs_ih) == 0)) {
    x <- as.data.frame(x)
    names(x) <- paste0("x", 1:ncol(x))
  }

  if(!complete_enroll) {
    obs <- complete_obs %>% filter(enroll == 1)
  } else {
    obs <- complete_obs
  }
  
  out <- list(obs = obs,
              complete_obs = complete_obs)
  
  if(!all(c(covs_eh, covs_ih) == 0)) {
    out <- append(out,
                  list(x = x))
  }

  return(out)
}


sim_siir <- function(eh_prob = 0.01, ih_prob = 0.05, n_hh = 100,
                     hh_size = 1:5, tmax = 100, gamma = c(1/5, 1/30), split = c(0.7, 0.3),
                     covs_eh = c(0, 0), covs_ih = c(0, 0),
                     obs_prob = list(c(0.05, 0.95, 0.95, 0.05),
                                     c(0.01, 0.01, 0.01, 0.8)),
                     start_prob = c(1, 0, 0, 0),
                     complete_enroll = TRUE) {

  if(length(ih_prob) == 1) {
    ih_prob <- rep(ih_prob, 2)
  }

  epsilon <- 1e-10

  hh_size <- sample(hh_size, n_hh, replace = TRUE) # household sizes
  enroll_per_hh <- numeric(n_hh) # number of participants enrolled per HH

  x <- matrix(nrow = sum(hh_size), ncol = length(covs_ih))

  for(i in 1:length(covs_ih)) {
    x[,i] <- rbinom(sum(hh_size), 1, 0.4)
  }

  # Create participant IDs
  part_ids <- list()
  for(i in 1:n_hh) {
    part_ids[[i]] <- 1:(hh_size[i])
    if(hh_size[i] == 1) {
      enroll_per_hh[i] <- 1
    } else {
      enroll_per_hh[i] <- sample(1:hh_size[i], 1)
    }
  }

  epsilon <- 1e-10
  enroll_ids <- list()

  # Infection state for all household members on all time steps
  complete_obs <- data.frame(t = numeric(),
                             part_id = numeric(),
                             enroll = numeric(),
                             state = numeric(),
                             hh_size = numeric(),
                             hh_id = numeric())

  last_x <- 0

  for(i in 1:length(hh_size)) {

    # Move HH members through SIR states
    for(d in 1:tmax) {
      wk <- base::ceiling(d/7)
      if(d == 1) {
        new_obs <- bind_rows(data.frame(t = rep(d, hh_size[i]),
                                        part_id = part_ids[[i]],
                                        enroll = c(rep(1, enroll_per_hh[i]), rep(0, hh_size[i]-enroll_per_hh[i])),
                                        state = sample(1:4, hh_size[i], replace = T, prob = start_prob),
                                        hh_size = rep(hh_size[i], hh_size[i]),
                                        hh_id = rep(i, hh_size[i])))

        complete_obs <- complete_obs %>%
          bind_rows(new_obs)

        prior <- new_obs$state

      } else {
        prior_inf <- c(sum(prior == 2), sum(prior == 3))
        new_states <- rep(0, hh_size[i])
        for(part in 1:hh_size[i]) {
          if(prior[part] == 1) {
            eh_prob_x <- inv_logit(logit(eh_prob)+sum(x[last_x+part,]*covs_eh))
            ih_prob_x <- inv_logit(logit(ih_prob)+sum(x[last_x+part,]*covs_ih))
            no_inf_prob <- (1-eh_prob_x)*prod((1-ih_prob_x)^prior_inf)
            new_states[part] = sample(x = c(1, 2, 3, 4),
                                      size = 1,
                                      prob = c(no_inf_prob,
                                               (1-no_inf_prob)*split[1],
                                               (1-no_inf_prob)*split[2],
                                               0))
          } else if(prior[part] == 2) {
            new_states[part] = sample(x = c(1, 2, 3, 4),
                                      size = 1,
                                      prob = c(0,
                                               1-gamma[1],
                                               0,
                                               gamma[1]))
          } else if(prior[part] == 3) {
            new_states[part] = sample(x = c(1, 2, 3, 4),
                                      size = 1,
                                      prob = c(0,
                                               0,
                                               1-gamma[2],
                                               gamma[2]))
          } else {
            new_states[part] = sample(x = c(1, 2, 3, 4),
                                      size = 1,
                                      prob = c(0,
                                               0,
                                               0,
                                               1))
          }
        }

        new_obs <- bind_rows(data.frame(t = rep(d, hh_size[i]),
                                        part_id = part_ids[[i]],
                                        enroll = c(rep(1, enroll_per_hh[i]), rep(0, hh_size[i]-enroll_per_hh[i])),
                                        state = new_states,
                                        hh_size = rep(hh_size[i], hh_size[i]),
                                        hh_id = rep(i, hh_size[i])))

        complete_obs <- complete_obs %>%
          bind_rows(new_obs)

        prior <- new_states
      }
    }
    last_x <- last_x + hh_size[i]
  }

  complete_obs <- complete_obs %>%
    arrange(hh_id, t, part_id)

  outcome <- matrix(nrow = nrow(complete_obs), ncol = length(obs_prob))
  outcome_names <- paste0("y", 1:length(obs_prob))
  colnames(outcome) <- outcome_names
  for(i in 1:nrow(complete_obs)) {
    for(j in 1:length(obs_prob)) {
      p1 <- obs_prob[[j]][complete_obs$state[i]]
      outcome[i,j] <- sample(c(0,1), 1, prob = c(1-p1, p1))
    }
  }

  complete_obs <- complete_obs %>%
    bind_cols(as.data.frame(outcome))

  if(!complete_enroll) {
    obs <- complete_obs %>% filter(enroll == 1)
  } else {
    obs <- complete_obs
  }

  return(list(obs = obs,
              complete_obs = complete_obs))
}

sim_seir <- function(eh_prob = 0.01, ih_prob = 0.05, n_hh = 100,
                    hh_size = 1:5, tmax = 100, sigma = 1/2, gamma = 1/5,
                    covs_eh = c(0, 0), covs_ih = c(0, 0),
                    obs_prob = list(c(0.05, 0.05, 0.95, 0.05),
                                    c(0.01, 0.01, 0.1, 0.8)),
                    start_prob = c(1, 0, 0, 0),
                    complete_enroll = TRUE) {

  epsilon <- 1e-10

  hh_size <- sample(hh_size, n_hh, replace = TRUE) # household sizes
  enroll_per_hh <- numeric(n_hh) # number of participants enrolled per HH

  x <- matrix(nrow = sum(hh_size), ncol = length(covs_ih))

  for(i in 1:length(covs_ih)) {
    x[,i] <- rbinom(sum(hh_size), 1, 0.4)
  }

  # Create participant IDs
  part_ids <- list()
  for(i in 1:n_hh) {
    part_ids[[i]] <- 1:(hh_size[i])
    if(hh_size[i] == 1) {
      enroll_per_hh[i] <- 1
    } else {
      enroll_per_hh[i] <- sample(1:hh_size[i], 1)
    }
  }

  epsilon <- 1e-10
  enroll_ids <- list()

  # Infection state for all household members on all time steps
  complete_obs <- data.frame(t = numeric(),
                             part_id = numeric(),
                             enroll = numeric(),
                             state = numeric(),
                             hh_size = numeric(),
                             hh_id = numeric())

  last_x <- 0

  for(i in 1:length(hh_size)) {

    # Move HH members through SIR states
    for(d in 1:tmax) {
      wk <- base::ceiling(d/7) # @ Claire is this doing anything?
      if(d == 1) {
        new_obs <- bind_rows(data.frame(t = rep(d, hh_size[i]),
                                        part_id = part_ids[[i]],
                                        enroll = c(rep(1, enroll_per_hh[i]), rep(0, hh_size[i]-enroll_per_hh[i])),
                                        state = sample(1:4, hh_size[i], replace = T, prob = start_prob),
                                        hh_size = rep(hh_size[i], hh_size[i]),
                                        hh_id = rep(i, hh_size[i])))

        complete_obs <- complete_obs %>%
          bind_rows(new_obs)

        prior <- new_obs$state

      } else {
        prior_inf <- sum(prior == 3)
        new_states <- rep(0, hh_size[i])
        for(part in 1:hh_size[i]) {
          if(prior[part] == 1) {
            eh_prob_x <- inv_logit(logit(eh_prob)+sum(x[last_x+part,]*covs_eh))
            ih_prob_x <- inv_logit(logit(ih_prob)+sum(x[last_x+part,]*covs_ih))
            no_inf_prob <- (1-eh_prob_x)*(1-ih_prob_x)^prior_inf
            new_states[part] = sample(x = c(1, 2, 3, 4),
                                      size = 1,
                                      prob = c(no_inf_prob,
                                               (1-no_inf_prob),
                                               0,
                                               0))
          } else if(prior[part] == 2) {
            new_states[part] = sample(x = c(1, 2, 3, 4),
                                      size = 1,
                                      prob = c(0,
                                               1-sigma,
                                               sigma,
                                               0))
          } else if(prior[part] == 3) {
            new_states[part] = sample(x = c(1, 2, 3, 4),
                                      size = 1,
                                      prob = c(0,
                                               0,
                                               1-gamma,
                                               gamma))
          } else {
            new_states[part] = sample(x = c(1, 2, 3, 4),
                                      size = 1,
                                      prob = c(0,
                                               0,
                                               0,
                                               1))
          }
        }

        new_obs <- bind_rows(data.frame(t = rep(d, hh_size[i]),
                                        part_id = part_ids[[i]],
                                        enroll = c(rep(1, enroll_per_hh[i]), rep(0, hh_size[i]-enroll_per_hh[i])),
                                        state = new_states,
                                        hh_size = rep(hh_size[i], hh_size[i]),
                                        hh_id = rep(i, hh_size[i])))

        complete_obs <- complete_obs %>%
          bind_rows(new_obs)

        prior <- new_states
      }
    }
    last_x <- last_x + hh_size[i]
  }

  complete_obs <- complete_obs %>%
    arrange(hh_id, t, part_id)

  outcome <- matrix(nrow = nrow(complete_obs), ncol = length(obs_prob))
  outcome_names <- paste0("y", 1:length(obs_prob))
  colnames(outcome) <- outcome_names
  for(i in 1:nrow(complete_obs)) {
    for(j in 1:length(obs_prob)) {
      p1 <- obs_prob[[j]][complete_obs$state[i]]
      outcome[i,j] <- sample(c(0,1), 1, prob = c(1-p1, p1))
    }
  }

  complete_obs <- complete_obs %>%
    bind_cols(as.data.frame(outcome))

  if(!complete_enroll) {
    obs <- complete_obs %>% filter(enroll == 1)
  } else {
    obs <- complete_obs
  }

  return(list(obs = obs,
              complete_obs = complete_obs))
}

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
### SIS model ###@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

sim_sis <- function(eh_prob = 0.05, 
                    ih_prob = 0.05, 
                    n_hh = 100,
                    hh_size = 1:5, 
                    tmax = 100, 
                    gamma = 1/3,
                    covs_eh = c(0, 0), 
                    covs_ih = c(0, 0),
                    covs_gamma = NULL,            # recovery effects on logit(gamma), length == length(covs_eh)
                    cat_covs = NULL,               # list of categorical covariate specs
                    obs_prob = list(c(0.05, 0.95), # S&+ , I&+
                                    c(0.01, 0.8)), # S&+ , I&+
                    start_prob = c(0.8, 0.2),
                    complete_enroll = TRUE,
                    sigma_pid = 0,                   # standard deviation by individual on logit scale
                    season_type = c("fourier", "spline"),
                    season_k = 0,                    # number of Fourier harmonics
                    season_period = 365,             # period in days
                    season_df = 6,                   # spline basis dimension (mgcv::s k)
                    season_knots = NULL,             # optional knots list for mgcv::s (e.g., list(day = c(0.5, 365.5)))
                    season_coefs_eh = NULL,          # length 2*season_k (sin1, cos1, sin2, cos2, ...)
                    season_coefs_ih = NULL,          # length 2*season_k (sin1, cos1, sin2, cos2, ...)
                    viral_acq = NULL,                # length 1 or viral_period: prob of viral acquisition by day-of-year
                    viral_period = season_period,    # period for viral acquisition probabilities
                    viral_eps = 1/7,                 # recovery probability per day for viral infection
                    viral_init_prob = 0,             # initial viral prevalence at day 1
                    viral_trans_boost = 0) {         # log-odds boost to ih_prob for viral-infected transmitters
  
  epsilon <- 1e-10 # tiny value to ensure no true zeros

  if (!is.numeric(gamma) || length(gamma) != 1 || is.na(gamma) || gamma < 0 || gamma > 1) {
    stop("gamma must be a single probability in [0, 1].")
  }
  
  if (!is.null(viral_acq)) {
    if (!is.numeric(viral_acq)) {
      stop("viral_acq must be numeric (probabilities).")
    }
    if (any(is.na(viral_acq)) || any(viral_acq < 0 | viral_acq > 1)) {
      stop("viral_acq values must be in [0, 1] with no NAs.")
    }
    if (!(length(viral_acq) == 1 || length(viral_acq) == viral_period)) {
      stop("viral_acq must have length 1 (fixed risk) or length of the viral_period for time variation.")
    }
  }
  if (!is.numeric(viral_eps) || length(viral_eps) != 1 || is.na(viral_eps) || viral_eps < 0 || viral_eps > 1) {
    stop("viral_eps must be a single probability in [0, 1].")
  }
  if (!is.numeric(viral_init_prob) || length(viral_init_prob) != 1 || is.na(viral_init_prob) ||
      viral_init_prob < 0 || viral_init_prob > 1) {
    stop("viral_init_prob must be a single probability in [0, 1].")
  }

  season_type <- match.arg(season_type)

  if (season_type == "fourier") {
    season_effect <- function(d, coefs, period) {
      if (is.null(coefs) || length(coefs) == 0) return(0)
      if ((length(coefs) %% 2) != 0) stop("season_coefs must have even length (sin, cos pairs).")
      k <- length(coefs) / 2
      harm <- seq_len(k)
      sin_terms <- sin(2 * pi * harm * d / period)
      cos_terms <- cos(2 * pi * harm * d / period)
      sum(coefs[seq(1, length(coefs), by = 2)] * sin_terms +
            coefs[seq(2, length(coefs), by = 2)] * cos_terms)
    }
    if (season_k > 0) {
      if (is.null(season_coefs_eh)) season_coefs_eh <- rep(0, 2 * season_k)
      if (is.null(season_coefs_ih)) season_coefs_ih <- rep(0, 2 * season_k)
      if (!is.null(season_coefs_eh) && length(season_coefs_eh) != 2 * season_k) {
        stop("season_coefs_eh must have length 2*season_k.")
      }
      if (!is.null(season_coefs_ih) && length(season_coefs_ih) != 2 * season_k) {
        stop("season_coefs_ih must have length 2*season_k.")
      }
    }
  }

  if (season_type == "spline") {
    day_seq <- seq_len(season_period)
    knots <- season_knots
    if (is.null(knots)) {
      knots <- list(day = c(0.5, season_period + 0.5))
    }
    smooth <- mgcv::smoothCon(
      mgcv::s(day, bs = "cc", k = season_df),
      data = data.frame(day = day_seq),
      knots = knots,
      absorb.cons = TRUE
    )
    season_basis <- smooth[[1]]$X
    if (!is.matrix(season_basis)) season_basis <- as.matrix(season_basis)
    n_basis <- ncol(season_basis)

    if (is.null(season_coefs_eh)) season_coefs_eh <- rep(0, n_basis)
    if (is.null(season_coefs_ih)) season_coefs_ih <- rep(0, n_basis)
    if (length(season_coefs_eh) != n_basis) {
      stop("season_coefs_eh must have length n_basis (from mgcv spline basis).")
    }
    if (length(season_coefs_ih) != n_basis) {
      stop("season_coefs_ih must have length n_basis (from mgcv spline basis).")
    }

    season_effect <- function(d, coefs, period) {
      if (is.null(coefs) || length(coefs) == 0) return(0)
      doy <- ((d - 1) %% period) + 1
      sum(season_basis[doy, ] * coefs)
    }
  }
  
  hh_size <- sample(hh_size, n_hh, replace = TRUE) # household sizes
  enroll_per_hh <- numeric(n_hh) # number of participants enrolled per HH ??
  
  n_people <- sum(hh_size)

  if (length(covs_eh) != length(covs_ih)) {
    stop("covs_eh and covs_ih must have the same length.")
  }
  if (is.null(covs_gamma)) covs_gamma <- rep(0, length(covs_eh))
  if (length(covs_gamma) != length(covs_eh)) {
    stop("covs_gamma must have the same length as covs_eh.")
  }

  build_cat_covs <- function(n, cat_covs_spec) {
    if (is.null(cat_covs_spec) || length(cat_covs_spec) == 0) {
      return(list(
        x = NULL,
        colnames = character(0),
        coefs_eh = numeric(0),
        coefs_ih = numeric(0),
        coefs_gamma = numeric(0),
        raw = NULL
      ))
    }
    if (!is.list(cat_covs_spec)) {
      stop("cat_covs must be a list of covariate specs.")
    }

    x_list <- list()
    colnames_out <- character(0)
    coefs_eh_out <- numeric(0)
    coefs_ih_out <- numeric(0)
    coefs_gamma_out <- numeric(0)
    # Pre-size to avoid "replacement has n rows, data has 0" when adding columns.
    raw_df <- as.data.frame(matrix(nrow = n, ncol = 0))

    for (idx in seq_along(cat_covs_spec)) {
      spec <- cat_covs_spec[[idx]]
      name <- spec$name
      if (is.null(name) || !nzchar(name)) {
        name <- paste0("cat", idx)
      }

      levels <- spec$levels
      if (is.null(levels) || length(levels) < 2) {
        stop("Each categorical covariate must define at least 2 levels.")
      }
      levels <- as.character(levels)

      probs <- spec$probs
      if (is.null(probs)) {
        probs <- rep(1 / length(levels), length(levels))
      }
      if (length(probs) != length(levels)) {
        stop("cat_covs$probs length must match levels length.")
      }
      probs <- probs / sum(probs)

      baseline <- spec$baseline
      if (is.null(baseline)) baseline <- levels[1]
      baseline <- as.character(baseline)
      if (!baseline %in% levels) {
        stop("cat_covs$baseline must be one of levels.")
      }

      cat_vals <- sample(levels, n, replace = TRUE, prob = probs)
      raw_df[[name]] <- cat_vals

      lev_no_base <- levels[levels != baseline]
      if (length(lev_no_base) < 1) {
        stop("Categorical covariate must have at least one non-baseline level.")
      }

      mat <- vapply(lev_no_base, function(lv) as.integer(cat_vals == lv), integer(n))
      if (is.null(dim(mat))) {
        mat <- matrix(mat, ncol = 1)
      }
      colnames(mat) <- paste0(name, "_", lev_no_base)

      ce <- spec$coefs_eh
      ci <- spec$coefs_ih
      cg <- spec$coefs_gamma
      n_coef <- length(lev_no_base)
      if (is.null(ce)) ce <- rep(0, n_coef)
      if (is.null(ci)) ci <- rep(0, n_coef)
      if (is.null(cg)) cg <- rep(0, n_coef)
      if (length(ce) != n_coef) {
        stop("cat_covs$coefs_eh length must match non-baseline levels.")
      }
      if (length(ci) != n_coef) {
        stop("cat_covs$coefs_ih length must match non-baseline levels.")
      }
      if (length(cg) != n_coef) {
        stop("cat_covs$coefs_gamma length must match non-baseline levels.")
      }

      x_list[[length(x_list) + 1]] <- mat
      colnames_out <- c(colnames_out, colnames(mat))
      coefs_eh_out <- c(coefs_eh_out, ce)
      coefs_ih_out <- c(coefs_ih_out, ci)
      coefs_gamma_out <- c(coefs_gamma_out, cg)
    }

    x_mat <- if (length(x_list) == 0) NULL else do.call(cbind, x_list)
    list(
      x = x_mat,
      colnames = colnames_out,
      coefs_eh = coefs_eh_out,
      coefs_ih = coefs_ih_out,
      coefs_gamma = coefs_gamma_out,
      raw = raw_df
    )
  }

  x_bin <- matrix(nrow = n_people, ncol = length(covs_ih))
  if (length(covs_ih) > 0) {
    for (i in 1:length(covs_ih)) {
      x_bin[, i] <- rbinom(n_people, 1, 0.4) # 4/10 get 1, 0
    }
  }

  cat_info <- build_cat_covs(n_people, cat_covs)
  x_cat <- cat_info$x

  x_parts <- list()
  x_colnames <- character(0)
  if (ncol(x_bin) > 0) {
    x_parts[[length(x_parts) + 1]] <- x_bin
    x_colnames <- c(x_colnames, paste0("x", 1:ncol(x_bin)))
  }
  if (!is.null(x_cat)) {
    x_parts[[length(x_parts) + 1]] <- x_cat
    x_colnames <- c(x_colnames, cat_info$colnames)
  }
  if (length(x_parts) == 0) {
    x <- matrix(nrow = n_people, ncol = 0)
  } else {
    x <- do.call(cbind, x_parts)
  }

  covs_eh_all <- c(covs_eh, cat_info$coefs_eh)
  covs_ih_all <- c(covs_ih, cat_info$coefs_ih)
  covs_gamma_all <- c(covs_gamma, cat_info$coefs_gamma)

  if (ncol(x) != length(covs_eh_all) || ncol(x) != length(covs_ih_all) || ncol(x) != length(covs_gamma_all)) {
    stop("Covariate matrix columns must match covs_eh / covs_ih lengths.")
  }

  u_pid <- rnorm(n_people, mean = 0, sd = sigma_pid)
  
  # Create participant IDs
  part_ids <- list()
  for(i in 1:n_hh) {o0§
    part_ids[[i]] <- 1:(hh_size[i])
    if(hh_size[i] == 1) {
      enroll_per_hh[i] <- 1
    } else {
      enroll_per_hh[i] <- sample(1:hh_size[i], 1)
    }
  }
  
  enroll_ids <- list()
  
  # Infection state for all household members on all time steps
  complete_obs <- data.frame(t = numeric(),
                             part_id = numeric(),
                             pid_global = numeric(),
                             enroll = numeric(),
                             state = numeric(),
                             viral = numeric(),
                             hh_size = numeric(),
                             hh_id = numeric(),
                             u_pid = numeric())
  
  last_x <- 0
  
  for(i in 1:length(hh_size)) {
    pid_global_vec <- last_x + seq_len(hh_size[i])
    
    # Move HH members through SIS states
    for(d in 1:tmax) {
      wk <- base::ceiling(d/7)
      if(d == 1) {

        viral_prior <- rbinom(hh_size[i], 1, viral_init_prob)

        new_obs <- bind_rows(data.frame(t = rep(d, hh_size[i]),
                                        part_id = part_ids[[i]],
                                        pid_global = pid_global_vec,
                                        enroll = c(rep(1, enroll_per_hh[i]), rep(0, hh_size[i]-enroll_per_hh[i])),
                                        state = sample(1:2, hh_size[i], replace = T, prob = start_prob),
                                        viral = viral_prior,
                                        hh_size = rep(hh_size[i], hh_size[i]),
                                        hh_id = rep(i, hh_size[i]),
                                        u_pid = u_pid[pid_global_vec]
                                        ))
        complete_obs <- complete_obs %>%
          bind_rows(new_obs)
        
        prior <- new_obs$state
        
      } else {
        prior_inf <- sum(prior == 2)
        prior_inf_viral <- sum(prior == 2 & viral_prior == 1)
        prior_inf_nonviral <- prior_inf - prior_inf_viral
        new_states <- integer(hh_size[i])
        viral_new <- integer(hh_size[i])
        if (!is.null(viral_acq)) {
          if (length(viral_acq) == 1) {
            p_acq <- viral_acq
          } else {
            doy_v <- ((d - 1) %% viral_period) + 1
            p_acq <- viral_acq[doy_v]
          }
        } else {
          p_acq <- 0
        }
        
        for (part in 1:hh_size[i]) {
          if (prior[part] == 1) {  # S state

            eh_prob_x <- plogis(qlogis(eh_prob) +
                                  sum(x[last_x + part, ] * covs_eh_all) +
                                  u_pid[last_x + part] +                            # individual level variability
                                  season_effect(d, season_coefs_eh, season_period)) # seasonality variability 

            ih_prob_x <- plogis(qlogis(ih_prob) +
                                  sum(x[last_x + part, ] * covs_ih_all) +
                                  u_pid[last_x + part] +                             # individual level variability
                                  season_effect(d, season_coefs_ih, season_period))  # seasonality variability

            ih_prob_viral <- plogis(qlogis(ih_prob_x) + viral_trans_boost)
            
            no_inf_prob <- (1 - eh_prob_x) *
              (1 - ih_prob_x)^prior_inf_nonviral *
              (1 - ih_prob_viral)^prior_inf_viral
            
            new_states[part] <- sample(x = c(1, 2), size = 1,
                                       prob = c(no_inf_prob, 1 - no_inf_prob))
          } else {                 # I state
            gamma_x <- plogis(qlogis(gamma) + sum(x[last_x + part, ] * covs_gamma_all))
            new_states[part] <- sample(x = c(1, 2), size = 1,
                                       prob = c(gamma_x, 1 - gamma_x))
          }

          if (viral_prior[part] == 1) {
            viral_new[part] <- sample(x = c(0, 1), size = 1, prob = c(viral_eps, 1 - viral_eps))
          } else {
            viral_new[part] <- rbinom(1, 1, p_acq)
          }
        }
        
        new_obs <- data.frame(
          t = rep(d, hh_size[i]),
          part_id = part_ids[[i]],
          pid_global = pid_global_vec,
          enroll = c(rep(1, enroll_per_hh[i]), rep(0, hh_size[i] - enroll_per_hh[i])),
          state = new_states,
          viral = viral_new,
          hh_size = rep(hh_size[i], hh_size[i]),
          hh_id = rep(i, hh_size[i]),
          u_pid = u_pid[pid_global_vec]
                             )
        
        complete_obs <- bind_rows(complete_obs, new_obs)
        prior <- new_states
        viral_prior <- viral_new
      }
    }
    last_x <- last_x + hh_size[i]
  }
  
  complete_obs <- complete_obs %>%
    arrange(hh_id, t, part_id)
  
  outcome <- matrix(nrow = nrow(complete_obs), ncol = length(obs_prob))
  outcome_names <- paste0("y", 1:length(obs_prob))
  colnames(outcome) <- outcome_names
  for(i in 1:nrow(complete_obs)) {
    for(j in 1:length(obs_prob)) {
      p1 <- obs_prob[[j]][complete_obs$state[i]]
      outcome[i,j] <- sample(c(0,1), 1, prob = c(1-p1, p1))
    }
  }
  
  complete_obs <- complete_obs %>%
    bind_cols(as.data.frame(outcome))
  
  if (ncol(x) > 0) {
    x <- as.data.frame(x)
    if (length(x_colnames) == ncol(x)) {
      names(x) <- x_colnames
    } else {
      names(x) <- paste0("x", 1:ncol(x))
    }
    x <- cbind(pid_global = seq_len(n_people), x)
  }
  
  if(!complete_enroll) {
    obs <- complete_obs %>% filter(enroll == 1)
  } else {
    obs <- complete_obs
  }
  
  out <- list(obs = obs,
              complete_obs = complete_obs)
  
  if (ncol(x) > 0) {
    out <- append(out,
                  list(x = x))
  }
  
  return(out)
}
