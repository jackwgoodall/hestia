
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
#@@@@@@@@@@@@@@@@@@@ SIS model ###@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

sim_sis <- function(eh_prob = 0.05,                  # Extra-household risk of infection per unit time
                    ih_prob = 0.05,                  # Intra-household risk of infection per infected household member per unit time
                    n_hh = 100,                      # Household numbers
                    hh_size = 1:5,                   # Range of household sizes
                    tmax = 100,                      # Time range
                    gamma = 1/3,                     # Recovery rate from the I state
                    covs_eh = c(0, 0), 
                    covs_ih = c(0, 0),
                    covs_prev = 0.4,                 # prevalence for binary covariates in covs_eh/covs_ih (scalar or length(covs_ih))
                    covs_gamma = NULL,               # recovery effects on logit(gamma), length == length(covs_eh)
                    cat_covs = NULL,                 # list of categorical covariate specifications 
                    obs_prob = list(c(0.05, 0.95),   # Specificity, Sensitivity (test 1)
                                    c(0.01, 0.8)),   # Specificity, Sensitivity (test 2)
                    start_prob = c(0.8, 0.2),        # Starting states (must sum to 1)
                    complete_enroll = TRUE,          
                    sigma_pid = 0,                   # standard deviation by individual on logit scale
                    season_type = c("fourier", "spline"),
                    season_k = 0,                    # number of Fourier harmonics
                    season_period = 365,             # period in days
                    season_df = 6,                   # spline basis dimension (mgcv::s k)
                    season_knots = NULL,             # optional knots list for mgcv::s (e.g., list(day = c(0.5, 365.5)))
                    season_coefs_eh = NULL,          # length 2*season_k (sin1, cos1, sin2, cos2, ...)
                    season_coefs_ih = NULL) {        # length 2*season_k (sin1, cos1, sin2, cos2, ...)
  
  epsilon <- 1e-10 # tiny value to ensure no true zeros

  if (!is.numeric(gamma) || length(gamma) != 1 || is.na(gamma) || gamma < 0 || gamma > 1) {
    stop("gamma must be a single probability in [0, 1].")
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
  if (!is.numeric(covs_prev) || any(is.na(covs_prev)) || any(covs_prev < 0 | covs_prev > 1)) {
    stop("covs_prev must contain probabilities in [0, 1].")
  }
  if (length(covs_ih) == 0) {
    if (!(length(covs_prev) %in% c(0, 1))) {
      stop("covs_prev must be length 1 (or empty) when no binary covariates are specified.")
    }
    covs_prev <- numeric(0)
  } else if (length(covs_prev) == 1) {
    covs_prev <- rep(covs_prev, length(covs_ih))
  } else if (length(covs_prev) != length(covs_ih)) {
    stop("covs_prev must have length 1 or length(covs_ih).")
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
      x_bin[, i] <- rbinom(n_people, 1, covs_prev[i])
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
  for(i in 1:n_hh) {
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
                             hh_size = numeric(),
                             hh_id = numeric(),
                             u_pid = numeric())
  
  last_x <- 0
  
  for(i in 1:length(hh_size)) {
    pid_global_vec <- last_x + seq_len(hh_size[i])
    
    # Move HH members through SIS states
    for(d in 1:tmax) {
      
      if(d == 1) {
        new_obs <- bind_rows(data.frame(t = rep(d, hh_size[i]),
                                        part_id = part_ids[[i]],
                                        pid_global = pid_global_vec,
                                        enroll = c(rep(1, enroll_per_hh[i]), rep(0, hh_size[i]-enroll_per_hh[i])),
                                        state = sample(1:2, hh_size[i], replace = T, prob = start_prob),
                                        hh_size = rep(hh_size[i], hh_size[i]),
                                        hh_id = rep(i, hh_size[i]),
                                        u_pid = u_pid[pid_global_vec]
                                        ))
        complete_obs <- complete_obs %>%
          bind_rows(new_obs)
        
        prior <- new_obs$state
        
      } else {
        prior_inf <- sum(prior == 2)
        new_states <- integer(hh_size[i])
        
        for (part in 1:hh_size[i]) {
          if (prior[part] == 1) {  

        # S state ---> 
            eh_prob_x <- plogis(qlogis(eh_prob) +
                                  sum(x[last_x + part, ] * covs_eh_all) +           # combined effect of all categorical covariates
                                  u_pid[last_x + part] +                            # individual level variability
                                  season_effect(d, season_coefs_eh, season_period)) # seasonality variability 

            ih_prob_x <- plogis(qlogis(ih_prob) +
                                  sum(x[last_x + part, ] * covs_ih_all) +            # combined effect of all categorical covariates
                                  u_pid[last_x + part] +                             # individual level variability
                                  season_effect(d, season_coefs_ih, season_period))  # seasonality variability

            no_inf_prob <- (1 - eh_prob_x) * (1 - ih_prob_x)^prior_inf
            
            new_states[part] <- sample(x = c(1, 2), size = 1,
                                       prob = c(no_inf_prob, 1 - no_inf_prob))
            
       # I state --->
          } else {                 
            gamma_x <- plogis(qlogis(gamma) + sum(x[last_x + part, ] * covs_gamma_all))
            new_states[part] <- sample(x = c(1, 2), size = 1,
                                       prob = c(gamma_x, 1 - gamma_x))
          }
        }
        
        new_obs <- data.frame(
          t = rep(d, hh_size[i]),
          part_id = part_ids[[i]],
          pid_global = pid_global_vec,
          enroll = c(rep(1, enroll_per_hh[i]), rep(0, hh_size[i] - enroll_per_hh[i])),
          state = new_states,
          hh_size = rep(hh_size[i], hh_size[i]),
          hh_id = rep(i, hh_size[i]),
          u_pid = u_pid[pid_global_vec]
                             )
        
        complete_obs <- bind_rows(complete_obs, new_obs)
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

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@   SIRS model   @@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

sim_sirs <- function(eh_prob = 0.05,
                     ih_prob = 0.05,
                     n_hh = 100,
                     hh_size = 1:5,
                     tmax = 100,
                     gamma = 1/3,
                     rho = 1/60,                          # Geometric waning immunity 
                     covs_eh = c(0, 0),
                     covs_ih = c(0, 0),
                     covs_prev = 0.4,                     # prevalence for binary covariates in covs_eh/covs_ih (scalar or length(covs_ih))
                     covs_gamma = NULL,
                     covs_rho = NULL,                     
                     cat_covs = NULL,
                     obs_prob = list(c(0.05, 0.95, 0.05),  # Test 1: S, I, R
                                     c(0.01, 0.01, 0.8)),  # Test 2: S, I, R
                     start_prob = c(0.75, 0.2, 0.05),      # S, I, R
                     complete_enroll = TRUE,
                     sigma_pid = 0,
                     season_type = c("fourier", "spline"),
                     season_k = 0,
                     season_period = 365,
                     season_df = 6,
                     season_knots = NULL,
                     season_coefs_eh = NULL,
                     season_coefs_ih = NULL) {

  if (!is.numeric(gamma) || length(gamma) != 1 || is.na(gamma) || gamma < 0 || gamma > 1) {
    stop("gamma must be a single probability in [0, 1].")
  }
  if (!is.numeric(rho) || length(rho) != 1 || is.na(rho) || rho < 0 || rho > 1) {
    stop("rho must be a single probability in [0, 1].")
  }
  if (length(start_prob) != 3 || any(start_prob < 0) || abs(sum(start_prob) - 1) > 1e-8) {
    stop("start_prob must have length 3 and sum to 1 for SIRS.")
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
      if (length(season_coefs_eh) != 2 * season_k) stop("season_coefs_eh must have length 2*season_k.")
      if (length(season_coefs_ih) != 2 * season_k) stop("season_coefs_ih must have length 2*season_k.")
    }
  } else {
    day_seq <- seq_len(season_period)
    knots <- season_knots
    if (is.null(knots)) knots <- list(day = c(0.5, season_period + 0.5))
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
    if (length(season_coefs_eh) != n_basis) stop("season_coefs_eh must match spline basis dimension.")
    if (length(season_coefs_ih) != n_basis) stop("season_coefs_ih must match spline basis dimension.")
    season_effect <- function(d, coefs, period) {
      if (is.null(coefs) || length(coefs) == 0) return(0)
      doy <- ((d - 1) %% period) + 1
      sum(season_basis[doy, ] * coefs)
    }
  }

  if (length(covs_eh) != length(covs_ih)) {
    stop("covs_eh and covs_ih must have the same length.")
  }
  if (is.null(covs_gamma)) covs_gamma <- rep(0, length(covs_eh))
  if (is.null(covs_rho)) covs_rho <- rep(0, length(covs_eh))
  if (length(covs_gamma) != length(covs_eh) || length(covs_rho) != length(covs_eh)) {
    stop("covs_gamma and covs_rho must have the same length as covs_eh.")
  }
  if (!is.numeric(covs_prev) || any(is.na(covs_prev)) || any(covs_prev < 0 | covs_prev > 1)) {
    stop("covs_prev must contain probabilities in [0, 1].")
  }
  if (length(covs_ih) == 0) {
    if (!(length(covs_prev) %in% c(0, 1))) {
      stop("covs_prev must be length 1 (or empty) when no binary covariates are specified.")
    }
    covs_prev <- numeric(0)
  } else if (length(covs_prev) == 1) {
    covs_prev <- rep(covs_prev, length(covs_ih))
  } else if (length(covs_prev) != length(covs_ih)) {
    stop("covs_prev must have length 1 or length(covs_ih).")
  }

  build_cat_covs <- function(n, cat_covs_spec) {
    if (is.null(cat_covs_spec) || length(cat_covs_spec) == 0) {
      return(list(
        x = NULL,
        colnames = character(0),
        coefs_eh = numeric(0),
        coefs_ih = numeric(0),
        coefs_gamma = numeric(0),
        coefs_rho = numeric(0)
      ))
    }
    x_list <- list()
    colnames_out <- character(0)
    coefs_eh_out <- numeric(0)
    coefs_ih_out <- numeric(0)
    coefs_gamma_out <- numeric(0)
    coefs_rho_out <- numeric(0)
    for (idx in seq_along(cat_covs_spec)) {
      spec <- cat_covs_spec[[idx]]
      name <- spec$name
      if (is.null(name) || !nzchar(name)) name <- paste0("cat", idx)
      levels <- as.character(spec$levels)
      if (length(levels) < 2) stop("Each categorical covariate must define at least 2 levels.")
      probs <- spec$probs
      if (is.null(probs)) probs <- rep(1 / length(levels), length(levels))
      if (length(probs) != length(levels)) stop("cat_covs$probs length must match levels length.")
      probs <- probs / sum(probs)
      baseline <- spec$baseline
      if (is.null(baseline)) baseline <- levels[1]
      baseline <- as.character(baseline)
      if (!baseline %in% levels) stop("cat_covs$baseline must be one of levels.")
      vals <- sample(levels, n, replace = TRUE, prob = probs)
      lev_no_base <- levels[levels != baseline]
      mat <- vapply(lev_no_base, function(lv) as.integer(vals == lv), integer(n))
      if (is.null(dim(mat))) mat <- matrix(mat, ncol = 1)
      colnames(mat) <- paste0(name, "_", lev_no_base)
      n_coef <- length(lev_no_base)
      ce <- spec$coefs_eh; if (is.null(ce)) ce <- rep(0, n_coef)
      ci <- spec$coefs_ih; if (is.null(ci)) ci <- rep(0, n_coef)
      cg <- spec$coefs_gamma; if (is.null(cg)) cg <- rep(0, n_coef)
      cr <- spec$coefs_rho; if (is.null(cr)) cr <- rep(0, n_coef)
      if (length(ce) != n_coef || length(ci) != n_coef || length(cg) != n_coef || length(cr) != n_coef) {
        stop("Categorical coefficient lengths must match non-baseline levels.")
      }
      x_list[[length(x_list) + 1]] <- mat
      colnames_out <- c(colnames_out, colnames(mat))
      coefs_eh_out <- c(coefs_eh_out, ce)
      coefs_ih_out <- c(coefs_ih_out, ci)
      coefs_gamma_out <- c(coefs_gamma_out, cg)
      coefs_rho_out <- c(coefs_rho_out, cr)
    }
    x_mat <- if (length(x_list) == 0) NULL else do.call(cbind, x_list)
    list(
      x = x_mat,
      colnames = colnames_out,
      coefs_eh = coefs_eh_out,
      coefs_ih = coefs_ih_out,
      coefs_gamma = coefs_gamma_out,
      coefs_rho = coefs_rho_out
    )
  }

  hh_size <- sample(hh_size, n_hh, replace = TRUE)
  n_people <- sum(hh_size)
  enroll_per_hh <- numeric(n_hh)
  x_bin <- matrix(nrow = n_people, ncol = length(covs_ih))
  if (length(covs_ih) > 0) {
    for (i in 1:length(covs_ih)) x_bin[, i] <- rbinom(n_people, 1, covs_prev[i])
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
  x <- if (length(x_parts) == 0) matrix(nrow = n_people, ncol = 0) else do.call(cbind, x_parts)

  covs_eh_all <- c(covs_eh, cat_info$coefs_eh)
  covs_ih_all <- c(covs_ih, cat_info$coefs_ih)
  covs_gamma_all <- c(covs_gamma, cat_info$coefs_gamma)
  covs_rho_all <- c(covs_rho, cat_info$coefs_rho)
  if (ncol(x) != length(covs_eh_all) || ncol(x) != length(covs_ih_all) ||
      ncol(x) != length(covs_gamma_all) || ncol(x) != length(covs_rho_all)) {
    stop("Covariate matrix columns must match coefficient lengths.")
  }

  u_pid <- rnorm(n_people, mean = 0, sd = sigma_pid)

  part_ids <- list()
  for (i in 1:n_hh) {
    part_ids[[i]] <- 1:(hh_size[i])
    enroll_per_hh[i] <- if (hh_size[i] == 1) 1 else sample(1:hh_size[i], 1)
  }

  complete_obs <- data.frame(
    t = numeric(),
    part_id = numeric(),
    pid_global = numeric(),
    enroll = numeric(),
    state = numeric(),
    hh_size = numeric(),
    hh_id = numeric(),
    u_pid = numeric()
  )

  last_x <- 0
  for (i in seq_along(hh_size)) {
    pid_global_vec <- last_x + seq_len(hh_size[i])
    for (d in 1:tmax) {
      if (d == 1) {
        new_obs <- data.frame(
          t = rep(d, hh_size[i]),
          part_id = part_ids[[i]],
          pid_global = pid_global_vec,
          enroll = c(rep(1, enroll_per_hh[i]), rep(0, hh_size[i] - enroll_per_hh[i])),
          state = sample(1:3, hh_size[i], replace = TRUE, prob = start_prob),
          hh_size = rep(hh_size[i], hh_size[i]),
          hh_id = rep(i, hh_size[i]),
          u_pid = u_pid[pid_global_vec]
        )
        complete_obs <- bind_rows(complete_obs, new_obs)
        prior <- new_obs$state
      } else {
        prior_inf <- sum(prior == 2)
        new_states <- integer(hh_size[i])
        for (part in 1:hh_size[i]) {
          idx <- last_x + part
          if (prior[part] == 1) {
            eh_prob_x <- plogis(qlogis(eh_prob) +
                                  sum(x[idx, ] * covs_eh_all) +
                                  u_pid[idx] +
                                  season_effect(d, season_coefs_eh, season_period))
            ih_prob_x <- plogis(qlogis(ih_prob) +
                                  sum(x[idx, ] * covs_ih_all) +
                                  u_pid[idx] +
                                  season_effect(d, season_coefs_ih, season_period))
            no_inf_prob <- (1 - eh_prob_x) * (1 - ih_prob_x)^prior_inf
            new_states[part] <- sample(c(1, 2, 3), 1, prob = c(no_inf_prob, 1 - no_inf_prob, 0))
          } else if (prior[part] == 2) {
            gamma_x <- plogis(qlogis(gamma) + sum(x[idx, ] * covs_gamma_all))
            new_states[part] <- sample(c(1, 2, 3), 1, prob = c(0, 1 - gamma_x, gamma_x))
          } else {
            rho_x <- plogis(qlogis(rho) + sum(x[idx, ] * covs_rho_all))
            new_states[part] <- sample(c(1, 2, 3), 1, prob = c(rho_x, 0, 1 - rho_x))
          }
        }
        new_obs <- data.frame(
          t = rep(d, hh_size[i]),
          part_id = part_ids[[i]],
          pid_global = pid_global_vec,
          enroll = c(rep(1, enroll_per_hh[i]), rep(0, hh_size[i] - enroll_per_hh[i])),
          state = new_states,
          hh_size = rep(hh_size[i], hh_size[i]),
          hh_id = rep(i, hh_size[i]),
          u_pid = u_pid[pid_global_vec]
        )
        complete_obs <- bind_rows(complete_obs, new_obs)
        prior <- new_states
      }
    }
    last_x <- last_x + hh_size[i]
  }

  complete_obs <- complete_obs %>% arrange(hh_id, t, part_id)
  outcome <- matrix(nrow = nrow(complete_obs), ncol = length(obs_prob))
  colnames(outcome) <- paste0("y", 1:length(obs_prob))
  for (i in 1:nrow(complete_obs)) {
    for (j in 1:length(obs_prob)) {
      p1 <- obs_prob[[j]][complete_obs$state[i]]
      outcome[i, j] <- sample(c(0, 1), 1, prob = c(1 - p1, p1))
    }
  }
  complete_obs <- complete_obs %>% bind_cols(as.data.frame(outcome))
  if (!complete_enroll) {
    obs <- complete_obs %>% filter(enroll == 1)
  } else {
    obs <- complete_obs
  }

  out <- list(obs = obs, complete_obs = complete_obs)
  if (ncol(x) > 0) {
    x <- as.data.frame(x)
    names(x) <- if (length(x_colnames) == ncol(x)) x_colnames else paste0("x", 1:ncol(x))
    out$x <- cbind(pid_global = seq_len(n_people), x)
  }
  return(out)
}

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
### SIS two-pathogen helper function ###@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

sim_sis_from_existing <- function(base_complete_obs,
                                  base_x = NULL,
                                  a_complete_obs = NULL,
                                  eh_prob = 0.05,
                                  ih_prob = 0.05,
                                  gamma = 1/3,
                                  covs_eh = NULL,
                                  covs_ih = NULL,
                                  covs_gamma = NULL,
                                  obs_prob = list(c(0.05, 0.95),
                                                  c(0.01, 0.8)),
                                  start_prob = c(0.8, 0.2),
                                  complete_enroll = TRUE,
                                  season_type = c("fourier", "spline"),
                                  season_k = 0,
                                  season_period = 365,
                                  season_df = 6,
                                  season_knots = NULL,
                                  season_coefs_eh = NULL,
                                  season_coefs_ih = NULL,
                                  a_infectious_states = 2,      # Define which of A's states are active (default is 2 - i.e. the I state)
                                  cross_eh_coef = 0,            # A effect on B EH susceptibility (logit scale)
                                  cross_ih_susc_coef = 0,       # A effect on B IH susceptibility (logit scale)
                                  cross_ih_trans_coef = 0) {    # A effect on B transmissible to infected household members (logit scale)

  req_cols <- c("t", "part_id", "pid_global", "enroll", "hh_size", "hh_id")
  if (!all(req_cols %in% names(base_complete_obs))) {
    stop("base_complete_obs must contain: t, part_id, pid_global, enroll, hh_size, hh_id.")
  }
  if (!is.numeric(gamma) || length(gamma) != 1 || is.na(gamma) || gamma < 0 || gamma > 1) {
    stop("gamma must be a single probability in [0, 1].")
  }

  # Pull out the base design
  base_design <- base_complete_obs %>%
    distinct(hh_id, hh_size, part_id, pid_global, enroll) %>%
    arrange(hh_id, part_id)
  hh_tbl <- base_design %>% distinct(hh_id, hh_size) %>% arrange(hh_id)
  hh_size <- hh_tbl$hh_size
  n_hh <- nrow(hh_tbl)
  tmax <- max(base_complete_obs$t)
  n_people <- nrow(base_design)

  # And covariates (from the x subdataframe)
  if (is.null(base_x)) {
    x <- matrix(nrow = n_people, ncol = 0)
    x_colnames <- character(0)
  } else {
    if (!("pid_global" %in% names(base_x))) stop("base_x must include pid_global.")
    x_df <- base_x %>%
      distinct(pid_global, .keep_all = TRUE) %>%
      arrange(pid_global)
    if (nrow(x_df) != n_people) stop("base_x must have one row per pid_global.")
    x <- as.matrix(x_df[, setdiff(names(x_df), "pid_global"), drop = FALSE])
    x_colnames <- colnames(x)
  }
  kx <- ncol(x)
  resolve_covs <- function(covs, arg_name) {
    if (is.null(covs)) return(rep(0, kx))   # i.e. no effect of covariates added 
    if (kx == 0) {
      if (length(covs) == 0) return(numeric(0))
      stop(paste0(arg_name, " must be NULL/empty when base_x has no covariate columns."))
    }
    if (!is.null(names(covs)) && any(nzchar(names(covs)))) {
      out <- rep(0, kx)
      names(out) <- x_colnames
      bad <- setdiff(names(covs), x_colnames)
      if (length(bad) > 0) {
        stop(paste0(arg_name, " has unknown covariate names: ", paste(bad, collapse = ", ")))
      }
      out[names(covs)] <- as.numeric(covs)
      return(unname(out))
    }
    if (length(covs) != kx) {
      stop(paste0(arg_name, " must be either a full length-", kx, " vector or a named subset matching base_x columns."))
    }
    as.numeric(covs)
  }
  covs_eh <- resolve_covs(covs_eh, "covs_eh")
  covs_ih <- resolve_covs(covs_ih, "covs_ih")
  covs_gamma <- resolve_covs(covs_gamma, "covs_gamma")

  # from main SIS framework
  season_type <- match.arg(season_type)
  if (season_type == "fourier") {
    season_effect <- function(d, coefs, period) {
      if (is.null(coefs) || length(coefs) == 0) return(0)
      if ((length(coefs) %% 2) != 0) stop("season_coefs must have even length.")
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
      if (length(season_coefs_eh) != 2 * season_k) stop("season_coefs_eh must have length 2*season_k.")
      if (length(season_coefs_ih) != 2 * season_k) stop("season_coefs_ih must have length 2*season_k.")
    }
  } else {
    day_seq <- seq_len(season_period)
    knots <- season_knots
    if (is.null(knots)) knots <- list(day = c(0.5, season_period + 0.5))
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
    if (length(season_coefs_eh) != n_basis || length(season_coefs_ih) != n_basis) {
      stop("season coefs must match spline basis dimension.")
    }
    season_effect <- function(d, coefs, period) {
      if (is.null(coefs) || length(coefs) == 0) return(0)
      doy <- ((d - 1) %% period) + 1
      sum(season_basis[doy, ] * coefs)
    }
  }

  u_vec <- base_complete_obs %>%
    filter(t == min(t)) %>%
    arrange(pid_global) %>%
    pull(u_pid)
  if (length(u_vec) != n_people || any(is.na(u_vec))) {
    u_vec <- rep(0, n_people)
  }

  if (!is.null(a_complete_obs)) {
    if (!all(c("pid_global", "t", "state") %in% names(a_complete_obs))) {
      stop("a_complete_obs must include pid_global, t, and state.")
    }
    
    # Determine whether they were in the infectious state at last t (converts to 0/1)
    a_prev <- a_complete_obs %>%
      transmute(pid_global = pid_global, t = t, a_active = as.integer(state %in% a_infectious_states))
  } else {
    a_prev <- expand.grid(pid_global = seq_len(n_people), t = seq_len(tmax)) %>%    # <-This is arguably a bit futile as the whole point
      mutate(a_active = 0L)                                                         # here is to make a function that will take one state and 
  }                                                                                 # interact with another - but might at some point try and
  a_map <- a_prev %>%                                                               # make this universal...
    mutate(key = paste(pid_global, t, sep = "_")) %>%
    select(key, a_active)
  a_active_lookup <- setNames(a_map$a_active, a_map$key)

  complete_obs <- data.frame(
    t = numeric(),
    part_id = numeric(),
    pid_global = numeric(),
    enroll = numeric(),
    state = numeric(),
    hh_size = numeric(),
    hh_id = numeric(),
    u_pid = numeric()
  )

  last_x <- 0
  for (i in seq_len(n_hh)) {
    hh_n <- hh_size[i]
    hh_members <- base_design %>%
      filter(hh_id == i) %>%
      arrange(part_id)
    pid_vec <- hh_members$pid_global
    enroll_vec <- hh_members$enroll
    part_vec <- hh_members$part_id

    for (d in seq_len(tmax)) {
      if (d == 1) {
        new_obs <- data.frame(
          t = rep(d, hh_n),
          part_id = part_vec,
          pid_global = pid_vec,
          enroll = enroll_vec,
          state = sample(1:2, hh_n, replace = TRUE, prob = start_prob),
          hh_size = rep(hh_n, hh_n),
          hh_id = rep(i, hh_n),
          u_pid = u_vec[pid_vec]
        )
        complete_obs <- bind_rows(complete_obs, new_obs)
        prior <- new_obs$state
      } else {
        new_states <- integer(hh_n) # vectors of 0s with length of hh_n to be filled later                   
        a_prev_vec <- sapply(pid_vec, function(pid) {
          val <- a_active_lookup[[paste(pid, d - 1, sep = "_")]]      # find the last t's state
          if (is.null(val)) 0L else val
        })
        prior_inf_active <- sum(prior == 2 & a_prev_vec == 1)       # previously a+ and new (b)+
        prior_inf_inactive <- sum(prior == 2 & a_prev_vec == 0)     # previously a- and new (b)+
        for (part in seq_len(hh_n)) {
          pid <- pid_vec[part]
          xrow <- if (kx > 0) x[pid, ] else numeric(0)
          if (prior[part] == 1) {
                    eh_lp <- qlogis(eh_prob) +
                      (if (kx > 0) sum(xrow * covs_eh) else 0) +
                      u_vec[pid] +
                      season_effect(d, season_coefs_eh, season_period) +
                      cross_eh_coef * a_prev_vec[part]                  # turns on/off the a risk 
                    ih_lp <- qlogis(ih_prob) +
                      (if (kx > 0) sum(xrow * covs_ih) else 0) +
                      u_vec[pid] +
                      season_effect(d, season_coefs_ih, season_period) +
                      cross_ih_susc_coef * a_prev_vec[part]
            eh_prob_x <- plogis(eh_lp)
            ih_prob_x <- plogis(ih_lp)
            ih_prob_x_active <- plogis(ih_lp + cross_ih_trans_coef)
            no_inf_prob <- (1 - eh_prob_x) *
              (1 - ih_prob_x)^prior_inf_inactive *
              (1 - ih_prob_x_active)^prior_inf_active
            new_states[part] <- sample(c(1, 2), 1, prob = c(no_inf_prob, 1 - no_inf_prob))
          } else {
            gamma_x <- plogis(qlogis(gamma) + if (kx > 0) sum(xrow * covs_gamma) else 0)
            new_states[part] <- sample(c(1, 2), 1, prob = c(gamma_x, 1 - gamma_x))
          }
        }
        new_obs <- data.frame(
          t = rep(d, hh_n),
          part_id = part_vec,
          pid_global = pid_vec,
          enroll = enroll_vec,
          state = new_states,
          hh_size = rep(hh_n, hh_n),
          hh_id = rep(i, hh_n),
          u_pid = u_vec[pid_vec]
        )
        complete_obs <- bind_rows(complete_obs, new_obs)
        prior <- new_states
      }
    }
    last_x <- last_x + hh_n
  }

  complete_obs <- complete_obs %>% arrange(hh_id, t, part_id)
  outcome <- matrix(nrow = nrow(complete_obs), ncol = length(obs_prob))
  colnames(outcome) <- paste0("y", 1:length(obs_prob))
  for (i in 1:nrow(complete_obs)) {
    for (j in 1:length(obs_prob)) {
      p1 <- obs_prob[[j]][complete_obs$state[i]]
      outcome[i, j] <- sample(c(0, 1), 1, prob = c(1 - p1, p1))
    }
  }
  complete_obs <- complete_obs %>% bind_cols(as.data.frame(outcome))
  if (!complete_enroll) {
    obs <- complete_obs %>% filter(enroll == 1)
  } else {
    obs <- complete_obs
  }

  out <- list(obs = obs, complete_obs = complete_obs)
  if (!is.null(base_x)) out$x <- base_x
  return(out)
}


sim_coinfection_ab <- function(a_args = list(),
                               b_args = list()) {
  sim_a <- do.call(sim_sirs, a_args)
  b_defaults <- list(
    base_complete_obs = sim_a$complete_obs,
    base_x = sim_a$x,
    a_complete_obs = sim_a$complete_obs
  )
  sim_b <- do.call(sim_sis_from_existing, utils::modifyList(b_defaults, b_args))

  merged_complete <- sim_a$complete_obs %>%
    select(t, part_id, pid_global, enroll, hh_size, hh_id, state_a = state) %>%
    left_join(
      sim_b$complete_obs %>%
        select(t, part_id, pid_global, state_b = state),
      by = c("t", "part_id", "pid_global")
    )

  merged_obs <- merged_complete
  if (all(c("obs", "complete_obs") %in% names(sim_a)) && all(c("obs", "complete_obs") %in% names(sim_b))) {
    if (!all(sim_a$obs$enroll == 1)) {
      merged_obs <- merged_complete %>% filter(enroll == 1)
    }
  }

  return(list(
    infection_a = sim_a,
    infection_b = sim_b,
    complete_obs = merged_complete,
    obs = merged_obs
  ))
}

complete_obs <- a$complete_obs
x <- a$x

plot_sim <- function(complete_obs,
                     x,
                     state = 2) { 
  
  indi_covs <- colnames(a$x)[grepl("^x", colnames(a$x))]
  other_covs <- setdiff(colnames(a$x), c(indi_covs, "pid_global"))

  ## Make co-variate plots 
  if(length(indi_covs) + length(other_covs) > 0) {
    complete_obs <- complete_obs %>%
      left_join(x, by = "pid_global")
    
    x %>%
      select(indi_covs) %>%
      pivot_longer(cols = indi_covs) %>%
      group_by(name, value) %>%
      summarise(n = n()) %>%
      ggplot(aes(x = name, y = n, fill = factor(value))) + 
      geom_bar(stat = "identity") + 
      labs(x = "Covariate",
           fill = "Present",
           y = "Count")
  }
  covs_plot <- 
    
  
  
}
