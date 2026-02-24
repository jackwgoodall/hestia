
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

# SIS model 
sim_sis <- function(eh_prob = 0.05, 
                    ih_prob = 0.05, 
                    n_hh = 100,
                    hh_size = 1:5, 
                    tmax = 100, 
                    gamma = 1/3,
                    covs_eh = c(0, 0), 
                    covs_ih = c(0, 0),
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
                    season_coefs_ih = NULL) {         # length 2*season_k (sin1, cos1, sin2, cos2, ...)
  
  epsilon <- 1e-10

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
      data = list(day = day_seq),
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
  
  x <- matrix(nrow = sum(hh_size), ncol = length(covs_ih))
  
  u_pid <- rnorm(sum(hh_size), mean = 0, sd = sigma_pid)
  
  for(i in 1:length(covs_ih)) {
    x[,i] <- rbinom(sum(hh_size), 1, 0.4) # 4/10 get 1, 0
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
  
  enroll_ids <- list()
  
  # Infection state for all household members on all time steps
  complete_obs <- data.frame(t = numeric(),
                             part_id = numeric(),
                             enroll = numeric(),
                             state = numeric(),
                             hh_size = numeric(),
                             hh_id = numeric(),
                             u_pid = numeric())
  
  last_x <- 0
  
  for(i in 1:length(hh_size)) {
    
    # Move HH members through SIS states
    for(d in 1:tmax) {
      wk <- base::ceiling(d/7)
      if(d == 1) {
        
        pid_global_vec <- last_x + part_ids[[i]]
        
        new_obs <- bind_rows(data.frame(t = rep(d, hh_size[i]),
                                        part_id = part_ids[[i]],
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
          if (prior[part] == 1) {  # S state

            eh_prob_x <- plogis(qlogis(eh_prob) +
                                  sum(x[last_x + part, ] * covs_eh) +
                                  u_pid[last_x + part] +                            # individual level variability
                                  season_effect(d, season_coefs_eh, season_period)) # seasonality variability 

            ih_prob_x <- plogis(qlogis(ih_prob) +
                                  sum(x[last_x + part, ] * covs_ih) +
                                  u_pid[last_x + part] +                             # individual level variability
                                  season_effect(d, season_coefs_ih, season_period))  # seasonality variability
            
            no_inf_prob <- (1 - eh_prob_x) * (1 - ih_prob_x)^prior_inf
            
            new_states[part] <- sample(x = c(1, 2), size = 1,
                                       prob = c(no_inf_prob, 1 - no_inf_prob))
          } else {                 # I state
            new_states[part] <- sample(x = c(1, 2), size = 1,
                                       prob = c(gamma, 1 - gamma))
          }
        }
        
        new_obs <- data.frame(
          t = rep(d, hh_size[i]),
          part_id = part_ids[[i]],
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
