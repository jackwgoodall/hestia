
#' @title Logit Transformation
#' @param x numeric value between 0 and 1
#' @export
logit <- function(x) qlogis(x)

#' @title Inverse Logit Transformation
#' @param x numeric value
#' @export
inv_logit <- function(x) plogis(x)


#' @title Checking split parameter specification
#' @description Utility function for checking whether split parameter is properly defined
#'
#' @param to character vector giving the name(s) of the destination compartment
#' @param split numeric or character vector indicating how people moving out of the starting compartment are split between the destination compartments
#'
split_check <- function(to, split) {
  # Make sure splits are valid
  if(!(length(split == 1) & sum(is.na(split)) == 1)) { # single NA OK
    if(sum(is.na(split)) >= 1) { # otherwise NA not OK
      stop("Improper specification for split - cannot contain NA")
    } else if(is.numeric(split)) { # numeric
      if(length(to) == 1) { # single destination compartment
        if(length(split != 1)) {
          stop("Must have length(split) == 1 if length(to) == 1.")
        }
      } else { # multiple destination compartments
        if(length(split) == length(to)) {
          if(sum(split) != 1) {
            stop("Values for split must sum to 1 if providing for all destination compartments.")
          }
        } else { # length(split) != length(to)
          if(length(split) == length(to) - 1) {
            if(sum(split) > 1) {
              stop("split cannot sum to greater than 1")
            }
          } else {
            stop("length(split) must be equal to length(to) or length(to)-1")
          }
        }
      }
    } else { # not numeric
      if(length(to) == 1) {
        if(length(split) != 1) {
          stop("Must have length(split) == 1 if length(to) == 1.")
        }
      } else { # multiple destination compartments
        if(length(split) != length(to)-1) {
          stop("Number of parameter names for split must be one less than the number of destination compartments.")
        }
      }
    }
  }
}


#' @title Define a non-infection transition
#'
#' @description
#' Defines a state transition in the infection process model which does not
#' represent an transmission (infection) event
#'
#' @param from string giving name of origin compartment
#' @param to string or vector of strings giving the name of the destination compartment
#' @param split an optional character or numeric vector indicating what proportion of individuals
#' transition into each of the `to` compartments
#'
#' @export
progress <- function(from, to, split = NA, ...) {

  .dots <- unlist(list(...))

  split_check(to, split)

  out <- list()
  for(i in 1:length(to)) {
    out[[i]] <- data.frame(from = from,
                           to = to[i],
                           source = NA,
                           rate_name = NA,
                           rate_value = NA,
                           split_name = NA,
                           split_value = NA)

    out[[i]]$rate_name <- names(.dots)
    out[[i]]$rate_value <- .dots

    if(i == 1) {
      if(is.numeric(split[i])) {
        out[[i]]$split_value <- split[i]
      } else {
        out[[i]]$split_name <- split[i]
      }
    } else {
      if(is.numeric(split)) {
        if(length(split) >= i) {
          out[[i]]$split_value <- split[i]
        } else {
          out[[i]]$split_value <- 1-sum(split[1:(i-1)])
        }
      } else {
        if(length(split) >= i) {
          out[[i]]$split_name <- split[i]
        } else {
          out[[i]]$split_name <- paste0("1-", paste(split[1:(i-1)], sep = "-", collapse = "-"))
        }
      }
    }

  }

  return(bind_rows(out))

}

#' @title Define a infection transition
#'
#' @description
#' Defines a state transition in the infection process model which is
#' the result of a transmission (infection) event
#'
#' @param from string giving name of origin compartment
#' @param to string giving the name of the destination compartment
#' @param source string (or vector of strings) designating which compartments are infectious.
#' If NULL, the destination compartment is presumed to be the infectious compartment.
#' @param split an optional character or numeric vector indicating what proportion of individuals
#' transition into each of the `to` compartments
#'
#' @export
transmit <- function(from, to, source = NA, split = NA) {

  split_check(to, split)

  out <- list()

  for(i in 1:length(to)) {

    out[[i]] <- data.frame(from = from,
                           to = to[i],
                           rate_name = NA,
                           rate_value = NA,
                           split_name = NA,
                           split_value = NA,
                           mult_inf_prob = NA)

    out[[i]]$source <- ifelse(sum(is.na(source))>0, to, list(source))

    if(i == 1) {
      if(is.numeric(split[i])) {
        out[[i]]$split_value <- split[i]
      } else {
        out[[i]]$split_name <- split[i]
      }
    } else {
      if(is.numeric(split)) {
        if(length(split) >= i) {
          out[[i]]$split_value <- split[i]
        } else {
          out[[i]]$split_value <- 1-sum(split[1:(i-1)])
        }
      } else {
        if(length(split) >= i) {
          out[[i]]$split_name <- split[i]
        } else {
          out[[i]]$split_name <- paste0("1-", paste(split[1:(i-1)], sep = "-", collapse = "-"))
        }
      }
    }
  }

  return(bind_rows(out))

}


#' @title Build Infection Process Model
#'
#' @param ... a series of \link{progress} or \link{transmit} function calls
#' @param mult_inf_probs If FALSE then all infection probabilities are shared
#' across infectious compartments
#'
#' @export
make_infection_model <- function(..., mult_inf_probs = FALSE) {
  .dots <- list(...)

  out  <- dplyr::bind_rows(.dots)
  out$mult_inf_probs <- mult_inf_probs

  return(out)

}

#' @title Create Transmission Probability Matrix
#'
#' @param inf_model infection process model object yielded by make_infection_model()
#'
#' @export
get_transmission_details <- function(inf_model) {
  states <- unique(c(inf_model$from, inf_model$to))
  trans <- matrix(1e-10, nrow = length(states), ncol = length(states))
  rownames(trans) <- paste("to", states, sep = "_")
  colnames(trans) <- paste("from", states, sep = "_")
  mult <- matrix(1, nrow = length(states), ncol = length(states))
  rownames(mult) <- paste("to", states, sep = "_")
  colnames(mult) <- paste("from", states, sep = "_")

  trans_to_fit <- data.frame(from = character(),
                             to = character(),
                             trans_row = numeric(),
                             trans_col = numeric(),
                             source = list(),
                             rate_name = character(),
                             param = numeric())

  mult_to_fit <- data.frame(from = character(),
                            to = character(),
                            mult_row = numeric(),
                            mult_col = numeric(),
                            mult_name = character(),
                            param = numeric())

  for(i in 1:nrow(inf_model)) {

    # Transition probabilities
    if(!is.na(inf_model$rate_value[i])) {
      trans[states == inf_model$to[i],states == inf_model$from[i]] <- inf_model$rate_value[i]
    } else {
      temp <- data.frame(from = inf_model$from[i],
                         to = inf_model$to[i],
                         trans_row = which(states == inf_model$to[i]),
                         trans_col = which(states == inf_model$from[i]),
                         rate_name = inf_model$rate_name[i],
                         param = NA)
      temp$source <- ifelse(is.null(inf_model$source[i][[1]]), list(0) , list(which(states %in% inf_model$source[i][[1]])))
      trans_to_fit <- bind_rows(trans_to_fit, temp)
    }

    # Multipliers
    if(!is.na(inf_model$split_value[i])) {
      mult[states == inf_model$to[i],states == inf_model$from[i]] <- inf_model$split_value[i]
    } else if(!is.na(inf_model$split_name[i])) {
      temp <- data.frame(from = inf_model$from[i],
                         to = inf_model$to[i],
                         mult_row = which(states == inf_model$to[i]),
                         mult_col = which(states == inf_model$from[i]),
                         mult_name = inf_model$split_name[i],
                         param = NA)
      mult_to_fit <- bind_rows(mult_to_fit, temp)
    }
  }

  # Identify unique parameters to fit - transitions
  if(sum(!is.na(trans_to_fit$rate_name)) > 0) {
    fac_levels <- unique(trans_to_fit$rate_name[!is.na(trans_to_fit$rate_name)])
    trans_to_fit$param <- as.numeric(factor(trans_to_fit$rate_name, levels = fac_levels))
  }
  trans_to_fit <- trans_to_fit %>%
    mutate(param = replace_na(param, 0))

  # Identify unique parameters to fit - multipliers
  if(sum(!is.na(mult_to_fit$mult_name)) > 0) {
    fac_levels <- unique(mult_to_fit$mult_name[!is.na(mult_to_fit$mult_name) & !str_detect(mult_to_fit$mult_name, "1-")])
    mult_to_fit$param <- as.numeric(factor(mult_to_fit$mult_name, levels = fac_levels))
  }
  mult_to_fit <- mult_to_fit %>%
    mutate(param = replace_na(param, 0))

  if(nrow(mult_to_fit) > 0) {
    for(i in 1:nrow(mult_to_fit)) {
      if(grepl("1-", mult_to_fit$mult_name[i])) {
        to_match <- strsplit(mult_to_fit$mult_name[i], "-")[[1]]
        to_match <- to_match[to_match != "1"]
        matches <- mult_to_fit$param[mult_to_fit$mult_name %in% to_match]
        mult_to_fit$param[i] <- list(-1*matches)
      }
    }
  }

  return(list(states = states,
              trans_matrix = trans,
              mult_matrix = mult,
              trans_to_fit = trans_to_fit,
              mult_to_fit = mult_to_fit,
              inf_states = unique(unlist(trans_to_fit$source[is.na(trans_to_fit$rate_name)])),
              mult_inf_probs = inf_model$mult_inf_probs[1]))
}

#' @title Make Observation Process Model
#'
#' @description
#' Specifies Beta distribution priors on the probability of a positive
#' observation for each test type and compartment combination. Used with
#' \code{hmm_tv_cov_reduce_sum_obs_prior.stan}, which estimates these
#' probabilities rather than treating them as fixed.
#'
#' @param ... A series of named lists. Each list corresponds to an observation
#'   type (e.g. \code{pcr}, \code{igg}). Each element of the list is named for
#'   a compartment and holds a length-2 numeric vector \code{c(alpha, beta)}
#'   giving the Beta prior hyperparameters for
#'   P(positive observation | that compartment).
#'   The prior mean is \code{alpha / (alpha + beta)}.
#'
#' @return A named list, one element per observation type. Each element is
#'   itself a list with named numeric vectors \code{alpha} and \code{beta}
#'   (compartment names preserved).
#'
#' @examples
#' make_observation_model(
#'   pcr = list("S" = c(1, 19), "I" = c(19, 1), "R" = c(1, 19)),
#'   igg = list("S" = c(1, 99), "I" = c(1, 99), "R" = c(16, 4))
#' )
#'
#' @export
make_observation_model <- function(...) {
  .dots <- list(...)

  ops <- list()
  for (i in seq_along(.dots)) {
    test_priors <- .dots[[i]]
    if (!is.list(test_priors)) {
      stop(
        "Each observation type must be a named list of c(alpha, beta) vectors. ",
        "Got a non-list for '", names(.dots)[i], "'. ",
        "Example: list(S = c(1, 19), I = c(19, 1))"
      )
    }
    alpha_vec <- sapply(test_priors, `[`, 1)
    beta_vec  <- sapply(test_priors, `[`, 2)
    names(alpha_vec) <- names(test_priors)
    names(beta_vec)  <- names(test_priors)

    # Derive per-state bounds from the prior mean to break HMM
    # label-switching symmetry.  Mean >= 0.5 (positive state) → [0.5, 1];
    # otherwise (negative state) → [0, 0.5].
    mean_vec <- alpha_vec / (alpha_vec + beta_vec)
    lb_vec   <- ifelse(mean_vec >= 0.5, 0.5, 0)
    ub_vec   <- ifelse(mean_vec >= 0.5, 1,   0.5)
    names(lb_vec) <- names(test_priors)
    names(ub_vec) <- names(test_priors)

    ops[[i]] <- list(alpha = alpha_vec, beta = beta_vec,
                     lb = lb_vec, ub = ub_vec)
  }
  names(ops) <- names(.dots)
  return(ops)
}


#' @title Create Data for Stan Model
#'
#' @param inf_model Infection process model object generated by \link{make_infection_model}
#' @param obs_model Observation process model object generated by \link{make_observation_model}
#' @param data data frame with columns....
#' @param init_probs vector of initial probabilities for each infection process model state
#' @param epsilon very small number to use for zero probability transitions
#' @param ih_cov NULL for run without covariates, otherwise data frame with intra-household covariates for each participant
#' @param eh_cov NULL for run without covariates, otherwise data frame with extra-household covariates for each participant
#' 
make_stan_data <- function(inf_model, obs_model, data, init_probs, epsilon = 1e-10,
                           ih_cov = NULL, eh_cov = NULL) {


  inf_details <- get_transmission_details(inf_model)
  dat <- data %>%
    arrange(hh_id, t, part_id)

  dat$row_id <- 1:nrow(dat)
  hh_sum <- dat %>%
    group_by(hh_id, hh_size) %>%
    summarize(hh_start_ind = min(row_id),
              hh_end_ind = max(row_id),
              hh_tmin = min(t),
              hh_tmax = max(t),
              obs_per_hh = n()) %>%
    ungroup()

  source_state_matrix <- matrix(0, nrow = nrow(inf_details$trans_to_fit), ncol = length(inf_details$states))
  for(i in 1:nrow(inf_details$trans_to_fit)) {
    if(any(inf_details$trans_to_fit$source[[i]] == 0)) {
      next
    } else {
      source_state_matrix[i,inf_details$trans_to_fit$source[[i]]] <- 1
    }
  }

  obs_prob_alpha <- matrix(nrow = length(obs_model), ncol = length(inf_details$states))
  obs_prob_beta  <- matrix(nrow = length(obs_model), ncol = length(inf_details$states))
  obs_lb         <- matrix(nrow = length(obs_model), ncol = length(inf_details$states))
  obs_ub         <- matrix(nrow = length(obs_model), ncol = length(inf_details$states))
  for (i in seq_along(obs_model)) {
    sn <- inf_details$states
    obs_prob_alpha[i, ] <- obs_model[[i]]$alpha[sn]
    obs_prob_beta[i, ]  <- obs_model[[i]]$beta[sn]
    obs_lb[i, ]         <- obs_model[[i]]$lb[sn]
    obs_ub[i, ]         <- obs_model[[i]]$ub[sn]
  }

  # Expand multipliers if needed
  if(is.list(inf_details$mult_to_fit$param)) {
    mult_info <- data.frame()
    for(i in 1:nrow(inf_details$mult_to_fit)) {
      if(length(inf_details$mult_to_fit$param[i][[1]]) == 1) {
        mult_info <- bind_rows(mult_info, inf_details$mult_to_fit[i,])
      } else {
        for(j in 1:length(inf_details$mult_to_fit$param[i][[1]])) {
          temp <- inf_details$mult_to_fit[i,] %>%
            unnest(param)
          mult_info <- bind_rows(mult_info, temp)
        }
      }
    }
  } else {
    mult_info <- inf_details$mult_to_fit
  }

  # TODO: deal with missing observations (change NA to -1)

  dat_stan <- list(n_states = length(inf_details$states),
                   trans = inf_details$trans_matrix,
                   n_inf_states = length(inf_details$inf_states),
                   inf_states = array(inf_details$inf_states),
                   n_trans_fit = nrow(inf_details$trans_to_fit),
                   param_index = array(inf_details$trans_to_fit$param),
                   trans_index = inf_details$trans_to_fit %>% select(trans_row, trans_col),
                   source_states = source_state_matrix,
                   transition_multiplier = inf_details$mult_matrix,
                   n_mult_fit = nrow(mult_info),
                   n_mult_params = length(unique(abs(unlist(mult_info$param)))),
                   mult_param_index = unlist(mult_info$param),
                   mult_index = mult_info %>% select(mult_row, mult_col),
                   n_params = length(unique(inf_details$trans_to_fit$param[inf_details$trans_to_fit$param != 0])),
                   n_hh = max(dat$hh_id),
                   hh_size = hh_sum$hh_size,
                   n_obs = nrow(dat),
                   n_obs_type = length(obs_model),
                   n_unique_obs = 2, #TODO: allow multi-level outcomes
                   y = dat %>% select(names(obs_model)) + 1,
                   part_id = dat$part_id,
                   t_day = dat$t,
                   obs_per_hh = hh_sum$obs_per_hh,
                   hh_start_ind = hh_sum$hh_start_ind,
                   hh_end_ind = hh_sum$hh_end_ind,
                   hh_tmin = hh_sum$hh_tmin,
                   hh_tmax = hh_sum$hh_tmax,
                   obs_prob_alpha = obs_prob_alpha,
                   obs_prob_beta  = obs_prob_beta,
                   obs_lb         = obs_lb,
                   obs_ub         = obs_ub,
                   init_probs = init_probs, #TODO: Toggle to fit
                   epsilon = epsilon,
                   n_inf_prob = ifelse(inf_details$mult_inf_probs, length(inf_details$inf_states), 1))

  # Get covariate information if applicable
  if(!(is.null(eh_cov) & is.null(ih_cov))) {
    if(!is.null(eh_cov) & !is.null(ih_cov)) {

      k_ih <- ncol(ih_cov)
      x_ih <- ih_cov

      k_eh <- ncol(eh_cov)
      x_eh <- eh_cov

      dat_stan <- append(dat_stan,
                         list(k_ih = k_ih,
                              x_ih = x_ih,
                              k_eh = k_eh,
                              x_eh = x_eh))

    } else {
      stop("Currently only support having covariates on both intra- and extra-household infeciton probabilities.")
    }
  }

  return(dat_stan)

}

#' @title Create Data for Stan Model (Time-Varying Covariates)
#'
#' @description
#' Variant of \link{make_stan_data} for when intra- and extra-household
#' covariates vary by both person and day. Adds \code{T_global} to the data
#' list and passes the covariate arrays in the 3-D format expected by
#' \code{hmm_cov_viral.stan}.
#'
#' @param inf_model Infection process model object generated by \link{make_infection_model}
#' @param obs_model Observation process model object generated by \link{make_observation_model}
#' @param data data frame with columns....
#' @param init_probs vector of initial probabilities for each infection process model state
#' @param epsilon very small number to use for zero probability transitions
#' @param ih_cov 3-D array \code{[T_global, N_people, k_ih]} of intra-household
#'   covariates. \code{ih_cov[t, p, ]} gives person \code{p}'s covariate vector
#'   on day \code{t}.
#' @param eh_cov 3-D array \code{[T_global, N_people, k_eh]} of extra-household
#'   covariates with the same indexing as \code{ih_cov}.
#'
make_stan_data_tv <- function(inf_model, obs_model, data, init_probs, epsilon = 1e-10,
                              ih_cov = NULL, eh_cov = NULL) {

  # Build the base data list (identical to make_stan_data up to the covariate block)
  dat_stan <- make_stan_data(inf_model, obs_model, data, init_probs, epsilon,
                             ih_cov = NULL, eh_cov = NULL)

  if (!(is.null(eh_cov) & is.null(ih_cov))) {
    if (!is.null(eh_cov) & !is.null(ih_cov)) {

      # Validate array dimensions
      if (length(dim(ih_cov)) != 3) stop("ih_cov must be a 3-D array [T_global, N_people, k_ih].")
      if (length(dim(eh_cov)) != 3) stop("eh_cov must be a 3-D array [T_global, N_people, k_eh].")

      T_global <- dim(ih_cov)[1]

      if (dim(eh_cov)[1] != T_global)
        stop("ih_cov and eh_cov must have the same first dimension (T_global).")
      if (T_global < max(dat_stan$hh_tmax))
        stop("T_global (dim(ih_cov)[1]) must be >= max(hh_tmax).")
      if (dim(ih_cov)[2] != sum(dat_stan$hh_size))
        stop("dim(ih_cov)[2] must equal sum(hh_size) (total number of people).")
      if (dim(eh_cov)[2] != sum(dat_stan$hh_size))
        stop("dim(eh_cov)[2] must equal sum(hh_size) (total number of people).")

      k_ih <- dim(ih_cov)[3]
      k_eh <- dim(eh_cov)[3]

      dat_stan <- append(dat_stan,
                         list(T_global = T_global,
                              k_ih     = k_ih,
                              x_ih     = ih_cov,
                              k_eh     = k_eh,
                              x_eh     = eh_cov))

    } else {
      stop("Currently only support having covariates on both intra- and extra-household infection probabilities.")
    }
  }

  return(dat_stan)
}


#' @title Run Household Transmission Model
#' 
#' @param inf_model Infection process model object generated by \link{make_infection_model}
#' @param obs_model Observation process model object generated by \link{make_observation_model}
#' @param data data frame with columns....
#' @param init_probs vector of initial probabilities for each infection process model state
#' @param epsilon very small number to use for zero probability transitions
#' @param ih_cov NULL for run without covariates. For \code{time_varying = FALSE},
#'   a matrix of intra-household covariates (one row per person). For
#'   \code{time_varying = TRUE}, a 3-D array \code{[T_global, N_people, k_ih]}.
#' @param eh_cov NULL for run without covariates. Same shape convention as
#'   \code{ih_cov} but for extra-household covariates.
#' @param time_varying logical. If \code{TRUE}, covariates are treated as
#'   varying by person and day and \code{make_stan_data_tv} is used to build
#'   the Stan data list. Defaults to \code{FALSE}.
#' @param file file path to stan model
#' @param iter number of MCMC iterations (total - split equally into warmup and
#'   sampling for cmdstanr)
#' @param chains number of MCMC chains
#' @param cores number of cores for parallelization
#' @param init initial conditions for MCMC chains
#' @param save_chains indicator for whether to save MCMC chains
#' @param save_states indicator for whether to save state probabilities
#' @param backend character string, either \code{"rstan"} or \code{"cmdstanr"}.
#'   Defaults to \code{"rstan"}.
#'
#' @export
run_model <- function(inf_model, 
                      obs_model, 
                      data, 
                      init_probs, 
                      epsilon    = 1e-10,
                      ih_cov     = NULL, 
                      eh_cov     = NULL, 
                      time_varying = FALSE,
                      file       = "stan/hmm.stan", 
                      iter       = 2000, 
                      chains     = 4,
                      parallel_chains = getOption("mc.cores", 1L), 
                      threads_per_chain = 4,    # 16 cores total for default
                      adapt_delta       = 0.9,
                      max_treedepth     = 12,
                      init       = NULL,
                      save_chains = TRUE, 
                      save_states = TRUE,
                      backend    = c("rstan", "cmdstanr")) {
  
  # Validate backend argument
  backend <- match.arg(backend)
  
  # Check the requested backend is installed
  if (backend == "rstan" && !requireNamespace("rstan", quietly = TRUE)) {
    stop("rstan is not installed. Install it with install.packages('rstan') or use backend = 'cmdstanr'.")
  }
  if (backend == "cmdstanr" && !requireNamespace("cmdstanr", quietly = TRUE)) {
    stop("cmdstanr is not installed. See https://mc-stan.org/cmdstanr/ or use backend = 'rstan'.")
  }
  
  # Build Stan data list
  if (time_varying) {
    dat_stan <- make_stan_data_tv(inf_model, obs_model, data, init_probs, epsilon, ih_cov, eh_cov)
  } else {
    dat_stan <- make_stan_data(inf_model, obs_model, data, init_probs, epsilon, ih_cov, eh_cov)
  }
  
  # Build init list if not supplied
  if (is.null(init)) {
    if (!is.null(eh_cov) & !is.null(ih_cov)) {
      init <- rep(list(list(
        logit_params       = array(rep(logit(0.5), dat_stan$n_params)),
        logit_mult_params  = array(rep(logit(0.5), dat_stan$n_mult_params)),
        beta_eh            = rep(0, dat_stan$k_eh),
        beta_ih            = rep(0, dat_stan$k_ih),
        beta0_eh           = logit(0.02),
        beta0_ih           = array(rep(logit(0.02), dat_stan$n_inf_prob))
      )), chains)
    } else {
      init <- rep(list(list(
        logit_params       = array(rep(logit(0.5), dat_stan$n_params)),
        logit_mult_params  = array(rep(logit(0.5), dat_stan$n_mult_params)),
        beta_eh            = logit(0.02),
        beta_ih            = array(rep(logit(0.02), dat_stan$n_inf_prob))
      )), chains)
    }
  } else {
    init <- rep(list(init), chains)
  }
  
  # rstan -------
  if (backend == "rstan") {
    
    library(rstan)
    
    if (save_states) {
      stan_fit <- rstan::stan(
        file   = file,
        data   = dat_stan,
        iter   = iter,
        chains = chains,
        cores  = cores,
        init   = init
      )
    } else {
      stan_fit <- rstan::stan(
        file    = file,
        data    = dat_stan,
        iter    = iter,
        chains  = chains,
        cores   = cores,
        init    = init,
        pars    = c("logalpha", "ih_prob", "eh_prob"),
        include = FALSE
      )
    }
  
    # cmdstanr -------  
  } else if (backend == "cmdstanr") {
    
    library(cmdstanr)
    
    # Helper to coerce data frames to matrices in the stan data list
    # cmdstanr is stricter than rstan about data types
    prep_stan_data_cmdstanr <- function(dat_stan) {
      lapply(dat_stan, function(x) {
        if (is.data.frame(x)) as.matrix(x) else x
      })
    }
    
    dat_stan <- prep_stan_data_cmdstanr(dat_stan)
    
    mod <- cmdstanr::cmdstan_model(file,
                                   cpp_options = list(stan_threads = TRUE))
    
    if (save_states) {
      stan_fit <- mod$sample(
        data            = dat_stan,
        iter_warmup     = iter / 2,
        iter_sampling   = iter / 2,
        chains          = chains,
        parallel_chains = cores,
        init            = init
      )
      
    } else {
      
      stan_fit <- mod$sample(
        data            = dat_stan,
        iter_warmup     = iter / 2,
        iter_sampling   = iter / 2,
        chains          = chains,
        parallel_chains = parallel_chains,
        threads_per_chain = threads_per_chain,
        adapt_delta = adapt_delta,
        max_treedepth = max_treedepth,
        init            = init
      )
    }
  }
  
  return(stan_fit)
  
}

