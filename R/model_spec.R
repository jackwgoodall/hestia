
library(tidyverse)

#' Defines a state transition in the infection process model which does not
#' represent an transmission (infection) event
#' 
#' @param from string giving name of origin compartment
#' @param to string giving the name of the destination compartment
#'
transit <- function(from, to, ...) {
  
  .dots <- unlist(list(...))
  
  out <- data.frame(from = from,
                    to = to,
                    source = NA,
                    rate_name = NA,
                    rate_value = NA)
  
  out$rate_name <- names(.dots)
  out$rate_value <- .dots
  
  return(out)
  
}

#' Defines a state transition in the infection process model which is
#' the result of a transmission (infection) event
#' 
#' @param from string giving name of origin compartment
#' @param to string giving the name of the destination compartment
#' @param source string (or vector of strings) designating which compartments are infectious. If NULL, the destination compartment is presumed to be the infectious compartment.
#'
transmit <- function(from, to, source = NA) {
  
  out <- data.frame(from = from,
                    to = to,
                    rate_name = NA,
                    rate_value = NA)
  out$source <- ifelse(sum(is.na(source))>0, to, list(source))
  
  return(out)
  
}


#' Builds infection process model.
#' 
#' @param ... a series of transit or transmit function calls
#' @param ih_cov indicator for whether intra-household infection risk is a function of covariates
#' @param eh_cov indicator for whether extra-household infection risk is a function of covariates
#' 
make_infection_model <- function(..., ih_cov = FALSE, eh_cov = FALSE) {
  .dots <- list(...)
   
  out  <- dplyr::bind_rows(.dots)
  
  return(out)
  
}

#' Creates transmission probability matrix, 
#' 
#' @param inf_model infection process model object yielded by make_infection_model()
#'
get_tranmission_details <- function(inf_model) {
  states <- unique(c(inf_model$from, inf_model$to))
  trans <- matrix(1e-10, nrow = length(states), ncol = length(states))
  rownames(trans) <- paste("to", states, sep = "_")
  colnames(trans) <- paste("from", states, sep = "_")
  
  trans_to_fit <- data.frame(from = character(),
                             to = character(),
                             trans_row = numeric(),
                             trans_col = numeric(),
                             source = list(),
                             rate_name = character(),
                             param = numeric())
  
  for(i in 1:nrow(inf_model)) {
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
  }

  if(sum(!is.na(trans_to_fit$rate_name)) > 0) {
    fac_levels <- unique(trans_to_fit$rate_name[!is.na(trans_to_fit$rate_name)])
    trans_to_fit$param <- as.numeric(factor(trans_to_fit$rate_name, levels = fac_levels))
  }
  trans_to_fit <- trans_to_fit %>%
    mutate(param = replace_na(param, 0))
  
  return(list(states = states,
              trans_matrix = trans,
              trans_to_fit = trans_to_fit,
              inf_states = unique(unlist(trans_to_fit$source[is.na(trans_to_fit$rate_name)]))))
}


make_observation_model <- function(...) {
  .dots <- list(...)
  
  ops <- list()
  for(i in 1:length(.dots)) {
    op <- matrix(nrow = 2, ncol = length(.dots[[i]]))
    op[1,] <- 1-.dots[[i]]
    op[2,] <- .dots[[i]]
    rownames(op) <- c("neg_obs", "pos_obs")
    colnames(op) <- names(.dots[[i]])
    ops[[i]] <- op
  }
  names(ops) <- names(.dots)
  return(ops)
}

make_stan_data <- function(inf_model, obs_model, data, init_probs, epsilon = 1e-10) {
  
  inf_details <- get_tranmission_details(inf_model)
  dat <- data %>%
    arrange(hh_id, t, part_id)
  
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
  
  obs_array <- array(dim = c(length(obs_model),2, length(inf_details$states)))
  for(i in 1:length(obs_model)) {
    obs_array[i,,] <- obs_process[[i]]
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
                   multiplier = array(rep(1, nrow(inf_details$trans_to_fit))), #TODO: update user input to support multipliers that aren't 1
                   n_params = length(unique(inf_details$trans_to_fit$param[inf_details$trans_to_fit$param != 0])),
                   n_hh = max(dat$hh_id),
                   hh_size = hh_sum$hh_size,
                   n_obs = nrow(dat),
                   n_obs_type = length(obs_process), 
                   n_unique_obs = 2, #TODO: allow multi-level outcomes
                   y = dat %>% select(names(obs_model)) + 1,
                   part_id = dat$part_id,
                   t_day = dat$t,
                   obs_per_hh = hh_sum$obs_per_hh,
                   hh_start_ind = hh_sum$hh_start_ind,
                   hh_end_ind = hh_sum$hh_end_ind,
                   hh_tmin = hh_sum$hh_tmin,
                   hh_tmax = hh_sum$hh_tmax,
                   obs_prob = obs_array,
                   init_probs = init_probs, #TODO: Toggle to fit
                   epsilon = epsilon)
  
  return(dat_stan)
  
}

# TODO: option to save state probabilities
# TODO: create processed model results
run_model <- function(inf_model, obs_model, data, init_probs, epsilon = 1e-10,
                      file = "stan/hmm.stan",
                      iter = 2000, chains = 4, cores = getOption("mc.cores", 1L),
                      save_chains = FALSE) {
  dat_stan <- make_stan_data(inf_model, obs_model, data, init_probs, epsilon)
  
  stan_fit <- stan(file = file,
                   data = dat_stan,
                   iter = iter,
                   chains = chains,
                   cores = cores)
  
  if(save_chains) {
    ch <- rstan::extract(stan_fit)
  } else {
    ch <- NULL
  }
  
}





