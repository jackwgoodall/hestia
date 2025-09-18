
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

