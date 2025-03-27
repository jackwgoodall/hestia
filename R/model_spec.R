
library(tidyverse)

#' Defines a state transition in the infection process model which does not
#' represent an transmission (infection) event
#' 
#' @param from string giving name of origin compartment
#' @param to string giving the name of the destination compartment
#'
transit <- function(from, to, rate, ...) {
  
  .dots <- unlist(list(...))
  
  out <- data.frame(from = from,
                    to = to,
                    source = NA,
                    rate_name = NA,
                    rate_value = NA)
  
  if(is.numeric(.dots) == 1) {
    out$rate_name <- names(.dots)
    out$rate_value <- .dots
  } else {
    out$rate_name <- deparse(substitute(rate))
  }
  
  return(out)
  
}

#' Defines a state transition in the infection process model which is
#' the result of a transmission (infection) event
#' 
#' @param from string giving name of origin compartment
#' @param to string giving the name of the destination compartment
#' @param source string (or vector of strings) designating which compartments are infectious. If NULL, the destination compartment is presumed to be the infectious compartment.
#'
transmit <- function(from, to, source = NULL, ...) {
  
  
  
  out <- data.frame(from = from,
                    to = to,
                    source = ifelse(is.null(source), to, source),
                    rate_name = NA,
                    rate_value = NA)
  
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


get_tranmission_details <- function(inf_model) {
  states <- unique(c(inf_model$from, inf_model$to))
  trans <- matrix(1e-10, nrow = length(states), ncol = length(states))
  rownames(trans) <- paste("to", states, sep = "_")
  colnames(trans) <- paste("from", states, sep = "_")
  
  trans_to_fit <- data.frame(from = character(),
                              to = character(),
                              trans_row = numeric(),
                              trans_col = numeric(),
                              source = numeric(),
                              rate_name = character())
  
  for(i in 1:nrow(inf_model)) {
    if(!is.na(inf_model$rate_value[i])) {
      trans[states == inf_model$to[i],states == inf_model$from[i]] <- inf_model$rate_value[i]
    } else {
      trans_to_fit <- bind_rows(trans_to_fit,
                                 data.frame(from = inf_model$from[i],
                                            to = inf_model$to[i],
                                            trans_row = which(states == inf_model$to[i]),
                                            trans_col = which(states == inf_model$from[i]),
                                            source = ifelse(is.na(inf_model$source[i]), 0 , which(states == inf_model$source[i])),
                                            rate_name = inf_model$rate_name[i]))
    }
  }

  
  
  return(list(trans_matrix = trans,
              trans_to_fit = trans_to_fit,
              inf_states = trans_to_fit$source[is.na(trans_to_fit$rate_name)]))
}

