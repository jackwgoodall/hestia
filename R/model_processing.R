logit <- function(x) {
  log(x/(1-x))
}

inv_logit <- function(x) {
  exp(x)/(1+exp(x))
}



