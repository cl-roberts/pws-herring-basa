#' Set initial TMB parameter values
#'
#' Provides a set of initial parameter values for BASA's MCMC chains. Upper
#' and lower bounds for parameter values are provided as input and random
#' error is added for initializing each markov chain.
#'
#' @param chains The number of markov chains used in the model.
#' @param lower Vector of lower bounds for parameter initialization
#' @param upper Vector of upper bounds for parameter initialization
#'
#' @returns List of named parameter values to initialize each markov chain in BASA.

init_tmb_params <- function(chains, lower, upper) {

  inits <- vector(mode = "list", length = chains)

  for(c in 1:chains) {
    inits[[c]] <- runif(length(lower), min = lower, max = upper)
    names(inits[[c]]) <- names(lower)
  }

  return(inits)
  
}
