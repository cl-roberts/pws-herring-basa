#' Set initial TMB parameter values
#'
#' Provides a set of initial parameter values for BASA's MCMC chains. Upper
#' and lower bounds for parameter values are provided as input and random
#' error is added for initializing each markov chain.
#'
#' @param start Vector of starting values for parameters
#' @param lower Vector of lower bounds for parameter initialization
#' @param upper Vector of upper bounds for parameter initialization
#' @param chains Integer number of markov chains for which initial values are desired
#' @param seed Integer seed for generating random parameter starting values
#'
#' @returns List of named parameter values to initialize each markov chain in BASA.

init_tmb_params <- function(start, lower, upper, chains, seed) {

    set.seed(seed)

    par_names <- names(lower) 
    sd <- abs(.75*start)

    inits <- list()
    
    for(j in 1:chains){

        inits[[j]] <- rnorm(length(start), start, sd = sd)

        while (any(inits[[j]] < lower) | any(inits[[j]] > upper)) {
            inits[[j]] <- rnorm(length(start), start, sd = sd)
        }
    
        names(inits[[j]]) <- par_names
    }

    return(inits)
}
