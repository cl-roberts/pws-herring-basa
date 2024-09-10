#' Posterior predictive distribution (negative binomial)
#'
#' Calculates posterior predictive distributions for negative-binomially distributed survey data.
#'
#' @param model.pred MCMC results for model-predicted values of annual survey data
#' @param sd Model-estimated standard deviation of survey data for each markov chain
#'
#' @returns Matrix of posterior predictive distributions (rows) for annual survey data (columns)

post.pred.negbin <- function(model.pred, sd){
  PP <- matrix(NA,nrow=nrow(model.pred),ncol=ncol(model.pred))
  # set.seed(100)
  for(i in 1:nrow(model.pred)){
    fits <- as.numeric(model.pred[i,]) # Take the true value of the survey estimate
    dispersion <- as.numeric(sd[i]) # Take the error based on the true value of the survey estimate

    E <- fits
    SD <- exp(dispersion)
    PP[i,] <- rnbinom(n=length(fits),size=SD,prob=SD/(SD+E))
    # MDM.PP[i,] <- rnorm(n=length(fits),mean=E,sd=SD) # Generate a random variable based on the expected value (predictive) using the assumed
  }
  return(PP)
}
