#' Posterior predictive distribution (lognormal)
#'
#' Calculates posterior predictive distributions for log-normally distributed survey data.
#' A reparameterization explained
#' \href{https://msalganik.wordpress.com/2017/01/21/making-sense-of-the-rlnorm-function-in-r/}{here}
#' is used.
#'
#' @param model.pred MCMC results for model-predicted values of annual survey data
#' @param sd Model-estimated standard deviation of survey data for each markov chain
#'
#' @returns Matrix of posterior predictive distributions (rows) for annual survey data (columns)

post.pred.lognorm <- function(model.pred, sd){
  PP <- matrix(NA,nrow=nrow(model.pred),ncol=ncol(model.pred))
  for(i in 1:nrow(model.pred)){
    fits <- as.numeric(model.pred[i,]) # Take the true value of the survey estimate
    errors <- as.numeric(sd[i]*fits) # Take the error based on the true value of the survey estimate

    # Need to reparameterize to calculate the log-normal mean and sd for the survey data
    # From: https://msalganik.wordpress.com/2017/01/21/making-sense-of-the-rlnorm-function-in-r/
    E <- log(fits^2/sqrt(errors^2+fits^2))
    SD <- sqrt(log(1+(errors^2/fits^2)))
    PP[i,] <- rlnorm(n=length(fits),meanlog=E,sdlog=SD)
    # MDM.PP[i,] <- rnorm(n=length(fits),mean=E,sd=SD) # Generate a random variable based on the expected value (predictive) using the assumed
  }
  return(PP)
}
