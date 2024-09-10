#' Posterior predictive distribution (multinomial)
#'
#' Calculates posterior predictive distributions for multinomially distributed survey data.
#'
#' @param model.pred MCMC results for model-predicted values of annual survey data
#' @param ess vector of effective sample sizes
#' @param nage scalar of the number of age classes
#'
#' @returns Matrix of posterior predictive distributions (rows) for annual survey data (columns)

post.pred.multinom <- function(model.pred, ess, nage){
  pp <- matrix(0, nrow(model.pred), ncol(model.pred)) # CHANGE
  for(i in 1:nrow(model.pred)){ # Loop through the MCMC draws
    for(j in 1:length(ess)){ # Loop through each year
      if(all(model.pred[i, (j*nage-(nage-1)):(j*nage)] == 0)){
        pp[i, (j*nage-(nage-1)):(j*nage)] <- rep(0, times=nage)
      }else{
        pp[i, (j*nage-(nage-1)):(j*nage) ] <- t(rmultinom(1, ess[j],
            model.pred[i, (j*nage-(nage-1)):(j*nage)])) / ess[j] *100
      }
    }
  }
  return(pp)
}
