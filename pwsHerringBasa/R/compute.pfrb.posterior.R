#' Compute Pre-fishery Biomass Posterior
#'
#' Estimate pre-fishery biomass (metric tons) posterior distribution from MCMC runs
#'
#' @param dir_mcmc_out String giving directory containing mcmc outputs
#' @param nyr Total number of years in model time series
#'
#' @returns A list of length 3. The first (`biomass.df`) is a single-column
#' `data.frame` giving posterior distribution. The second (`biomass.quants`)
#' is a numeric vector giving the .025, .5, and .975 quantiles of posterior
#' distribution. The third (`prob.below.threshold`) is a scalar giving the probability
#' that pre-fishery biomass lies below 20,000 metric tons (an approximation of the
#' 22,000 short-ton harvest threshold).
#'

compute.pfrb.posterior <- function(dir_mcmc_out, nyr){

  biomass <- read.table(here::here(dir_mcmc_out, "PFRBiomass.csv"), header = FALSE, sep = ",", dec=".")[,nyr]

  biomass.df <- data.frame(biomass=biomass)
  prob.below.threshold <- round(sum(biomass < 20000)/length(biomass), 2)

  biomass.quants <- round(as.vector(apply(as.matrix(biomass), 2, quantile, c(0.025, 0.5, 0.975)))/1000, 2)

  return(listN(biomass.df, biomass.quants, prob.below.threshold))
}
