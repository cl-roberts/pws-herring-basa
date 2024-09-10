#' Read in Biomass Estimates
#'
#' Read in model-estimated biomass over full time-series
#'
#' @param dir_mcmc_out String giving directory containing mcmc outputs
#' @param nyr Total number of years in model time series
#'
#' @returns A data frame containing MCMC results (rows) for annual biomass
#' (columns, in metric tons)
#'

read.biomass.estimates <- function(dir_mcmc_out, nyr = NA) {
  fname <- here::here(dir_mcmc_out, "PFRBiomass.csv")
  biomass.data <- read.table(fname, header = FALSE, sep = ",", dec=".")

  nyr <- ifelse(is.na(nyr), ncol(biomass.data)-1, nyr)
  years <- 1980:(1980+nyr)

  colnames(biomass.data) <- years

  return(biomass.data[,1:nyr])

}
