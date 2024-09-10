#' Compute Biomass Trajectory
#'
#' Calculate model-estimated biomass credible intervals over full time-series
#'
#' @param dir_mcmc_out String giving directory containing mcmc outputs
#' @param nyr Total number of years in model time series
#'
#' @returns A \link[tibble]{tibble} containing biomass (metric tons) over the full
#' model time series, as well as credible intervals and the probability that the
#' stock is below the harvest threshold (22000 tons or approximately 19958 metric tons).
#'

compute.biomass.traj <- function(dir_mcmc_out, nyr) {

  year <- NULL
  start.year <- 1980
  curr.year <- start.year+nyr

  pfrb <- compute.pfrb.posterior(dir_mcmc_out, nyr+1)

  biomass <- read.biomass.estimates(dir_mcmc_out, nyr) |>
    mutate(forecast = unlist(pfrb$biomass.df))

  biomass.df <- tibble::as_tibble(biomass) |>
    tidyr::pivot_longer(dplyr::everything(), names_to="year", values_to="biomass") |>
    dplyr::group_by(year) |>
    ggdist::median_qi(biomass, .width=c(0.10, 0.50, 0.95)) |>
    print(n=10)

  prob.below.threshhold <- apply(biomass, 2, function(x) sum(x < (22000*0.907185)))/nrow(biomass)

  biomass.df$prob <- as.vector(rep(prob.below.threshhold, 3))

  return(biomass.df)
}
