#' Compute Recruitment Time Series
#'
#' Calculate model-estimated age-3 recruits (in millions of fish) with credible
#' intervals.
#'
#' @param dir_mcmc_out String giving directory containing mcmc outputs
#' @param nyr Total number of years in model time series
#' @param years Integer vector of years in model
#'
#' @returns A \link[tibble]{tibble} containing model-estimated age-3 recruits
#' (millions of fish) with credible intervals for full time series in model.
#'

compute.recruitment <- function(dir_mcmc_out, nyr, years){
  year <- recruits <- NULL

  age.3.recruits <- read.table(here::here(dir_mcmc_out, "Age3.csv"), header = FALSE, sep = ",", dec=".")[,1:nyr]
  colnames(age.3.recruits) <- years

  age.3.recruits.df <- tibble::as_tibble(age.3.recruits) %>%
    tidyr::pivot_longer(dplyr::everything(), names_to="year", values_to="recruits") %>%
    dplyr::group_by(year) %>%
    ggdist::median_qi(recruits, .width=c(0.50, 0.95)) %>%
    print(n=10)

  return(age.3.recruits.df)
}
