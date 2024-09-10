#' Compute Catch Biomass
#'
#' Estimate total catch biomass across different fleets using fishery dependent
#' data and ASL survey over full model time series.
#'
#' @param data List of model data read-in from `PWS_ASA.dat`
#' @param model.dir String giving model directory (if `data` not supplied, this
#' will be used to read in data)
#'
#' @returns A numeric vector giving catch biomass time series calculated from
#' fishery and ASL survey data.
#'

compute.catch.biomass <- function(data, model.dir){

  if((!exists("data")) & (exists("model.dir"))){
    data <- read.data.files(model.dir)$PWS_ASA.dat
  } else if((!exists("data")) & (!exists("model.dir"))){
    stop("Must supply data or model.dir!")
  }

  nyr <- data$nyr

  weight.at.age <- data$waa[1:nyr,]

  fb.nya <- data$foodbait_catch[1:nyr,]
  pound.nya <- data$pound_catch[1:nyr,]
  gillnet.nya <- data$gillnet_catch[1:nyr,]
  seine.yield  <- data$seine_yield[1:nyr]

  fb.nya.biomass <- weight.at.age * fb.nya
  pound.nya.biomass <- weight.at.age * pound.nya
  gillnet.nya.biomass <- weight.at.age * gillnet.nya

  fb.biomass.annual       <- rowSums(fb.nya.biomass) # Now sum the biomass over all age classes for each year
  pound.biomass.annual    <- rowSums(pound.nya.biomass)
  gillnet.biomass.annual  <- rowSums(gillnet.nya.biomass)

  fb.biomass.annual       <- replace(fb.biomass.annual, fb.biomass.annual == 0, NA)
  pound.biomass.annual    <- replace(pound.biomass.annual, pound.biomass.annual == 0, NA)
  gillnet.biomass.annual  <- replace(gillnet.biomass.annual, gillnet.biomass.annual == 0, NA)
  seine.yield             <- replace(seine.yield, seine.yield == 0, NA)

  # Matrix of catches by gear type in mt
  total.catch <- cbind(fb.biomass.annual, pound.biomass.annual, gillnet.biomass.annual, seine.yield)
  total.catch[is.na(total.catch)] <- 0
  total.catch.biomass <- rowSums(total.catch) # total catches by year in mt

  return(total.catch.biomass)

}
