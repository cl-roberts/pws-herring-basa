% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute.catch.biomass.R
\name{compute.catch.biomass}
\alias{compute.catch.biomass}
\title{Compute Catch Biomass}
\usage{
compute.catch.biomass(data, model.dir)
}
\arguments{
\item{data}{List of model data read-in from \code{PWS_ASA.dat}}

\item{model.dir}{String giving model directory (if \code{data} not supplied, this
will be used to read in data)}
}
\value{
A numeric vector giving catch biomass time series calculated from
fishery and ASL survey data.
}
\description{
Estimate total catch biomass across different fleets using fishery dependent
data and ASL survey over full model time series.
}
