% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_biomass_trajectory.R
\name{plot_biomass_trajectory}
\alias{plot_biomass_trajectory}
\title{Plot Biomass Trajectory}
\usage{
plot_biomass_trajectory(df, years, legend = TRUE, show.probs = TRUE)
}
\arguments{
\item{df}{A data frame (or \link[tibble]{tibble}) returned by
\link[pwsHerringBasa]{compute.biomass.traj} containing model-estimated biomass
time series with credible intervals}

\item{years}{Integer vector of years in time series}

\item{legend}{Logical; if TRUE display legend}

\item{show.probs}{Logical; if TRUE display biomass credible intervals}
}
\value{
A \link[ggplot2]{ggplot} object showing model-estimated biomass
(metric tons) over the full model time series, as well as credible intervals
}
\description{
Wrapper function for plotting the model-estimated biomass trajectory with
credible intervals.
}
