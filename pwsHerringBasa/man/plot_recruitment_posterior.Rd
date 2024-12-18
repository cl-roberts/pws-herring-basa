% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_recruitment_posterior.R
\name{plot_recruitment_posterior}
\alias{plot_recruitment_posterior}
\title{Plot Age-3 Recruitments}
\usage{
plot_recruitment_posterior(df, years, legend = TRUE)
}
\arguments{
\item{df}{A data frame (or \link[tibble]{tibble}) returned by
\link[pwsHerringBasa]{compute.recruitment} containing model-estimated
age-3 recruits with credible intervals}

\item{years}{A integer vector giving the years in model time series}

\item{legend}{Logical; displays plot legend if TRUE}
}
\value{
A \link[ggplot2]{ggplot} object showing showing model-estimated
time series of age-3 recruits with credible intervals
}
\description{
Wrapper function for plotting the model-estimated age-3 recruits (in millions
of fish) for the full model time series.
}
