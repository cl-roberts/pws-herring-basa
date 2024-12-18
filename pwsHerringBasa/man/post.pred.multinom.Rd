% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/post.pred.multinom.R
\name{post.pred.multinom}
\alias{post.pred.multinom}
\title{Posterior predictive distribution (multinomial)}
\usage{
post.pred.multinom(model.pred, ess, nage)
}
\arguments{
\item{model.pred}{MCMC results for model-predicted values of annual survey data}

\item{ess}{vector of effective sample sizes}

\item{nage}{scalar of the number of age classes}
}
\value{
Matrix of posterior predictive distributions (rows) for annual survey data (columns)
}
\description{
Calculates posterior predictive distributions for multinomially distributed survey data.
}
