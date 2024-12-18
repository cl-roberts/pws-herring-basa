% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate.colors.R
\name{generate.colors}
\alias{generate.colors}
\title{Generate Color Vector for Age Composition Plot}
\usage{
generate.colors(nyr, nages = 7, color.options, plus.group = "grey")
}
\arguments{
\item{nyr}{Number of years for which age compositions are modeled}

\item{nages}{Number of age classes available to seine or survey gear (i.e.
7 age classes for herring that recruit at age 3 and using a plus-group of
9+)}

\item{color.options}{Comprehensive list of colors that will be used (and
recycled as necessary) to illustrate cohorts}

\item{plus.group}{Color used for plus group}
}
\value{
Vector of colors recycled from \code{color.options} and \code{plus.group} with
a length equal \code{nyr}*\code{nages}
}
\description{
This utility function generates a vector for coloring age composition plot to
enable ease of tracking cohorts through time
}
