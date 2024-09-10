#' Calculate Buckland CV
#'
#' Calculates 95% confidence intervals from survey CV
#'
#' @param cv Coefficient of variation from survey
#'
#' @returns A scalar used to calculate a 95% confidence interval for survey data.
#'
#' @references \insertRef{buckland1992}{pwsHerringBasa}

calc.buck.cv <- function(cv) {
  return(exp(1.96*sqrt(log(1+(cv^2)))))
}
