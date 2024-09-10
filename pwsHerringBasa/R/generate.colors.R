#' Generate Color Vector for Age Composition Plot
#'
#' This utility function generates a vector for coloring age composition plot to
#' enable ease of tracking cohorts through time
#'
#' @param nyr Number of years for which age compositions are modeled
#' @param nages Number of age classes available to seine or survey gear (i.e.
#' 7 age classes for herring that recruit at age 3 and using a plus-group of
#' 9+)
#' @param color.options Comprehensive list of colors that will be used (and
#' recycled as necessary) to illustrate cohorts
#' @param plus.group Color used for plus group
#'
#' @returns Vector of colors recycled from `color.options` and `plus.group` with
#' a length equal `nyr`*`nages`
#'


generate.colors <- function(nyr, nages = 7, color.options, plus.group = "grey"){

  n <- nyr*nages

  shifter <- function(x, n = 1) {
    return(if (n == 0) x else c(tail(x, -n), head(x, n)))
  }

  colors <- rep(NA, n)
  i=1
  j=0
  cs <- color.options
  while(i < n){
    cs <- shifter(cs, length(color.options)-1)
    for(c in cs){
      colors[i] <- c
      colors[i] <- c
      i = i+1
    }
    colors[i] <- "gray"
    colors[i] <- "gray"
    i = i+1
    j = j+1
  }
  return(colors)
}
