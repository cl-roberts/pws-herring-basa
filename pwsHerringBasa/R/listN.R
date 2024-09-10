#' Create named list
#'
#' A wrapper function for \link[base]{list} which names list elements using their
#' object name
#'
#' @param ... Objects to include in named list
#'
#' @returns A named list


listN <- function(...){
  anonList <- list(...)
  names(anonList) <- as.character(substitute(list(...)))[-1]
  anonList
}
