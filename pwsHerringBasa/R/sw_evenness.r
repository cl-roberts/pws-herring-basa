#' Shannon-Weiner Evenness 
#'
#' Calculates the Shannon-Weiner evenness index score for age composition model 
#' estimates.
#'
#' @param naa A model-estimated numbers-at-age vector for a single year
#'
#' @returns A scalar, J, indicating the Shannon-Weiner evenness index score for
#' a single year's estimated age composition
#'
#' @references \insertRef{shannonweaver1949}{pwsHerringBasa}

sw_evenness <- function(naa){

    if(sum(naa) != 1.0){
        naa <- naa / sum(naa)
    }

    as.log <- log(naa)
    as.log[is.infinite(as.log)] <- 0

    H <- -sum(naa*as.log)
    J <- H / log(length(naa))
    
    return(J)

}