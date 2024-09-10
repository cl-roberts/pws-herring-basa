#' Compute retrospective matrix
#'
#' Summarizes retrospective BASA model runs into a matrix for plotting in the
#' `plot_retrospective.r` script.
#'
#' @param n.peels Integer giving the maximum number of years peeled in 
#' retrospective analysis
#' @param dir_mcmc_out String giving directory containing base model mcmc outputs
#'
#' @returns A matrix of median pre-fishery biomass for each year (rows) and 
#' retrospective peel (columns) 
#'

compute.retrospective.matrix <- function(n.peels, dir_mcmc_out) {

    base <- readr::read_csv(here::here(dir_mcmc_out, "PFRBiomass.csv"), 
                            col_names=FALSE, show_col_types = FALSE) |>
                summarise(across(everything(), median)) |>
                as.numeric()

    nyr <- length(base)

    # Aggregate biomass estimates for each retrospective run
    annual.estimate <- matrix(NA, nyr, n.peels+1)
    annual.estimate[1:nyr, 1] <- base

    for(i in n.peels:1){

        print(i)
        dir_retro_mcmc_out <- here::here("retrospectives", paste0("basa-", i), "mcmc_out")
        est <- readr::read_csv(here::here(dir_retro_mcmc_out, "PFRBiomass.csv"),     
                               col_names=FALSE, show_col_types = FALSE) |>
                    summarise(across(everything(), median)) |>
                    t()        
        rownames(est) <- NULL

        for(j in 1:length(est)){
            print(j)
            annual.estimate[j, i+1] <- est[j, 1]
        }

    }

    rownames(annual.estimate) <- 1980:(1980+nyr-1)
    colnames(annual.estimate) <- c("base", paste0("-", 1:n.peels))
    
    return(annual.estimate)
}