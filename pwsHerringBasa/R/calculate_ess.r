#' Calculate effective sample sizes for age composition data
#'
#' Estimates effective sample sizes for age composition data in BASA using a 
#' modified iterative reweighting procedure.
#'
#' estimate age comps --> calculate ESS's --> estimate age comps --> repeat until convergence
#' 
#' @param data List of model data passed to `TMB::MakeADFun`
#' @param max_iter The maximum number of iterations in reweighting procedure
#' @param ... Other parameters to be passed to `TMB::MakeADFun` and `stats::nlminb`
#'
#' @returns List of effective sample sizes for age composition data and a convergence exit message.
#'
#' @references \insertRef{muradian2017}{pwsHerringBasa}
#' @references \insertRef{stewart2014}{pwsHerringBasa}

calculate_ess <- function(data, max_iter = 10, ...) {

    args <- list(...)

    model_data <- data
    seine_samp_size <- model_data$seine_sample_size
    spawn_samp_size <- model_data$spawn_sample_size

    convergence <- FALSE 
    its <- 0

    while (!convergence) {

        its <- its + 1

        # stop while loop if ESS's haven't converged
        if (its > max_iter) {
            
            exit_message <- paste(
                "ESS calculation didn't converged after", max_iter, "iterations.\n",
                "Try increasing max_iter!"
            )

            break 

        }

        # fit model using agecomp sample sizes of current iteration 
        model_iteration <- TMB::MakeADFun(
            data = model_data, 
            parameters = args$parameters, map = args$map, DLL = args$DLL, 
            silent = TRUE, hessian = FALSE
        )

        # optimize model 
        fit_iteration <- stats::nlminb(
            start = model_iteration$par, 
            objective = model_iteration$fn, 
            gradient = model_iteration$gr,  
            lower = args$lower, upper = args$upper, 
            control = args$control
        )

        # obtain estimated seine and spawner age comps
        seine_ac_est <- model_iteration$report(model_iteration$env$last.par.best)$seine_age_comp_est
        spawn_ac_est <- model_iteration$report(model_iteration$env$last.par.best)$spawn_age_comp_est

        # calculate the ESS's
        seine_ess_numer <- rowSums(seine_ac_est[!seine_missing,]*(1-seine_ac_est[!seine_missing,])) 
        seine_ess_denom <- rowSums((seine_ac[!seine_missing,]-seine_ac_est[!seine_missing,])^2) 
        seine_ess <- seine_ess_numer / seine_ess_denom
            
        spawn_ess_numer <- rowSums(spawn_ac_est[!spawn_missing,]*(1-spawn_ac_est[!spawn_missing,])) 
        spawn_ess_denom <- rowSums((spawn_ac[!spawn_missing,]-spawn_ac_est[!spawn_missing,])^2) 
        spawn_ess <- spawn_ess_numer / spawn_ess_denom

        # calculate the ratio of ESS's to raw sample sizes
        seine_ratio <- seine_ess/seine_samp_size[!seine_missing,]
        spawn_ratio <- spawn_ess/spawn_samp_size[!spawn_missing,]

        # calculate harmonic means (see Muradian et al. 2017 and Stewart & Hamel 2014)
        seine_hm <- 1 / mean(1/seine_ratio)
        spawn_hm <- 1 / mean(1/spawn_ratio)

        # compare this harmonic mean to the previous using a convergence criteria
        if(its > 1) {

            seine_test <- 100*abs(seine_hm_last - seine_hm) / seine_hm_last
            spawn_test <- 100*abs(spawn_hm_last - spawn_hm) / spawn_hm_last

            # this criteria was arbitrarily chosen (0.1% change)
            convergence <- all(seine_test < 0.09, spawn_test < 0.09) 

        }

        seine_hm_last <- seine_hm
        spawn_hm_last <- spawn_hm

        # multiply the current harmonic mean by the sample size to get the new ESS
        seine_ess <- round(seine_hm*seine_samp_size)
        spawn_ess <- round(spawn_hm*spawn_samp_size)

        # identify the missing years
        seine_ess[seine_missing,] <- -9
        spawn_ess[spawn_missing,] <- -9

        # save new ess's to model data
        model_data$seine_ess <- seine_ess
        model_data$spawn_ess <- spawn_ess

        exit_message <- paste("ESS calculations converged after", its, "iterations")
    
    }

    if(its <= max_iter) {
        message(exit_message)
    } else {
        warning(exit_message)
    }

    return(list("seine_ess" = seine_ess, "spawn_ess" = spawn_ess, message = exit_message))

}