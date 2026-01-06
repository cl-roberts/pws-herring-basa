calculate_ess_spatialmodel <- function(data, max_iter = 10, ...) {

    args <- list(...)

    model_data <- data
    seine_samp_size <- model_data$seine_sample_size
    spawn_samp_size <- model_data$spawn_sample_size
    KI_spawn_samp_size <- model_data$KI_spawn_sample_size
    seine_ac <- model_data$seine_age_comp
    spawn_ac <- model_data$spawn_age_comp
    KI_spawn_ac <- model_data$KI_spawn_age_comp
    seine_missing <- apply(seine_ac, MARGIN = 1, FUN = \(x) any(x == -9))
    spawn_missing <- apply(spawn_ac, MARGIN = 1, FUN = \(x) any(x == -9))
    KI_spawn_missing <- apply(KI_spawn_ac, MARGIN = 1, FUN = \(x) any(x == -9))

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
            control = list(eval.max = 10000, iter.max = 1000, rel.tol = 1e-10)
        )

        cat("obj: ", fit_iteration$objective, "\n")

        # obtain estimated seine and spawner age comps
        seine_ac_est <- model_iteration$report(model_iteration$env$last.par.best)$seine_age_comp_est
        spawn_ac_est <- model_iteration$report(model_iteration$env$last.par.best)$spawn_age_comp_est
        KI_spawn_ac_est <- model_iteration$report(model_iteration$env$last.par.best)$KI_spawn_age_comp_est

        # calculate the ESS's
        seine_ess_numer <- rowSums(seine_ac_est[!seine_missing,]*(1-seine_ac_est[!seine_missing,])) 
        seine_ess_denom <- rowSums((seine_ac[!seine_missing,]-seine_ac_est[!seine_missing,])^2) 
        seine_ess <- seine_ess_numer / seine_ess_denom
            
        spawn_ess_numer <- rowSums(spawn_ac_est[!spawn_missing,]*(1-spawn_ac_est[!spawn_missing,])) 
        spawn_ess_denom <- rowSums((spawn_ac[!spawn_missing,]-spawn_ac_est[!spawn_missing,])^2) 
        spawn_ess <- spawn_ess_numer / spawn_ess_denom

        KI_spawn_ess_numer <- rowSums(KI_spawn_ac_est[!KI_spawn_missing,]*(1-KI_spawn_ac_est[!KI_spawn_missing,])) 
        KI_spawn_ess_denom <- rowSums((KI_spawn_ac[!KI_spawn_missing,]-KI_spawn_ac_est[!KI_spawn_missing,])^2) 
        KI_spawn_ess <- KI_spawn_ess_numer / KI_spawn_ess_denom

        # calculate the ratio of ESS's to raw sample sizes
        seine_ratio <- seine_ess/seine_samp_size[!seine_missing,]
        spawn_ratio <- spawn_ess/spawn_samp_size[!spawn_missing,]
        KI_spawn_ratio <- KI_spawn_ess/KI_spawn_samp_size[!KI_spawn_missing,]


        # calculate harmonic means (see Muradian et al. 2017 and Stewart & Hamel 2014)
        seine_hm <- 1 / mean(1/seine_ratio)
        spawn_hm <- 1 / mean(1/spawn_ratio)
        KI_spawn_hm <- 1 / mean(1/KI_spawn_ratio)

        # compare this harmonic mean to the previous using a convergence criteria
        if(its > 1) {

            seine_test <- 100*abs(seine_hm_last - seine_hm) / seine_hm_last
            spawn_test <- 100*abs(spawn_hm_last - spawn_hm) / spawn_hm_last
            KI_spawn_test <- 100*abs(KI_spawn_hm_last - KI_spawn_hm) / KI_spawn_hm_last

            # this criteria was arbitrarily chosen (0.09% change)
            convergence <- all(seine_test < 0.09, spawn_test < 0.09, KI_spawn_test < 0.09) 

        }

        seine_hm_last <- seine_hm
        spawn_hm_last <- spawn_hm
        KI_spawn_hm_last <- KI_spawn_hm

        # multiply the current harmonic mean by the sample size to get the new ESS
        seine_ess <- round(seine_hm*seine_samp_size)
        spawn_ess <- round(spawn_hm*spawn_samp_size)
        KI_spawn_ess <- round(KI_spawn_hm*KI_spawn_samp_size)

        # identify the missing years
        seine_ess[seine_missing,] <- -9
        spawn_ess[spawn_missing,] <- -9
        KI_spawn_ess[KI_spawn_missing,] <- -9

        # save new ess's to model data
        model_data$seine_ess <- seine_ess
        model_data$spawn_ess <- spawn_ess
        model_data$KI_spawn_ess <- KI_spawn_ess

        exit_message <- paste("ESS calculations converged after", its, "iterations")
    
    }

    if(its <= max_iter) {
        message(exit_message)
    } else {
        warning(exit_message)
    }

    return(
        list(
            "seine_ess" = seine_ess, "spawn_ess" = spawn_ess, 
            "KI_spawn_ess" = KI_spawn_ess, message = exit_message)
    )

}
