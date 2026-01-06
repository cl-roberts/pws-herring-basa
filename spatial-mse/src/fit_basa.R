fit_basa <- function(data, par, map, DLL, do_mcmc = TRUE, lower, upper, chains, seed, iter, control, warmup) {

    # make estimation model object  
    model <- TMB::MakeADFun(
        data, par, map = map, DLL = DLL, 
        silent = TRUE, hessian = TRUE
    )

    new_seed <- seed

    # fit with ML
    converged <- 0
    while(!converged) {
        start <- pwsHerringBasa::init_tmb_params(
                chains = chains, seed = new_seed,
                start = model$par, 
                lower = lower, upper = upper
            ) |>
            unlist()    
        names(start) <- names(model$par)
        fit_ML <- nlminb(
                start = start, objective = model$fn, gradient = model$gr, 
                lower = lower, upper = upper,
                control = list(eval.max = 10000, iter.max = 1000, rel.tol = 1e-3)
            )
        converged <- fit_ML$convergence == 0
        if (!converged) {
            new_seed <- new_seed + 1
            cat("Initial model fit failed, trying again with new initial parameter values...\n")
        }
    }

    if (do_mcmc) {

        new_seed <- seed
        success <- 0
        mcmc_start_time <- Sys.time()
        while(!success) {
            # randomize initial values for MCMC chains
            inits <- pwsHerringBasa::init_tmb_params(
                chains = chains, seed = new_seed,
                start = model$env$last.par.best, 
                lower = lower, upper = upper
            )
            # run NUTS to set GHL
            fit <- tryCatch({
                tmbstan::tmbstan(
                    model, chains = chains, lower = lower, upper = upper, 
                    cores = chains, init = inits, iter = iter, 
                    control = control, warmup = warmup, seed = new_seed, algorithm = "NUTS", 
                    silent = TRUE
                )
            })
            mcmc_end_time <- Sys.time()
            success <- fit@mode == 0
            if (!success) {
                new_seed <- new_seed + 1
                if (new_seed-seed > 100) stop("MSE Iteration failed")
                cat("seed", new_seed, "\n")
                cat("MCMC failed, trying again with new initial parameter values...\n")
            }
        }
        mcmc_time <- mcmc_end_time - mcmc_start_time
        message(paste("MCMC time:", round(mcmc_time, 4), units(mcmc_time)))

        # save parameter posteriors
        mcmc_results <- as.data.frame(fit) |> 
            dplyr::select(!lp__)

        # calculate forecast
        Btilde_forecast <- apply(mcmc_results, MARGIN = 1, FUN = \(x) model$report(x)$Btilde_forecast) |>
            median()

        # recruits
        annual_age0devs <- mcmc_results |>
            dplyr::select(dplyr::contains("annual_age0devs")) |>
            apply(MARGIN = 2, FUN = median)

        return(
            list(mcmc_results = mcmc_results, 
                Btilde_forecast = Btilde_forecast,
                annual_age0devs = annual_age0devs,
                model = model)
        )

    }

    if (!do_mcmc) return(list(model = model))

}