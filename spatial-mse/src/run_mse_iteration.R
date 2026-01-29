run_mse_iteration <- function (num_sim_years, new_R_deviates, om_data, om_map, om_parameters, om_DLL, em_data, em_map, em_parameters, em_DLL, lower, upper, ...) {

    # num_sim_years = num_sim_years
    # new_R_deviates = new_R_deviates[[run]]
    # om_data = om_data
    # om_parameters = om_parameters
    # om_map = om_map
    # om_DLL = om_DLL 
    # em_data = em_data
    # em_parameters = em_parameters
    # em_map = em_map
    # em_DLL = em_DLL 
    # lower = em_lower
    # upper = em_upper
    # args <- list(
    #     chains = chains, seed = seed, iter = iter, control = control, warmup = warmup,
    #     timeout = timeout
    # )
    # y <- 1

    args <- list(...)

    ## a few objects to store MSE results 
    simulated_mdm_data <- c()
    simulated_KI_mdm_data <- c()
    mean_ghl <- 0
    num_KI_fisheries <- 0
    KI_ghl <- 0

    for (y in 1:num_sim_years) {

        message("Simulating year ", y, " ----\n")

        ## Estimation model ----

        # fit estimation model
        fit_em <- R.utils::withTimeout(
            fit_basa(
                data = em_data, par = em_parameters, map = em_map, DLL = em_DLL, 
                lower = lower, upper = upper, chains = args$chains, seed = args$seed, 
                iter = ifelse(y == num_sim_years, 2*args$iter, args$iter), 
                warmup = ifelse(y == num_sim_years, 2*args$warmup, args$warmup),
                control = args$control),
            timeout = args$timeout)
        

        # set PWS ghl based on estimation model
        ghl <- fit_em$Btilde_forecast*hcr(fit_em$Btilde_forecast)
        mean_ghl <- mean_ghl + ghl/num_sim_years

        # set KI ghl based on milt survey
        if (y > 1) {
            KI_ghl <- KI_hcr(KI_new_mdm)
            num_KI_fisheries <- sum(num_KI_fisheries, KI_ghl > 0)
        }

        ## Operating model ----

        # simulate new recruitments
        new_R_deviate <- new_R_deviates[y,]

        # PWS recruitment     
        om_parameters$annual_age0devs <- c(
            om_parameters$annual_age0devs, 
            annual_age0devs = new_R_deviate[1]
        )

        # KI recruitment     
        om_parameters$KI_annual_age0devs <- c(
            om_parameters$KI_annual_age0devs, 
            KI_annual_age0devs = new_R_deviate[2]
        )

        # fixed weight at age        
        new_waa <- om_data$new_waa

        # fixed percent female spawners        
        new_perc_female <- om_data$new_perc_female

        # fixed disease covariates for mortality
        new_disease_covs <- c(0, 0, 0)

        # ess's for simulating new agecomp data
        new_seine_ess <- ifelse(ghl > 0, sample(20:100, 1), -9)
        new_spawn_ess <- sample(20:100, 1)
        new_KI_spawn_ess <- sample(10:70, 1)
        new_KI_seine_ess <- ifelse(KI_ghl > 0, sample(10:70, 1), 0)

        # raw sample sizes (for EM with dirichlet-multinomial likelihood)
        new_seine_sample_size <- ifelse(ghl > 0, sample(1500:3000, 1), 0)
        new_spawn_sample_size <- sample(1500:3000, 1)

        # update operating model data
        om_data <- update_om_data(
            data = om_data, ghl = ghl, new_waa = new_waa, 
            new_perc_female = new_perc_female,
            new_disease_covs = new_disease_covs, 
            new_seine_ess = new_seine_ess, new_spawn_ess = new_spawn_ess,
            KI_ghl = KI_ghl, new_KI_seine_ess = new_KI_seine_ess, 
            new_KI_spawn_ess = new_KI_spawn_ess
        )

        # make om object and simulate survey data
        om_map$KI_annual_age0devs <- c(om_map$KI_annual_age0devs, factor(NA))        
        om <- TMB::MakeADFun(
            om_data, om_parameters, map = om_map, DLL = om_DLL, 
            silent = TRUE, hessian = FALSE
        )

        # save simulated survey data
        simulated_data <- om$simulate()
        new_mdm <- simulated_data$new_mdm
        KI_new_mdm <- simulated_data$KI_new_mdm
        new_juv <- simulated_data$new_juv
        new_seine_agecomp <- simulated_data$new_seine_agecomp
        new_spawn_agecomp <- simulated_data$new_spawn_agecomp

        if (y < num_sim_years) {
            simulated_mdm_data[y] <- new_mdm
            simulated_KI_mdm_data[y] <- KI_new_mdm
        }

        # update om data with simulated KI catch age comps
        new_KI_seine_age_comp <- om$report()$new_KI_seine_agecomp
        new_KI_spawn_age_comp <- simulated_data$KI_new_spawn_agecomp
        om_data$KI_seine_age_comp[Y+y,] <- new_KI_seine_age_comp

        # update estimation model data for next simulation year
        em_data <- update_em_data(
            data = em_data, ghl = ghl, new_waa = new_waa, 
            new_perc_female = new_perc_female, new_disease_covs = new_disease_covs, 
            new_seine_ess = new_seine_ess, new_spawn_ess = new_spawn_ess, 
            new_mdm = new_mdm, new_juv = new_juv, 
            new_seine_agecomp = new_seine_agecomp, new_spawn_agecomp = new_spawn_agecomp,
            new_seine_sample_size = new_seine_sample_size, new_spawn_sample_size = new_spawn_sample_size
        )

        # if spatial model, update em with KI catch and survey information
        if (em_data$em == "spatialmodel") {
            em_data$KI_mdm <- rbind(em_data$KI_mdm, KI_new_mdm)
            em_data$KI_seine_yield <- rbind(em_data$KI_seine_yield, KI_ghl)
            em_data$KI_seine_ess <- rbind(em_data$KI_seine_ess, new_KI_seine_ess)
            em_data$KI_seine_age_comp <- rbind(em_data$KI_seine_age_comp, new_KI_seine_age_comp)
            em_data$KI_spawn_ess <- rbind(em_data$KI_spawn_ess, new_KI_spawn_ess)
            em_data$KI_spawn_age_comp <- rbind(em_data$KI_spawn_age_comp, new_KI_spawn_age_comp)
            em_parameters$KI_annual_age0devs <- c(em_parameters$KI_annual_age0devs, 0)
            deviates_index <- max(which(grepl("KI_annual_age0devs", names(lower))))
            lower <- append(lower, setNames(-10, paste0("KI_annual_age0devs", Y+y-1)), after = deviates_index)
            upper <- append(upper, setNames(10, paste0("KI_annual_age0devs", Y+y-1)), after = deviates_index)
        }
        
        # update parameter dimensionality for next simulation year
        em_parameters$annual_age0devs <- c(em_parameters$annual_age0devs, 0)
        deviates_index <- max(which(grepl("annual_age0devs", names(lower)) & !grepl("KI_annual_age0devs", names(lower))))
        lower <- append(lower, setNames(-10, paste0("annual_age0devs", Y+y-1)), after = deviates_index)
        upper <- append(upper, setNames(10, paste0("annual_age0devs", Y+y-1)), after = deviates_index)

    }

    return(
        list(
            om = om, 
            fit_em = fit_em,
            om_data_sim = om_data,
            mean_ghl = mean_ghl,
            num_KI_fisheries = num_KI_fisheries,
            simulated_mdm_data = simulated_mdm_data,
            simulated_KI_mdm_data = simulated_KI_mdm_data
        )
    )

}