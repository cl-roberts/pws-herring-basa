################################################################################

# Run MSE 

# CL Roberts

# This script simulates from the operating model (om) and fits simulated data 
# using estimation model (em)

# The following script options can be used to control the MSE: 

# em: string, child directory (within spatial_mse/src/em) containing estimation model used in MSE assessments
# num_mse_iters: integer, number of mse iterations to run
# num_sim_years: integer, number of years to simulate in each MSE iteration
# KI_recruitment_factor: numeric, controls relative size of KI to PWS log mean recruits
# R_corr: numeric, correlation between KI and PWS recruitments
# KI_milt_add_var: variance for KI milt index
# KI_lambda: numeric, proportion that move from KI to PWS
# PWS_lambda: numeric, proportion that move from PWS to KI
# lambda_ages: integer vector, ones denote ages that move between PWS and KI. Model uses new spawners if all zero
# KI_hcr: function, arbitrary knife-edged harvest control rule for kayak island
# probability_near: function to calculate the probability that assess biomass is near simulated biomass
# write_outputs: logical, write outputs to file if true

################################################################################

# start Rscript processes with high priority
# shell('cmd.exe /c start "" /high "C:\Program Files\R\R-4.4.1\bin\Rscript.exe"')

rm(list = ls())

start_time <- Sys.time()

#### script controls ####

#### mse controls ----

em <- "spatialmodel"
# num_mse_iters <- 55
# num_sim_years <- 10
num_mse_iters <- 10
num_sim_years <- 10
KI_recruitment_factor <- 3/4
R_corr <- 0.25
KI_milt_add_var <- 0.8
lambda_states <- c(0, 0.01, 0.15) 
# lambda_states <- 0.15 
lambda_ages_states <- list(
    all_ages = rep(1, 10),
    young_ages = c(rep(1, 5), rep(0, 5)),
    new_spawners = rep(0, 10)
)
# lambda_ages_states <- list(
#     all_ages = rep(1, 10)
# )
# KI_lambda <- .2
# PWS_lambda <- KI_lambda
# lambda_ages <- rep(1, 10)

KI_hcr <- function(T, l = 40, ghl = 2000) {
    out <- ifelse(T > l, ghl, 0)
    return(out)
}

probability_near <- function(x, vec, percent = 30) {
    p <- percent/100
    prob <- sum((vec > (x-p*x)) & (vec < (x+p*x))) / length(vec)
    return(prob)
} 

write_outputs <- TRUE
write_log <- TRUE
test_em <- FALSE

n_cores <- 10

#### model controls ----

# declare which parameters to fix
fix <- c("pk", "egg_add", "Z_0_8", "log_MeanAge0") 
# fix <- c("pk", "egg_add", "Z_0_8") 

# MCMC controls
seed <- 406
set.seed(seed)
chains <- 1
iter <- 700
warmup <- 300
control <- list(adapt_delta = 0.85, max_treedepth = 12)
timeout <- 10000

# forecast controls
forecast_controls <- list(
    recruitment_average_years = 10,
    waa_average_years = 10, 
    disease_cov_average_years = 10, 
    expected_spring_harvest = 0,
    perc_female_forecast_years = 10    
)

# ------------------------------------------------------------------------------ 

#### set up ####

# attach packages
suppressPackageStartupMessages({
    library(TMB)
    library(tmbstan)
    library(rstan)
    library(pwsHerringBasa)
    library(dplyr)
    library(here)
    library(ggplot2)
    library(tidyr)
    library(mvtnorm)
    library(foreach)
    library(doParallel)
})

theme_set(theme_bw())


# directory handling
dir_mse <- here("spatial-mse")
dir_src <- here(dir_mse, "src")
dir_em_base <- here(dir_src, "em", "base")
dir_em <- here(dir_src, "em", em)
dir_om <- here(dir_src, "om")
dir_out <- here(dir_mse, "mse_out")
if (!dir.exists(dir_out)) dir.create(dir_out)
# dir_results <- here(dir_out, "results")
# if (!dir.exists(dir_results)) dir.create(dir_results)
# dir_figures <- here(dir_out, "figures")
# if (!dir.exists(dir_figures)) dir.create(dir_figures)

# source functions
lapply(here(dir_src, list.files(path = here(dir_src), pattern = ".R")), source)


#### MSE start up ----

est_runtime <- ceiling((70*(num_sim_years + 1)/60))

if(write_log) {
    log_file <- here(dir_out, em, "mse-log.txt")
    sink(log_file)
    cat("\n\n+", rep("-", 60), "+\n")
    cat("Starting up MSE analysis ------")
    # cat("Estimated runtime: ", est_runtime, " minutes plus compile time")
    cat("\nStartup time:", as.character(Sys.time()))
    cat("\nEstimation model:", em)
    cat("\nNumber of states of nature:", 1+(length(lambda_ages_states)*(length(lambda_states)-1)))
    cat("\nMSE iters:", num_mse_iters)
    cat("\nNumber of simulation years:", num_sim_years)    
    cat("\n+", rep("-", 60), "+\n")
    sink()
}

# ------------------------------------------------------------------------------ 

#### compile models ####

em_DLL <- paste0("PWS_ASA_em_", em)
em_base_DLL <- "PWS_ASA_em_base"
om_DLL <- "PWS_ASA_om"

## estimation model

## base estimation model
if(em_base_DLL %in% names(getLoadedDLLs())) {
    dyn.unload(dynlib(here(dir_em_base, em_base_DLL)))
}
compile(here(dir_em_base, paste0(em_base_DLL, ".cpp")))
dyn.load(dynlib(here(dir_em_base, em_base_DLL)))

## estimation model
if(em_DLL %in% names(getLoadedDLLs())) {
    dyn.unload(dynlib(here(dir_em, em_DLL)))
}
compile(here(dir_em, paste0(em_DLL, ".cpp")))
dyn.load(dynlib(here(dir_em, em_DLL)))

## operating model 
if(om_DLL %in% names(getLoadedDLLs())) {
    dyn.unload(dynlib(here(dir_om, om_DLL)))
}
compile(here(dir_om, paste0(om_DLL, ".cpp")))
dyn.load(dynlib(here(dir_om, om_DLL)))

# ------------------------------------------------------------------------------ 

#### model data ####

#### estimation model data ----

# model data
PWS_ASA <- read.data.files(dir_em)$"PWS_ASA.dat"  

# recruitment and natural mortality deviate formulations
PWS_ASA_covariate <- read.data.files(dir_em)$"PWS_ASA_covariate.ctl"

# disease data
PWS_ASA_disease <- read.data.files(dir_em)$"PWS_ASA_disease.dat"

# raw agecomp sample sizes   
agecomp_samp_sizes <- read.data.files(dir_em)$"agecomp_samp_sizes.txt"

# identify missing data in raw sample sizes with -9's
seine_ac <- PWS_ASA$seine_age_comp
seine_missing <- apply(seine_ac, MARGIN = 1, FUN = \(x) any(x == -9))
spawn_ac <- PWS_ASA$spawn_age_comp
spawn_missing <- apply(spawn_ac, MARGIN = 1, FUN = \(x) any(x == -9))
agecomp_samp_sizes$seine_sample_size[seine_missing,] <- -9
agecomp_samp_sizes$spawn_sample_size[spawn_missing,] <- -9

# KI data
KI_data <- list(
    KI_mdm = as.matrix(read.csv(here(dir_mse, "data", "KI_mdm.csv"))$KI_mdm),
    KI_spawn_ess = as.matrix(read.csv(here(dir_mse, "data", "KI_spawn_ess.csv"))$KI_spawn_ess),
    KI_spawn_age_comp = as.matrix(read.csv(here(dir_mse, "data", "KI_spawn_age_comp.csv"))[,-1])
)

# will later calculate ESS's from raw sample sizes
PWS_ASA_ESS <- agecomp_samp_sizes
names(PWS_ASA_ESS) <- c("seine_ess", "spawn_ess", "vhsv_ess", "ich_ess")

# collate data into single list to pass to model
model_data <- c(
    PWS_ASA, PWS_ASA_covariate, PWS_ASA_disease, 
    agecomp_samp_sizes, PWS_ASA_ESS, KI_data,
    forecast_controls
)

#### global variables ----

Y <- model_data$nyr
start.year <- 1980
curr.year <- start.year+Y
A <- model_data$nage

# ------------------------------------------------------------------------------ 

#### model parameters ####

## estimation model parameters

fixed_pars <- list(
    pk = 0.75, egg_add = 0.4, Z_0_8 = 0.25, 
    log_MeanAge0 = 6.2 
    # sigma_age0devs = 0
)

phase1_pars <- list(
    # log_MeanAge0 = 6.2, 
    annual_age0devs = rep(0, Y-1), 
    log_juvenile_q = 4.22, 
    loginit_pop = c(6.35,  5.66,  5.92,  6.74,  4.74)
)
phase2_pars <- list(
    Z_9 = 0.93, beta_mortality = rep(0.2, 3), 
    logmdm_c = 5.87, adfg_hydro_q = -0.38, pwssc_hydro_q = -0.21
)
phase3_pars <- list(
    VHSV_age3_4_mort_93 = 0.08, ICH_age5_8_mort_93 = 0.22,
    mat_age3 = 0.60, mat_age4 = 0.99
)
phase4_pars <- list(
    seine_selex_alpha = 3.66, seine_selex_beta = 2.83
)
phase5_pars <- list(
    milt_add_var = 0.33,
    adfg_hydro_add_var = 0.30, pwssc_hydro_add_var = 0.32,
    juvenile_overdispersion = 1.00
)

model_parameters <- c(
    fixed_pars, phase1_pars, phase2_pars, 
    phase3_pars, phase4_pars, phase5_pars
)


# ------------------------------------------------------------------------------ 

#### other model attributes ####

## the map object is used to fix parameters
map <- as.list(rep(factor(NA), length(fix)))
names(map) <- fix

## parameter lower bounds
lower <- c(
    log_MeanAge0 = -10, 
    annual_age0devs = rep(-10, Y-1), 
    log_juvenile_q = -5,
    loginit_pop = rep(3, 5), 
    Z_9 = 0.3,  
    beta_mortality = rep(-30, 3), 
    logmdm_c = 2.3,  adfg_hydro_q = -5, pwssc_hydro_q = -5,     
    VHSV_age3_4_mort_93 = 0, ICH_age5_8_mort_93 = 0, 
    mat_age3 = 0.01, mat_age4 = 0.3,     
    seine_selex_alpha = 3, seine_selex_beta = 1,
    milt_add_var = 0.01, 
    adfg_hydro_add_var = 0.01, pwssc_hydro_add_var = 0.01,
    juvenile_overdispersion = 0.01
) 

lower <- lower[!(names(lower) %in% fix)]

## parameter upper bounds
upper <- c(
    log_MeanAge0 = 10, 
    annual_age0devs = rep(10, Y-1), 
    log_juvenile_q = 8,
    loginit_pop = rep(8, 5), 
    Z_9 = 1.6, 
    beta_mortality = rep(30, 3), 
    logmdm_c = 9, adfg_hydro_q = 5, pwssc_hydro_q = 5,     
    VHSV_age3_4_mort_93 = 5, ICH_age5_8_mort_93 = 5, 
    mat_age3 = 0.9, mat_age4 = 1, 
    seine_selex_alpha = 5, seine_selex_beta = 7,
    milt_add_var = 8,  
    adfg_hydro_add_var = 0.7, pwssc_hydro_add_var = 0.6,
    juvenile_overdispersion = 4
)

upper <- upper[!(names(upper) %in% fix)]


# ------------------------------------------------------------------------------

#### non-base estimation model assets ####

em_data <- model_data 
em_data$em <- em
em_parameters <- model_parameters 
em_lower <- lower 
em_upper <- upper 
em_map <- map

if (em == "gengamma") {
    em_parameters <- c(em_parameters, Q = .33)
    em_lower <- c(lower, Q = -2)
    em_upper <- c(upper, Q = 12)
}

if (em == "mvtweedie") {
    em_parameters <- c(em_parameters, seine_log_phi = log(3), seine_psi_star = 0, spawn_log_phi = log(3), spawn_psi_star = 0)
    em_lower <- c(lower, seine_log_phi = -5, seine_psi_star = -5, spawn_log_phi = -5, spawn_psi_star = -5)
    em_upper <- c(upper, seine_log_phi = 5, seine_psi_star = 5, spawn_log_phi = 5, spawn_psi_star = 5)
}

if (em == "spatialmodel") {

    em_data$lambda_ages <- lambda_ages_states[[1]]

    em_data$KI_mdm <- KI_data$KI_mdm
    em_data$KI_seine_yield <- matrix(0, nrow = Y, ncol = 1)          
    em_data$KI_seine_age_comp <- matrix(-9, nrow = Y, ncol = A)      
    em_data$KI_seine_ess <- matrix(0, nrow = Y, ncol = 1)            

    # em_data$KI_spawn_sample_size <- as.matrix(read.csv(here(dir_em, "KI_spawn_sample_size.csv"))$KI_spawn_sample_size)
    em_data$KI_spawn_age_comp <- KI_data$KI_spawn_age_comp
    em_data$KI_spawn_ess <- KI_data$KI_spawn_ess
    # em_parameters$KI_loginit_pop <- rep(1.5, 5)
    # em_parameters$KI_log_MeanAge0 <-  log(exp(model_parameters$log_MeanAge0)*KI_recruitment_factor)
    em_parameters$KI_annual_age0devs <- rep(0, Y-1) 
    # em_map$KI_annual_age0devs <- rep(factor(NA), Y-1)
    # em_parameters$logit_lambda <- -3
    # em_map$logit_lambda <- factor(NA) 
    em_parameters$lambda <- 0
    em_map$lambda <- factor(NA) 
    em_parameters$KI_milt_add_var <- KI_milt_add_var
    em_map$KI_milt_add_var <- factor(NA) 
    # em_map$log_MeanAge0 <- factor(NA)
    # em_map$KI_log_MeanAge0 <- factor(NA)
    em_parameters$sigma_age0devs <- 2
    em_map$sigma_age0devs <- factor(NA)
    em_parameters$KI_recruitment_factor <- KI_recruitment_factor
    em_map$KI_recruitment_factor <- factor(NA)
    em_parameters$R_corr <- R_corr
    em_map$R_corr <- factor(NA)
    em_lower <- c(
        lower, 
        # KI_loginit_pop = rep(1.0, 5), 
        # KI_log_MeanAge0 = -5, 
        KI_annual_age0devs = rep(-10, Y-1)
        # logit_lambda = -10
        # KI_milt_add_var = 0.01
    )
    em_upper <- c(
        upper, 
        # KI_loginit_pop = rep(7.0, 5), 
        # KI_log_MeanAge0 = 5, 
        KI_annual_age0devs = rep(10, Y-1)
        # logit_lambda = 10
        # KI_milt_add_var = 8
    )
}

# ------------------------------------------------------------------------------


#### run MSE simulation ####

#### Pseudocode ----

# calculate effective sample sizes
# fit om with MCMC to generate draws from posterior to simulate from
# simulate recruitment devs for all runs and years
# for each state of nature:
#    for run in number of mse iterations:
#      for year y in simulation years:
#          - fit estimation model
#          - apply harvest strategy and calculate GHL
#          - update operating model data
#              + fit OM to true data if first simulation year
#              + use parameter estimates as true parameters in simulations
#          - update model population dynamics
#          - simulate survey data
#          - sample new ess's for simulate agecomp data 
#          - update estimation model data
#          - update estimation model parameters
#          - y = y+1
#      run = run + 1

#### calculate ess's ----

ess <- calculate_ess(
    data = model_data, max_iter = 10, 
    parameters = model_parameters, map = map, DLL = em_base_DLL,
    lower = lower, upper = upper, 
    control = list(eval.max = 10000, iter.max = 1000, rel.tol = 1e-10)
)
# if (em == "spatialmodel") {
#     ess <- calculate_ess_spatialmodel(
#         data = em_data, max_iter = 10, 
#         parameters = em_parameters, map = em_map, DLL = em_DLL,
#         lower = em_lower, upper = em_upper, 
#         control = list(eval.max = 10000, iter.max = 1000, rel.tol = 1e-10)
#     )
# } else {
#     ess <- calculate_ess(
#         data = model_data, max_iter = 10, 
#         parameters = model_parameters, map = map, DLL = em_base_DLL,
#         lower = lower, upper = upper, 
#         control = list(eval.max = 10000, iter.max = 1000, rel.tol = 1e-10)
#     )
# }

# save ESS's to model data for fitting BASA
model_data$seine_ess <- ess$seine_ess
model_data$spawn_ess <- ess$spawn_ess
em_data$seine_ess <- ess$seine_ess
em_data$spawn_ess <- ess$spawn_ess

# if (em == "spatialmodel") {
#     model_data$KI_spawn_ess <- ess$KI_spawn_ess
#     em_data$KI_spawn_ess <- ess$KI_spawn_ess
# }

state_of_nature <- 0

for (cc in seq(lambda_ages_states)) {

    lambda_ages <- lambda_ages_states[[cc]]
    if (em == "spatialmodel") {
        em_data$lambda_ages <- lambda_ages_states[[cc]]
    }
    movement_state <- names(lambda_ages_states)[[cc]] 

    for (ll in seq(lambda_states)) {

        if (all(lambda_states[ll] == 0, cc > 1)) next 

        if (em == "spatialmodel") {
            em_parameters$lambda <- lambda_states[ll]
        }

        PWS_lambda <- lambda_states[ll]
        KI_lambda <- PWS_lambda
        
        state_start_time <- Sys.time()
        
        state_of_nature <- state_of_nature + 1

        if(write_log) {
            sink(log_file, append = TRUE)
            cat("\nStarting MSE", state_of_nature, "with the following state of nature:")
            cat("\n  ")
            cat("-KI recruitment factor:", KI_recruitment_factor)
            cat("\n  ")
            cat("-KI recruitment correlation:", R_corr)
            cat("\n  ")
            cat("-KI lambda:", KI_lambda)
            cat("\n  ")
            cat("-PWS lambda:", PWS_lambda)
            cat("\n  ")
            cat("-Ages that move:", paste(lambda_ages, collapse = " "))
            sink()
        }

        dir_results <- here(dir_out, em, paste("lambda", PWS_lambda, movement_state, sep = "_"))
        if (!dir.exists(dir_results)) dir.create(dir_results)


        #### operating model data ----

        new_waa <- apply(model_data$waa, 2, mean)
        new_perc_female <- mean(model_data$perc_female)

        simulation_inputs <- list(
            new_waa = new_waa,                                       # waa used in simulations
            new_perc_female = new_perc_female,                       # perc female in simulations
            lambda_ages = lambda_ages,                               # ages that move
            KI_seine_yield = matrix(0, nrow = Y, ncol = 1),          # catch for KI (none historically)
            KI_seine_age_comp = matrix(-9, nrow = Y, ncol = A),      # KI catch age comp
            KI_seine_ess = matrix(0, nrow = Y, ncol = 1),            # ess for KI catch agecomp
            do_simulation = 0                                        # boolean for controlling simulation
        )

        if (em == "spatialmodel") {
            simulation_inputs$KI_spawn_ess <- em_data$KI_spawn_ess
        }

        om_data <- c(
            PWS_ASA, PWS_ASA_covariate, PWS_ASA_disease, 
            agecomp_samp_sizes, PWS_ASA_ESS, KI_data,
            simulation_inputs
        )

        #### for testing em ----
        if(test_em) {

            stop("estimation model ready for testing")

            model_iteration <- TMB::MakeADFun(
                data = em_data, 
                parameters = em_parameters, map = em_map, DLL = em_DLL, 
                silent = FALSE, hessian = FALSE
            )

            model_iteration$report()$res

            fit_em <- fit_basa(
                data = em_data, par = em_parameters, map = em_map, DLL = em_DLL, 
                lower = em_lower, upper = em_upper, chains = chains, seed = seed, 
                iter = 2000, control = list(adapt_delta = .85, max_treedepth = 10), 
                warmup = 700
            )

            colMeans(fit_em$mcmc_results)

            ggplot(fit_em$mcmc_results) +
                geom_histogram(aes(x = logit_lambda))

            Btilde_posterior <- matrix(nrow = Y, ncol = (iter-warmup))
            spawn_agecomp_posterior <- list(length = (iter-warmup))
            for (i in 1:(iter-warmup)) {
                Btilde_posterior[,i] <- fit_em$model$report(fit_em$mcmc_results[i,])$Btilde_y
                spawn_agecomp_posterior[[i]] <- matrix(fit_em$model$report(fit_em$mcmc_results[i,])$spawn_age_comp, nrow = Y, ncol = A)
            }

            # biomass
            Btilde_summary <- Btilde_posterior |>
                apply(MARGIN = 1, FUN = quantile, probs = c(.025, .25, .5, .75, .975)) |>
                t() |>
                as.data.frame() |>
                cbind(data.frame(Year = 1980:2024))

            ggplot(Btilde_summary, aes(x = Year)) +
                geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), alpha = .25) +
                geom_ribbon(aes(ymin = `25%`, ymax = `75%`), alpha = .25) +
                geom_line(aes(y = `50%`)) 

            # spawn age comp
            spawn_years_index <- which(em_data$spawn_ess[-em_data$nyr,] != -9)
            spawn_years <- 1980+spawn_years_index-1
            spawn_agecomp_median <- spawn_agecomp_posterior |>
                lapply(FUN = \(x) as.matrix(x[spawn_years_index,4:10])) |>
                simplify2array() |>
                apply(1:2, median) |>
                as.data.frame() |>
                setNames(paste0("age_", 3:9)) |>
                cbind(data.frame(Year = spawn_years)) |>
                pivot_longer(cols = !Year, names_to = "age", values_to = "median") 
            ggplot(spawn_agecomp_median) +
                geom_col(aes(x = age, y = median)) +
                facet_wrap(~ Year, dir = "v") 

        }

        #### fit OM with MCMC to generate states of nature to simulate from ----

        # turn off simulation to fit model
        om_data$do_simulation <- 0

        # ess's for operating model data
        om_data$seine_ess <- ess$seine_ess
        om_data$spawn_ess <- ess$spawn_ess

        # add spatial model parameters to operating model
        om_parameters <- model_parameters
        om_parameters$KI_loginit_pop <- log(exp(model_parameters$loginit_pop)*KI_recruitment_factor)
        om_parameters$KI_log_MeanAge0 <- log(exp(model_parameters$log_MeanAge0)*KI_recruitment_factor)
        om_parameters$KI_annual_age0devs <- rep(0, Y-1)
        om_parameters$KI_lambda <- KI_lambda                       
        om_parameters$PWS_lambda <- PWS_lambda
        om_parameters$sigma_age0devs <- 2
        om_parameters$R_corr <- R_corr
        om_parameters$KI_milt_add_var <- KI_milt_add_var

        # don't estimate spatial model parameters
        om_map <- map
        om_map$KI_loginit_pop <- rep(factor(NA), length(om_parameters$KI_loginit_pop))
        om_map$KI_log_MeanAge0 <- factor(NA)
        om_map$KI_annual_age0devs <- rep(factor(NA), length(om_parameters$KI_annual_age0devs))
        om_map$KI_lambda <- factor(NA)                       
        om_map$PWS_lambda <- factor(NA)
        om_map$sigma_age0devs <- factor(NA)
        om_map$R_corr <- factor(NA)
        om_map$KI_milt_add_var <- factor(NA)

        # fit operating model
        fit_om <- fit_basa(
            data = om_data, par = om_parameters, map = om_map, DLL = om_DLL, 
            do_mcmc = num_mse_iters>=1,
            lower = lower, upper = upper, chains = chains, seed = seed, 
            iter = iter, 
            warmup = warmup,
            control = list(adapt_delta = 0.95, max_treedepth = 16)
        )

        if (num_mse_iters < 1) stop("Quitting without running MSE simulations")

        # OM states of nature
        om_sample <- fit_om$mcmc_results |>
            slice_sample(n = num_mse_iters, replace = FALSE)

        om_states <- om_sample |>
            split(seq(num_mse_iters)) |>
            lapply(FUN = \(x) setNames(x, gsub("\\[\\d+\\]", "", colnames(x)))) |>
            lapply(FUN = \(x) c(fixed_pars, split(unlist(x), colnames(x)))[names(model_parameters)]) |>
            lapply(FUN = \(x) c(x, om_parameters[!(names(om_parameters) %in% names(x))]))

        # assumed recruitment parameters used in simulations
        low_regime_index <- 1993:(curr.year-1) - 1980
        low_regime_Rbar <- mean(fit_om$annual_age0devs[low_regime_index])
        low_regime_sigmaR <- sd(fit_om$annual_age0devs[low_regime_index])   

        # create mean vector and variance-covariance matrix for simulating recruitments
        R_mu <- c(low_regime_Rbar, low_regime_Rbar)
        R_vcov <- diag(2)*low_regime_sigmaR
        R_vcov[lower.tri(R_vcov) | upper.tri(R_vcov)] <- R_corr * low_regime_sigmaR^2

        # turn simulation back on
        om_data$do_simulation <- 1


        #### pre-create recruitment devs for all simulation runs ----

        # recruitment deviates in simulation
        new_R_deviates <- vector(mode = "list", length = num_mse_iters) |>
            lapply(FUN = \(x) rmvnorm(n = num_sim_years, mean = R_mu, sigma = R_vcov))

        # KI initial population sizes in simulation
        KI_loginit_pops <- om_states |>
            lapply(FUN = \(x) x$KI_log_MeanAge0+rnorm(5, low_regime_Rbar, low_regime_sigmaR)-(1:5)*x$Z_0_8)


        # generate hindcast recruit devs for KI in simulations
        # such that they are correlated with estimated PWS recruit devs
        hindcast_KI_R_deviates <- om_states |>
            lapply(FUN = \(x) {
                    PWS_devs <- as.matrix(scale(x$annual_age0devs[1:(Y-1)]))
                    y <- rnorm(Y-1)
                    e <- residuals(lm(y ~ PWS_devs))
                    z <- R_corr * PWS_devs[,1] + sqrt(1 - R_corr^2) * scale(e)[,1]
                    return(low_regime_Rbar + low_regime_sigmaR * z)
                    }
                )

        #### Run MSE with parallel processing ----

        # names of objects to return
        names_to_return <- c(
            "biomass_sim", "recruits_sim", "mdm_sim",  
            "KI_biomass_sim", "KI_recruits_sim", "KI_mdm_sim", 
            "biomass_em", "recruits_em", "mdm_em",
            "simulated_mdm_data", "simulated_KI_mdm_data",
            "performance_metrics", "successful_iterations"
        )

        cl <- makeCluster(n_cores)
        registerDoParallel(cl)

        mse <- foreach (run = 1:num_mse_iters, .errorhandling = "remove") %dopar% {

            #### MSE iteration start up -----
            iter_time_start <- Sys.time()
            message("Starting MSE iteration ", run)

            ## load model dlls
            dyn.load(TMB::dynlib(here::here(dir_em, em_DLL)))
            dyn.load(TMB::dynlib(here::here(dir_om, om_DLL)))

            ## OM simulation state is a draw from posterior
            om_parameters <- om_states[[run]]

            # simulate KI initial state and recruitment hindcast
            om_parameters$KI_loginit_pop <- KI_loginit_pops[[run]]
            om_parameters$KI_annual_age0devs <- hindcast_KI_R_deviates[[run]]

            # project forward one MSE iteration
            if (em == "spatialmodel") message(paste0("em: lambda=", em_parameters$lambda))
            suppressWarnings({
                mse_iteration <- tryCatch(
                    run_mse_iteration(
                        num_sim_years = num_sim_years, new_R_deviates = new_R_deviates[[run]],
                        om_data = om_data, om_parameters = om_parameters, om_map = om_map, om_DLL = om_DLL, 
                        em_data = em_data, em_parameters = em_parameters, em_map = em_map, em_DLL = em_DLL, 
                        lower = em_lower, upper = em_upper,
                        chains = chains, seed = seed, iter = iter, control = control, warmup = warmup,
                        timeout = timeout
                    )
                )
            })


            om <- mse_iteration$om
            fit_em <- mse_iteration$fit_em
            om_data_sim <- mse_iteration$om_data
            mean_ghl <- mse_iteration$mean_ghl
            num_KI_fisheries <- mse_iteration$num_KI_fisheries
            simulated_mdm_data <- mse_iteration$simulated_mdm_data
            simulated_KI_mdm_data <- mse_iteration$simulated_KI_mdm_data

            #### save simulation run results -----

            ## simulated dynamics

            om_fits <- om$report()

            # PWS
            biomass_sim <- om_fits$Btilde_y[-(Y+num_sim_years)]
            recruits_sim <- om_fits$age_0[-(Y+num_sim_years)]
            mdm_sim <- om_fits$That_y[-(Y+num_sim_years)]

            # KI
            KI_biomass_sim <- om_fits$KI_Btilde_y[-(Y+num_sim_years)]
            KI_recruits_sim <- om_fits$KI_age_0[-(Y+num_sim_years)]
            KI_mdm_sim <- om_fits$KI_That_y[-(Y+num_sim_years)]

            ## fit biomass
            biomass_em_posterior <- fit_em$mcmc_results |>
                apply(MARGIN = 1, FUN = \(x) fit_em$model$report(x)$Btilde_y) |>
                as.data.frame() |>
                dplyr::slice_tail(n = 1) |>
                unlist()

            biomass_em <- fit_em$mcmc_results |>
                apply(MARGIN = 1, FUN = \(x) fit_em$model$report(x)$Btilde_y) |>
                apply(MARGIN = 1, FUN = median) 

            ## fit recruits
            recruits_em <- fit_em$mcmc_results |>
                apply(MARGIN = 1, FUN = \(x) fit_em$model$report(x)$age_0) |>
                apply(MARGIN = 1, FUN = median) 

            ## fit mdm
            mdm_em <- fit_em$mcmc_results |>
                apply(MARGIN = 1, FUN = \(x) fit_em$model$simulate(x)$That_y_pp) |>
                apply(MARGIN = 1, FUN = median)
            
            ## performance metrics 
            true_biomass <- biomass_sim[Y:(Y+num_sim_years-1)]
            assessed_biomass <- median(biomass_em_posterior)
            Btilde_ratio <- assessed_biomass / rev(true_biomass)[1]
            Btilde_prob <- probability_near(rev(true_biomass)[1], biomass_em_posterior, percent = 30) 
            yield_difference <- true_biomass*hcr(true_biomass) - assessed_biomass*hcr(assessed_biomass)
            mean_yield_difference <- mean(yield_difference)
            mean_lost_yield <- sum(yield_difference[assessed_biomass < true_biomass])/num_sim_years

            performance_metrics <- data.frame(
                Btilde_ratio = Btilde_ratio,
                Btilde_prob = Btilde_prob,
                mean_ghl = mean_ghl,
                mean_yield_difference = mean_yield_difference,
                mean_lost_yield = mean_lost_yield,
                num_KI_fisheries = num_KI_fisheries
            )

            #### end of MSE iteration -----

            iter_time_end <- Sys.time()
            iter_time <- iter_time_end - iter_time_start
            message("MSE iteration ", run, " runtime: ", round(iter_time, 4), " ", units(iter_time))
            message(rep(" ", 10), "- Btilde_ratio: ", round(Btilde_ratio, 3))
            message(rep(" ", 10), "- Btilde_prob: ", round(Btilde_prob, 3))
            message(rep(" ", 10), "- mean_ghl: ", round(mean_ghl))
            message(rep(" ", 10), "- mean_yield_difference: ", round(mean_yield_difference), "\n")

            to_return <- list(
                    biomass_sim = biomass_sim, recruits_sim = recruits_sim, mdm_sim = mdm_sim,
                    KI_biomass_sim = KI_biomass_sim, KI_recruits_sim = KI_recruits_sim, KI_mdm_sim = KI_mdm_sim,
                    biomass_em = biomass_em, recruits_em = recruits_em, mdm_em = mdm_em,
                    simulated_mdm_data = simulated_mdm_data, simulated_KI_mdm_data = simulated_KI_mdm_data, 
                    performance_metrics = t(performance_metrics), successful_iterations = as.matrix(run)
                ) |>
                setNames(paste(names_to_return, run, sep = "_"))

            return(to_return)

        }

        stopCluster(cl)

        mse_unlisted <- mse |>
            unlist(recursive = FALSE)

        mse_results <- split(mse_unlisted, sub("_[0-9]+$", "", names(mse_unlisted))) |>
            lapply(FUN = \(x) do.call(cbind, x)) 

        mse_results$performance_metrics <- t(mse_results$performance_metrics)

        # mse_results$performance_metrics[,"Btilde_ratio"] |> mean()

        #### summarize MSE results ####

        years_df <- data.frame(Year = 1979 + 1:(Y+num_sim_years-1))

        summarize_simulations <- function(df, label) {
            quants <- c(.025, .25, .5, .75, .95)
            quant_names <- c("lower_95", "lower_50", "median", "upper_50", "upper_95")
            out <- df |>
                apply(MARGIN = 1, FUN = quantile, probs = quants, na.rm = TRUE) |>
                t() |>
                as.data.frame() |>
                setNames(paste(label, quant_names, sep = "_")) 
            return(out)
        }

        simulation_results <- years_df |>
            cbind(
                summarize_simulations(mse_results$biomass_sim, "Btilde"),
                summarize_simulations(mse_results$recruits_sim, "age0"),
                summarize_simulations(mse_results$mdm_sim, "T"),
                summarize_simulations(mse_results$KI_biomass_sim, "KI_Btilde"),
                summarize_simulations(mse_results$KI_recruits_sim, "KI_age0"),
                summarize_simulations(mse_results$KI_mdm_sim, "KI_T")
            )

        assessment_results <- years_df |>
            cbind(mse_results$biomass_em, mse_results$recruits_em, mse_results$mdm_em) |>
            pivot_longer(cols = !Year, names_to = c(".value", "iteration"), names_pattern = "(.*)_em_(.*)") |>
            mutate(iteration = as.integer(sub("[^0-9]+", "", iteration)))

        simulated_data <- years_df |>
            filter(Year >= curr.year) |>
            cbind(
                rbind(mse_results$simulated_mdm_data)
                ) |>
            cbind(
                rbind(mse_results$simulated_KI_mdm_data)
                ) 

        if (num_sim_years > 1) {
            simulated_data <- simulated_data |>
                    pivot_longer(cols = !Year, names_to = c(".value", "iteration"), names_pattern = "simulated_(.*)_data_(.*)") |>
                    mutate(iteration = as.integer(sub("[^0-9]+", "", iteration)))
        }

        #### plot MSE results ####

        #### biomass ----
        simulation_biomass_plot <- ggplot(simulation_results, aes(x = Year)) +
            geom_ribbon(aes(ymin = Btilde_lower_95, ymax = Btilde_upper_95), alpha = .25) +
            geom_ribbon(aes(ymin = Btilde_lower_50, ymax = Btilde_upper_50), alpha = .25) +
            geom_line(aes(y = KI_Btilde_median), color = "blue") +
            geom_vline(aes(xintercept = curr.year-1)) +
            geom_line(data = assessment_results, aes(x = Year, y = biomass, group = iteration), color = "red") 

        #### recruitment ---
        simulation_recruits_plot <- ggplot(simulation_results, aes(x = Year)) +
            geom_ribbon(aes(ymin = age0_lower_95, ymax = age0_upper_95), alpha = .25) +
            geom_ribbon(aes(ymin = age0_lower_50, ymax = age0_upper_50), alpha = .25) +
            geom_line(aes(y = KI_age0_median), color = "blue") +
            geom_vline(aes(xintercept = curr.year-1)) +
            geom_line(data = assessment_results, aes(x = Year, y = recruits, group = iteration), color = "red") 

        #### mdm ----
        simulation_mdm_plot <- ggplot(simulation_results, aes(x = Year)) +
            geom_ribbon(aes(ymin = T_lower_95, ymax = T_upper_95), alpha = .25) +
            geom_ribbon(aes(ymin = T_lower_50, ymax = T_upper_50), alpha = .25) +
            geom_line(aes(y = KI_T_median), color = "blue") +
            geom_vline(aes(xintercept = curr.year-1)) +
            geom_line(data = assessment_results, aes(x = Year, y = mdm, group = iteration), color = "red")

        if (num_sim_years > 1) {    
            simulation_mdm_plot <- simulation_mdm_plot +
                geom_point(data = simulated_data, aes(x = Year, y = mdm), color = "red", na.rm = TRUE) +
                geom_point(data = simulated_data, aes(x = Year, y = KI_mdm), color = "blue", na.rm = TRUE) 
        } 

        # kayak island biomass
        simulation_KI_biomass_plot <- ggplot(simulation_results, aes(x = Year)) +
            geom_ribbon(aes(ymin = KI_Btilde_lower_95, ymax = KI_Btilde_upper_95), alpha = .25) +
            geom_ribbon(aes(ymin = KI_Btilde_lower_50, ymax = KI_Btilde_upper_50), alpha = .25) +
            geom_line(aes(y = KI_Btilde_median)) +
            geom_vline(aes(xintercept = curr.year-1))

        # kayak island recruits
        simulation_KI_recruits_plot <- ggplot(simulation_results, aes(x = Year)) +
            geom_ribbon(aes(ymin = KI_age0_lower_95, ymax = KI_age0_upper_95), alpha = .25) +
            geom_ribbon(aes(ymin = KI_age0_lower_50, ymax = KI_age0_upper_50), alpha = .25) +
            geom_line(aes(y = KI_age0_median)) +
            geom_vline(aes(xintercept = curr.year-1))

        # kayak island milt
        simulation_KI_mdm_plot <- ggplot(simulation_results, aes(x = Year)) +
            geom_ribbon(aes(ymin = KI_T_lower_95, ymax = KI_T_upper_95), alpha = .25) +
            geom_ribbon(aes(ymin = KI_T_lower_50, ymax = KI_T_upper_50), alpha = .25) +
            geom_line(aes(y = KI_T_median)) +
            geom_vline(aes(xintercept = curr.year-1))

        state_end_time <- Sys.time()
        mse_state_time <- state_end_time - state_start_time

        if (write_log) {
            sink(log_file, append = TRUE)
            cat("\n\n****MSE", state_of_nature, "finished****") 
            cat("\nRuntime: ") 
            cat(round(mse_state_time, 4), " ", units(mse_state_time))
            cat("\nSuccessful runs: ") 
            cat(paste(mse_results$successful_iterations, collapse = ", "))
            cat("\nPerformance metrics:\n")
            print(mse_results$performance_metrics)
            cat("\n\n+", rep("-", 60), "+\n")
            sink()
        }

        if (write_outputs) {
            write.csv(mse_results$performance_metrics, here(dir_results, "performance-metrics.csv"), row.names = FALSE)
            write.csv(t(mse_results$successful_iterations), here(dir_results, "successful_iterations.csv"), row.names = FALSE)
            write.csv(simulation_results, here(dir_results, "simulation-results.csv"), row.names = FALSE)
            write.csv(assessment_results, here(dir_results, "assessment-results.csv"), row.names = FALSE)
            write.csv(simulated_data, here(dir_results, "simulated-data.csv"), row.names = FALSE) 
            ggsave(here(dir_results, "biomass.png"), simulation_biomass_plot, width = 5.5, height = 4)
            ggsave(here(dir_results, "recruits.png"), simulation_recruits_plot, width = 5.5, height = 4)
            ggsave(here(dir_results, "mdm.png"), simulation_mdm_plot, width = 5.5, height = 4)
            ggsave(here(dir_results, "KI_biomass.png"), simulation_KI_biomass_plot, width = 5.5, height = 4)
            ggsave(here(dir_results, "KI_recruits.png"), simulation_KI_recruits_plot, width = 5.5, height = 4)
            ggsave(here(dir_results, "KI_mdm.png"), simulation_KI_mdm_plot, width = 5.5, height = 4)
        }

    }
}

# ------------------------------------------------------------------------------

#### end of file ####

end_time <- Sys.time()
mse_time <- end_time - start_time
message("\n\n+", rep("-", 60), "+\n")
message(rep(" ", 15), "Total MSE runtime: ", round(mse_time, 4), " ", units(mse_time))
message("\n+", rep("-", 60), "+\n\n")

