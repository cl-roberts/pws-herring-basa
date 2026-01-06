################################################################################

# Run MSE 

# CL Roberts

# This script simulates from the operating model (om) and fits simulated data 
# using estimation model (em)

################################################################################

start_time <- Sys.time()

#### mse controls ####

num_mse_iters <- 3
num_sim_years <- 5
KI_recruitment_factor <- 1/2               # controls relative size of KI to PWS log mean recruits
R_corr <- 0                                # correlation between KI and PWS recruitments
KI_lambda <- 0                             # proportion that move from KI to PWS
PWS_lambda <- 0                            # proportion that move from PWS to KI
lambda_ages <- rep(1, 10)                  # ages that move between PWS and KI. uses new spawners if all zero

# arbitrary knife-edged harvest control rule for kayak island
KI_hcr <- function(T, l = 50, ghl = 2000) {
    out <- ifelse(T > l, ghl, 0)
    return(out)
}

write_outputs <- FALSE

cat("Starting up MSE analysis ------\n")
cat("Estimated runtime: ")
cat(ceiling((25*(num_mse_iters*num_sim_years + 1)/60)))
cat(" minutes plus compile time\n\n")

#### model controls ####

## declare which parameters to fix
fix <- c("pk", "egg_add", "Z_0_8", "log_MeanAge0") 

## MCMC controls
seed <- 406
set.seed(seed)
chains <- 1
iter <- 2000
warmup <- 700
control <- list(adapt_delta = 0.95)

## forecast controls
forecast_controls <- list(
    recruitment_average_years = 10,
    waa_average_years = 10, 
    disease_cov_average_years = 10, 
    expected_spring_harvest = 0,
    perc_female_forecast_years = 10    
)


# ------------------------------------------------------------------------------ 

#### set up ####

## attach packages

library(TMB)
library(tmbstan)
library(rstan)
library(pwsHerringBasa)
library(dplyr)
library(here)
library(ggplot2)
theme_set(theme_bw())
library(tidyr)
library(mvtnorm)

## directory handling

dir_model <- here("model")
dir_mse <- here("spatial-mse")
dir_src <- here(dir_mse, "src")
dir_om <- here(dir_mse, "om")
dir_out <- here(dir_mse, "mse_out")
dir_results <- here(dir_out, "results")
dir_figures <- here(dir_out, "figures")

if (!dir.exists(dir_out)) dir.create(dir_out)
if (!dir.exists(dir_results)) dir.create(dir_results)
if (!dir.exists(dir_figures)) dir.create(dir_figures)


## source functions

lapply(here(dir_src, list.files(path = here(dir_src), pattern = ".R")), source)


# ------------------------------------------------------------------------------ 

#### compile models ####

## estimation model
if("PWS_ASA" %in% names(getLoadedDLLs())) {
    dyn.unload(dynlib(here(dir_model, "PWS_ASA")))
}
compile(here(dir_model, "PWS_ASA.cpp"))
dyn.load(dynlib(here(dir_model, "PWS_ASA")))

## operating model 

if("PWS_ASA_om" %in% names(getLoadedDLLs())) {
    dyn.unload(dynlib(here(dir_om, "PWS_ASA_om")))
}
compile(here(dir_om, "PWS_ASA_om.cpp"))
dyn.load(dynlib(here(dir_om, "PWS_ASA_om")))

# ------------------------------------------------------------------------------ 

#### model data ####

## estimation model data

# model data
PWS_ASA <- read.data.files(dir_model)$"PWS_ASA.dat"  

# recruitment and natural mortality deviate formulations
PWS_ASA_covariate <- read.data.files(dir_model)$"PWS_ASA_covariate.ctl"

# disease data
PWS_ASA_disease <- read.data.files(dir_model)$"PWS_ASA_disease.dat"

# raw agecomp sample sizes   
agecomp_samp_sizes <- read.data.files(dir_model)$"agecomp_samp_sizes.txt"

# identify missing data in raw sample sizes with -9's
seine_ac <- PWS_ASA$seine_age_comp
seine_missing <- apply(seine_ac, MARGIN = 1, FUN = \(x) any(x == -9))
spawn_ac <- PWS_ASA$spawn_age_comp
spawn_missing <- apply(spawn_ac, MARGIN = 1, FUN = \(x) any(x == -9))
agecomp_samp_sizes$seine_sample_size[seine_missing,] <- -9
agecomp_samp_sizes$spawn_sample_size[spawn_missing,] <- -9

# will later calculate ESS's from raw sample sizes
PWS_ASA_ESS <- agecomp_samp_sizes
names(PWS_ASA_ESS) <- c("seine_ess", "spawn_ess", "vhsv_ess", "ich_ess")

# collate data into single list to pass to model
model_data <- c(
    PWS_ASA, PWS_ASA_covariate, PWS_ASA_disease, 
    agecomp_samp_sizes, PWS_ASA_ESS,
    forecast_controls
)

## global variables

Y <- model_data$nyr
start.year <- 1980
curr.year <- start.year+Y
A <- model_data$nage
quants <- c(.025, .25, .5, .75, .95)
quant_names <- c("lower_95", "lower_50", "median", "upper_50", "upper_95")

## operating model data

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

om_data <- c(
    PWS_ASA, PWS_ASA_covariate, PWS_ASA_disease, 
    agecomp_samp_sizes, PWS_ASA_ESS,
    simulation_inputs
)

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
    log_MeanAge0 = 2, 
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
    log_MeanAge0 = 20, 
    annual_age0devs = rep(10, Y-1), 
    log_juvenile_q = 8,
    loginit_pop = rep(8, 5), 
    Z_9 = 1.6, 
    beta_mortality = rep(30, 3), 
    logmdm_c = 9, adfg_hydro_q = 5, pwssc_hydro_q = 5,     
    VHSV_age3_4_mort_93 = 5, ICH_age5_8_mort_93 = 5, 
    mat_age3 = 0.9, mat_age4 = 1, 
    seine_selex_alpha = 5, seine_selex_beta = 7,
    milt_add_var = 0.9,  
    adfg_hydro_add_var = 0.7, pwssc_hydro_add_var = 0.6,
    juvenile_overdispersion = 4
)

upper <- upper[!(names(upper) %in% fix)]


# ------------------------------------------------------------------------------

#### run MSE simulation ####

## Pseudocode

# calculate effective sample sizes
# fit om with MCMC to generate draws from posterior to simulate from
# for run in number of states
#   for year y in simulation years:
#       - fit estimation model
#       - apply harvest strategy and calculate GHL
#       - update operating model data
#           + fit OM to true data if first simulation year
#           + use parameter estimates as true parameters in simulations
#       - simulate operating model population dynamics
#       - simulate survey data
#       - sample new ess's for simulate agecomp data 
#       - update estimation model data
#       - update estimation model parameters
#       - y = y+1
#   run = run + 1

# calculate ess's
ess <- calculate_ess(
    data = model_data, max_iter = 10, 
    parameters = model_parameters, map = map, DLL = "PWS_ASA",
    lower = lower, upper = upper, 
    control = list(eval.max = 10000, iter.max = 1000, rel.tol = 1e-10)
)

# save ESS's to model data for fitting BASA
model_data$seine_ess <- ess$seine_ess
model_data$spawn_ess <- ess$spawn_ess


# fit OM with MCMC to generate states of nature to simulate from

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


# don't estimate spatial model parameters
om_map <- map
om_map$KI_loginit_pop <- rep(factor(NA), length(om_parameters$loginit_pop))
om_map$KI_log_MeanAge0 <- factor(NA)
om_map$KI_annual_age0devs <- rep(factor(NA), length(om_parameters$KI_annual_age0devs))
om_map$KI_lambda <- factor(NA)                       
om_map$PWS_lambda <- factor(NA)

# fit operating model
fit_om <- fit_basa(
    data = om_data, par = om_parameters, map = om_map, DLL = "PWS_ASA_om", 
    do_mcmc = num_mse_iters>=1,
    lower = lower, upper = upper, chains = chains, seed = seed, iter = iter, 
    control = control, warmup = warmup
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

#### set up objects to store MSE results ####

# simulation
biomass_sim <- matrix(nrow = Y+num_sim_years-1, ncol = num_mse_iters)
recruits_sim <- matrix(nrow = Y+num_sim_years-1, ncol = num_mse_iters)
mdm_sim <- rbind(
    apply(om_sample, MARGIN = 1, FUN = \(x) fit_om$model$report(x)$That_y),
    matrix(NA, nrow = num_sim_years-1, ncol = num_mse_iters)
)
juv_sim <- rbind(
    apply(om_sample, MARGIN = 1, FUN = \(x) fit_om$model$report(x)$Jhat_y),
    matrix(NA, nrow = num_sim_years-1, ncol = num_mse_iters)
)
seine_agecomp_sim <- apply(om_sample, MARGIN = 1, FUN = \(x) data.frame(seine_agecomp = fit_om$model$report(x)$seine_age_comp_est))
spawn_agecomp_sim <- apply(om_sample, MARGIN = 1, FUN = \(x) data.frame(spawn_agecomp = fit_om$model$report(x)$spawn_age_comp_est))

KI_biomass_sim <- matrix(nrow = Y+num_sim_years-1, ncol = num_mse_iters)
KI_recruits_sim <- matrix(nrow = Y+num_sim_years-1, ncol = num_mse_iters)
KI_mdm_sim <- rbind(
    apply(om_sample, MARGIN = 1, FUN = \(x) fit_om$model$report(x)$KI_That_y),
    matrix(NA, nrow = num_sim_years-1, ncol = num_mse_iters)
)

# assessment at end of each simulation
biomass_em <- matrix(nrow = Y+num_sim_years-1, ncol = num_mse_iters)
recruits_em <- matrix(nrow = Y+num_sim_years-1, ncol = num_mse_iters)
mdm_em <- matrix(nrow = Y+num_sim_years-1, ncol = num_mse_iters)
juv_em <- matrix(nrow = Y+num_sim_years-1, ncol = num_mse_iters)
seine_agecomp_em <- vector(mode = "list", length = num_mse_iters)
spawn_agecomp_em <- vector(mode = "list", length = num_mse_iters)

# performance metrics
Btilde_ratio <- c()
mean_ghl <- rep(0, num_mse_iters)
mean_yield_difference <- c()
mean_lost_yield <- c()

#### pre-create recruitment devs for all simulation runs ####

# recruitment deviates in simulation
new_R_deviates <- vector(mode = "list", length = num_mse_iters) |>
    lapply(FUN = \(x) rmvnorm(n = num_sim_years, mean = R_mu, sigma = R_vcov))

# KI initial population sizes in simulation
KI_loginit_pops <- om_states |>
    lapply(FUN = \(x) log(exp(x$KI_log_MeanAge0+rnorm(5, low_regime_Rbar, low_regime_sigmaR))*exp(x$Z_0_8)^(1:5)))


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

#### start MSE ####

for (run in 1:num_mse_iters){

    ## simulation run start up

    cat("\n\n+", rep("-", 20), "+\n")
    cat(rep(" ", 5), "MSE iteration", run)
    cat("\n+", rep("-", 20), "+\n\n")

    # create em data and parameter objects for mse run
    em_data <- model_data
    em_parameters <- model_parameters
    em_lower <- lower 
    em_upper <- upper

    # select om parameters based on list of states
    om_parameters <- om_states[[run]]
    om_data_sim <- om_data

    # recreate map object
    om_map_sim <- om_map

    # simulate KI initial state and recruitment hindcast
    om_parameters$KI_loginit_pop <- KI_loginit_pops[[run]]
    om_parameters$KI_annual_age0devs <- hindcast_KI_R_deviates[[run]]

    for (y in 1:num_sim_years) {

        cat("\n\n---- Simulating year", y, "----\n\n")

        ## Estimation model ----

        # fit estimation model
        fit_em <- fit_basa(
            data = em_data, par = em_parameters, map = map, DLL = "PWS_ASA", 
            lower = em_lower, upper = em_upper, chains = chains, seed = seed, iter = iter, 
            control = control, warmup = warmup
        )

        # set PWS ghl based on estimation model
        ghl <- fit_em$Btilde_forecast*hcr(fit_em$Btilde_forecast)
        mean_ghl[run] <- mean_ghl[run] + ghl/num_sim_years

        # set KI ghl based on milt survey
        KI_ghl <- KI_hcr(KI_mdm_sim[Y+y-1, run])

        ## Operating model ----

        # simulate new recruitments
        new_R_deviate <- new_R_deviates[[run]][y,]

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

        # simulate disease covariates for mortality
        new_disease_covs <- c(0, 0, 0)

        # ess's for simulating new agecomp data
        new_seine_ess <- sample(20:100, 1)
        new_spawn_ess <- sample(20:100, 1)
        new_KI_seine_ess <- sample(20:100, 1)
        
        # update operating model data
        om_data_sim <- update_om_data(
            om_data_sim, A = 10, new_waa = new_waa, 
            ghl = ghl, new_perc_female = new_perc_female,
            new_disease_covs = new_disease_covs, 
            new_seine_ess = new_seine_ess, new_spawn_ess = new_spawn_ess,
            KI_ghl = KI_ghl, new_KI_seine_ess = new_KI_seine_ess
        )

        # make om object and simulate survey data
        om_map_sim$KI_annual_age0devs <- c(om_map_sim$KI_annual_age0devs, factor(NA))        
        om <- MakeADFun(
            om_data_sim, om_parameters, map = om_map_sim, DLL = "PWS_ASA_om", 
            silent = TRUE, hessian = FALSE
        )

        # save simulated survey data
        simulated_data <- om$simulate()

        new_mdm <- simulated_data$new_mdm
        KI_new_mdm <- simulated_data$KI_new_mdm
        new_juv <- simulated_data$new_juv
        new_seine_agecomp <- simulated_data$new_seine_agecomp
        new_spawn_agecomp <- simulated_data$new_spawn_agecomp

        new_KI_seine_age_comp <- om$report()$new_KI_seine_agecomp

        # update om data with simulated KI catch age comps
        om_data_sim$KI_seine_age_comp[Y+y,] <- new_KI_seine_age_comp
        
        if (y < num_sim_years) {
            mdm_sim[Y+y, run] <- new_mdm
            KI_mdm_sim[Y+y, run] <- KI_new_mdm
            juv_sim[Y+y, run] <- new_juv
            seine_agecomp_sim[[run]] <- rbind(seine_agecomp_sim[[run]], new_seine_agecomp)
            spawn_agecomp_sim[[run]] <- rbind(spawn_agecomp_sim[[run]], new_spawn_agecomp)
        }

        # update estimation model data for next simulation year
        em_data <- update_em_data(
            data = em_data, A = 10, ghl = ghl, new_waa = new_waa, 
            new_perc_female = new_perc_female, new_disease_covs = new_disease_covs, 
            new_seine_ess = new_seine_ess, new_spawn_ess = new_spawn_ess, 
            new_mdm = new_mdm, new_juv = new_juv, 
            new_seine_agecomp = new_seine_agecomp, new_spawn_agecomp = new_spawn_agecomp
        )

        # update parameter dimensionality for next simulation year
        em_parameters$annual_age0devs <- c(em_parameters$annual_age0devs, 0)
        deviates_index <- max(which(grepl("annual_age0devs", names(em_lower))))
        em_lower <- append(em_lower, setNames(-10, paste0("annual_age0devs", Y+y-1)), after = deviates_index)
        em_upper <- append(em_upper, setNames(10, paste0("annual_age0devs", Y+y-1)), after = deviates_index)

    }

    # save simulation run results -----

    ## simulated dynamics

    # PWS
    biomass_sim[,run] <- om$report()$Btilde_y[-(om_data$nyr+y)]
    recruits_sim[,run] <- om$report()$age_0[-(om_data$nyr+y)]

    # KI
    KI_biomass_sim[,run] <- om$report()$KI_Btilde_y[-(om_data$nyr+y)]
    KI_recruits_sim[,run] <- om$report()$KI_age_0[-(om_data$nyr+y)]

    ## fit biomass
    biomass_em[,run] <- fit_em$mcmc_results |>
        apply(MARGIN = 1, FUN = \(x) fit_em$model$report(x)$Btilde_y) |>
        apply(MARGIN = 1, FUN = median) 

    ## fit recruits
    recruits_em[,run] <- fit_em$mcmc_results |>
        apply(MARGIN = 1, FUN = \(x) fit_em$model$report(x)$age_0) |>
        apply(MARGIN = 1, FUN = median) 

    ## fit mdm
    mdm_em[,run] <- fit_em$mcmc_results |>
        apply(MARGIN = 1, FUN = \(x) fit_em$model$simulate(x)$That_y_pp) |>
        apply(MARGIN = 1, FUN = median)
    
    ## fit juvenile schools
    juv_em[,run] <- fit_em$mcmc_results |>
        apply(MARGIN = 1, FUN = \(x) fit_em$model$simulate(x)$Jhat_y_pp) |>
        apply(MARGIN = 1, FUN = median) 

    ## fit seine age comps
    seine_years_index <- which(em_data$seine_ess[-em_data$nyr,] != -9)
    seine_years <- 1980+seine_years_index-1

    seine_agecomp_em[[run]] <- fit_em$mcmc_results |>
        apply(MARGIN = 1, FUN = \(x) fit_em$model$simulate(x)$seine_agecomp_pp) |>
        as.data.frame() |>
        lapply(FUN = \(x) matrix(x, nrow = em_data$nyr-1, ncol = A)[seine_years_index,4:10]) |>
        simplify2array() |>
        apply(1:2, median) |>
        as.data.frame() |>
        setNames(paste0("age_", 3:9)) |>
        cbind(data.frame(Year = seine_years))


    ## fit spawn age comps
    spawn_years_index <- which(em_data$spawn_ess[-em_data$nyr,] != -9)
    spawn_years <- 1980+spawn_years_index-1

    spawn_agecomp_em[[run]] <- fit_em$mcmc_results |>
        apply(MARGIN = 1, FUN = \(x) fit_em$model$simulate(x)$spawn_agecomp_pp) |>
        as.data.frame() |>
        lapply(FUN = \(x) matrix(x, nrow = em_data$nyr-1, ncol = A)[spawn_years_index,4:10]) |>
        simplify2array() |>
        apply(1:2, median) |>
        as.data.frame() |>
        setNames(paste0("age_", 3:9)) |>
        cbind(data.frame(Year = spawn_years))

    ## performance metrics 
    true_biomass <- biomass_sim[Y:(om_data_sim$nyr-1),run]
    assessed_biomass <- biomass_em[Y:(om_data_sim$nyr-1),run]
    Btilde_ratio[run] <- rev(assessed_biomass)[1] / rev(true_biomass)[1]
    yield_difference <- true_biomass*hcr(true_biomass) - assessed_biomass*hcr(assessed_biomass)
    mean_yield_difference[run] <- mean(yield_difference)
    mean_lost_yield[run] <- sum(yield_difference[assessed_biomass < true_biomass])/num_sim_years

}

performance_metrics <- data.frame(
    Btilde_ratio = Btilde_ratio,
    mean_ghl = mean_ghl,
    mean_yield_difference = mean_yield_difference,
    mean_lost_yield = mean_lost_yield
)

# ------------------------------------------------------------------------------

#### plot simulated data ####

## biomass
biomass_sim_df <- biomass_sim |>
    apply(MARGIN = 1, FUN = quantile, probs = c(.025, .25, .5, .75, .975)) |>
    t() |>
    as.data.frame() |>
    setNames(paste("Btilde", quant_names, sep = "_")) |>
    cbind(
        data.frame(Year = 1980:(1980+om_data_sim$nyr-2), 
                   Btilde_em = biomass_em)
        )

KI_biomass_sim_df <- KI_biomass_sim |>
    apply(MARGIN = 1, FUN = quantile, probs = c(.025, .25, .5, .75, .975)) |>
    t() |>
    as.data.frame() |>
    setNames(paste("Btilde", quant_names, sep = "_")) |>
    cbind(
        data.frame(Year = 1980:(1980+om_data_sim$nyr-2))
        )

simulation_biomass_plot <- ggplot(biomass_sim_df, aes(x = Year)) +
    geom_ribbon(aes(ymin = Btilde_lower_95, ymax = Btilde_upper_95), alpha = .25) +
    geom_ribbon(aes(ymin = Btilde_lower_50, ymax = Btilde_upper_50), alpha = .25) +
    geom_line(data = KI_biomass_sim_df, aes(y = Btilde_median), color = "blue") +
    geom_vline(aes(xintercept = curr.year)) +
    geom_line(aes(x = Year, y = Btilde_em.1), color = "red") +
    geom_line(aes(x = Year, y = Btilde_em.2), color = "red") + 
    geom_line(aes(x = Year, y = Btilde_em.3), color = "red") 

## recruitment
recruits_sim_df <- recruits_sim |>
    apply(MARGIN = 1, FUN = quantile, probs = c(.025, .25, .5, .75, .975)) |>
    t() |>
    as.data.frame() |>
    setNames(paste("age0", quant_names, sep = "_")) |>
    cbind(
        data.frame(Year = 1980:(1980+om_data_sim$nyr-2), 
                   age0_em = recruits_em)
        )

KI_recruits_sim_df <- KI_recruits_sim |>
    apply(MARGIN = 1, FUN = quantile, probs = c(.025, .25, .5, .75, .975)) |>
    t() |>
    as.data.frame() |>
    setNames(paste("age0", quant_names, sep = "_")) |>
    cbind(
        data.frame(Year = 1980:(1980+om_data_sim$nyr-2))
        )


simulation_recruits_plot <- ggplot(recruits_sim_df, aes(x = Year)) +
    geom_ribbon(aes(ymin = age0_lower_95, ymax = age0_upper_95), alpha = .25) +
    geom_ribbon(aes(ymin = age0_lower_50, ymax = age0_upper_50), alpha = .25) +
    geom_line(data = KI_recruits_sim_df, aes(y = age0_median), color = "blue") +
    geom_vline(aes(xintercept = curr.year)) +
    geom_line(aes(x = Year, y = age0_em.1), color = "red") +
    geom_line(aes(x = Year, y = age0_em.2), color = "red") + 
    geom_line(aes(x = Year, y = age0_em.3), color = "red") 

## mdm
mdm_sim_df <- mdm_sim |>
    apply(MARGIN = 1, FUN = quantile, probs = c(.025, .25, .5, .75, .975)) |>
    t() |>
    as.data.frame() |>
    setNames(paste("T", quant_names, sep = "_")) |>
    cbind(
        data.frame(Year = 1980:(1980+om_data_sim$nyr-2), 
                   That_em = mdm_em)
        )

KI_mdm_sim_df <- KI_mdm_sim |>
    apply(MARGIN = 1, FUN = quantile, probs = c(.025, .25, .5, .75, .975)) |>
    t() |>
    as.data.frame() |>
    setNames(paste("T", quant_names, sep = "_")) |>
    cbind(
        data.frame(Year = 1980:(1980+om_data_sim$nyr-2))
        )

simulation_mdm_plot <- ggplot(mdm_sim_df, aes(x = Year)) +
    geom_ribbon(aes(ymin = T_lower_95, ymax = T_upper_95), alpha = .25) +
    geom_ribbon(aes(ymin = T_lower_50, ymax = T_upper_50), alpha = .25) +
    geom_line(data = KI_mdm_sim_df, aes(y = T_median), color = "blue") +
    geom_vline(aes(xintercept = curr.year)) +
    geom_line(aes(x = Year, y = That_em.1), color = "red") +
    geom_line(aes(x = Year, y = That_em.2), color = "red") + 
    geom_line(aes(x = Year, y = That_em.3), color = "red") 


## juvenile schools
juv_sim_df <- juv_sim |>
    apply(MARGIN = 1, FUN = quantile, probs = c(.025, .25, .5, .75, .975)) |>
    t() |>
    as.data.frame() |>
    setNames(paste("J", quant_names, sep = "_")) |>
    cbind(
        data.frame(Year = 1980:(1980+om_data_sim$nyr-2), 
                   Jhat_em = juv_em)
        )

simulation_juv_plot <- ggplot(juv_sim_df, aes(x = Year)) +
    geom_ribbon(aes(ymin = J_lower_95, ymax = J_upper_95), alpha = .25) +
    geom_ribbon(aes(ymin = J_lower_50, ymax = J_upper_50), alpha = .25) +
    # geom_line(aes(y = J_median)) +
    geom_vline(aes(xintercept = curr.year)) +
    geom_line(aes(x = Year, y = Jhat_em.1), color = "red") +
    geom_line(aes(x = Year, y = Jhat_em.2), color = "red") + 
    geom_line(aes(x = Year, y = Jhat_em.3), color = "red") 


# kayak island biomass
ggplot(KI_biomass_sim_df, aes(x = Year)) +
    geom_ribbon(aes(ymin = Btilde_lower_95, ymax = Btilde_upper_95), alpha = .25) +
    geom_ribbon(aes(ymin = Btilde_lower_50, ymax = Btilde_upper_50), alpha = .25) +
    geom_line(aes(y = Btilde_median)) +
    geom_vline(aes(xintercept = curr.year))

# kayak island recruits
ggplot(KI_recruits_sim_df, aes(x = Year)) +
    geom_ribbon(aes(ymin = age0_lower_95, ymax = age0_upper_95), alpha = .25) +
    geom_ribbon(aes(ymin = age0_lower_50, ymax = age0_upper_50), alpha = .25) +
    geom_line(aes(y = age0_median)) +
    geom_vline(aes(xintercept = curr.year))

# kayak island milt
ggplot(KI_mdm_sim_df, aes(x = Year)) +
    geom_ribbon(aes(ymin = T_lower_95, ymax = T_upper_95), alpha = .25) +
    geom_ribbon(aes(ymin = T_lower_50, ymax = T_upper_50), alpha = .25) +
    geom_line(aes(y = T_median)) +
    geom_vline(aes(xintercept = curr.year))

# ------------------------------------------------------------------------------

#### plot simulated agecomps ####

## seine age comps
# the following is an array with third dimension being summary statistics
seine_agecomp_summary <- seine_agecomp_sim |>
    lapply(FUN = \(x) as.matrix(x[seine_years_index,4:10])) |>
    simplify2array() |>
    apply(1:2, quantile, prob = quants)

dimnames(seine_agecomp_summary)[[1]] <- quant_names
dimnames(seine_agecomp_summary)[[2]] <- seine_years
dimnames(seine_agecomp_summary)[[3]] <- 3:9

viridis_palette <- "H" 
color.options <- scales::viridis_pal(option = viridis_palette)(6)
seine_sim_agecomp_df <- seine_agecomp_em |> 
    sapply(FUN = \(x) select(pivot_longer(x, cols = !Year, names_to = "Age", values_to = "seine_agecomp_em"), seine_agecomp_em)) |>
    do.call(cbind, args = _) |>
    as.data.frame() |>
    cbind(
        data.frame(Year = rep(seine_years, times = rep(length(3:9), length(seine_years))),
                   Age = 3:9, 
                   fill_colors = generate.colors(nyr = length(seine_years), color.options = color.options))
    )

colnames(seine_sim_agecomp_df)[colnames(seine_sim_agecomp_df) == "seine_agecomp_em"] <- paste0("seine_agecomp_em.", 1:num_mse_iters)

seine_sim_agecomp_df$Proportion <- pivot_longer(
        as.data.frame(seine_agecomp_summary[3,,]), 
        cols = everything(), 
        names_to = "Age", 
        values_to = "Proportion"
        ) |>
    select(Proportion) |>
    unlist()

seine_sim_agecomp_df$Lower_95 <- pivot_longer(
        as.data.frame(seine_agecomp_summary[1,,]), 
        cols = everything(), 
        names_to = "Age", 
        values_to = "Proportion"
        ) |>
    select(Proportion) |>
    unlist()

seine_sim_agecomp_df$Lower_50 <- pivot_longer(
        as.data.frame(seine_agecomp_summary[2,,]), 
        cols = everything(), 
        names_to = "Age", 
        values_to = "Proportion"
        ) |>
    select(Proportion) |>
    unlist()

seine_sim_agecomp_df$Upper_50 <- pivot_longer(
        as.data.frame(seine_agecomp_summary[4,,]), 
        cols = everything(), 
        names_to = "Age", 
        values_to = "Proportion"
        ) |>
    select(Proportion) |>
    unlist()

seine_sim_agecomp_df$Upper_95 <- pivot_longer(
        as.data.frame(seine_agecomp_summary[5,,]), 
        cols = everything(), 
        names_to = "Age", 
        values_to = "Proportion"
        ) |>
    select(Proportion) |>
    unlist()

seine_sim_agecomp_df_sim <- seine_sim_agecomp_df |>
    filter(Year >= curr.year)

simulation_seine_agecomp_plot <- ggplot(seine_sim_agecomp_df_sim, aes(x = Age)) +
    geom_linerange(aes(ymin = Lower_95, ymax = Upper_95)) +
    geom_rect(aes(ymin = Lower_50, ymax = Upper_50, xmin = as.numeric(Age) - 0.4, xmax = as.numeric(Age) + 0.4, fill = fill_colors), color = "black", alpha = .5) +
    geom_segment(aes(y = Proportion, yend = Proportion, x = as.numeric(Age) - 0.4, xend = as.numeric(Age) + 0.4)) +
    geom_point(aes(x = jitter(as.numeric(Age)), y = seine_agecomp_em.1)) +
    geom_point(aes(x = jitter(as.numeric(Age)), y = seine_agecomp_em.2)) +
    geom_point(aes(x = jitter(as.numeric(Age)), y = seine_agecomp_em.3)) +
    facet_wrap(~ Year, dir = "v") +
    theme(legend.position = "none")


## spawn age comps
# the following is an array with third dimension being summary statistics
spawn_agecomp_summary <- spawn_agecomp_sim |>
    lapply(FUN = \(x) as.matrix(x[spawn_years_index,4:10])) |>
    simplify2array() |>
    apply(1:2, quantile, prob = quants)

dimnames(spawn_agecomp_summary)[[1]] <- quant_names
dimnames(spawn_agecomp_summary)[[2]] <- spawn_years
dimnames(spawn_agecomp_summary)[[3]] <- 3:9

spawn_sim_agecomp_df <- spawn_agecomp_em |> 
    sapply(FUN = \(x) select(pivot_longer(x, cols = !Year, names_to = "Age", values_to = "spawn_agecomp_em"), spawn_agecomp_em)) |>
    do.call(cbind, args = _) |>
    as.data.frame() |>
    cbind(
        data.frame(Year = rep(spawn_years, times = rep(length(3:9), length(spawn_years))),
                   Age = 3:9, 
                   fill_colors = generate.colors(nyr = length(spawn_years), color.options = color.options))
    )

colnames(spawn_sim_agecomp_df)[colnames(spawn_sim_agecomp_df) == "spawn_agecomp_em"] <- paste0("spawn_agecomp_em.", 1:num_mse_iters)

spawn_sim_agecomp_df$Proportion <- pivot_longer(
        as.data.frame(spawn_agecomp_summary[3,,]), 
        cols = everything(), 
        names_to = "Age", 
        values_to = "Proportion"
        ) |>
    select(Proportion) |>
    unlist()

spawn_sim_agecomp_df$Lower_95 <- pivot_longer(
        as.data.frame(spawn_agecomp_summary[1,,]), 
        cols = everything(), 
        names_to = "Age", 
        values_to = "Proportion"
        ) |>
    select(Proportion) |>
    unlist()

spawn_sim_agecomp_df$Lower_50 <- pivot_longer(
        as.data.frame(spawn_agecomp_summary[2,,]), 
        cols = everything(), 
        names_to = "Age", 
        values_to = "Proportion"
        ) |>
    select(Proportion) |>
    unlist()

spawn_sim_agecomp_df$Upper_50 <- pivot_longer(
        as.data.frame(spawn_agecomp_summary[4,,]), 
        cols = everything(), 
        names_to = "Age", 
        values_to = "Proportion"
        ) |>
    select(Proportion) |>
    unlist()

spawn_sim_agecomp_df$Upper_95 <- pivot_longer(
        as.data.frame(spawn_agecomp_summary[5,,]), 
        cols = everything(), 
        names_to = "Age", 
        values_to = "Proportion"
        ) |>
    select(Proportion) |>
    unlist()

spawn_sim_agecomp_df_sim <- spawn_sim_agecomp_df |>
    filter(Year >= curr.year)

simulation_spawn_agecomp_plot <- ggplot(spawn_sim_agecomp_df_sim, aes(x = Age)) +
    geom_linerange(aes(ymin = Lower_95, ymax = Upper_95)) +
    geom_rect(aes(ymin = Lower_50, ymax = Upper_50, xmin = as.numeric(Age) - 0.4, xmax = as.numeric(Age) + 0.4, fill = fill_colors), color = "black", alpha = .5) +
    geom_segment(aes(y = Proportion, yend = Proportion, x = as.numeric(Age) - 0.4, xend = as.numeric(Age) + 0.4)) +
    geom_point(aes(x = jitter(as.numeric(Age)), y = spawn_agecomp_em.1)) +
    geom_point(aes(x = jitter(as.numeric(Age)), y = spawn_agecomp_em.2)) +
    geom_point(aes(x = jitter(as.numeric(Age)), y = spawn_agecomp_em.3)) +
    facet_wrap(~ Year, dir = "v") +
    theme(legend.position = "none")


# ------------------------------------------------------------------------------

#### write outputs ####

if (write_outputs) {

    ## results
    write.csv(biomass_sim_df, here(dir_results, "simulated-biomass.csv"), row.names = FALSE)
    write.csv(mdm_sim_df, here(dir_results, "simulated-mdm.csv"), row.names = FALSE)
    write.csv(juv_sim_df, here(dir_results, "simulated-juv.csv"), row.names = FALSE)
    write.csv(seine_sim_agecomp_df, here(dir_results, "simulated-seine-agecomp.csv"), row.names = FALSE)
    write.csv(spawn_sim_agecomp_df, here(dir_results, "simulated-spawn-agecomp.csv"), row.names = FALSE)

    # performance metrics
    write.csv(performance_metrics, here(dir_results, "performance-metrics.csv"), row.names = FALSE)

    ## figures
    ggsave(here(dir_figures, "simulated-biomass.png"), simulation_biomass_plot, width = 6.5, height = 4)
    ggsave(here(dir_figures, "simulated-mdm.png"), simulation_mdm_plot, width = 6.5, height = 4)
    ggsave(here(dir_figures, "simulated-juv.png"), simulation_juv_plot, width = 6.5, height = 4)
    ggsave(here(dir_figures, "simulated-seine-agecomp.png"), simulation_seine_agecomp_plot, width = 10, height = 8.5)
    ggsave(here(dir_figures, "simulated-spawn-agecomp.png"), simulation_spawn_agecomp_plot, width = 12, height = 10)

}

# ------------------------------------------------------------------------------

#### end of file ####

end_time <- Sys.time()
mse_time <- end_time - start_time
cat("\n", paste("MSE runtime:", round(mse_time, 4), units(mse_time)), "\n")
