################################################################################

# Run BASA model and obtain parameter posteriors

# CL Roberts

# This script compiles the BASA TMB model, calculates effective sample sizes,
# and  

################################################################################

#### set up ####

start_time <- Sys.time()

## model controls

# declare which parameters to fix
fix <- c(
    "pk", "egg_add", "Z_0_8", "log_MeanAge0"
    # "sigma_age0devs"                                   
    # ------------------------------------------------------------------ 
    # "Z_9", 
    # "VHSV_age3_4_mort_93", "ICH_age5_8_mort_93" 
    # "beta_mortality",
    # ------------------------------------------------------------------
    # "mat_age3", "mat_age4",                  
    # ------------------------------------------------------------------
    # "seine_selex_alpha", "seine_selex_beta", 
    # ------------------------------------------------------------------
    # "logmdm_c",
    # "milt_add_var"
    # "adfg_hydro_q", "adfg_hydro_add_var", "pwssc_hydro_q", "pwssc_hydro_add_var",
    # "log_juvenile_q", "juvenile_overdispersion"
    # ------------------------------------------------------------------
    # "annual_age0devs"
    # ------------------------------------------------------------------
    # "loginit_pop"
) 

# MCMC controls
seed <- 406
set.seed(seed)
chains <- 4
iter <- 2000
warmup <- 700
control <- list(adapt_delta = 0.95)

# run retrospective analysis?
run_retro <- FALSE
n_peels <- 5
retro_iter <- 2000
retro_warmup <- 700
retro_chains <- 1

# populate app data?
app_data <- TRUE

# forecast controls
forecast_controls <- list(
    recruitment_average_years = 10,
    waa_average_years = 10, 
    disease_cov_average_years = 1, 
    expected_spring_harvest = 0,
    perc_female_forecast_years = 10    
)


## attach packages

library(TMB)
library(tmbstan)
library(rstan)
library(pwsHerringBasa)
library(dplyr)
library(here)

## directory handling

dir_model <- here("model")

dir_mcmc <- here(dir_model, "mcmc_out")
if (!dir.exists(dir_mcmc)) {
  dir.create(dir_mcmc)
}
dir_rep <- here(dir_model, "rep_out")
if (!dir.exists(dir_rep)) {
  dir.create(dir_rep)
}
dir_figures <- here("figures")
if (!dir.exists(dir_figures)) {
  dir.create(dir_figures)
}
dir_outputs <- here("data_outputs")
if (!dir.exists(dir_outputs)) {
  dir.create(dir_outputs)
}
if (run_retro) {
    dir_retro <- here("retrospectives")
    if (!dir.exists(dir_retro)) {
        dir.create(dir_retro)
    }
}
if (app_data) {
    dir_app_data <- here("mid-year-management-app", "data")
    if (!dir.exists(dir_app_data)) {
        dir.create(dir_app_data)
    }
}


## compile model
if("PWS_ASA" %in% names(getLoadedDLLs())) {
    dyn.unload(dynlib(here(dir_model, "PWS_ASA")))
}
compile(here(dir_model, "PWS_ASA.cpp"))
dyn.load(dynlib(here(dir_model, "PWS_ASA")))

# ------------------------------------------------------------------------------ 

#### data section ####

# model data
PWS_ASA <- read.data.files(dir_model)$"PWS_ASA.dat"  

# recruitment and natural mortality deviate formulations
PWS_ASA_covariate <- read.data.files(dir_model)$"PWS_ASA_covariate.ctl"

# disease data
PWS_ASA_disease <- read.data.files(dir_model)$"PWS_ASA_disease.dat"

# age composition sample sizes   
agecomp_samp_sizes <- read.data.files(dir_model)$"agecomp_samp_sizes.txt"

# effective sample size
# initially fit model using raw sample sizes
# iteratively compute ESS's below
PWS_ASA_ESS <- agecomp_samp_sizes

seine_ac <- PWS_ASA$seine_age_comp
seine_missing <- apply(seine_ac, MARGIN = 1, FUN = \(x) any(x == -9))

# identify missing data in raw sample sizes with -9's
spawn_ac <- PWS_ASA$spawn_age_comp
spawn_missing <- apply(spawn_ac, MARGIN = 1, FUN = \(x) any(x == -9))

PWS_ASA_ESS$seine_sample_size[seine_missing,] <- -9
PWS_ASA_ESS$spawn_sample_size[spawn_missing,] <- -9

names(PWS_ASA_ESS) <- c("seine_ess", "spawn_ess", "vhsv_ess", "ich_ess")

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

## sections of report files

llik <- c(paste0("L", 1:7), "penCount", "priors") 
derived <- c("Btilde_y", "Btilde_post_y", "N_y_a", "Ntilde_y_a", 
             "winter_survival", "summer_survival", "maturity")
survey <- c("seine_age_comp_est", "spawn_age_comp_est", "Ehat_y", 
            "That_y", "Hhat_adfg_y", "Hhat_pwssc_y", "Jhat_y")
forecast <- c("Btilde_forecast", "waa_forecast", "mean_log_rec", "N_a_forecast",
              "Btilde_agecomp_forecast", "Ntilde_agecomp_forecast", "mdm_forecast")

# ------------------------------------------------------------------------------ 

#### Parameter section ####

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

parameters <- c(
    fixed_pars, phase1_pars, phase2_pars, 
    phase3_pars, phase4_pars, phase5_pars
) 

# ------------------------------------------------------------------------------

#### set up other model attributes ####

# the map object is used to fix parameters
map <- list(
    pk = factor(NA), egg_add = factor(NA), Z_0_8 = factor(NA),       # fixed pars
    log_MeanAge0 = factor(NA),
    sigma_age0devs = factor(NA),                                   
    # ------------------------------------------------------------------ 
    Z_9 = factor(NA), VHSV_age3_4_mort_93 = factor(NA),             
    ICH_age5_8_mort_93 = factor(NA),                                # mort pars
    beta_mortality = rep(factor(NA), 3),
    # ------------------------------------------------------------------
    mat_age3 = factor(NA), mat_age4 = factor(NA),                   # maturity pars
    # ------------------------------------------------------------------
    seine_selex_alpha = factor(NA), seine_selex_beta = factor(NA),  # selectivity pars
    # ------------------------------------------------------------------
    logmdm_c = factor(NA), milt_add_var = factor(NA),
    adfg_hydro_q = factor(NA), adfg_hydro_add_var = factor(NA),     # survey pars
    pwssc_hydro_q = factor(NA), pwssc_hydro_add_var = factor(NA),
    log_juvenile_q = factor(NA), 
    juvenile_overdispersion = factor(NA),
    # ------------------------------------------------------------------
    annual_age0devs = rep(factor(NA), Y-1),                           # recruit deviates
    # ------------------------------------------------------------------
    loginit_pop = rep(factor(NA), 5)                                # init pop size (ages 1-5)
) 

map <- map[names(map) %in% fix]

# estimated parameter lower bounds
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
# lower <- lower[!grepl("beta_mortality", names(lower))]

# estimated parameter upper bounds
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
# upper <- upper[!grepl("beta_mortality", names(upper))]


# make model object
# model <- MakeADFun(model_data, parameters, map = map, DLL = "PWS_ASA", 
#                    silent = TRUE, random = "annual_age0devs")
# model <- MakeADFun(model_data, parameters, map = map, DLL = "PWS_ASA", 
#                    silent = TRUE)

# ------------------------------------------------------------------------------


#### calculate effective sample sizes ####

# calculate ESS's
ess <- calculate_ess(
    data = model_data, max_iter = 10, 
    parameters = parameters, map = map, DLL = "PWS_ASA",
    lower = lower, upper = upper, 
    control = list(eval.max = 10000, iter.max = 1000, rel.tol = 1e-10)
)

# write ESS's to file
write(
    c(ess$message, 
      "\n#seine_ess", ess$seine_ess, 
      "\n#spawn_ess", ess$spawn_ess),
    file = here(dir_model, "agecomp_effective_sample_size.txt"),
    sep = "\n"
)

# save ESS's to model data for fitting BASA
model_data$seine_ess <- ess$seine_ess
model_data$spawn_ess <- ess$spawn_ess


# ------------------------------------------------------------------------------

# compile model
# if("PWS_ASA" %in% names(getLoadedDLLs())) {
#     dyn.unload(dynlib(here(dir_model, "PWS_ASA")))
# }
# compile(here(dir_model, "PWS_ASA.cpp"))
# dyn.load(dynlib(here(dir_model, "PWS_ASA")))


#### create model object and write pre-optimization report ####

# fit model with iteratively re-weighted effective sample sizes 
model <- MakeADFun(
    model_data, parameters, map = map, DLL = "PWS_ASA", 
    silent = TRUE, hessian = TRUE
)

# write pre-optimization report
write_report(
    obj = model, par = model$par, 
    file = here(dir_rep, "pre-optim-report.txt"),
    # dummy = c("dummy", "dummy_vector"),
    llik = llik, derived = derived, survey = survey, forecast = forecast
)

# ------------------------------------------------------------------------------

#### fit model with ML ####

# fit with ML
fit_ML <- nlminb(
    start = model$par, objective = model$fn, gradient = model$gr, 
    lower = lower, upper = upper,
    control = list(eval.max = 10000, iter.max = 1000, rel.tol = 1e-10)
)

# write maximum likelihood report
write_report(
    obj = model, par = fit_ML$par, 
    file = here(dir_rep, "ml-report.txt"),
    llik = llik, derived = derived, survey = survey, forecast = forecast
)

# ------------------------------------------------------------------------------

#### fit model with markov chain monte carlo ####

# randomize initial values for MCMC chains
inits <- init_tmb_params(
    chains = chains, seed = seed,
    start = model$env$last.par.best, 
    lower = lower, upper = upper
)

# run NUTS
mcmc_start_time <- Sys.time()
fit <- tmbstan(
    model, chains = chains, lower = lower, upper = upper, 
    cores = chains, init = inits, iter = iter, 
    control = control, warmup = warmup, seed = seed, algorithm = "NUTS", 
    silent = FALSE
)
mcmc_end_time <- Sys.time()
mcmc_time <- mcmc_end_time - mcmc_start_time
message(paste("MCMC time:", round(mcmc_time, 4), units(mcmc_time)))


# save parameter posteriors
mcmc_results <- as.data.frame(fit) |> select(!lp__)

# parameter estimates are taken to be posterior medians
fit_mcmc_par <- apply(mcmc_results, MARGIN = 2, FUN = median)

# this bit is so that the write_report function does not insert a line break after 
# each element in vector parameters. to be improved another day
names(fit_mcmc_par)[grep("annual_age0devs", names(fit_mcmc_par))] <- "annual_age0devs"
names(fit_mcmc_par)[grep("loginit_pop", names(fit_mcmc_par))] <- "loginit_pop"
names(fit_mcmc_par)[grep("beta_mortality", names(fit_mcmc_par))] <- "beta_mortality"

# write mcmc report
write_report(
    obj = model, par = fit_mcmc_par, 
    file = here(dir_rep, "mcmc-report.txt"),
    llik = llik, derived = derived, survey = survey, forecast = forecast
)

# cool shiny app for mcmc chain diagnostics
# shinystan::launch_shinystan(fit)


# ------------------------------------------------------------------------------ 

#### write results ####

## set up objects for saving mcmc results ----

n_iters <- nrow(mcmc_results)

# survey quantities
SpAC <- matrix(NA, nrow = n_iters, ncol = A*Y)
SeAC <- matrix(NA, nrow = n_iters, ncol = A*Y)
HYD_ADFG <- matrix(NA, nrow = n_iters, ncol = Y)
HYD_PWSSC <- matrix(NA, nrow = n_iters, ncol = Y)
EGG <- matrix(NA, nrow = n_iters, ncol = Y)
MDM <- matrix(NA, nrow = n_iters, ncol = Y)
juv_schools <- matrix(NA, nrow = n_iters, ncol = Y)

# population dynamics
N <- matrix(NA, nrow = n_iters, ncol = A*Y)
Ntilde <- matrix(NA, nrow = n_iters, ncol = A*Y)
Btilde_y <- vector(mode = "list", length = n_iters)
Btilde_post_y <- vector(mode = "list", length = n_iters)
age_3 <- vector(mode = "list", length = n_iters)
summer_survival <- matrix(NA, nrow = n_iters, ncol = A*Y)
winter_survival <- matrix(NA, nrow = n_iters, ncol = A*Y)
VARSreport <- matrix(NA, nrow = n_iters, ncol = 4)

# forecast quantities
mean_log_rec <- matrix(NA, nrow = n_iters, ncol = 1)
N_a_forecast <- matrix(NA, nrow = n_iters, ncol = A)
Ntilde_a_forecast <- matrix(NA, nrow = n_iters, ncol = A)
Ntilde_agecomp_forecast <- matrix(NA, nrow = n_iters, ncol = A)
Btilde_forecast <- matrix(NA, nrow = n_iters, ncol = 1)
Btilde_a_forecast <- matrix(NA, nrow = n_iters, ncol = A)
Btilde_agecomp_forecast <- matrix(NA, nrow = n_iters, ncol = A)
winter_survival_forecast <- matrix(NA, nrow = n_iters, ncol = A)
mdm_forecast <- matrix(NA, nrow = n_iters, ncol = 1)

# objective function quantities
likelihoods <- matrix(NA, nrow = n_iters, ncol = 8)
penalties <- matrix(NA, nrow = n_iters, ncol = 5)
priors <- matrix(NA, nrow = n_iters, ncol = 1)

## extract raw mcmc results ----

for(i in 1:n_iters){

    other_posteriors <- model$report(mcmc_results[i,])

    # survey quantities
    SpAC[i,] <- other_posteriors$spawn_age_comp_est |> t() |> c()
    SeAC[i,] <- other_posteriors$seine_age_comp_est |> t() |> c()
    HYD_ADFG[i,] <- other_posteriors$Hhat_adfg_y
    HYD_PWSSC[i,] <- other_posteriors$Hhat_pwssc_y
    EGG[i,] <- other_posteriors$Ehat_y
    MDM[i,] <- other_posteriors$That_y
    juv_schools[i,] <- other_posteriors$Jhat_y

    # population dynamics
    N[i, ] <- other_posteriors$N_y_a |> t() |> c()
    Ntilde[i, ] <- other_posteriors$Ntilde_y_a |> t() |> c()
    Btilde_y[[i]] <- other_posteriors$Btilde_y
    Btilde_post_y[[i]] <- other_posteriors$Btilde_post_y
    age_3[[i]] <- other_posteriors$N_y_a[,4]
    summer_survival[i,] <- other_posteriors$summer_survival |> t() |> c()
    winter_survival[i,] <- other_posteriors$winter_survival |> t() |> c()
    VARSreport[i,] <- unlist(
        c(mcmc_results[i, c("milt_add_var")], fixed_pars["egg_add"],
          mcmc_results[i, c("adfg_hydro_add_var", "pwssc_hydro_add_var")])
    )

    # forecast quantities
    mean_log_rec[i,] <- other_posteriors$mean_log_rec
    N_a_forecast[i,] <- other_posteriors$N_a_forecast
    Ntilde_a_forecast[i,] <- other_posteriors$Ntilde_a_forecast
    Ntilde_agecomp_forecast[i,] <- other_posteriors$Ntilde_agecomp_forecast
    Btilde_forecast[i,] <- other_posteriors$Btilde_forecast
    Btilde_a_forecast[i,] <- other_posteriors$Btilde_a_forecast
    Btilde_agecomp_forecast[i,] <- other_posteriors$Btilde_agecomp_forecast
    winter_survival_forecast[i,] <- other_posteriors$winter_survival_forecast
    mdm_forecast[i,] <- other_posteriors$mdm_forecast

    # objective function quantities
    likelihoods[i,1] <- other_posteriors$L1
    likelihoods[i,2] <- other_posteriors$L2
    likelihoods[i,3] <- other_posteriors$L3
    likelihoods[i,4] <- other_posteriors$L4
    likelihoods[i,5] <- other_posteriors$L5
    likelihoods[i,6] <- other_posteriors$L6
    likelihoods[i,7] <- other_posteriors$L7
    likelihoods[i,8] <- other_posteriors$negLogLik
    penalties[i,] <- c(
        other_posteriors$penCount, other_posteriors$naa_pen, other_posteriors$ntilde_pen,
        other_posteriors$winter_surv_pen, other_posteriors$summer_surv_pen
    )
    priors[i,] <- other_posteriors$priors

}


## write results to .csv files ----

## estimated parameters

# all mcmc iterations 
write.csv(
    mcmc_results, here(dir_mcmc, "iterations.csv"), row.names = FALSE
)

## survey quantities

# estimated spawn age comps 
write.table(
    SpAC, sep = ",", here(dir_mcmc, "SpAC.csv"), 
    row.names = FALSE, col.names = FALSE
)

# estimated seine age comps
write.table(
    SeAC, sep = ",", here(dir_mcmc, "SeAC.csv"), 
    row.names = FALSE, col.names = FALSE
)

# adfg hydro fit 
write.table(
    HYD_ADFG, sep = ",", here(dir_mcmc, "HYD_ADFG.csv"), 
    row.names = FALSE, col.names = FALSE
)

# pwssc hydro fit 
write.table(
    HYD_PWSSC, sep = ",", here(dir_mcmc, "HYD_PWSSC.csv"), 
    row.names = FALSE, col.names = FALSE
)

# egg fit 
write.table(
    EGG, sep = ",", here(dir_mcmc, "EGG.csv"), 
    row.names = FALSE, col.names = FALSE
)

# mdm fit 
write.table(
    MDM, sep = ",", here(dir_mcmc, "MDM.csv"), 
    row.names = FALSE, col.names = FALSE
)

# juv aerial fit 
write.table(
    juv_schools, sep = ",", here(dir_mcmc, "juv_schools.csv"), 
    row.names = FALSE, col.names = FALSE
)

## population dynamics

# estimated numbers-at-age 
write.csv(
    N, here(dir_mcmc, "Num_at_age.csv"), 
    row.names = FALSE
)

# estimated pre fishery spawning biomass 
write.table(
    cbind(t(do.call(cbind, Btilde_y)), Btilde_forecast), 
    sep = ",", here(dir_mcmc, "PFRBiomass.csv"), 
    row.names = FALSE, col.names = FALSE
)

# estimated post fishery spawning biomass
write.csv(
    do.call(cbind, Btilde_post_y), 
    here(dir_mcmc, "post-fishery-spawning-biomass.csv"), 
    row.names = FALSE
)

# estimated age-3 herring 
write.table(
    t(do.call(cbind, age_3)), sep = ",", here(dir_mcmc, "Age3.csv"), 
    row.names = FALSE, col.names = FALSE
)   

# summer survival 
write.csv(
    summer_survival, 
    here(dir_mcmc, "adult_survival_effects_summer.csv"), 
    row.names = FALSE
)

# winter survival
write.csv(
    cbind(winter_survival, winter_survival_forecast), 
    here(dir_mcmc, "adult_survival_effects_winter.csv"), 
    row.names = FALSE
)

# estimated and fixed variances 
write.table(
    VARSreport, sep = ",", here(dir_mcmc, "VARSreport.csv"), 
    row.names = FALSE, col.names = FALSE
)

## forecast quantities

# numbers-at-age forecast
write.csv(
    N_a_forecast, here(dir_mcmc, "N_a_forecast.csv"), 
    row.names = FALSE
)

# mature numbers-at-age forecast
write.csv(
    Ntilde_a_forecast, here(dir_mcmc, "Ntilde_a_forecast.csv"), 
    row.names = FALSE
)

# biomass-at-age forecast
write.csv(
    Btilde_a_forecast, here(dir_mcmc, "Btilde_a_forecast.csv"), 
    row.names = FALSE
)


# biomass agecomp forecast
write.csv(
    Btilde_agecomp_forecast, here(dir_mcmc, "Btilde_agecomp_forecast.csv"), 
    row.names = FALSE
)

# mdm forecast
write.csv(
    mdm_forecast, here(dir_mcmc, "mdm_forecast.csv"), 
    row.names = FALSE
)


## objective function quantities

# likelihood components 
write.csv(
    likelihoods, here(dir_mcmc, "likelihoods.csv"), 
    row.names = FALSE
)

# penalties 
write.csv(
    penalties, here(dir_mcmc, "penalties.csv"), 
    row.names = FALSE
)


## estimates for mid-year management calculations 

mid_year_management <- data.frame(
    Btilde_forecast, Btilde_a_forecast, 
    mcmc_results$logmdm_c, mdm_forecast,
    N_a_forecast, mcmc_results$mat_age3, mcmc_results$mat_age4, 
    winter_survival_forecast, mean_log_rec,
    Ntilde_agecomp_forecast, Btilde_agecomp_forecast
)

colnames(mid_year_management) <- c(
    "Btilde_forecast", paste0("Btilde_age", 0:9, "_forecast"), 
    "logmdm_c", "mdm_forecast", 
    paste0("N_a_forecast_age", 0:9), "mat_age3", "mat_age4",
    paste0("survival_forecast_age", 0:9), 
    "mean_log_rec", 
    paste0("Ntilde_agecomp_forecast_age", 0:9),
    paste0("Btilde_agecomp_forecast_age", 0:9)
)

mid_year_management <- mid_year_management |>
    arrange(Btilde_forecast)

write.csv(
    mid_year_management, here(dir_mcmc, "mid-year-management.csv"), 
    row.names = FALSE
)

# move data files to app dir
if (app_data) {

    file.copy(
        from = here(dir_mcmc, "mid-year-management.csv"), 
        to = here(dir_app_data, "mid-year-management.csv")
    )
    file.copy(
        from = here(dir_model, "PWS_ASA.dat"), 
        to = here(dir_app_data, "PWS_ASA.dat")
    )

}

## print diagnostics ----

mon <- monitor(fit)
end_time <- Sys.time()
time <- end_time - start_time

message(paste("maximum Rhat:", max(mon$Rhat)))
message(paste("minimum bulk ESS:", min(mon$Bulk_ESS)))
message(paste("minimum tail ESS:", min(mon$Tail_ESS)))
check_treedepth(fit)
check_divergences(fit)
message(paste("Total model run time:", round(time, 4), units(time)))


## run retrospective analysis ----

if (run_retro) {

    retro_start_time <- Sys.time()

    for(i in 1:n_peels){

        dir_retro_i <- here(dir_retro, paste0("basa-", i))
        dir_mcmc_i <- here(dir_retro_i, "mcmc_out")
        if(!dir.exists(dir_retro_i)) dir.create(dir_retro_i)
        if(!dir.exists(dir_mcmc_i)) dir.create(dir_mcmc_i)
        
        # incrementally remove data to fit model
        retro_data <- lapply(model_data, FUN = \(x) {
            if (all(!is.null(dim(x)), nrow(x)==Y)) {
                x[1:(Y-i), , drop = FALSE]
            } else {
                x
            }
        })
        retro_data$nyr <- Y-i

        # fit model
        model_retro <- MakeADFun(
            retro_data, parameters, map = map, DLL = "PWS_ASA", 
            silent = TRUE, hessian = FALSE
        )

        # ML optimization
        retro_fit_ML <- nlminb(
            start = model_retro$par, objective = model_retro$fn, gradient = model_retro$gr, 
            lower = lower, upper = upper,
            control = list(eval.max = 10000, iter.max = 1000, rel.tol = 1e-10)
        )

        # run NUTS
        retro_inits <- init_tmb_params(
            chains = retro_chains, seed = seed,
            start = model_retro$env$last.par.best, 
            lower = lower, upper = upper
        )
        retro_fits <- vector(mode = "list", length = n_peels)
        retro_fits[[i]] <- tmbstan(
            model, chains = retro_chains, lower = lower, upper = upper, 
            cores = retro_chains, init = retro_inits, iter = retro_iter, 
            control = control, warmup = retro_warmup, seed = seed, algorithm = "NUTS", 
            silent = FALSE
        )

        # write results
        retro_mcmc_results <- as.data.frame(retro_fits[[i]]) |> 
            select(!lp__)
        
        n_retro_iter <- nrow(retro_mcmc_results)

        retro_Btilde_y <- vector(mode = "list", length = n_retro_iter)
        retro_Btilde_forecast <- matrix(NA, nrow = n_retro_iter, ncol = 1)
        for(iter in 1:n_retro_iter) {
            retro_Btilde_y[[iter]] <- model_retro$report(retro_mcmc_results[iter,])$Btilde_y 
            retro_Btilde_forecast[iter] <- model_retro$report(retro_mcmc_results[iter,])$Btilde_forecast
        }

        # estimated pre fishery spawning biomass 
        write.table(
            cbind(t(do.call(cbind, retro_Btilde_y)), retro_Btilde_forecast), 
            sep = ",", here(dir_mcmc_i, "PFRBiomass.csv"), 
            row.names = FALSE, col.names = FALSE
        )
        saveRDS(retro_fits[[i]], here(dir_retro_i, "fit.rds"))

    }

    retro_end_time <- Sys.time()
    retro_time <- retro_end_time - retro_start_time
    message(paste("Retrospective analysis run time:", round(retro_time, 4), units(retro_time)))

}

