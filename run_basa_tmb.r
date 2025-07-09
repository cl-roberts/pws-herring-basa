################################################################################

# Run BASA model and obtain parameter posteriors

# CL Roberts

# This script compiles the BASA TMB model 

################################################################################

#### front matter ####

## MCMC controls

seed <- 406
set.seed(seed)
chains <- 1
iter <- 2000
warmup <- 700
control <- list(adapt_delta = 0.95)
# control <- list(adapt_delta = 0.9, max_treedepth = 12)

## attach packages

library(TMB)
library(tmbstan)
library(rstan)
library(pwsHerringBasa)
library(dplyr)
library(here)

## directory handling

dir_model <- here("model")

# compile model
if("PWS_ASA_tmb" %in% names(getLoadedDLLs())) {
    dyn.unload(dynlib(here(dir_model, "PWS_ASA_tmb")))
}
compile(here(dir_model, "PWS_ASA_tmb.cpp"))
dyn.load(dynlib(here(dir_model, "PWS_ASA_tmb")))

## read data

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

# forecast controls
forecast_controls <- list(
    recruitment_average_years = 10,
    waa_average_years = 5
)

model_data <- c(
    PWS_ASA, PWS_ASA_covariate, PWS_ASA_disease, 
    agecomp_samp_sizes, PWS_ASA_ESS,
    forecast_controls
)

## local vars

Y <- model_data$nyr
start.year <- 1980
curr.year <- start.year+Y

A <- model_data$nage


# ------------------------------------------------------------------------------ 

#### Compile Model and load DLL ####

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
    Z_9 = 0.93, beta_mortality = rep(.2, 3), 
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

#### Create model object ####

# fix some parameters
map <- list(
    pk = factor(NA), egg_add = factor(NA), Z_0_8 = factor(NA),       # fixed pars
    log_MeanAge0 = factor(NA)
    # sigma_age0devs = factor(NA)                                   
    # # ------------------------------------------------------------------ 
    # Z_9 = factor(NA), VHSV_age3_4_mort_93 = factor(NA),             
    # ICH_age5_8_mort_93 = factor(NA),                                # mort pars
    # beta_mortality = rep(factor(NA), 3)
    # # ------------------------------------------------------------------
    # mat_age3 = factor(NA), mat_age4 = factor(NA),                   # maturity pars
    # # ------------------------------------------------------------------
    # seine_selex_alpha = factor(NA), seine_selex_beta = factor(NA),  # selectivity pars
    # # ------------------------------------------------------------------
    # logmdm_c = factor(NA), milt_add_var = factor(NA),
    # adfg_hydro_q = factor(NA), adfg_hydro_add_var = factor(NA),     # survey pars
    # pwssc_hydro_q = factor(NA), pwssc_hydro_add_var = factor(NA),
    # log_juvenile_q = factor(NA), 
    # juvenile_overdispersion = factor(NA),   
    # # ------------------------------------------------------------------
    # annual_age0devs = rep(factor(NA), Y-1)                           # recruit deviates
    # # ------------------------------------------------------------------
    # loginit_pop = rep(factor(NA), 5)                                # init pop size (ages 1-5)
) 

# parameter bounds
lower <- c(
    # log_MeanAge0 = 2, 
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
upper <- c(
    # log_MeanAge0 = 20, 
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

est_parameters <- unlist(
    c(phase1_pars, phase2_pars, phase3_pars, phase4_pars, phase5_pars)
) 
                

# make model object
# model <- MakeADFun(model_data, parameters, map = map, DLL = "PWS_ASA_tmb", 
#                    silent = TRUE, random = "annual_age0devs")
# model <- MakeADFun(model_data, parameters, map = map, DLL = "PWS_ASA_tmb", 
#                    silent = FALSE)

# before.optim <- c(
#     last.par = model$env$last.par,
#     par = model$par,
#     fn = model$fn(),
#     gr = model$gr()
# )


# ------------------------------------------------------------------------------


#### calculate effective sample sizes ####

# estimate age comps --> calculate ESS's --> estimate age comps --> repeat until convergence

convergence <- FALSE 

seine_samp_size <- model_data$seine_sample_size
spawn_samp_size <- model_data$spawn_sample_size

its <- 0
max_iter <- 10

while (!convergence) {

    its <- its+1

    # stop while loop if ESS's haven't converged
    if (its > max_iter) {
        stop(
            paste("ESS calculation hasn't converged after", max_iter, "iterations")
        )
    }

    # fit model using agecomp sample sizes of current iteration 
    model_iteration <- MakeADFun(
        model_data, parameters, map = map, DLL = "PWS_ASA_tmb", 
        silent = TRUE, hessian = FALSE
    )

    # optimize model 
    fit_iteration <- nlminb(
        start = model_iteration$par, 
        objective = model_iteration$fn, 
        gradient = model_iteration$gr, 
        lower = lower, upper = upper, 
        control = list(eval.max = 10000, iter.max = 1000, rel.tol = 1e-10)
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

    model_data$seine_ess <- seine_ess
    model_data$spawn_ess <- spawn_ess

}

exit_message <- paste("ESS calculations converged after", its, "iterations")

print(exit_message)

write(
    c(exit_message, 
      "\n#seine_ess", model_data$seine_ess, 
      "\n#spawn_ess", model_data$spawn_ess),
    file = here(dir_model, "agecomp_effective_sample_size.txt"),
    sep = "\n"
)

# ------------------------------------------------------------------------------

#### fit model with ML using final effective sample sizes ####

# fit model using agecomp sample sizes of current iteration 
model <- MakeADFun(
    model_data, parameters, map = map, DLL = "PWS_ASA_tmb", 
    silent = TRUE, hessian = TRUE
)

# fit with ML
fit_ML <- nlminb(
    start = model$par, objective = model$fn, gradient = model$gr, 
    lower = lower, upper = upper,
    control = list(eval.max = 10000, iter.max = 1000, rel.tol = 1e-10)
)

sdreport(model)
rep <- model$report(fit_ML$par)


# calculate initial values for mcmc chains
init_tmb_params <- function(chains, start, lower, upper, seed){

    set.seed(seed)

    par_names <- names(lower) 
    sd <- abs(.75*start)

    inits <- list()
    
    for(j in 1:chains){

        inits[[j]] <- rnorm(length(start), start, sd = sd)

        while (any(inits[[j]] < lower) | any(inits[[j]] > upper)) {
            inits[[j]] <- rnorm(length(start), start, sd = sd)
        }
    
        names(inits[[j]]) <- par_names
    }

    return(inits)
}

# randomize initial values for all chains
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
time <- mcmc_end_time - mcmc_start_time


# save parameter posteriors
mcmc_results <- as.data.frame(fit) |> select(!lp__)
apply(mcmc_results, MARGIN = 2, FUN = median)

# cool shiny app for mcmc chain diagnostics
# shinystan::launch_shinystan(fit)


# ------------------------------------------------------------------------------ 

#### write results ####

## extract raw mcmc results

n_iters <- nrow(mcmc_results)

MDM <- matrix(NA, nrow = n_iters, ncol = Y)
EGG <- matrix(NA, nrow = n_iters, ncol = Y)
HYD_ADFG <- matrix(NA, nrow = n_iters, ncol = Y)
HYD_PWSSC <- matrix(NA, nrow = n_iters, ncol = Y)
juv_schools <- matrix(NA, nrow = n_iters, ncol = Y)

SpAC <- matrix(NA, nrow = n_iters, ncol = A*Y)
SeAC <- matrix(NA, nrow = n_iters, ncol = A*Y)

N <- matrix(NA, nrow = n_iters, ncol = A*Y)
Ntilde <- matrix(NA, nrow = n_iters, ncol = A*Y)

Btilde_y <- vector(mode = "list", length = n_iters)
Btilde_post_y <- vector(mode = "list", length = n_iters)
age_3 <- vector(mode = "list", length = n_iters)

mean_log_rec <- matrix(NA, nrow = n_iters, ncol = 1)
Btilde_forecast <- matrix(NA, nrow = n_iters, ncol = 1)

winter_survival <- vector(mode = "list", length = n_iters)
summer_survival <- vector(mode = "list", length = n_iters)

VARSreport <- matrix(NA, nrow = n_iters, ncol = 4)

likelihoods <- matrix(NA, nrow = n_iters, ncol = 8)

penalties <- matrix(NA, nrow = n_iters, ncol = 5)

priors <- matrix(NA, nrow = n_iters, ncol = 1)

for(i in 1:n_iters){

    other_posteriors <- model$report(mcmc_results[i,])

    MDM[i,] <- other_posteriors$That_y
    EGG[i,] <- other_posteriors$Ehat_y
    HYD_ADFG[i,] <- other_posteriors$Hhat_adfg_y
    HYD_PWSSC[i,] <- other_posteriors$Hhat_pwssc_y
    juv_schools[i,] <- other_posteriors$Jhat_y

    SpAC[i,] <- other_posteriors$spawn_age_comp_est |> t() |> c()
    SeAC[i,] <- other_posteriors$seine_age_comp_est |> t() |> c()

    N[i, ] <- other_posteriors$N_y_a |> t() |> c()
    Ntilde[i, ] <- other_posteriors$Ntilde_y_a |> t() |> c()

    Btilde_y[[i]] <- other_posteriors$Btilde_y
    Btilde_post_y[[i]] <- other_posteriors$Btilde_post_y
    age_3[[i]] <- other_posteriors$N_y_a[,4]

    mean_log_rec[i,] <- other_posteriors$mean_log_rec
    Btilde_forecast[i,] <- other_posteriors$Btilde_forecast

    winter_survival[[i]] <- other_posteriors$winter_survival
    summer_survival[[i]] <- other_posteriors$summer_survival

    VARSreport[i,] <- unlist(c(mcmc_results[i, c("milt_add_var")],
                        fixed_pars["egg_add"],
                        mcmc_results[i, c("adfg_hydro_add_var", "pwssc_hydro_add_var")]
                        ))

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

# ---------------------------

dir_mcmc_out <- here("model/mcmc_out")

# investigate likelihoods

# tmb_lllk <- likelihoods
# admb_llik <- read.csv(here(dir_mcmc_out, "llikcomponents.csv"))

# # 3, 4, 5, 6, 7, 8/10
# llik <- data.frame(admb = admb_llik[,3], tmb = tmb_llik[-1,3])

# ggplot(llik) +
#     geom_histogram(aes(x = admb, fill = "admb"), alpha = .25) +
#     geom_histogram(aes(x = tmb, fill = "tmb"), alpha = .25) +
#     geom_vline(aes(xintercept = median(admb)), color = "firebrick") +
#     geom_vline(aes(xintercept = median(tmb)), color = "#404080") +
#     scale_fill_manual(values=c("firebrick", "#404080")) +
#     xlab("Total Negative Log-Likelihood")

# ---------------------------

# investigate priors


# tmb_priors <- priors 
# admb_priors_tmp <- read.csv(here(dir_mcmc_out, "priordensities.csv"), header = FALSE) |>
#     as.matrix() |>
#     t() |>
#     na.omit()
# admb_priors <- matrix(admb_priors_tmp, ncol = 45, byrow = TRUE) |>
#     apply(MARGIN = 1, FUN = sum)

# compare_priors <- data.frame(admb = admb_priors, tmb = tmb_priors)

# ggplot(compare_priors) +
#     geom_histogram(aes(x = admb, fill = "admb"), alpha = .25) +
#     geom_histogram(aes(x = tmb, fill = "tmb"), alpha = .25) +
#     geom_vline(aes(xintercept = median(admb)), color = "firebrick") +
#     geom_vline(aes(xintercept = median(tmb)), color = "#404080") +
#     scale_fill_manual(values=c("firebrick", "#404080")) +
#     xlab("Prior densities")

# ---------------------------

# check to see if penalties worked

sapply(age_3, FUN = median)

sapply(winter_survival, \(x) all((x>0) & (x<1))) |>
    all()

sapply(summer_survival, \(x) all((x>0) & (x<1))) |>
    all()

sapply(SpAC, \(x) all((x>=0) & (x<=1))) |>
    all()

sapply(SeAC, \(x) all((x>=0) & (x<=1))) |>
    all()

sapply(N, \(x) all(x >= 0)) |>
    all()

sapply(Ntilde, \(x) all(x >= 0)) |>
    all()

sapply(Btilde_y, \(x) all(x >= 0)) |>
    all()

sapply(Btilde_post_y, \(x) all(x >= 0)) |>
    all()

# pens_admb <- read.csv(here(dir_mcmc_out, "penalties.csv"))[,1:3]
# apply(pens_admb, MARGIN = 2, sum)
apply(penalties, MARGIN = 2, sum)

# age3_admb <- read.csv(here(dir_mcmc_out, "Age3.csv"))

# dim(age3_admb)


# ---------------------------


## create directories for tmb outputs

dir_mcmc_tmb <- here::here(dir_model, "mcmc_out_tmb")
if(!dir.exists(dir_mcmc_tmb)) {
  dir.create(dir_mcmc_tmb)
}
dir_figures <- here::here("figures/tmb")
if(!dir.exists(dir_figures)) {
  dir.create(dir_figures)
}
dir_outputs <- here::here("data_outputs/tmb")
if(!dir.exists(dir_outputs)) {
  dir.create(dir_outputs)
}

## write results to .csv files

# mdm fit ----
write.table(MDM, sep = ",",
          here::here(dir_mcmc_tmb, "MDM.csv"), 
          row.names = FALSE, col.names = FALSE)

# egg fit ----
write.table(EGG, sep = ",",
          here::here(dir_mcmc_tmb, "EGG.csv"), 
          row.names = FALSE, col.names = FALSE)

# adfg hydro fit ----
write.table(HYD_ADFG, sep = ",",
          here::here(dir_mcmc_tmb, "HYD_ADFG.csv"), 
          row.names = FALSE, col.names = FALSE)

# pwssc hydro fit ----
write.table(HYD_PWSSC, sep = ",",
          here::here(dir_mcmc_tmb, "HYD_PWSSC.csv"), 
          row.names = FALSE, col.names = FALSE)

# juv aerial fit ----
write.table(juv_schools, sep = ",",
          here::here(dir_mcmc_tmb, "juv_schools.csv"), 
          row.names = FALSE, col.names = FALSE)

# estimated spawn age comps ----
write.table(SpAC, sep = ",",
          here::here(dir_mcmc_tmb, "SpAC.csv"), 
          row.names = FALSE, col.names = FALSE)

# estimated seine age comps ----
write.table(SeAC, sep = ",", 
          here::here(dir_mcmc_tmb, "SeAC.csv"), 
          row.names = FALSE, col.names = FALSE)

# estimated numbers-at-age ----
write.table(N, sep = ",", 
          here::here(dir_mcmc_tmb, "Num_at_age.csv"), 
          row.names = FALSE, col.names = FALSE)

# estimated age-3 herring ----
write.table(t(do.call(cbind, age_3)), sep = ",", 
          here::here(dir_mcmc_tmb, "Age3.csv"), 
          row.names = FALSE, col.names = FALSE)

# estimated pre fishery spawning biomass ----
write.table(cbind(t(do.call(cbind, Btilde_y)), Btilde_forecast), sep = ",", 
          here::here(dir_mcmc_tmb, "PFRBiomass.csv"), 
          row.names = FALSE, col.names = FALSE)

# estimated post fishery spawning biomass ----
write.csv(do.call(cbind, Btilde_post_y), 
          here::here(dir_mcmc_tmb, "post-fishery-spawning-biomass.csv"), 
          row.names = FALSE)

# estimated and fixed variances ----
write.table(VARSreport, sep = ",", 
          here::here(dir_mcmc_tmb, "VARSreport.csv"), 
          row.names = FALSE, col.names = FALSE)

# likelihood components ----
write.csv(likelihoods, 
          here::here(dir_mcmc_tmb, "likelihoods.csv"), 
          row.names = FALSE)

# iterations ----
write.csv(mcmc_results, 
          here::here(dir_mcmc_tmb, "iterations.csv"), 
          row.names = FALSE)

## print diagnostics

mon <- monitor(fit)
print(max(mon$Rhat))
print(min(mon$Bulk_ESS))
print(min(mon$Tail_ESS))
check_treedepth(fit)
check_divergences(fit)


