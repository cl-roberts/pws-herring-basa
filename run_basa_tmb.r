################################################################################

# Run BASA model and obtain parameter posteriors

# CL Roberts

# This script compiles the BASA TMB model 

################################################################################

#### front matter ####

## controls

chains <- 4
set.seed(8558)
seeds <- sample(1:1e4, size = chains)
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

## directory handling

dir_model <- here::here("model")

## TODO - calculate effective sample sizes

# IMPORTANT: insert code here to calculate effective sample sizes once model
# .cpp file is finished. This procedure is an iterative reweighting procedure 
# which calculates ESS's based on a ratio of observed and estimated age 
# compositions, I.e.
#
# estimate age comps --> calculate ESS's --> estimate age comps --> repeat until convergence
#
# Thus, model implementation must be finished in order estimate age comps for 
# calculating ESS's.
#
# For now, we will use the values in PWS_ASA(ESS).ctl calculated in the ADMB
# version of the model

PWS_ASA_ESS <- read.data.files(dir_model)$"PWS_ASA_ESS.ctl" 

## read data

# model data
PWS_ASA <- read.data.files(dir_model)$"PWS_ASA.dat"  
# recruitment and natural mortality deviate formulations
PWS_ASA_covariate <- read.data.files(dir_model)$"PWS_ASA_covariate.ctl"
# age composition sample sizes   
agecomp_samp_sizes <- read.data.files(dir_model)$"agecomp_samp_sizes.txt"
# disease data
PWS_ASA_disease <- read.data.files(dir_model)$"PWS_ASA_disease.dat"

model_data <- c(PWS_ASA_ESS, PWS_ASA, PWS_ASA_covariate, 
                agecomp_samp_sizes, PWS_ASA_disease)

## local vars

Y <- model_data$nyr
start.year <- 1980
curr.year <- start.year+Y

A <- model_data$nage


# ------------------------------------------------------------------------------ 

#### Compile Model and load DLL ####

fixed_pars <- list(pk = 0.75, egg_add = 0.4, Z_0_8 = 0.25, 
                    log_MeanAge0 = 6.2, 
                    sigma_age0devs = 0)

phase1_pars <- list(
                    annual_age0devs = rep(0, Y-1), 
                    log_juvenile_q = 4.22, 
                    loginit_pop = c(6.35,  5.66,  5.92,  6.74,  4.74)
                    )
phase2_pars <- list(Z_9 = 0.93, beta_mortality = rep(.2, 3), 
                    logmdm_c = 5.87, adfg_hydro_q = -0.38, pwssc_hydro_q = -0.21)
phase3_pars <- list(VHSV_age3_4_mort_93 = 0.08, ICH_age5_8_mort_93 = 0.22,
                    mat_age3 = 0.60, mat_age4 = 0.99)
phase4_pars <- list(seine_selex_alpha = 3.66, seine_selex_beta = 2.83)
phase5_pars <- list(milt_add_var = 0.33,
                    adfg_hydro_add_var = 0.30, pwssc_hydro_add_var = 0.32,
                    juvenile_overdispersion = 1.00)

parameters <- c(fixed_pars, phase1_pars, phase2_pars, phase3_pars, 
                phase4_pars, phase5_pars) 

if("PWS_ASA_tmb" %in% names(getLoadedDLLs())) {
    dyn.unload(dynlib(here::here(dir_model, "PWS_ASA_tmb")))
}
compile(here::here(dir_model, "PWS_ASA_tmb.cpp"))
dyn.load(dynlib(here::here(dir_model, "PWS_ASA_tmb")))

# ------------------------------------------------------------------------------

#### Fit and optimize model without phases ####

# fix some parameters
map <- list(pk = factor(NA), egg_add = factor(NA), Z_0_8 = factor(NA),       # fixed pars
            log_MeanAge0 = factor(NA), 
            sigma_age0devs = factor(NA)                                   
            # # ------------------------------------------------------------------ 
            # Z_9 = factor(NA), VHSV_age3_4_mort_93 = factor(NA),             
            # ICH_age5_8_mort_93 = factor(NA),                                # mort pars
            # beta_mortality = rep(factor(NA), 3)
            # # ------------------------------------------------------------------
            # mat_age3 = factor(NA), mat_age4 = factor(NA),                   # maturity pars
            # # ------------------------------------------------------------------
            # seine_selex_alpha = factor(NA), seine_selex_beta = factor(NA),  # select pars
            # # ------------------------------------------------------------------
            # logmdm_c = factor(NA), milt_add_var = factor(NA),
            # adfg_hydro_q = factor(NA), adfg_hydro_add_var = factor(NA),     # survey pars
            # pwssc_hydro_q = factor(NA), pwssc_hydro_add_var = factor(NA),
            # log_juvenile_q = factor(NA), 
            # juvenile_overdispersion = factor(NA),   
            # # ------------------------------------------------------------------
            # annual_age0devs = rep(factor(NA), Y-1)                           # recruit pars
            # # ------------------------------------------------------------------
            # loginit_pop = rep(factor(NA), 5)                                # init pop size (ages 1-5)
            ) 

# parameter bounds
lower <- c(annual_age0devs = rep(-5, Y-1), 
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
            juvenile_overdispersion = 0.01) 
upper <- c(annual_age0devs = rep(5, Y-1), 
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
            juvenile_overdispersion = 4)

est_parameters <- unlist(c(phase1_pars, phase2_pars, phase3_pars, 
                phase4_pars, phase5_pars)) 
                

# fit model
model <- MakeADFun(model_data, parameters, map = map, DLL = "PWS_ASA_tmb", hessian = TRUE, silent = TRUE)

# generate initial values for MCMC chains
# lower_init <- lower
# upper_init <- upper

# lower_init[grep("annual_age0devs", names(lower))] <- -0.2
# upper_init[grep("annual_age0devs", names(upper))] <- 0.2
# lower_init["VHSV_age3_4_mort_93"] <- .01
# upper_init["VHSV_age3_4_mort_93"] <- .5
# lower_init["ICH_age5_8_mort_93"] <- .01
# upper_init["ICH_age5_8_mort_93"] <- .5

lower_init <- c(annual_age0devs = rep(-.1, Y-1), 
				log_juvenile_q = -1,
				loginit_pop = rep(4.5, 5), 
				Z_9 = 0.4,  
				beta_mortality = rep(-5, 3), 
				logmdm_c = 3,  adfg_hydro_q = -1, pwssc_hydro_q = -1,     
				VHSV_age3_4_mort_93 = 0.15, ICH_age5_8_mort_93 = 0.1, 
				mat_age3 = 0.4, mat_age4 = 0.8,     
				seine_selex_alpha = 3, seine_selex_beta = 4,
				milt_add_var = 0.3, 
				adfg_hydro_add_var = 0.1, pwssc_hydro_add_var = 0.2,
				juvenile_overdispersion = 0.1) 
upper_init <- c(annual_age0devs = rep(.1, Y-1), 
				log_juvenile_q = 3,
				loginit_pop = rep(6, 5), 
				Z_9 = 1.1, 
				beta_mortality = rep(5, 3), 
				logmdm_c = 6, adfg_hydro_q = 1, pwssc_hydro_q = 1,     
				VHSV_age3_4_mort_93 = .25, ICH_age5_8_mort_93 = .4, 
				mat_age3 = 0.6, mat_age4 = .9, 
				seine_selex_alpha = 5, seine_selex_beta = 5,
				milt_add_var = 0.6,  
				adfg_hydro_add_var = 0.3, pwssc_hydro_add_var = 0.3,
				juvenile_overdispersion = 2)

inits <- init_tmb_params(chains = chains, lower = lower_init, upper = upper_init)

for (i in 1:length(inits)) {
	inits[[i]] <- est_parameters + rnorm(length(est_parameters), mean = 0, sd = .01)
}

# pre-MCMC model optimization 
# fit_ML <- nlminb(start = model$par, objective = model$fn, gradient = model$gr, 
#                  hessian = model$he, lower = lower, upper = upper,
#               control = list(eval.max = 10000, iter.max = 1000, rel.tol = 1e-10))

# mass_matrix <- solve(numDeriv::hessian(func = model$fn, x = fit_ML$par))


# run NUTS
mcmc_start_time <- Sys.time()
fit <- tmbstan(model, chains = chains, cores = chains, init = inits, iter = iter, control = control,
                warmup = warmup, seed = seeds[1], algorithm = "NUTS", silent = FALSE,
                lower = lower, upper = upper)
mcmc_end_time <- Sys.time()
time <- mcmc_end_time - mcmc_start_time

# save parameter posteriors
mcmc_results <- as.data.frame(fit) |> select(!lp__)
apply(mcmc_results, MARGIN = 2, FUN = median)


# ------------------------------------------------------------------------------ 

#### Model diagnostics ####

mon <- monitor(fit)
print(max(mon$Rhat))
print(min(mon$Tail_ESS))


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

Btilde_y <- vector(mode = "list", length = nrow(mcmc_results))
Btilde_post_y <- vector(mode = "list", length = nrow(mcmc_results))
age_3 <- vector(mode = "list", length = nrow(mcmc_results))
Btilde_forecast <- matrix(NA, nrow = n_iters, ncol = 1)

VARSreport <- matrix(NA, nrow = n_iters, ncol = 4)

likelihoods <- matrix(NA, nrow = nrow(mcmc_results), ncol = 8)

for(i in 1:nrow(mcmc_results)){

    other_posteriors <- model$report(mcmc_results[i,])

    MDM[i,] <- other_posteriors$That_y
    EGG[i,] <- other_posteriors$Ehat_y
    HYD_ADFG[i,] <- other_posteriors$Hhat_adfg_y
    HYD_PWSSC[i,] <- other_posteriors$Hhat_pwssc_y
    juv_schools[i,] <- other_posteriors$Jhat_y

    SpAC[i,] <- other_posteriors$spawn_age_comp_est |> t() |> c()
    SeAC[i,] <- other_posteriors$seine_age_comp_est |> t() |> c()
    N[i, ] <- other_posteriors$N_y_a |> t() |> c()

    Btilde_y[[i]] <- other_posteriors$Btilde_y
    Btilde_post_y[[i]] <- other_posteriors$Btilde_post_y
    age_3[[i]] <- other_posteriors$N_y_a[,4]
    Btilde_forecast[i,] <- other_posteriors$Btilde_forecast

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

}

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

## write diagnostics
