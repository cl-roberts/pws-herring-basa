################################################################################

# Fit BASA iteratively omitting years of survey data

# CL Roberts
# 05/29/2025

# Summary:
# This script runs the Bayesian ASA model using the No-U-turn (NUTS) MCMC sampler 
# to obtain posteriors. 

# repeats model fitting procedure iteratively omitting each year in each of the 
# following surveys:
#
# Milt
# Juvenile schools
# Egg deposition
# Hydroacoustics
# Survey age comp 
#

################################################################################

#### front matter ####

# load packages

library(data.table)
library(adnuts)
library(snowfall)
library(rstan)
library(r4ss)
library(pwsHerringBasa)
library(here)

# set up directories
# each survey gets its own directory

survey <- "spawn_age_comp"

dir_model <- here("cost-of-info/model")
dir_survey <- here("cost-of-info", survey) 

if (!dir.exists(dir_survey)) {
    dir.create(dir_survey)
}

# mcmc options

n.samples <- 1000 
n.warmup <- 400 
n.time <- 6
n.chains <- 4

#-------------------------------------------------------------------------------

#### compile model, calculate ESS's execute NUTS for parameter posteriors ####
# see ?pwsHerringBasa::run.basa() for more details

# identify non-missing years in survey data set

dat_file <- read.data.files(dir_model)
all_years <- 1980:2024

if (survey != "base") {
    survey_years <- which(dat_file$PWS_ASA.dat[[survey]][,1] != -9)
} else {
    survey_years <- 1
    to_omit <- "base"
    dir_survey_year <- here("cost-of-info", "base") 
}

# loop through all non-missing years in survey dataset 

for (y in survey_years) {

    if (survey != "base") {

        # read data to start from full dataset
        dat_file_omitted <- read.data.files(dir_model)

        # delete year y
        dat_file_omitted$PWS_ASA.dat[[survey]][y,] <- -9

        if (survey == "egg") {
            dat_file_omitted$PWS_ASA.dat[["egg_se"]][y,] <- -9
        }

        if (survey == "pwssc_hydro") {
            dat_file_omitted$PWS_ASA.dat[["pwssc_hydro_se"]][y,] <- -9
        }

        # identify survey and year to omit
        to_omit <- paste(survey, all_years[y], sep = "_")

        # make directory for omitted survey and year
        dir_survey_year <- here("cost-of-info", survey, to_omit) 

        if (!dir.exists(dir_survey_year)) {
            dir.create(dir_survey_year)
        }

        # write dat file with omitted year
        dir_dat_file <- here(dir_model, "PWS_ASA.dat")
        write(paste("#", names(dat_file_omitted$PWS_ASA.dat[1])), file = dir_dat_file)
        write(dat_file$PWS_ASA.dat[[1]], file = dir_dat_file, append = TRUE)

        for (i in 2:length(dat_file_omitted$PWS_ASA.dat)) {

            write(paste0("\n#", " (", names(dat_file_omitted$PWS_ASA.dat[i]), ")"), 
                file = dir_dat_file, append = TRUE)
            to_write <- dat_file_omitted$PWS_ASA.dat[[i]]
            write(t(to_write), 
                ncolumns = ifelse(is.null(dim(to_write)), 1, ncol(to_write)),
                file = dir_dat_file, append = TRUE)

        }

    }

    # run model
    run <- run.basa(dir_model, n.samples = n.samples, n.warmup = n.warmup,
                    n.time = n.time, n.chains = n.chains)

    #### extract model run metadata and diagnostics ####

    # Extracts NUTS stats (energy, leapfrog transitions,etc)
    mon <- monitor(run$fit1$samples, warmup=run$fit1$warmup, print=FALSE)
    x <- extract_sampler_params(run$fit1)

    # Quick check for divergences & Gelman-Ruben statistic
    n.divergences <- sum(x$divergent__)/nrow(x)
    r.hat <- max(mon[, "Rhat"])<=1.1

    # summarize NUTS/MCMC diagnostics
    sum.dia <- data.frame(divergences.from.extract.function=sum(x$divergent__)/nrow(x),
                        min.ESS=min(mon[, "n_eff"]),
                        which.min.ESS=names(which.min(mon[, "n_eff"])),
                        max.Rhat=max(mon[, "Rhat"]),
                        which.max.Rhat=names(which.max(mon[, "Rhat"])),
                        time.elapsed=run$time)

    #### save information from model run ####

    # Write summary of parameter posteriors (medians, percentiles, etc)
    write.csv(as.data.frame(mon), 
            file = here(dir_survey_year, paste0("posterior_summary_", to_omit, ".csv")))

    # Write all MCMC samples of the parameters
    mcmc.samps <- data.frame(matrix(run$fit1$samples, ncol=dim(run$fit1$samples)[3], byrow=FALSE))
    names(mcmc.samps) <- run$fit1$par_names
    write.csv(mcmc.samps, 
            file = here(dir_survey_year, paste0("iterations_", to_omit, ".csv")),
            row.names=FALSE)

    # write summary file of NUTS/MCMC diagnostics
    write.table(sum.dia, 
                file = here(dir_survey_year, paste0("convergence_diagnostics_", to_omit, ".csv")), 
                sep=",", row.names=FALSE)
    saveRDS(run$fit1, file = here(dir_survey_year, paste0("NUTS_fit_", to_omit, ".RDS")))

    # copy biomass posterior samples to dir_survey_year
    file.copy(from = here(dir_model, "mcmc_out/PFRBiomass.csv"),
            to = here(dir_survey_year, paste0("PFRBiomass_", to_omit, ".csv")),
            overwrite = TRUE)

    # copy back unaltered data file
    file.copy(from = here("model/PWS_ASA.dat"),
            to = here(dir_model, "PWS_ASA.dat"),
            overwrite = TRUE)

}

