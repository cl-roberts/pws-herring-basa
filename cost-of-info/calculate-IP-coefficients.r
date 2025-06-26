#------------------------------------------------------------------------------#

# calculate value of information coefficients for integer programming model

# CL Roberts

# for SEFS 540

# 5/30/2025

#------------------------------------------------------------------------------#


## set up ----

# pkgs
library(here)
library(ggplot2)
theme_set(theme_bw())

# dir 
dir_coi <- here("cost-of-info")

# functions

mpe <- function(x, y) {

    if (!all(c(is.numeric(x), is.numeric(y)))) {
        stop("x and y must be numeric vectors")
    }
    if (length(x) != length(y)) {
        stop("x and y must be vectors of the same length")
    }

    n <- length(x)
    out <- (100 / n) * sum(((x - y) / x))    

    return(out)
}

rmse <- function(x, y) {

    if (!all(c(is.numeric(x), is.numeric(y)))) {
        stop("x and y must be numeric vectors")
    }
    if (length(x) != length(y)) {
        stop("x and y must be vectors of the same length")
    }

    n <- length(x)
    out <- sqrt(sum((x-y)^2)/n)    

    return(out)
}

read_pfrb <- function(file, years) {

    pfrb <- read.csv(file = file) |>
        apply(MARGIN = 2, FUN = median)
    out <- data.frame(year = years, pfrb = pfrb)

    return(out)
}

# local vars

all_years <- 1980:2025

# base model run

base_pfrb <- read_pfrb(here(dir_coi, "base", "PFRBiomass_base.csv"), all_years)
base_par <- read.csv(here(dir_coi, "base", "posterior_summary_base.csv"))[-68,]


## calculate coefficients and plot surveys ----

## mdm

# read data

dirs_mdm <- list.files(here(dir_coi, "mdm"))

n_mdm <- length(dirs_mdm)

mdm_list_pfrb <- mdm_list_par <- vector(mode = "list", length = n_mdm)

for (i in 1:n_mdm) {
    mdm_year <- dirs_mdm[i]
    year <- unlist(stringr::str_split(mdm_year, "_"))[2] |>
        as.integer()
    dir_file <- here(dir_coi, "mdm", mdm_year, paste0("PFRBiomass_mdm_", year, ".csv"))
    mdm_list_pfrb[[i]] <- read_pfrb(dir_file, all_years)
    mdm_list_pfrb[[i]]$survey_year <- mdm_year
    dir_file <- here(dir_coi, "mdm", mdm_year, paste0("posterior_summary_mdm_", year, ".csv"))
    mdm_list_par[[i]] <- read.csv(dir_file)[-68,]
}

mdm_pfrb <- do.call(rbind, mdm_list_pfrb)
mdm_par <- do.call(rbind, mdm_list_par)


# calculate rmse
mdm_pfrb_mpe <- mpe(mdm_pfrb$pfrb, rep(base_pfrb$pfrb, n_mdm))
mdm_par_rmse <- rmse(mdm_par$mean, rep(base_par$mean, n_mdm))


# plot PFRB with missing data against model with all data

ggplot() +
    geom_line(data = base_pfrb, aes(x = year, y = pfrb)) +
    geom_line(data = mdm_pfrb, aes(x = year, y = pfrb, color = survey_year)) +
    theme(legend.position = "none")

## egg

# read data

dirs_egg <- list.files(here(dir_coi, "egg"))

n_egg <- length(dirs_egg)

egg_list_pfrb <- egg_list_par <- vector(mode = "list", length = n_egg)

for (i in 1:n_egg) {
    egg_year <- dirs_egg[i]
    year <- unlist(stringr::str_split(egg_year, "_"))[2] |>
        as.integer()
    dir_file <- here(dir_coi, "egg", egg_year, paste0("PFRBiomass_egg_", year, ".csv"))
    egg_list_pfrb[[i]] <- read_pfrb(dir_file, all_years)
    egg_list_pfrb[[i]]$survey_year <- egg_year
    dir_file <- here(dir_coi, "egg", egg_year, paste0("posterior_summary_egg_", year, ".csv"))
    egg_list_par[[i]] <- read.csv(dir_file)[-68,]
}

egg_pfrb <- do.call(rbind, egg_list_pfrb)
egg_par <- do.call(rbind, egg_list_par)


# calculate rmse
egg_pfrb_mpe <- mpe(egg_pfrb$pfrb, rep(base_pfrb$pfrb, n_egg))
egg_par_rmse <- rmse(egg_par$mean, rep(base_par$mean, n_egg))


# plot PFRB with missing data against model with all data

ggplot() +
    geom_line(data = egg_pfrb, aes(x = year, y = pfrb, color = survey_year)) +
    geom_line(data = base_pfrb, aes(x = year, y = pfrb)) +
    theme(legend.position = "none")


## juvenile_survey

# read data

dirs_juv <- list.files(here(dir_coi, "juvenile_survey"))

n_juv <- length(dirs_juv)

juv_list_pfrb <- juv_list_par <- vector(mode = "list", length = n_juv)

for (i in 1:n_juv) {
    juv_year <- dirs_juv[i]
    year <- unlist(stringr::str_split(juv_year, "_"))[3] |>
        as.integer()
    dir_file <- here(dir_coi, "juvenile_survey", juv_year, paste0("PFRBiomass_juvenile_survey_", year, ".csv"))
    juv_list_pfrb[[i]] <- read_pfrb(dir_file, all_years)
    juv_list_pfrb[[i]]$survey_year <- juv_year
    dir_file <- here(dir_coi, "juvenile_survey", juv_year, paste0("posterior_summary_juvenile_survey_", year, ".csv"))
    if (i < 13) {
        juv_list_par[[i]] <- read.csv(dir_file)[-68,]
    }
    if (i >= 13) {
        juv_list_par[[i]] <- read.csv(dir_file)[-66,]
    }
    print(nrow(juv_list_par[[i]]))
}

juv_pfrb <- do.call(rbind, juv_list_pfrb)
juv_par <- do.call(rbind, juv_list_par)


# calculate rmse
juv_pfrb_mpe <- mpe(juv_pfrb$pfrb, rep(base_pfrb$pfrb, n_juv))
juv_par_rmse <- rmse(juv_par$mean, c(rep(base_par$mean, n_juv-2), rep(base_par$mean[-63:-64], 2)))


# plot PFRB with missing data against model with all data

ggplot() +
    geom_line(data = juv_pfrb, aes(x = year, y = pfrb, color = survey_year)) +
    geom_line(data = base_pfrb, aes(x = year, y = pfrb)) +
    theme(legend.position = "none")



## pwssc_hydro

# read data

dirs_hydro <- list.files(here(dir_coi, "pwssc_hydro"))

n_hydro <- length(dirs_hydro)

hydro_list_pfrb <- hydro_list_par <- vector(mode = "list", length = n_hydro)

for (i in 1:n_hydro) {
    hydro_year <- dirs_hydro[i]
    year <- unlist(stringr::str_split(hydro_year, "_"))[3] |>
        as.integer()
    dir_file <- here(dir_coi, "pwssc_hydro", hydro_year, paste0("PFRBiomass_pwssc_hydro_", year, ".csv"))
    hydro_list_pfrb[[i]] <- read_pfrb(dir_file, all_years)
    hydro_list_pfrb[[i]]$survey_year <- hydro_year
    dir_file <- here(dir_coi, "pwssc_hydro", hydro_year, paste0("posterior_summary_pwssc_hydro_", year, ".csv"))
    hydro_list_par[[i]] <- read.csv(dir_file)[-68,]
}

hydro_pfrb <- do.call(rbind, hydro_list_pfrb)
hydro_par <- do.call(rbind, hydro_list_par)


# calculate rmse
hydro_pfrb_mpe <- mpe(hydro_pfrb$pfrb, rep(base_pfrb$pfrb, n_hydro))
hydro_par_rmse <- rmse(hydro_par$mean, rep(base_par$mean, n_hydro))


# plot PFRB with missing data against model with all data

ggplot() +
    geom_line(data = hydro_pfrb, aes(x = year, y = pfrb, color = survey_year)) +
    geom_line(data = base_pfrb, aes(x = year, y = pfrb)) +
    theme(legend.position = "none")


## disease_covs

# read data

dirs_disease <- list.files(here(dir_coi, "disease_covs"))

n_disease <- length(dirs_disease)

disease_list_pfrb <- disease_list_par <- vector(mode = "list", length = n_disease)

for (i in 1:n_disease) {
    disease_year <- dirs_disease[i]
    year <- unlist(stringr::str_split(disease_year, "_"))[3] |>
        as.integer()
    dir_file <- here(dir_coi, "disease_covs", disease_year, paste0("PFRBiomass_disease_covs_", year, ".csv"))
    disease_list_pfrb[[i]] <- read_pfrb(dir_file, all_years)
    disease_list_pfrb[[i]]$survey_year <- disease_year
    dir_file <- here(dir_coi, "disease_covs", disease_year, paste0("posterior_summary_disease_covs_", year, ".csv"))
    disease_list_par[[i]] <- read.csv(dir_file)[-68,]
}

disease_pfrb <- do.call(rbind, disease_list_pfrb)
disease_par <- do.call(rbind, disease_list_par)


# calculate rmse
disease_pfrb_mpe <- mpe(disease_pfrb$pfrb, rep(base_pfrb$pfrb, n_disease))
disease_par_rmse <- rmse(disease_par$mean, rep(base_par$mean, n_disease))


# plot PFRB with missing data against model with all data

ggplot() +
    geom_line(data = disease_pfrb, aes(x = year, y = pfrb, color = survey_year)) +
    geom_line(data = base_pfrb, aes(x = year, y = pfrb)) +
    theme(legend.position = "none")

## spawn_age_comp

# read data

dirs_agecomp <- list.files(here(dir_coi, "spawn_age_comp"))

n_agecomp <- length(dirs_agecomp)

agecomp_list_pfrb <- agecomp_list_par <- vector(mode = "list", length = n_agecomp)

for (i in 1:n_agecomp) {
    agecomp_year <- dirs_agecomp[i]
    year <- unlist(stringr::str_split(agecomp_year, "_"))[4] |>
        as.integer()
    dir_file <- here(dir_coi, "spawn_age_comp", agecomp_year, 
                     paste0("PFRBiomass_spawn_age_comp_", year, ".csv"))
    agecomp_list_pfrb[[i]] <- read_pfrb(dir_file, all_years)
    agecomp_list_pfrb[[i]]$survey_year <- agecomp_year
    dir_file <- here(dir_coi, "spawn_age_comp", agecomp_year, 
                     paste0("posterior_summary_spawn_age_comp_", year, ".csv"))
    agecomp_list_par[[i]] <- read.csv(dir_file)[-68,]
}

agecomp_pfrb <- do.call(rbind, agecomp_list_pfrb)
agecomp_par <- do.call(rbind, agecomp_list_par)


# calculate rmse
agecomp_pfrb_mpe <- mpe(agecomp_pfrb$pfrb, rep(base_pfrb$pfrb, n_agecomp))
agecomp_par_rmse <- rmse(agecomp_par$mean, rep(base_par$mean, n_agecomp))


# plot PFRB with missing data against model with all data

ggplot() +
    geom_line(data = agecomp_pfrb, aes(x = year, y = pfrb, color = survey_year)) +
    geom_line(data = base_pfrb, aes(x = year, y = pfrb)) +
    theme(legend.position = "none")



## write coefficients ---

ip_coefs <- data.frame(
    mdm = c(mdm_pfrb_mpe, mdm_par_rmse),
    egg = c(egg_pfrb_mpe, egg_par_rmse),
    juv = c(juv_pfrb_mpe, juv_par_rmse),
    hydro = c(hydro_pfrb_mpe, hydro_par_rmse),
    disease = c(disease_pfrb_mpe, disease_par_rmse),
    agecomp = c(agecomp_pfrb_mpe, agecomp_par_rmse),
    row.names = c("PFRB MPE", "Parameter RMSE")
)

saveRDS(ip_coefs, file = here(dir_coi, "ip_coefs.rds"))
