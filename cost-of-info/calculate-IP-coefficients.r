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

base <- read_pfrb(here(dir_coi, "base", "PFRBiomass_base.csv"), all_years)

## calculate coefficients and plot surveys ----

## mdm

# read data

dirs_mdm <- list.files(here(dir_coi, "mdm"))

n_mdm <- length(dirs_mdm)

mdm_list <- vector(mode = "list", length = n_mdm)

for (i in 1:n_mdm) {
    mdm_year <- dirs_mdm[i]
    year <- unlist(stringr::str_split(mdm_year, "_"))[2] |>
        as.integer()
    dir_file <- here(dir_coi, "mdm", mdm_year, paste0("PFRBiomass_mdm_", year, ".csv"))
    mdm_list[[i]] <- read_pfrb(dir_file, all_years)
    mdm_list[[i]]$survey_year <- mdm_year
}

mdm <- do.call(rbind, mdm_list)

# calculate rmse
rmse(mdm$pfrb, rep(base$pfrb, n_mdm))

# plot PFRB with missing data against model with all data

ggplot() +
    geom_line(data = base, aes(x = year, y = pfrb)) +
    geom_line(data = mdm, aes(x = year, y = pfrb, color = survey_year)) +
    theme(legend.position = "none")
