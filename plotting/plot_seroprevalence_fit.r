################################################################################

# plot seroprevalence fit

# plots BASA model fits to seroprevalence data for PWS herring 1980-present

# authors: John Trochta, Joshua Zahner, CL Roberts

# inputs: BASA model inputs (model/PWS_ASA.dat) and outputs (from mcmc_out/)

# outputs: 
#   - plots: 
#   - tables: 

################################################################################


#### front matter ####

# choose TMB or ADMB
software <- "ADMB"

# attach packages

library(pwsHerringBasa)
library(dplyr)
library(ggplot2)
theme_set(theme_bw())
library(ggpubr)
library(tidyr)

# directory handling

dir_model <- here::here("model")

if (software == "ADMB") {
    dir_mcmc_out <- here::here(dir_model, "mcmc_out")
    dir_figures <- here::here("figures")
    dir_outputs <- here::here("data_outputs")
} else if (software == "TMB") {
    dir_mcmc_out <- here::here(dir_model, "mcmc_out_tmb")
    dir_figures <- here::here("figures/tmb")
    dir_outputs <- here::here("data_outputs/tmb")
} else {
    stop("choose valid software")
}

if (!dir.exists(dir_figures)) {
    dir.create(dir_figures)
}

if (!dir.exists(dir_outputs)) {
    dir.create(dir_outputs)
}

# load data

raw.data <- read.data.files(dir_model)$PWS_ASA.dat
disease.data <- read.data.files(dir_model)$PWS_ASA_disease.dat


# local vars 

nyr <- raw.data$nyr
start.year <- 1980
curr.year <- start.year+nyr
nyr.sim <- 0
years <- seq(start.year, curr.year+nyr.sim-1)

#### data ####

# observed seroprevalence by age

sero_raw <- disease.data$vhsv_age_prevalence
sero_raw[sero_raw == -9] <- NA

sero_year_ind <- apply(!is.na(sero_raw), MARGIN = 1, FUN = all)
sero_years <- years[sero_year_ind]

sero_obs_pos <- sero_raw[sero_year_ind, seq(1, 19, by = 2)]
sero_obs_neg <- sero_raw[sero_year_ind, seq(2, 20, by = 2)]
colnames(sero_obs_pos) <- colnames(sero_obs_neg) <- paste0("age", 1:10) 
sero_obs_pos <- data.frame(year = sero_years, sero_obs_pos, vhsv = "pos")
sero_obs_neg <- data.frame(year = sero_years, sero_obs_neg, vhsv = "neg")

sero_obs <- rbind(
    pivot_longer(sero_obs_pos, cols = paste0("age", 1:10), 
                 names_to = "age", values_to = "seroprevalence"),
    pivot_longer(sero_obs_neg, cols = paste0("age", 1:10), 
                 names_to = "age", values_to = "seroprevalence")
    ) 

sero_obs$age <- gsub("age", "", sero_obs$age) |>
    as.numeric()

# predicted seroprevalence by age

sero_pred_raw <- read.csv(here::here(dir_mcmc_out, "Sero_pred.csv")) |>
    apply(MARGIN = 2, FUN = \(x) quantile(x, probs = c(.025, .5, .975))) |>
    t() |>
    as.data.frame() |>
    lapply(\(x) data.frame(
            year = sero_years,
            matrix(x, ncol = 20, byrow = TRUE)[sero_year_ind,]
        )
    ) 

sero_pred_raw[[1]]$prediction <- "lower"
sero_pred_raw[[2]]$prediction <- "estimate"
sero_pred_raw[[3]]$prediction <- "upper"

sero_pred_raw <- do.call(rbind, sero_pred_raw)

sero_pred_pos <- sero_pred_raw[, c(1, seq(2, 20, by = 2), 22)]
sero_pred_neg <- sero_pred_raw[, c(1, seq(3, 21, by = 2), 22)]
colnames(sero_pred_pos) <- colnames(sero_pred_neg) <- c("year", paste0("age", 1:10), "prediction") 
sero_pred_pos$vhsv <- "pos"
sero_pred_neg$vhsv <- "neg"

sero_pred <- rbind(
    pivot_longer(sero_pred_pos, cols = c(paste0("age", 1:10)), 
                 names_to = "age", values_to = "estimate"),
    pivot_longer(sero_pred_neg, cols = paste0("age", 1:10), 
                 names_to = "age", values_to = "estimate")
    ) |>
    pivot_wider(id_cols = c(year, vhsv, age), 
                names_from = prediction, 
                values_from = estimate)

sero_pred$age <- gsub("age", "", sero_pred$age) |>
    as.numeric()

sero <- cbind(sero_obs, select(sero_pred, estimate, lower, upper))

#### make plot ####

ggplot(filter(sero, vhsv == "pos"), aes(x = age)) +
    geom_line(aes(y = seroprevalence)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .25) +
    geom_line(aes(y = estimate), linetype = 2) +
    facet_wrap(~ year)

ggplot(filter(sero, vhsv == "neg"), aes(x = age)) +
    geom_line(aes(y = seroprevalence)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .25) +
    geom_line(aes(y = estimate), linetype = 2) +
    facet_wrap(~ year)

# ggplot() +
#     geom_col(data = sero_obs, aes(x = age, y = seroprevalence, fill = vhsv, group = vhsv)) +
#     facet_wrap(~ year)
