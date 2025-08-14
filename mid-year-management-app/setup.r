#---- read and wrangle data ---#

# pkgs

library(dplyr)
library(stringr)
library(here)
library(ggplot2)
theme_set(theme_bw())
library(tidyr)
library(tibble)


# dir

if(interactive()) {

    app_dir <- grep("\\bapp.r\\b", list.files(recursive = TRUE), value = TRUE) |>
        here() |>
        dirname()

} else {

    app_dir <- "mid-year-management-app"

}

# scripts 

source(here(app_dir, "helper.r"))

# colors 

theme_navy <- "#072f49"
theme_carolina <- "#3F78A7"
theme_argentina <- "#b5d8f2"
theme_gold <- "#edbd03"
theme_offwhite <- "#f2f2f2"
theme_blood <- "#ba141a"
theme_transparent <- "#ffffff00"

# forecast controls

waa_average_years <- 10
perc_female_forecast_years <- 10    

# read model input data

model_data <- data.reader(here(app_dir, "data", "PWS_ASA.dat"))

nyr <- model_data[[1]]
curr_year <- 1980 + nyr
perc_females <- model_data[[11]]
waa <- model_data[[4]]
pk <- model_data[[7]]

# averaged quantities used in forecast

colnames(waa) <- paste("Age", c(0:8, "9+"))
waa <- waa[,colnames(waa) %in% paste("Age", c(3:8, "9+"))]
waa_forecast <- window_average(waa, waa_average_years)

perc_female_forecast <- mean(perc_females[(nyr-perc_female_forecast_years+1):nyr,])
threshold <- 22000

# read mcmc results

mcmc_results <- read.csv(here(app_dir, "data", "mid-year-management.csv"))

n_iters <- nrow(mcmc_results)

mean_log_rec <- mcmc_results$mean_log_rec
Btilde_forecast <- mcmc_results$Btilde_forecast
logmdm_c <- mcmc_results$logmdm_c
N_a_forecast <- mcmc_results[,grepl("N_a_forecast_age", names(mcmc_results))]
maturity <- cbind(mcmc_results$mat_age3, mcmc_results$mat_age4, 1, 1, 1, 1, 1)
survival_forecast <- mcmc_results[,grepl("survival_forecast_age", names(mcmc_results))]
Btilde_agecomp_forecast <- mcmc_results[,grepl("Btilde_agecomp", names(mcmc_results))]
Ntilde_agecomp_forecast <- mcmc_results[,grepl("Ntilde_agecomp", names(mcmc_results))]


# calculate one-year forecast quantities

Btilde_forecast_tons <- Btilde_forecast*1.10231

Btilde_forecast_med <- median(Btilde_forecast_tons) |>
    round() |>
    prettyNum(big.mark = ",")
Btilde_forecast_ci <- quantile(Btilde_forecast_tons, c(.025, .975)) |>
    round() |>
    prettyNum(big.mark = ",") |>
    paste(collapse = ", ")
prob_below_threshold <- paste0(
    round(100*sum(Btilde_forecast_tons < threshold)/n_iters, 1),
    "%"
)

# plot one-year forecast 

biomass_forecast_plot <- ggplot(data.frame(biomass = Btilde_forecast_tons)) +
    geom_histogram(aes(x = biomass, y = after_stat(count / sum(count))), 
                    fill = "NA", color = "black", bins = 50) +
    xlab("Mature Biomass (tons)") +
    ylab("Frequency") +
    labs(title = paste(curr_year, "Biomass Forecast Posterior Distribution"), 
         color = NULL) +
    geom_vline(aes(xintercept = threshold, color = "22,000-ton Threshold")) +
    scale_color_manual(values = theme_gold) +
    theme(legend.position.inside = c(.8, .8), legend.position = "inside")


# make agecomp forecast tables

colnames(Ntilde_agecomp_forecast) <- paste("Age", c(0:8, "9+"))
colnames(Btilde_agecomp_forecast) <- paste("Age", c(0:8, "9+"))

Ntilde_agecomp_tbl <- apply(
    Ntilde_agecomp_forecast, MARGIN = 2, FUN = \(x) quantile(x, c(.5, .025, .975))
    ) |>
    as.data.frame() |>
    select(contains(paste("Age", 3:9)))

Btilde_agecomp_tbl <- apply(
    Btilde_agecomp_forecast, MARGIN = 2, FUN = \(x) quantile(x, c(.5, .025, .975))
    ) |>
    as.data.frame() |>
    select(contains(paste("Age", 3:9)))

Ntilde_agecomp_tbl <- cbind(
    data.frame(Percentile = rownames(Ntilde_agecomp_tbl),
               Type = "Numbers"),
    Ntilde_agecomp_tbl
)
Btilde_agecomp_tbl <- cbind(
    data.frame(Percentile = rownames(Btilde_agecomp_tbl),
               Type = "Biomass"),
    Btilde_agecomp_tbl
)

agecomp_tbl <- rbind(Ntilde_agecomp_tbl, Btilde_agecomp_tbl)


# make agecomp forecast plots

agecomp_posteriors <- rbind(
        cbind(
            data.frame(Type = "Numbers"), 
            Ntilde_agecomp_forecast
        ), 
        cbind(
            data.frame(Type = "Biomass"), 
            Btilde_agecomp_forecast
        )
    ) |>
    select(Type, contains(paste("Age", 3:9))) |>
    pivot_longer(cols = !Type, names_to = "Age", values_to = "Proportion")

agecomp_posteriors$Type <- factor(agecomp_posteriors$Type, c("Numbers", "Biomass"))

naa_agecomp_fig <- ggplot(filter(agecomp_posteriors, Type == "Numbers")) +
    geom_boxplot(aes(x = Age, y = Proportion), outliers = FALSE) +
    # facet_wrap(~ Type) +
    xlab(NULL) +
    labs(color = NULL, subtitle = "Numbers-at-age Composition") 

baa_agecomp_fig <- ggplot(filter(agecomp_posteriors, Type == "Biomass")) +
    geom_boxplot(aes(x = Age, y = Proportion), outliers = FALSE) +
    # facet_wrap(~ Type) +
    xlab(NULL) +
    labs(color = NULL, subtitle = "Biomass Age Composition") 

# calculate numbers-at-age forecast

colnames(N_a_forecast) <- paste("Age", c(0:8, "9+")) 
N_a_forecast <- N_a_forecast |>
    select(contains(paste("Age", 3:9))) 

# calculate two-year biomass forecast quantities

colnames(maturity) <- paste("Age", c(3:8, "9+"))

colnames(survival_forecast) <- paste("Age", c(0:8, "9+"))
survival_forecast_age2 <- survival_forecast |>
    select(contains(paste("Age", 2))) 
survival_forecast <- survival_forecast |>
    select(contains(paste("Age", 3:9))) 

# get table ready for download

forecastTbl <- tibble(
    Quantity = "Mature Biomass Forecast", 
    Year = curr_year, 
    Age = "3+",
    Units = "Tons",
    `50%` = median(Btilde_forecast_tons), 
    `2.5%` = quantile(Btilde_forecast_tons, .025),
    `97.5%` = quantile(Btilde_forecast_tons, .975)
)

agecompTbl <- agecomp_tbl |>
    rename("Quantity" = Type) |>
    mutate(Year = curr_year, Units = "Proportion") |>
    pivot_longer(cols = contains("Age "), names_to = "Age") |>
    pivot_wider(names_from = Percentile, values_from = value)

agecompTbl$Quantity <- paste(forecastTbl$Quantity, "agecomp")
agecompTbl$Age <- gsub("Age", "", agecompTbl$Age)
