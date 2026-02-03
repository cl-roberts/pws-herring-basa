################################################################################

# Plot catch

# Creates a bar plot for most recent fishery catch

# authors: CL Roberts

# inputs: BASA inputs 

# outputs: bar plot of catch-at-ages

################################################################################


## set up ----

# pkgs

library(here)
library(ggplot2)
theme_set(theme_bw())
library(dplyr)
library(tidyr)
library(pwsHerringBasa)

# dir

dir_model <- here::here("model")
dir_data <- here::here("data")

# dir_mcmc_out <- here::here(dir_model, "mcmc_out")
dir_figures <- here::here("figures")
dir_outputs <- here::here("data_outputs")

if (!dir.exists(dir_figures)) {
    dir.create(dir_figures)
}


## data ----

data_files <- read.data.files(dir_model)
nyr <- data_files$PWS_ASA.dat$nyr

fall_2024_catch <- data_files$PWS_ASA.dat$foodbait_catch[nyr-1,]
fall_2025_catch <- data_files$PWS_ASA.dat$foodbait_catch[nyr,]

waa_fall_2024 <- read.csv(here(dir_data, "2024", "2024-food-bait-waa.csv")) |>
    filter(age < 18) |>
    mutate(age = ifelse(age >= 9, 9, age))

waa_fall_2024[waa_fall_2024$age == 9,]

mean_waa_fall_2024 <- tapply(waa_fall_2024$weight, waa_fall_2024$age, mean)

waa_fall_2025 <- read.csv(here(dir_data, "2025", "2025-food-bait-waa.csv")) |>
    filter(age < 18) |>
    mutate(age = ifelse(age >= 9, 9, age)) 

mean_waa_fall_2025 <- tapply(waa_fall_2025$weight, waa_fall_2025$age, mean)

catch <- data.frame(
        age = c(0:8, "9+"),
        fall_2024 = fall_2024_catch, #* c(0, mean_waa_fall_2024) 
        fall_2025 = fall_2025_catch #* c(0, 0, mean_waa_fall_2025, 0)
    ) |>
    pivot_longer(cols = !age, names_to = "fishery", values_to = "catch")

sum(fall_2024_catch * c(0, mean_waa_fall_2024))
sum(fall_2025_catch * c(0, 0, mean_waa_fall_2025, 0))

## plot ----

fishery_catch_plot <- ggplot(catch) +
    geom_col(aes(x = age, y = catch, fill = fishery), color = "black", position = "dodge") +
    scale_fill_manual(
        labels = c("fall_2024" = "Fall 2024 Food/bait", "fall_2025" = "Fall 2025 Food/bait"),
        values = c("fall_2024" = "grey50", "fall_2025" = "grey90")
    ) +
    xlab("Age") +
    ylab("Catch (millions)") +
    labs(title = "Commercial food/bait fishery catch-at-age", fill = NULL) +
    # scale_y_continuous(breaks = seq(0, 3, by = 0.5), limits = c(0, 3.25), expand=c(0, 0)) +
    theme_bw(base_size = 8) +
    theme(
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "bottom",
        axis.title = element_text(size=12),
        axis.text = element_text(size=12),
        plot.title = element_text(size=14)
    )
ggsave("catch-at-age.png", path = dir_figures, height = 4.5, width = 7)

