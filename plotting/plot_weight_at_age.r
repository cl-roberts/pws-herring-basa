################################################################################

# Plot weight-at-ages

# Creates a box plot for most recent fishery and survey weight-at-age samples

# authors: CL Roberts

# inputs: BASA inputs 

# outputs: box plot of weight-at-ages

################################################################################


## set up ----

# pkgs

library(here)
library(ggplot2)
theme_set(theme_bw())
library(dplyr)
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

# data

waa_fall_2024 <- read.csv(here(dir_data, "2024", "2024-food-bait-waa.csv")) |>
    filter(age < 18) |>
    mutate(season = "Fall 2024", age = ifelse(age >= 9, 9, age)) |>
    filter(age != 9)

nrow(waa_fall_2024)

waa_spring_2025 <- read.csv(here(dir_data, "2025", "2025-asl-waa.csv")) |>
    # filter(age < 18) |>
    mutate(season = "Spring 2025")

waa_fall_2025 <- read.csv(here(dir_data, "2025", "2025-food-bait-waa.csv")) |>
    filter(age < 18) |>
    mutate(season = "Fall 2025", age = ifelse(age >= 9, 9, age)) |>
    filter(age != 9)

nrow(waa_fall_2025)

waa_spring_2025$age[waa_spring_2025$age == "9+"] <- 9

waa_spring_2025 <- waa_spring_2025 |>
    mutate(age = ifelse(as.numeric(age) >= 9, 9, as.numeric(age))) |>
    filter(age != 9) |>
    select(!gcon)

waa_data <- rbind(
    waa_fall_2024,
    waa_spring_2025,
    waa_fall_2025
)

waa_data$season <- factor(waa_data$season, levels = c("Fall 2024", "Spring 2025", "Fall 2025"))


## table ----

get_box_stats <- function(y, upper_limit = max(waa_data$weight) * 1.15) {
  return(data.frame(
    y = 0.95 * upper_limit,
    label = paste0(
        round(mean(y), 1), "\n",
        "n=", length(y)
    )
  ))
}

waa_plot <- ggplot(filter(waa_data, age >= 3), aes(x = factor(age), y = weight, fill = season)) +
    geom_boxplot() +
    stat_summary(
        fun.data = get_box_stats, 
        geom = "text", 
        # aes(label = round(..y.., 2)), 
        # hjust = 0.5, vjust = 0.9, 
        size = 2, fontface = "italic",
        position = position_dodge(width = 0.75)
    ) +
    scale_fill_manual(
        labels = c("Fall 2024" = "Fall 2024 Food/bait", "Spring 2025" = "Spring 2025 ASL", "Fall 2025" = "Fall 2025 Food/bait"),
        values = c("Fall 2024" = "grey50", "Spring 2025" = "grey70", "Fall 2025" = "grey90")
    ) +
    labs(x = "Age", y = "Weight (g)", title = "2024/25 Weight-at-age", fill = NULL) +
    scale_y_continuous(breaks = seq(0, 250, by = 50), limits = c(0, 275), expand=c(0, 0)) +
    theme_bw(base_size = 8) +
    theme(
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "bottom",
        axis.title = element_text(size=12),
        axis.text = element_text(size=12),
        plot.title = element_text(size=14)
    )
ggsave("weight-at-age.png", path = dir_figures, height = 4.5, width = 7)
