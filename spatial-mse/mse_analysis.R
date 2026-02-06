################################################################################

# Analyze MSE outputs

# CL Roberts


################################################################################


#### set up ####

# attach packages
library(here)
library(ggplot2)
theme_set(theme_bw(base_size = 14))
library(dplyr)
library(tidyr)
library(ggh4x)  # for nested facet plots
library(gt)


# directory handling
dir_mse <- here("spatial-mse")
dir_out <- here(dir_mse, "mse_out")
dir_figures <- here(dir_mse, "figures")
if (!dir.exists(dir_figures)) dir.create(dir_figures)
dir_tables <- here(dir_mse, "tables")
if (!dir.exists(dir_tables)) dir.create(dir_tables)
# ---------------------------------------------------------------------------- #

#### read outputs ####

# read successful MSE iterations and only compare iterations in which all ems were successful
# list.files(pattern = "successful_iterations.csv", recursive = TRUE) |>
#     lapply(\(x) unlist(read.csv(x))) |>
#     Reduce(f = intersect)

ems <- list.dirs(dir_out, full.names = FALSE, recursive = FALSE)

# read performance metrics
metrics_list <- vector(mode = "list")
list_index <- 0
for (p in seq(ems)) {

    metrics_files <- list.files(here(dir_out, ems[p]), pattern = "performance-metrics.csv", recursive = TRUE) 

    for (f in seq(metrics_files)) {
        list_index <- list_index + 1

        metrics_list[[list_index]] <- read.csv(here(dir_out, ems[p], metrics_files[f]))
        metrics_list[[list_index]]$em <- ems[p]
        metrics_list[[list_index]]$state <- dirname(metrics_files[f])

    }

}

metrics <- do.call(rbind, metrics_list)

# ---------------------------------------------------------------------------- #

#### plot performance metrics ####

common_theme <- theme(
    axis.text.x = element_text(angle = -45, vjust = -.25, hjust = .5)
)

metrics <- metrics |>
    mutate(
      lambda = paste(as.numeric(gsub("[^0-9.]", "", state))),
      c = gsub("^([^_]+_){2}", "", state)
    ) |>
    mutate(
      c = factor(gsub("_", " ", c), c("all ages", "young ages", "new spawners")),
      em = factor(em, levels = c("base", "mvtweedie", "weibull", "spatialmodel"))  
    ) 

metrics$lambda <- recode_factor(metrics$lambda, `0` = "no\nmovement", `0.01` = "low\nmovement", `0.15` = "high\nmovement")

boxplot_summary <- function(x) {
    q <- quantile(x, probs = c(0.10, 0.25, 0.50, 0.75, 0.90))
    names(q) <- c("ymin", "lower", "middle", "upper", "ymax")
    return(data.frame(as.list(q)))
}

Btilde_ratio_plot <- ggplot(metrics, aes(x = em, y = Btilde_ratio)) +
    stat_summary(fun.data = boxplot_summary, geom = "boxplot", width = 0.75) +
    geom_hline(aes(yintercept = 1), color = "red") +
    facet_grid2(vars(lambda), vars(c), scales = "free", independent = "y", render_empty = FALSE) +
    xlab(NULL) + 
    ylab("Assessment accuracy ratio") +
    scale_y_continuous(breaks = seq(0, 1.5, by = 0.5)) +
    # ylab(bquote(frac(hat(B)[Y+10], B[Y+10]))) +
    common_theme 

Btilde_prob_plot <- ggplot(metrics, aes(x = em, y = Btilde_prob)) +
    stat_summary(fun.data = boxplot_summary, geom = "boxplot", width = 0.75) +
    geom_hline(aes(yintercept = 1), color = "red") +
    facet_grid2(vars(lambda), vars(c), scales = "free_y", independent = "y", render_empty = FALSE) +
    xlab(NULL) + 
    ylab("Assessment precision probability") +
    scale_y_continuous(breaks = seq(0, 1.5, by = 0.5)) +
    # ylab(bquote(P("|"^frac(B[Y+10] - hat(B)[Y+10], B[Y+10])~"|" < 0.3))) +
    common_theme

yield_difference_plot <- ggplot(metrics, aes(x = em, y = mean_yield_difference)) +
    stat_summary(fun.data = boxplot_summary, geom = "boxplot", width = 0.75) +
    geom_hline(aes(yintercept = 0), color = "red") +
    facet_grid2(vars(lambda), vars(c), scales = "free_y", independent = "y", render_empty = FALSE) +
    xlab(NULL) + 
    ylab("Mean difference in catch from HCR applied to true biomass") +
    common_theme

lost_yield_plot <- ggplot(metrics, aes(x = em, y = mean_lost_yield)) +
    stat_summary(fun.data = boxplot_summary, geom = "boxplot", width = 0.75) +
    geom_hline(aes(yintercept = 0), color = "red") +
    facet_grid2(vars(lambda), vars(c), scales = "free_y", independent = "y", render_empty = FALSE) +
    common_theme


# ---------------------------------------------------------------------------- #

#### make risk tables ####

head(metrics)

metrics_summary <- metrics |>
  group_by(em, state, lambda, c) |>
  summarize(
    Btilde_ratio = median(Btilde_ratio),
    Btilde_prob = median(Btilde_prob),
    mean_yield_difference = median(mean_yield_difference), 
    .groups = "drop"
  )

metrics_summary$c[metrics_summary$lambda == "\u039B = 0"] <- NA

Btilde_ratio_tbl <- metrics_summary |>
    select(em, lambda, c, Btilde_ratio) |>
    pivot_wider(names_from = c(lambda, c), values_from = Btilde_ratio)

colnames(Btilde_ratio_tbl) <- gsub("_NA", "", colnames(Btilde_ratio_tbl))
Btilde_ratio_tbl <- Btilde_ratio_tbl[,c(1,8,2,4,3,5,7,6)]

Btilde_prob_tbl <- metrics_summary |>
    select(em, lambda, c, Btilde_prob) |>
    pivot_wider(names_from = c(lambda, c), values_from = Btilde_prob)

colnames(Btilde_prob_tbl) <- gsub("_NA", "", colnames(Btilde_prob_tbl))
Btilde_prob_tbl <- Btilde_prob_tbl[,c(1,8,2,4,3,5,7,6)]

mean_yield_difference_tbl <- metrics_summary |>
    select(em, lambda, c, mean_yield_difference) |>
    pivot_wider(names_from = c(lambda, c), values_from = mean_yield_difference)

colnames(mean_yield_difference_tbl) <- gsub("_NA", "", colnames(mean_yield_difference_tbl))
mean_yield_difference_tbl <- mean_yield_difference_tbl[,c(1,8,2,4,3,5,7,6)]


gt(Btilde_ratio_tbl, rowname_col = "Name") |> 
  tab_spanner_delim(
    delim = "_", reverse = TRUE
  ) |>
  fmt_number(decimals = 2) |>
  cols_label(em = html("Estimation<br>model")) |>
  tab_style(
    style = cell_fill(color = "yellow"), locations = cells_body(columns = 2, rows = 1)
  ) |>
  tab_style(
    style = cell_fill(color = "yellow"), locations = cells_body(columns = 3, rows = 2)
  ) |>
  tab_style(
    style = cell_fill(color = "yellow"), locations = cells_body(columns = 4, rows = 2)
  ) |>
  tab_style(
    style = cell_fill(color = "yellow"), locations = cells_body(columns = 5, rows = 3)
  ) |>
  tab_style(
    style = cell_fill(color = "yellow"), locations = cells_body(columns = 6, rows = 3)
  ) |>
  tab_style(
    style = cell_fill(color = "yellow"), locations = cells_body(columns = 7, rows = 3)
  ) |>
  tab_style(
    style = cell_fill(color = "yellow"), locations = cells_body(columns = 8, rows = 2)
  ) 

gt(Btilde_prob_tbl, rowname_col = "Name") |> 
  tab_spanner_delim(
    delim = "_", reverse = TRUE
  ) |>
  fmt_number(decimals = 2) |>
  cols_label(em = html("Estimation<br>model")) |>
  tab_style(
    style = cell_fill(color = "yellow"), locations = cells_body(columns = 2, rows = 3)
  ) |>
  tab_style(
    style = cell_fill(color = "yellow"), locations = cells_body(columns = 3, rows = 2)
  ) |>
  tab_style(
    style = cell_fill(color = "yellow"), locations = cells_body(columns = 4, rows = 2)
  ) |>
  tab_style(
    style = cell_fill(color = "yellow"), locations = cells_body(columns = 5, rows = 3)
  ) |>
  tab_style(
    style = cell_fill(color = "yellow"), locations = cells_body(columns = 6, rows = 3)
  ) |>
  tab_style(
    style = cell_fill(color = "yellow"), locations = cells_body(columns = 7, rows = 3)
  ) |>
  tab_style(
    style = cell_fill(color = "yellow"), locations = cells_body(columns = 8, rows = 2)
  ) 

gt(mean_yield_difference_tbl, rowname_col = "Name") |> 
  tab_spanner_delim(
    delim = "_", reverse = TRUE
  ) |>
  fmt_number(decimals = 0) |>
  cols_label(em = html("Estimation<br>model")) |>
  tab_style(
    style = cell_fill(color = "yellow"), locations = cells_body(columns = 2, rows = 1)
  ) |>
  tab_style(
    style = cell_fill(color = "yellow"), locations = cells_body(columns = 3, rows = 3)
  ) |>
  tab_style(
    style = cell_fill(color = "yellow"), locations = cells_body(columns = 4, rows = 1)
  ) |>
  tab_style(
    style = cell_fill(color = "yellow"), locations = cells_body(columns = 5, rows = 1)
  ) |>
  tab_style(
    style = cell_fill(color = "yellow"), locations = cells_body(columns = 6, rows = 1)
  ) |>
  tab_style(
    style = cell_fill(color = "yellow"), locations = cells_body(columns = 7, rows = 1)
  ) |>
  tab_style(
    style = cell_fill(color = "yellow"), locations = cells_body(columns = 8, rows = 2)
  ) 

# ---------------------------------------------------------------------------- #

#### save figures and tables ####

ggsave(here(dir_figures, "Btilde-ratio-plot.png"), Btilde_ratio_plot, width = 7.5, height = 4)
ggsave(here(dir_figures, "Btilde-prob-plot.png"), Btilde_prob_plot, width = 7.5, height = 4)
ggsave(here(dir_figures, "yield-difference-plot.png"), yield_difference_plot, width = 7.5, height = 4)
ggsave(here(dir_figures, "lost-yield-plot.png"), lost_yield_plot, width = 7.5, height = 4)

write.csv(Btilde_ratio_tbl, here(dir_tables, "Btilde-ratio-tbl.csv"), row.names = FALSE)
write.csv(Btilde_prob_tbl, here(dir_tables, "Btilde-prob-tbl.csv"), row.names = FALSE)
write.csv(mean_yield_difference_tbl, here(dir_tables, "mean-yield-difference-tbl.csv"), row.names = FALSE)


# ---------------------------------------------------------------------------- #

#### plot mse examples ####

# fit to real data

  # mutate(
  #   em = "base", 
  #   Btilde_lower_95 = if_else(Year <= 2024, spatialmodel_om$Btilde_lower_95, Btilde_lower_95),
  #   Btilde_lower_50 = if_else(Year <= 2024, spatialmodel_om$Btilde_lower_50, Btilde_lower_50),
  #   Btilde_upper_50 = if_else(Year <= 2024, spatialmodel_om$Btilde_upper_50, Btilde_upper_50),
  #   Btilde_upper_95 = if_else(Year <= 2024, spatialmodel_om$Btilde_upper_95, Btilde_upper_95)
  # )

# no movement

spatialmodel_om <- read.csv(here(dir_out, "spatialmodel", "lambda_0_all_ages", "simulation-results.csv")) |>
  mutate(em = "spatialmodel")
spatialmodel_em <- read.csv(here(dir_out, "spatialmodel", "lambda_0_all_ages", "assessment-results.csv")) |>
  mutate(em = "spatialmodel")

base_om <- read.csv(here(dir_out, "base", "lambda_0_all_ages", "simulation-results.csv")) |>
    mutate(
      em = "base", 
      Btilde_lower_95 = spatialmodel_om$Btilde_lower_95,
      Btilde_lower_50 = spatialmodel_om$Btilde_lower_50,
      Btilde_upper_50 = spatialmodel_om$Btilde_upper_50,
      Btilde_upper_95 = spatialmodel_om$Btilde_upper_95
  )

base_em <- read.csv(here(dir_out, "base", "lambda_0_all_ages", "assessment-results.csv")) |>
  mutate(em = "base")

mvtweedie_om <- read.csv(here(dir_out, "mvtweedie", "lambda_0_all_ages", "simulation-results.csv")) |>
  mutate(
    em = "mvtweedie", 
      Btilde_lower_95 = spatialmodel_om$Btilde_lower_95,
      Btilde_lower_50 = spatialmodel_om$Btilde_lower_50,
      Btilde_upper_50 = spatialmodel_om$Btilde_upper_50,
      Btilde_upper_95 = spatialmodel_om$Btilde_upper_95
  )
mvtweedie_em <- read.csv(here(dir_out, "mvtweedie", "lambda_0_all_ages", "assessment-results.csv")) |>
  mutate(em = "mvtweedie")

weibull_om <- read.csv(here(dir_out, "weibull", "lambda_0_all_ages", "simulation-results.csv")) |>
  mutate(
    em = "weibull", 
      Btilde_lower_95 = spatialmodel_om$Btilde_lower_95,
      Btilde_lower_50 = spatialmodel_om$Btilde_lower_50,
      Btilde_upper_50 = spatialmodel_om$Btilde_upper_50,
      Btilde_upper_95 = spatialmodel_om$Btilde_upper_95
  )
weibull_em <- read.csv(here(dir_out, "weibull", "lambda_0_all_ages", "assessment-results.csv")) |>
  mutate(em = "weibull")


all_om <- rbind(base_om, mvtweedie_om, weibull_om, spatialmodel_om)
all_em <- rbind(base_em, mvtweedie_em, weibull_em, spatialmodel_em)

all_om$em <- factor(factor(all_om$em, levels = c("base", "mvtweedie", "weibull", "spatialmodel")))
all_em$em <- factor(factor(all_em$em, levels = c("base", "mvtweedie", "weibull", "spatialmodel")))

ggplot(all_om, aes(x = Year)) +
  geom_vline(aes(xintercept = 2024, linetype = "Simulation start")) +
  geom_ribbon(aes(ymin = Btilde_lower_95/1000, ymax = Btilde_upper_95/1000, fill = "True biomass"), alpha = .25) +
  geom_ribbon(aes(ymin = Btilde_lower_50/1000, ymax = Btilde_upper_50/1000, fill = "True biomass"), alpha = .25) +
  geom_line(data = filter(all_em, iteration <= 10), aes(y = biomass/1000, group = iteration, color = "Assessed biomass"), linewidth = .25) +
  facet_wrap(~ em, scales = "fixed") +
  scale_color_manual(values = "black") +
  scale_fill_manual(values = "black") +
  scale_linetype_manual(values = "dotted") +
  labs(color = NULL, fill = NULL, linetype = NULL, y = "Spring mature biomass (1000 mt)") +
  theme(legend.position = "bottom")
ggsave(here(dir_figures, "mse-example-lambda_0_all_ages.png"), width = 7.5, height = 4)


# low movement

spatialmodel_om <- read.csv(here(dir_out, "spatialmodel", "lambda_0.01_young_ages", "simulation-results.csv")) |>
  mutate(em = "spatialmodel")
spatialmodel_em <- read.csv(here(dir_out, "spatialmodel", "lambda_0.01_young_ages", "assessment-results.csv")) |>
  mutate(em = "spatialmodel")

base_om <- read.csv(here(dir_out, "base", "lambda_0.01_young_ages", "simulation-results.csv")) |>
  mutate(
    em = "base", 
      Btilde_lower_95 = spatialmodel_om$Btilde_lower_95,
      Btilde_lower_50 = spatialmodel_om$Btilde_lower_50,
      Btilde_upper_50 = spatialmodel_om$Btilde_upper_50,
      Btilde_upper_95 = spatialmodel_om$Btilde_upper_95
  )
base_em <- read.csv(here(dir_out, "base", "lambda_0.01_young_ages", "assessment-results.csv")) |>
  mutate(em = "base")

mvtweedie_om <- read.csv(here(dir_out, "mvtweedie", "lambda_0.01_young_ages", "simulation-results.csv")) |>
  mutate(
    em = "mvtweedie", 
      Btilde_lower_95 = spatialmodel_om$Btilde_lower_95,
      Btilde_lower_50 = spatialmodel_om$Btilde_lower_50,
      Btilde_upper_50 = spatialmodel_om$Btilde_upper_50,
      Btilde_upper_95 = spatialmodel_om$Btilde_upper_95
  )
mvtweedie_em <- read.csv(here(dir_out, "mvtweedie", "lambda_0.01_young_ages", "assessment-results.csv")) |>
  mutate(em = "mvtweedie")

weibull_om <- read.csv(here(dir_out, "weibull", "lambda_0.01_young_ages", "simulation-results.csv")) |>
  mutate(
    em = "weibull", 
      Btilde_lower_95 = spatialmodel_om$Btilde_lower_95,
      Btilde_lower_50 = spatialmodel_om$Btilde_lower_50,
      Btilde_upper_50 = spatialmodel_om$Btilde_upper_50,
      Btilde_upper_95 = spatialmodel_om$Btilde_upper_95
  )
weibull_em <- read.csv(here(dir_out, "weibull", "lambda_0.01_young_ages", "assessment-results.csv")) |>
  mutate(em = "weibull")

all_om <- rbind(base_om, mvtweedie_om, weibull_om, spatialmodel_om)
all_em <- rbind(base_em, mvtweedie_em, weibull_em, spatialmodel_em)

all_om$em <- factor(factor(all_om$em, levels = c("base", "mvtweedie", "weibull", "spatialmodel")))
all_em$em <- factor(factor(all_em$em, levels = c("base", "mvtweedie", "weibull", "spatialmodel")))

ggplot(all_om, aes(x = Year)) +
  geom_vline(aes(xintercept = 2024, linetype = "Simulation start")) +
  geom_ribbon(aes(ymin = Btilde_lower_95/1000, ymax = Btilde_upper_95/1000, fill = "True biomass"), alpha = .25) +
  geom_ribbon(aes(ymin = Btilde_lower_50/1000, ymax = Btilde_upper_50/1000, fill = "True biomass"), alpha = .25) +
  geom_line(data = filter(all_em, iteration <= 10), aes(y = biomass/1000, group = iteration, color = "Assessed biomass"), linewidth = .25) +
  facet_wrap(~ em, scales = "fixed") +
  scale_color_manual(values = "black") +
  scale_fill_manual(values = "black") +
  scale_linetype_manual(values = "dotted") +
  labs(color = NULL, fill = NULL, linetype = NULL, y = "Spring mature biomass (1000 mt)") +
  theme(legend.position = "bottom")
ggsave(here(dir_figures, "mse-example-lambda_0.01_young_ages.png"), width = 7.5, height = 4)



# high movement

spatialmodel_om <- read.csv(here(dir_out, "spatialmodel", "lambda_0.15_all_ages", "simulation-results.csv")) |>
  mutate(em = "spatialmodel")
spatialmodel_em <- read.csv(here(dir_out, "spatialmodel", "lambda_0.15_all_ages", "assessment-results.csv")) |>
  mutate(em = "spatialmodel")

base_om <- read.csv(here(dir_out, "base", "lambda_0.15_all_ages", "simulation-results.csv")) |>
  mutate(
    em = "base", 
      Btilde_lower_95 = spatialmodel_om$Btilde_lower_95,
      Btilde_lower_50 = spatialmodel_om$Btilde_lower_50,
      Btilde_upper_50 = spatialmodel_om$Btilde_upper_50,
      Btilde_upper_95 = spatialmodel_om$Btilde_upper_95
  )
base_em <- read.csv(here(dir_out, "base", "lambda_0.15_all_ages", "assessment-results.csv")) |>
  mutate(em = "base")

mvtweedie_om <- read.csv(here(dir_out, "mvtweedie", "lambda_0.15_all_ages", "simulation-results.csv")) |>
  mutate(
    em = "mvtweedie", 
      Btilde_lower_95 = spatialmodel_om$Btilde_lower_95,
      Btilde_lower_50 = spatialmodel_om$Btilde_lower_50,
      Btilde_upper_50 = spatialmodel_om$Btilde_upper_50,
      Btilde_upper_95 = spatialmodel_om$Btilde_upper_95
  )
mvtweedie_em <- read.csv(here(dir_out, "mvtweedie", "lambda_0.15_all_ages", "assessment-results.csv")) |>
  mutate(em = "mvtweedie")

weibull_om <- read.csv(here(dir_out, "weibull", "lambda_0.15_all_ages", "simulation-results.csv")) |>
  mutate(
    em = "weibull", 
      Btilde_lower_95 = spatialmodel_om$Btilde_lower_95,
      Btilde_lower_50 = spatialmodel_om$Btilde_lower_50,
      Btilde_upper_50 = spatialmodel_om$Btilde_upper_50,
      Btilde_upper_95 = spatialmodel_om$Btilde_upper_95
  )
weibull_em <- read.csv(here(dir_out, "weibull", "lambda_0.15_all_ages", "assessment-results.csv")) |>
  mutate(em = "weibull")


all_om <- rbind(base_om, mvtweedie_om, weibull_om, spatialmodel_om)
all_em <- rbind(base_em, mvtweedie_em, weibull_em, spatialmodel_em)

all_om$em <- factor(factor(all_om$em, levels = c("base", "mvtweedie", "weibull", "spatialmodel")))
all_em$em <- factor(factor(all_em$em, levels = c("base", "mvtweedie", "weibull", "spatialmodel")))

ggplot(all_om, aes(x = Year)) +
  geom_vline(aes(xintercept = 2024, linetype = "Simulation start")) +
  geom_ribbon(aes(ymin = Btilde_lower_95/1000, ymax = Btilde_upper_95/1000, fill = "True biomass"), alpha = .25) +
  geom_ribbon(aes(ymin = Btilde_lower_50/1000, ymax = Btilde_upper_50/1000, fill = "True biomass"), alpha = .25) +
  geom_line(data = filter(all_em, iteration <= 10), aes(y = biomass/1000, group = iteration, color = "Assessed biomass"), linewidth = .25) +
  # geom_line(data = all_om, aes(y = KI_Btilde_median/1000, color = "KI")) +
  # geom_errorbar(data = all_om, aes(ymin = KI_Btilde_lower_95/1000, ymax = KI_Btilde_upper_95/1000, color = "KI")) +
  facet_wrap(~ em, scales = "free_y") +
  scale_color_manual(values = "black") +
  scale_fill_manual(values = "black") +
  scale_linetype_manual(values = "dotted") +
  labs(color = NULL, fill = NULL, linetype = NULL, y = "Spring mature biomass (1000 mt)") +
  theme(legend.position = "bottom")
ggsave(here(dir_figures, "mse-example-lambda_0.15_all_ages.png"), width = 7.5, height = 4)


# ---------------------------------------------------------------------------- #

# boxplot for annual evostc report 


head(metrics)

metrics_base <- metrics |>
  filter(em == "base") |>
  select(Btilde_ratio, Btilde_prob, lambda, c) |>
  pivot_longer(cols = c(Btilde_ratio, Btilde_prob), names_to = "metric") 

metrics_base$metric <- factor(metrics_base$metric, levels = c("Btilde_ratio", "Btilde_prob"))

ggplot(metrics_base, aes(x = metric, y = value)) +
    stat_summary(fun.data = boxplot_summary, geom = "boxplot", width = 0.75) +
    geom_hline(aes(yintercept = 1), color = "red") +
    facet_grid2(vars(lambda), vars(c), scales = "fixed", render_empty = FALSE) +
    xlab(NULL) + 
    ylab("Performance metric") +
    scale_x_discrete(labels = c("Btilde_ratio" = "Assessment\naccuracy\nratio", Btilde_prob = "Assessment\nprecision\nprobability")) +
    scale_y_continuous(breaks = seq(0, 1.5, by = 0.5)) 
ggsave(here::here(dir_figures, "boxplot-base-model.png"), width = 7.5, height = 5)
