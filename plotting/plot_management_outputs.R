################################################################################

# plot management outputs

# plots BASA-estimated and derived time series parameters of interest to management

# authors: John Trochta, Joshua Zahner, CL Roberts

# inputs: BASA model inputs (model/PWS_ASA.dat) and outputs (from mcmc_out/)

# outputs: 
#   - plots: 
#       - recruitment (figures/recruitment_trajectory.pdf)
#       - spawning biomass (figures/biomass_trajectory.pdf)
#       - exploitation rate (figures/exploitation_rate.pdf)
#       - pre-fishery biomass posterior (figures/pfrb_posterior.pdf)
#       - a 4-panel composite plot of all outputs (figures/management_outputs.pdf)
#   - tables: 
#       - a single .csv file (data_outputs/outputs_for_management.csv) giving 
#         estimated recruitment, PFRB, catch, exploitation rate, and probability 
#         below harvest threshold with 95% credible intervals where relevant

################################################################################


#### front matter ####

# load packages

library(ggplot2)
library(ggpubr)
library(pwsHerringBasa)
library(dplyr)

# directory handling

dir_model <- here::here("model")

dir_mcmc_out <- here::here(dir_model, "mcmc_out")
dir_figures <- here::here("figures")
dir_outputs <- here::here("data_outputs")
dir_data <- here::here("data")

if (!dir.exists(dir_figures)) {
    dir.create(dir_figures)
}

# read model input data

raw.data <- read.data.files(dir_model)

# save global variables
nyr <- raw.data$PWS_ASA.dat$nyr
start.year <- 1980
curr.year <- start.year+nyr
nyr.sim <- 0
years <- seq(start.year, curr.year+nyr.sim-1)


#-------------------------------------------------------------------------------

#### calculate and plot outputs ####

# the following functions read in, compute posteriors, and plot each management 
# output of interest see ?pwsHerringBasa::compute.* and ?pwsHerringBasa::plot_* 
# for more details

# recruitment ----

recruit.df <- compute.recruitment(dir_mcmc_out, nyr, years) 

recruit_df <- recruit.df |>
    filter(.width == 0.5) |>
    select(year, recruits, .lower, .upper) |>
    mutate(year = as.integer(year)) |>
    rename(lower_50 = .lower, upper_50 = .upper)

recruit_df$lower_95 <- filter(recruit.df, .width == 0.95) |>
    select(.lower) |>
    unlist()

recruit_df$upper_95 <- filter(recruit.df, .width == 0.95) |>
    select(.upper) |>
    unlist()

recruit_forecast <- read.csv(here::here(dir_mcmc_out, "N_a_forecast.csv"))[,4] |>
    quantile(probs = c(0.5, 0.25, 0.75, 0.025, 0.975)) |>
    setNames(c("recruits", "lower_50", "upper_50", "lower_95", "upper_95"))

recruit_df <- recruit_df |>
    add_row(
        year = curr.year, recruits = recruit_forecast["recruits"], 
        lower_50 = recruit_forecast["lower_50"], upper_50 = recruit_forecast["upper_50"],
        lower_95 = recruit_forecast["lower_95"], upper_95 = recruit_forecast["upper_95"]
    )

recruit_df$forecast <- recruit_df$year == curr.year

recruit.plot <- plot_recruitment_posterior(recruit_df, years, legend = FALSE)

recruit.plot$layers$geom_line <- NULL

recruit.plot <- recruit.plot +
    ylim(c(0, 2200)) +
    geom_line(aes(x = year, y = recruits, linetype = forecast)) +
    geom_point(aes(x = year, y = recruits, shape = forecast), size = 2) +
    scale_linetype_manual(values = c(`TRUE` = NULL, `FALSE` = 1)) +
    scale_shape_manual(values = c(`TRUE` = 8, `FALSE` = NULL))



# biomass ----

biomass.df <- compute.biomass.traj(dir_mcmc_out, nyr)

biomass_df <- biomass.df |>
    filter(.width == 0.5) |>
    select(year, prob, biomass, .lower, .upper) |>
    mutate(year = c(years, curr.year)) |>
    rename(lower_50 = .lower, upper_50 = .upper)

biomass_df$lower_95 <- filter(biomass.df, .width == 0.95) |>
    select(.lower) |>
    unlist()

biomass_df$upper_95 <- filter(biomass.df, .width == 0.95) |>
    select(.upper) |>
    unlist()

biomass_df$forecast <- biomass_df$year == curr.year

biomass.plot <- plot_biomass_trajectory(biomass_df, c(years, curr.year), legend=FALSE)

biomass.plot$layers$geom_line <- NULL

biomass.plot <- biomass.plot +
    geom_line(aes(x = year, y = biomass/1000, linetype = forecast)) +
    geom_point(aes(x = year, y = biomass/1000, shape = forecast), size = 2) +
    scale_linetype_manual(values = c(`TRUE` = NULL, `FALSE` = 1)) +
    scale_shape_manual(values = c(`TRUE` = 8, `FALSE` = NULL))


# exploitation ----

exploit.df <- compute.exploit.rate(dir_mcmc_out, nyr)

exploit_df <- exploit.df$exploit.rate.df |>
    select(year, exploit, .lower, .upper) |>
    mutate(year = as.numeric(year)) |>
    rename(lower_95 = .lower, upper_95 = .upper)

exploit.rate.plot <- plot_exploit_rate(exploit_df, years)

# exploitation ----

pfrb.posterior <- compute.pfrb.posterior(dir_mcmc_out, nyr+1)
pfrb.posterior.plot <- plot_pfrb_posterior(pfrb.posterior$biomass.df, 
                                           pfrb.posterior$biomass.quants, 
                                           pfrb.posterior$prob.below.threshold,
                                           curr.year,
                                           font.size=14)

#-------------------------------------------------------------------------------

#### sensitivity analysis #### 

biomass_low_mdm <- read.csv(here::here(dir_mcmc_out, "PFRBiomass_low_mdm.csv"), header = FALSE) 
biomass_high_mdm <- read.csv(here::here(dir_mcmc_out, "PFRBiomass_high_mdm.csv"), header = FALSE)

biomass_sensitivity <- data.frame(
        year = 1980:curr.year, 
        biomass_low_mdm = apply(biomass_low_mdm, MARGIN = 2, FUN = median), 
        biomass_high_mdm = apply(biomass_high_mdm, MARGIN = 2, FUN = median)
    ) |>
    mutate(is_forecast = year == curr.year)

biomass_sensitivity_plot <- ggplot(biomass_df) +
    geom_ribbon(aes(x=year, ymin=lower_50/1000, ymax=upper_50/1000), alpha = .2) +
    geom_ribbon(aes(x=year, ymin=lower_95/1000, ymax=upper_95/1000), alpha = .2) +
    # geom_hline(yintercept=c(22, 42.5), linetype="dashed") +
    geom_line(aes(x = year, y = biomass/1000, linetype = "missing_mdm")) +
    geom_line(data = filter(biomass_sensitivity, !is_forecast), aes(x=year, y = biomass_low_mdm/1000, color = "low_mdm"), linewidth=1, alpha = .75) +
    geom_point(data = filter(biomass_sensitivity, is_forecast), aes(x=year, y = biomass_low_mdm/1000, color = "low_mdm"), size = 2, shape = 16, alpha = .75) +
    geom_line(data = filter(biomass_sensitivity, !is_forecast), aes(x=year, y = biomass_high_mdm/1000, color = "high_mdm"), linewidth=1, alpha = .75) +
    geom_point(data = filter(biomass_sensitivity, is_forecast), aes(x=year, y = biomass_high_mdm/1000, color = "high_mdm"), size = 2, shape = 16, alpha = .75) +
    scale_x_continuous("Year", breaks=seq(min(years), max(years), by=5), expand=c(.05, .05)) +
    ylab("Spring pre-fishery mature biomass (1000 tons)") +
    scale_linetype_manual(
        labels = c("missing_mdm" = "Missing MDM"),
        breaks = "missing_mdm",
        values = c("missing_mdm" = "dotted")
    ) +
    scale_color_manual(
        labels = c("low_mdm" = "25 MDM", "high_mdm" = "30 MDM"),
        breaks = c("low_mdm", "high_mdm"),
        values = c("low_mdm" = "#377EB8FF", "high_mdm" = "#E41A1CFF")
    ) +
    labs(title = "Biomass sensitivity to missing milt data", color = NULL, linetype = NULL) +
    coord_cartesian(clip="off") +
    theme_bw(base_size = 14) +
    theme(
        axis.title = element_text(size=12),
        axis.text = element_text(size=12),
        plot.title = element_text(size=14),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "bottom",
        legend.background = element_rect(fill = "grey90") 
    )

library(patchwork)
library(grid)

low_mdm_pct_diffs <- 100 * (biomass_df$biomass - biomass_sensitivity$biomass_low_mdm) / biomass_df$biomass
high_mdm_pct_diffs <- 100 * (biomass_df$biomass - biomass_sensitivity$biomass_high_mdm) / biomass_df$biomass

low_mdm_ci_width <- diff(quantile(biomass_low_mdm[,nyr+1], probs = c(.025, .975)))
high_mdm_ci_width <- diff(quantile(biomass_high_mdm[,nyr+1], probs = c(.025, .975)))
missing_mdm_ci_width <- filter(biomass_df, year == 2026) |>
    select(lower_95, upper_95) |>
    unlist() |>
    diff()

low_mdm_pct_diffs_ci <- 100 * (missing_mdm_ci_width - low_mdm_ci_width) / missing_mdm_ci_width
high_mdm_pct_diffs_ci <- 100 * (missing_mdm_ci_width - high_mdm_ci_width) / missing_mdm_ci_width

sink(here::here(dir_outputs, "sensitivity-report.txt"))
cat("low mdm average pct difference: ")
cat(mean(low_mdm_pct_diffs), "%")
cat("\nhigh mdm average pct difference: ")
cat(mean(high_mdm_pct_diffs), "%")
cat("\n\nlow mdm forecast pct difference: ")
cat(rev(low_mdm_pct_diffs)[1], "%")
cat("\nhigh mdm forecast pct difference: ")
cat(rev(high_mdm_pct_diffs)[1], "%")
cat("\n\nlow mdm prob under threshold: ")
cat(sum(biomass_low_mdm[,nyr+1] < 19958) / nrow(biomass_low_mdm), "%")
cat("\nhigh mdm prob under threshold: ")
cat(sum(biomass_high_mdm[,nyr+1] < 19958) / nrow(biomass_low_mdm), "%")
cat("\n\nlow mdm percent difference in CI width: ")
cat(low_mdm_pct_diffs_ci, "%")
cat("\nhigh mdm percent difference in CI width: ")
cat(high_mdm_pct_diffs_ci, "%")
sink()

biomass_sensitivity_plot_inset <- biomass_sensitivity_plot +
    labs(title = NULL, x = NULL, y = NULL) +
    xlim(c(2022, curr.year)) +
    ylim(c(0, 55)) +
    scale_y_continuous(breaks=c(0, 20, 40, 60, 80), limits = c(0, 80), expand = c(0, 0)) +
    # annotation_custom(biomass_sensitivity_plot_text_layer) +
    theme(
        legend.position = "none"
    )

biomass_sensitivity_plot_combined <- biomass_sensitivity_plot +
    inset_element(
        biomass_sensitivity_plot_inset,
        left = .4, right = .95, bottom = .4, top = .95
    )


#-------------------------------------------------------------------------------

#### save output files ####

# save individual plots

ggsave(here::here(dir_figures, "recruitment_trajectory.pdf"), plot = recruit.plot, height=8.5, width=11) 
ggsave(here::here(dir_figures, "biomass_trajectory.pdf"), plot = biomass.plot, height=8.5, width=11) 
ggsave(here::here(dir_figures, "exploitation_rate.pdf"), plot = exploit.rate.plot, height=8.5, width=11) 
ggsave(here::here(dir_figures, "pfrb_posterior.pdf"), plot = pfrb.posterior.plot, height=8.5, width=11) 
ggsave(here::here(dir_figures, "biomass-sensitivity.png"), plot = biomass_sensitivity_plot_combined, width = 7, height = 4.5)


# save 4-panel composite plot

management_outputs <- ggarrange(
    recruit.plot, biomass.plot, exploit.rate.plot, pfrb.posterior.plot,
    nrow=2,
    ncol=2
)
ggsave(here::here(dir_figures, "management_outputs.pdf"), management_outputs, 
       height=8.5, width=11)
ggsave(here::here(dir_figures, "management_outputs.png"), management_outputs, 
       height=8.5, width=11)


# save outputs-for-management.csv table

biomass.matrix <- biomass.df |> 
    filter(.width==0.95, year != "forecast") |> 
    select(biomass, .lower, .upper) |> 
    as.matrix() |> 
    round(2)
recruit.matrix <- recruit.df |> 
    filter(.width==0.95) |> 
    select(recruits, .lower, .upper) |> 
    as.matrix() |> 
    round(2)
exploit.matrix <- exploit.df$exploit.rate.df |> 
    filter(.width==0.95) |> 
    select(exploit, .lower, .upper) |> 
    as.matrix() |> 
    round(2)
prob.matrix <- biomass.df |> 
    filter(.width==0.95, year != "forecast") |> 
    select(prob) |> 
    as.matrix() |> 
    round(2)

total.catch.biomass <- compute.catch.biomass(data = raw.data$PWS_ASA.dat) |> 
    round(2) 


# FIX HOW FOOD/BAIT CATCH IS CALCULATED
# ADD FALL FISHERY WEIGHT-AT-AGE TO DATA FILE
waa_fall_2024 <- read.csv(here::here(dir_data, "2024", "2024-food-bait-waa.csv")) |>
    filter(age < 18) |>
    mutate(age = ifelse(age >= 9, 9, age))

mean_waa_fall_2024 <- c(0, tapply(waa_fall_2024$weight, waa_fall_2024$age, mean))

waa_fall_2025 <- read.csv(here::here(dir_data, "2025", "2025-food-bait-waa.csv")) |>
    filter(age < 18) |>
    mutate(age = ifelse(age >= 9, 9, age)) 

mean_waa_fall_2025 <- c(0, 0, tapply(waa_fall_2025$weight, waa_fall_2025$age, mean), 0)

total.catch.biomass[nyr-1] <- sum(mean_waa_fall_2024*raw.data$PWS_ASA.dat$foodbait[nyr-1,])
total.catch.biomass[nyr] <- sum(mean_waa_fall_2025*raw.data$PWS_ASA.dat$foodbait[nyr,])

biomass_forecast <- filter(biomass_df, year == curr.year) |>
    unlist() |>
    setNames(c("year", "prob", "biomass", "lower_50", "upper_50", "lower_95", "upper_95"))

final.table <- data.frame(
        years, recruit.matrix, biomass.matrix/1000, 
        total.catch.biomass, exploit.matrix, prob.matrix
    ) |>
    rbind(
        data.frame(
            years = curr.year, recruits = recruit_forecast["recruits"],
            .lower = recruit_forecast["lower_95"], .upper = recruit_forecast["upper_95"],
            biomass = biomass_forecast["biomass"]/1000, .lower.1 = biomass_forecast["lower_95"]/1000, 
            .upper.1 = biomass_forecast["upper_95"]/1000, total.catch.biomass = NA,
            exploit = NA, .lower.2 = NA, .upper.2 = NA, prob = pfrb.posterior$prob.below.threshold
        )
    )

names(final.table) <- c("Years",
                        "Median Age 3 (in millions)",
                        "Lower 95th Age 3 (in millions)",
                        "Upper 95th Age 3 (in millions)",
                        "Median Pre-fishery biomass (in 1000s metric tons)",
                        "Lower 95th Biomass (in 1000s metric tons)",
                        "Upper 95th Biomass (in 1000s metric tons)",
                        "Catch (metric tons)",
                        "Median Exploitation rate",
                        "Lower 95th ER",
                        "Upper 95th ER",
                        "Probability B<20K")
write.csv(final.table, here::here(dir_outputs, "outputs-for-management.csv"), row.names=FALSE)