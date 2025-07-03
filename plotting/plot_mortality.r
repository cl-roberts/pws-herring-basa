################################################################################

# plot mortality time series and posteriors

# plots mortality time series estimates and parameter posteriors

# authors: CL Roberts

# inputs: BASA model inputs (model/PWS_ASA.dat) and outputs (from mcmc_out/)

# outputs: 
#   - plots: 
#   - tables: 

################################################################################


#### front matter ####

# choose TMB or ADMB
software <- "ADMB"

# load packages

library(ggplot2)
theme_set(theme_bw())
library(ggpubr)
library(pwsHerringBasa)
library(dplyr)
library(tidyr)
library(here)

# directory handling

dir_model <- here::here("model")

if (software == "ADMB") {
    dir_mcmc_out <- here(dir_model, "mcmc_out")
    dir_figures <- here("figures")
    dir_outputs <- here("data_outputs")
} else if (software == "TMB") {
    dir_mcmc_out <- here(dir_model, "mcmc_out_tmb")
    dir_figures <- here("figures/tmb")
    dir_outputs <- here("data_outputs/tmb")
} else {
    stop("choose valid software")
}

if (!dir.exists(dir_figures)) {
    dir.create(dir_figures)
}

# read model input data

raw.data <- read.data.files(dir_model)

# save global variables
nyr <- raw.data$PWS_ASA.dat$nyr
nage <- raw.data$PWS_ASA.dat$nage
start.year <- 1980
curr.year <- start.year+nyr
nyr.sim <- 0
years <- seq(start.year, curr.year+nyr.sim-1)

# read model output data
summer_survival <- read.csv(here(dir_mcmc_out, "adult_survival_effects_summer.csv"))
winter_survival <- read.csv(here(dir_mcmc_out, "adult_survival_effects_winter.csv"))

beta_posteriors <- read.csv(here(dir_mcmc_out, "iterations.csv")) |>
    select(contains("beta_mortality"))


#-------------------------------------------------------------------------------

#### plot mortality time series ####

# calculate mortalities
summer_mortality <- -log(summer_survival)
winter_mortality <- -log(winter_survival)

# combine mortalities into a single matrix for annual mortality rate 
annual_mortality_post <- summer_mortality + winter_mortality 


# calculate 50% and 95% intervals from posterior
annual_mortality <- expand.grid(1:nage-1, years)

quants <- c(.025, .25, .5, .75, .975)
for (q in seq(quants)) {
    mortality_quantiles <- apply(annual_mortality_post, MARGIN = 2, FUN = \(x) quantile(x, probs = quants[q]))    
    annual_mortality <- cbind(annual_mortality, mortality_quantiles)
}

colnames(annual_mortality) <- c("age", "year", quants)

# filter by the only relevant ages
mortality <- annual_mortality |>
    filter(age %in% c(3, 5)) 
mortality$age_range <- recode_factor(mortality$age, "3" = "Ages 3-4", "5" = "Ages 5-8")

# add raw disease data to mortality df for secondary axis
disease_covs <- raw.data$PWS_ASA_covariate.ctl$disease_covs 
disease_covs[disease_covs[,1] == -9,1] <- disease_covs[disease_covs[,1] == -9,2]

vhs <- data.frame(age = 3, year = years, disease = disease_covs[,3])
ich <- data.frame(age = 5, year = years, disease = disease_covs[,1])

disease <- rbind(ich, vhs)
disease$disease <- ifelse(disease$disease == -9, NA, disease$disease)

mortality <- merge(mortality, disease)

# plot

vhs_scalar <- 10*max(mortality$`0.975`)
vhs_mortality_ts <- ggplot(data = filter(mortality, age == 3), aes(x = year)) +
    geom_ribbon(aes(ymin = `0.025`, ymax = `0.975`), fill = "grey5", alpha = .25) +
    geom_ribbon(aes(ymin = `0.25`, ymax = `0.75`), fill = "grey5", alpha = .25) +
    geom_line(aes(y = `0.5`)) +
    geom_line(aes(y = disease*vhs_scalar)) +
    xlab("Year") + ylab("Instantaneous Mortality") +
    labs(title = "BASA Mortality Time Series", subtitle = "Ages 3-4") +
    scale_y_continuous(
        breaks = seq(0, max(mortality$`0.975`), by = .25),
        sec.axis = sec_axis(~./vhs_scalar, name="VHSV Prevalence")
    )

ich_scalar <- max(mortality$`0.975`)
ich_mortality_ts <- ggplot(data = filter(mortality, age == 5), aes(x = year)) +
    geom_ribbon(aes(ymin = `0.025`, ymax = `0.975`), fill = "grey5", alpha = .25) +
    geom_ribbon(aes(ymin = `0.25`, ymax = `0.75`), fill = "grey5", alpha = .25) +
    geom_line(aes(y = `0.5`)) +
    geom_line(aes(y = disease*ich_scalar)) +
    xlab("Year") + ylab("Instantaneous Mortality") +
    labs(subtitle = "Ages 5-8") +
    scale_y_continuous(
        breaks = seq(0, max(mortality$`0.975`), by = .25),
        sec.axis = sec_axis(~./ich_scalar, name="Ich. Prevalence")
    )

mortality_ts <- ggarrange(vhs_mortality_ts, ich_mortality_ts, ncol = 1, heights = c(1,.875))
# ggsave(here::here(dir_figures, "mortality_ts.pdf"), mortality_ts, height=8.5, width=11)
ggsave(here::here(dir_figures, "mortality_ts.png"), mortality_ts, height=5, width=7)


#-------------------------------------------------------------------------------

#### plot beta mortality parameters posteriors ####

beta_quants <- apply(beta_posteriors, MARGIN = 2, FUN = \(x) quantile(x, quants))

beta_hist_1 <- ggplot(data = beta_posteriors) +
    geom_histogram(aes(beta_mortality.1.), bins = 50, fill = "white", color = "black") +
    geom_vline(aes(xintercept = beta_quants[3,1]), color = "red") +
    labs(subtitle = paste0("Ich. Beta Mortality posterior (pre-2007)"))

beta_hist_2 <- ggplot(data = beta_posteriors) +
    geom_histogram(aes(beta_mortality.2.), bins = 50, fill = "white", color = "black") +
    geom_vline(aes(xintercept = beta_quants[3,2]), color = "red") +
    labs(subtitle = paste0("Ich. Beta Mortality posterior (post-2007)")) 

# beta_hist_3 <- ggplot(data = beta_posteriors) +
#     geom_histogram(aes(beta_mortality.3.), bins = 50, fill = "white", color = "black") +
#     geom_vline(aes(xintercept = beta_quants[3,3]), color = "red") +
#     labs(subtitle = paste0("VHSV Beta Mortality posterior"))

beta_mortality_posteriors <- ggarrange(
    beta_hist_1, beta_hist_2, 
    # beta_hist_3, 
    # nrow = 2, ncol = 2
    ncol = 2
)
# ggsave(here::here(dir_figures, "beta_mortality_posteriors.pdf"), beta_mortality_posteriors, height=8.5, width=11)
# ggsave(here::here(dir_figures, "beta_mortality_posteriors.png"), beta_mortality_posteriors, height=5, width=7)
ggsave(here::here(dir_figures, "beta_mortality_posteriors.png"), beta_mortality_posteriors, height=2.5, width=7)
