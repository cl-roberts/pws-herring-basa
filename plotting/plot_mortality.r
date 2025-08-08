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
software <- "TMB"

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

dim(summer_survival)
dim(winter_survival) # 46 years

beta_posteriors <- read.csv(here(dir_mcmc_out, "iterations.csv")) |>
    select(contains("beta_mortality"))


#-------------------------------------------------------------------------------

#### plot mortality time series ####

# calculate mortalities

summer_mortality_age_3 <- -log(summer_survival)[seq(4, ncol(summer_survival), by = 10)]
winter_mortality_age_3 <- -log(winter_survival)[seq(4, ncol(winter_survival), by = 10)]

summer_mortality_age_5 <- -log(summer_survival)[seq(6, ncol(summer_survival), by = 10)]
winter_mortality_age_5 <- -log(winter_survival)[seq(6, ncol(winter_survival), by = 10)]

# combine mortalities into a single matrix for annual mortality rate 

annual_mortality_post_age_3 <- cbind(
    summer_mortality_age_3[,1:(nyr-1)] + winter_mortality_age_3[,2:nyr], 
    summer_mortality_age_3[,nyr] + winter_mortality_age_3[,(nyr+1)]
)
annual_mortality_post_age_5 <- cbind(
    summer_mortality_age_5[,1:(nyr-1)] + winter_mortality_age_5[,2:nyr], 
    summer_mortality_age_5[,nyr] + winter_mortality_age_5[,(nyr+1)]
)

# calculate 50% and 95% intervals from posterior
annual_mortality <- expand.grid(years, c(3, 5))

quants <- c(.025, .25, .5, .75, .975)
for (q in seq(quants)) {
    quantiles_age_3 <- apply(annual_mortality_post_age_3, MARGIN = 2, FUN = \(x) quantile(x, probs = quants[q]))    
    quantiles_age_5 <- apply(annual_mortality_post_age_5, MARGIN = 2, FUN = \(x) quantile(x, probs = quants[q]))    
    annual_mortality <- cbind(annual_mortality, c(quantiles_age_3, quantiles_age_5))
}

colnames(annual_mortality) <- c("year", "age", quants)

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
rownames(mortality) <- 1:nrow(mortality)

vhs_scalar <- 15*max(mortality$`0.975`)
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

ich_scalar <- .75*max(mortality$`0.975`)
ich_mortality_ts <- ggplot(data = filter(mortality, age == 5), aes(x = year)) +
    geom_ribbon(aes(ymin = `0.025`, ymax = `0.975`), fill = "grey5", alpha = .25) +
    geom_ribbon(aes(ymin = `0.25`, ymax = `0.75`), fill = "grey5", alpha = .25) +
    geom_line(aes(y = `0.5`)) +
    geom_line(aes(y = disease*ich_scalar)) +
    geom_vline(aes(xintercept = 2007), color = "red") +
    xlab("Year") + ylab("Instantaneous Mortality") +
    labs(subtitle = "Ages 5-8") +
    scale_y_continuous(
        breaks = seq(0, max(mortality$`0.975`), by = .25),
        sec.axis = sec_axis(~./ich_scalar, name="Ich. Prevalence")
    )

mortality_ts <- ggarrange(vhs_mortality_ts, ich_mortality_ts, ncol = 1, heights = c(1,.875))
ggsave(here::here(dir_figures, "mortality_ts.pdf"), mortality_ts, height=8.5, width=11)
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
ggsave(here::here(dir_figures, "beta_mortality_posteriors.png"), beta_mortality_posteriors, height=5, width=7)
# ggsave(here::here(dir_figures, "beta_mortality_posteriors.png"), beta_mortality_posteriors, height=2.5, width=7)