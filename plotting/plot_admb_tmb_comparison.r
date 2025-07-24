################################################################################

# plot ADMB/TMB comparison

# plot comparison between outputs of ADMB and TMB versions of BASA

# authors: CL Roberts

# inputs: BASA model inputs (model/PWS_ASA.dat) and outputs (from mcmc_out/)

# outputs: 
#   - plots: 
#   - tables: 

################################################################################

## controls ----

method <- "MCMC"

## set up ----

# attach packages

library(pwsHerringBasa)
library(dplyr)
library(here)
library(ggplot2)
theme_set(theme_bw())
library(ggpubr)

# dir

dir_model <- here("model")

dir_mcmc <- here(dir_model, "mcmc_out")
dir_mcmc_tmb <- here(dir_model, "mcmc_out_tmb")

dir_rep <- here(dir_model, "rep_out")
dir_rep_tmb <- here(dir_model, "rep_out_tmb")

dir_outputs <- here("data_outputs")
dir_fig_tmb <- here("figures/tmb")


# model data
PWS_ASA <- read.data.files(dir_model)$"PWS_ASA.dat"  

# local vars

Y <- PWS_ASA$nyr
start.year <- 1980
curr.year <- start.year+Y
years <- start.year:(start.year + Y - 1)

## Read data ----

# admb <- read.csv(here(dir_outputs, "outputs-for-management.csv"))

# admb_biomass <- 1000*admb$Median.Pre.fishery.biomass..in.1000s.metric.tons.
# admb_recruitment <- admb$Median.Age.3..in.millions.

# nrow(admb_biomass_raw)
# nrow(tmb_biomass_raw)
# nrow(admb_recruitment_raw)
# nrow(tmb_recruitment_raw)

if (method == "MCMC") {

    admb_biomass_raw <- read.csv(here(dir_mcmc, "PFRBiomass.csv"))[-(Y+1)]
    admb_recruitment_raw <- read.csv(here(dir_mcmc, "Age3.csv"))

    admb_biomass <- admb_biomass_raw |>
        apply(MARGIN = 2, FUN = median)
    admb_recruitment <- admb_recruitment_raw |>
        apply(MARGIN = 2, FUN = median)

    admb_biomass_ci <- admb_biomass_raw |>
        apply(MARGIN = 2, FUN = \(x) quantile(x, .975) - quantile(x, .025))
    admb_recruitment_ci <- admb_recruitment_raw |>
        apply(MARGIN = 2, FUN = \(x) quantile(x, .975) - quantile(x, .025))

    tmb_biomass_raw <- read.csv(here(dir_mcmc_tmb, "PFRBiomass.csv"))[-(Y+1)]
    tmb_recruitment_raw <- read.csv(here(dir_mcmc_tmb, "Age3.csv"))

    tmb_biomass <- tmb_biomass_raw |>
        apply(MARGIN = 2, FUN = median)
    tmb_recruitment <- tmb_recruitment_raw |>
        apply(MARGIN = 2, FUN = median)

    tmb_biomass_ci <- tmb_biomass_raw |>
        apply(MARGIN = 2, FUN = \(x) quantile(x, .975) - quantile(x, .025))
    tmb_recruitment_ci <- tmb_recruitment_raw |>
        apply(MARGIN = 2, FUN = \(x) quantile(x, .975) - quantile(x, .025))

}

if (method == "ML") {

    admb_biomass <- data.reader(here(dir_rep, "PWS_ASA.rep"))[[37]]
    admb_recruitment <- data.reader(here(dir_rep, "PWS_ASA.rep"))[[42]] |>
        unlist()

    tmb_biomass <- data.reader(here(dir_rep_tmb, "ml-report.txt"))[[28]] |>
        unlist()
    tmb_recruitment <- data.reader(here(dir_rep_tmb, "ml-report.txt"))[[33]][,4] |>
        unlist()

}

## Make plots ----

# plot percent differences between biomass and recruitment time series


# comparison <- rbind(
#     data.frame(year = years,
#                biomass = admb_biomass,
#                recruitment = admb_recruitment,
#                software = "admb"),
#     data.frame(year = years,
#                biomass = tmb_biomass,
#                recruitment = tmb_recruitment,
#                software = "tmb")
# )

# ggplot(comparison) +
#     geom_line(aes(x = year, y = biomass, color = software)) +
#     ylab("Pre-fishery spawning biomass") +
#     labs(title = "ADMB vs. TMB model comparison")

# ggplot(comparison) +
#     geom_line(aes(x = year, y = recruitment, color = software)) +
#     ylab("Age-3 fish") +
#     labs(title = "ADMB vs. TMB model comparison")

comparison <- data.frame(
    year = years, 
    biomass = 100 * (admb_biomass - tmb_biomass) / admb_biomass,
    recruitment = 100 * (admb_recruitment - tmb_recruitment) / admb_recruitment,
    biomass_ci = 100 * (admb_biomass_ci - tmb_biomass_ci) / admb_biomass_ci,
    recruitment_ci = 100 * (admb_recruitment_ci - tmb_recruitment_ci) / admb_recruitment_ci
)

if (method == "MCMC") {

    comparison_median <- ggplot(comparison) +
        geom_line(aes(x = year, y = biomass, color = "biomass")) +
        geom_line(aes(x = year, y = recruitment, color = "recruitment")) +
        ylab("percent difference") +
        labs(title = "ADMB vs. TMB model comparison", subtitle = "Posterior median")

    comparison_ci <- ggplot(comparison) +
        geom_line(aes(x = year, y = biomass_ci, color = "biomass")) +
        geom_line(aes(x = year, y = recruitment_ci, color = "recruitment")) +
        ylab("percent difference") +
        labs(title = "", subtitle = "95% credible interval width")

    ggarrange(comparison_median, comparison_ci, ncol = 2, common.legend = TRUE, legend = "bottom")
    ggsave(here(dir_fig_tmb, "admb_tmb_comparison_mcmc.png"), width = 8, height = 3.5)

}

if (method == "ML") {

    comparison_mean <- ggplot(comparison) +
        geom_line(aes(x = year, y = biomass, color = "biomass")) +
        geom_line(aes(x = year, y = recruitment, color = "recruitment")) +
        ylab("percent difference") +
        labs(title = "ADMB vs. TMB model comparison", subtitle = "Maximum Likelihood")

    ggsave(here(dir_fig_tmb, "admb_tmb_comparison_ml.png"), width = 6, height = 3.5)

}


# plot likelihood histograms

admb_llik <- read.csv(here(dir_mcmc, "llikcomponents.csv"), header = FALSE)
tmb_llik <- read.csv(here(dir_mcmc_tmb, "likelihoods.csv"))

iters <- min(nrow(admb_llik), nrow(tmb_llik))
admb_llik <- admb_llik[1:iters,]
tmb_llik <- tmb_llik[1:iters,]

likelihoods <- data.frame(admb = admb_llik$V14, tmb = tmb_llik$V8)

ggplot(likelihoods) +
    geom_histogram(aes(x = admb, fill = "admb"), alpha = .25, bins = 50) +
    geom_histogram(aes(x = tmb, fill = "tmb"), alpha = .25, bins = 50) +
    geom_vline(aes(xintercept = median(admb)), color = "firebrick") +
    geom_vline(aes(xintercept = median(tmb)), color = "#404080") +
    scale_fill_manual(values=c("firebrick", "#404080")) +
    xlab("Total Negative Log-Likelihood")
ggsave(here(dir_fig_tmb, "admb_tmb_likelihoods.png"))



# plot penalties histograms

admb_pen <- admb_llik$V12
tmb_pen <- read.csv(here(dir_mcmc_tmb, "penalties.csv"))[,4]

admb_pen <- admb_pen[1:iters]
tmb_pen <- tmb_pen[1:iters]

sum(tmb_pen > 0)

penalties <- data.frame(admb = admb_pen, tmb = tmb_pen)

sum(penalties$admb)

ggplot(penalties) +
    geom_histogram(aes(x = admb, fill = "admb"), alpha = .25) +
    geom_histogram(aes(x = tmb, fill = "tmb"), alpha = .25) +
    geom_vline(aes(xintercept = median(admb)), color = "firebrick") +
    geom_vline(aes(xintercept = median(tmb)), color = "#404080") +
    scale_fill_manual(values=c("firebrick", "#404080")) +
    xlab("Penalties") + 
    coord_trans(y = "sqrt")

# admb_survival <- read.csv(here(dir_mcmc, "adult_survival_effects_summer.csv"))
# tmb_survival <- read.csv(here(dir_mcmc_tmb, "summer_survival.csv"))

# sum(admb_survival > .99)
# sum(tmb_survival > .99)
