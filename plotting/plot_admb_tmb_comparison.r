################################################################################

# plot management outputs

# plot comparison between outputs of ADMB and TMB versions of BASA

# authors: CL Roberts

# inputs: BASA model inputs (model/PWS_ASA.dat) and outputs (from mcmc_out/)

# outputs: 
#   - plots: 
#   - tables: 

################################################################################

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
dir_outputs <- here("data_outputs")
dir_fig_tmb <- here("figures/tmb")


# model data
PWS_ASA <- read.data.files(dir_model)$"PWS_ASA.dat"  

# local vars

Y <- PWS_ASA$nyr
start.year <- 1980
curr.year <- start.year+Y
years <- start.year:(start.year + Y - 1)


## Make plots ----

# plot percent differences between biomass and recruitment time series

admb <- read.csv(here(dir_outputs, "outputs-for-management.csv"))

admb_biomass <- 1000*admb$Median.Pre.fishery.biomass..in.1000s.metric.tons.
admb_recruitment <- admb$Median.Age.3..in.millions.

admb_biomass_ci <- 1000*admb$Upper.95th.Biomass..in.1000s.metric.tons. - 
    1000*admb$Lower.95th.Biomass..in.1000s.metric.tons.
admb_recruitment_ci <- admb$Upper.95th.Age.3..in.millions. -
    admb$Lower.95th.Age.3..in.millions. 

tmb_biomass <- apply(do.call(cbind, Btilde_y), MARGIN = 1, FUN = median)
tmb_recruitment <- apply(do.call(cbind, age_3), MARGIN = 1, FUN = median)

tmb_biomass_ci <- apply(
    do.call(cbind, Btilde_y), MARGIN = 1, 
    FUN = \(x) quantile(x, .975) - quantile(x, .025)
)
tmb_recruitment_ci <- apply(
    do.call(cbind, age_3), MARGIN = 1, 
    FUN = \(x) quantile(x, .975) - quantile(x, .025)
)

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
ggsave(here(dir_fig_tmb, "admb_tmb_comparison.png"), width = 8, height = 3.5)


# plot likelihood histograms

admb_llik <- read.csv(here(dir_mcmc, "llikcomponents.csv"), header = FALSE)
tmb_llik <- read.csv(here(dir_mcmc_tmb, "likelihoods.csv"))

iters <- min(nrow(admb_llik), nrow(tmb_llik))
admb_llik <- admb_llik[1:iters,]
tmb_llik <- tmb_llik[1:iters,]

likelihoods <- data.frame(admb = admb_llik$V10, tmb = tmb_llik$V8)

ggplot(likelihoods) +
    geom_histogram(aes(x = admb, fill = "admb"), alpha = .25) +
    geom_histogram(aes(x = tmb, fill = "tmb"), alpha = .25) +
    geom_vline(aes(xintercept = median(admb)), color = "firebrick") +
    geom_vline(aes(xintercept = median(tmb)), color = "#404080") +
    scale_fill_manual(values=c("firebrick", "#404080")) +
    xlab("Total Negative Log-Likelihood")
ggsave(here(dir_fig_tmb, "admb_tmb_likelihoods.png"))


