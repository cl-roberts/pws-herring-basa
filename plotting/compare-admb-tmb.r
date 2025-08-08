################################################################################

# compare ADMB and TMB model outputs

# generates plots and statistics for comparing BASA model outputs implemented in
# ADMB and TMB

# authors: CL Roberts

# inputs: BASA model outputs (from data_outputs/ and model/mcmc_out_tmb/)

# outputs: 
#   - plots: 
#       - compare_admb_tmb.png: a time series plot showing the percent difference
#                               in recruits and pre-fishery biomass between the
#                               ADMB and TMB models 

################################################################################

#### front matter ####

## attach packages

library(ggplot2)
library(dplyr)
library(pwsHerringBasa)
library(data.table)

## directory handling

dir_model <- here::here("model")
dir_outputs <- here::here("data_outputs")
dir_outputs_tmb <- here::here(dir_outputs, "tmb")
dir_mcmc_tmb <- here::here(dir_model, "mcmc_out_tmb")
dir_figures_tmb <- here::here("figures/tmb")

# local vars

Y <- read.data.files(dir_model)$"PWS_ASA.dat"$nyr
start.year <- 1980
curr.year <- start.year+Y

n_iters <- read.csv(here::here(dir_mcmc_tmb, "likelihoods.csv")) |>
    nrow()

## load data

# admb

admb_outputs <- read.csv(here::here(dir_outputs, "outputs-for-management.csv"))

biomass_admb <- admb_outputs[["Median.Pre.fishery.biomass..in.1000s.metric.tons."]]*1000 
recruits_admb <- admb_outputs[["Median.Age.3..in.millions."]] 

# tmb

tmb_outputs <- read.csv(here::here(dir_outputs_tmb, "outputs-for-management.csv"))

biomass_tmb <- tmb_outputs[["Median.Pre.fishery.biomass..in.1000s.metric.tons."]]*1000 
recruits_tmb <- tmb_outputs[["Median.Age.3..in.millions."]] 


## plot comparison

admb_tmb_compare <- data.frame(year = 1980:(curr.year-1))
admb_tmb_compare$biomass_pct <- 100*((biomass_tmb-biomass_admb) / biomass_tmb)
admb_tmb_compare$recruits_pct <- 100*((recruits_tmb-recruits_admb) / recruits_tmb)

admb_tmb_melt <- admb_tmb_compare |>
  as.data.table() |>
  melt(measure.vars = c("biomass_pct", "recruits_pct"))

compare_admb_tmb <- ggplot(admb_tmb_melt) + 
  geom_line(aes(x = year, y = value, color = variable)) +
  theme_bw()
ggsave(here::here(dir_figures_tmb, "compare_admb_tmb.png"), compare_admb_tmb)

