
################################################################################

# plot ADMB/TMB posteriors

# plot posteriors of parameters of ADMB and TMB versions of BASA

# authors: CL Roberts

# inputs: BASA model inputs (model/PWS_ASA.dat) and outputs (from mcmc_out/)

# outputs: 
#   - plots: 
#   - tables: 

################################################################################


## controls ----

# mcmc controls

chains <- 4
iter <- 2000
warmup <- 700

## set up ----

# attach packages

library(pwsHerringBasa)
library(dplyr)
library(ggplot2)
theme_set(theme_bw())
library(here)

# dir

dir_model <- here("model")

dir_mcmc <- here(dir_model, "mcmc_out")
dir_mcmc_tmb <- here(dir_model, "mcmc_out_tmb")

# read posteriors

admb_post <- read.csv(here(dir_mcmc, "iterations.csv"))[-(1:(warmup*chains)),]
tmb_post <- read.csv(here(dir_mcmc_tmb, "iterations.csv"))

tmp <- read.csv(here(dir_mcmc, "adult_survival_effects_summer.csv"))
all(tmp >= 0)

nrow(admb_post)
nrow(tmb_post)

param_post <- data.frame(
    software = c(rep("admb", (iter-warmup)*chains), rep("tmb", (iter-warmup)*chains))
)

pars <- names(tmb_post)[match(names(admb_post), names(tmb_post))] |>
    na.omit()

for (i in seq(pars)) {
    par <- pars[i]
    param_post[,i+1] <- c(admb_post[,par], tmb_post[,par])
}

names(param_post)[2:(length(pars)+1)] <- pars

param_post$seine_selex_alpha <- c(admb_post[,"seine_vuln_alpha"], tmb_post[,"seine_selex_alpha"])
param_post$seine_selex_beta <- c(admb_post[,"seine_vuln_beta"], tmb_post[,"seine_selex_beta"])


ggplot(param_post) +
    geom_point(aes(x = seine_selex_beta, y = beta_mortality.3., color = software), size = .25)

admb_param <- admb_post[,"beta_mortality.3."]
tmb_param <- tmb_post[,"beta_mortality.3."]

ggplot() +
    geom_histogram(aes(x = admb_param, fill = "admb"), alpha = .25) +
    geom_histogram(aes(x = tmb_param, fill = "tmb"), alpha = .25) +
    geom_vline(aes(xintercept = median(admb_param)), color = "firebrick") +
    geom_vline(aes(xintercept = median(tmb_param)), color = "#404080") +
    scale_fill_manual(values=c("firebrick", "#404080")) +
    xlab("parameter")
