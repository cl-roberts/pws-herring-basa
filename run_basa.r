# basa_run.r
# Created by John Trochta
# Date created:  06/08/2019
# Summary:
# This script runs the Bayesian ASA model using the No-U-turn (NUTS) MCMC sampler to obtain posteriors.
# The adnuts package developed by Cole Monnahan to run NUTS with ADMB (hence adnuts) is used.
# I highly recommend users should read the following for more details on NUTS application to stock assessments:
#   Cole C Monnahan, Trevor A Branch, James T Thorson, Ian J Stewart, Cody S Szuwalski, 
#         Overcoming long Bayesian run times in integrated fisheries stock assessments, 
#         ICES Journal of Marine Science, Volume 76, Issue 6, November-December 2019, 
#         Pages 1477–1488, https://doi-org.eres.qnl.qa/10.1093/icesjms/fsz059

# Some annotated guidance on diagnostic checks is provided below:
# 1) If divergent transitions>0, increase target acceptance rate to reduce step size
#    e.g. control = list(adapt_delta = 0.95)
# 2) IF extreme global correlations, pass dense matrix estimated from previous run
#    e.g. control-list(metric=M) where M is matrix in untransformed space
#    for ADMB models, use MLE covairance with control=list(metric="mle")

# With this adnuts, most important diagnostics are the:
#   1) ESS (accounts for autocorrelation)-500 ESS is sufficient for most quantities
#   2) Potential Scale reduction (R hat)-R hat fails if >1.1
#   3) No max tree depths exceeded (<12)
#   4) 0% divergences

# For divergence diagnoses and resolutions:
# https://discourse.mc-stan.org/t/divergent-transitions-a-primer/17099

#################################################################

packages <- c("data.table", "tidyverse", "adnuts", "rstan", "here", "r4ss")
if(length(packages[which(packages %in% rownames(installed.packages()) == FALSE )]) > 0){
  install.packages(packages[which(packages %in% rownames(installed.packages()) == FALSE)])
}
lapply(packages, library, character.only = TRUE)

# library(data.table)
# library(tidyverse)
# library(adnuts)
# library(snowfall)
# library(rstan)
# library(r4ss)
source(paste0(here::here("functions/"), "/run_basa.r"))

run <- run.basa(here::here("model"))
# Extracts NUTS stats (energy, leapfrog transitions,etc)
mon <- run$fit1$monitor 
x <- extract_sampler_params(run$fit1)

# Quick check for divergences, Gelman-Ruben statistic, and tail ESS
n.divergences <- sum(x$divergent__)/nrow(x)
r.hat <- max(mon[, "Rhat"]) <= 1.1
ess <- min(mon[, "Tail_ESS"]) >= 500

n.divergences
r.hat
ess

# If this returns TRUE, diagnostic convergence checks pass
ifelse(
  n.divergences < 0.001 & r.hat & ess, 
  print("Diagnostics pass. Convergence likely."), 
  print("One or more diagnostic checks failed.")
)

# Write summary of parameter posteriors (medians, percentiles, etc)
write.csv(as.data.frame(mon), file="mcmc_out/table_par_posterior_summary.csv")

# Write all MCMC samples of the parameters
mcmc.samps <- data.frame(matrix(run$fit1$samples, ncol=dim(run$fit1$samples)[3], byrow=FALSE))
names(mcmc.samps) <- run$fit1$par_names
write.csv(mcmc.samps, file="mcmc_out/iterations.csv", row.names=FALSE)

## Examine the slowest mixing parameters
# slow <- names(sort(mon[,"n_eff"]))[1:8]
# pairs_admb(fit=fit, pars=slow)

#Create summary file of NUTS/MCMC diagnostics
sum.dia <- data.frame(
    divergences.from.extract.function=n.divergences,
    min.ESS=min(mon[, "n_eff"]),
    which.min.ESS=names(which.min(mon[, "n_eff"])),
    max.Rhat=max(mon[, "Rhat"]),
    which.max.Rhat=names(which.max(mon[, "Rhat"])),
    time.elapsed=run$time
)

write.table(sum.dia, file="mcmc_out/table_convergence_diagnostics.csv", sep=",", row.names=FALSE)
saveRDS(run$fit1, file="mcmc_out/NUTS_fit.RDS")

# Launch Shiny App to check diagnostics online
# launch_shinyadmb(fit.1)

rm(list = ls(all.names = TRUE))

# plot_marginals(fit.1)
# plot_sampler_params(fit.1)
# plot_uncertainties(fit.1)

# pairs_admb(fit.1,pars=1:20)
