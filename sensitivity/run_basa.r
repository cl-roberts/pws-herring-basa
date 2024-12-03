################################################################################

# Run BASA model and obtain parameter posteriors

# Created by John Trochta
# Date created:  06/08/2019

# Summary:
# This script runs the Bayesian ASA model using the No-U-turn (NUTS) MCMC sampler 
# to obtain posteriors.

# The adnuts package developed by Cole Monnahan to run NUTS with ADMB (hence adnuts) is used.

# I highly recommend users should read the following for more details on NUTS 
# application to stock assessments:
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

################################################################################

#### front matter ####

# load packages

library(data.table)
library(adnuts)
library(snowfall)
library(rstan)
library(r4ss)
library(pwsHerringBasa)

#-------------------------------------------------------------------------------

#### compile model, calculate ESS's execute NUTS for parameter posteriors ####

dir_model <- here::here("sensitivity/model")

template.files <- here::here(dir_model)
print(template.files)
setwd(template.files)

  inits <- init.admb.params(1)

if(dir.exists("mcmc_out")){
    system("rm mcmc_out/*.csv")
}

start.time <- Sys.time()
fit.1 <- adnuts::sample_nuts(model='./PWS_ASA', 
                    path=template.files,
                    iter=1000,
                    warmup=500,
                    #warmup=100,
                    duration = 2,
                    init=inits, seeds=1867, chains=1, cores=1,
                    mceval=TRUE,
                    control=list(
                        adapt_delta=0.9,
                        #max_treedepth=16,
                        metric="mle"
                    )
)
end.time <- Sys.time()
total.time <- end.time - start.time

run <- list(fit1=fit.1, time=total.time)

# see ?pwsHerringBasa::run.basa() for more details

# run <- run.basa(here::here("sensitivity/model"), 
#                 n.samples = 1000, n.warmup = 500, 
#                 n.chains = 1, n.time = 2)

#-------------------------------------------------------------------------------

#### extract model run metadata and diagnostics ####

# Extracts NUTS stats (energy, leapfrog transitions,etc)
mon <- monitor(run$fit1$samples, warmup=run$fit1$warmup, print=FALSE)
x <- extract_sampler_params(run$fit1)

# Quick check for divergences & Gelman-Ruben statistic
n.divergences <- sum(x$divergent__)/nrow(x)
r.hat <- max(mon[, "Rhat"])<=1.1

n.divergences
r.hat

## Examine the slowest mixing parameters
# slow <- names(sort(mon[,"n_eff"]))[1:8]
# pairs_admb(fit=fit, pars=slow)

# summarize NUTS/MCMC diagnostics
sum.dia <- data.frame(divergences.from.extract.function=sum(x$divergent__)/nrow(x),
                      min.ESS=min(mon[, "n_eff"]),
                      which.min.ESS=names(which.min(mon[, "n_eff"])),
                      max.Rhat=max(mon[, "Rhat"]),
                      which.max.Rhat=names(which.max(mon[, "Rhat"])),
                      time.elapsed=run$time)

#-------------------------------------------------------------------------------

#### save information from model run ####

# Write summary of parameter posteriors (medians, percentiles, etc)
write.csv(as.data.frame(mon), file="mcmc_out/table_par_posterior_summary.csv")

# Write all MCMC samples of the parameters
mcmc.samps <- data.frame(matrix(run$fit1$samples, ncol=dim(run$fit1$samples)[3], byrow=FALSE))
names(mcmc.samps) <- run$fit1$par_names
write.csv(mcmc.samps, file="mcmc_out/iterations.csv", row.names=FALSE)

# write summary file of NUTS/MCMC diagnostics
write.table(sum.dia, file="mcmc_out/table_convergence_diagnostics.csv", sep=",", row.names=FALSE)
saveRDS(run$fit1, file="mcmc_out/NUTS_fit.RDS")

# Launch Shiny App to check diagnostics online
# launch_shinyadmb(fit.1)

rm(list = ls(all.names = TRUE))

# plot_marginals(fit.1)
# plot_sampler_params(fit.1)
# plot_uncertainties(fit.1)

# pairs_admb(fit.1,pars=1:20)
