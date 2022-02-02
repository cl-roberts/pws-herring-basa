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
library(data.table)
library(tidyverse)
library(adnuts)
library(snowfall)
library(rstan)
library(r4ss)

#################################################################
# Are you running this on a PC or a Mac
OS <- "MAC"

# BE SURE TO CHECK YOUR DIRECTORY
if(!exists("model")){
  template.files <- here::here("model")
  #template_files <- here::here("admb/vhs_sero_byage copy/")
  #template_files <- here::here("admb/ich_prev_byage/")
  #template_files <- here::here("admb/both_byage/")
}
function.dir <- here::here("functions/")
setwd(template.files)

source(file=paste0(function.dir, "data_reader.R"))
source(file=paste0(function.dir, "data_header_reader.R"))
source(file=paste0(function.dir, "read_admb_dat_files.R"))
source(file=paste0(function.dir, "calculate_ess.R"))
source(file=paste0(function.dir, "init_admb_params.R"))

system("admb -s PWS_ASA")

model.data <- read.admb.files()

nyr.fit <- model.data$nyr.fit # 42 for data up through 2021
nyr.tot <- model.data$nyr.tot # 42 for data up through 2021
nage <- 10                    # number of age classes (age 0 to 9+) (only care about 3-9+)

# Read in measured age comps
seine.age.comp <- model.data$seine.age.comp
spawn.age.comp <- model.data$spawn.age.comp
vhsv.age.comp <- model.data$vhsv.age.comp
ich.age.comp <- model.data$ich.age.comp

# Read in the actual sample sizes
# The following commands skips over lines of the file to read in tables 
# so change the number skipped if file is modified  
seine.samp.size <- read.table("agecomp_samp_sizes.txt", header=FALSE, skip=4,                 nrows=nyr.tot)
spawn.samp.size <- read.table("agecomp_samp_sizes.txt", header=FALSE, skip=4+1*(nyr.tot+1),   nrows=nyr.tot)
vhsv.samp.size  <- read.table("agecomp_samp_sizes.txt", header=FALSE, skip=4+2*(nyr.tot+1),   nrows=nyr.tot)
ich.samp.size   <- read.table("agecomp_samp_sizes.txt", header=FALSE, skip=4+3*(nyr.tot+1)+1, nrows=nyr.tot)


# Create empty matrices to fill estimated ESS and age comps
seine.ess.its <- matrix(0, nyr.tot, 1)  # Matrix to hold all iterations of the routine
seine.ess.its <- seine.samp.size        # fill in the first column with the recorded sample size 

spawn.ess.its <- matrix(0, nyr.tot, 1)
spawn.ess.its <- spawn.samp.size

vhsv.ess.its  <- matrix(0, nyr.tot, 1)
vhsv.ess.its  <- vhsv.samp.size
ich.ess.its   <- matrix(0, nyr.tot, 1)
ich.ess.its   <- ich.samp.size

# Change phases of the ESS in phases file to use the PWS_ASA(ESS_estimate)
phases <- readLines("PWS_ASA(ESS).ctl", -1)
phases[5] <- 1
writeLines(phases, "PWS_ASA(ESS).ctl")

# LOOP THROUGH AND ITERATIVELY CALCULATE ESS
convergence <- 0

seine.ess <- seine.samp.size
spawn.ess <- spawn.samp.size
vhsv.ess  <- vhsv.samp.size
ich.ess   <- ich.samp.size
its <- 1

age.comps <- list(seine=seine.age.comp, spawn=spawn.age.comp, vshv=vhsv.age.comp, ich=ich.age.comp)
start.ess <- list(seine=seine.ess, spawn=spawn.ess, vhvs=vhsv.ess, ich=ich.ess)
samp.size <- list(seine=seine.samp.size, spawn=spawn.samp.size, vhsv=vhsv.samp.size, ich=ich.samp.size)

# calc.ess <- calculate.ess(age.comps, start.ess, samp.size, nyr.fit)

for(i in 1:2){
  # Create "PWS_ASA(ESS_estimate).ctl" with sample sizes (the original sample sizes on the first iteration)
  write.table(
        rbind("# PWS age comp effective sample sizes", 
            "# Seine ESS",      seine.ess,  " ", 
            "# Spawn ESS",      spawn.ess,  " ", 
            "# VHS sero ESS",   vhsv.ess,   " ", 
            "# Ich prev ESS",   ich.ess
        ),
        file = "PWS_ASA(ESS_estimate).ctl", 
        append = F, 
        sep = " ",
        row.names=FALSE, 
        col.names=FALSE, 
        quote=F
    )
  
  # Compile and Run PWS_ASA
  if(OS=="MAC"){
    # system("./PWS_ASA -pinwrite -nohess")
    system("./PWS_ASA -pinwrite")
  }else if(OS=="PC"){
    shell("PWS_ASA  -pinwrite -nohess")
  }
  
  
  # Read in the estimated seine and spawner age comps
  seine.age.comp.est <- read.table("rep_out/SeAC_pd.rep",     header = FALSE) 
  spawn.age.comp.est <- read.table("rep_out/SpAC_pd.rep",     header = FALSE)
  vhsv.age.comp.est  <- read.table("rep_out/vhscomp_pd.rep",  header = FALSE)
  ich.age.comp.est   <- read.table("rep_out/ichcomp_pd.rep",  header = FALSE)
  
  # Calculate the ESS
  seine.ess <- rowSums(seine.age.comp.est*(1-seine.age.comp.est))/rowSums((seine.age.comp[1:nyr.fit, ]-seine.age.comp.est)^2)
  spawn.ess <- rowSums(spawn.age.comp.est*(1-spawn.age.comp.est))/rowSums((spawn.age.comp[1:nyr.fit, ]-spawn.age.comp.est)^2)
  vhsv.ess  <- rowSums(vhsv.age.comp.est*(1-vhsv.age.comp.est))/rowSums((vhsv.age.comp[1:nyr.fit, ]-vhsv.age.comp.est)^2)
  ich.ess   <- rowSums(ich.age.comp.est*(1-ich.age.comp.est))/rowSums((ich.age.comp[1:nyr.fit, ]-ich.age.comp.est)^2)
  
  # Remove the missing years of age comps
  seine.ess.rem <- seine.ess[!(seine.age.comp[1:nyr.fit, 1]==-9)]
  spawn.ess.rem <- spawn.ess[!(spawn.age.comp[1:nyr.fit, 1]==-9)]
  vhsv.ess.rem <- vhsv.ess[!(vhsv.age.comp[1:nyr.fit, 1]==-9)]
  ich.ess.rem <- ich.ess[!(ich.age.comp[1:nyr.fit, 1]==-9)]
  
  # Calculate the ratio of ESS to original sample sizes
  seine.ratio <- seine.ess.rem/seine.samp.size[1:nyr.fit, 1][seine.age.comp[1:nyr.fit, 1]!=-9]
  spawn.ratio <- spawn.ess.rem/spawn.samp.size[1:nyr.fit, 1][spawn.age.comp[1:nyr.fit, 1]!=-9]
  vhsv.ratio  <- vhsv.ess.rem/vhsv.samp.size[1:nyr.fit, 1][vhsv.age.comp[1:nyr.fit, 1]!=-9]
  ich.ratio   <- ich.ess.rem/ich.samp.size[1:nyr.fit, 1][ich.age.comp[1:nyr.fit, 1]!=-9]
  
  # Calculate the harmonic means (see Muradian et al. 2017 and Stewart & Hamel 2014)
  seine.hm <- 1/mean(1/seine.ratio)
  spawn.hm <- 1/mean(1/spawn.ratio)
  vhsv.hm <- 1/mean(1/vhsv.ratio)
  ich.hm <- 1/mean(1/ich.ratio)
  
  # Compare this harmonic mean to the previous using a convergence criteria (WHAT AM I CONVERGING!!!!)
  if(its==1) {
    convergence <- 0
    seine.hmS <- seine.hm
    spawn.hmS <- spawn.hm
    vhsv.hmS  <- vhsv.hm
    ich.hmS   <- ich.hm
  }else{
    seine.test <- abs(seine.hm - seine.hmS[its-1])/seine.hmS[its-1]*100
    spawn.test <- abs(spawn.hm - spawn.hmS[its-1])/spawn.hmS[its-1]*100
    vhsv.test  <- abs(vhsv.hm - vhsv.hmS[its-1])/vhsv.hmS[its-1]*100
    ich.test   <- abs(ich.hm - ich.hmS[its-1])/ich.hmS[its-1]*100

    convergence <- (seine.test<0.1 & spawn.test<0.1 & vhsv.test<0.1) # This criteria was arbitrarily chosen (0.1% change)
    seine.hmS <- rbind(seine.hmS, seine.hm)
    spawn.hmS <- rbind(spawn.hmS, spawn.hm) 
    vhsv.hmS  <- rbind(vhsv.hmS, vhsv.hm) 
    ich.hmS   <- rbind(ich.hmS, ich.hm) 
  }
  
  # Now multiply the harmonic mean by the sample size to get the new ESS 
  seine.ess <- round(seine.hm*seine.samp.size, 0)
  spawn.ess <- round(spawn.hm*spawn.samp.size, 0)
  vhsv.ess  <- round(vhsv.hm*vhsv.samp.size,   0)
  ich.ess   <- round(ich.hm*ich.samp.size,     0)
  
  # Use the average ESS for all years (each years obs weighted equally)
  # seine.ess[seine.ess>0] <- round(mean(seine.ess[seine.ess>0]), digits=0)
  # spawn.ess[spawn.ess>0] <- round(mean(spawn.ess[spawn.ess>0]), digits=0)
  
  # Denote the missing values
  seine.ess[(seine.age.comp[, 1] == -9), 1] <- -9
  spawn.ess[(spawn.age.comp[, 1] == -9), 1] <- -9
  vhsv.ess[(vhsv.age.comp[, 1]   == -9), 1] <- -9
  ich.ess[(ich.age.comp[, 1]     == -9), 1] <- -9
  
  vhsv.ess <- vhsv.samp.size
  ich.ess <- ich.samp.size
  
  # Fill in this iteration"s ESS
  seine.ess.its <- cbind(seine.ess.its, round(seine.ess, 0))
  spawn.ess.its <- cbind(spawn.ess.its, round(spawn.ess, 0))
  vhsv.ess.its <- cbind(vhsv.ess.its, round(vhsv.ess, 0))
  ich.ess.its <- cbind(ich.ess.its, round(ich.ess, 0))
  
  # Cease iterations if convergence hasn"t happened after so many...
  if(its==10) {
    break
  }
  its <- its+1
}

# Turn of the phases for the ESS calculation so it no longer recalculates ESS in future model runs
# Now write the converged ESS to a ctl file to be used for model runs
write.table(
  rbind(
    "# PWS age comp effective sample sizes", paste0("# (", date(), ")"), " ",
    "# Determines which ctl file for the age comps and ESS to use (1 uses ESS control to be iteratively estimated)", -1, " ",
    "# Seine ESS", seine.ess, " ", 
    "# Spawn ESS", spawn.ess, " ", 
    "# Sero ESS", vhsv.ess, " ", 
    "# Ich prev ESS", ich.ess
  ),
  file = "PWS_ASA(ESS).ctl", 
  append = F, 
  sep = " ",
  row.names=FALSE, 
  col.names=FALSE, 
  quote=F
)
  
######################################################
# Create reps x starting par vectors, and run NUTS
setwd(template.files)
reps <- 4
set.seed(8558)
seeds <- sample(1:1e4, size=reps)
#system("admb -s PWS_ASA")
system("./PWS_ASA -pinwrite -hbf 1")

inits <- init.admb.params(reps)

# ADMB command for running nuts
# PWS_ASA -nox -noest -nohess -maxfn 0 -nuts -mcmc 2000 -warmup 500 -chain 1 -mcseed 8682524 -max_treedepth 12 -adapt_delta 0.8 -adapt_mass -mcpin init.pin

# Pilot run to check 

if(dir.exists("mcmc_out")){
    system("rm mcmc_out/*.csv")
}

start.time <- Sys.time()
fit.1 <- sample_nuts(model='./PWS_ASA',path=template.files,
                        #iter=500,
                        #warmup=700,
                        #warmup=100,
                        duration = 5,
                        init=inits, seeds=seeds, chains=reps,cores=reps,
                        mceval=TRUE,
                        control=list(
                            adapt_delta=0.9,
                            #max_treedepth=16,
                            metric="mle"
                        )
                    )
end.time <- Sys.time()
end.time - start.time

# Extracts NUTS stats (energy, leapfrog transitions,etc)
mon <- monitor(fit.1$samples, warmup=fit.1$warmup, print=FALSE)
x <- extract_sampler_params(fit.1)

# Quick check for divergences & Gelman-Ruben statistic
n.divergences <- sum(x$divergent__)/nrow(x)
r.hat <- max(mon[, "Rhat"])<=1.1

n.divergences
r.hat

# Write summary of parameter posteriors (medians, percentiles, etc)
write.csv(mon, file="mcmc_out/table_par_posterior_summary.csv")

# Write all MCMC samples of the parameters
mcmc.samps <- data.frame(matrix(fit.1$samples, nrow=dim(fit.1$samples)[3], byrow=TRUE))
names(mcmc.samps) <- fit.1$par_names
write.csv(mcmc.samps, file="mcmc_out/iterations.csv", row.names=FALSE)

## Examine the slowest mixing parameters
# slow <- names(sort(mon[,"n_eff"]))[1:8]
# pairs_admb(fit=fit, pars=slow)

# Create summary file of NUTS/MCMC diagnostics
sum.dia <- data.frame(divergences.from.extract.function=sum(x$divergent__)/nrow(x),
                      min.ESS=min(mon[, "n_eff"]),
                      which.min.ESS=names(which.min(mon[, "n_eff"])),
                      max.Rhat=max(mon[, "Rhat"]),
                      which.max.Rhat=names(which.max(mon[, "Rhat"])),
                      time.elapsed=end.time-start.time)
write.table(sum.dia, file="mcmc_out/table_convergence_diagnostics.csv", sep=",", row.names=FALSE)
saveRDS(fit.1, file="mcmc_out/NUTS_fit.RDS")

# Launch Shiny App to check diagnostics online
# launch_shinyadmb(fit.1)

rm(list = ls(all.names = TRUE))

# plot_marginals(fit.1)
# plot_sampler_params(fit.1)
# plot_uncertainties(fit.1)

# pairs_admb(fit.1,pars=1:20)
