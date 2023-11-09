# Script for performing retrospective analyses on BASA
# Joshua Zahner | 05/17/2022 
#
# Creates a new subdirectory ('retrospectives/') and runs
# BASA on successivley fewer years of data (up to a provided
# number of 'peels'). Final biomass estimates and recrtuiment
# deviates are aggregated and analysed for bias patterns relative
# to the most recent (current year) full BASA run.
library(here)
library(doParallel)
library(tidyverse)
library(icesAdvice)
setwd(here::here())
source(file=file.path(here::here(), "functions/fun_read_dat.R"))
source(file=file.path(here::here(), "functions/fun_write_dat.R"))
source(file=file.path(here::here(), "functions/run_basa.R"))
source(file=file.path(here::here(), "plotting/compute_plot_products.R"))

n.peels <- 5
main.basa.directory <- paste0(here::here("model/"))

run.parallel <- FALSE
total.cores <- parallel::detectCores()
parallel.runs <- (total.cores-1) %/% 4      # BASA runs using 4 nodes

run.retropspective <- function(i, ...){
    new.dir.name <- paste0("basa-", i)
    if(!dir.exists(new.dir.name)){
        dir.create(new.dir.name)
    }
    setwd(new.dir.name)
    
    if(!file.exists("mcmc_out/PFRBiomass.csv")){
        print(paste("mcmc_out directory does not exist in", new.dir.name))
        file.symlink(paste0(main.basa.directory, "PWS_ASA.TPL"), ".")
        file.symlink(paste0(main.basa.directory, "PWS_ASA(par).ctl"), ".")
        file.symlink(paste0(main.basa.directory, "PWS_ASA(phases).ctl"), ".")
        file.symlink(paste0(main.basa.directory, "PWS_ASA(sim_settings).ctl"), ".")
        dat.files <- read.data.files(main.basa.directory)
        fun_write_dat(dat.files, i)
        run.basa(paste0("retrospectives/", new.dir.name))    
    }else{
        print("mcmc_out directory already exists")
    }
    setwd("..")
}

# Set up retrospectives directory
if(!dir.exists("retrospectives/")){
    dir.create("retrospectives/")
}
setwd("retrospectives/")

# Need to copy and then modify all of input files so they only reflect the data
# that was available in that year. Then actually run the assessment on all of
# the old datasets. Do this in parallel to speed things up.
if(run.parallel){
    cluster <- parallel::makeCluster(parallel.runs, type="FORK")
    registerDoParallel(cluster)
    foreach(i=1:n.peels) %dopar% {
        run.retropspective(i)
    }
    stopCluster(cluster)
}else{
    for(i in 1:5){
        run.retropspective(i=i)
    }
}
