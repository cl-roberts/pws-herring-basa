library(here)
library(tidyverse)
source(paste0(here::here("functions/"), "fun_read_dat.R"))



compute.exploit.rate <- function(model.dir, nyr, years){
    exploit.rate <- read.exploit.rates(model.dir, nyr)

    exploit.rate.df <- as_tibble(exploit.rate) %>%
                    pivot_longer(everything(), names_to="year", values_to="exploit") %>%
                    group_by(year) %>%
                    median_qi(exploit, .width=c(0.95)) %>%
                    print(n=10)

    exploit.zeros <- exploit.rate.df[exploit.rate.df$exploit == 0, ]

    return(listN(exploit.rate.df, exploit.zeros))
}

compute.biomass.traj <- function(model.dir, nyr, years){
    biomass <- read.biomass.estimates(model.dir, nyr)

    biomass.df <- as_tibble(biomass) %>%
                            pivot_longer(everything(), names_to="year", values_to="biomass") %>%
                            group_by(year) %>%
                            median_qi(biomass, .width=c(0.10, 0.50, 0.95)) %>%
                            print(n=10)

    prob.below.threshhold <- apply(biomass, 2, function(x) sum(x < 19958))/nrow(biomass)

    biomass.df$prob <- as.vector(rep(prob.below.threshhold, 3))

    return(biomass.df)
}

compute.pfrb.posterior <- function(model.dir, nyr, years){
    biomass <- read.table(paste0(model.dir, "mcmc_out/PFRBiomass.csv"), header = FALSE, sep = ",", dec=".")[,nyr]

    biomass.df <- data.frame(biomass=biomass)
    prob.below.threshold <- round(sum(biomass < 20000)/length(biomass), 2)

    biomass.quants <- round(as.vector(apply(as.matrix(biomass), 2, quantile, c(0.025, 0.5, 0.975)))/1000, 2)

    return(listN(biomass.df, biomass.quants, prob.below.threshold))
}

compute.recruitment <- function(model.dir, nyr, years){
    age.3.recruits <- read.table(paste0(model.dir, "mcmc_out/Age3.csv"), header = FALSE, sep = ",", dec=".")[,1:nyr]
    colnames(age.3.recruits) <- years

    age.3.recruits.df <- as_tibble(age.3.recruits) %>%
                            pivot_longer(everything(), names_to="year", values_to="recruits") %>%
                            group_by(year) %>%
                            median_qi(recruits, .width=c(0.50, 0.95)) %>%
                            print(n=10)

    return(age.3.recruits.df)
}

compute.catch.biomass <- function(model.dir, nyr, years){

    data <- read.data.files(model.dir)

    weight.at.age <- data$waa[1:nyr,]

    fb.nya <- data$foodbait_catch[1:nyr,]
    pound.nya <- data$pound_catch[1:nyr,]
    gillnet.nya <- data$gillnet_catch[1:nyr,]
    seine.yield  <- data$seine_yield[1:nyr]
    
    fb.nya.biomass <- weight.at.age * fb.nya 
    pound.nya.biomass <- weight.at.age * pound.nya
    gillnet.nya.biomass <- weight.at.age * gillnet.nya 

    fb.biomass.annual       <- rowSums(fb.nya.biomass) # Now sum the biomass over all age classes for each year
    pound.biomass.annual    <- rowSums(pound.nya.biomass)
    gillnet.biomass.annual  <- rowSums(gillnet.nya.biomass)

    fb.biomass.annual       <- replace(fb.biomass.annual, fb.biomass.annual == 0, NA)
    pound.biomass.annual    <- replace(pound.biomass.annual, pound.biomass.annual == 0, NA)
    gillnet.biomass.annual  <- replace(gillnet.biomass.annual, gillnet.biomass.annual == 0, NA)
    seine.yield             <- replace(seine.yield, seine.yield == 0, NA)
    
    # Matrix of catches by gear type in mt
    total.catch <- cbind(fb.biomass.annual, pound.biomass.annual, gillnet.biomass.annual, seine.yield)
    total.catch[is.na(total.catch)] <- 0 
    total.catch.biomass <- rowSums(total.catch) # total catches by year in mt

    return(total.catch.biomass)

}