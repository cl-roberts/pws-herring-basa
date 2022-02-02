# results_outputs_for_management.r
# Created by Melissa Muradian
# Edited by John Trochta
# Date updated:  04/07/2021
# Summary:
# This script outputs plots of posteriors of select derived quantities from the Bayesian ASA
# including the recruitment trajectory, spawning biomass trajectory, exploitation rates by year, and the estimated spawning biomass in the most recent year
###################################

model.path <- here::here("model/")
rfunctions.path <- here::here("functions/")
source(file=here::here("functions/data_reader.R"))
source(file=here::here("functions/data_header_reader.R"))

# Years of the model
years<-1980:2021 # CHANGE with each update
nYr<-length(years)+1
n.burn <- 1 # This is the burn-in period - be sure to adjust accordingly with you specified chain length and thinning rate

# load("~/PWS_herring/Model/catches_data.RData")
read.in.files <- function(model.path,n.burn,nYr){
    ############### Calculate F (fishing mortality)using catch data and estimated biomass ########
    # Reads in the catches from different fisheries ("catches_data.RData")
    # The variables are named as follows:
    # a:  Food and Bait catches (millions of fish), where first column is year, X3 is age 3's, X4 is age 4's, and so on (the following columns)
    # b:  Pound fishery catches (millions of fish), where first column is year, X3 is age 3's, X4 is age 4's, and so on (the following columns)
    # c:  Gill Net catches (millions of fish), where first column is year, X3 is age 3's, X4 is age 4's, and so on (the following columns)
    # seineYield: Purse Seine catches in metric tonnes
    # weight: Weight-at-age of purse seine catches, with first column as year, X3 is age 3's, X4 is age 4's, and so on
    
    # Store the file names from which data is available
    filename <- vector(length=2)
    filename[1] = "PWS_ASA.dat"
    filename[2] = "PWS_ASA(ESS).ctl"
    #filename[3]="PWS_ASA(age_comps).ctl"
    #filename[4]="PWS_ASA(surveys).ctl"
    
    PWS_ASA.dat <- data.reader(filename=filename[1])
    PWS_ASA.dat <- PWS_ASA.dat[-2]
    
    foodbait <- rowSums(PWS_ASA.dat[[7]]) # Removes the first column (years) from summations
    pound <- rowSums(PWS_ASA.dat[[5]])
    gillnet <- rowSums(PWS_ASA.dat[[8]])
    seineYield <- PWS_ASA.dat[[9]]
    seineYield <- replace(seineYield, seineYield == 0, NA) # If an element of seineYield is 0, replace with NA
    
    weight <- PWS_ASA.dat[[3]]
    fb.nya <- PWS_ASA.dat[[7]]
    fb.nya.biomass <- weight*fb.nya # Now calculate the tonnage (e.g. biomass) of each age class
    pound.nya <- PWS_ASA.dat[[5]]
    pound.nya.biomass <- weight*pound.nya
    gillnet.nya <- PWS_ASA.dat[[8]]
    gillnet.nya.biomass <- weight*gillnet.nya 

    fb.biomass.annual <- rowSums(fb.nya.biomass) # Now sum the biomass over all age classes for each year
    fb.biomass.annual <- replace(fb.biomass.annual, fb.biomass.annual == 0, NA)

    pound.biomass.annual <- rowSums(pound.nya.biomass)
    pound.biomass.annual <- replace(pound.biomass.annual, pound.biomass.annual == 0, NA)

    gillnet.biomass.annual <- rowSums(gillnet.nya.biomass)
    gillnet.biomass.annual <- replace(gillnet.biomass.annual, gillnet.biomass.annual == 0, NA)
    
    # Matrix of catches by gear type in mt
    total.catch <- cbind(fb.biomass.annual, pound.biomass.annual, gillnet.biomass.annual, seineYield)
    total.catch[is.na(total.catch)] <- 0 
    total.catch.mass <- rowSums(total.catch) # total catches by year in mt
    
    ###########################################################################  
    # Read in the MCMC draws of recruitment (millions of age 3's) and 
    # pre-fishery spawning biomass (metric tons)
    
    setwd(model.path)
    Age3 <- read.table(paste0(model.path, "mcmc_out/Age3.csv"), header = FALSE, sep = ",", dec=".")
    Age3 <- Age3[-c(1:n.burn), ]

    PFRB <- read.table(paste0(model.path, "mcmc_out/PFRBiomass.csv"), header = FALSE, sep = ",", dec=".") 
    PFRB <- PFRB[-c(1:n.burn), ]

    #quantile(PFRB[,nYr]/1000, probs=c(.025, .5, .975)) # Calculate the 95% credibility interval of PFRB converted to 1000's metric tons
    final.PFRB <- PFRB[, nYr] # year 2018
    #allIts <- length(final.PFRB) # number of thinned iterations after burn-in
    
    exploit.rate.mat<-total.catch.mass[1:(nYr-1)]/t(PFRB[,1:(nYr-1)]) # Calculate the exploitation rate
    exploit.rate.mat<-rbind(exploit.rate.mat,rep(0,ncol(exploit.rate.mat)))
    
    fMedian<-apply(exploit.rate.mat, 1, median) # And its median
    fMedian<-replace(fMedian, fMedian == 0, NA)
    
    # lower .95 int is row 1, upper is row 2
    exploit.rate.ci <- apply(exploit.rate.mat, 1, quantile, probs=c(.025, .975)) # And its credibility interval

    return(
      list(
        Age3 = Age3,
        PFRB = PFRB,
        final.PFRB = final.PFRB,
        exploit.rate.mat = exploit.rate.mat,
        f.median = fMedian,
        exploit.rate.ci = exploit.rate.ci,
        catch = total.catch.mass
      )
    )

}
 
## PLOTS ===================================================

## Plot Age-3 Rec
source(file=here::here("plotting/plot_recruitment_posterior.R"))
## Plot histogram of projected final year biomass in mt
# 95% CI with median
source(file=here::here("plotting/plot_pfrb_posterior.R"))
## Plot biomass
source(file=here::here("plotting/plot_ssb_trajectory.R"))
## Exploitation Rate plot - from Catches.R
source(file=here::here("plotting/plot_exploitation_rate.R"))

###############################################################
## Plot management quantities for a single model
setwd(model.path)
data <- read.in.files(model.path, n.burn, nYr)

#pdf(here::here("figures/outputs-for-management.pdf"), height=8, width=12)
#par(mfcol=c(2, 2), mar=c(4, 5, 2, 2), oma=c(1, 1, 0, 3))

ylabel <- TRUE
labels.y <- TRUE
ylabel.2 <- TRUE
labels.y2 <- TRUE
xlabel <- TRUE
labels.x <- TRUE
projected.recruit <- apply(data$Age3, 1, FUN=function(x) sum(x[(length(x)-9):length(x)])/10)
#Age3 <- Age3[,1:(nYr-1)]
#Age3 <- cbind(Age3,projected.recruit)
ylim <- 2000
myMGP=c(3, 0.8, 0) 

pdf(here::here("figures/recruitment_and_ssb.pdf"), height=4, width=12)
par(mfcol=c(1, 2), mar=c(4, 5, 2, 2), oma=c(1, 1, 0, 3))
rec.posterior.figure <- rec_posterior(years=years, nYr-1, data$Age3, myMGP=myMGP, modelPath=rfunctions.path, 
                                      ylim=ylim, ylabel, xlabel, labels.x, labels.y)

ylim <- 180
ssb.traj.figure <- ssb_trajectory(data$PFRB, years=c(years, tail(years,1)), nYr-1, myMGP=myMGP, modelPath=rfunctions.path, 
                                  ylim=ylim, ylabel, xlabel, labels.x, labels.y, ylabel.2, labels.y2)
dev.off()

pdf(here::here("figures/exploitation_rate.pdf"), height=4, width=6)
par(mfcol=c(1, 1), mar=c(4, 5, 2, 2), oma=c(1, 1, 0, 3))
exploit.rate.figure <- exploit_rate(data$f.median[1:nYr-1], data$exploit.rate.ci[, 1:nYr-1], years=c(years, tail(years,1)), nYr-1, 
                                    myMGP=myMGP, ylabel, xlabel, labels.x, labels.y)
dev.off()


xlabel <- TRUE
xlim <- 40
ylim <- 0.13
labels.y <- TRUE

pdf(here::here("figures/post_fishery_biomass_posterior.pdf"), height=4, width=6)
par(mfcol=c(1, 1), mar=c(4, 5, 2, 2), oma=c(1, 1, 0, 3))
pfrb.posterior.figure <- pfrb_posterior(data$PFRB, years, nYr-1, myMGP=myMGP, xlim=xlim, ylim=ylim, ylabel, xlabel, labels.x, labels.y)
dev.off()

###############################################################
## Save csv of values displayed in table above (with 95% credibility intervals) 
round.df <- function(x, digits) {
  # round all numeric variables
  # x: data frame 
  # digits: number of digits to round
  numeric.columns <- sapply(x, mode) == "numeric"
  x[numeric.columns] <-  round(x[numeric.columns], digits)
  return(x)
}
recruitment.est <-  as.matrix(t(apply(data$Age3[, 1:(nYr-1)], 2, quantile,
                                     probs=c(0.5, 0.025, 0.975), name=F, na.rm=T)))
ssb.est <- as.matrix(t(apply(data$PFRB[, 1:(nYr-1)], 2, quantile,
                             probs=c(0.5, 0.025, 0.975), name=F, na.rm=T)))
exploit.rate.est <- as.matrix(t(apply(data$exploit.rate.mat[1:(nYr-1), ], 1, quantile,
                            probs=c(0.5, 0.025, 0.975), name=F, na.rm=T)))
final.table <- round.df(data.frame(years, recruitment.est, ssb.est/1000, data$catch, exploit.rate.est, round.df(ssb.traj.figure, 2)), 2)
names(final.table) <- c("Years",
                        "Median Age 3 (in millions)",
                        "Lower 95th Age 3 (in millions)",
                        "Upper 95th Age 3 (in millions)",
                        "Median Pre-fishery biomass (in 1000s metric tons)",
                        "Lower 95th Biomass (in 1000s metric tons)",
                        "Upper 95th Biomass (in 1000s metric tons)",
                        "Median Exploitation rate",
                        "Lower 95th ER",
                        "Upper 95th ER",
                        "Probability B<20K")
write.csv(final.table, here::here("data_outputs/outputs-for-management.csv"), row.names=FALSE)
