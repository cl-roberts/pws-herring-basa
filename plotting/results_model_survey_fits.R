# model_fits_to_data.r
# Created by Melissa Muradian
# Edited by John Trochta
# Date updated:  04/07/2021
# Summary:
#   Plot Bayesian ASA model fits to the survey data

###############################################################################
source(file=here::here("functions/data_reader.R"))
source(file=here::here("functions/data_header_reader.R"))

# Calculates 95% confidence intervals from survey CV (from Buckland 1992 in References)
calc.buck.cv <- function(cv=cv) {
    buck.c <- NULL
    cv <- cv
    buck.c <- exp(1.96*sqrt(log(1+(cv^2))))
    return(buck.c)
}
# Function for creating shaded polygons as credibility intervals around model estimates
add.polygon2 <- function(x, y, z, alpha.level, alpha.min=0, alpha.max=1){
    ## Pass this a series of years (x) and a matrix of trajectories (df) and it
    ## adds a polygon to the current plot with whose area contains 1-alpha.level
    alpha.level <- sort(alpha.level)
    for(i in 1:length(alpha.level)){
        alpha <- alpha.level[i]
        alpha.col <- alpha.min+alpha*(alpha.max-alpha.min)
        col.poly <- rgb(1-alpha.col, 1-alpha.col, 1-alpha.col, alpha=1)
        quantiles.temp <-  as.matrix(t(apply(z, 2, quantile, probs=c(alpha/2, 1-alpha/2), name=F, na.rm=T)))
        # Next part accounts for missing values - plots separate polygons over observed sections of time series
        if(any(is.na(quantiles.temp[, 2]))){
        indices <- which(!is.na(quantiles.temp[, 2]))
        indices_start <- indices[c(TRUE, diff(indices)!=1)]
        indices_end <- indices[c(diff(indices)!=1, TRUE)]
        for(j in 1:length(indices_start)){
            x_1 <- x[indices_start[j]:indices_end[j]]
            x_1 <- c(head(x_1, 1)-0.5, x_1, tail(x_1, 1)+0.5)
            quantiles_1 <- quantiles.temp[indices_start[j]:indices_end[j], 1]
            quantiles_1 <- c(head(quantiles_1, 1), quantiles_1, tail(quantiles_1, 1))
            quantiles_2 <- quantiles.temp[indices_start[j]:indices_end[j], 2]
            quantiles_2 <- c(head(quantiles_2, 1), quantiles_2, tail(quantiles_2, 1))
            polygon(x=c(x_1, rev(x_1)), y=c(quantiles_1, rev(quantiles_2)),
                    col=col.poly, border=NA)
        }
        }else{
            polygon(x=c(x, rev(x)), y=c(quantiles.temp[, 1], rev(quantiles.temp[, 2])), col=col.poly, border=NA)
        }
    }
}
# Read in both input to and output files from BASA
read.in.files <- function(modelPath,nburn,sYr,nYr){
    # Store the file names from which data is available
    filename <- vector(length=2)
    filename[1]="PWS_ASA.dat"
    filename[2]="PWS_ASA(ESS).ctl"
    
    PWS_ASA.dat <- data.reader(filename=filename[1])
    PWS_ASA.dat <- PWS_ASA.dat[-2]
    
    #  ADFGHydData:     Hydroacoustic estimates from ADF&G of herring biomass, set #2
    #  PWSSCHydData:    Hydroacoustic estimates from PWSSC of herring biomass
    #  mdm:             Mile-days of milt index
    #  egg:             Egg deposition indices
    #  egg.se:          Standard error for the egg deposition survey
    #  PWSSC.se:        Standard error for the PWSSC hydroacoustic survey
    #  age1schools:     Aerial surveys of age-1 schools
    mdm <- PWS_ASA.dat[[11]]
    egg <- PWS_ASA.dat[[12]]
    egg.se <- PWS_ASA.dat[[13]]
    ADFGHydData <- PWS_ASA.dat[[15]]
    PWSSCHydData <- PWS_ASA.dat[[17]]
    PWSSC.se <- PWS_ASA.dat[[18]]
    age1schools <- PWS_ASA.dat[[21]]
    
    mdm[mdm==-9] <- NA
    egg[egg==-9] <- NA
    egg.se[egg.se==-9] <- NA
    ADFGHydData[ADFGHydData==-9] <- NA
    PWSSCHydData[PWSSCHydData==-9] <- NA
    PWSSC.se[PWSSC.se==-9] <- NA
    age1schools[age1schools==-9] <- NA
    
    # The following reads in the model estimates (draws from MCMC chain of N length) for the different survey indices
    HYD_ADFG      <- read.table(paste0(modelPath, "mcmc_out/HYD_ADFG.csv"),    header = FALSE, sep = ",", dec=".")  
    HYD_PWSSC     <- read.table(paste0(modelPath, "mcmc_out/HYD_PWSSC.csv"),   header = FALSE, sep = ",", dec=".")
    MDM           <- read.table(paste0(modelPath, "mcmc_out/MDM.csv"),         header = FALSE, sep = ",", dec=".")
    EGG           <- read.table(paste0(modelPath, "mcmc_out/EGG.csv"),         header = FALSE, sep = ",", dec=".")
    VARSReport    <- read.table(paste0(modelPath, "mcmc_out/VarsReport.csv"),  header = FALSE, sep = ",", dec=".")
    AGE1SCHOOlS   <- read.table(paste0(modelPath, "mcmc_out/juv_schools.csv"), header = FALSE, sep = ",", dec=".")
    age1.overdisp <- read.table(paste0(modelPath, "mcmc_out/iterations.csv"),  header = FALSE, sep = ",", dec=".")
    
    # Get rid of the burn-in draws
    HYD_ADFG        <- HYD_ADFG[-c(1:nburn), ] 
    HYD_PWSSC       <- HYD_PWSSC[-c(1:nburn), ]
    MDM             <- MDM[-c(1:nburn), ]
    EGG             <- EGG[-c(1:nburn), ]
    VARSReport      <- VARSReport[-c(1:nburn), ]
    AGE1SCHOOlS     <- AGE1SCHOOlS[-c(1:nburn), ]
    age1.overdisp   <- age1.overdisp[-c(1:nburn), 28]

    # Replace years with no data (read in as zero to these output files) with NA so that they do not appear on the plot
    EGG[EGG==0] <- NA 
    AGE1SCHOOlS[, is.na(age1schools)] <- NA
    
    # Order is m_add, egg_add, hydADFG_add, hydPWSSC_add.
    # m_add:        additional error on the mile-days of milt survey observations (estimated within model)
    # egg_add:      additional error on the egg deposition survey observations (fixed within model)
    # hydADFG_add:  additional error on the ADF&G hydroacoustic survey observations (estimated within model)
    # hydPWSSC_add: additional error on the PWSSC hydroacoustic survey observations (estimated within model)
    variance <- VARSReport
    names(variance) <- c("m_add", "egg_add", "hydADFG_add", "hydPWSSC_add")
    # Calculate CV's for each survey type, using both estimated and provided variances (where appropriate)
    mediansAdd  <- apply(variance, 2, median); mediansAdd
    mdm.cv      <- rep(mediansAdd[1], nYr)
    egg.cv      <- as.vector(sqrt(egg.se^2 + mediansAdd[2]^2))
    ADFG.cv     <- rep(mediansAdd[3], nYr) 
    PWSSC.cv    <- sqrt((PWSSC.se^2) + (mediansAdd[4]^2))
    
    return(
        list(
            mdm              = mdm,
            egg              = egg,
            ADFG.hydro       = ADFGHydData,
            PWSSC.hydro      = PWSSCHydData,
            age1.schools     = age1schools,
            egg.se           = egg.se,
            PWSSC.se         = PWSSC.se,

            ADFG.hydro.est   = HYD_ADFG,
            PWSSC.hydro.est  = HYD_PWSSC,
            mdm.est          = MDM,
            egg.est          = EGG,
            age1.schools.est = AGE1SCHOOlS,

            variance         = variance,
            mdm.cv           = mdm.cv,
            egg.cv           = egg.cv,
            ADFG.cv          = ADFG.cv,
            PWSSC.cv         = PWSSC.cv,
            age1.overdisp    = age1.overdisp
        )
    )
}
# Computes posterior predictive intervals for log normal
post.pred.lognorm <- function(model.pred, sd, obs, cv){
    PP <- matrix(NA,nrow=NROW(model.pred),ncol=NCOL(model.pred))
    for(i in 1:NROW(model.pred)){
        fits <- as.numeric(model.pred[i,]) # Take the true value of the survey estimate
        errors <- as.numeric(sd[i]*fits) # Take the error based on the true value of the survey estimate
        
        # Need to reparameterize to calculate the log-normal mean and sd for the survey data
        # From: https://msalganik.wordpress.com/2017/01/21/making-sense-of-the-rlnorm-function-in-r/
        E <- log(fits^2/sqrt(errors^2+fits^2))
        SD <- sqrt(log(1+(errors^2/fits^2)))
        PP[i,] <- rlnorm(n=length(fits),meanlog=E,sdlog=SD)
        # MDM.PP[i,] <- rnorm(n=length(fits),mean=E,sd=SD) # Generate a random variable based on the expected value (predictive) using the assumed 
    }
    return(PP)
}
# Computes posterior predictive intervals for a negative binomial on juvenile school counts
post.pred.negbin <- function(PP, sd, obs){
    PP <- matrix(NA,nrow=NROW(model.pred),ncol=NCOL(model.pred))
    for(i in 1:NROW(model.pred)){
        fits <- as.numeric(model.pred[i,]) # Take the true value of the survey estimate
        dispersion <- as.numeric(sd[i]) # Take the error based on the true value of the survey estimate
        
        E <- fits
        SD <- exp(dispersion)
        PP[i,] <- rnbinom(n=length(fits),size=SD,prob=SD/(SD+E))
        # MDM.PP[i,] <- rnorm(n=length(fits),mean=E,sd=SD) # Generate a random variable based on the expected value (predictive) using the assumed 
    }
    return(PP)
}

# Plots log normal posterior predictive intervals
pp.plotter.ln <- function(PP, sd, obs, cv, sYr, nYr, dataYears, years, ylabel, xlabel, myMGP, labels.x, labels.y){
    # First set up the axes
    ymax.1 <- max(t(obs[sYr:nYr,1])*calc.buck.cv(cv=cv)[sYr:nYr],na.rm=T)
    ymax.2 <- max(apply(PP[,sYr:nYr], 2, quantile,probs=0.99,name=F, na.rm=T),na.rm=T)
    ylimit <- signif(max(c(ymax.1,ymax.2),na.rm=T),1)
    plot(dataYears[sYr:nYr], rep(NA,(nYr-sYr+1)), col=4, ylab='', xlab='',axes=F, xaxs="i", yaxs="i",ylim=c(0,ylimit),
        pch=16, cex.main=1.6, cex.axis=2, xlim=c(years[sYr],years[nYr]+0.3), las=1,xpd=NA, main=ylabel, cex=2.5)
    # Then calculate posterior predictive intervals, specificied at the inner 5% and 95% intervals around the median
    add.polygon2(x=dataYears[sYr:nYr], z=PP[,sYr:nYr], alpha.level=c(.05,.5,.95), alpha.min=.2, alpha.max=.8)
    # Plot the survey observations
    points(dataYears[sYr:nYr], t(obs[sYr:nYr,1]), pch=16,xpd=NA) 
    lines(dataYears[sYr:nYr], t(obs[sYr:nYr,1]), type="l",xpd=NA)
    # Calculates confidence intervals around observations based on the survey CV
    segments(x0=dataYears[sYr:nYr],y0=(t(obs[sYr:nYr,1])/calc.buck.cv(cv=cv)[sYr:nYr]),x1=dataYears[sYr:nYr],y1=(t(obs[sYr:nYr,1])*calc.buck.cv(cv=cv)[sYr:nYr]),lwd=3,xpd=NA)
    # Add axes labels and ticks
    axis(1, at=years, labels=F, tcl=-.15, cex.axis=1.4)
    axis(1,at= seq(years[sYr], years[nYr], 5), labels=labels.x, tcl=-.35, lwd.ticks=2, mgp=myMGP, cex.axis=1.4)
    #axis(2, at=seq(0,500,100), labels=T, las=1, mgp=myMGP)
    axis(2, at=seq(0,ylimit,length.out=5), labels=labels.y, las=1, mgp=myMGP, cex.axis=1.4)
    
    #text(x=max(dataYears)-1,y=ylimit*0.9,labels=ylabel,cex=2.5,adj=c(1,1))
    
    # Plot the posterior predictive intervals, specificied at the inner 5% and 95% intervals around the median
    # add.polygon2 was not used because only a select number of years had observations
    # points(dataYears[sYr:nYr], apply(EGG.PP, 2, median), pch=16, col="grey60", cex=2)
    # segments(x0=dataYears[sYr:nYr],y0=apply(EGG.PP, 2, quantile, probs=0.025, na.rm=T),x1=dataYears[sYr:nYr],
    #          y1=apply(EGG.PP, 2, quantile, probs=0.975, na.rm=T),col="grey60",lwd=8)
}
# Plots posterior predictive intervals for a negative binomial
pp.plotter.nb <- function(PP, sd, obs, cv, sYr, nYr, dataYears, years, ylabel, xlabel, myMGP, labels.x, labels.y){
    # First set up the axes
    # ymax.1 <- max(t(obs[sYr:nYr,1])*calc.buck.cv(cv=cv)[sYr:nYr],na.rm=T)
    ymax.2 <- max(apply(PP[,sYr:nYr], 2, quantile,probs=0.99,name=F, na.rm=T),na.rm=T)
    #ylimit <- signif(ymax.2,1)
    ylimit <- 20000
    plot(dataYears[sYr:nYr], rep(NA,(nYr-sYr+1)), col=4, ylab='', xlab=xlabel,axes=F, xaxs="i", yaxs="i",ylim=c(0,ylimit),
        pch=16, cex.main=1.6, cex.lab=1.8, xlim=c(years[sYr],years[nYr]+0.3), las=1,xpd=NA, main=ylabel, cex=2.5)
    # Then calculate posterior predictive intervals, specificied at the inner 5% and 95% intervals around the median
    add.polygon2(x=dataYears[sYr:nYr], z=PP[,sYr:nYr], alpha.level=c(.05,.5,.95), alpha.min=.2, alpha.max=.8)
    # Plot the survey observations
    points(dataYears[sYr:nYr], t(obs[sYr:nYr,1]), pch=16,xpd=NA) 
    lines(dataYears[sYr:nYr], t(obs[sYr:nYr,1]), type="l",xpd=NA)
        # Add axes labels and ticks
    axis(1, at=years, labels=F, tcl=-.15, cex.axis=1.5)
    axis(1, at= seq(years[sYr], years[nYr], 5), labels=labels.x, tcl=-.35, lwd.ticks=2, mgp=myMGP, cex.axis=1.5)
    #axis(2, at=seq(0,500,100), labels=T, las=1, mgp=myMGP)
    axis(2, at=seq(0,ylimit,length.out=5), labels=labels.y, las=1, mgp=myMGP, cex.axis=1.5)
    #text(x=min(dataYears)+5,y=ylimit*0.9,labels=ylabel,cex=2.5,adj=c(0,1))
}

###################################
library(dplyr)

years <- 1980:2021 # Model fitting year - CHANGE!
dataYears <- 1980:2021 # Years the data spans - CHANGE#!
nYr <- length(dataYears) # Number of years in current data set
sYr <- 1
#sYr <- 29

# Defines the burn-in the number of initial draws to remove
# Remnant from old code using Metropolis Hastings MCMC (& not year updated here...)
nburn <- 1

modelPath<-here::here("model/")
setwd(modelPath) 

###############################################################################
# Plot the Mile-days of Milt, Egg deposition, ADF&G acoustic, and the PWSSC acoustic biomass predictive posteriors
# Predictive posteriors are model fits (estimate in a given MCMC draw) generated with random error 

data <- read.in.files(modelPath,nburn,sYr,nYr)

#pdf(paste0(last(years)," Model herring estimates (All years).pdf"),height=7.5,width=4,family="Times")
pdf(here::here("figures/predicted-survey-values.pdf"),height=10,width=12)

# fill <- c(1:4,5,5)
# footPrint <- matrix(fill, nrow=3, byrow=T) ## an 18 by 8 matrix
# layout(footPrint)

par(mfrow=c(3,2), mar=c(5,7,2,2))
myMGP <- c(3,1,0) # Plot formatting spcification
X1 <- list(data$mdm.est, data$egg.est, data$ADFG.hydro.est/1000, data$PWSSC.hydro.est/1000, data$age1.schools.est)
X2 <- list(data$variance$m_add, data$variance$egg_add, data$variance$hydADFG_add, data$variance$hydPWSSC_add, data$age1.overdisp)
X3 <- list(data$mdm.cv, data$egg.cv, data$ADFG.cv, data$PWSSC.cv, data$age1.overdisp)
# X4 <- c("Mile-days\nmilt",
#         "Eggs\ndeposition\n(trillions)",
#         "ADF&G\nacoustic\nbiomass\n(1000s tons)",
#         "PWSSC\nacoustic\nbiomass\n(1000s tons)",
#         "Age 1 schools")
X4 = c("Mile-days of Milt", "Egg Deposition (trillions)", "ADF&G Acoustic Biomass (1000s tons)", "PWSSC Acoustic Biomass (1000s tons)", "Age 1 Schools")
X5 <- list(data$mdm, data$egg, data$ADFG.hydro/1000, data$PWSSC.hydro/1000, data$age1.schools)
post.pred <- list()
for(i in 1:(length(X1))){
  model.pred <- X1[[i]]
  sd <- X2[[i]]
  cv <- X3[[i]]
  ylabel <- X4[i]
  #xlabel <- NA
  xlabel <- "Year"
  labels.x <- TRUE
  labels.y <- TRUE
  obs <- X5[[i]]
  if(i==5){
    post.pred[[i]] <- post.pred.negbin(PP,sd,obs)
    pp.plotter.nb(PP=post.pred[[i]],sd,obs,cv,sYr,nYr,dataYears,years,ylabel,xlabel,myMGP,labels.x,labels.y)
  }else{
    post.pred[[i]] <- post.pred.lognorm(model.pred,sd,obs,cv)
    pp.plotter.ln(PP=post.pred[[i]],sd,obs,cv,sYr,nYr,dataYears,years,ylabel,xlabel,myMGP,labels.x,labels.y)
  }
}

dev.off()

###############################################################
## Save csv of values displayed in table above (with 95% posterior predictive intervals) 
round_df <- function(x, digits) {
  # round all numeric variables
  # x: data frame 
  # digits: number of digits to round
  numeric_columns <- sapply(x, mode) == 'numeric'
  x[numeric_columns] <-  round(x[numeric_columns], digits)
  x
}


sink(here::here("data_outputs/predicted-survey-values.csv"))
cat('Milt (mile-days) posterior predictions\n')
write.table(cbind(years,t(round_df(apply(post.pred[[1]],2,quantile,probs=c(0.5,0.025,0.975),na.rm=TRUE),2))),row.names=FALSE,col.names=c('Year','Median Posterior Prediction','Lower 95th','Upper 95th'),sep=",")
cat('\nEggs deposited (trillions) posterior predictions\n')
write.table(cbind(years,t(round_df(apply(post.pred[[2]],2,quantile,probs=c(0.5,0.025,0.975),na.rm=TRUE),3))),row.names=FALSE,col.names=c('Year','Median Posterior Prediction','Lower 95th','Upper 95th'),sep=",")
cat('\nADF&G acoustic biomass (thousands of tons) posterior predictions\n')
write.table(cbind(years,t(round_df(apply(post.pred[[3]],2,quantile,probs=c(0.5,0.025,0.975),na.rm=TRUE),2))),row.names=FALSE,col.names=c('Year','Median Posterior Prediction','Lower 95th','Upper 95th'),sep=",")
cat('\nPWSSC acoustic biomass (thousands of tons) posterior predictions\n')
write.table(cbind(years,t(round_df(apply(post.pred[[4]],2,quantile,probs=c(0.5,0.025,0.975),na.rm=TRUE),2))),row.names=FALSE,col.names=c('Year','Median Posterior Prediction','Lower 95th','Upper 95th'),sep=",")
cat('\nAge-1 Aerial survey (# of small schools) posterior predictions\n')
write.table(cbind(years,t(round_df(apply(post.pred[[5]],2,quantile,probs=c(0.5,0.025,0.975),na.rm=TRUE),2))),row.names=FALSE,col.names=c('Year','Median Posterior Prediction','Lower 95th','Upper 95th'),sep=",")
sink()

