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

create.retrospective.matrix <- function(n.peels, fname="PFRBiomass.csv"){

    base <- as.numeric(read_csv(paste0(here::here(), "/model/mcmc_out/", fname), col_names=FALSE, show_col_types = FALSE) %>%
                                            summarise(
                                                across(
                                                everything(),
                                                median
                                                )
                                            ))
    nyr <- length(base)

    # Aggregate biomass estimates for each retrospective run
    annual.estimate <- matrix(NA, nyr, n.peels+1)
    annual.estimate[1:nyr, 1] <- base

    for(i in n.peels:1){
        print(i)
        est <- read_csv(paste0(here::here(), "/retrospectives/", "basa-", i,"/mcmc_out/", fname), col_names=FALSE, show_col_types = FALSE) %>%
                summarise(
                    across(
                    everything(),
                    median
                    )
                ) 
        
        est <- t(as.matrix(est))
        rownames(est) <- NULL
        for(j in 1:length(est)){
            print(j)
            annual.estimate[j, i+1] <- est[j, 1]
        }
    }

    rownames(annual.estimate) <- 1980:(1980+nyr-1)
    colnames(annual.estimate) <- c("base", paste0("-", 1:n.peels))
    
    return(annual.estimate)
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

########################## Analaysis and Plotting Code ###########################

annual.biomass.estimate <- create.retrospective.matrix(5)

mohns.rho <- icesAdvice::mohn(annual.biomass.estimate, details=TRUE)
rel.bias <- c(mohns.rho$compare[,"relbias"], NA)
rho <- mohns.rho$rho

# Compute the peel estimates for each run
#final.year.est <- annual.biomass.estimate[38:nrow(annual.biomass.estimate), 1]
peel.estimates <- apply(as.matrix(annual.biomass.estimate), 2, function(x) x[max(which(!is.na(x)))])
peel.estimates <- rev(peel.estimates)


# Plot retrospective plot (can also be donee with icesAdvice::mohn(..., plot=TRUE))
colors <- palette(rainbow(n.peels+1))

yr = 2024

plot(1980:yr, rep(NA, length(1980:yr)), ylim=c(0, 150000), xlab="Year", ylab="SSB (mt)", main="Retrospective Analaysis of SSB")
for(p in 1:n.peels){
    data <- annual.biomass.estimate[,1+p]
    lines(1980:yr, data, col=colors[p], lwd=2)
    points(yr-p, rev(peel.estimates)[1+p], col=colors[p], pch=19, lwd=10)
}
lines(1980:yr, annual.biomass.estimate[,1], col="black", lwd=2.3)
text(2019.5, 110000, paste0("Mohn's Rho: ", round(rho, 3)))
abline(h=19958, lty=2)
abline(h=39885, lty=2)

# Fancy thing for adding the relative biases to the legend
legs <- rep(NA, length((yr-n.peels):yr))
for(i in 1:length(legs)){
    legs[i] <- paste0(yr-n.peels-1+i, " (bias = ", round(rel.bias[i], 3), ")")
}
legend(x=2014.5, y=150000, legend=rev(legs), col=c("black", colors), lwd=2)

# Make squid plot
# plot(1:n.peels, rep(NA, n.peels), xlim=c(0, n.peels+1), ylim=c(0, 150000), xlab="")
# colors <- palette(rainbow(n.peels+1))
# for(i in 1:(n.peels+1)){
#     y <- 2024-1-n.peels+i
#     data <- rev(annual.biomass.estimate[as.character(y),])
#     non.na.index <- which(!is.na(data))[1]
#     data <- data[non.na.index:length(data)]
#     lines(0:length(data), c(0, data), col=colors[i], lwd=3)
#     text(x=length(data), y=data[length(data)]+2500, y)
# }

# Make reverse squid plot
# plot(1:n.peels, rep(NA, n.peels), xlim=c(0.5, n.peels+1), ylim=c(-2, 2), xlab="")
# colors <- palette(rainbow(n.peels+1))
# for(i in 1:(n.peels)){
#     y <- 2024-1-n.peels+i
#     data <- rev(annual.biomass.estimate[as.character(y),])
#     data <- (data - data[length(data)])/data[length(data)]
#     non.na.index <- which(!is.na(data))[1]
#     data <- data[non.na.index:length(data)]
#     lines(1:length(data), c(data), col=colors[i], lwd=3)
#     text(x=0.75, y=data[1], y)
# }
# abline(h=0, lty=2)


retro <- as_tibble(annual.biomass.estimate) %>% 
    pivot_longer(everything(), names_to="lag", values_to="biomass") %>%
    arrange(lag) %>%
    mutate(year=rep(1980:yr, length.out=n())) %>%
    mutate(
        lag = factor(lag, levels=c("base", paste0("-", 1:n.peels)), labels = c(yr:(yr-n.peels)))
    )

biomass.df <- compute.biomass.traj(here::here("model"), length(1980:(yr-1)), 1980:(yr-1)) %>% mutate(year=as.numeric(year)) %>% filter(.width >= 0.5)

ggplot(retro) +
    geom_line(aes(x=year, y=biomass, color=lag), size=1)+
    geom_lineribbon(data=biomass.df, aes(x=year, y=biomass, ymin=.lower, ymax=.upper), size=0, alpha=0.35)+
    geom_hline(yintercept=c(20000, 40000), linetype="dashed")+
    scale_fill_grey(start=0.8, end=0.4)+
    scale_y_continuous(limits=c(0, 200000), breaks=c(0, 20000, 40000, 50000, 100000, 150000, 200000), labels=scales::comma)+
    coord_cartesian(expand=0)+
    labs(x="Year", y="Spawning Biomass (mt)", color="Data Lag", title=paste0(n.peels, "-year Retrospective Pattern"))+
    theme_classic()+
    theme(
        axis.title = element_text(size=16),
        axis.text = element_text(size=12),
        plot.title = element_text(size=18)
    )


##########
#########
##########

create.retrospective.matrix.probs <- function(n.peels, fname="PFRBiomass.csv"){

    nyr.global <- yr-1980

    base <- compute.biomass.traj(here::here("model"), nyr.global, 1980:yr)
    prob <- base$prob[1:nyr.global]

    # Aggregate biomass estimates for each retrospective run
    annual.estimate <- matrix(NA, nyr.global, n.peels+1)
    annual.estimate[1:nyr.global, 1] <- prob

    for(i in n.peels:1){
        nyr <- nyr.global-i
        est <- compute.biomass.traj(file.path(here::here("retrospectives"), paste0("basa-", i)), nyr, years=1980:(yr-i))
        prob <- est$prob[1:nyr]
        #annual.estimate[,i] <- 
        #prob <- t(as.matrix(prob))
        #rownames(prob) <- NULL
        for(j in 1:length(prob)){
           annual.estimate[j, i+1] <- prob[j]
        }
    }

    rownames(annual.estimate) <- 1980:(1980+nyr.global-1)
    colnames(annual.estimate) <- c("base", paste0("-", 1:n.peels))
    
    return(annual.estimate)
}

retro.prob <- as_tibble(create.retrospective.matrix.probs(n.peels)) %>% 
    pivot_longer(everything(), names_to="lag", values_to="prob") %>%
    arrange(lag) %>%
    mutate(year=rep(1980:(yr-1), length.out=n())) %>%
    mutate(
        lag = factor(lag, levels=c("base", paste0("-", 1:n.peels)), labels = c(yr:(yr-n.peels)))
    )

ggplot(retro.prob) +
    geom_line(aes(x=year, y=prob, color=lag), size=1)+
    geom_point(aes(x=year, y=prob, color=lag), size=4)+
    #geom_lineribbon(data=biomass.df, aes(x=year, y=biomass, ymin=.lower, ymax=.upper), size=0, alpha=0.35)+
    geom_hline(yintercept=c(0.50), linetype="dashed")+
    #scale_fill_grey(start=0.8, end=0.4)+
    #scale_y_continuous(limits=c(0, 200000), breaks=c(0, 20000, 40000, 50000, 100000, 150000, 200000), labels=scales::comma)+
    coord_cartesian(expand=0)+
    labs(x="Year", y="Probability", color="Data Lag", title=paste0(n.peels, "-year Retrospective Pattern"))+
    theme_classic()+
    theme(
        axis.title = element_text(size=16),
        axis.text = element_text(size=12),
        plot.title = element_text(size=18)
    )

m = icesAdvice::mohn(create.retrospective.matrix.probs(n.peels))
