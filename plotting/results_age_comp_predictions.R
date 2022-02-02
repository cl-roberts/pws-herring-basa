# results_age_comp_predictions.R
# Created by Melissa Muradian
# Modified by John Trochta
#   This code produces and plots posterior predictions 
#   of the age composition data

# Estimated numbers-at-age is also output from this file
# Last updated:  04/07/2021

modelPath<-here::here("model")
source(file=here::here("functions/data_reader.R"))
source(file=here::here("functions/data_header_reader.R"))

setwd(modelPath)

nburn <- 1
final.year <- 2021 # CHANGE with each update

#Color-reordering function: 
reorder<-function(vec,pos){
    x<-NULL
    for(i in 1:pos){
        x[i]<-vec[length(vec)-pos+i]
    }
    for(i in (pos+1):length(vec)){
        x[i]<-vec[i-pos]
    }
    return(x)
}

xx<-seq(1, by=1.2, length.out=7) #secret x-axis formula from D. Scott Rinnan

# Formerly using colors<-rich.colors(7) but needed to replace the yellow, 
#  so developed the pallette below
colors <- c("#13F24AFF", "#FFBA00FF", "firebrick1", "darkmagenta", "blue", "#00A4FFFF")

read.ADMB.files <- function(){
  # Parameters that are fixed and thus included as data
  #Z_3_8 <- 0.25
  egg_add <- 0.4
  
  # Store the file names from which data is available
  filename <- vector(length=2)
  filename[1]="PWS_ASA.dat"
  filename[2]="PWS_ASA(ESS).ctl"
  
  PWS_ASA.dat <- data.reader(filename=filename[1])
  PWS_ASA.dat <- PWS_ASA.dat[-2]
  nyr <- PWS_ASA.dat[[1]]
  nage <- PWS_ASA.dat[[2]]
  seine <- PWS_ASA.dat[[19]]
  spac <-  PWS_ASA.dat[[20]]
  
  PWS_ASA_ESS.ctl <- data.reader(filename=filename[2])
  ESS_Se <- PWS_ASA_ESS.ctl[[2]]
  ESS_Sp <- PWS_ASA_ESS.ctl[[3]]
  
  
  seine_indices <- which(rowSums(seine)>0)
  spawnsurvey_indices <- which(rowSums(spac)>0)
  
  model.data <- list(nyr=nyr,
                     nage=nage,
                     ESS_Se=ESS_Se,
                     ESS_Sp=ESS_Sp,
                     seine=seine,
                     spac=spac,
                     seine_indices=seine_indices,
                     spawnsurvey_indices=spawnsurvey_indices)
  return(model.data)
}
model.data <- read.ADMB.files()

nyr <- model.data$nyr 
ncol <- model.data$nage
seine.data <- model.data$seine
spawn.data <- model.data$spac

seine.data<-seine.data*100
spawn.data<-spawn.data*100

#comp<-cacomp[,2:8]*100
short.years<-seq(1980, 1992)
long.years<-seq(1980, final.year)   

#comp<-comp[-c(14),]
seine.data[10, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
seine.data[14:17, ] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
seine.data[20:length(long.years), ] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

# read in the ESS
seine.ESS <- model.data$ESS_Se
spawn.ESS <- model.data$ESS_Sp

seine.ESS[seine.ESS==-9] <- 0
spawn.ESS[spawn.ESS==-9] <- 0

# read in the model estimates for spawners
## spawn.age.comp is a # of saved draws X 231 matrix (33yrs * 7 age classes)
spawn.age.comp<-read.csv("mcmc_out/SpAC.csv", header = FALSE, dec=".") 
spawn.age.comp<-spawn.age.comp[-c(1:nburn), ]*100 # Just 'cause it's easier to work with this scale
###################
## read in the model estimates for seine catch
seine.age.comp<-read.csv("mcmc_out/SeAC.csv", header = FALSE, dec=".") 
seine.age.comp<-seine.age.comp[-c(1:nburn), ]*100 # Just 'cause it's easier to work with this scale


##############  Posterior Predictive Intervals #####################################################
## SEINE
pp.seine.age.comp<-matrix(0, NROW(seine.age.comp), NCOL(seine.age.comp)) # CHANGE
pp.spawn.age.comp<-matrix(0, NROW(spawn.age.comp), NCOL(spawn.age.comp))
# sample using the MCMC draws 
for(i in 1:length(seine.age.comp[, 1])){ # Loop through the MCMC draws
  for(j in 1:length(seine.ESS)){ # Loop through each year
    pp.seine.age.comp[i, (j*ncol-(ncol-1)):(j*ncol) ] <- t(rmultinom(1,size=as.integer(seine.ESS[j]),seine.age.comp[i, (j*ncol-(ncol-1)):(j*ncol)]))/seine.ESS[j]*100
    if(all(spawn.age.comp[i, (j*ncol-(ncol-1)):(j*ncol)] == 0)){
      pp.spawn.age.comp[i, (j*ncol-(ncol-1)):(j*ncol)] <- rep(0, times=ncol)
    }else{
      pp.spawn.age.comp[i, (j*ncol-(ncol-1)):(j*ncol) ] <- t(rmultinom(1,spawn.ESS[j], spawn.age.comp[i, (j*ncol-(ncol-1)):(j*ncol)])) / spawn.ESS[j] *100
    }
  }
}

## seine.age.comp.interval is a 2 X (# of ages in age comp * years) matrix; e.g. 7 ages (ages 3-9+) * 35 years (from 1980-2014)=35
seine.age.comp.interval <- matrix(rep(0, 2*length(seine.age.comp[1, ])), nrow=2, ncol=length(seine.age.comp[1, ]))
seine.age.comp.interval[1:2,] <- apply(pp.seine.age.comp, 2, quantile, probs=c(.025,.975), na.rm=T) # Fills first row with 0.025 percentile and second with 0.975 percentile
seine.age.comp.median <- apply(pp.seine.age.comp, 2, median)
## Now cut each up into years
seine.age.comp.median.mat <- matrix(seine.age.comp.median, nrow=length(seine.ESS[,1]), ncol=ncol, byrow=T) # nyr x nages matrix filled with posterior predictive median
seine.age.comp.lower.mat <- matrix(seine.age.comp.interval[1, ], nrow=length(seine.ESS[, 1]), ncol=ncol, byrow=T)
seine.age.comp.upper.mat <- matrix(seine.age.comp.interval[2, ], nrow=length(seine.ESS[, 1]), ncol=ncol, byrow=T)

## seine.age.comp.interval is a 2 X (# of ages in age comp * years) matrix; e.g. 7 ages (ages 3-9+) * 35 years (from 1980-2014)=35
spawn.age.comp.interval <- matrix(rep(0, 2*length(spawn.age.comp[1, ])), nrow=2, ncol=length(spawn.age.comp[1, ])) 
spawn.age.comp.interval[1:2, ] <-apply(pp.spawn.age.comp, 2, quantile, probs=c(0.025, 0.975), na.rm=T)
spawn.age.comp.median <- apply(pp.spawn.age.comp, 2, median)
## Now cut each up into years
spawn.age.comp.median.mat <- matrix(spawn.age.comp.median, nrow=length(spawn.ESS[, 1]), ncol=ncol, byrow=T)
spawn.age.comp.lower.mat <- matrix(spawn.age.comp.interval[1, ], nrow=length(spawn.ESS[, 1]), ncol=ncol, byrow=T)
spawn.age.comp.upper.mat <- matrix(spawn.age.comp.interval[2, ], nrow=length(spawn.ESS[, 1]), ncol=ncol, byrow=T)

names.arg <- c(3:8, "9+")

# Create the underlying foot print of the order I want to plot the 88 cells:
#fill <- c(1:11, 34:44, 56:66, 12:22, 45:55, 67:77, 23:33, 78:88)
#foot.print <- matrix(fill, nrow=11, byrow=F) ## an 11 by 8 matrix

# Below is from when I tried without the year columns for version 2, not using...
fig.rows <- ceiling(length(long.years)/3) # Split into 3 columns, according to Melissa's original layout
fill <- c(1:fig.rows, 
          (2*fig.rows+1):(3*fig.rows), 
          (fig.rows+1):(2*fig.rows), 
          (3*fig.rows+1):(4*fig.rows),
          (4*fig.rows+1):(5*fig.rows))
foot.print <- matrix(fill, nrow=14, byrow=F) ## an 18 by 8 matrix
my.bg <- "grey91"

##############   Plot it all !!!   #####################################################
#pdf(paste0(nyr,"_PP_age_comps.pdf"),height=7,width=6.6)
pdf(here::here("figures/predicted-age-comps.pdf"), height=8, width=6.6)
layout(foot.print)
#layout(foot.print, widths=c(1,3,3, 1,3,3, 1,3))
par(oma=c(4, 4, 4, 2.75), mar=rep(0, 4), bg="white")
#layout.show(55)

# Reduce seine and spawner comp matrices to correspond to ages 3-9+
seine.data <- seine.data[, -(1:3)]
seine.age.comp.median.mat <- seine.age.comp.median.mat[, -(1:3)]
seine.age.comp.lower.mat <- seine.age.comp.lower.mat[, -(1:3)]
seine.age.comp.upper.mat <- seine.age.comp.upper.mat[, -(1:3)]

spawn.data <- spawn.data[, -(1:3)]
spawn.age.comp.median.mat <- spawn.age.comp.median.mat[, -(1:3)]
spawn.age.comp.lower.mat <- spawn.age.comp.lower.mat[, -(1:3)]
spawn.age.comp.upper.mat <- spawn.age.comp.upper.mat[, -(1:3)]


## plot SEINE comps - only two columns
# Find missing
missing.seine <- which(apply(seine.data, MARGIN=1, FUN=function(x) all(x<=0)))
available.seine <- which(apply(seine.data, MARGIN=1, FUN=function(x) any(x>0)))
par(mar=c(0, 1, 0, .82))
for(i in 1:foot.print[length(foot.print[, 3]), 3]){
    if(i %in% missing.seine){ # Finds missing years in the purse-seine catch
        barplot(rep(NA, 7), col=c(reorder(colors, (i-1)%%6), "grey45"), 
                ylim=c(-8, 100), axes=F, border=NA)
        rect(par("usr")[1],
             par("usr")[3]-par("mar")[1],
             par("usr")[2]-par("mar")[4]+5,
             par("usr")[4]-par("mar")[3],
             col=my.bg, xpd=NA, border=NA)
        #text(x=8.5, y=75, labels=long.years[i], xpd=T, cex=1.2)
        text(x=8, y=80, labels=long.years[i], xpd=T, cex=1.2)
    }else{
        if( i == 1 | i == foot.print[length(foot.print[,1]),1] | i == foot.print[1,3] ){ # Finds those years at the top of the first, and bottoms of the first and second columns
            barplot(as.numeric(seine.data[i, 1:7]), col=c(reorder(colors, (i-1)%%6), "grey45"), 
                    ylim=c(-8, 100), yaxt="n", border=NA, mgp=c(2, 0.2, 0),  names.arg=c(3:8, "9+"))
            rect(par("usr")[1],
                 par("usr")[3]-par("mar")[1],
                 par("usr")[2]-par("mar")[4]+5,
                 par("usr")[4]-par("mar")[3],
                 col=my.bg, xpd=NA, border=NA)
            barplot(as.numeric(seine.data[i, 1:7]), col=c(reorder(colors,(i-1)%%6), "grey45"), 
                    ylim=c(-8, 100), yaxt="n", border=NA, mgp=c(2, 0.2, 0),  names.arg=c(3:8, "9+"), add=T)
            axis(2, at=c(0, 50, 100), labels=c("0.0", "0.5", "1.0"), las=2, cex.axis=.7, mgp=c(1, 0.4, 0), tcl=-.25)
            #axis(4, at=c(0,50,90) , labels=F, las=2, cex.axis=.7, mgp=c(1, .4, 0), tcl=-.25)
            #text(x=8.5, y=75, labels=long.years[i], xpd=T, cex=1.2)
            text(x=8, y=80, labels=long.years[i], xpd=T, cex=1.2)
            points(x=xx-0.3, y=seine.age.comp.median.mat[i, ], pch=20, cex=0.7, xpd=T)
            segments(xx-0.3, seine.age.comp.lower.mat[i, ], xx-0.3, seine.age.comp.upper.mat[i, ], lwd=1.5, lend=1) 
        }else{
            barplot(as.numeric(seine.data[i, 1:7]), col=c(reorder(colors, (i-1)%%6), "grey45"), 
                    ylim=c(-8, 100), axes=F, border=NA)
            rect(par("usr")[1],
                 par("usr")[3]-par("mar")[1],
                 par("usr")[2]-par("mar")[4]+5,
                 par("usr")[4]-par("mar")[3],
                 col=my.bg, xpd=NA, border=NA)
            barplot(as.numeric(seine.data[i, 1:7]), col=c(reorder(colors, (i-1)%%6), "grey45"), 
                    ylim=c(-8, 100), axes=F, border=NA, add=T)
            #axis(2, at=c(0,50,90) , labels=c("0.0","0.5","0.9"), las=2, cex.axis=.7, mgp=c(1, .4, 0), tcl=-.25)
            #text(x=8.5, y=75, labels=long.years[i], xpd=T, cex=1.2)
            text(x=8, y=80, labels=long.years[i], xpd=T, cex=1.2)
            points(x=xx-0.3, y=seine.age.comp.median.mat[i, ], pch=20, cex=0.7, xpd=T)
            segments(xx-0.3, seine.age.comp.lower.mat[i, ], xx-0.3, seine.age.comp.upper.mat[i, ], lwd=1.5, lend=1) 
        }
    }
}
axis(2, at=c(0, 50, 100), labels=c("0.0", "0.5", "1.0"), las=2, cex.axis=0.7, mgp=c(1, .04, 0), tcl=-0.25)

## plot SPAWN comps - three columns
par(mar=c(0, 1, 0, .5))
for(i in 1:2){ 
    barplot(rep(NA, 7), col=c(reorder(colors, (i-1)%%6), "grey45"), 
            ylim=c(-8, 100), axes=F, border=NA)
    #plot(0:1,0:1, type= "n", axes=F, ann=F)
    rect(par("usr")[1],
         par("usr")[3]-par("mar")[1],
         par("usr")[2]-par("mar")[4]+.5,
         par("usr")[4]-par("mar")[3],
         col=my.bg, xpd=NA, border=NA)
    #mtext(at=c(-.18), text=long.years[i], xpd=T, cex=.85, padj=1)
    # text(x=-1.35, y=80, labels=long.years[i], xpd=T, cex=1.2)
} # Plots empty gray spaces in first two years (1980-1981) of spawn survey

first.column.bottom <- foot.print[length(foot.print[, 1]), 1]
second.column.bottom <- foot.print[length(foot.print[, 3]), 3]
for(i in 3:second.column.bottom){
    if(i == first.column.bottom | i== second.column.bottom){
        barplot(as.numeric(spawn.data[i, 1:7]), col=c(reorder(colors, (i-1)%%6), "grey45"), ylim=c(-8, 100), 
                yaxt="n", border=NA, mgp=c(2, 0.2, 0),  names.arg=c(3:8, "9+"))
        rect(par("usr")[1],
             par("usr")[3]-par("mar")[1],
             par("usr")[2]-par("mar")[4]+.5,
             par("usr")[4]-par("mar")[3],
             col=my.bg, xpd=NA, border=NA)
        barplot(as.numeric(spawn.data[i, 1:7]), col=c(reorder(colors, (i-1)%%6), "grey45"), ylim=c(-8, 100), 
                yaxt="n", border=NA, mgp=c(2, 0.2, 0),  names.arg=c(3:8, "9+"), add=T)
        #axis(2, at=c(0,50,90) , labels=c("0.0","0.5","0.9"), las=2, cex.axis=.7, mgp=c(1, .4, 0), tcl=-.25)
        #text(x=8.5, y=75, labels=long.years[i], xpd=T, cex=1.2)
        #mtext(at=c(-1.25), text=long.years[i], xpd=T, cex=.85, padj=1)
        #text(x=-1.35, y=80, labels=long.years[i], xpd=T, cex=1.2)
        points(x=xx-.3, y=seine.age.comp.median.mat[i, ], pch=20, cex=.7, xpd=T)
        segments(xx-.3, spawn.age.comp.lower.mat[i, ], xx-.3, spawn.age.comp.upper.mat[i, ], lwd=1.5, lend=1)
    }
    else{
        barplot(as.numeric(spawn.data[i, 1:7]), col=c(reorder(colors, (i-1)%%6), "grey45"), ylim=c(-8, 100), axes=F, 
                border=NA)
        rect(par("usr")[1],
             par("usr")[3]-par("mar")[1],
             par("usr")[2]-par("mar")[4]+.5,
             par("usr")[4]-par("mar")[3],
             col=my.bg, xpd=NA, border=NA)
        barplot(as.numeric(spawn.data[i, 1:7]), col=c(reorder(colors, (i-1)%%6), "grey45"), ylim=c(-8, 100), axes=F, 
                border=NA, add=T)
        #axis(2, at=c(0,50,90) , labels=c("0.0","0.5","0.9"), las=2, cex.axis=.7, mgp=c(1, .4, 0), tcl=-.25)
        #text(x=8.5, y=75, labels=long.years[i], xpd=T, cex=1.2)
        #mtext(at=c(-1.25), text=long.years[i], xpd=T, cex=.85, padj=1)
        points(x=xx-0.3, y=seine.age.comp.median.mat[i, ], pch=20, cex=.7, xpd=T)
        segments(xx-0.3, spawn.age.comp.lower.mat[i, ], xx-0.3, spawn.age.comp.upper.mat[i, ], lwd=1.5, lend=1)
        #text(x=-1.35, y=80, labels=long.years[i], xpd=T, cex=1.2)
    }
}
par(mar=c(0, 1, 0, .5))
for(i in (second.column.bottom+1):length(long.years)){
    if(i == (second.column.bottom+1) | i == length(long.years)){
        barplot(as.numeric(spawn.data[i, 1:7]), col=c(reorder(colors, (i-1)%%6), "grey45"), ylim=c(-8, 100), 
                yaxt="n", border=NA, mgp=c(2, 0.2, 0),  names.arg=c(3:8, "9+"))
        rect(par("usr")[1],
             par("usr")[3]-par("mar")[1],
             par("usr")[2]-par("mar")[4]+.5,
             par("usr")[4]-par("mar")[3],
             col=my.bg, xpd=NA, border=NA)
        barplot(as.numeric(spawn.data[i, 1:7]), col=c(reorder(colors, (i-1)%%6), "grey45"), ylim=c(-8, 100), 
                yaxt="n", border=NA, mgp=c(2, 0.2, 0),  names.arg=c(3:8, "9+"), add=T)
        axis(2, at=c(0, 50, 100), labels=c("0.0", "0.5", "1.0"), las=2, cex.axis=.7, mgp=c(1, 0.4, 0), tcl=-.25)
        text(x=5, y=80, labels=long.years[i], xpd=T, cex=1.2)
        #mtext(at=c(-1.25), text=long.years[i], xpd=T, cex=.85, padj=1)
        points(x=xx-0.3, y=seine.age.comp.median.mat[i, ], pch=20, cex=0.7, xpd=T)
        segments(xx-0.3, spawn.age.comp.lower.mat[i, ], xx-0.3, spawn.age.comp.upper.mat[i, ], lwd=1.5, lend=1)  
    }
    else{
        barplot(as.numeric(spawn.data[i, 1:7]), col=c(reorder(colors, (i-1)%%6), "grey45"), ylim=c(-8, 100), axes=F, 
                border=NA)
        rect(par("usr")[1],
             par("usr")[3]-par("mar")[1],
             par("usr")[2]-par("mar")[4]+.5,
             par("usr")[4]-par("mar")[3],
             col=my.bg, xpd=NA, border=NA)
        barplot(as.numeric(spawn.data[i, 1:7]), col=c(reorder(colors, (i-1)%%6), "grey45"), ylim=c(-8, 100), axes=F, 
                border=NA, add=T)
        #axis(2, at=c(0,50,90) , labels=c("0.0","0.5","0.9"), las=2, cex.axis=.7, mgp=c(1, .4, 0), tcl=-.25)
        text(x=5, y=80, labels=long.years[i], xpd=T, cex=1.2)
        points(x=xx-0.3, y=seine.age.comp.median.mat[i, ], pch=20, cex=0.7, xpd=T)
        segments(xx-0.3, spawn.age.comp.lower.mat[i, ], xx-0.3, spawn.age.comp.upper.mat[i, ], lwd=1.5, lend=1)
    }
}

# Label the individual columns
# NOTE: spaces needed since the text is xx-justified using 'adj' argument
mtext(side=3, text="Purse-seine \n    catch",        adj=0.05, outer=T, cex=0.9) # far left-justified  
mtext(side=3, text="Herring-spawn \n   survey",      adj=0.27, outer=T, cex=0.9) 
mtext(side=3, text="Purse-seine \n catch   ",        adj=0.51, outer=T, cex=0.9) # centered-ish
mtext(side=3, text="Herring-spawn \n survey       ", adj=0.74, outer=T, cex=0.9) # right-justified
mtext(side=3, text="Herring-spawn \n survey       ", adj=1.00, outer=T, cex=0.9) # far right-justified
## Main axis labels
mtext(side=1, text="Age Class", outer=T, line=1.25, cex=1)
#mtext(side=2, text="Year", line=0, cex=1, outer=T)
mtext(side=2, text="Proportion", line=0, cex=1, outer=T)
#mtext(side=4, text="Proportion", line=.95, cex=1, outer=T)
dev.off()

###############################################################
## Save csv of values displayed in table above (with 95% posterior predictive intervals) 
round_df <- function(x, digits) {
    # round all numeric variables
    # x: data frame 
    # digits: number of digits to round
    numeric.columns <- sapply(x, mode) == "numeric"
    x[numeric.columns] <-  round(x[numeric.columns], digits)
    return(x)
}

sink(here::here("data_outputs/predicted-age-comps.csv"))
cat("# Seine Fishery Age Composition (proportion of each age): Median Posterior Prediction\n")
write.table(
    cbind(long.years, round_df(seine.age.comp.median.mat/100, 3)), 
    row.names=FALSE, col.names=c("# Year", 3:8, "9+"), sep=","
)

cat("\n# Seine Fishery Age Composition (proportion of each age): Lower 95% Posterior Predictive Interval\n")
write.table(
    cbind(long.years, round_df(seine.age.comp.lower.mat/100, 3)),  
    row.names=FALSE, col.names=c("# Year", 3:8, "9+"), sep=","
)

cat("\n# Seine Fishery Age Composition (proportion of each age): Upper 95% Posterior Predictive Interval\n")
write.table(
    cbind(long.years, round_df(seine.age.comp.upper.mat/100, 3)),  
    row.names=FALSE, col.names=c("# Year", 3:8, "9+"), sep=","
)

cat("\n\n# Spawner Survey Age Composition (proportion of each age): Median Posterior Prediction\n")
write.table(
    cbind(long.years, round_df(seine.age.comp.median.mat/100, 3)), 
    row.names=FALSE, col.names=c("# Year", 3:8, "9+"), sep=","
)
cat("\n# Spawner Survey Age Composition (proportion of each age): Lower 95% Posterior Predictive Interval\n")
write.table(
    cbind(long.years, round_df(spawn.age.comp.lower.mat/100, 3)),  
    row.names=FALSE, col.names=c("# Year", 3:8, "9+"), sep=","
)
cat("\n# Spawner Survey Age Composition (proportion of each age): Upper 95% Posterior Predictive Interval\n")
write.table(
    cbind(long.years, round_df(spawn.age.comp.upper.mat/100, 3)),  
    row.names=FALSE, col.names=c("# Year", 3:8, "9+"), sep=","
)
sink()

###############################################################
## Read-in Numbers-at-age and output into nice tables (.csv) 

n.y.a<-read.csv(here::here("model/mcmc_out/Num_at_age.csv"), header = FALSE, dec=".") 
n.y.a<-n.y.a[-c(1:nburn), ] # Just 'cause it's easier to work with this scale

## seine.age.comp.interval is a 2 X (# of ages in age comp * years) matrix; e.g. 10 ages (ages 0-9+) * 35 years (from 1980-2014)=350
n.y.a.interval <- matrix(rep(0, 2*length(n.y.a[1, ])), nrow=2, ncol=length(n.y.a[1, ]))
n.y.a.interval[1:2,] <- apply(n.y.a, 2, quantile, probs=c(0.025, 0.975), na.rm=T) # Fills first row with 0.025 percentile and second with 0.975 percentile
n.y.a.median <- apply(n.y.a, 2, median)
## Now cut each up into years
n.y.a.median.mat <- matrix(n.y.a.median, nrow=nyr, ncol=ncol, byrow=T) # nyr x nages matrix filled with posterior predictive median
n.y.a.lower.mat <- matrix(n.y.a.interval[1, ], nrow=nyr, ncol=ncol, byrow=T)
n.y.a.upper.mat <- matrix(n.y.a.interval[2, ], nrow=nyr, ncol=ncol, byrow=T)

sink(here::here("data_outputs/numbers-at-age.csv"))
cat("# Numbers-at-age (millions): Median estimate\n")
write.table(
    cbind(long.years, round_df(n.y.a.median.mat, 3)), 
    row.names=FALSE, col.names=c("# Year", 0:8, "9+"), sep=","
)

cat("\n# Numbers-at-age (millions): Lower 95% Credibility Interval\n")
write.table(
    cbind(long.years, round_df(n.y.a.lower.mat, 3)),  
    row.names=FALSE, col.names=c("# Year", 0:8, "9+"), sep=","
)

cat("\n# Numbers-at-age (millions): Upper 95% Credibility Interval\n")
write.table(
    cbind(long.years, round_df(n.y.a.upper.mat, 3)),  
    row.names=FALSE, col.names=c("# Year", 0:8, "9+"), sep=","
)

sink()
