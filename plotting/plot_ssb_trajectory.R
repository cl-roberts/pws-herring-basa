# ssb_trajectory.R
# Created by Melissa Muradian
# Date updated:  12/14/2016
# Summary:
# This function plots two series: 
# 1) The pre-fishery biomass in each year, with the posterior median and 95% credibility interval
# 2) The probability that biomass falls below the management threshold (20,000 tons) in each year based on the posterior

## Plot biomass
ssb_trajectory <- function(PFRB,years,nYr,myMGP,modelPath,ylabel,xlabel,labels.x,labels.y,ylabel.2,labels.y2,ylim) {
    font.size <- 1.5
    #par(mar=c(3,3.5,1,3.5))
    finalPFRB<-PFRB[,nYr] # year 2015
    allIts<-length(finalPFRB) # number of thinned iterations after burn-in
    
    if(exists("ylim")){
      vert.max <- ylim
    }else{
      vert.max <- signif(max(apply(PFRB[,1:(nYr)]/1000,2,quantile,probs=.95)),2)
    }

    if(xlabel){
        xlabel <- "Year"
    }else{
        xlabel <- NA
    }
    
    if(ylabel){
        ylabel <- expression(Pre-fishery~biomass~(10^3~mt))
    }else{
        ylabel <- NA
    }
    
    plot(years[1:(nYr)], apply(PFRB[,1:(nYr)]/1000,2,FUN=quantile,probs=0.5),#rep(0,nYr), # CHANGE year
         type="l", lwd=2, las=1, xaxs="i", yaxs="i",
         xlab=xlabel,mgp=myMGP, ylab=ylabel,
         ylim=c(0,vert.max), axes=F, cex.lab=font.size, xpd=NA)
    axis(1, at=years, labels=F, tcl=-.15, col="grey60")
    axis(1,at= seq(years[1], years[length(years)], 5), labels=labels.x, tcl=-.3, lwd.ticks=2, cex.axis=font.size-0.6, mgp=myMGP, col="grey60")
    ## using alpha.level=c(.05,.5,.95) will plot the central 95%, 50% and 5%, in that order.
    source(file=paste0(modelPath,"function_add.polygon2.R"))
    add.polygon2(x=years[1:(nYr)], z=PFRB[,1:(nYr)]/1000, alpha.level=c(.05,.5,.95), alpha.min=.2, alpha.max=.9)
    axis(2, at=c(0,20,40,seq(100,vert.max,length.out = 2)), labels=labels.y, las=1, mgp=myMGP, cex.axis=font.size-0.6, col="grey60")
    abline(h=20, lty=2) 
    abline(h=40, lty=2)
    #lines(years, prob_below_threshold*200, lwd=2)
    
    count<-rep(0,nYr)
    prob_below_threshold<-rep(0,nYr)
    for(j in 1:(nYr)){ # CHANGE BACK TO JUST nYr 10/15/2015
    #   for(i in 1:allIts){
    #     #i<-1; j<-1
    #     if(PFRB[i,j]<20000)
    #       count[j]<-count[j]+1
    #   }
    #   print(count)
    #   prob_below_threshold[j]<-count[j]/allIts
        prob_below_threshold[j]<-sum(PFRB[,j] < 20000)/length(finalPFRB)
    }
    
    if(labels.y2){
      labels.y2 <- c(0,.25,.5,.75,1)
    }
    axis(4, at=seq(0,vert.max,vert.max/4), labels=labels.y2, las=1, mgp=myMGP, cex.axis=font.size-0.6, col="grey60")
    points(years[1:(nYr)], prob_below_threshold[1:(nYr)]*vert.max, pch=16, xpd=T) # CHANGE 10/15/2015
    lines(years[1:(nYr)], prob_below_threshold[1:(nYr)]*vert.max) # CHANGE 10/15/2015
    max.y <- 1.05*max(PFRB[,1:(nYr)]/1000)
    if(ylabel.2){
      mtext("Probability below threshold",side=4,line=3.5,xpd=NA, cex=font.size-0.2, srt=270)
    }
    return(prob_below_threshold)
}
