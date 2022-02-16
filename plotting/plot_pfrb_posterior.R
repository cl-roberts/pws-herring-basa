# pfrb_posterior.R
# Created by Melissa Muradian
# Date updated:  12/14/2016
# Summary:
# This function plots pre-fishery biomass posterior distribution (forecast), with the median and 95% credibility interval

## Plot histogram of projected final year biomass in mt
# 95% CI with median
pfrb_posterior <- function(PFRB, years, nyr, myMGP, xlim, ylabel, xlabel, labels.x, labels.y, ylim) {
    #par(mar=c(2.7,3,1,3))
    
    PFRB.quantiles <- matrix(NA, nyr, 3)
    PFRB.quantiles[, 1:3]<- t(apply(PFRB[, 1:nyr]/1000, 2, quantile, probs=c(0.025, 0.5, 0.975))) 
    PFRB.quantiles[, 1:3] <- apply(PFRB.quantiles, 2, round, 2)
    PFRB.quantiles[nyr, ]
    final.PFRB<-PFRB[, nyr] # year of projection

    h <- hist(final.PFRB/1000, breaks=60, plot=FALSE)
    #h.cuts <- cut(h$breaks, c(-Inf, PFRB.quantiles[nyr, 1], PFRB.quantiles[nyr, 3], Inf))

    #final.PFRB_2 <- final.PFRB[(final.PFRB/1000)>=PFRB.quantiles[nyr,1] & (final.PFRB/1000)<=PFRB.quantiles[nyr,3]]
    #pdf("FinalPFRB.pdf", height=4, width=6, family="Times")
    
    if(!exists("ylim")){
      ylim <- signif(max(h$density)*1.2, 1)
    }
    
    
    if(xlabel){
      xlabel <- paste0(nyr+1980-1, " pre-fishery biomass (1000s mt)")
    }else{
      xlabel <- NA
    }
    
    if(ylabel){
      ylabel <- " Probability Density"
    }else{
      ylabel <- NA
    }
    
    font.size <- 1

    plot(
        h,
        col="grey70",
        xlim=c(0, xlim), ylim=c(0, ylim), 
        xlab=xlabel, 
        ylab=ylabel, main="", yaxs="i",
        xaxs="i", axes=F, mgp=myMGP, lwd=4,
        cex.lab=font.size, freq=F, border="grey70", xpd=NA
    )
    
    axis(1, at=seq(0, xlim, by=10), mgp=myMGP, cex.axis=font.size-0.2, pos=0, col="grey60", labels=labels.x) 
    axis(2, at=c(0, ylim), las=2, mgp=myMGP, cex.axis=font.size-0.2, pos=0, col="grey60", labels=labels.y)
    segments(x0=PFRB.quantiles[nyr, ], x1=PFRB.quantiles[nyr, ],
            y0=c(0, 0, 0), y1=c(ylim*0.5, ylim, ylim*0.5), 
            lty=c(2, 1, 2), lwd=c(2, 3.5, 2))
    # legend("topright",legend="Management \nthreshold \n(20,000 mt)",lty=2,lwd=2,bty="n",col=1)
    
    # text(0.98*x.limit,y.limit*0.85,
    #      paste0("95% interval:\n(", PFRB.quantiles[nyr,1], ", ",
    #             PFRB.quantiles[nyr,3], ")"), cex=font.size-0.4,pos=2)
    
    # Probability below regulatory threshold
    count<-0
    allIts<-length(final.PFRB)
    for(i in 1:allIts){
      #i<-1; j<-1
      if(final.PFRB[i]<19958)
        count<-count+1
    }
    prob_below_threshold <- count/allIts
    #print(count)
    #prob_below_threshold <- sum(final.PFRB < 20000)/length(final.PFRB)
    
    text(0.85*xlim, ylim*0.94, 
        paste0("Median: ", PFRB.quantiles[nyr, 2], "\n\n",
               "95% interval:\n(", PFRB.quantiles[nyr, 1], ", ", PFRB.quantiles[nyr, 3], ")", "\n\n",
               "Probability below\nthreshold: ", format(prob_below_threshold, digits=3)), cex=font.size-0.2, pos=1)
    #text(0.98*x.limit,y.limit*0.69,paste0("Probability below\nthreshold: ",format(prob_below_threshold,digits=3)), cex=font.size-0.4,pos=2)
    #text(0.98*x.limit,y.limit*0.63,paste0("threshold: ",format(prob_below_threshold,digits=3)), cex=font.size-0.4,pos=2)
    #dev.off()
}

