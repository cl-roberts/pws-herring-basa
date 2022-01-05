# rec_posterior.R
# Created by Melissa Muradian
# Date updated:  12/14/2016
# Summary:
# This function plots the recruitment trajectory (millions of Age-3's) posteriors, as the median and 95% credibility intervals in each year

## Plot Age-3 Rec
rec_posterior <- function(years, nYr, Age3, myMGP, modelPath, ylim, ylabel, xlabel, labels.x, labels.y) {
    source(file=paste0(modelPath, "function_add.polygon2.R"))

    font.size <- 1.8
    
    if(xlabel){
        xlabel <- "Year"
    }else{
        xlabel <- NA
    }
    if(ylabel){
        ylabel <- expression(Age-3~recruits~(10^6))
    }else{
        ylabel <- NA
    }

    plot(years, rep(0, nYr), type="l", lwd=2, las=1, xaxs="i",
        yaxs="i", axes=F, ylab=ylabel, xlab=xlabel, 
        mgp=myMGP,  ylim=c(0,ylim), cex.lab=font.size, xpd=NA)
    add.polygon2(x=years, z=Age3, alpha.level=c(.05,.5,.95),
                alpha.min=.2, alpha.max=.9)
    axis(1, at=years, labels=F, tcl=-.15, col="grey60")
    axis(1, at=seq(years[1], years[length(years)], 5),
        labels=labels.x, tcl=-.35, lwd.ticks=2, cex.axis=font.size-0.6, mgp=myMGP, col="grey60")
    vert.max <- ylim
    axis(2, at=seq(0,vert.max,500), labels=labels.y, las=1, mgp=myMGP,
        cex.axis=font.size-0.6, col="grey60")
}
