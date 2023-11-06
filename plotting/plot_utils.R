library(here)
library(ggplot2)
library(tidyverse)
library(ggdist)

color.palette <- "Blues"

present.theme <- theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = ifelse(legend, "right", "none"),
    axis.title = element_text(size=16, face="bold"),
    axis.text = element_text(size=14),
    plot.title = element_text(size=20, face="bold")
)

plot.biomass.trajectory <- function(df, years, legend=TRUE, show.probs=TRUE, new.theme=NA){
  
    # the.theme <- ifelse(
    #     is.na(new.theme),
    #     theme(
    #         panel.grid.minor = element_blank(),
    #         panel.grid.major = element_blank(),
    #         panel.background = element_blank(),
    #         axis.line = element_line(colour = "black"),
    #         legend.position = ifelse(legend, "right", "none")
    #     ),
    #     new.theme
    # )

    plot <- ggplot(df ) +
      geom_lineribbon(aes(x=year, y=biomass/1000, ymin=.lower/1000, ymax=.upper/1000, group=1), size=0.25)+
      scale_fill_grey(start=0.8, end=0.4)+
      #scale_fill_brewer(palette = color.palette)+
      scale_x_discrete("Year", breaks=seq(min(years), max(years), by=5), expand=c(0, 0))+
      geom_hline(yintercept=c(20, 40), linetype="dashed")+
      geom_vline(xintercept = length(years)-1, linetype="dashed")+
      ggtitle("Biomass Trajectory")+
      scale_y_continuous(
        "Pre-fishery biomass (1000 mt)", 
        breaks=c(0, 20, 40, 50, 100, 150, 200),
        expand=c(0, 0)
      )+
      coord_cartesian(clip="off")+
      theme(
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = ifelse(legend, "right", "none"),
        axis.title = element_text(size=16, face="bold"),
        axis.text = element_text(size=14),
        plot.title = element_text(size=20, face="bold")
      )
    
    if(show.probs){
      plot <- plot + geom_point(aes(x=year, y=prob*200), size=2)+
              geom_line(aes(x=year, y=prob*200, group=1), size=0.8)+
              scale_y_continuous(
                "Pre-fishery biomass (1000 mt)", 
                breaks=c(0, 20, 40, 50, 100, 150, 200),
                expand=c(0, 0),
                sec.axis = sec_axis(trans=~.*1/200, name="Probability below 20k metric tons")
              )
    }
  
    return(
        plot
    )
}

plot.exploit.rate <- function(df, zeros, years){
    return(
        ggplot(df, aes(x=year, y=exploit, ymin=.lower, ymax=.upper))+
            geom_pointinterval(size=0.25)+
            geom_pointinterval(data=zeros, color="black", shape=4)+
            scale_x_discrete("Year", breaks=seq(min(years), max(years), 5), expand=c(0, 0))+
            scale_y_continuous("Exploitation rate", breaks=seq(0, 0.3, 0.05), expand=c(0, 0))+
            ggtitle("Exploitation Rate")+
            coord_cartesian(clip="off")+
            theme(
                panel.grid.minor = element_blank(),
                panel.grid.major = element_blank(),
                panel.background = element_blank(),
                axis.line = element_line(colour = "black"),
                axis.title = element_text(size=16, face="bold"),
                axis.text = element_text(size=14),
                plot.title = element_text(size=20, face="bold")
            )
    )
}

plot.pfrb.posterior <- function(df, quants, prob, curr.year, font.size=1){

    extra <- data.frame(q1=quants[1], q2=quants[2], q3=quants[3], prob=prob)

    return(
        ggplot(df)+
            geom_histogram(aes(x=biomass/1000, y= ..density..), bins=60)+
            scale_fill_grey(start=0.8, end=0.6)+
            geom_vline(xintercept = quants, linetype=c("dashed", "solid", "dashed"), size=c(0.5, 1, 0.5))+
            geom_text(data=extra, aes(x=58, y=0.14, label=paste("Median:", q2)), size=font.size, hjust=1)+
            geom_text(data=extra, aes(x=58, y=0.12, label=paste0("95% interval:\n", "(", q1, ", ", q2, ")")), size=font.size, hjust=1)+
            geom_text(data=extra, aes(x=58, y=0.09, label=paste("Probability below\nthreshold:", prob)), size=font.size, hjust=1)+
            scale_x_continuous(paste(curr.year-1, "Pre-Fishery Biomass (mt)"), breaks=seq(0, 60, 5), expand=c(0, 0))+
            scale_y_continuous("Probability density", breaks=seq(0, 0.15, 0.05), expand=c(0, 0))+
            coord_cartesian(ylim=c(0, 0.15), xlim=c(0, 60))+
            ggtitle(paste(curr.year-1, "Pre-fishery Biomass Posterior Probability Density"))+
            theme(
                panel.grid.minor = element_blank(),
                panel.grid.major = element_blank(),
                panel.background = element_blank(),
                axis.line = element_line(colour = "black"),
                axis.title = element_text(size=16, face="bold"),
                axis.text = element_text(size=14),
                plot.title = element_text(size=20, face="bold")
            )
    )
}

plot.recruitment.posterior <- function(df, years, legend=TRUE){
     return(
        ggplot(df) +
            geom_lineribbon(aes(x=year, y=recruits, ymin=.lower, ymax=.upper, group=1), size=0.25)+
            scale_fill_grey(start=0.8, end=0.6)+
            scale_x_discrete("Year", breaks=seq(min(years), max(years), by=5), expand=c(0, 0))+
            scale_y_continuous("Age-3 recruits (millions)", breaks=seq(0, 2000, by=500), limits=c(0, 2000), expand=c(0, 0))+
            ggtitle("Age-3 Recruitment")+
            theme(
                panel.grid.minor = element_blank(),
                panel.grid.major = element_blank(),
                panel.background = element_blank(),
                axis.line = element_line(colour = "black"),
                legend.position = ifelse(legend, "right", "none"),
                axis.title = element_text(size=16, face="bold"),
                axis.text = element_text(size=14),
                plot.title = element_text(size=20, face="bold")
            )
    )
}

