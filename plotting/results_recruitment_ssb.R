# plotting_outputs.R

# BASA_vhsv_outputs.r
# Created by John Trochta
# Date updated:  11/17/2020
# Summary:
###################################

library(dplyr)
library(ggplot2)

model.path <- here::here("model") 

models <- c("base/")
model.name <- c("Base model")

# Years of the model
years<-1980:2021
nyr<-length(years)+1

dat.plot.2 <- data.frame()

for(j in 1:length(models)){
  #setwd(paste0(modelPath,models[j]))
  setwd(model.path)
  # Read in key outputs for VHSV outbreak severity
  ssb <- read.table("mcmc_out/PFRBiomass.csv", 
                    header = FALSE, 
                    sep = ",", 
                    dec=".",
                    col.names = paste0(c(years, max(years)+1)),
                    check.names = FALSE, 
                    stringsAsFactors=FALSE
                    )
  age3 <- read.table("mcmc_out/Age3.csv", 
                     header = FALSE, 
                     sep = ",", 
                     dec=".",
                     col.names = paste0(years),
                     check.names = FALSE, 
                     stringsAsFactors=FALSE
                     )
  
  # Combine quantities for output
  dat.plot <- rbind(
    data.frame(Variable="SSB (metric tons)", reshape2::melt(ssb, variable.name=c("Year"), factorsAsStrings = FALSE)),
    data.frame(Variable="Age 3 (millions)", reshape2::melt(age3, variable.name=c("Year"), factorsAsStrings = FALSE))
  )
  
  #dat.plot <- dat.plot %>% filter(Year %in% 2008:(max(years)-1))
  
  # Calc 95% quantiles for each year
  dat.plot <- dat.plot %>% group_by(Variable, Year) %>%  
      summarize(
          Q.025=quantile(value,probs=0.025, na.rm=TRUE),
          Q.250=quantile(value,probs=0.250, na.rm=TRUE),
          Q.500=quantile(value,probs=0.500, na.rm=TRUE),
          Q.750=quantile(value,probs=0.750, na.rm=TRUE),
          Q.975=quantile(value,probs=0.975, na.rm=TRUE)
      )

  dat.plot.2 <- bind_rows(dat.plot.2, data.frame(model=model.name[j], dat.plot))
}

# Estimate stats for paper
sum.stat <- dat.plot.2 %>% group_by(model, Variable) %>% 
    summarise(
        avg.rates=median(Q.500),
        min.rates=min(Q.500),
        year.min=Year[which.min(Q.500)],
        max.rates=max(Q.500),
        year.max=Year[which.max(Q.500)]
    )

sum.stat.2 <- dat.plot.2 %>% group_by(Variable, Year) %>% summarise(max.diff=max(Q.500)-min(Q.500))

max.bounds <- dat.plot.2 %>% group_by(Variable) %>% summarise(max.val=max(Q.975))
sub_panel_labs <- data.frame(Year=0.5,
                             Variable=max.bounds$Variable,
                             value=max.bounds$max.val*0.95,
                             label=LETTERS[1:2])


min.fish.threshhold <- c(rep(NA, length(years)), rep(20000, length(years)+1))
max.fish.threshhold <- c(rep(NA, length(years)), rep(40000, length(years)+1))

dat.plot.2$min.threshhold <- min.fish.threshhold
dat.plot.2$max.threshhold <- max.fish.threshhold

font.size <- 11
ggplot(data=dat.plot.2, aes(x=Year, y=Q.500, group=Variable)) + 
  geom_ribbon(aes(ymin=Q.025, ymax=Q.975), fill="grey70")+
  geom_ribbon(aes(ymin=Q.250, ymax=Q.750), fill="grey85")+
  geom_hline(aes(yintercept=min.threshhold), color="red")+
  geom_hline(aes(yintercept=max.threshhold), color="red")+
  geom_line(size=1)+
  scale_x_discrete(breaks=c(seq(1980, 2021, 3)))+
  scale_y_continuous(limits=c(0,NA))+
  facet_grid(Variable~model, switch="y", scales="free_y")+
  #scale_x_discrete(limits=c(0,50),expand=c(0,0))+
  geom_text(data=sub_panel_labs, aes(x=Year, y=value, label=label,
                                    group=Variable, size=font.size+3), hjust=-1)+
  theme_classic()+
  ggtitle("Age 3 and Spawning Stock Biomass of PWS Herring")+
  theme(strip.background = element_blank(),
        panel.border=element_rect(fill=NA),
        # plot.title = element_text(hjust = 0.5),
        #plot.title = element_blank(),
        # strip.text.x = element_text(size=font.size),
        strip.text.x = element_blank(),
        strip.text.y = element_text(size=font.size),
        strip.placement = "outside",
        axis.text.y = element_text(size=font.size-1),
        axis.text.x = element_text(size=font.size-1, angle=40, hjust=1),
        axis.ticks.x= element_line(color="black"),
        # axis.title.y = element_text(size=font.size+2),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=font.size),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        panel.spacing = unit(0.5, "lines"),
        legend.position ="none")

ggsave(filename="~/Desktop/rec_ssb_2021.png",
       width=6, height=8.5, units="in",dpi=600) 

# ggsave(filename=here::here(paste0(modelPath,models[1],"Figure_BASA ssb and rec.png")),
#        width=6, height=8.5, units="in",dpi=600) 