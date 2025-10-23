################################################################################

# plot mortality time series and posteriors

# plots mortality time series estimates and parameter posteriors

# authors: CL Roberts

# inputs: BASA model inputs (model/PWS_ASA.dat) and outputs (from mcmc_out/)

# outputs: 
#   - plots: 
#   - tables: 

################################################################################

#### front matter ####

# load packages

library(ggplot2)
theme_set(theme_bw())
# library(ggpubr)
library(pwsHerringBasa)
library(dplyr)
# library(tidyr)
library(here)
library(stringr)

# directory handling

dir_model <- here::here("model")

dir_mcmc_out <- here::here(dir_model, "mcmc_out")
dir_figures <- here::here("figures")
dir_outputs <- here::here("data_outputs")

if (!dir.exists(dir_figures)) {
    dir.create(dir_figures)
}

# read model input data

raw.data <- read.data.files(dir_model)

# save global variables
nyr <- raw.data$PWS_ASA.dat$nyr
nage <- raw.data$PWS_ASA.dat$nage
start.year <- 1980
curr.year <- start.year+nyr
years <- seq(start.year, curr.year-1)

# read model output data
incidence <- read.table(
    here(dir_mcmc_out, "VHSV_incidence.csv"), header = FALSE, 
    sep = ",", dec = ".", col.names = paste0(years), 
    check.names = FALSE, stringsAsFactors = FALSE
)

fatality <- read.table(
    here(dir_mcmc_out, "VHSV_fatalities.csv"), header = FALSE, 
    sep = ",", dec = ".", col.names = paste0(years), 
    check.names = FALSE, stringsAsFactors = FALSE
)

seroprevalence <- read.table(
    here(dir_mcmc_out, "VSHV_seroprevalence.csv"), header = FALSE, 
    sep = ",", dec = ".", col.names = paste0(years),
    check.names = FALSE, stringsAsFactors = FALSE
)

ssb <- read.table(
    here(dir_mcmc_out, "PFRBiomass.csv"), header = FALSE, 
    sep = ",", dec = ".", col.names = paste0(c(years,max(years)+1)),
    check.names = FALSE, stringsAsFactors = FALSE
)

age3 <- read.table(
    here(dir_mcmc_out, "Age3.csv"), header = FALSE, 
    sep = ",", dec = ".", col.names = paste0(years),
    check.names = FALSE, stringsAsFactors = FALSE
)


#### sereoprevalence survey fits

dat.plot.2 <- data.frame()

# Combine quantities for output
# USE PIVOT_LONGER INSTEAD
dat.plot <- rbind(
    data.frame(Variable="VHSV infection rate\n(Prop. of Spawners)",reshape2::melt(incidence,variable.name=c('Year'),factorsAsStrings = FALSE)),
    data.frame(Variable="VHSV mortality rate\n(Prop. of Spawners)",reshape2::melt(fatality,variable.name=c('Year'),factorsAsStrings = FALSE)),
    data.frame(Variable="Population immunity\n(Prop. of Spawners)",reshape2::melt(seroprevalence,variable.name=c('Year'),factorsAsStrings = FALSE)),
    data.frame(Variable="SSB (metric tons)",reshape2::melt(ssb,variable.name=c('Year'),factorsAsStrings = FALSE)),
    data.frame(Variable="Age 3 (millions)",reshape2::melt(age3,variable.name=c('Year'),factorsAsStrings = FALSE))
)

# Calc 95% quantiles for each year
dat.plot <- dat.plot |> 
    filter(Year %in% 2008:(max(years)-1)) |> 
    group_by(Variable,Year) |> 
    summarize(Q.025 = quantile(value,probs=0.025,na.rm=TRUE),
              Q.25 = quantile(value,probs=0.25,na.rm=TRUE),
              Q.50 = quantile(value,probs=0.5,na.rm=TRUE),
              Q.75 = quantile(value,probs=0.75,na.rm=TRUE),
              Q.975 = quantile(value,probs=0.975,na.rm=TRUE))

dat.plot.2 <- bind_rows(
    dat.plot.2, 
    data.frame(model = "seroprevalence", dat.plot)
)

# Estimate stats for paper
sum.stat <- dat.plot.2 |> 
    group_by(model, Variable) |> 
    summarise(avg.rates=median(Q.50),
              min.rates=min(Q.50),
              year.min=Year[which.min(Q.50)],
              max.rates=max(Q.50),
              year.max=Year[which.max(Q.50)])

trial <- dat.plot.2 |> 
    group_by(model,Year) |> 
    summarise(mor_rate_2=
              Q.50[Variable=='VHSV mortality rate\n(Prop. of Spawners)']/
              Q.50[Variable=='VHSV infection rate\n(Prop. of Spawners)'])
median(trial$mor_rate_2)

sum.stat.2 <- dat.plot.2 |> 
    group_by(Variable,Year) |> 
    summarise(max.diff=max(Q.50)-min(Q.50))

max.bounds <- dat.plot.2 |> 
    group_by(Variable) |> 
    summarise(max.val=max(Q.975))

sub_panel_labs <- data.frame(
    Year=0.5,
    Variable=max.bounds$Variable,
    value=max.bounds$max.val*0.95,
    label=LETTERS[1:5]
)
  
font.size = 11
seroprevalence_estimates_plot <- ggplot(data=dat.plot.2,aes(x=Year,y=Q.50,group=Variable)) + 
  geom_ribbon(aes(ymin=Q.025,ymax=Q.975),fill="grey70")+
  geom_ribbon(aes(ymin=Q.25,ymax=Q.75),fill="grey85")+
  geom_line(size=1)+
  scale_y_continuous(limits=c(0,NA))+
  facet_grid(Variable~model,switch="y",scales="free_y")+
  #scale_x_discrete(limits=c(0,50),expand=c(0,0))+
  geom_text(data=sub_panel_labs,aes(x=Year,y=value,label=label,
                                    group=Variable,size=font.size+3),hjust=0)+
  theme(strip.background = element_blank(),
        panel.border=element_rect(fill=NA),
        # plot.title = element_text(hjust = 0.5),
        plot.title = element_blank(),
        # strip.text.x = element_text(size=font.size),
        strip.text.x = element_blank(),
        strip.text.y = element_text(size=font.size),
        strip.placement = "outside",
        axis.text.y = element_text(size=font.size-1),
        axis.text.x = element_text(size=font.size-1,angle=40,hjust=1),
        axis.ticks.x= element_line(color="black"),
        # axis.title.y = element_text(size=font.size+2),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=font.size),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        panel.spacing = unit(0.5, "lines"),
        legend.position ="none")

ggsave(filename = here(dir_figures, "vsvs_outbreak_charcteristics.pdf"),
       width=6, height=8.5, units="in",dpi=600) 


###########################################################################
# Now assess posterior predictive fits to the seroprevalence data
###########################################################################  

rbinom_sero <- function(x){
  x <- as.numeric(x)
  x1 <- x[-(1:2)]
  pos_indices <- seq(1,length(x)-2,by=2)
  
  # X is a vector where each indice is a proportion of each age with or without antibodies
  sizes <- round(colSums(matrix(x1,nrow=2))*x[2])
  probs <- x1[pos_indices]/(x1[pos_indices] + x1[pos_indices+1])
  
  pos_sam <- apply(cbind(sizes,probs),1,function(x) rbinom(1,x[1],x[2]))
  neg_sam <- sizes-pos_sam
  
  y <- rep(NA,length(x1))
  y[pos_indices] <- pos_sam
  y[pos_indices+1] <- neg_sam
  return(y)
}

# Read in seroprevalence sample data
sero_obs <- raw.data$PWS_ASA_disease$vhsv_obs
sero_obs[sero_obs<0] <- NA
sample_size <- raw.data$agecomp_samp_sizes.txt$vhsv_sample_size
sample_size <- data.frame(years=years, ess=sample_size[,1])

# Create data frame for sero_obs
row.names(sero_obs) <- years
colnames(sero_obs) <- paste0(rep(c(paste0(0:8),"9+"),each=2),"_",rep(c("Present","Absent"),10))
sero_obs <- reshape2::melt(sero_obs) |> filter(!is.na(value))
names(sero_obs) <- c('years','age_spec_sero','obs')
sero_obs <- left_join(sero_obs,sample_size,by='years')


# Obtain effective sample size
ESS <- data.reader(here(dir_model, "agecomp_effective_sample_size.txt"))[[3]]
sero_samp_size <- data.frame(year=years,ss=ESS[,1]) |> 
    mutate(ss=ifelse(ss<0,NA,ss))

# Now check fits of seroprevalence to data
samp.pred.sero <- read.table(
    here(dir_mcmc_out, "VHSV.csv"), header = FALSE, 
    sep = ",", dec = ".", stringsAsFactors = FALSE
)

# Create key for the age-specific seroprevalence data
var.key <- expand.grid(1:nrow(samp.pred.sero),c("Present","Absent"),c(paste0(0:8),"9+"),years)
names(var.key) <- c('draw','antibody','age','years')

# Melt predicted sample seroprevalence & combine with key (year, age, antibodies, draw)
# CHANGE TO DPLYR::PIVOT_LONGER
samp.pred.sero.2 <- cbind(var.key,reshape2::melt(samp.pred.sero,factorsAsStrings = FALSE))

# Filter out years with no data
samp.pred.sero.2 <- samp.pred.sero.2 |> 
    filter(years>=2012) |>
    group_by(years) |>
    mutate(samp_size=sero_samp_size$ss[sero_samp_size$year%in%years])

# Make posterior predictive draws
samp.sero.post.pred <- NULL
for(i in min(samp.pred.sero.2$years):max(samp.pred.sero.2$years)) {

    indices <- which(samp.pred.sero.2$years==i)
    temp.df <- reshape2::dcast(samp.pred.sero.2[indices,],draw + samp_size ~ age + antibody)

    temp.df.2 <- reshape2::melt(t(apply(temp.df,1,function(x) rbinom_sero(x)/x[2])),id=c(1,2))

    temp.df.2 <- data.frame(samp.pred.sero.2[indices,-6],pp.value=temp.df.2$value) |> 
        group_by(draw,age,years) |>
        summarise(seroprev=pp.value[1]/(pp.value[1] + pp.value[2]),
                agecomp=pp.value[1] + pp.value[2])

    samp.sero.post.pred <- rbind(samp.sero.post.pred,temp.df.2)

    print(paste0(i," done"))

}

# Check for nan's
samp.sero.post.pred$seroprev[is.nan(samp.sero.post.pred$seroprev)] <- 0.0

# Melt so we just see the age comp and seroprevalence
samp.sero.post.pred <- reshape2::melt(samp.sero.post.pred,id.vars=1:3)

# Calculate 95% posterior predictive intervals
sero.fits <- samp.sero.post.pred |> 
    group_by(years,age,variable) |> 
    summarize(Q.025=quantile(value,probs=0.025,na.rm=TRUE),
                Q.50=quantile(value,probs=0.5,na.rm=TRUE),
                MEAN=mean(value,na.rm=TRUE),
                Q.975=quantile(value,probs=0.975,na.rm=TRUE))

sero.fits.2 <- sero.fits

# Create indices for cohorts
cohort.indices <- NULL
for(y in (min(sero.fits.2$years)-9):max(sero.fits.2$years)){
  cohort.indices <- rbind(cohort.indices,data.frame(Cohort=y,Age=c(paste0(0:8),"9+"),Year=y:(y+10-1)))
}

sero_obs$age <- str_extract_all(sero_obs$age_spec_sero, "\\d+|\\+") |>
    lapply(\(x) paste(x, collapse = "")) |>
    unlist()

sero_obs$antibody <- str_extract(sero_obs$age_spec_sero, "Absent|Present")

# Add separate columns for age and seropositivity
sero_obs_2 <- sero_obs |> 
    group_by(years,age) |>
    summarise(seroprev=obs[1]/(obs[1] + obs[2]),
                agecomp=obs[1] + obs[2],
                age_sampsize=round(unique(ess*agecomp)))
# This is quick & dirty for obtaining sample sizes, but I inspected and there
# are no rounding errors for raw sample sizes

# Check for nan's
sero_obs_2$seroprev[is.nan(sero_obs_2$seroprev)] <- 0.0

# Melt so we just see the age comp and seroprevalence
sero_obs_2 <- reshape2::melt(sero_obs_2,id.vars=c(1:2,5))


# Now add columns for cohort year to add color coding to plot
sero_obs_plot <- sero_obs_2 |> 
    filter(age!='0') |> 
    group_by(years, age) |> 
    # group_by(years, age, antibody) |> 
    mutate(cohort=cohort.indices$Cohort[cohort.indices$Year==years & cohort.indices$Age==age],
            variable=ifelse(variable=="seroprev","Seroprevalence","Age composition"))

sero_fits_plot <- sero.fits.2 |> 
    filter(age!='0') |> 
    group_by(years, age) |> 
    # group_by(years, age, antibody) |> 
    mutate(cohort=cohort.indices$Cohort[cohort.indices$Year==years & cohort.indices$Age==age],
            variable=ifelse(variable=="seroprev","Seroprevalence","Age composition"))

font.size <- 15

sero_fits_plot_2 <- sero_fits_plot |>
    filter(variable == "Seroprevalence", years %in% sero_obs_plot$years)
sero_obs_plot <- filter(sero_obs_plot,variable=="Seroprevalence")

year.labels <- data.frame(
    age = 1, 
    value = 0.85,
    years = sort(unique(sero_obs_plot$years)),
    variable = 'Seroprevalence'
)

# geom_text(data = sub.labs, mapping = aes(label=label.1, fill=NULL, color=NULL),x=max(mod_sel.df.2$y.axis.rank)-4, y=Inf,hjust=2,vjust=2,size=4.5) +
#   geom_text(aes(label = value), vjust = 0.5, hjust=-0.1,size=3.5)+

seroprevalence_fits_plot <- ggplot() + 
  geom_col(data=sero_obs_plot,aes(x=age, y=value, fill=factor(cohort)),position = 'identity')+
  # geom_pointrange(data=sero_fits_plot,aes(x=age,y=Q.50,ymin=Q.025,ymax=Q.975,color=model),
  #                 position=position_dodge(width=0.6),size=0.8,fatten=0.8) +
  geom_pointrange(data=sero_fits_plot_2,aes(x=age,y=MEAN,ymin=Q.025,ymax=Q.975,),
                  position=position_dodge(width=0.7),size=0.9,fatten=0.8) +
  # geom_col(aes(x=age, y=seropositive,group=year),fill="grey75")+
  # facet_wrap(year~.,dir="v",nrow=10)+
  # scale_fill_manual(values = c('age_comp'="black",'seropositive'='grey75'),
  #                   labels = c('- antibodies','+ antibodies')) + 
  scale_color_manual(values = c("black","grey55","grey35","grey75")) +
  scale_fill_manual(values = rep(c("olivedrab2","lightgoldenrod", "orange", "firebrick2", "darkmagenta", "blue", "cyan3"),length.out=diff(range(sero_obs_plot$cohort))+1)) +
  facet_grid(years~variable)+
  labs(y="Proportions")+
  xlab("Age")+
  geom_text(data=sero_obs_plot,aes(x=age, y=0.6, label = age_sampsize),
            hjust = -0.3, size=3.5,colour="grey40",fontface="italic")+
  geom_text(data=year.labels,aes(x=age,y=value,label=years,size=font.size))+
  scale_y_continuous(limits=c(0,1),breaks=c(0,0.5,1),labels=c(0,0.5,""))+
  guides(fill=FALSE,size=FALSE)+
  theme(strip.background = element_blank(),
        #strip.text.x = element_blank(),
        axis.line.x = element_blank(),
        strip.text.x =  element_text(size=font.size),
        strip.text.y = element_blank(),
        #panel.border=element_rect(fill=NA),
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_text(size=font.size),
        axis.title.x = element_text(size=font.size),
        axis.text.y = element_text(size=font.size-5,hjust=0.5,angle=90),
        #axis.text.y = element_blank(),
        axis.text.x = element_text(size=font.size-4),
        axis.ticks.x = element_blank(),
        #axis.ticks.x= element_line(color="black"),
        #axis.title.y = element_blank(),
        legend.position ="bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=font.size-4),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        panel.background = element_rect(fill="grey95"),
        panel.spacing.x = unit(1, "lines"),
        panel.spacing.y = unit(0, "lines"))

ggsave(filename = here(dir_figures, "seroprevalence_fit.pdf"),
       width = 6, height = 9, units = "in", dpi = 600)
ggsave(filename = here(dir_figures, "seroprevalence_fit.png"),
       width = 6, height = 9, units = "in")

# ggsave(filename=here::here(paste0("results/figures/Figure_basa seroprevalence estimates_comp with binom.png")),
#        width=9, height=9, units="in",dpi=600)

# ggsave(filename=here::here(paste0("results/figures/Figure_basa seroprevalence estimates.png")),
#        width=6.5, height=9, units="in",dpi=600)
# 
# ggsave(filename=here::here(paste0("results/figures/Figure_basa seroprevalence estimates_ESS 50.png")),
#        width=6.5, height=9, units="in",dpi=600)
# 
# ggsave(filename=here::here(paste0("results/figures/Figure_basa seroprevalence estimates_multinom vs binom.png")),
#        width=6.5, height=9, units="in",dpi=600)


###########################################################################
# Alternative visualization of seroprevalence data (as antibody + and - of total composition)
###########################################################################  
# Read in seroprevalence sample data

sero_obs <- raw.data$PWS_ASA_disease$vhsv_obs
sero_obs[sero_obs<0] <- NA

# Create data frame for sero_obs
row.names(sero_obs) <- years
colnames(sero_obs) <- paste0(rep(c(paste0(0:8),"9+"),each=2),"_",rep(c("Present","Absent"),10))
sero_obs <- reshape2::melt(sero_obs) |> filter(!is.na(value))
names(sero_obs) <- c('years','age_spec_sero','obs')

sero.fits.2 <- data.frame()

# Obtain effective sample size
ESS <- data.reader(here(dir_model, "agecomp_effective_sample_size.txt"))[[3]]
sero_samp_size <- data.frame(year=years,ss=ESS[,1]) |> 
    mutate(ss=ifelse(ss<0,NA,ss))

# Now check fits of seroprevalence to data
samp.pred.sero <- read.table(
    here(dir_mcmc_out, "VHSV.csv"), header = FALSE, 
    sep = ",", dec = ".", stringsAsFactors = FALSE
)

# Create key for the age-specific seroprevalence data
var.key <- expand.grid(1:nrow(samp.pred.sero),c("Present","Absent"),c(paste0(0:8),"9+"),years)
names(var.key) <- c('draw','antibody','age','years')

# Melt predicted sample seroprevalence & combine with key (year, age, antibodies, draw)
samp.pred.sero.2 <- cbind(var.key,reshape2::melt(samp.pred.sero,factorsAsStrings = FALSE))

# Filter out years with no data
samp.pred.sero.2 <- samp.pred.sero.2 |> 
    filter(years>=2012) |>
    group_by(years) |>
    mutate(samp_size=sero_samp_size$ss[sero_samp_size$year%in%years])

# Make posterior predictive draws
samp.sero.post.pred <- NULL
for(i in min(samp.pred.sero.2$years):max(samp.pred.sero.2$years)){
    indices <- which(samp.pred.sero.2$years==i)
    temp.df <- reshape2::dcast(samp.pred.sero.2[indices,],draw + samp_size ~ age + antibody)
    temp.df.2 <- reshape2::melt(t(apply(temp.df,1,function(x) rbinom_sero(x)/x[2])),id=c(1,2))
    samp.sero.post.pred <- rbind(samp.sero.post.pred,data.frame(samp.pred.sero.2[indices,1:5],pp.value=temp.df.2$value))

    print(paste0(i," done"))
}

sero.fits <- samp.sero.post.pred |> 
    group_by(years,age,antibody) |> 
    summarize(Q.025=quantile(pp.value,probs=0.025,na.rm=TRUE),
            Q.50=quantile(pp.value,probs=0.5,na.rm=TRUE),
            Q.975=quantile(pp.value,probs=0.975,na.rm=TRUE))

sero.fits <- sero.fits |> 
    mutate(age_spec_sero=paste0(age,"_",antibody))

sero.fits.2 <- sero.fits

# Create indices for cohorts
cohort.indices <- NULL
for(y in (min(sero.fits.2$years)-9):max(sero.fits.2$years)){
  cohort.indices <- rbind(cohort.indices,data.frame(Cohort=y,Age=c(paste0(0:8),"9+"),Year=y:(y+10-1)))
}

# Add separate columns for age and seropositivity
sero_obs <- sero_obs |> 
    group_by(years,age_spec_sero) |> 
    mutate(age=sero.fits.2$age[sero.fits.2$age_spec_sero==age_spec_sero][1],
            antibody=sero.fits.2$antibody[sero.fits.2$age_spec_sero==age_spec_sero][1])

# Now add columns for cohort year to add color coding to plot
sero_obs_plot <- sero_obs |> 
    filter(age!='0') |> 
    group_by(years,age) |> 
    mutate(cohort=cohort.indices$Cohort[cohort.indices$Year==years & cohort.indices$Age==age])

sero_fits_plot <- sero.fits.2 |> 
    filter(age!='0', years %in% sero_obs_plot$years) |> 
    group_by(years,age) |> 
    mutate(cohort=cohort.indices$Cohort[cohort.indices$Year==years & cohort.indices$Age==age])

font.size <- 15

year.labels <- data.frame(
    age = 8.5,
    obs = 0.85,
    years = sort(unique(sero_obs_plot$years)),
    antibody = 'Absent'
)

seroprevalence_fits_plot2 <- ggplot() + 
  geom_col(data=sero_obs_plot,aes(x=age, y=obs, fill=factor(cohort)),position = 'identity')+
  geom_pointrange(data=sero_fits_plot,aes(x=age,y=Q.50,ymin=Q.025,ymax=Q.975),
                  position=position_dodge(width=0.6),size=0.8,fatten=0.8) +
  # geom_col(aes(x=age, y=seropositive,group=year),fill="grey75")+
  # facet_wrap(year~.,dir="v",nrow=10)+
  # scale_fill_manual(values = c('age_comp'="black",'seropositive'='grey75'),
  #                   labels = c('- antibodies','+ antibodies')) + 
  scale_color_manual(values = c("black","grey70","grey45")) +
  scale_fill_manual(values = rep(c("olivedrab2","lightgoldenrod", "orange", "firebrick2", "darkmagenta", "blue", "cyan3"),length.out=diff(range(sero_obs_plot$cohort))+1)) +
  facet_grid(years~antibody)+
  labs(y="Age-specific antibody prevalence")+
  xlab("Age")+
  geom_text(data=year.labels,aes(x=age,y=obs,label=years,size=font.size))+
  scale_y_continuous(limits=c(0,1),breaks=c(0,0.5,1),labels=c(0,0.5,1))+
  theme(strip.background = element_blank(),
        #strip.text.x = element_blank(),
        axis.line.x = element_blank(),
        strip.text.x =  element_text(size=font.size),
        strip.text.y = element_blank(),
        #panel.border=element_rect(fill=NA),
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_text(size=font.size),
        axis.title.x = element_text(size=font.size),
        axis.text.y = element_text(size=font.size-5,hjust=0.5,angle=90),
        #axis.text.y = element_blank(),
        axis.text.x = element_text(size=font.size-4),
        axis.ticks.x = element_blank(),
        #axis.ticks.x= element_line(color="black"),
        #axis.title.y = element_blank(),
        legend.position ="none",
        legend.title = element_blank(),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        panel.background = element_rect(fill="grey95"),
        panel.spacing.x = unit(1, "lines"),
        panel.spacing.y = unit(0.5, "lines"))

ggsave(filename=here(dir_figures, "seroprevalence_fit2.pdf"),
       width=6.5, height=9, units="in",dpi=600)  
ggsave(filename=here(dir_figures, "seroprevalence_fit2.png"),
       width=6.5, height=9, units="in")  
