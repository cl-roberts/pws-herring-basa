library(here)
library(ggplot2)
library(tidyverse)
library(ggdist)

source(file=paste0(here::here("functions/"), "fun_read_dat.R"))

generate.post.pred <- function(seine.age.comp, spawn.age.comp, seine.ess, spawn.ess, ncol=10){
    pp.seine.age.comp<-matrix(0, nrow(seine.age.comp), ncol(seine.age.comp)) # CHANGE
    pp.spawn.age.comp<-matrix(0, nrow(spawn.age.comp), ncol(spawn.age.comp))
    # sample using the MCMC draws 
    #set.seed(100)
    for(i in 1:nrow(seine.age.comp)){ # Loop through the MCMC draws
        for(j in 1:length(seine.ess)){ # Loop through each year
            pp.seine.age.comp[i, (j*ncol-(ncol-1)):(j*ncol) ] <- t(rmultinom(1, size=as.integer(seine.ess[j]),seine.age.comp[i, (j*ncol-(ncol-1)):(j*ncol)])) / seine.ess[j]*100
            if(all(spawn.age.comp[i, (j*ncol-(ncol-1)):(j*ncol)] == 0)){
                pp.spawn.age.comp[i, (j*ncol-(ncol-1)):(j*ncol)] <- rep(0, times=ncol)
            }else{
                pp.spawn.age.comp[i, (j*ncol-(ncol-1)):(j*ncol) ] <- t(rmultinom(1, spawn.ess[j], spawn.age.comp[i, (j*ncol-(ncol-1)):(j*ncol)])) / spawn.ess[j] *100
            }
        }
    }
    return(list(
        seine.pp = pp.seine.age.comp,
        spawn.pp = pp.spawn.age.comp
    ))
}

color.options <- c("#13F24AFF", "#FFBA00FF", "firebrick1", "darkmagenta", "blue", "#00A4FFFF")
generate.colors <- function(n){

    shifter <- function(x, n = 1) {
        return(if (n == 0) x else c(tail(x, -n), head(x, n)))
    }

    colors <- rep(NA, n)
    i=1
    j=0
    cs <- color.options
    while(i < n){
        cs <- shifter(cs, length(color.options)-1)
        for(c in cs){
            colors[i] <- c
            colors[i] <- c
            i = i+1
        }
        colors[i] <- "gray"
        colors[i] <- "gray"
        i = i+1
        j = j+1
    }
    return(colors)
}

start.year <- 1980
curr.year <- 2024
nyr.sim <- 0
years <- seq(start.year, curr.year+nyr.sim-1)
nyr <- length(years)

model.dir <- here::here("model/")

raw.data <- read.data.files(model.dir)
model.data <- list(nyr=raw.data$PWS_ASA.dat$nyr,
                   nage=raw.data$PWS_ASA.dat$nage,
                   ess.seine=raw.data$PWS_ASA_ESS.ctl$seine_ess[1:nyr],
                   ess.spawn=raw.data$PWS_ASA_ESS.ctl$spawn_ess[1:nyr],
                   seac=raw.data$PWS_ASA.dat$seine_age_comp[1:nyr,]*100,
                   spac=raw.data$PWS_ASA.dat$spawn_age_comp[1:nyr,]*100,
                   seine_indices=which(rowSums(raw.data$PWS_ASA.dat$seine_age_comp[1:nyr,])>0),
                   spawn_indices=which(rowSums(raw.data$PWS_ASA.dat$spawn_age_comp[1:nyr,])>0)
                   )
model.data$seac[model.data$seac == -900] <- 0
model.data$spac[model.data$spac == -900] <- 0
model.data$ess.seine[model.data$ess.seine == -9] <- 0
model.data$ess.spawn[model.data$ess.spawn == -9] <- 0

colnames(model.data$seac) <- seq(0, 9, 1)
colnames(model.data$spac) <- seq(0, 9, 1)

rownames(model.data$seac) <- years
rownames(model.data$spac) <- years

colors <- generate.colors(n=nyr*7)

seac.raw.df <- as_tibble(model.data$seac) %>%
                mutate(year=years) %>%
                pivot_longer(matches("[0-9]"), names_to="age", values_to="val") %>%
                filter(age > 2) %>%
                mutate(
                    fill.color=colors,
                    type="seine"
                ) %>%
                print(n=30)
spac.raw.df <- as_tibble(model.data$spac) %>%
                mutate(year=years) %>%
                pivot_longer(matches("[0-9]"), names_to="age", values_to="val") %>%
                filter(age > 2) %>%
                mutate(
                    fill.color=colors,
                    type="spawn"
                ) %>%
                print(n=30)
raw.df <- rbind(seac.raw.df, spac.raw.df)

spawn.age.comp<-read.csv(paste0(model.dir, "mcmc_out/SpAC.csv"), header = FALSE, dec=".") 
spawn.age.comp<-spawn.age.comp[-c(1:nburn), 1:(nyr*10)]*100

seine.age.comp<-read.csv(paste0(model.dir, "mcmc_out/SeAC.csv"), header = FALSE, dec=".") 
seine.age.comp<-seine.age.comp[-c(1:nburn), 1:(nyr*10)]*100

pps <- generate.post.pred(seine.age.comp, spawn.age.comp, model.data$ess.seine, model.data$ess.spawn)

seine.age.comp.quants <- apply(pps$seine.pp, 2, quantile, probs=c(0.025,0.5,0.975), na.rm=T)
spawn.age.comp.quants <- apply(pps$spawn.pp, 2, quantile, probs=c(0.025,0.5,0.975), na.rm=T)

cnames <- rep(NA, ncol(seine.age.comp))
i=1
for(y in years){
    for(a in 0:9){
        cnames[i] <- paste0(y, "_", a)
        i = i+1
    }
}
colnames(seine.age.comp.quants) <- cnames
colnames(spawn.age.comp.quants) <- cnames

seac.df <- as_tibble(seine.age.comp.quants) %>%
            pivot_longer(everything(), names_to="year_age", values_to="val") %>%
            separate(year_age, c("year", "age"), sep="_") %>%
            mutate(percentile=rep(c("2.5%", "50%", "97.5%"), each=(10*nyr))) %>%
            pivot_wider(names_from=percentile, values_from=val) %>%
            filter(age > 2) %>%
            mutate(
                type="seine"
            ) %>%
            print(n=30)

spac.df <- as_tibble(spawn.age.comp.quants) %>%
            pivot_longer(everything(), names_to="year_age", values_to="val") %>%
            separate(year_age, c("year", "age"), sep="_") %>%
            mutate(percentile=rep(c("2.5%", "50%", "97.5%"), each=(10*nyr))) %>%
            pivot_wider(names_from=percentile, values_from=val) %>%
            filter(age > 2) %>%
            mutate(
                type="spawn"
            ) %>%
            print(n=30)

age.comp.df <- rbind(seac.df, spac.df)
age.comp.df$age <- factor(age.comp.df$age, levels=sort(unique(age.comp.df$age)), ordered=TRUE)
age.comp.df$type <- factor(age.comp.df$type)

age.class.df <- data.frame(year=raw.df$year, age=rep(c("3", "4", "5", "6", "7", "8", "9"), nrow(raw.df)), color=rep(c("black", "white", "black", "white", "black", "white", "black"), nrow(raw.df)), type=rep(c("seine", "spawn"), each=nrow(raw.df)/2))
year.df <- data.frame(year=years, year_str=years)

age.struct.plot <- ggplot(raw.df)+
    geom_col(aes(x=type, y=val/100, color=age, fill=fill.color), position=position_dodge(0.9), size=0.0)+
    scale_fill_manual(values=c(color.options, "grey")) + 
    geom_pointinterval(data=age.comp.df, aes(x=type, y=`50%`/100, ymin=`2.5%`/100, ymax=`97.5%`/100, color=age), position=position_dodge(0.9)) +
    geom_text(data=year.df, aes(x=2.3, y=1, label=year_str), size=4)+
    geom_text(data=age.class.df, aes(x=type, y=-0.25, label=age, group=age), position=position_dodge(0.9), size=3)+
    scale_color_manual(values=rep("black", 7))+
    scale_y_continuous("Proportion", breaks=c(0.0, 0.5, 1.0))+
    scale_x_discrete(labels=c("Seine Catch", "Spawner Survey"))+
    facet_wrap(~year, nrow=12, dir="v")+
    coord_cartesian(clip="off", ylim=c(0, 1.15))+
    labs(x="Age Class", y="Proportion", title="Proportional Age Structure")+
    theme(
        legend.position = "none",
        panel.grid = element_blank(),
        panel.spacing.y = unit(0.1, "line"),
        panel.spacing.x = unit(0.25, "line"),
        axis.text.x = element_text(margin=margin(t=12)),
        axis.ticks.x = element_blank(),
        axis.title.x = element_text(margin=margin(t=10)),
        strip.background = element_blank(),
        strip.text.x = element_blank()
    )
age.struct.plot

ggsave(paste0(here::here("figures/"), "age_compositions.pdf"), height = 11, width=8.5)

###############################################################
## Save csv of values displayed in table above (with 95% posterior predictive intervals) 

age.comp.median.wide <- age.comp.df %>%
                        select(year, age, `50%`, type) %>%
                        mutate(`50%` = round(`50%`/100, 3)) %>%
                        pivot_wider(
                            everything(),
                            names_from=age,
                            values_from=`50%`
                        ) %>%
                        print(n=80)

age.comp.lower.wide <- age.comp.df %>%
                        select(year, age, `2.5%`, type) %>%
                        mutate(`2.5%` = round(`2.5%`/100, 3)) %>%
                        pivot_wider(
                            everything(),
                            names_from=age,
                            values_from=`2.5%`
                        ) %>%
                        print(n=80)

age.comp.upper.wide <- age.comp.df %>%
                        select(year, age, `97.5%`, type) %>%
                        mutate(`97.5%` = round(`97.5%`/100, 3)) %>%
                        pivot_wider(
                            everything(),
                            names_from=age,
                            values_from=`97.5%`
                        ) %>%
                        print(n=80)

sink(here::here("data_outputs/predicted-age-comps.csv"))

cat("# Seine Fishery Age Composition (proportion of each age): Median Posterior Prediction\n")
write.table(
    age.comp.median.wide %>% filter(type == "seine") %>% select(-c(type)), 
    row.names=FALSE, col.names=c("# Year", 3:8, "9+"), sep=","
)

cat("\n# Seine Fishery Age Composition (proportion of each age): Lower 95% Posterior Predictive Interval\n")
write.table(
    age.comp.lower.wide %>% filter(type == "seine") %>% select(-c(type)),  
    row.names=FALSE, col.names=c("# Year", 3:8, "9+"), sep=","
)

cat("\n# Seine Fishery Age Composition (proportion of each age): Upper 95% Posterior Predictive Interval\n")
write.table(
    age.comp.upper.wide %>% filter(type == "seine") %>% select(-c(type)),  
    row.names=FALSE, col.names=c("# Year", 3:8, "9+"), sep=","
)

cat("\n\n# Spawner Survey Age Composition (proportion of each age): Median Posterior Prediction\n")
write.table(
    age.comp.median.wide %>% filter(type == "spawn") %>% select(-c(type)), 
    row.names=FALSE, col.names=c("# Year", 3:8, "9+"), sep=","
)
cat("\n# Spawner Survey Age Composition (proportion of each age): Lower 95% Posterior Predictive Interval\n")
write.table(
    age.comp.lower.wide %>% filter(type == "spawn") %>% select(-c(type)),  
    row.names=FALSE, col.names=c("# Year", 3:8, "9+"), sep=","
)
cat("\n# Spawner Survey Age Composition (proportion of each age): Upper 95% Posterior Predictive Interval\n")
write.table(
    age.comp.upper.wide %>% filter(type == "spawn") %>% select(-c(type)),  
    row.names=FALSE, col.names=c("# Year", 3:8, "9+"), sep=","
)
sink()

# TODO: tidyverse-ify this
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
