<<<<<<< HEAD
################################################################################

# plot age structure data

# plots BASA-estimated numbers-at-age and age compositions fit to seine fishery 
# and spawn survey data for PWS herring 1980-present

# authors: John Trochta, Joshua Zahner, CL Roberts

# inputs: BASA model inputs (model/PWS_ASA.dat) and outputs (from mcmc_out/)

# outputs: 
#   - plots: 
#       - age compositions 1980-present (figures/age_compositions.pdf)
#       - numbers-at-age (figures/numbers_at_age.pdf)
#   - tables: 
#       - age compositions 1980-present (data_outputs/predicted-age-comps.csv)
#       - numbers-at-age (data_outputs/numbers-at-age.csv)

################################################################################

#### front matter ####

# choose TMB or ADMB
software <- "ADMB"

# load packages

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggdist)
library(pwsHerringBasa)
library(data.table)
library(ggridges)

# directory handling

dir_model <- here::here("model")

if (software == "ADMB") {
    dir_mcmc_out <- here::here(dir_model, "mcmc_out")
    dir_figures <- here::here("figures")
    dir_outputs <- here::here("data_outputs")
} else if (software == "TMB") {
    dir_mcmc_out <- here::here(dir_model, "mcmc_out_tmb")
    dir_figures <- here::here("figures/tmb")
    dir_outputs <- here::here("data_outputs/tmb")
} else {
    stop("choose valid software")
}

if (!dir.exists(dir_figures)) {
    dir.create(dir_figures)
}

# read in model input data 
raw.data <- read.data.files(dir_model)

# save global variables

nyr <- raw.data$PWS_ASA.dat$nyr
start.year <- 1980
curr.year <- start.year+nyr
nyr.sim <- 0
years <- seq(start.year, curr.year+nyr.sim-1)
nburn <- 1  # CLR: added this line
ncol <- 10   # CLR: added this line
# color.options <- RColorBrewer::brewer.pal(7, "Set2")
# color.options <- RColorBrewer::brewer.pal(7, "Dark2")
color.options <- gplots::rich.colors(n = 7, alpha = .5)
colors <- generate.colors(nyr = nyr, color.options = color.options[1:6])

#-------------------------------------------------------------------------------

#### format data ####


# format model input data

model.data <- list(nyr=raw.data$PWS_ASA.dat$nyr,
                   nage=raw.data$PWS_ASA.dat$nage,
                   ess.seine=raw.data$PWS_ASA_ESS.ctl$seine_ess[1:nyr],
                   ess.spawn=raw.data$PWS_ASA_ESS.ctl$spawn_ess[1:nyr],
                   seac=raw.data$PWS_ASA.dat$seine_age_comp[1:nyr,]*100,
                   spac=raw.data$PWS_ASA.dat$spawn_age_comp[1:nyr,]*100,
                   seine_indices=which(rowSums(raw.data$PWS_ASA.dat$seine_age_comp[1:nyr,])>0),
                   spawn_indices=which(rowSums(raw.data$PWS_ASA.dat$spawn_age_comp[1:nyr,])>0)
                   )

# replace missing data with 0's

model.data$seac[model.data$seac == -900] <- 0
model.data$spac[model.data$spac == -900] <- 0
model.data$ess.seine[model.data$ess.seine == -9] <- 0  
model.data$ess.spawn[model.data$ess.spawn == -9] <- 0

colnames(model.data$seac) <- seq(0, 9, 1)
colnames(model.data$spac) <- seq(0, 9, 1)
rownames(model.data$seac) <- years
rownames(model.data$spac) <- years

# tidy-ify data

seac.raw.df <- as_tibble(model.data$seac) |>
                mutate(year=years) |>
                pivot_longer(matches("[0-9]"), names_to="age", values_to="val") |>
                filter(age > 2) |>
                mutate(
                    fill.color=colors,
                    type="seine"
                ) |>
                print(n=30)
spac.raw.df <- as_tibble(model.data$spac) |>
                mutate(year=years) |>
                pivot_longer(matches("[0-9]"), names_to="age", values_to="val") |>
                filter(age > 2) |>
                mutate(
                    fill.color=colors,
                    type="spawn"
                ) |>
                print(n=30)
raw.df <- rbind(seac.raw.df, spac.raw.df)


#-------------------------------------------------------------------------------

#### read in model fit to age comps, calculate posterior intervals ####

# save locations for model age comp fits 

# CLR: added / to file paths
spawn.age.comp.fname <- here::here(dir_mcmc_out, "SpAC.csv")
seine.age.comp.fname <- here::here(dir_mcmc_out, "SeAC.csv")

# save fits and calculate posterior predictive intervals
# see ?pwsHerringBasa::generate.post.pred() for more info

pp.spawn.age.comp <- generate.post.pred(spawn.age.comp.fname, ess = model.data$ess.spawn, 
                                        years = years, nburn = nburn, dist = "multinom")
pp.seine.age.comp <- generate.post.pred(seine.age.comp.fname, ess = model.data$ess.seine, 
                                        years = years, nburn = nburn, dist = "multinom")

# save and format median and 95% credible intervals

seine.age.comp.quants <- apply(pp.seine.age.comp, 2, quantile, probs=c(0.025,0.5,0.975), na.rm=T)
spawn.age.comp.quants <- apply(pp.spawn.age.comp, 2, quantile, probs=c(0.025,0.5,0.975), na.rm=T)

cnames <- rep(NA, ncol(pp.spawn.age.comp))
i=1
for(y in years){
    for(a in 0:9){
        cnames[i] <- paste0(y, "_", a)
        i = i+1
    }
}
colnames(seine.age.comp.quants) <- cnames
colnames(spawn.age.comp.quants) <- cnames

# tidy-ify median age comps and credible intervals for plotting

seac.df <- as_tibble(seine.age.comp.quants) |>
            pivot_longer(everything(), names_to="year_age", values_to="val") |>
            separate(year_age, c("year", "age"), sep="_") |>
            mutate(percentile=rep(c("2.5%", "50%", "97.5%"), each=(10*nyr))) |>
            pivot_wider(names_from=percentile, values_from=val) |>
            filter(age > 2) |>
            mutate(
                type="seine"
            ) |>
            print(n=30)

spac.df <- as_tibble(spawn.age.comp.quants) |>
            pivot_longer(everything(), names_to="year_age", values_to="val") |>
            separate(year_age, c("year", "age"), sep="_") |>
            mutate(percentile=rep(c("2.5%", "50%", "97.5%"), each=(10*nyr))) |>
            pivot_wider(names_from=percentile, values_from=val) |>
            filter(age > 2) |>
            mutate(
                type="spawn"
            ) |>
            print(n=30)

age.comp.df <- rbind(seac.df, spac.df)
age.comp.df$age <- factor(age.comp.df$age, levels=sort(unique(age.comp.df$age)), ordered=TRUE)
age.comp.df$type <- factor(age.comp.df$type)

age.class.df <- data.frame(year=raw.df$year, age=rep(c("3", "4", "5", "6", "7", "8", "9"), nrow(raw.df)), color=rep(c("black", "white", "black", "white", "black", "white", "black"), nrow(raw.df)), type=rep(c("seine", "spawn"), each=nrow(raw.df)/2))
year.df <- data.frame(year=years, year_str=years)


#-------------------------------------------------------------------------------

#### read in format model-estimated numbers-at-age ####

# read-in Numbers-at-age and output into nice tables (.csv) 

n.y.a<-read.csv(here::here(dir_mcmc_out, "Num_at_age.csv"), header = FALSE, dec=".") 
n.y.a<-n.y.a[-c(1:nburn), ] # Just 'cause it's easier to work with this scale

## seine.age.comp.interval is a 2 X (# of ages in age comp * years) matrix; e.g. 10 ages (ages 0-9+) * 35 years (from 1980-2014)=350
n.y.a.interval <- matrix(rep(0, 2*length(n.y.a[1, ])), nrow=2, ncol=length(n.y.a[1, ]))
n.y.a.interval[1:2,] <- apply(n.y.a, 2, quantile, probs=c(0.025, 0.975), na.rm=T) # Fills first row with 0.025 percentile and second with 0.975 percentile
n.y.a.median <- apply(n.y.a, 2, median)
## Now cut each up into years

n.y.a.median.mat <- matrix(n.y.a.median, nrow=nyr, ncol=ncol, byrow=T) # nyr x nages matrix filled with posterior predictive median
n.y.a.lower.mat <- matrix(n.y.a.interval[1, ], nrow=nyr, ncol=ncol, byrow=T)
n.y.a.upper.mat <- matrix(n.y.a.interval[2, ], nrow=nyr, ncol=ncol, byrow=T)

# calculate shannon weiner evenness index (J) for spawner age comps

spac.50 <- spac.df |> 
            select(year, age, `50%`) |>
            filter(!is.na(`50%`)) |>
            pivot_wider(names_from=age, values_from=`50%`) |>
            select(c(`3`, `4`, `5`, `6`, `7`, `8`, `9`)) |>
            as.matrix()
evenness.50 <- apply(spac.50, 1, sw_evenness)

# format numbers-at-age for ggridges plot

naa_ridgedat <- n.y.a.median.mat |>
    as.data.table() |>
    mutate(Year = years) |>
    melt(id.vars = "Year", variable.name = "Age", value.name = "Numbers") |>
    mutate(Age = as.numeric(Age), Numbers = round(Numbers)) |> 
    uncount(Numbers)


#-------------------------------------------------------------------------------

#### make plots ####

# make age composition plot

raw.df <- raw.df |> 
    left_join(data.frame(year = 1982:(curr.year+nyr.sim-1), J = round(evenness.50, 2)))

age.struct.plot <- ggplot(raw.df)+
    geom_col(aes(x=type, y=val/100, color=age, fill=fill.color), position=position_dodge(0.9), linewidth=0)+
    scale_fill_manual(values=color.options) + 
    geom_pointinterval(data=age.comp.df, aes(x=type, y=`50%`/100, ymin=`2.5%`/100, ymax=`97.5%`/100, color=age), position=position_dodge(0.9)) +
    geom_text(data=year.df, aes(x=0.65, y=1, label=year), size=4)+
    geom_text(aes(x=2.25, y=1, label=ifelse(is.na(J), NA, paste0("J=", J))), size=4)+
    geom_text(data=age.class.df, aes(x=type, y=-0.25, label=age, group=age), position=position_dodge(0.9), size=3)+
    scale_color_manual(values=rep("black", 7))+
    scale_y_continuous("Proportion", breaks=c(0.0, 0.5, 1.0))+
    scale_x_discrete(labels=c("Seine Catch", "Spawner Survey"))+
    facet_wrap(~year, nrow=12, dir="v")+
    coord_cartesian(clip="off", ylim=c(0, 1.15))+
    labs(x="Age Class", y="Proportion", title="Proportional Age Structure")+
    theme_bw() +
    theme(
        legend.position = "none",
        panel.border = element_rect(linewidth = 0),
        panel.grid = element_blank(),
        panel.spacing.y = unit(0.0, "line"),
        panel.spacing.x = unit(0.25, "line"),
        axis.text.x = element_text(margin=margin(t=12)),
        axis.ticks.x = element_blank(),
        axis.title.x = element_text(margin=margin(t=10)),
        strip.background = element_blank(),
        strip.text.x = element_blank()
    )

# make numbers-at-age plot

# naa_plot <- ggplot(naa_ridgedat, aes(x = Age, y = Year, fill = Year, group = Year)) +
#     geom_density_ridges2() +
#     scale_fill_viridis_c(option = viridis_palette, trans = "reverse") + 
#     scale_y_reverse() +
#     xlim(c(0, 12)) +
#     labs(title = "Estimated Numbers-at-age", subtitle = "Prince William Sound Herring, 1980-Present") +
#     theme_bw()


#-------------------------------------------------------------------------------

#### save outputs ####

# save age comps plot

ggsave(here::here(dir_figures, "age_compositions.pdf"), plot = age.struct.plot, height = 11, width=8.5)

# save numbers-at-age plot

# ggsave(here::here(dir_figures, "numbers_at_age.pdf"), plot = naa_plot, height = 11, width=8.5)

# save csv table of median age comps with 95% posterior predictive intervals

age.comp.median.wide <- age.comp.df |>
                        select(year, age, `50%`, type) |>
                        mutate(`50%` = round(`50%`/100, 3)) |>
                        pivot_wider(
                            everything(),
                            names_from=age,
                            values_from=`50%`
                        ) |>
                        print(n=80)

age.comp.lower.wide <- age.comp.df |>
                        select(year, age, `2.5%`, type) |>
                        mutate(`2.5%` = round(`2.5%`/100, 3)) |>
                        pivot_wider(
                            everything(),
                            names_from=age,
                            values_from=`2.5%`
                        ) |>
                        print(n=80)

age.comp.upper.wide <- age.comp.df |>
                        select(year, age, `97.5%`, type) |>
                        mutate(`97.5%` = round(`97.5%`/100, 3)) |>
                        pivot_wider(
                            everything(),
                            names_from=age,
                            values_from=`97.5%`
                        ) |>
                        print(n=80)

sink(here::here(dir_outputs, "predicted-age-comps.csv"))

cat("# Seine Fishery Age Composition (proportion of each age): Median Posterior Prediction\n")
write.table(
    age.comp.median.wide |> filter(type == "seine") |> select(-c(type)), 
    row.names=FALSE, col.names=c("# Year", 3:8, "9+"), sep=","
)

cat("\n# Seine Fishery Age Composition (proportion of each age): Lower 95% Posterior Predictive Interval\n")
write.table(
    age.comp.lower.wide |> filter(type == "seine") |> select(-c(type)),  
    row.names=FALSE, col.names=c("# Year", 3:8, "9+"), sep=","
)

cat("\n# Seine Fishery Age Composition (proportion of each age): Upper 95% Posterior Predictive Interval\n")
write.table(
    age.comp.upper.wide |> filter(type == "seine") |> select(-c(type)),  
    row.names=FALSE, col.names=c("# Year", 3:8, "9+"), sep=","
)

cat("\n\n# Spawner Survey Age Composition (proportion of each age): Median Posterior Prediction\n")
write.table(
    age.comp.median.wide |> filter(type == "spawn") |> select(-c(type)), 
    row.names=FALSE, col.names=c("# Year", 3:8, "9+"), sep=","
)
cat("\n# Spawner Survey Age Composition (proportion of each age): Lower 95% Posterior Predictive Interval\n")
write.table(
    age.comp.lower.wide |> filter(type == "spawn") |> select(-c(type)),  
    row.names=FALSE, col.names=c("# Year", 3:8, "9+"), sep=","
)
cat("\n# Spawner Survey Age Composition (proportion of each age): Upper 95% Posterior Predictive Interval\n")
write.table(
    age.comp.upper.wide |> filter(type == "spawn") |> select(-c(type)),  
    row.names=FALSE, col.names=c("# Year", 3:8, "9+"), sep=","
)
sink()


# save csv table of numbers-at-age

# CLR: changed long.years to years and round_df to round in the following write.table() calls

sink(here::here(dir_outputs, "numbers-at-age.csv"))
cat("# Numbers-at-age (millions): Median estimate\n")
write.table(
    cbind(years, round(n.y.a.median.mat, 3), round(apply(n.y.a.median.mat, 1, sw_evenness), 2)), 
    row.names=FALSE, col.names=c("# Year", 0:8, "9+", "J"), sep=","
)

cat("\n# Numbers-at-age (millions): Lower 95% Credibility Interval\n")
write.table(
    cbind(years, round(n.y.a.lower.mat, 3), round(apply(n.y.a.lower.mat, 1, sw_evenness), 2)),  
    row.names=FALSE, col.names=c("# Year", 0:8, "9+", "J"), sep=","
)

cat("\n# Numbers-at-age (millions): Upper 95% Credibility Interval\n")
write.table(
    cbind(years, round(n.y.a.upper.mat, 3), round(apply(n.y.a.upper.mat, 1, sw_evenness), 2)),  
    row.names=FALSE, col.names=c("# Year", 0:8, "9+", "J"), sep=","
)

sink()

=======
################################################################################

# plot age structure data

# plots BASA-estimated numbers-at-age and age compositions fit to seine fishery 
# and spawn survey data for PWS herring 1980-present

# authors: John Trochta, Joshua Zahner, CL Roberts

# inputs: BASA model inputs (model/PWS_ASA.dat) and outputs (from mcmc_out/)

# outputs: 
#   - plots: 
#       - age compositions 1980-present (figures/age_compositions.pdf)
#       - numbers-at-age (figures/numbers_at_age.pdf)
#   - tables: 
#       - age compositions 1980-present (data_outputs/predicted-age-comps.csv)
#       - numbers-at-age (data_outputs/numbers-at-age.csv)

################################################################################

#### front matter ####

# choose TMB or ADMB
software <- "ADMB"

# load packages

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggdist)
library(pwsHerringBasa)
library(data.table)
library(ggridges)
library(gplots)

# directory handling

dir_model <- here::here("model")

if (software == "ADMB") {
    dir_mcmc_out <- here::here(dir_model, "mcmc_out")
    dir_figures <- here::here("figures")
    dir_outputs <- here::here("data_outputs")
} else if (software == "TMB") {
    dir_mcmc_out <- here::here(dir_model, "mcmc_out_tmb")
    dir_figures <- here::here("figures/tmb")
    dir_outputs <- here::here("data_outputs/tmb")
} else {
    stop("choose valid software")
}

if (!dir.exists(dir_figures)) {
    dir.create(dir_figures)
}

if (!dir.exists(dir_outputs)) {
    dir.create(dir_outputs)
}

# read in model input data 
raw.data <- read.data.files(dir_model)

# save global variables

nyr <- raw.data$PWS_ASA.dat$nyr
start.year <- 1980
curr.year <- start.year+nyr
nyr.sim <- 0
years <- seq(start.year, curr.year+nyr.sim-1)
nburn <- 1  
ncol <- 10   
color.options <- rich.colors(n = 7, alpha = .5)
colors <- generate.colors(nyr = nyr, color.options = color.options[1:6])

#-------------------------------------------------------------------------------

#### format data ####


# format model input data

model.data <- list(nyr=raw.data$PWS_ASA.dat$nyr,
                   nage=raw.data$PWS_ASA.dat$nage,
                   ess.seine=raw.data$PWS_ASA_ESS.ctl$seine_ess[1:nyr],
                   ess.spawn=raw.data$PWS_ASA_ESS.ctl$spawn_ess[1:nyr],
                   seac=raw.data$PWS_ASA.dat$seine_age_comp[1:nyr,]*100,
                   spac=raw.data$PWS_ASA.dat$spawn_age_comp[1:nyr,]*100,
                   seine_indices=which(rowSums(raw.data$PWS_ASA.dat$seine_age_comp[1:nyr,])>0),
                   spawn_indices=which(rowSums(raw.data$PWS_ASA.dat$spawn_age_comp[1:nyr,])>0)
                   )

# replace missing data with 0's

model.data$seac[model.data$seac == -900] <- 0
model.data$spac[model.data$spac == -900] <- 0
model.data$ess.seine[model.data$ess.seine == -9] <- 0  
model.data$ess.spawn[model.data$ess.spawn == -9] <- 0

colnames(model.data$seac) <- seq(0, 9, 1)
colnames(model.data$spac) <- seq(0, 9, 1)
rownames(model.data$seac) <- years
rownames(model.data$spac) <- years

# tidy-ify data

seac.raw.df <- as_tibble(model.data$seac) |>
                mutate(year=years) |>
                pivot_longer(matches("[0-9]"), names_to="age", values_to="val") |>
                filter(age > 2) |>
                mutate(
                    fill.color=colors,
                    type="seine"
                ) |>
                print(n=30)
spac.raw.df <- as_tibble(model.data$spac) |>
                mutate(year=years) |>
                pivot_longer(matches("[0-9]"), names_to="age", values_to="val") |>
                filter(age > 2) |>
                mutate(
                    fill.color=colors,
                    type="spawn"
                ) |>
                print(n=30)
raw.df <- rbind(seac.raw.df, spac.raw.df)


#-------------------------------------------------------------------------------

#### read in model fit to age comps, calculate posterior intervals ####

# save locations for model age comp fits 

# CLR: added / to file paths
spawn.age.comp.fname <- here::here(dir_mcmc_out, "SpAC.csv")
seine.age.comp.fname <- here::here(dir_mcmc_out, "SeAC.csv")

# save fits and calculate posterior predictive intervals
# see ?pwsHerringBasa::generate.post.pred() for more info

pp.spawn.age.comp <- generate.post.pred(spawn.age.comp.fname, ess = model.data$ess.spawn, 
                                        years = years, nburn = nburn, dist = "multinom")
pp.seine.age.comp <- generate.post.pred(seine.age.comp.fname, ess = model.data$ess.seine, 
                                        years = years, nburn = nburn, dist = "multinom")

# save and format median and 95% credible intervals

seine.age.comp.quants <- apply(pp.seine.age.comp, 2, quantile, probs=c(0.025,0.5,0.975), na.rm=T)
spawn.age.comp.quants <- apply(pp.spawn.age.comp, 2, quantile, probs=c(0.025,0.5,0.975), na.rm=T)

cnames <- rep(NA, ncol(pp.spawn.age.comp))
i=1
for(y in years){
    for(a in 0:9){
        cnames[i] <- paste0(y, "_", a)
        i = i+1
    }
}
colnames(seine.age.comp.quants) <- cnames
colnames(spawn.age.comp.quants) <- cnames

# tidy-ify median age comps and credible intervals for plotting

seac.df <- as_tibble(seine.age.comp.quants) |>
            pivot_longer(everything(), names_to="year_age", values_to="val") |>
            separate(year_age, c("year", "age"), sep="_") |>
            mutate(percentile=rep(c("2.5%", "50%", "97.5%"), each=(10*nyr))) |>
            pivot_wider(names_from=percentile, values_from=val) |>
            filter(age > 2) |>
            mutate(
                type="seine"
            ) |>
            print(n=30)

spac.df <- as_tibble(spawn.age.comp.quants) |>
            pivot_longer(everything(), names_to="year_age", values_to="val") |>
            separate(year_age, c("year", "age"), sep="_") |>
            mutate(percentile=rep(c("2.5%", "50%", "97.5%"), each=(10*nyr))) |>
            pivot_wider(names_from=percentile, values_from=val) |>
            filter(age > 2) |>
            mutate(
                type="spawn"
            ) |>
            print(n=30)

age.comp.df <- rbind(seac.df, spac.df)
age.comp.df$age <- factor(age.comp.df$age, levels=sort(unique(age.comp.df$age)), ordered=TRUE)
age.comp.df$type <- factor(age.comp.df$type)

age.class.df <- data.frame(year=raw.df$year, age=rep(c("3", "4", "5", "6", "7", "8", "9"), nrow(raw.df)), color=rep(c("black", "white", "black", "white", "black", "white", "black"), nrow(raw.df)), type=rep(c("seine", "spawn"), each=nrow(raw.df)/2))
year.df <- data.frame(year=years, year_str=years)


#-------------------------------------------------------------------------------

#### read in format model-estimated numbers-at-age ####

# read-in Numbers-at-age and output into nice tables (.csv) 

n.y.a<-read.csv(here::here(dir_mcmc_out, "Num_at_age.csv"), header = FALSE, dec=".") 
n.y.a<-n.y.a[-c(1:nburn), ] # Just 'cause it's easier to work with this scale

## seine.age.comp.interval is a 2 X (# of ages in age comp * years) matrix; e.g. 10 ages (ages 0-9+) * 35 years (from 1980-2014)=350
n.y.a.interval <- matrix(rep(0, 2*length(n.y.a[1, ])), nrow=2, ncol=length(n.y.a[1, ]))
n.y.a.interval[1:2,] <- apply(n.y.a, 2, quantile, probs=c(0.025, 0.975), na.rm=T) # Fills first row with 0.025 percentile and second with 0.975 percentile
n.y.a.median <- apply(n.y.a, 2, median)
## Now cut each up into years

n.y.a.median.mat <- matrix(n.y.a.median, nrow=nyr, ncol=ncol, byrow=T) # nyr x nages matrix filled with posterior predictive median
n.y.a.lower.mat <- matrix(n.y.a.interval[1, ], nrow=nyr, ncol=ncol, byrow=T)
n.y.a.upper.mat <- matrix(n.y.a.interval[2, ], nrow=nyr, ncol=ncol, byrow=T)

# calculate shannon weiner evenness index (J) for spawner age comps

spac.50 <- spac.df |> 
            select(year, age, `50%`) |>
            filter(!is.na(`50%`)) |>
            pivot_wider(names_from=age, values_from=`50%`) |>
            select(c(`3`, `4`, `5`, `6`, `7`, `8`, `9`)) |>
            as.matrix()
evenness.50 <- apply(spac.50, 1, sw_evenness)

# format numbers-at-age for ggridges plot

naa_ridgedat <- n.y.a.median.mat |>
    as.data.table() |>
    mutate(Year = years) |>
    melt(id.vars = "Year", variable.name = "Age", value.name = "Numbers") |>
    mutate(Age = as.numeric(Age), Numbers = round(Numbers)) |> 
    uncount(Numbers)


#-------------------------------------------------------------------------------

#### make plots ####

# make age composition plot

raw.df <- raw.df |> 
    left_join(data.frame(year = 1982:(curr.year+nyr.sim-1), J = round(evenness.50, 2)))

year.df$J <- raw.df |> 
    as.data.frame() |>
    slice(match(year.df$year, raw.df$year)) |>
    select(J) |>
    unlist()

age.struct.plot <- ggplot(raw.df)+
    geom_col(aes(x=type, y=val/100, color=age, fill=fill.color), position=position_dodge(0.9), linewidth=0)+
    # scale_fill_manual(values=c(color.options, "grey")) + 
    scale_fill_manual(values=color.options) + 
    geom_pointinterval(data=age.comp.df, aes(x=type, y=`50%`/100, ymin=`2.5%`/100, ymax=`97.5%`/100, color=age), position=position_dodge(0.9)) +
    geom_text(data=year.df, aes(x=0.65, y=1, label=year), size=4)+
    geom_text(data=year.df, aes(x=2.25, y=1, label=ifelse(is.na(J), NA, paste0("J=", J))), size=4)+
    geom_text(data=age.class.df, aes(x=type, y=-0.25, label=age, group=age), position=position_dodge(0.9), size=3)+
    scale_color_manual(values=rep("black", 7))+
    scale_y_continuous("Proportion", breaks=c(0.0, 0.5, 1.0))+
    scale_x_discrete(labels=c("Seine Catch", "Spawner Survey"))+
    facet_wrap(~year, nrow=12, dir="v")+
    coord_cartesian(clip="off", ylim=c(0, 1.15))+
    labs(x="Age class", y="Proportion", title="Proportional age structure")+
    theme_bw() +
    theme(
        legend.position = "none",
        # panel.background = element_rect(fill = "grey95"),
        panel.border = element_rect(linewidth = 0),
        panel.grid = element_blank(),
        panel.spacing.y = unit(0.0, "line"),
        panel.spacing.x = unit(0.25, "line"),
        axis.text.x = element_text(margin=margin(t=12)),
        axis.ticks.x = element_blank(),
        axis.title.x = element_text(margin=margin(t=10)),
        strip.background = element_blank(),
        strip.text.x = element_blank()
    )

# make numbers-at-age plot

# naa_plot <- ggplot(naa_ridgedat, aes(x = Age, y = Year, fill = Year, group = Year)) +
#     geom_density_ridges2() +
#     scale_fill_viridis_c(option = viridis_palette, trans = "reverse") + 
#     scale_y_reverse() +
#     xlim(c(0, 12)) +
#     labs(title = "Estimated Numbers-at-age", subtitle = "Prince William Sound Herring, 1980-Present") +
#     theme_bw()


#-------------------------------------------------------------------------------

#### save outputs ####

# save age comps plot

ggsave(here::here(dir_figures, "age_compositions.pdf"), plot = age.struct.plot, height = 11, width=8.5)

# save numbers-at-age plot

# ggsave(here::here(dir_figures, "numbers_at_age.pdf"), plot = naa_plot, height = 11, width=8.5)

# save csv table of median age comps with 95% posterior predictive intervals

age.comp.median.wide <- age.comp.df |>
                        select(year, age, `50%`, type) |>
                        mutate(`50%` = round(`50%`/100, 3)) |>
                        pivot_wider(
                            everything(),
                            names_from=age,
                            values_from=`50%`
                        ) |>
                        print(n=80)

age.comp.lower.wide <- age.comp.df |>
                        select(year, age, `2.5%`, type) |>
                        mutate(`2.5%` = round(`2.5%`/100, 3)) |>
                        pivot_wider(
                            everything(),
                            names_from=age,
                            values_from=`2.5%`
                        ) |>
                        print(n=80)

age.comp.upper.wide <- age.comp.df |>
                        select(year, age, `97.5%`, type) |>
                        mutate(`97.5%` = round(`97.5%`/100, 3)) |>
                        pivot_wider(
                            everything(),
                            names_from=age,
                            values_from=`97.5%`
                        ) |>
                        print(n=80)

sink(here::here(dir_outputs, "predicted-age-comps.csv"))

cat("# Seine Fishery Age Composition (proportion of each age): Median Posterior Prediction\n")
write.table(
    age.comp.median.wide |> filter(type == "seine") |> select(-c(type)), 
    row.names=FALSE, col.names=c("# Year", 3:8, "9+"), sep=","
)

cat("\n# Seine Fishery Age Composition (proportion of each age): Lower 95% Posterior Predictive Interval\n")
write.table(
    age.comp.lower.wide |> filter(type == "seine") |> select(-c(type)),  
    row.names=FALSE, col.names=c("# Year", 3:8, "9+"), sep=","
)

cat("\n# Seine Fishery Age Composition (proportion of each age): Upper 95% Posterior Predictive Interval\n")
write.table(
    age.comp.upper.wide |> filter(type == "seine") |> select(-c(type)),  
    row.names=FALSE, col.names=c("# Year", 3:8, "9+"), sep=","
)

cat("\n\n# Spawner Survey Age Composition (proportion of each age): Median Posterior Prediction\n")
write.table(
    age.comp.median.wide |> filter(type == "spawn") |> select(-c(type)), 
    row.names=FALSE, col.names=c("# Year", 3:8, "9+"), sep=","
)
cat("\n# Spawner Survey Age Composition (proportion of each age): Lower 95% Posterior Predictive Interval\n")
write.table(
    age.comp.lower.wide |> filter(type == "spawn") |> select(-c(type)),  
    row.names=FALSE, col.names=c("# Year", 3:8, "9+"), sep=","
)
cat("\n# Spawner Survey Age Composition (proportion of each age): Upper 95% Posterior Predictive Interval\n")
write.table(
    age.comp.upper.wide |> filter(type == "spawn") |> select(-c(type)),  
    row.names=FALSE, col.names=c("# Year", 3:8, "9+"), sep=","
)
sink()


# save csv table of numbers-at-age

# CLR: changed long.years to years and round_df to round in the following write.table() calls

sink(here::here(dir_outputs, "numbers-at-age.csv"))
cat("# Numbers-at-age (millions): Median estimate\n")
write.table(
    cbind(years, round(n.y.a.median.mat, 3), round(apply(n.y.a.median.mat, 1, sw_evenness), 2)), 
    row.names=FALSE, col.names=c("# Year", 0:8, "9+", "J"), sep=","
)

cat("\n# Numbers-at-age (millions): Lower 95% Credibility Interval\n")
write.table(
    cbind(years, round(n.y.a.lower.mat, 3), round(apply(n.y.a.lower.mat, 1, sw_evenness), 2)),  
    row.names=FALSE, col.names=c("# Year", 0:8, "9+", "J"), sep=","
)

cat("\n# Numbers-at-age (millions): Upper 95% Credibility Interval\n")
write.table(
    cbind(years, round(n.y.a.upper.mat, 3), round(apply(n.y.a.upper.mat, 1, sw_evenness), 2)),  
    row.names=FALSE, col.names=c("# Year", 0:8, "9+", "J"), sep=","
)

sink()

>>>>>>> fa4cd14 (incorporate trevors report edits)
