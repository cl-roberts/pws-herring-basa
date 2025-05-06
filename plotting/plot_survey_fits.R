################################################################################

# plot survey fits

# plots BASA model fits to survey data for all five indices collected on PWS 
# herring 1980-present

# authors: John Trochta, Joshua Zahner, CL Roberts

# inputs: BASA model inputs (model/PWS_ASA.dat) and outputs (from mcmc_out/)

# outputs: 
#   - plots: five individuals plots to all survey fits, and a combined plot
#   - tables: a single csv file ("data_outputs/predicted-survey-values.csv") 
#             containing all survey fits

################################################################################


#### front matter ####

# choose TMB or ADMB
software <- "ADMB"

# attach packages

library(pwsHerringBasa)
library(dplyr)
library(ggpubr)

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

# load data

raw.data <- read.data.files(dir_model)$PWS_ASA.dat

# local vars 

nyr <- raw.data$nyr
start.year <- 1980
curr.year <- start.year+nyr
nyr.sim <- 0
years <- seq(start.year, curr.year+nyr.sim-1)

#-------------------------------------------------------------------------------

#### read-in model outputs, calculate posteriors, and plot ####

# model estimated variances

variances <- read.table(here::here(dir_mcmc_out, "VARSReport.csv"), 
                        sep=",", header=FALSE)[-1,]
names(variances) <- c("mdm", "egg", "adfg.hydro", "pwssc.hydro")

add.vars <- apply(variances, 2, median)
cvs <- list(
    mdm         = as.vector(rep(add.vars[1], nyr)),
    egg         = as.vector(sqrt(raw.data$egg_se[1:nyr]^2 + add.vars[2]^2)),
    adfg.hydro  = as.vector(rep(add.vars[3], nyr)),
    pwssc.hydro = as.vector(sqrt((raw.data$pwssc_hydro_se[1:nyr]^2) + (add.vars[4]^2)))
)

# juv.overdisp <- read.table(here::here(dir_mcmc_out, "iterations.csv"),  
#                            header = TRUE, sep = ",")[, "juvenile_overdispersion"]
# juv.schools.mask <- c(raw.data$juvenile_survey == -9)[1:nyr]

# save file paths containing model outputs

mdm.fname <- here::here(dir_mcmc_out, "MDM.csv")
egg.fname <- here::here(dir_mcmc_out, "EGG.csv")
adfg.hydro.fname <- here::here(dir_mcmc_out, "HYD_ADFG.csv")
pwssc.hydro.fname <- here::here(dir_mcmc_out, "HYD_PWSSC.csv")
# juv.schools.fname <- here::here(dir_mcmc_out, "juv_schools.csv")

# calculate posterior predictive intervals for each survey fit
# see ?pwsHerringBasa::generate.post.pred() for more info

mdm.pp <- generate.post.pred(fname = mdm.fname, var = variances$mdm, years = years)
egg.pp <- generate.post.pred(fname = egg.fname, var = variances$egg, years = years)
adfg.hydro.pp <- generate.post.pred(fname = adfg.hydro.fname, 
                                    var = variances$adfg.hydro, years = years)
pwssc.hydro.pp <- generate.post.pred(fname = pwssc.hydro.fname, 
                                     var = variances$pwssc.hydro, years = years)
# juv.schools.pp <- generate.post.pred(fname = juv.schools.fname, var = juv.overdisp, 
#                                      years = years, dist="negbin", mask=juv.schools.mask)


# save and format raw survey data for each annual survey

mdm.data            <- data.frame(year=as.character(years), data=raw.data$mdm[1:length(years)], lower=raw.data$mdm[1:length(years)]/calc.buck.cv(cvs$mdm), upper=raw.data$mdm[1:length(years)]*calc.buck.cv(cvs$mdm))
egg.data            <- data.frame(year=as.character(years), data=raw.data$egg[1:length(years)], lower=raw.data$egg[1:length(years)]/calc.buck.cv(cvs$egg), upper=raw.data$egg[1:length(years)]*calc.buck.cv(cvs$egg))
adfg.hydro.data     <- data.frame(year=as.character(years), data=raw.data$adfg_hydro[1:length(years)], lower=raw.data$adfg_hydro[1:length(years)]/calc.buck.cv(cvs$adfg.hydro), upper=raw.data$adfg_hydro[1:length(years)]*calc.buck.cv(cvs$adfg.hydro))
pwssc.hydro.data    <- data.frame(year=as.character(years), data=raw.data$pwssc_hydro[1:length(years)], lower=raw.data$pwssc_hydro[1:length(years)]/calc.buck.cv(cvs$pwssc.hydro), upper=raw.data$pwssc_hydro[1:length(years)]*calc.buck.cv(cvs$pwssc.hydro))
# aer.juvenile.data   <- data.frame(year=as.character(years), data=raw.data$juvenile_survey[1:length(years)])

data <- list(
    mdm = mdm.data,
    egg = egg.data,
    adfg.hydro = adfg.hydro.data,
    pwssc.hydro = pwssc.hydro.data
    # juvenile = aer.juvenile.data
)
for(i in 1:length(data)){
    missing <- data[[i]]$data == -9
    data[[i]]$data[missing] <- NA
    data[[i]]$lower[missing] <- NA  # replaces missing data (-9's) with NAs
    data[[i]]$upper[missing] <- NA
}

# plot survey fits 
# see pwsHerringBasa::plot_survey_fits()

mdm.fit.plot <- plot_survey_fits(mdm.pp, data$mdm, y.max=500, title="Mile days of milt")
egg.fit.plot <- plot_survey_fits(egg.pp, data$egg, y.max=21, title="Egg deposition (trillions)")
adfg.hydro.fit.plot <- plot_survey_fits(adfg.hydro.pp, data$adfg.hydro, y.max=175000, title="ADF&G hydroacoustic biomass (1000s tons)", scale=1000)
pwssc.hydro.fit.plot <- plot_survey_fits(pwssc.hydro.pp, data$pwssc.hydro, y.max=205000, title="PWSSC hydroacoustic biomass (1000s tons)", scale=1000)
# aer.juvenile.fit.plot <- plot_survey_fits(juv.schools.pp, data$juvenile, y.max=40000, title="Age-1 schools", cvs=FALSE)

#-------------------------------------------------------------------------------


#### save outputs ####

# save individual plots

ggsave(here::here(dir_figures, "mdm_fit.pdf"), plot = mdm.fit.plot, width=12, height=8, dpi=300)
ggsave(here::here(dir_figures, "egg_fit.pdf"), plot = egg.fit.plot, width=12, height=8, dpi=300)
ggsave(here::here(dir_figures, "adfg_hydro_fit.pdf"), plot = adfg.hydro.fit.plot, width=12, height=8, dpi=300)
ggsave(here::here(dir_figures, "pwssc_hydro_fit.pdf"), plot = pwssc.hydro.fit.plot, width=12, height=8, dpi=300)
# ggsave(here::here(dir_figures, "aer_juvenile_fit.pdf"), plot = aer.juvenile.fit.plot, width=12, height=8, dpi=300)

# save 5-panel plot

survey_fits <- ggarrange(
    mdm.fit.plot, egg.fit.plot, adfg.hydro.fit.plot, pwssc.hydro.fit.plot, 
    nrow=3,
    ncol=2,
    common.legend = TRUE, legend="right"
)

# CLR: added / to file path
ggsave(here::here(dir_figures, "survey_fits.pdf"), survey_fits, width=12, height=8, dpi=300)


# save .csv table with clean survey fits

sink(here::here(dir_outputs, "predicted-survey-values.csv"))
cat('Milt (mile-days) posterior predictions\n')
write.table(
    mdm.pp %>% filter(.width==0.95) %>% select(-c(.point, .interval, .width)) %>% rename(median=data) %>% mutate(median=round(median, 3), .lower=round(.lower, 3), .upper=round(.upper, 3)),
    row.names=FALSE,col.names=c('Year','Median Posterior Prediction','Lower 95th','Upper 95th'),sep=","
)

cat('\nEggs deposited (trillions) posterior predictions\n')
write.table(
    egg.pp %>% filter(.width==0.95) %>% select(-c(.point, .interval, .width)) %>% rename(median=data) %>% mutate(median=round(median, 3), .lower=round(.lower, 3), .upper=round(.upper, 3)),
    row.names=FALSE,col.names=c('Year','Median Posterior Prediction','Lower 95th','Upper 95th'),sep=","
)

cat('\nADF&G acoustic biomass (thousands of tons) posterior predictions\n')
write.table(
    adfg.hydro.pp %>% filter(.width==0.95) %>% select(-c(.point, .interval, .width)) %>% rename(median=data) %>% mutate(median=round(median, 3), .lower=round(.lower, 3), .upper=round(.upper, 3)),
    row.names=FALSE,col.names=c('Year','Median Posterior Prediction','Lower 95th','Upper 95th'),sep=","
)

cat('\nPWSSC acoustic biomass (thousands of tons) posterior predictions\n')
write.table(
    pwssc.hydro.pp %>% filter(.width==0.95) %>% select(-c(.point, .interval, .width)) %>% rename(median=data) %>% mutate(median=round(median, 3), .lower=round(.lower, 3), .upper=round(.upper, 3)),
    row.names=FALSE,col.names=c('Year','Median Posterior Prediction','Lower 95th','Upper 95th'),sep=","
)

# cat('\nAge-1 Aerial survey (# of small schools) posterior predictions\n')
# write.table(
#     juv.schools.pp %>% filter(.width==0.95) %>% select(-c(.point, .interval, .width)) %>% rename(median=data) %>% mutate(median=round(median, 3), .lower=round(.lower, 3), .upper=round(.upper, 3)),
#     row.names=FALSE,col.names=c('Year','Median Posterior Prediction','Lower 95th','Upper 95th'),sep=","
# )

sink()
