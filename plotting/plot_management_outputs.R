<<<<<<< HEAD
################################################################################

# plot management outputs

# plots BASA-estimated and derived time series parameters of interest to management

# authors: John Trochta, Joshua Zahner, CL Roberts

# inputs: BASA model inputs (model/PWS_ASA.dat) and outputs (from mcmc_out/)

# outputs: 
#   - plots: 
#       - recruitment (figures/recruitment_trajectory.pdf)
#       - spawning biomass (figures/biomass_trajectory.pdf)
#       - exploitation rate (figures/exploitation_rate.pdf)
#       - pre-fishery biomass posterior (figures/pfrb_posterior.pdf)
#       - a 4-panel composite plot of all outputs (figures/management_outputs.pdf)
#   - tables: 
#       - a single .csv file (data_outputs/outputs_for_management.csv) giving 
#         estimated recruitment, PFRB, catch, exploitation rate, and probability 
#         below harvest threshold with 95% credible intervals where relevant

################################################################################


#### front matter ####

# choose TMB or ADMB
software <- "ADMB"

# load packages

library(ggplot2)
#library(ggpubr)
library(pwsHerringBasa)
library(dplyr)

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

# read model input data

raw.data <- read.data.files(dir_model)

# save global variables
nyr <- raw.data$PWS_ASA.dat$nyr
start.year <- 1980
curr.year <- start.year+nyr
nyr.sim <- 0
years <- seq(start.year, curr.year+nyr.sim-1)


#-------------------------------------------------------------------------------

#### calculate and plot outputs ####

# the following functions read in, compute posteriors, and plot each management 
# output of interest see ?pwsHerringBasa::compute.* and ?pwsHerringBasa::plot_* 
# for more details

recruit.df <- compute.recruitment(dir_mcmc_out, nyr, years)
recruit.plot <- plot_recruitment_posterior(recruit.df, years, legend=FALSE)

biomass.df <- compute.biomass.traj(dir_mcmc_out, nyr)
biomass.plot <- plot_biomass_trajectory(biomass.df, c(years, curr.year), legend=FALSE)

exploit.df <- compute.exploit.rate(dir_mcmc_out, nyr)
exploit.rate.plot <- plot_exploit_rate(exploit.df$exploit.rate.df,
                                       exploit.df$exploit.zeros,
                                       years)

pfrb.posterior <- compute.pfrb.posterior(dir_mcmc_out, nyr+1)
pfrb.posterior.plot <- plot_pfrb_posterior(pfrb.posterior$biomass.df, 
                                           pfrb.posterior$biomass.quants, 
                                           pfrb.posterior$prob.below.threshold,
                                           curr.year,
                                           font.size=5)

#-------------------------------------------------------------------------------

#### save output files ####

# save individual plots

ggsave(here::here(dir_figures, "recruitment_trajectory.pdf"), plot = recruit.plot, height=8.5, width=11) 
ggsave(here::here(dir_figures, "biomass_trajectory.pdf"), plot = biomass.plot, height=8.5, width=11) 
ggsave(here::here(dir_figures, "exploitation_rate.pdf"), plot = exploit.rate.plot, height=8.5, width=11) 
ggsave(here::here(dir_figures, "pfrb_posterior.pdf"), plot = pfrb.posterior.plot, height=8.5, width=11) 


# save 4-panel composite plot

#management_outputs <- ggarrange(
#    recruit.plot, biomass.plot, exploit.rate.plot, pfrb.posterior.plot,
#    nrow=2,
#    ncol=2
# )
# ggsave(here::here(dir_figures, "management_outputs.pdf"), management_outputs, 
#        height=8.5, width=11)


# save outputs-for-management.csv table

biomass.matrix <- biomass.df |> 
    filter(.width==0.95, year != "forecast") |> 
    select(biomass, .lower, .upper) |> 
    as.matrix() |> 
    round(2)
recruit.matrix <- recruit.df |> 
    filter(.width==0.95) |> 
    select(recruits, .lower, .upper) |> 
    as.matrix() |> 
    round(2)
exploit.matrix <- exploit.df$exploit.rate.df |> 
    filter(.width==0.95) |> 
    select(exploit, .lower, .upper) |> 
    as.matrix() |> 
    round(2)
prob.matrix <- biomass.df |> 
    filter(.width==0.95, year != "forecast") |> 
    select(prob) |> 
    as.matrix() |> 
    round(2)

total.catch.biomass <- compute.catch.biomass(data = raw.data$PWS_ASA.dat) |> 
    round(2) 

final.table <- data.frame(years, recruit.matrix, biomass.matrix/1000, 
                          total.catch.biomass, exploit.matrix, prob.matrix) 
names(final.table) <- c("Years",
                        "Median Age 3 (in millions)",
                        "Lower 95th Age 3 (in millions)",
                        "Upper 95th Age 3 (in millions)",
                        "Median Pre-fishery biomass (in 1000s metric tons)",
                        "Lower 95th Biomass (in 1000s metric tons)",
                        "Upper 95th Biomass (in 1000s metric tons)",
                        "Catch (metric tons)",
                        "Median Exploitation rate",
                        "Lower 95th ER",
                        "Upper 95th ER",
                        "Probability B<20K")
write.csv(final.table, here::here(dir_outputs, "outputs-for-management.csv"), row.names=FALSE)
=======
################################################################################

# plot management outputs

# plots BASA-estimated and derived time series parameters of interest to management

# authors: John Trochta, Joshua Zahner, CL Roberts

# inputs: BASA model inputs (model/PWS_ASA.dat) and outputs (from mcmc_out/)

# outputs: 
#   - plots: 
#       - recruitment (figures/recruitment_trajectory.pdf)
#       - spawning biomass (figures/biomass_trajectory.pdf)
#       - exploitation rate (figures/exploitation_rate.pdf)
#       - pre-fishery biomass posterior (figures/pfrb_posterior.pdf)
#       - a 4-panel composite plot of all outputs (figures/management_outputs.pdf)
#   - tables: 
#       - a single .csv file (data_outputs/outputs_for_management.csv) giving 
#         estimated recruitment, PFRB, catch, exploitation rate, and probability 
#         below harvest threshold with 95% credible intervals where relevant

################################################################################


#### front matter ####

# choose TMB or ADMB
software <- "ADMB"

# load packages

library(ggplot2)
library(ggpubr)
library(pwsHerringBasa)
library(dplyr)

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

# read model input data

raw.data <- read.data.files(dir_model)

# save global variables
nyr <- raw.data$PWS_ASA.dat$nyr
start.year <- 1980
curr.year <- start.year+nyr
nyr.sim <- 0
years <- seq(start.year, curr.year+nyr.sim-1)


#-------------------------------------------------------------------------------

#### calculate and plot outputs ####

# the following functions read in, compute posteriors, and plot each management 
# output of interest see ?pwsHerringBasa::compute.* and ?pwsHerringBasa::plot_* 
# for more details

recruit.df <- compute.recruitment(dir_mcmc_out, nyr, years)
recruit.plot <- plot_recruitment_posterior(recruit.df, years, legend=FALSE) +
    labs(main = "Age-3 recruitment") +
    theme(axis.title = element_text(size=12),
          title = element_text(size=14)) 

biomass.df <- compute.biomass.traj(dir_mcmc_out, nyr)
biomass.plot <- plot_biomass_trajectory(biomass.df, c(years, curr.year), legend=FALSE) +
    labs(main = "Biomass trajectory") +
    theme(axis.title = element_text(size=12),
          title = element_text(size=14))

exploit.df <- compute.exploit.rate(dir_mcmc_out, nyr)
exploit.rate.plot <- plot_exploit_rate(exploit.df$exploit.rate.df,
                                       exploit.df$exploit.zeros,
                                       years) +
    labs(main = "Exploitation rate") +
    theme(axis.title = element_text(size=12),
          title = element_text(size=14))

pfrb.posterior <- compute.pfrb.posterior(dir_mcmc_out, nyr+1)
pfrb.posterior.plot <- plot_pfrb_posterior(pfrb.posterior$biomass.df, 
                                           pfrb.posterior$biomass.quants, 
                                           pfrb.posterior$prob.below.threshold,
                                           curr.year,
                                           font.size=5) +
    labs(main=paste0(curr.year+1, "pre-fishery biomass posterior probability density"),
         x=paste0(curr.year+1, "pre-fishery biomass (mt)")) +
    theme(axis.title = element_text(size=12),
          title = element_text(size=14))

#-------------------------------------------------------------------------------

#### save output files ####

# save individual plots

ggsave(here::here(dir_figures, "recruitment_trajectory.pdf"), plot = recruit.plot, height=8.5, width=11) 
ggsave(here::here(dir_figures, "biomass_trajectory.pdf"), plot = biomass.plot, height=8.5, width=11) 
ggsave(here::here(dir_figures, "exploitation_rate.pdf"), plot = exploit.rate.plot, height=8.5, width=11) 
ggsave(here::here(dir_figures, "pfrb_posterior.pdf"), plot = pfrb.posterior.plot, height=8.5, width=11) 


# save 4-panel composite plot

management_outputs <- ggarrange(
    recruit.plot, biomass.plot, exploit.rate.plot, pfrb.posterior.plot,
    nrow=2,
    ncol=2
)
ggsave(here::here(dir_figures, "management_outputs.pdf"), management_outputs, 
       height=8.5, width=11)


# save outputs-for-management.csv table

biomass.matrix <- biomass.df |> 
    filter(.width==0.95, year != "forecast") |> 
    select(biomass, .lower, .upper) |> 
    as.matrix() |> 
    round(2)
recruit.matrix <- recruit.df |> 
    filter(.width==0.95) |> 
    select(recruits, .lower, .upper) |> 
    as.matrix() |> 
    round(2)
exploit.matrix <- exploit.df$exploit.rate.df |> 
    filter(.width==0.95) |> 
    select(exploit, .lower, .upper) |> 
    as.matrix() |> 
    round(2)
prob.matrix <- biomass.df |> 
    filter(.width==0.95, year != "forecast") |> 
    select(prob) |> 
    as.matrix() |> 
    round(2)

total.catch.biomass <- compute.catch.biomass(data = raw.data$PWS_ASA.dat) |> 
    round(2) 

final.table <- data.frame(years, recruit.matrix, biomass.matrix/1000, 
                          total.catch.biomass, exploit.matrix, prob.matrix) 
names(final.table) <- c("Years",
                        "Median Age 3 (in millions)",
                        "Lower 95th Age 3 (in millions)",
                        "Upper 95th Age 3 (in millions)",
                        "Median Pre-fishery biomass (in 1000s metric tons)",
                        "Lower 95th Biomass (in 1000s metric tons)",
                        "Upper 95th Biomass (in 1000s metric tons)",
                        "Catch (metric tons)",
                        "Median Exploitation rate",
                        "Lower 95th ER",
                        "Upper 95th ER",
                        "Probability B<20K")
write.csv(final.table, here::here(dir_outputs, "outputs-for-management.csv"), row.names=FALSE)
>>>>>>> fa4cd14 (incorporate trevors report edits)
