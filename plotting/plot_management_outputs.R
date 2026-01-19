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

# load packages

library(ggplot2)
library(ggpubr)
library(pwsHerringBasa)
library(dplyr)

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
start.year <- 1980
curr.year <- start.year+nyr
nyr.sim <- 0
years <- seq(start.year, curr.year+nyr.sim-1)


#-------------------------------------------------------------------------------

#### calculate and plot outputs ####

# the following functions read in, compute posteriors, and plot each management 
# output of interest see ?pwsHerringBasa::compute.* and ?pwsHerringBasa::plot_* 
# for more details

# recruitment ----

recruit.df <- compute.recruitment(dir_mcmc_out, nyr, years) 

recruit_df <- recruit.df |>
    filter(.width == 0.5) |>
    select(year, recruits, .lower, .upper) |>
    mutate(year = as.integer(year)) |>
    rename(lower_50 = .lower, upper_50 = .upper)

recruit_df$lower_95 <- filter(recruit.df, .width == 0.95) |>
    select(.lower) |>
    unlist()

recruit_df$upper_95 <- filter(recruit.df, .width == 0.95) |>
    select(.upper) |>
    unlist()

recruit.plot <- plot_recruitment_posterior(recruit_df, years, legend = FALSE)

# biomass ----

biomass.df <- compute.biomass.traj(dir_mcmc_out, nyr)

biomass_df <- biomass.df |>
    filter(.width == 0.5) |>
    select(year, prob, biomass, .lower, .upper) |>
    mutate(year = c(years, curr.year)) |>
    rename(lower_50 = .lower, upper_50 = .upper)

biomass_df$lower_95 <- filter(biomass.df, .width == 0.95) |>
    select(.lower) |>
    unlist()

biomass_df$upper_95 <- filter(biomass.df, .width == 0.95) |>
    select(.upper) |>
    unlist()

biomass.plot <- plot_biomass_trajectory(biomass_df, c(years, curr.year), legend=FALSE)

# exploitation ----

exploit.df <- compute.exploit.rate(dir_mcmc_out, nyr)

exploit_df <- exploit.df$exploit.rate.df |>
    select(year, exploit, .lower, .upper) |>
    mutate(year = as.numeric(year)) |>
    rename(lower_95 = .lower, upper_95 = .upper)

exploit.rate.plot <- plot_exploit_rate(exploit_df, years)

# exploitation ----

pfrb.posterior <- compute.pfrb.posterior(dir_mcmc_out, nyr+1)
pfrb.posterior.plot <- plot_pfrb_posterior(pfrb.posterior$biomass.df, 
                                           pfrb.posterior$biomass.quants, 
                                           pfrb.posterior$prob.below.threshold,
                                           curr.year,
                                           font.size=14)

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