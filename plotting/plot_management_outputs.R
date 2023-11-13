library(ggplot2)
library(ggpubr)
library(ggdist)

source(file=paste0(here::here("plotting/"), "compute_plot_products.R"))
source(file=paste0(here::here("plotting/"), "plot_utils.R"))
source(paste0(here::here("functions/"), "fun_read_dat.R"))

start.year <- 1980
curr.year <- 2024
nyr.sim <- 0
years <- seq(start.year, curr.year+nyr.sim)
nyr <- length(years)

model.dir <- here::here("model/")

biomass.df <- compute.biomass.traj(model.dir, nyr, years)
biomass.plot <- plot.biomass.trajectory(biomass.df, years, legend=FALSE, new.theme=present.theme)
# biomass.plot+present.theme
ggsave(paste0(here::here("figures/"), "biomass_trajectory.pdf"), height=8.5, width=11)

exploit.df <- compute.exploit.rate(model.dir, nyr-1, years[1:(nyr-1)])
exploit.rate.plot <- plot.exploit.rate(exploit.df$exploit.rate.df,
                                       exploit.df$exploit.zeros,
                                       years)
ggsave(paste0(here::here("figures/"), "exploitation_rate.pdf"), height=8.5, width=11)

pfrb.posterior <- compute.pfrb.posterior(model.dir, nyr, years)
pfrb.posterior.plot <- plot.pfrb.posterior(pfrb.posterior$biomass.df, 
                                           pfrb.posterior$biomass.quants, 
                                           pfrb.posterior$prob.below.threshold,
                                           curr.year,
                                           font.size=5)
ggsave(paste0(here::here("figures/"), "pfrb_posterior.pdf"), height=8.5, width=11)

recruit.df <- compute.recruitment(model.dir, nyr-1, years[1:(nyr-1)])
recruit.plot <- plot.recruitment.posterior(recruit.df, years, legend=FALSE)
ggsave(paste0(here::here("figures/"), "recrtuitment_trajectory.pdf"), height=8.5, width=11)

ggarrange(
    recruit.plot, biomass.plot, exploit.rate.plot, pfrb.posterior.plot,
    nrow=2,
    ncol=2
)

ggsave(paste0(here::here("figures/"), "management_outputs.pdf"), height=8.5, width=11)

biomass.matrix <- biomass.df %>% filter(.width==0.95) %>% select(biomass, .lower, .upper) %>% head(nyr-1) %>% as.matrix %>% round(., 2)
recruit.matrix <- recruit.df %>% filter(.width==0.95) %>% select(recruits, .lower, .upper) %>% as.matrix %>% round(., 2)
exploit.matrix <- exploit.df$exploit.rate.df %>% filter(.width==0.95) %>% select(exploit, .lower, .upper) %>% as.matrix %>% round(., 2)
prob.matrix    <- biomass.df %>% filter(.width==0.95) %>% select(prob) %>% head(nyr-1) %>% as.matrix %>% round(., 2)


raw.data <- read.data.files(model.dir)
total.catch.biomass <- compute.catch.biomass(raw.data$PWS_ASA.dat, nyr-1) %>% round(., 2)

final.table <- data.frame(years[1:(nyr-1)], recruit.matrix, biomass.matrix/1000, total.catch.biomass, exploit.matrix, prob.matrix) 
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
write.csv(final.table, here::here("data_outputs/outputs-for-management.csv"), row.names=FALSE)
