source(file=paste0(here::here("plotting/", "compute_plot_products.R")))
source(file=paste0(here::here("plotting/", "plot_utils.R")))

start.year <- 1980
curr.year <- 2022
nyr.sim <- 0
years <- seq(start.year, curr.year+nyr.sim-1)
nyr <- length(years)

model.dir <- here::here("model/")

biomass.df <- compute.biomass.traj(model.dir, nyr, years)
plot <- plot.biomass.trajectory(biomass.df, years, show.probs=FALSE, legend=FALSE)
plot

ggsave(paste0(here::here("figures/"), "biomass_trajectory.pdf"))
