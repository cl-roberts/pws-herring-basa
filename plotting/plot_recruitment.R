source(file=paste0(here::here("plotting/", "compute_plot_products.R")))
source(file=paste0(here::here("plotting/", "plot_utils.R")))

start.year <- 1980
curr.year <- 2023
nyr.sim <- 0
years <- seq(start.year, curr.year+nyr.sim-1)
nyr <- length(years)

model.dir <- here::here("model/")

recruit.df <- compute.recruitment(model.dir, nyr, years)
plot <- plot.recruitment.posterior(recruit.df, years)
plot

ggsave(paste0(here::here("figures/"), "recruitment_trajectory.pdf"))
