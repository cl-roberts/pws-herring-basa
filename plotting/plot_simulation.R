library(tidyverse)
library(ggplot2)
library(ggdist)

simulation.fname <- file.path(here::here(), "data_outputs", "simulation_biomass.csv")

sim.data.ci <- read_csv(simulation.fname) %>% na.omit() %>%
                  `colnames<-`(1:ncol(.)) %>%
                  pivot_longer(everything(), names_to="year", values_to="biomass") %>%
                  mutate(year=as.numeric(year)) %>%
                  group_by(year) %>%
                  median_qi(biomass, .width=c(0.5, 0.95))


regimes <- rep(rep(rep(c("high", "low"), each=regime.length), length.out = max(sim.data.ci$year)), 2)

sim.data.ci$regime = regimes

plot.regimes <- data.frame(regime.change <- seq(0, nyr, regime.length))
ggplot(sim.data.ci, aes(x=year, y=biomass, ymin=.lower, ymax=.upper)) +
  geom_lineribbon(size=0.5)+
  geom_hline(aes(yintercept=40000), linetype="dashed")+
  geom_hline(aes(yintercept=20000), linetype="dashed")+
  geom_vline(data=plot.regimes, aes(xintercept=regime.change), color="grey", linetype="dashed")+
  scale_fill_brewer(palette="Blues")+
  scale_y_continuous("Spawning Biomass (mt)", labels=scales::comma, breaks=c(0, 20000, 40000, 50000, 100000, 150000, 200000))+
  labs(
    x = "Simulation Year",
    title = "Simulated Spawning Biomass",
    fill = "Uncertainty"
  ) + 
  coord_cartesian(ylim=c(0, 2e5), expand=0)+
  theme_classic()+
  theme(
    axis.title = element_text(size=16),
    axis.text = element_text(size=12),
    plot.title = element_text(size=18)
  )
