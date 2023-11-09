library(here)
library(tidyverse)
library(icesAdvice)
source(file=file.path(here::here(), "plotting/compute_plot_products.R"))

n.peels <- 5
yr <- 2024

annual.biomass.estimate <- compute.retrospective.matrix(n.peels)

mohns.rho <- icesAdvice::mohn(annual.biomass.estimate, details=TRUE)
rel.bias <- c(mohns.rho$compare[,"relbias"], NA)
rho <- mohns.rho$rho

# Compute the peel estimates for each run
#final.year.est <- annual.biomass.estimate[38:nrow(annual.biomass.estimate), 1]
peel.estimates <- apply(as.matrix(annual.biomass.estimate), 2, function(x) x[max(which(!is.na(x)))])
peel.estimates <- rev(peel.estimates)

retro <- as_tibble(annual.biomass.estimate) %>% 
    pivot_longer(everything(), names_to="lag", values_to="biomass") %>%
    arrange(lag) %>%
    mutate(year=rep(1980:yr, length.out=n())) %>%
    mutate(
        lag = factor(lag, levels=c("base", paste0("-", 1:n.peels)), labels = c(yr:(yr-n.peels)))
    )

final.year <- retro %>% mutate(year.lag=yr-as.numeric(lag)) %>% filter(if_all(c(year.lag, year), ~ year == .x))

biomass.df <- compute.biomass.traj(here::here("model"), length(1980:(yr)), 1980:(yr)) %>% mutate(year=as.numeric(year)) %>% filter(.width >= 0.5)

ggplot(retro) +
    geom_lineribbon(data=biomass.df, aes(x=year, y=biomass, ymin=.lower, ymax=.upper), size=0, alpha=0.35)+
    geom_line(aes(x=year, y=biomass, color=lag), size=1)+
    geom_point(data=final.year, aes(x=year, y=biomass, color=lag), size=3)+
    geom_hline(yintercept=c(20000, 40000), linetype="dashed")+
    geom_text(data=data.frame(), aes(x=2020, y=180000, label=paste("Mohn's Rho:", round(rho, 3))))+
    scale_fill_grey(start=0.8, end=0.4)+
    scale_y_continuous(limits=c(0, 200000), breaks=c(0, 20000, 40000, 50000, 100000, 150000, 200000), labels=scales::comma)+
    coord_cartesian(expand=0)+
    labs(x="Year", y="Spawning Biomass (mt)", color="Data Lag", title=paste0(n.peels, "-year Retrospective Pattern"))+
    theme_classic()+
    theme(
        axis.title = element_text(size=16),
        axis.text = element_text(size=12),
        plot.title = element_text(size=18)
    )

ggsave(file.path(here::here(), "figures", "retrospective.pdf"))
