################################################################################

# Plot retrospectives

# Creates a line plot for all retrospective model fits. Retrospective analysis
# is performed in run_retrospectives.r.

# authors: Joshua Zahner, CL Roberts

# inputs: BASA outputs for each retrospective peel. Each peeled run is saved to a 
#          subdirectory of 'retrospectives/'.

# outputs: line plot of retrospective model fits

################################################################################

#### front matter ####

# load packages

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggdist)
library(icesAdvice)
library(pwsHerringBasa)

# directories

dir_model <- here::here("model")
dir_mcmc_out <- here::here(dir_model, "mcmc_out")
dir_retro <- here::here("retrospectives")
dir_figures <- here::here("figures")

if (!dir.exists(dir_figures)) {
    dir.create(dir_figures)
}

if (!dir.exists(dir_retro)) {
    stop("Retrospective analysis has not yet been performed! Execute run_retrospectives.r ...")
}

# save local variables

nyr <- read.data.files(dir_model)$PWS_ASA.dat$nyr
start.year <- 1980
curr.year <- start.year+nyr

n.peels <- list.dirs(dir_retro, recursive = FALSE, full.names = FALSE) |>
    grep(pattern = "basa-") |>
    length()

#### analysis ####

# calc mohn's rho

annual.biomass.estimate <- compute.retrospective.matrix(n.peels, dir_mcmc_out)

mohns.rho <- icesAdvice::mohn(annual.biomass.estimate, details=TRUE)
rel.bias <- c(mohns.rho$compare[,"relbias"], NA)
rho <- mohns.rho$rho

saveRDS(rho, file = here::here("report/mohns_rho.rds"))

# Compute the peel estimates for each run
# final.year.est <- annual.biomass.estimate[38:nrow(annual.biomass.estimate), 1]
peel.estimates <- apply(as.matrix(annual.biomass.estimate), 2, function(x) x[max(which(!is.na(x)))])
peel.estimates <- rev(peel.estimates)

retro <- as_tibble(annual.biomass.estimate) |> 
    pivot_longer(everything(), names_to="lag", values_to="biomass") |>
    arrange(lag) |>
    mutate(year=rep(1980:curr.year, length.out=n())) |>
    mutate(
        lag = factor(lag, levels=c("base", paste0("-", 1:n.peels)), labels = c(curr.year:(curr.year-n.peels)))
    )

final.year <- retro |> 
    mutate(year.lag=curr.year-as.numeric(lag)) |> 
    filter(if_all(c(year.lag, year), ~ year == .x))

biomass.df <- compute.biomass.traj(dir_mcmc_out, nyr) |> 
                mutate(year = as.numeric(year)) |> 
                filter(.width >= 0.5)

ggplot(retro) +
    geom_lineribbon(data=biomass.df, aes(x=year, y=biomass, ymin=.lower, ymax=.upper), size=0, alpha=0.35)+
    geom_line(aes(x=year, y=biomass, color=lag), size=1)+
    geom_point(data=final.year, aes(x=year, y=biomass, color=lag), size=3)+
    geom_hline(yintercept=c(20000, 40000), linetype="dashed")+
    geom_text(data=data.frame(), aes(x=2016, y=180000, label=paste("Mohn's rho:", round(rho, 3))))+
    scale_fill_grey(start=0.8, end=0.4)+
    scale_y_continuous(limits=c(0, 200000), breaks=c(0, 20000, 40000, 50000, 100000, 150000, 200000), labels=scales::comma)+
    coord_cartesian(expand=0)+
    labs(x="Year", y="Spawning biomass (mt)", color="Data lag", title=paste0(n.peels, "-year retrospective pattern"))+
    theme_classic()+
    theme(
        axis.title = element_text(size=16),
        axis.text = element_text(size=12),
        plot.title = element_text(size=18)
    )

ggsave(here::here(dir_figures, "retrospective.pdf"))
