#' Plot Biomass Trajectory
#'
#' Wrapper function for plotting the model-estimated biomass trajectory with
#' credible intervals.
#'
#' @param df A data frame (or \link[tibble]{tibble}) returned by
#' \link[pwsHerringBasa]{compute.biomass.traj} containing model-estimated biomass
#' time series with credible intervals
#' @param years Integer vector of years in time series
#' @param legend Logical; if TRUE display legend
#' @param show.probs Logical; if TRUE display biomass credible intervals
#'
#' @returns A \link[ggplot2]{ggplot} object showing model-estimated biomass
#' (metric tons) over the full model time series, as well as credible intervals
#'


plot_biomass_trajectory <- function(df, years, legend=TRUE, show.probs=TRUE){

  out <- ggplot(df) +
    geom_ribbon(aes(x=year, ymin=lower_50/1000, ymax=upper_50/1000), alpha = .25) +
    geom_ribbon(aes(x=year, ymin=lower_95/1000, ymax=upper_95/1000), alpha = .25) +
    geom_line(aes(x=year, y = biomass/1000), linewidth=0.25) +
    scale_x_continuous("Year", breaks=seq(min(years), max(years), by=5), expand=c(.05, .05))+
    geom_hline(yintercept=c(20, 40), linetype="dashed")+
    ggtitle("Biomass trajectory")+
    scale_y_continuous(
      "Pre-fishery biomass (1000 mt)",
      breaks=c(0, 20, 40, 50, 100, 150, 200),
      expand=c(0, 0))+
    coord_cartesian(clip="off") +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      legend.position = ifelse(legend, "right", "none")
    )

  if(show.probs){
    out <- out + geom_point(aes(x=year, y=prob*200))+
      geom_line(aes(x=year, y=prob*200))+
      scale_y_continuous(
        "Pre-fishery biomass (1000 mt)",
        breaks=c(0, 20, 40, 50, 100, 150, 200),
        expand=c(0, 0),
        sec.axis = sec_axis(transform=~.*1/200, name="Probability below 20k metric tons")
      )
  }

  return(
    out
  )
}

