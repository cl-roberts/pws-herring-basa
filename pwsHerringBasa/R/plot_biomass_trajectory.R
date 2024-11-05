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
    ggdist::geom_lineribbon(aes(x=.data$year, y=.data$biomass/1000,
                                ymin=.data$.lower/1000, ymax=.data$.upper/1000, group=1),
                            size=0.25)+
    scale_fill_grey(start=0.8, end=0.4)+
    scale_x_discrete("Year", breaks=seq(min(years), max(years), by=5), expand=c(0, 0))+
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
    out <- out + geom_point(aes(x=.data$year, y=.data$prob*200))+
      geom_line(aes(x=.data$year, y=.data$prob*200, group=1))+
      scale_y_continuous(
        "Pre-fishery biomass (1000 mt)",
        breaks=c(0, 20, 40, 50, 100, 150, 200),
        expand=c(0, 0),
        sec.axis = sec_axis(trans=~.*1/200, name="Probability below 20k metric tons")
      )
  }

  return(
    out
  )
}
