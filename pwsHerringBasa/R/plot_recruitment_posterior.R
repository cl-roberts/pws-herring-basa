#' Plot Age-3 Recruitments
#'
#' Wrapper function for plotting the model-estimated age-3 recruits (in millions
#' of fish) for the full model time series.
#'
#' @param df A data frame (or \link[tibble]{tibble}) returned by
#' \link[pwsHerringBasa]{compute.recruitment} containing model-estimated
#' age-3 recruits with credible intervals
#' @param years A integer vector giving the years in model time series
#' @param legend Logical; displays plot legend if TRUE
#'
#' @returns A \link[ggplot2]{ggplot} object showing showing model-estimated
#' time series of age-3 recruits with credible intervals
#'

plot_recruitment_posterior <- function(df, years, legend=TRUE){
  out <- ggplot(df) +
    ggdist::geom_lineribbon(aes(x=.data$year, y=.data$recruits, ymin=.data$.lower, ymax=.data$.upper, group=1), size=0.25)+
    scale_fill_grey(start=0.8, end=0.6)+
    scale_x_discrete("Year", breaks=seq(min(years), max(years), by=5), expand=c(0, 0))+
    scale_y_continuous("Age-3 recruits (millions)", breaks=seq(0, 2000, by=500), limits=c(0, 2000), expand=c(0, 0))+
    ggtitle("Age-3 Recruitment")+
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      legend.position = ifelse(legend, "right", "none")
    )

  return(out)
}
