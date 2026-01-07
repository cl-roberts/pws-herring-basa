#' Plot survey fits
#'
#' Plots BASA model fits to survey values with posterior predictive interval with
#' \link[ggplot2]{ggplot}
#'
#' @param fits Model fits of survey data
#' @param survey.data Observed values from surveys
#' @param y.max Upper limit of y-axis, passed to \link[ggplot2]{ylim}
#' @param title Plot title, passed to \link[ggplot2]{labs}
#' @param ylabel y-axis label (usually name of survey), passed to
#' \link[ggplot2]{labs}
#' @param scale Scaling factor for dependent variable (survey data)
#' @param cvs Logical indicating whether survey CV's are included in plot
#'
#' @returns A \link[ggplot2]{ggplot} object
#'


plot_survey_fits <- function(fits, survey.data, y.max, title, ylabel="",
                             scale=1, cvs=TRUE){

  if(cvs){
    points <- ggdist::geom_pointinterval(data=survey.data,
                                         aes(x=.data$year, y=data/scale,
                                             ymin=.data$lower/scale, ymax=.data$upper/scale))
  }else{
    points <- geom_point(data=survey.data, aes(x=.data$year, y=data/scale))
  }

  return(
    ggplot(fits) +
      ggdist::geom_lineribbon(aes(x=.data$year, y=data/scale, ymin=.data$.lower/scale, ymax=.data$.upper/scale, group=1), linewidth=0.25)+
      scale_fill_grey(start=0.8, end=0.6) +
      points +
      geom_line(data=survey.data, aes(x=.data$year, y=data/scale, group=1))+
      scale_x_discrete("Year", breaks=seq(min(survey.data$year), max(survey.data$year), by=5))+
      scale_y_continuous(expand=c(0, 0))+
      labs(y=ylabel, title=title)+
      coord_cartesian(xlim=c(1, length(survey.data$year)), ylim=c(0, y.max/scale))+
      theme_bw(base_size = 12) +
      theme(
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
      )
  )
}
