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

    survey.data$year <- as.numeric(survey.data$year)  

    out <- ggplot(fits) +
      geom_ribbon(aes(x=year, ymin=lower_50/scale, ymax=upper_50/scale), alpha = .25) +
      geom_ribbon(aes(x=year, ymin=lower_95/scale, ymax=upper_95/scale), alpha = .25) +
      geom_line(aes(x=year, y = data/scale), linewidth=0.25) +
      geom_point(data=survey.data, aes(x=.data$year, y=data/scale), na.rm=TRUE) +
      geom_line(data=survey.data, aes(x=.data$year, y=data/scale, group=1), na.rm=TRUE)+
      scale_x_continuous("Year", breaks=seq(min(survey.data$year), max(survey.data$year), by=5))+
      scale_y_continuous(expand=c(0, 0))+
      labs(y=ylabel, title=title)+
      coord_cartesian(xlim=range(survey.data$year), ylim=c(0, y.max/scale))+
      theme_bw(base_size = 12) +
      theme(
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
      )

  if(cvs){
    out <- out + 
      geom_errorbar(data=survey.data, aes(x=year, ymin=lower/scale, ymax=upper/scale), width=0)
  }

  return(out)
}
