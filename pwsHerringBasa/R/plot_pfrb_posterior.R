#' Plot Pre-fishery biomass posterior
#'
#' Wrapper function for plotting the model-fit pre-fishery biomass (metric tons)
#' posterior distribution
#'
#' @param df A data frame returned by
#' \link[pwsHerringBasa]{compute.pfrb.posterior} containing model-fit
#' pre-fishery biomass posterior distribution
#' @param quants A numeric vector returned by
#' \link[pwsHerringBasa]{compute.pfrb.posterior} containing .025, .5, and .975
#' quantiles.
#' @param prob A scalar returned by
#' \link[pwsHerringBasa]{compute.pfrb.posterior} giving the probability
#' that pre-fishery biomass lies below the 22,000 short-ton (19,958 metric ton)
#' harvest threshold.
#' @param curr.year Integer giving the current year in the model run
#' @param font.size Plotting parameter passed to \link[ggplot2]{ggplot}
#'
#' @returns A \link[ggplot2]{ggplot} object showing histogram of model-fit
#' pre-fishery biomass (metric tons)
#'

plot_pfrb_posterior <- function(df, quants, prob, curr.year, font.size=1){
  q3 <- q2 <- q1 <- NULL

  extra <- data.frame(q1=quants[1], q2=quants[2], q3=quants[3], prob=prob)
  out <- ggplot(df)+
    geom_histogram(aes(x=.data$biomass/1000, y= after_stat(density)), bins=60)+
    # CLR: changed ..density.. to after_stat(density) to address deprecation warning
    scale_fill_grey(start=0.8, end=0.6)+
    geom_vline(xintercept = quants, linetype=c("dashed", "solid", "dashed"), linewidth=c(0.5, 1, 0.5))+
    geom_text(data=extra, aes(x=75, y=0.05, label=paste("Median:", q2)), size=font.size, hjust=1)+
    geom_text(data=extra, aes(x=75, y=0.035, label=paste0("95% interval:\n", "(", q1, ", ", q3, ")")), size=font.size, hjust=1)+
    geom_text(data=extra, aes(x=75, y=0.02, label=paste("Probability below\nthreshold:", prob)), size=font.size, hjust=1)+
    # scale_x_continuous(paste(curr.year, "spawning biomass (mt)"), breaks=seq(0, 60, 5), expand=c(0, 0))+
    # scale_y_continuous("Probability density", breaks=seq(0, 0.10, 0.025), expand=c(0, 0))+
    # coord_cartesian(ylim=c(0, 0.15), xlim=c(0, 60))+
    xlab(paste(curr.year, "spawning biomass (mt)")) + 
    ylab("Probability density") +
    ggtitle(paste(curr.year, "spawning biomass"),
            subtitle = "Posterior probability density")+
    theme_bw(base_size = 12) +    
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
    )

  return(out)
}
