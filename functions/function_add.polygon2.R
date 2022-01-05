# add.polygon2
# Created by Melissa Muradian
# Date updated:  12/14/2015
# Summary:
# This function plots shaded quantiles on a plot of a time series of posterior distributions

add.polygon2 <- function(x, y, z, alpha.level, alpha.min=0, alpha.max=1){
  ## Pass this a series of years (x) and a matrix of trajectories (df) and it
  ## adds a polygon to the current plot with whose area contains 1-alpha.level
  alpha.level <- sort(alpha.level)
  for(i in 1:length(alpha.level)){
    alpha <- alpha.level[i]
    alpha.col <- alpha.min+alpha*(alpha.max-alpha.min)
    col.poly <- rgb(1-alpha.col,1-alpha.col,1-alpha.col, alpha=1)
    quantiles.temp <-  as.matrix(t(apply(z, 2, quantile,
                                         probs=c(alpha/2,1-alpha/2),name=F, na.rm=T)))
    polygon(x=c(x, rev(x)), y=c(quantiles.temp[,1], rev(quantiles.temp[,2])),
            col=col.poly, border=NA)
  }
}

