# like tabPlot.bootstrap, but add arguments "new" and ...

tabPlot.resamp <- 
function(obj, plotHist = F, plotQQ = F, plot.both = F, new=T, ..., main = NULL,
	 subset.statistic = seq(length=p))
{
  p <- ncol(obj$replicates)
  if(!length(subset.statistic)){
    if(plotHist || plotQQ || plot.both)
      warning("statistic is length 0, skipping plots")
    return(invisible(obj))
  }
  if((plotHist || plotQQ || plot.both) && new)
    new.graphsheet()
  if(plotHist)
    plot(obj, grid.layout = F, ..., subset.statistic = subset.statistic)
  if(plotQQ)
    qqnorm(obj,  grid.layout = F, ..., subset.statistic = subset.statistic)
  if(plot.both){
    oldpar <- par(mfrow=c(1,2))
    on.exit(par(oldpar))
    haveMain <- (is.null(main) || main != "")
    if(haveMain){
      cexMain <- 1.5 * par("cex")
      par(oma = c(0, 0, .1 + cexMain, 0) + par("oma"))
      if(is.null(main))
	main <- obj$label
      if(is.null(main))
	main <- obj$defaultLabel
    }
    for(j in subset.statistic){
      plot(obj, grid.layout = F, ..., subset.statistic = j, main = "")
      qqnorm(obj, grid.layout = F, ..., subset.statistic = j, main = "")
      if(haveMain)
	mtext(side=3, outer=T, line = -2.1, cex = cexMain, text=main)
    }
  }
  invisible(obj)
}
# do not abbreviate grid.layout
# treat plot and qqnorm consistently for grid.layout
# add argument subset.statistic
# plot hist & qq together, when multiple statistics    
# single title when plotting hist and qq together
# skip plots if none to be produced

tabPlot.bootstrap <- function(obj, ...){
  # this function is deprecated, use tabPlot.resamp directly.
  tabPlot.resamp(obj, ...)
  invisible(obj)
}
tabPlot.jackknife <- tabPlot.bootstrap


"Done with tabplot.q" # for script file input
