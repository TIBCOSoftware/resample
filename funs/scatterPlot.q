# Copyright 2021. TIBCO Software Inc.
# This file is subject to the license terms contained
# in the license file that is distributed with this file.


##################################################
# scatterPlot
##################################################
scatterPlot <- function(x, ...){
  UseMethod("scatterPlot")
}

##################################################
# scatterPlot.default
##################################################
scatterPlot.default <- function(x, ...){
  # This default method for scatterPlot handles ordinary vectors
  plot(x, ...)
}

##################################################
# scatterPlot.resamp
##################################################
scatterPlot.resamp <- function(x, subset.statistic = 1:2,
			      observed = T,
			      par.observed = list(pch=2, col=5),
			      ...,
			      xlab = statisticNames[1],
			      ylab = statisticNames[2],
			      xlim = NULL, ylim = NULL,
			      main = NULL){
  # Do a scatterplot of two columns of x$replicates
  # If observed = T, show observed using abline
  # par.observed may be NULL, or list of arguments to abline
  statisticNames <- names(x$observed[subset.statistic])
  if(is.null(main)) main <- x$label
  if(is.null(main))
    main <- x$defaultLabel
  X <- x$replicates[ ,subset.statistic[1]]
  Y <- x$replicates[ ,subset.statistic[2]]
  obs <- x$observed[subset.statistic]
  if(is.null(xlim))
    xlim <- if(observed) range(X, obs[1]) else range(X)
  if(is.null(ylim))
    ylim <- if(observed) range(Y, obs[2]) else range(Y)
  plot(X, Y, 
       xlab = xlab, ylab = ylab, main = main, xlim = xlim, ylim = ylim, ...)
  if(observed)
    do.call("abline",
	    c(list(v = obs[1], h = obs[2]),
	      par.observed))
  invisible(NULL)
}

##################################################
# scatterPlot.data.frame
##################################################
scatterPlot.data.frame <- function(x, ...){
  # Find numerical columns.  If only two, then do scatterplot.
  # If more than two, then call pairs.
  x <- x[!sapply(x, is.factor)]
  x <- x[!sapply(x, is.character)]
  p <- ncol(x)
  if(p < 2){
    warning("Only ", p, " columns remaining after excluding character and factor data, cannot do scatter plot")
    return(invisible(NULL))
  }
  if(p == 2)
    plot(x[[1]], x[[2]],
	 ..., xlab = names(x)[1], ylab = names(x)[2])
  # Have ... before xlab in case ... includes xlab
  else
    pairs(x, ...)
  invisible(NULL)
}
