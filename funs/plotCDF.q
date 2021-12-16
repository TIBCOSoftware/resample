# In this file
#   plotCDF		(generic)
#   plotCDF.default
#   plotCDF.resamp


plotCDF <- function(x, ... )
  UseMethod("plotCDF")


##########################################################
# plotCDF.default
##########################################################
plotCDF.default <- function(x, weights = NULL, cumWeights, new=T, ... ) {
  # Plot cumulative distribution function
  # If cumWeights is supplied x must be sorted, and ignore weights
  n <- length(x)

  if(missing(cumWeights)){
    # Use weights, if supplied
    if(length(weights)){
      if(length(weights) != n )
	stop("length of weights does not match x in plotCDF")
      if(notSorted(x)){
	o <- order(x)
	x <- x[o]
	cumWeights <- cumsum(weights[o]) / sum(weights)
      }
      else
	cumWeights <- cumsum(weights) / sum(weights)
    }
    else {
      cumWeights <- 1:n/n
      if(notSorted(x))
	x <- sort(x)
    }
  }
  else {
    # cumWeights is supplied
    if(notSorted(x))
      stop("x must be sorted in plotCDF if cumWeights supplied")
    if(abs(cumWeights[n] -1) > 1E-6 )
      warning("last value of cumWeights is not 1 in plotCDF")
  }

  # These arguments are par() arguments
  dots <- list(...)
  dots <- dots[ is.element(names(dots), names(par())) ]

  # New plot, or add to plot
  if(new){
    plot(c(x[1], x), c(0, cumWeights), type="s", ..., xlab="", ylab="")
    if(is.null(list(...)$xlab))
      title(xlab=as.character(substitute(x)))
    if(is.null(list(...)$ylab))
      title(ylab="CDF")
  }
  else
    do.call("lines", c(list(c(x[1], x), c(0, cumWeights), type="s"), dots))
  # Add orizontal lines at both ends
  do.call("segments", c(list(par("usr")[1:2], 0:1, x[c(1,n)], 0:1), dots))
  invisible(list(x=x, y=cumWeights))
}

# x <- runif(10)
# plotCDF(x)
# plotCDF(x,weights=1:10)
# plotCDF(x,weights=10:1)
# cw <- c(sort(runif(9)), 1)
# plotCDF(sort(x),cw=cw)



##########################################################
# plotCDF.resamp
##########################################################
plotCDF.resamp <- function(x, weights = x$weights,
			   nrow = NULL, grid.layout = T,
			   subset.statistic = 1:p,
			   ..., main=NULL){
  # plotCDF plots of resample replicates
  p <- ncol(x$replicates)
  P <- length(subset.statistic)
  if(is.character(subset.statistic)){
    j <- match(subset.statistic, names(x$observed))
    if(anyMissing(j))
      stop("subset.statistic names not recognized:\n    ",
	   paste(subset.statistic[which.na(j)], collapse=", "))
    subset.statistic <- j
  }
  if(grid.layout && P != 1) {
    if(is.null(nrow)){
      if(length(x$dim.obs) == 2)
	old.par <- par(mfcol = x$dim.obs)
      else {
        nrow <- min(P, 2)
        old.par <- par(mfrow = c(nrow, ceiling(P/nrow)))
      }
    }
    else
      old.par <- par(mfrow = c(nrow, ceiling(P/nrow)))
    on.exit(par(old.par))
  }
  # Main title; if not supplied then use x$label or x$defaultLabel
  if(is.null(main))
    main <- x$label
  if(is.null(main))
    main <- x$defaultLabel

  for(i in subset.statistic) {
    xi <- x$replicates[, i]
    plotCDF(xi,
	   ylab = names(x$observed)[i], weights = weights, ...,
	   main=if(is.null(main)) names(x$observed)[i] else main)
    dots <- list(...)
    if(!is.element("axes", names(dots)) || dots$axes){
      if(length(dots)){
	do.call("box", dots[is.element(names(dots), names(par()))])
	# like box(...), but only pass par parameters like bty, col
      }
      else
	box()
    }
  }
  invisible(NULL)
}
# Copied from qqnorm.resamp, minor modifications
# Use box(...) instead of box(), to support bty


#set.seed(0)
#x <- runif(9)
#y <- cbind(x, y = runif(9))
#boot1 <- bootstrap(x, mean, trace=F, seed=1)
#boot2 <- bootstrap(y, colMeans, trace=F, seed=1)
#plotCDF(boot1)
#plotCDF(boot2)
