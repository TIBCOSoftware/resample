# quantile

# Features:
#   weights
#   alpha (1: narrowest definition of quantiles
#          0: widest),
#   rule   choice of rules for quantile tails out of the range of
#	   x (1: NA , 2: extreme x values, 3: extrapolation)
#   freq   frequencies

quantile.resamp <- 
function(x, probs= 0:4/4, na.rm= F, ..., weights = x$weights){
  # x is a resample object; use the $replicates
  # q may be a vector, or matrix with p columns, 
  # where p = ncol(x$replicates).
  # Result is matrix with p columns.
  #
  # The default value of alpha = 1 gives quantiles which are too narrow.
  # alpha = 0 is better, and is used by limits.percentile
  xr <- x$replicates
  if(ncol(xr) > 1){
    result <- apply(xr, 2, quantile,
		     probs = probs, na.rm = na.rm, weights = weights,
		     ...)
    dimnames(result)[[2]] <- dimnames(xr)[[2]]
    result
  }
  else
    quantile(xr, probs = probs, na.rm = na.rm, weights = weights, ...)
}
