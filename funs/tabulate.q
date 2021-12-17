# Copyright 2021. TIBCO Software Inc.
# This file is subject to the license terms contained
# in the license file that is distributed with this file.

# Have removed tabulate -- that version is now in S+.

colTabulate <- function(bin, nbins, weights = NULL)
{
  # bin should be a matrix (n rows), weights a vector of length n or 0
  # Result should be equivalent to apply(bin, 2, tabulate, ...)
  if(anyMissing(bin))
    stop("missing values in bin are not supported")
  if(!is.matrix(bin))
    dim(bin) <- c(length(bin), 1)
  n <- length(weights)
  if(n){
    if(nrow(bin) != n)
      stop("Different lengths for bin and weights")
    if(anyMissing(weights))
      stop("missing values in weights are not supported")
  }
  else
    n <- nrow(bin)
  p <- ncol(bin)
  if(!n)
    return(if(missing(nbins)) matrix(integer(0), 0, p)
	   else matrix(0, nbins, p))
  binRange <- range(bin, na.rm = T)
  if(binRange[1] < 1)
    stop("values of bin must be positive integers")
  if(missing(nbins))
    nbins <- binRange[2]
  else if(nbins < binRange[2])
    stop("nbins must be at least as large as max(bin)")
  # C code does not handle NA, but already removed, don't check again.
  if(length(weights)){
    result <- 
      .C("S_tabulate_weights",
	 as.integer(n * p),
	 as.integer(bin + nbins * rep(0:(p-1), each=n)),
	 rep(as.double(weights), p),
	 ans = double(nbins * p),
	 NAOK = T, specialsok = T,
	 COPY = c(F, F, F, T))$ans
  }
  else {
    # To avoid loop, add nbins * (j-1) to j'th column of bins
    result <- .C("S_tabulate",
	   as.integer(n*p),
	   as.integer(bin + nbins * rep(0:(p-1), each=n)),
	   ans = integer(nbins*p),
	   NAOK = T, specialsok = T,
	   COPY = c(F, F, T))$ans
  }
  dim(result) <- c(nbins, p)
  result
}
