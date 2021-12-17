# Copyright 2021. TIBCO Software Inc.
# This file is subject to the license terms contained
# in the license file that is distributed with this file.

# In this file:
#  cdf
#  cdf.resamp
#  cdf.default
#
# pdiscrete, ddiscrete, qdiscrete, rdiscrete are now in pdiscrete.q


cdf <- function(x, q, ...)
  UseMethod("cdf")

cdf.resamp <- function(x, q, weights = x$weights, ...){
  # x is a resample object; use the $replicates
  # q may be a vector, or matrix with p columns,
  # where p = ncol(x$replicates).
  # Result is matrix with p columns.
  xr <- x$replicates
  if(is.matrix(q) && ncol(q) == ncol(xr))
    result <- sapply(1:ncol(q),
		    function(j, x, q, ...) cdf(x[,j], q[,j], ...),
		    x = xr, q = q, weights = weights, ...)
  else
    result <- sapply(1:ncol(xr),
		    function(j, x, q, ...) cdf(x[,j], q, ...),
		    x = xr, q = q, weights = weights, ...)
  if(!is.matrix(result))
    result <- matrix(result, ncol = ncol(xr))
  dimnames(result) <- list(NULL, dimnames(xr)[[2]])
  result
}

cdf.default <- function(x, q, weights=NULL, na.rm=F, normalize=T) {
  # F(q), where F is the empirical distribution function of data x.
  # Result matches length of q.
  # If normalize=F, then how many x values are <= each value of q.
  # If normalize=F, then out[i] = (how many x values are <= q[i])
  if(!length(q))
    return(numeric(0))
  wna <- which.na(x)
  nWeights <- length(weights)
  if(nWeights){
    if(nWeights != length(x))
      stop("length of weights must match length of x")
    wna <- unique(c(wna, which.na(weights)))
  }
  if(length(wna)){
    if(!na.rm)
      stop("missing values in x and/or weights")
    if(nWeights)
      weights <- weights[-wna]
    x <- x[-wna]
  }
  if(nWeights){
    o <- sort.list(c(x, q))
    i <- rep(c(T,F), c(length(x), length(q)))[o]
    wi <- c(weights, rep(0, length(q)))[o]
    i <- cumsum(wi)[!i]
    if(normalize)
      i <- i / sum(weights)
  } else {
    i <- rep(c(T,F), c(length(x), length(q)))[sort.list(c(x, q))]
    i <- cumsum(i)[!i]
    if(normalize)
      i <- i/length(x)
  }
  i[sort.list(q)] <- i
  if(anyMissing(q))
    i[which.na(q)] <- NA
  i
}
