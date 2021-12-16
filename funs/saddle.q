#
# Functions saddlepointP
#           saddlepointD
#           saddlepointPSolve
#           pDiscreteMean
#           dDiscreteMean
#           qDiscreteMean
#
# For computing the saddlepoint approximations of the density or cumulative
# distribution function for the mean of `size' observations drawn from `L'.
#
# saddlepointP: given `tau', find the values of the cdf.
# saddlepointD: given `tau', find the values of the density.
# pDiscreteMean: given quantiles, find the values of the cdf.
# dDiscreteMean: given quantiles, find the values of the density.
# qDiscreteMean: given probabilities (cdf), find the quantiles.
# saddlepointPSolve: given probabilties (cdf), find the `tau'.
#

##################################################
# saddlepointP
##################################################
saddlepointP <- function(tau, L, size = NULL, weights = NULL, group=NULL,
			mean = T,
                        conv.factor = 0)
{
  # Calculate a saddlepoint estimate of lower tail probability,
  #    P( Lbar < tiltMean(tau, L))
  # where Lbar is the mean of "size" random observations from L.
  # tau is on the scale of L, not the sample mean.
  #
  # If multiple groups and mean=F, then estimate
  #    P( sum of n[g] Lbar[g] < sum of n[g] tiltMean(tau, L[group g] )
  #
  # If multiple groups and mean=T, then estimate
  #    P( sum of Lbar[g] < tiltMean(tau, L, group)
  #                      = sum of tiltMean(tau, L[group g] )
  # In this case, the tilting parameter tau / (n[g]/n) is used for group g.

  n <- length(L)
  nw <- length(weights)
  if(nw){
    if(nw != n)
      stop("length of weights must be same as length of L")
    if(min(weights) < 0)
      stop("weights must be nonnegative")
  }

  if(byGroup <- length(group)){
    if(byGroup != n)
      stop("length of group must be same as length of L")
    group.inds <- split(1:n, group)
    nGroups <- length(group.inds)
    groupSizes <- unlist(lapply(group.inds, length))
    if(notSorted(ord <- unlist(group.inds))){
      # Reorder so that groups are contiguous
      L <- L[ord]
      weights <- weights[ord]
    }
    if(is.null(size))
      size <- n
    else
      stop("size argument not supported with stratified data")
  }
  else{
    nGroups <- 1
    groupSizes <- n
    if(is.null(size))
      size <- n
    else
      if(size <= 0) stop("size must be a positive integer")
  }
  .C("S_saddlepointP",
     as.double(tau),
     as.integer(length(tau)),
     as.double(L),
     as.integer(n),
     as.integer(size),
     as.double(weights),
     as.integer(nw),
     as.integer(nGroups),
     as.integer(groupSizes),
     as.integer(mean),
     as.double(conv.factor),
     p = double(length(tau)),
     COPY = c(F,F,
       (nGroups > 1 && mean),  # L is modified by the C code
       F,F,F,F,F,F,F,T))$p
}
# add size argument
# add COPY argument to .C call
# add conv.factor argument
# change value of size to n in stratified case (for C code)

##################################################
# saddlepointD: given tau, find the density
##################################################
saddlepointD <- function(tau, L, size = NULL, weights = NULL, group = NULL,
                        mean = T,
			conv.factor = 0)
{
  # Calculate a saddlepoint density estimate for:
  #   if mean=T:  the mean of size observations from L
  #   if mean=F:  the sum  of size observations from L 
  # In either case, tau is on the scale of L, not the mean or sum.
  #
  # If multiple groups and mean=F, then estimate density for
  # the sum of random group sums (of size n[g] from group g) at
  #   sum of n[g] tiltMean(tau, L[group g]
  #
  # If multiple groups and mean=T, then estimate density for
  # the sum of group means (of size n[g] from group g) at
  #   tiltMean(tau, L, group) (equals sum of tiltMean(tau, L[group g] ))
  # In this case, the tilting parameter tau / (n[g]/n) is used for group g.

  n <- length(L)
  nw <- length(weights)
  if(nw){
    if(nw != n)
      stop("length of weights must be same as length of L")
    if(min(weights) < 0)
      stop("weights must be nonnegative")
  }

  if(byGroup <- length(group)){
    if(byGroup != n)
      stop("length of group must be same as length of L")
    if(!is.null(size))
      stop("size argument not supported with stratified data")
    else
      size <- n
    group.inds <- split(1:n, group)
    nGroups <- length(group.inds)
    groupSizes <- unlist(lapply(group.inds, length))
    if(notSorted(ord <- unlist(group.inds))){
      # Reorder so that groups are contiguous
      L <- L[ord]
      weights <- weights[ord]
    }
  }
  else{
    nGroups <- 1
    groupSizes <- n
    if(is.null(size))
      size <- n
    else
      if(size <= 0 ) stop("size must be a positive integer")
  }
  .C("S_saddlepointD",
     as.double(tau),
     as.integer(length(tau)),
     as.double(L),
     as.integer(n),
     as.integer(size),
     as.double(weights),
     as.integer(nw),
     as.integer(nGroups),
     as.integer(groupSizes),
     as.integer(mean),
     as.double(conv.factor),
     p = double(length(tau)),
     COPY = c(F,F,
       (nGroups > 1 && mean),  # L is modified by the C code
       F,F,F,F,F,F,F,T))$p
}


##################################################
# saddlepointPSolve
##################################################
saddlepointPSolve <-
  function(probs, L, size = NULL, weights = NULL, group = NULL,
	   mean = T,
           conv.factor = 0, initial, tol = 1E-6, tol.tau = tol,
           maxiter = 100)
{
  # Find tau that solves saddlepointP(tau, ...) = probs

  n <- length(L)
  nw <- length(weights)
  if(nw){
    if(nw != n)
      stop("length of weights must be same as length of L")
    if(min(weights) < 0)
      stop("weights must be nonnegative")
  }

  k <- length(probs)
  if(missing(initial) || length(initial) != k)
    initial <- qnorm(probs) / sqrt(var(as.vector(L), SumSquares=T))

  if(byGroup <- length(group)){
    if(byGroup != n)
      stop("length of group must be same as length of L")
    if(is.null(size))
      size <- n
    else
      stop("size argument not supported with stratified data")
    group.inds <- split(1:n, group)
    nGroups <- length(group.inds)
    groupSizes <- unlist(lapply(group.inds, length))
    if(notSorted(ord <- unlist(group.inds))){
      # Reorder so that groups are contiguous
      L <- L[ord]
      weights <- weights[ord]
    }
  }
  else{
    nGroups <- 1
    groupSizes <- n
    if(is.null(size))
      size <- n
    else
      if(size <= 0 ) stop("size must be a positive integer")
  }
  .C("S_saddlepointPSolve",
     as.double(probs),
     as.integer(k),
     as.double(L),
     as.integer(n),
     as.integer(size),
     as.double(weights),
     as.integer(nw),
     as.integer(nGroups),
     as.integer(groupSizes),
     as.integer(mean),
     as.double(conv.factor),
     as.double(initial),
     as.double(tol),
     as.double(tol.tau),
     as.integer(maxiter),
     tau = double(k),
     COPY = c(F,F,
       (nGroups > 1 && mean),  # L is modified by the C code
       F,F,F,F,F,F,F,F,F,F,F,T))$tau
}
# add size argument
# add COPY argument to .C call
# add conv.factor argument
# change value of `size' to n in stratified case (for C code)



##################################################
# pDiscreteMean: given quantiles, find the probabilities
##################################################
pDiscreteMean <- function(q, values, size = NULL, weights = NULL, group = NULL,
			 conv.factor = 0, ...){
  # note:  ... arguments are passed to tiltMeanSolve

  tau <- tiltMeanSolve(q = q, L = values, tilt = "exponential",
		       weights = weights,
                       group = group, ...)$tau.exp
  saddlepointP(tau = tau, L = values, size = size, weights = weights,
               group = group, conv.factor = conv.factor)
}
# Change name, pSaddlepoint to pDiscreteMean
# Change argument from L to values, for consistency with [pdqr]discrete

##################################################
# qDiscreteMean: given probabilities, find the quantiles
##################################################
qDiscreteMean <- function(p, values, size = NULL, weights = NULL, group = NULL,
			 conv.factor = 0, ...){
  # note:  ... arguments passed to saddlepointPSolve

  tau <- saddlepointPSolve(probs = p, L = values, size = size,
			   weights = weights,
                           group = group, conv.factor = conv.factor, ...)
  tiltMean(tau = tau, L = values, tilt = "exponential", weights = weights,
           group = group)$q
}
# Change name, qSaddlepoint to qDiscreteMean
# Change argument from L to values, for consistency with [pdqr]discrete

##################################################
# dDiscreteMean: given quantiles, find the density
##################################################
dDiscreteMean <- function(x, values, size = NULL, weights = NULL, group = NULL,
			 conv.factor = 0, ...){
  # note:  ... arguments passed to tiltMeanSolve

  tau <- tiltMeanSolve(q = x, L = values, tilt = "exponential",
		       weights = weights,
                       group = group, ...)$tau.exp
  saddlepointD(tau = tau, L = values, size = size, weights = weights,
               group = group, conv.factor = conv.factor)
}
# Change name, dSaddlepoint to dDiscreteMean
# Change argument from L to values, for consistency with [pdqr]discrete
