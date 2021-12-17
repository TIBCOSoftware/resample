# Copyright 2021. TIBCO Software Inc.
# This file is subject to the license terms contained
# in the license file that is distributed with this file.

# bootstrapT and its methods, and pivot functions

# resampPivotT
# resampPivotDiff
# resampPivotRatio
# bootstrapT
# bootstrapT.default
# bootstrapT.bootstrap


##########################################################
# Pivots
##########################################################
# bootstrapT expects that pivot is a list containing two functions,
#   * pivot()
#   * inverse()

# pivot() takes arguments replicates[B,k] and observed[k] where
# B is the number of bootstrap replications and k is typically 2 (est. & se)
# It returns a matrix[B,K] (or a vector[B]?).
# For example
#	replicates[,1] = bootstrap sample means "Xbar"
#	replicates[,2] = bootstrap sample stdevs "S"
#	observed = c(xbar, s)                          # observed data
#	pivot( ) = (Xbar - xbar) / S
# In this case k=2 and K=1.  K is the number of pivotal statistics,
# which should equal the number of parameters of interest.

# inverse() solves for the parameter(s) of interest.
# It takes arguments observed[k] and quantiles[K,r], where r=length(probs)
# and quantiles contains the (1-probs) quantiles of K pivotal statistics.
# It returns a matrix theta[K,r] that corresponds to solving
#	pivot(observed, theta+junk) = quantiles
# in terms of the true parameters of interest -- in this case mu but not sigma,
#	inverse(c(xbar,s), quantiles) = xbar - quantiles * s
# It should also add appropriate row names.

resampPivotT <-
  list(pivot = function(replicates, observed){
	 # assume estimates are in odd columns and stderrs in even columns
	 # Compute:  (estimate-observed)/stderr; 
	 K <- trunc(length(observed) / 2)
	 if(K < 1)
	   stop("length of statistic is only ", length(observed),
		".\nPerhaps you forgot to include both an estimate and standard error?")
	 odd <- seq(1, by=2, length=K)
	 (replicates[, odd, drop=F] -
	  rep(observed[odd], each=nrow(replicates))
	  )/replicates[, 1+odd, drop=F]
       },
       inverse = function(observed, quantiles){
	 # Inverse of the default pivot function.
	 # Given pivot = quantiles, solve for observed.
	 # Return matrix with r=length(quantiles) rows and K columns.
	 K <- nrow(quantiles)
	 odd <- seq(1, by=2, length=K)
	 # order is important in next line, so result is a matrix
	 - quantiles * observed[odd+1] + observed[odd]
       })


resampPivotDiff <-
  list(pivot = function(replicates, observed) scale(replicates, observed, F),
       inverse = function(observed, quantiles) (-quantiles) + observed)

resampPivotRatio <-
  list(pivot = function(replicates, observed) scale(replicates, F, observed),
       inverse = function(observed, quantiles) 1/quantiles * observed)

# Use unusual order of operations to work around bug #18178 x+y and y+x differ.


##########################################################
# Bootstrap t
##########################################################
# Support for bootstrapping pivots.
# First call bootstrap to get the statistic and standard error,
# then call this.
# Or call this directly; but then replicates aren't saved.


# ToDo (for bootstrap t):
#   allow cases without analytical se's.
#   incorporate within bootstrap? (have bootstrapT.default instead)
#   add a class and print method? (I prefer not to)


#####################
bootstrapT <-
function(data, ...)
UseMethod("bootstrapT")
# Change name of argument to data.  Similarly for methods.

#####################
# Do bootstrap t for a vector or data frame.
# Default statistic is the mean.
bootstrapT.default <-
function(data, statistic = NULL, ...,
	 pivot = resampPivotT,
	 probs = c(25, 50, 950, 975)/1000)
{
  # This default method performs bootstrap T calculations for raw data
  # First call bootstrap(data, statistic, ...)
  #   (with a default statistic that calls colMeans() and colStdevs())
  # Then call bootstrapT.bootstrap for pivot calculations
  m <- match.call()
  if(is.null(m$statistic))
    m$statistic <- function(x) rbind(mean = colMeans(x), stdev = colStdevs(x))
  m[[1]] <- as.name("bootstrap")
  m$pivot <- NULL
  m$probs <- NULL
  boot.obj <- eval(m, sys.parent())
  bootstrapT.bootstrap(boot.obj, pivot=pivot, probs=probs)
}
# Arguments pivot and probs rather than t.args
# Add argument statistic (so user can supply by position)


#####################
# Do bootstrap t based on a bootstrap object.
# For the default pivot, replicates should have 2p columns, odd columns
# are estimates and even columns are standard errors or std deviations

bootstrapT.bootstrap <-
function(data,
	 pivot = resampPivotT,
	 probs = c(25, 50, 950, 975)/1000)
{
  # Perform bootstrap T calculations pivot calculations, based
  # on a bootstrap object.
  u <- pivot$pivot(data$replicates, data$observed)
  if(!is.matrix(u))
    u <- as.matrix(u)
  r <- length(probs)
  K <- numCols(u)
  qprobs <- 1-probs  # upper quantiles of bootstrap t distn give lower CI's

  # Create matrix of quantiles, K rows and r columns
  quantiles <- matrix(NA, K, r)
  for(j in 1:K)
    quantiles[j,] <- quantile(u[,j], probs=qprobs, na.rm=T,
			 alpha=0, weights=data$weights)
  dimnames(quantiles) <- list(dimnames(u)[[2]], names(quantile(1:3, probs=qprobs)))

  CI <- pivot$inverse(data$observed, quantiles)

  # Use column names as formatted by quantile
  dimnames(CI) <- list(dimnames(CI)[[1]], names(quantile(1:3, probs)))

  list(pivot = quantiles,
       limits = CI)
}
# name of first return component now pivot, not t-quantiles
# pivot is now a list of two functions, pivot$pivot and pivot$inverse
# First argument is "data", for ease of use, even though it is a bootstrap obj



