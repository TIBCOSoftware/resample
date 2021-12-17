# Copyright 2021. TIBCO Software Inc.
# This file is subject to the license terms contained
# in the license file that is distributed with this file.

# Bootstrap samplers
#	samp.bootstrap
#	samp.MonteCarlo (same as samp.bootstrap)
#       samp.boot.mc	(same as samp.bootstrap)
#	samp.boot.bal
#	samp.permute
#	samp.permute.old
#	samp.combinations
#	samp.blockBootstrap
#	blockBootstrap
#       samp.bootknife
#       samp.finite
#       samp.half
#
# Other sampling functions
#	balancedSample

# The "samp." functions are bootstrap samplers, with arguments:
#	n  (now a scalar, earlier was a vector)
#	B
#	size = n
# Return [size,B] array of numbers from 1:n

# Common additional arguments:
#	reduceSize (e.g. 1 to reduce size in each stratum by 1)
#	prob
# Additional arguments currently present in one or more functions:
#	method
#	partial
#	blockLength
#       random
#	k
#
# Function balancedSample is not a bootstrap sampler, but is included 
# here because it shares functionality with samp.permute.


##########################################################

samp.bootstrap <-
  function(n, B, size = n - reduceSize, reduceSize = 0, prob = NULL)
{
# Simple Monte Carlo bootstrap.  Each column is one sample.
  if (length(prob))
    x <- sample(x = n, size = B*size, replace = T, prob = prob)
  else
    x <- sample(x = n, size = B*size, replace = T)
  dim(x) <- c(size, B)
  x
}
# allow weights argument
# changed `weights' to `prob'
# add argument reduceSize


##########################################################
# samp.boot.mc and samp.MonteCarlo retained for backwards compatibility
#
samp.MonteCarlo <- samp.bootstrap
samp.boot.mc <- samp.bootstrap


##########################################################

samp.boot.bal <-
  function(n, B, size = n - reduceSize, reduceSize = 0, method = "biased")
{
  # Balanced bootstrap.  Each column is one sample.
  method <- match.arg(method, c("biased", "unbiased", "semi"))
  if(method == "biased"){
    # Usual balanced bootstrap.  Bias of order O(1/B).
    # Bias is due to dependence between rows.
    if((size * B) %% n)
      stop("size * B must be a multiple of n for method = 'biased'")
    x <- sample(rep(1:n, size*B/n), replace = F)
    dim(x) <- c(size, B)
    return(x)
  }

  # Define a function for use below
  permuteCols <- function(x){
    y <- .C("S_samp_permute",
	   nrow(x),
	   ncol(x),
	   inds = as.integer(x),
	   COPY = c(F,F,T))$inds
    dim(y) <- dim(x)
    y
  }

  # If n divides B, then the remaining methods use the same algorithm,
  # which is unbiased and balanced:  place a balanced set of observations
  # in each column, reshuffle each column, then transpose.
  remainder <- B %% n
  if(!remainder){
    x <- rep(1:n, length = B * size)
    dim(x) <- c(B, size)
    return(t(permuteCols(x)))
  }

  if(method == "unbiased"){
    # Unbiased (independence between samples); not exactly balanced
    return(t(samp.permute(n, B = size, size = B)))
  }

  # method = "semi"
  # Balanced (if n divides size*B), biased O(k/B^2) where k=(B%%n)
  if(B == 1)
    return(samp.permute(n, B = 1, size = size)) 

  k <- B %/% n
  x1 <- rep(rep(1:n, k), size)
  dim(x1) <- c(n*k, size)
  # if n divides size*B, then n also divides remainder*size, and 
  # x2 will be balanced.
  x2 <- balancedSample(n = n, size = remainder*size)  
  dim(x2) <- c(remainder, size)
  t(permuteCols(rbind(x1, x2)))
}
# Added argument "method"; rm(samp.boot.bal.[unbiased,semi])
# Added argument "size"
# Rewrote to use samp.permute; rm(samp.permProbs)
# add argument reduceSize

# Need a tech report on balanced bootstrapping.

#| temp <- samp.boot.bal(5, 20, "unbiased")
#| temp
#| tabulate(temp, nbins=5)
#| apply(temp, 1, tabulate, nbins=5)
#| temp <- samp.boot.bal(5, 17, "unbiased")
#| temp
#| tabulate(temp, nbins=5)  # not balanced
#| apply(temp, 1, tabulate, nbins=5)
#|
#| temp <- samp.boot.bal(5, 20, "semi")
#| temp
#| tabulate(temp, nbins=5)
#| apply(temp, 1, tabulate, nbins=5)
#| temp <- samp.boot.bal(5, 17, "semi")
#| temp
#| tabulate(temp, nbins=5)
#| apply(temp, 1, tabulate, nbins=5)


##########################################################

samp.permute <-
function(n, B, size = n - reduceSize, reduceSize = 0, prob = NULL,
	 full.partition = "none", partial = NULL)
{
  # If only n and B are supplied, return n x B matrix of random
  # permutations of 1:n.  Each column is one permutation.
  # See help file for other options.
  # Argument partial is ignored, is here for compatibility.
  #
  # Easy case: no probabilities and n divides size, in which case 
  # generate B permutations of size/n copies of 1:n, 
  if(!(nprob <- length(prob)) && !(size %% n) )
    x <- .C("S_samp_permute", 
	   as.integer(size),
	   as.integer(B),
	   inds = as.integer(rep(1:n,B*size/n)),
	   COPY = c(F,F,T))$inds
  else{
    # Each column is a "balanced" sample of size `size'
    # from numbers 1:n, with probabilities proportional to `prob'.
    #
    # Let K = size * prob / sum(prob) = expected frequences (in each column)
    # "Balanced" means that obs j is included either trunc(K[j]) or
    # ceiling(K[j]) times.
    if(nprob){
      if(n != nprob)
	stop("n must match length(prob)")
      prob <- prob / sum(prob)
    }
    if(size >= n)
      full <- full.last <- F
    else{
      full.partition <- match.arg(full.partition, c("first", "last", "none"))
      if(full.partition == "first"){
        full <- T
        full.last <- F
      }
      else if(full.partition == "last"){
        full <- full.last <- T
      }
      else
        full <- full.last <- F
    }
    x <- .C("S_samp_perm_probs", 	
	   as.integer(size),
	   as.integer(n),
	   as.integer(nprob),
	   as.double(prob / sum(prob)),
           as.integer(full),
           as.integer(full.last),
	   as.integer(B),
	   index = integer(size*B),
	   COPY = c(F,F,F,F,F,F,F,T))$index
  }
    
  dim(x) <- c(size, B)
  x
}
# Revised version, address bug 18227.
# 	version 9 from ~timh/bootstrap/samp.permute.S
# 	I earlier checked in version 8 to the Splus tree.
# Revised to use simple transposition algorithm.  
# 	Argument partial now ignored, since the algorithm is already O(n).
# Functionality combined with samp.permProbs (rm(samp.permProbs)), 
#       added size argument
# NOTE: The case of unequal probs (here and in balancedSample) was written 
#       to replace "sample" in the case of replace=F and unequal probs 
#       (see bug #15647), and the previous version of balancedSample defined 
#       in crossValidation.q. (RT, 5/16/01) 
# add argument reduceSize


##########################################################

samp.permute.old <-
function(n, B, partial = NULL)
{
  # This function is provided in case you need to reproduce results
  # from the version of samp.permute in Splus6.0 and earlier.
  # The new version is faster.
  # This is deprecated, and may be removed in future versions.
  #
  # Random permutations of 1:n.  Each column is one permutation.
  f <- function(i, n, x = runif(n), partial = NULL)
  .Internal(sort.list(x, partial), "S_sort_list", T, 0)
  x <- unlist(lapply(1:B, f, n = n, partial = partial))
  dim(x) <- c(n, B)
  x
}

##########################################################

samp.combinations <-
function(n, B, k, both=T)
{
  # Return all combinations of k observations from 1:n,
  # one combination in each column.
  # This is useful for exhaustive calculations in two-sample permutation
  # applications.
  #
  # B must be equal to choose(n, k).
  # If both=F, then return an [k, B] matrix, where each sample
  # contains only the k indices selected.
  # If both=T, then return an [n, B] matrix, appending the n-k
  # indices omitted.
  if(B != choose(n, k))
    stop("B must be equal to choose(n, k) = ", choose(n, k))
  x <- combinations(n, k)
  if(both)
    rbind(x, apply(-x, 2, "[", x=1:n))  # (1:n)[-x[,j]]
  else
    x
}

##########################################################

samp.permutations <-
function(n, B)
{
  # Return all permutations of 1:n, one permutation in each column.
  # This is useful for exhaustive calculations in one-sample permutation
  # applications.
  #
  # B must be equal to factorial(n).
  if(B != factorial(n))
    stop("B must be equal to factorial(n) = ", factorial(n))
  permutations(n)
}


##########################################################

balancedSample <- 
function(n, size = n, prob = NULL, full.partition = "none")
{
  # Draw `size' observations from 1 to `n',
  # with probabilities proportional to `prob' (NULL means equal probs)
  # Let K = trunc(size * prob/sum(prob));
  # observation j is included either K[j] or K[j]+1 times
  if(lprob <- length(prob)){
    if(n != lprob)
      stop("n must match length(prob)")
    prob <- prob / sum(prob)
  }
  if(size >= n)
    full <- full.last <- F
  else{
    full.partition <- match.arg(full.partition, c("first", "last", "none"))
    if(full.partition == "first"){
      full <- T
      full.last <- F
    }
    else if(full.partition == "last"){
      full <- full.last <- T
    }
    else
      full <- full.last <- F
  }
  .C("S_samp_perm_probs",
     as.integer(size),
     as.integer(n),
     as.integer(length(prob)),
     as.double(prob),
     as.integer(full),
     as.integer(full.last),
     as.integer(1),  
     index = integer(size),
     COPY = c(F,F,F,F,F,F,F,T))$index
}
# Set default prob to NULL
# Do not define prob if not supplied (is not used in C code).
# Set default size to n.


##########################################################

samp.blockBootstrap <-
function(n, B, size = n - reduceSize, reduceSize = 0, blockLength)
{
  # Simple block bootstrap, overlapping blocks, no wrap-around,
  # no matching of ends.
  # To change blockLength use blockBootstrap()
  if(blockLength > n)
    stop("blockLength must be <= n")
  nblocks <- ceiling(size/blockLength)
  x <- sample(0:(n-blockLength), size = B*nblocks, replace = T)
  dim(x) <- c(nblocks, B)
  # x = random indices for block starts, -1
  apply(x, 2, function(y, L) 1:L + rep(y, each = L), 
	L = blockLength)[1:size,,drop = F]
}
# add size argument
# add drop = F to handle B = 1 case, add test for blockLength < n
# changed name from samp.boot.block to samp.blockBootstrap
# add argument reduceSize


##########################################################

blockBootstrap <- function(blockLength){
  # return a copy of samp.blockBootstrap with the desired blockLength
  f <- samp.blockBootstrap

  f$blockLength <- blockLength
  f
}
# bootstrap now has a sampler.args argument, but this is more convenient
# changed name from make.samp.boot.block to blockBootstrap

#| blockBootstrap(5)
#| samp.blockBootstrap(5, 10, blockLength = 3)
#| samp.blockBootstrap(6, 10, blockLength = 3)
#| samp.blockBootstrap(5, 10, 10, blockLength = 3)
#| samp.blockBootstrap(6, 10, 12, blockLength = 3)
#| bootstrap(1:25, mean)
#| bootstrap(1:25, mean, sampler=blockBootstrap(5))
#| bootstrap(1:25, mean, sampler=samp.blockBootstrap, sampler.args = list(blockLength = 5))
#|# the block version has a larger standard error, as it should


##########################################################

samp.bootknife <-
function(n, B, size = n - reduceSize, reduceSize = 0, njack = 1)
{
  # Bootknife Monte Carlo bootstrap.  Each column is one sample of 
  # size `size' chosen with replacement from a jackknife sample
  # (1:n with njack values omitted)
  #
  # If njack = 1, then omit each observation in B/n samples.
  # Otherwise omissions are not balanced.
  if(njack == 1){
    omit <- balancedSample(n, B)
    # Generate samples from 1:(n-1), then replace omitted with n
    x <- sample(n - 1, size = B * size, replace = T)
    x[x == rep(omit, each = size)] <- n
  }
  else{
    n2 <- n - njack
    keep <- samp.permute(n, B, size = n2)
    # Following is a vectorized way to do equivalent of:
    # inds <- samp.bootstrap(n2, B, size=size)
    # result[,j] = ( keep[,j] )[ inds[,j] ]
    x <- keep[ ceiling(n2 * runif(size*B))
	      + n2 * rep(0:(B-1), each = size) ]
  }
  dim(x) <- c(size, B)
  x
}
# Combined random and deterministic versions.
# Change meaning of random.  Now balanced to the extent possible,
#   and random refers to random reordering.
# Add size argument (7/10/01, RT)
# Use balancedSample; remove argument "random" (bad for stratified sampling)
# add argument reduceSize
# Add 'njack' argument, to allow omitting more than one value.


##########################################################

samp.finite <-
function(n, B, size = n - reduceSize, reduceSize = 0, N, bootknife=F)
{
  # Finite-population bootstrap sampling (without replacement)
  # data are a sample of size n from a population of size N
  # Create "superpopulation", repeating original data to make
  # sample of size N, then sample from that.
  # If N is not a multiple of size, some superpopulations will be larger
  # than N, others smaller.
  if(N < size) stop("N must be >= size")
  if(bootknife)
    m <- n-1
  else
    m <- n
  if(N %% m == 0){		# easy case, N is a multiple of m (n or n-1)
    k <- N/m			
    x1 <- rep(1:m, k)		# superpopulation
    result <- x1[samp.permute(N, B, size=size)]
    dim(result) <- c(size, B)
  }
  else {
    # Superpopulation sizes will vary
    k <- trunc(N/m)
    N1 <- k*m
    N2 <- N1+m
    lambda <- (1/N - 1/N2) / (1/N1 - 1/N2)
    x1 <- rep(1:m, k)		# smaller superpopulation
    x2 <- rep(1:m, k+1)		# larger  superpopulation
    B1 <- lambda * B
    B1 <- trunc(B1) + (runif(1) < (B1 %% 1)) # E[B1] = lambda*B
    B2 <- B - B1
    result <- c(x1[samp.permute(N1, B1, size=size)],
	       x2[samp.permute(N2, B2, size=size)])
    dim(result) <- c(size, B)
    # randomly permute columns, so small superpops aren't all first
    result <- result[, balancedSample(B)]
  }
  if(bootknife){
    omit <- balancedSample(n, B)
    result[result == rep(omit, each = size)] <- n
  }
  result
}
# add argument reduceSize


samp.half <-
function(n, B, size = n/2 - reduceSize, reduceSize = 0){
  # Return a matrix with ceiling(size) rows and B columns
  # In the usual case when size is not supplied,
  # each sample is half a random permutation
  # (and the following sample is the other half).
  #
  # If n is odd then half the samples are size floor(n/2), half ceiling(n/2).
  # The matrix has ceiling(n/2) rows, but alternate columns contain a 0.
  #
  # size and reduceSize may be integers or half-integers;
  # size is a half-integer, then alternate columns contain 0.
  #
  # If zeroes are needed, use them in the odd samples, to be conservative
  if(2*size > n)
    stop("size > n/2 is not supported; use samp.permute")
  x <- samp.permute(n, ceiling(B/2))
  if(2*size < n)
    x <- x[seq(length=2*size),,drop=F]
  if(size %% 1)
    x <- rbind(0, x)
  dim(x) <- dim(x)*c(1/2,2)
  if(B%%2)
    x[,1:B,drop=F]
  else
    x
}
