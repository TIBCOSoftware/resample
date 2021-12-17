# Copyright 2021. TIBCO Software Inc.
# This file is subject to the license terms contained
# in the license file that is distributed with this file.

findPvalues <- function(observed, replicates,
                       B, p, alternative, combine, combinationFunction){
  # Find p-values, for individual variables;
  # if length(combine), also find combined p-values.
  # In: observed[p], replicates[B,p], alternative[p] (char)
  #     combine[m] (list of m sets of indices to combine)
  # Out: Pvalues[p], combP[m]
  # This is a utility function called by permutation test functions,
  # and is not intended to be called directly by users.
  B1 <- sum(B) + 1
  if(!length(combine)){
    # No combination, individual p-values only
    Pvalues <- rep(NA, p)
    for(i in seq(length=p)){
      X <- replicates[,i]
      counts <- switch(alternative[i],
      		less    = sum(X <= observed[i]),
      		greater = sum(X >= observed[i]),
      		two.sided = 1 + 2 * min(sum(X <= observed[i]),
      		                        sum(X >= observed[i])))
      # Use (counts+1)/(B+1) rather than counts/B to be conservative,
      # and prevent zeros.  Add 1 in the two.sided case so
      # Pvalue = 2*min(one-sided Pvalues)
      Pvalues[i] <- min(1, (counts+1)/B1)
    }
    return(Pvalues)
  }
  # Find p-value(s) for observed & every replication, combine.
  rank <- function(x){
    # Resolve ties using the minimum (not mean), use na.last=F
    # (this is conservative, if the observed is one ofthe ties)
    ranks <- sort.list(sort.list(x, na.last = F))
    for(i in unique(x[duplicated(x)])) {
      which <- x == i & !is.na(x)
      ranks[which] <- min(ranks[which])
    }
    ranks
  }

  lambda <- matrix(as.double(NA), B1, p)
  for(i in seq(length=p)){
      X <- c(observed[i], replicates[,i])
      counts <- switch(alternative[i],
                      less = B1 - rank(-X),
                      greater = B1 - rank(X),
                      two.sided = 1 + 2 * (B1 - pmax(rank(X), rank(-X)))
                      )
      lambda[,i] <- pmin(1, (counts+1)/B1)
  }
  Pvalues <- lambda[1,]
  comblambda <- numeric(B1)
  combP <- sapply(combine,
      	     function(cols, lambda, f){
      	       comblambda <- f(lambda[ , cols, drop=F])
      	       sum(comblambda >= comblambda[1])
      	     },
      	     lambda=lambda, f=combinationFunction) / B1
  list(Pvalues = Pvalues,
       combP = combP)
}
# allow p=0
