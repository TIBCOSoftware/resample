# Copyright 2021. TIBCO Software Inc.
# This file is subject to the license terms contained
# in the license file that is distributed with this file.

# Combination functions.
# Each must take a matrix as input, one column for every variable
# and any number of rows, containing p-values for individual variables.
# These functions must be monotone non-increasing in every column,
# and approach their maximum value (possibly Inf) as any column approaches 0.

combinePValues.Fisher <- 
  # Nonparametric combination of p-values, Fisher's form
  function(p, ...){
  -rowSums(log(p))
}

combinePValues.Liptak <- 
  # Nonparametric combination of p-values, Liptak's form
  function(p, ...){
  -rowSums(qnorm(p), n = nrow(p))
}

combinePValues.chisquare <- 
  # Nonparametric combination of p-values,  Chi-square form
  function(p, ...){
  rowSums(pmax(qchisq(1-p, df=1), 0), n = nrow(p))
}

combinePValues.halfnormal <- 
  # Nonparametric combination of p-values,  Half-normal form
  function(p, ...){
  # Assume that two-sided normal problem are converted to absolute values
  -rowSums(qnorm(p/2), n = nrow(p))
}

combinePValues.Tippett <- 
  # Nonparametric combination of p-values,  Tippett form
  function(p, ...){
    1 - apply(p, 1, min)
}


# Later generalize these, so that the transformation may depend on
# other rows.
# May need a `First' argument indicating that the first
# row should not be used in defining the transformation.
# May want an "alternative argument", then handle two-sided and one-sided
# separately in the normal case.
# ... allows for varying argument lists
