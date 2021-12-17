# Copyright 2021. TIBCO Software Inc.
# This file is subject to the license terms contained
# in the license file that is distributed with this file.

resampSubDF <- function(x, i){
  # subscript rows of a data frame.
  # Intended for use in bootstrapping or other resampling.
  # i should be numeric or logical, not character
  ax <- attributes(x)
  f <- function(xj, i) if(length(dim(xj)) == 2) xj[i,,drop=F] else xj[i]
  x <- lapply(x, f, i=i)
  ax$row.names <- ax$row.names[i]
  attributes(x) <- ax
  x
}
