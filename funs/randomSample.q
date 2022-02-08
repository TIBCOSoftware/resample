# Copyright 2021. TIBCO Software Inc.
# This file is subject to the license terms contained
# in the license file that is distributed with this file.

"randomSample" <- 
function(x, size=numRows(x), prob=NULL, full.partition = "none")
{
  ## Draw one sample with replacement from a vector, matrix, or data frame.
  inds <- balancedSample(numRows(x), size=size, prob=prob,
                         full.partition=full.partition)
  if(is.null(dim(x)))
    x[inds, drop=F]
  else if(is.data.frame(x))
    resampSubDF(x, inds)
  else
    x[inds,, drop=F]
}
