# Copyright 2021. TIBCO Software Inc.
# This file is subject to the license terms contained
# in the license file that is distributed with this file.

# file contains one function:  (for S version 4)
# makeParametricBootstrap: auxiliary function for parametricBootstrap()
#           : calculates parametricBootstrap statistic, creates parametricBootstrap object.

makeParametricBootstrap <- function(replicates, observed, n, call, seed.start, seed.end, 
		       dim.obs = NULL, samples = NULL, rsampler, args.rsampler,
		       dsampler, ...){

  # Calculate parametricBootstrap statistics and create a "parametricBootstrap" object.

  if(anyMissing(replicates)) {
    warning("NA's encountered in replicates.  Mean and standard error are calculated using remaining observations")
    B.missing <- sum(rowSums(is.na(replicates)) > 0)
    boot.mean <- colMeans(replicates, na.rm = T)
    boot.se <- colStdevs(replicates, na.rm = T)
  }
  else {
    B.missing <- NULL
    boot.mean <- colMeans(replicates)
    boot.se <- colStdevs(replicates)
  }

  boot.bias <- boot.mean - observed
  est <- data.frame(Mean = boot.mean, Bias = boot.bias, SE = boot.se)
  B <- dim(replicates)[1]
  result <- list(call = call, 
		 observed = observed, 
		 replicates = replicates, 
		 estimate = est, 
		 B = B, 
		 n = n, 
		 dim.obs = dim.obs, 
		 seed.start = seed.start, 
		 seed.end = seed.end, 
		 rsampler = rsampler, 
		 args.rsampler = args.rsampler, 
		 dsampler = dsampler, 
		 B.missing = B.missing, 
		 samples = samples,
		 ...)
  # Components after rsampler may be NULL & hence omitted

  # Strip certain components if they are NULL
  for(i in c("args.rsampler", "dsampler", "B.missing", "samples",
	     "label")){
    if(is.null(result[[i]]))
      result[[i]] <- NULL
  }
  oldClass(result) <- c("parametricBootstrap", "resamp")
  result
}

# Add ... argument to accomodate any new components (for now, parent.frame)
# (RT, 7/16/01)
# Switch order of Mean & Bias for consistency with other functions.
