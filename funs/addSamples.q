# This file contains functions:
#   addSamples (generic function)
#   addSamples.bootstrap
#   addSamples.bootstrap2
#   addSamples.parametricBootstrap
#   addSamples.smoothedBootstrap
#   resampMakeWeightsRecip

# addSamples is used to add additional samples to a (p|s)bootstrap
# object.


##########################################################
addSamples <- function(object, B.add = 100, ...)
  UseMethod("addSamples")



##########################################################
# addSamples.bootstrap (was update.bootstrap)
##########################################################
addSamples.bootstrap <- function(object, B.add = 100, sampler.prob = NULL, 
				 ..., frame.eval = object$parent.frame)
{
  # Add B.add more resamples to an existing bootstrap object.
  if(is.null(frame.eval))
    frame.eval <- sys.parent(1)  

  # Create two calls; one to obtain new samples (evaluate this one),
  # the other to go in the output object.
  newCall <- object$call
  if(any(B.add < 100) && is.null(newCall$block.size))
    newCall$block.size <- min(100, eval(newCall$B, local = frame.eval))
  resultCall <- newCall
  newCall$B <- B.add
  newCall$seed <- object$seed.end
  newCall$sampler.prob <- sampler.prob
  newCall$L <- object$L
  
  save.indices <- !is.null(object$indices) || !is.null(object$compressedIndices)
  newCall$save.indices <- save.indices

  # Need to compute weights if there are weights in the original object
  # or if sampler.prob is specified.
  weights.old <- object$weights
  if(addWeights <- (!is.null(weights.old) || !is.null(sampler.prob)))
    newCall$save.indices <- T

  # Evaluate the new call
  New <- eval(newCall, local = frame.eval)

  if(addWeights){
    #
    # There will be sum(object$B) + sum(B.add) weights.  Instead of 
    # computing them from scratch, we take advantage of weights which 
    # have already been computed. 
    #
    # The first sum(object$B) weights are from the original ("old") 
    # resample indices, and include effects from the old sampler.probs
    # and new sampler.probs:
    #
    # 1/weights.old = k.old*(1/(weights.old.old) + k.new*(1/(weights.old.new)), 
    #
    # where k.old = sum(B.old)/(sum(B.old) + sum(B.add)), and 
    #       k.new = sum(B.add)/(sum(B.old) + sum(B.add))
    #
    # The final sum(B.add) weights are from the new indices, and also include 
    # effects from the old sampler.probs and new sampler.probs:
    #
    # 1/weights.new = k.old*(1/(weights.new.old) + k.new*(1/(weights.new.new)), 
    #
    # The final result is c(weights.old, weights.old).
    # weights.old.old and weights.new.new have already been computed. 


    # need stratification information for later
    if(!is.null(group <- resampGetArgument(object, "group"))){
      group.inds <- split(1:length(group), group)
      groupSizes <- sapply(group.inds, length)
    }
    else
      group.inds <- groupSizes <- NULL

    B.old <- resampGetArgument(object, "B")
    sumB.old <- sum(B.old)
    sumB.new <- sum(B.add)
    k.old <- sumB.old/(sumB.old + sumB.new)
    k.new <- sumB.new/(sumB.old + sumB.new)
    nB.old <- length(B.old)
    nB.new <- length(B.add)

    if(is.null(weights.old)){
      sampler.prob.old <- vector("list", nB.old)
      weightsRecip.old.old <- rep(1, sumB.old)
      weightsRecip.new.old <- rep(1, sumB.new)
    }
    else{
      sampler.prob.old <- resampGetArgument(object, "sampler.prob")
      if(!is.list(sampler.prob.old))
	sampler.prob.old <- rep(list(sampler.prob.old), nB.old)
      weightsRecip.old.old <- 1/weights.old
      weightsRecip.new.old <- 
	resampMakeWeightsRecip(B.old, sampler.prob.old, 
			       New$indices, group.inds, groupSizes)
    }
    if(is.null(sampler.prob)){
      sampler.prob.new <- vector("list", nB.new)
      weightsRecip.old.new <- rep(1, sumB.old) 
      weightsRecip.new.new <- rep(1, sumB.new) 
    }
    else{
      sampler.prob.new <- sampler.prob
      if(!is.list(sampler.prob.new))
	sampler.prob.new <- rep(list(sampler.prob.new), nB.new)
      weightsRecip.old.new <- 
	resampMakeWeightsRecip(B.add, sampler.prob.new, 
			       resampGetIndices(object), group.inds, groupSizes)
      weightsRecip.new.new <- 1/New$weights
    }

    weights <- c(1/(k.old*weightsRecip.old.old + k.new*weightsRecip.old.new), 
		 1/(k.old*weightsRecip.new.old + k.new*weightsRecip.new.new))
  }

  resultCall$B <- c(if(is.null(resultCall$B)) object$B
		   else eval(resultCall$B, local = frame.eval),
		   B.add)
  if(addWeights)
    resultCall$sampler.prob <- c(sampler.prob.old, sampler.prob.new)

  makeBootstrap(replicates = rbind(object$replicates, New$replicates),
	    observed = object$observed,
	    n = object$n,
	    call =  resultCall,
	    seed.start = object$seed.start,
	    seed.end = New$seed.end,
	    dim.obs = object$dim.obs,
	    group = object$group,
	    subject = object$subject,
	    parent.frame = frame.eval,
	    indices = if(!is.null(object$indices))
	      cbind(object$indices, New$indices),
	    compressedIndices = if(!is.null(object$compressedIndices))
	      cbind(object$compressedIndices, compressIndices(New$indices)),
	    weights = if(addWeights) weights, 
	    statistic.is.random = object$statistic.is.random,
	    label = object$label,
	    defaultLabel = object$defaultLabel,
            L = object$L,
            Lstar = if(length(object$Lstar) && length(New$Lstar))
	      rbind(object$Lstar, New$Lstar),
	    actual.calls = c(object$actual.calls,
	      list(object$call, match.call())))
}
# B can be vector (bug 13641), no update.B
# Correct call when B was not supplied originally #20239
# Use full components names, e.g. x$replicates instead of x$rep
# Change `x' to `object'
# Use named arguments.
# shorten comments at beginning
# Comment - B is total number of replications, call$B is a vector
# Changed name from update.bootstrap to addSamples.bootstrap (see addSamples.q)
#   (due to bug #14071, and to make update(bootstrap) more standard
# Restructured how calls are performed to avoid duplicate work and
#   make code more common across this and other methods
# Add parent.frame, remove group (may add group back in later) (like bootstrap)
# Added frame.eval argument (7/16/10, RT)
# Add actual.calls
# Support weights; (new argument sampler.prob is part of this)
# Add L from original object
# Support compressedIndices
# Support Lstar


# Comment - don't check for sys.parent() != object$parent.frame; would be rare.

# To do:
# * Need to support optional components to object, like weights, L


##########################################################

addSamples.bootstrap2 <- function(object, B.add = 100, ...)
  stop("bootstrap2 objects are not supported by addSamples")


##########################################################

addSamples.parametricBootstrap <- function(object, B.add = 100, ..., 
				  frame.eval = object$parent.frame)
{
  # Add B.add more resamples to an existing parametricBootstrap object.

  if(!is.null(object$weights))
    stop("weights detected, not yet supported")

  if(is.null(frame.eval))
    frame.eval <- sys.parent(1)  

  # Create two calls; one to obtain new samples (evaluate this one),
  # the other to go in the output object.
  newCall <- object$call
  if(any(B.add < 100) && is.null(newCall$block.size))
    newCall$block.size <- min(100, eval(newCall$B, local = frame.eval))
  resultCall <- newCall
  newCall$B <- B.add
  newCall$seed <- object$seed.end
  New <- eval(newCall, local = frame.eval)
  resultCall$B <- c(if(is.null(resultCall$B)) object$B
		   else eval(resultCall$B, local = frame.eval),
		   B.add)
  makeParametricBootstrap(replicates = rbind(object$replicates, New$replicates),
	     observed = object$observed,
	     n = object$n,
	     call = resultCall,
	     seed.start = object$seed.start,
	     seed.end = New$seed.end,
	     dim.obs = object$dim.obs,
	     parent.frame = frame.eval,
	     samples = c(object$samples, New$samples),
	     rsampler = object$rsampler,
	     args.rsampler = object$args.rsampler,
	     dsampler = object$dsampler,
	     defaultLabel = object$defaultLabel,
	     label = object$label,
	     actual.calls = c(object$actual.calls,
	       list(object$call, match.call())))

}
# Change log:
# Add actual.calls

##########################################################

addSamples.smoothedBootstrap <- function(object, B.add = 100, ...,
				  frame.eval = object$parent.frame)
{
  # Add B.add more resamples to an existing smoothedBootstrap object.

  if(!is.null(object$weights))
    stop("weights detected, not yet supported")

  if(is.null(frame.eval))
    frame.eval <- sys.parent(1)  

  # Create two calls; one to obtain new samples (evaluate this one),
  # the other to go in the output object.
  newCall <- object$call
  if(any(B.add < 100) && is.null(newCall$block.size))
    newCall$block.size <- min(100, eval(newCall$B, local = frame.eval))
  resultCall <- newCall
  newCall$B <- B.add
  newCall$seed <- object$seed.end
  New <- eval(newCall, local = frame.eval)
  resultCall$B <- c(if(is.null(resultCall$B)) object$B
		   else eval(resultCall$B, local = frame.eval),
		   B.add)
  result <- makeParametricBootstrap(replicates = rbind(object$replicates, New$replicates),
		       observed = object$observed,
		       n = object$n,
		       call = resultCall,
		       seed.start = object$seed.start,
		       seed.end = New$seed.end,
		       dim.obs = object$dim.obs,
		       parent.frame = frame.eval,
		       samples = c(object$samples, New$samples),
		       rsampler = object$rsampler,
		       args.rsampler = object$args.rsampler,
		       dsampler = object$dsampler,
		       defaultLabel = object$defaultLabel,
		       label = object$label,
		       actual.calls = c(object$actual.calls,
			 list(object$call, match.call())))
  result$rsampler <- result$args.rsampler <-
    result$dsampler <- NULL
  oldClass(result) <- c("smoothedBootstrap", "resamp")
  result
}
# Change log:
# Add actual.calls


##########################################################
# resampMakeWeightsRecip
##########################################################
resampMakeWeightsRecip <- function(B, sampler.prob, inds.mat, 
				   group.inds, groupSizes)
{
  # Compute the reciprocal of the weights for importance sampling. 
  # 
  # inds.mat may represent some fraction of the B resample indices, 
  # in which case return that fraction of the full reciprocal weights 
  # vector (length = numCols(inds.mat) := B.computed).
  #
  # sampler.prob is list[length(B)], elements are normed to sum to 1 (by group)
  # probRatio     is list[length(B)], normed to average to 1 (by group)
  # probRatioB    is [B.computed, length(B)] matrix of likelihood ratios

  dims <- dim(inds.mat)
  n <- dims[1]
  B.computed <- dims[2]
  withProb <- sapply(sampler.prob, length) > 0
  if(!is.null(group.inds)){# stratified data
    sampler.prob[withProb] <-
      lapply(sampler.prob[withProb],
	     function(x, groupl) {
	       for(i in groupl)
		 x[i] <- x[i]/sum(x[i])
	       x
	     },
	     groupl = group.inds)
    temp <- rep(groupSizes, groupSizes)
    temp[unlist(group.inds)] <- temp
    probRatio <- lapply(sampler.prob, get("*"), temp)
  }
  else {
    sampler.prob[withProb] <-
      lapply(sampler.prob[withProb], function(x) x/sum(x))
    probRatio <- lapply(sampler.prob, get("*"), n)
  }
  probRatio[!withProb] <- vector("list", sum(!withProb))
  probRatioB <- matrix(1, nrow =  B.computed, ncol = length(B))
  probRatioB[, withProb] <-
    sapply(probRatio[withProb],
	   indexProducts, indices = inds.mat)
  probRatioB %*% (B/sum(B))
}
