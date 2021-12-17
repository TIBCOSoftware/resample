# Copyright 2021. TIBCO Software Inc.
# This file is subject to the license terms contained
# in the license file that is distributed with this file.

# Functions in this file:
#   tiltBootProbs
#   limits.tilt
#   print.limits.tilt
# should add p-value calculations

tiltBootProbs <-
function(tau, boot.obj, tilt = "exponential",
	 subset.statistic = 1,
	 observed = boot.obj$observed[subset.statistic],
	 replicates = boot.obj$replicates[,subset.statistic],
	 L = resampGetL(boot.obj,
	   frame.eval=frame.eval)[,subset.statistic],
	 indices = resampGetIndices(boot.obj,
	   frame.eval=frame.eval),
	 group = resampGetArgument(boot.obj, "group",
	   frame.eval=frame.eval),
	 treatment = resampGetArgument(boot.obj, "treatment",
	   frame.eval=frame.eval),
	 weights = boot.obj$weights,
	 normalizeWeights = T,
	 tol = 1e-6, maxiter = 20, debug=F,
	 frame.eval = boot.obj$parent.frame)
{
  if(is.null(frame.eval))
    frame.eval <- sys.parent(1)

  if(length(observed) > 1)
    stop("tiltBoot only supports univariate statistics")

  n <- length(L)
  B <- length(replicates)
  nw <- length(weights)
  if(nw && nw != B)
    stop("length of weights must match number of replications")
  if(nw && normalizeWeights)
    weights <- weights/mean(weights)

  # If both treatment and group, stratify on combination.
  if(length(treatment)){
    if(length(group))
      group <- paste(match(group, unique(group)),
		    match(treatment, unique(treatment)))
    else
      group <- treatment  # so can refer just to "group" below
  }

  if(byGroup <- length(group)){
    if(byGroup != n)
      stop("length of group must be same as length of L")
    group.inds <- split(1:n, group)
    group.names <- names(group.inds)
    nGroups <- length(group.inds)
    groupSizes <- unlist(lapply(group.inds, length))
    if(notSorted(ord <- unlist(group.inds))){
      # Reorder so that groups are contiguous, and modify indices
      L <- L[ord]
      # modify indices, so they point to current positions for L values
      ord2 <- 1:n
      ord2[ord] <- ord2
      olddim <- dim(indices)
      indices <- ord2[indices]
      dim(indices) <- olddim
    }
  }
  else{
    nGroups <- 1
    groupSizes <- n
  }


  # What tilting method; in "ml" case subtract group means
  tilt <- casefold(tilt)
  if(pmatch(tilt, "mle", nomatch = 0)){	# "ml", or "mle"
    p <- .C("S_tiltBoot_ml",
	    as.integer(n),
	    as.integer(B),
	    as.integer(length(tau)),
	    as.double(observed),
	    as.double(replicates),
	    as.double(subtractMeans(L, group = rep(1:nGroups, groupSizes))),
	    as.integer(indices),
	    as.integer(nw),
	    as.double(weights),
	    as.integer(nGroups),
	    as.integer(groupSizes),
	    as.double(tol), # for solving for lambda
	    as.integer(maxiter), # for solving for lambda
	    as.integer(debug),
	    p = as.double(tau),
	    COPY = c(rep(F, 5), L=T, rep(F,8), p=T))$p
  }
  else if(pmatch(tilt, "exponential", nomatch = 0))
    p <- .C("S_tiltBoot_exp",
	    as.integer(n),
	    as.integer(B),
	    as.integer(length(tau)),
	    as.double(observed),
	    as.double(replicates),
	    as.double(L),  # L is modified when there are multiple groups.
	    as.integer(indices),
	    as.integer(nw),
	    as.double(weights),
	    as.integer(nGroups),
	    as.integer(groupSizes),
	    as.integer(debug),
	    p = as.double(tau),
	    COPY = c(rep(F, 5), L=T, rep(F, 6), p=T))$p
  else
    stop("\n Supplied tilting method ", tilt, " is not recognized.")

  list(tau = tau, p=p)
}
# Add argument subset.statistic
# Add argument group, treatment, weights, debug
# Support bootstrap2 objects (implicitly, using treatment)

# Technical comments:
#   L is modified in C code, if multiple groups. (divide by n[g]/n)


##################################################

"limits.tilt"<-
function(boot.obj, probs = c(25, 50, 950, 975)/1000,
	 L = resampGetL(boot.obj, frame.eval=frame.eval),
	 t.adjustment = F,
	 subset.statistic = 1:p,
         subjectDivide = F, modifiedStatistic = NULL,
	 initial, tol = 1e-6, tol.tau = tol, maxiter = 20,
	 indices = resampGetIndices(boot.obj,
	   frame.eval = frame.eval),
	 group = resampGetArgument(boot.obj, "group",
	   frame.eval = frame.eval),
	 treatment = resampGetArgument(boot.obj, "treatment",
	   frame.eval = frame.eval),
	 bootWeights = boot.obj$weights,
	 frame.eval = boot.obj$parent.frame, debug=F, ...)
{
  # Find bootstrap tilting interval
  # This version doesn't update L; call it twice to do that.
  # Does not support missing values in replicates.
  if(is.null(frame.eval))
    frame.eval <- sys.parent(1)
  if(any(probs <= 0 | probs >= 1))
    stop("illegal value for probs. Must be in (0,1)")

  message <- resampCheckIfOrderMatters(boot.obj$order.matters)
  if(!is.null(message)){
    warning("Cannot calculate tilting limits--", message)
    p <- length(boot.obj$observed)
    return(matrix(NA, nrow=length(subset.statistic), ncol=length(probs)))
  }

  ###### Setup for evaluating the statistic #####
  bootcall <- boot.obj$call
  if(is(boot.obj, "concomitants"))
    bootcall <- boot.obj$original$call

  # Get statistic: if modifiedStatistic is supplied, use that.
  # Otherwise, get it from the bootstrap call.
  if(!missing(modifiedStatistic)){
    # This is either an expression like
    # mean(data, weights = Splus.resamp.weights)
    # in which case we use substitute(), or name of object containing
    # a stored expression.
    substitute.stat <- substitute(modifiedStatistic)
    if(mode(substitute.stat) != "name")
      modifiedStatistic <- substitute.stat
    if(!is.element("Splus.resamp.weights", all.names(modifiedStatistic)))
      stop('Your modifiedStatistic does not use a weights vector with name "Splus.resamp.weights"')
    statistic <- modifiedStatistic
  }
  else{
    substitute.stat <- bootcall$statistic

    # If statistic isn't a function or name of function or expression,
    # store it as an expression to pass to fit.func.
    if(!is.element(mode(substitute.stat), c("name", "function")))
      statistic <- substitute.stat
    else
      statistic <- eval(substitute.stat, frame.eval)

    if(is.function(statistic)){
      if(!any(is.element(c("weights","..."), names(statistic))))
        stop("the bootstrap statistic must have a weights argument")
    }
    else{
      # create a modifiedStatistic, add a weights argument
      modifiedStatistic <- addWeightsInCall(statistic)
      if(identical(statistic, modifiedStatistic))
	stop("unable to find a function in the bootstrap statistic with a `weights' argument")
      statistic <- modifiedStatistic
    }
  }

  # Get args.stat
  args.stat <- eval(bootcall$args.stat, frame.eval)
  # Get name of data.
  data.name <- bootcall$data
  data <- eval(data.name, frame.eval)
  if(!is.name(data.name))
    data.name <- "data"
  # Handle case of bootstrap2(data, data2)
  if(!is.null(bootcall$data2)){
    data2.name <- bootcall$data2
    data2 <- eval(data2.name, frame.eval)
    if(length(dim(data)))
      data <- rbind(data, data2)
    else
      data <- c(data, data2)
    data.name <- "data"
  }
  # Assign data to frame 1 if needed
  assign.frame1 <- bootcall$assign.frame1
  if(is.null(assign.frame1))
    assign.frame1 <- F
  if(assign.frame1){
    on.exit(if(exists(data.name, frame = 1)) remove(data.name, frame = 1))
    assign(data.name, data, frame=1)
  }
  # Get function to evaluate the statistic given data and weights
  is.df.data <- is.data.frame(data)
  fit.func <-
    resampMakeFuncWeighted(statistic = statistic,
                           data.name = data.name,
                           is.df.data = is.df.data,
                           df.var.names = if(is.df.data) names(data),
                           is.null.args.stat = is.null(args.stat),
                           assign.frame1 = assign.frame1)

  n <- boot.obj$n
  B <- boot.obj$B
  K <- length(probs)
  observed  <- boot.obj$observed
  bySubject <- length(subject <- resampGetArgument(boot.obj, "subject"))
  nObs <- if(bySubject) bySubject else boot.obj$n
  p <- length(observed)
  nw <- length(bootWeights)

  # Create a function to manipulate weights, if necessary, to handle
  # dividing or replicating subject weight to observations.
  # Call with first argument only; set defaults for others.
  if(bySubject){
    subjectIndices <- split(1:nObs, subject)
    subjectSizes <- sapply(subjectIndices, length)
    subjectNames <- names(subjectIndices)
    subjectI <- match(subject, subjectNames)
    if(subjectDivide){
      # by subject, and divide subject weight among observations
      fixWeights <- function(weight, Order, div){(weight/div)[Order]}
      fixWeights$Order <- subjectI
      fixWeights$div <- subjectSizes
    }
    else {
      # by subject, and replicate subject weight to observations
      fixWeights <- function(weight, Order) {weight[Order]}
      fixWeights$Order <- subjectI
    }
  }
  else # do nothing in this case
    fixWeights <- function(weight) weight

  # Use subset.statistic
  if(!missing(subset.statistic)){
    observed <- observed[subset.statistic]
    if(is.logical(subset.statistic))
      subset.statistic <- (1:p)[subset.statistic]
    p <- length(observed)
    L <- L[,subset.statistic,drop=F]
    replicates <- boot.obj$replicates[,subset.statistic,drop=F]
  }
  else
    replicates <- boot.obj$replicates

  if(!K || !p)
    stop("no confidence limits requested")
  dnames <- list(names(observed), paste(100 * probs, "%", sep = ""))

  if(t.adjustment)
    probs <- pnorm(qt(probs, df=n-1))

  # Check for missing values
  if(anyMissing(replicates))
    warning(sum(rowSums(is.na(replicates)) > 0),
	    " bootstrap samples have missing values; this may cause bias")

  # Check that the observed probabilities are between 0 and 1
  # and that L has nonzero variance.
  observedProbs <- colMeans(replicates >= rep(observed, each=B), na.rm=T,
			   weights=bootWeights)
  Lstdevs <- colStdevs(L, SumSquares=T, na.rm=T)
  omit <- (is.element(observedProbs, 0:1) | Lstdevs == 0)
  omit[is.na(omit)] <- T
  if(any(omit)){
    warning("Cannot compute confidence intervals for variable(s):\n",
	    paste(names(observed)[omit], collapse=", "),
	    "\nbecause L has zero variance or the observed value is not within the range of the bootstrap distribution")
    original.dnames <- dnames
    observed <- observed[!omit]
    p <- length(observed)
    if(!p)
      stop("No variables remaining to do confidence intervals for")
    subset.statistic <- subset.statistic[!omit]
    L <- L[ , !omit, drop=F]
    Lstdevs <- Lstdevs[ !omit ]
    observedProbs <- observedProbs[ !omit ]
    replicates <- replicates[ , !omit, drop=F]
    dnames[[1]] <- dnames[[1]][ !omit ]
  }

  # If both treatment and group, stratify on both;
  # if only treatment, then call it "group"
  if(byTreatment <- length(treatment)){
    if(byTreatment != n)
      stop("length of treatment must match number of observations")
    if(length(group))
      group <- paste(match(group, unique(group)),
		    match(treatment, unique(treatment)))
    else
      group <- treatment  # so can refer just to "group" below
  }

  if(byGroup <- length(group)){
    if(byGroup != n)
      stop("length of group must match number of observations")
    group.inds <- split(1:n, group)
    nGroups <- length(group.inds)
    groupSizes <- unlist(lapply(group.inds, length))
    if(notSorted(ord <- unlist(group.inds))){
      # Reorder so that groups are contiguous, and modify indices
      L <- L[ord,,drop=F]
      # modify indices, so they point to current positions for L values
      ord2 <- 1:n
      ord2[ord] <- ord2
      olddim <- dim(indices)
      indices <- ord2[indices]
      dim(indices) <- olddim
    }
  }
  else{
    nGroups <- 1
    groupSizes <- n
  }

  L <- subtractMeans(L, group = rep(1:nGroups, groupSizes), na.rm=T)

  if(missing(initial) || length(initial) != K*p)
    initial <- (qnorm(probs) -
	       rep(qnorm(observedProbs), each=K)) /
		   (rep(Lstdevs, each=K) * nGroups)

  out <- .C("S_tiltBootSolve",
	    as.integer(n),
	    as.integer(B),
	    as.integer(p),
	    as.integer(K),
	    as.double(observed),
	    as.double(replicates),
	    as.double(probs),
	    as.double(L),
	    as.integer(indices),
	    as.integer(nw),
	    as.double(bootWeights),
	    as.integer(nGroups),
	    as.integer(groupSizes),
	    as.double(tol.tau),
	    as.double(tol),
	    as.integer(debug),
	    as.integer(maxiter),
	    exp.tau = as.double(matrix(initial, p, K, byrow = T)),
	    ml.tau = double(p * K),
	    COPY = c(rep(F, 17),
	      exp.tau = T, ml.tau = T))[c("exp.tau", "ml.tau")]
  out <- list(exp.tau = matrix(out$exp.tau, p, K, dimnames=dnames),
	     ml.tau  = matrix(out$ml.tau,  p, K, dimnames=dnames))
  exp.limits <- ml.limits <- NA * out$exp.tau

  for(j in 1:p) {
    W1 <- tiltWeights(out$exp.tau[j,], L[, j], tilt="exponential",
		     group=rep(1:nGroups, groupSizes))
    W2 <- tiltWeights(out$ml.tau[j,],  L[, j], tilt="ml",
		     group=rep(1:nGroups, groupSizes))
    if(byGroup && notSorted(ord)){ # reorder W's to match order of data
      W1 <- W1[ord2,,drop=F]
      W2 <- W2[ord2,,drop=F]
    }
    if(!byTreatment){  # Usual case, not bootstrap2 application
      for(k in 1:K) {
	exp.limits[j, k] <-
	  fit.func(fixWeights(W1[,k]), data, statistic,
		   args.stat)[subset.statistic[j]]
	ml.limits[j, k] <-
	  fit.func(fixWeights(W2[,k]), data, statistic,
		   args.stat)[subset.statistic[j]]
      }
    }
    else { # bootstrap2 application, have treatment
      treat1 <- treatment == treatment[1]
      if(is.matrix(data)){
	data1 <- data[treat1,,drop=F]
	data2 <- data[!treat1,,drop=F]
      } else {
	data1 <- data[treat1]
	data2 <- data[!treat1]
      }
      assign("%-%",  # either - or /
	     if(!is.null(boot.obj$ratio) && boot.obj$ratio) get("/")
	     else get("-"))
      for(k in 1:K) {
	exp.limits[j, k] <-
	  fit.func(fixWeights(W1[treat1,k]), data1, statistic,
		   args.stat)[subset.statistic[j]] %-%
	  fit.func(fixWeights(W1[!treat1,k]), data2, statistic,
		   args.stat)[subset.statistic[j]]
	ml.limits[j, k] <-
	  fit.func(fixWeights(W2[treat1,k]), data1, statistic,
		   args.stat)[subset.statistic[j]] %-%
	  fit.func(fixWeights(W2[!treat1,k]), data2, statistic,
		   args.stat)[subset.statistic[j]]
      }
    }
  }
  # Create a matrix, with the most conservative values.
  #   Add the two types of limits, tau values, etc. as attributes.
  x <- list(exp.limits = matrix(exp.limits, p, K, dimnames = dnames),
	   ml.limits = matrix(ml.limits, p, K, dimnames = dnames),
	   exp.tau = out$exp.tau, ml.tau = out$ml.tau,
	   probs = probs)
  result <- x$ml.limits
  result[,probs >  .5] <- pmax(x$ml.limits[,probs > .5],
			       x$exp.limits[,probs > .5])
  result[,probs <= .5] <- pmin(x$ml.limits[,probs <= .5],
			       x$exp.limits[,probs <= .5])
  # Check for non-monotonicity (unless probs is not monotone)
  if(length(probs) > 1 && !notSorted(probs) &&
     any(result[,-1] < result[,-K])){
    warning("Confidence limits are not monotone; will modify.")
    lower <- which(probs <= .5)
    upper <- which(probs > .5)
    if(length(upper) && length(lower)){
      middle <- result[, max(lower) + 0:1, drop=F]
      result[,max(lower)] <- pmin(middle[,1], middle[,2])
      result[,min(upper)] <- pmax(middle[,1], middle[,2])
    }
    if(length(lower) > 1)
      for(k in rev(lower)[-1])
	result[,k] <- pmin(result[,k], result[,k+1])
    if(length(upper) > 1)
      for(k in upper[-1])
	result[,k] <- pmax(result[,k], result[,k-1])
  }

  # If omitted some columns, restructure output matrices now.
  if(any(omit)){
    restructure <- function(x, original.dnames, omit){
      X <- matrix(NA, length(omit), ncol(x), dimnames = original.dnames)
      X[!omit,] <- x
      X
    }
    result <- restructure(result, original.dnames, omit)
    for(i in intersect(names(x),
		       c("exp.limits", "ml.limits", "exp.tau", "ml.tau")))
      x[[i]] <- restructure(x[[i]], original.dnames, omit)
  }

  attributes(result) <- c(attributes(result), x)
  oldClass(result) <- "limits.tilt"
  result
}
# Add quit.if.nf argument, used by summary.bootstrap
# Add ... argument.  Though extra args ignored, need for summary.bootstrap.
# Handle replicates with missing values
# Add checks for columns of L with zero variance, or where the
# observed is not inside the range of the replicates.
# Don't check functionality for statistics for which intervals aren't done.
# Remove checks for functionality.  Is less important than I thought, at
#  least for the common cases (stdev, var) because both boot$observed
#  and boot$replicates are subject to same effects.
# Add components probs and call.
# Stop if sampling by treatment.
# Return an object containing the most conservative limits, with the
#   individual limits as attributes.
# Remove the call component (did not have a frame.eval component, and
#   seems not important to have this here).
# For monotonicity.
# Change order of arguments.
# Add args:  indices, group, treatment, tol.tau, maxiter, bootWeights, debug
# support group, treatment, and importance sampling weights
# Support bootstrap2(data, data2)
# Support bootstrap2(..., ratio)
# Better missing data handling - don't omit rows with NA in any column.
# return NA if resampling residuals


##################################################

print.limits.tilt <- function(x, ...){
  # Print bootstrap tilting limits.
  if(is.list(x)){  # backward compatibility, just print the list
    NextMethod("print")
    return(invisible(x))
  }
  # Print just the limits, not various attributes
  print.matrix(x, ...)
  invisible(x)
}
# Print only one set of limits (the most extreme of the two sets).
# Changed to assume that the object is a matrix with attributes.


##################################################

summary.limits.tilt <- function(object, ...){
  cat("Combined limits:\n")
  print.matrix(object, ...)
  ax <- attributes(object)
  cat("\nIndividual limits:\n")
  both <- rbind(exp = ax$exp.limits,
	       ml  = ax$ml.limits)
  dimnames(both)[[1]] <- paste(rep(c("exp"," ml"), each=nrow(ax$exp.limits)),
			      dimnames(both)[[1]])
  print(both, ...)
  cat("\nTilting parameters (tau):\n")
  both <- rbind(exp = ax$exp.tau,
	       ml  = ax$ml.tau)
  dimnames(both)[[1]] <- paste(rep(c("exp"," ml"), each=nrow(ax$exp.limits)),
			      dimnames(both)[[1]])
  print(both, ...)
  invisible(x)
}
# New function
