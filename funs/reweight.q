# In this file:
#   reweight  (for a single tau, or weights)
#   tiltAfterBootstrap
#   print.tiltAfterBootstrap
#   plot.tiltAfterBootstrap


# Define a reweighting function, bootstrap sensitivity
reweight <-
  function(boot.obj, data.weights,
	   probs = NULL, L = NULL, tau = NULL, tilt="exponential",
	   subset.statistic = 1,
           subjectDivide = F, modifiedStatistic = NULL,
	   frame.eval=boot.obj$parent.frame)
{
  # Reweight a bootstrap object, to estimate results that would
  # be obtained by sampling from a different empirical distribution
  # (with the same values, but different probabilities.
  #
  # Result is of class c("reweight", "resamp")
  # Input is a bootstrap
  # Output has the usual components of a resamp object
  #  (call, observed, replicates, estimate, B, n, dim.obs)
  # as well as new components
  #  (weights, data.weights)
  # and optional components
  #  (indices)
  # If object has weights components, assume it was created using
  # importance sampling, and the output weights are multiplied by those
  # weights.
  #
  # Note that observed is modified - it is recomputed using data.weights
  #
  # The reweighting is controlled by one of:
  #  (data.weights)
  #  (L, tau, tilt, subset.statistic) -- exponential or ML tilting method

  if(is.null(tilt))
    tilt <- "exponential"
  else
    tilt <- casefold(tilt)
  if(is.null(frame.eval))
    frame.eval <- sys.parent(1)
  group <- resampGetArgument(boot.obj, "group", frame.eval = frame.eval)
  treatment <- resampGetArgument(boot.obj, "treatment", frame.eval = frame.eval)
  args.stat <- eval(boot.obj$call$args.stat, frame.eval)

  # If both treatment and group, stratify on both;
  # if only treatment, then call it "group".
  if(length(treatment)){
    if(length(group))
      group <- paste(match(group, unique(group)),
		    match(treatment, unique(treatment)))
    else
      group <- treatment
  }

  # Get data.weights, and normalize
  if(missing(data.weights)){
    if(is.null(L))
      L <- resampGetL(boot.obj, frame.eval = frame.eval)[,subset.statistic]
    if(is.matrix(L) && ncol(L) > 1)
      stop("L must be a vector, not a matrix")
    if(is.null(tau)){
      if(is.null(probs))
        stop("one of probs and tau must be supplied if data.weights is missing.")
      tau <- saddlepointPSolve(probs, L, group=group)
      # add choice of "mle" in case user mistypes "ml" to be "mle":
      if(pmatch(tilt, "mle", nomatch = 0))
	tau <- tiltMeanSolve(tiltMean(tau, L, group=group)$q,
			    L, tilt = "ml", group=group)$tau.ml
    }
    if(length(tau) > 1){
      warning("reweight currently supports only scalar tau or probs")
      tau <- tau[1]
    }
    data.weights <- tiltWeights(tau=tau, L=L, tilt=tilt, group=group)
  }

  # Get subject argument, if any, from bootstrap object, for special
  # handling of weights.
  subject <- resampGetArgument(boot.obj, "subject", frame.eval = frame.eval)
  if(length(subject)){
    data.weights <-
      resampSubjectWeights(subject, weights = data.weights,
                           weights.out = "both", subjectDivide = subjectDivide)
    subjWeights <- data.weights$subject
    obsWeights <- data.weights$observation
    equal.weights <-
      resampSubjectWeights(subject, weights.out = "obs")
      
  }
  else{
    subjWeights <- obsWeights <- data.weights
    equal.weights <- rep(1, boot.obj$n)
  }
                                         
  # Get the weights
  weights <- indexProducts(if(is.null(group)) subjWeights/mean(subjWeights)
			      else subjWeights/groupMeans(subjWeights, group=group, repeats=T),
                              resampGetIndices(boot.obj))
  if(!is.null(boot.obj$weights))
    weights <- weights * boot.obj$weights

  ###### Get the reweighted observed value #######
  bootcall <- boot.obj$call

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
        stop("statistic does not have a weights argument")
    }
    else{
      # create a modifiedStatistic, add a weights argument
      modifiedStatistic <- addWeightsInCall(statistic)
      if(identical(statistic, modifiedStatistic))
	stop("unable to find a function in your statistic with a `weights' argument")
      statistic <- modifiedStatistic
    }
  }

  # Get name of data.
  data.name <- bootcall$data
  data <- eval(data.name, frame.eval)
  if(!is.name(data.name))
    data.name <- "data"

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
  fit.func <- resampMakeFuncWeighted(statistic = statistic,
				     data.name = data.name,
				     is.df.data = is.df.data,
				     df.var.names = if(is.df.data) names(data),
				     is.null.args.stat =
				     is.null(bootcall$args.stat),
				     assign.frame1 = assign.frame1)

  #  Get parameter values for observed data.
  if(length(treatment)){  # dealing with a bootstrap2 objects
    treat1 <- treatment == treatment[1]
    if(is.matrix(data)){
      data1 <- data[treat1,,drop=F]
      data2 <- data[!treat1,,drop=F]
    } else {
      data1 <- data[treat1]
      data2 <- data[!treat1]
    }
    observed <- (fit.func(obsWeights[treat1], data1, statistic, args.stat) -
		 fit.func(obsWeights[!treat1], data2, statistic, args.stat))
  }
  else
    observed <- fit.func(obsWeights, data, statistic, args.stat)
  if(!is.atomic(observed))
    stop("Statistic must return numerical (non-list) data.")

  # Compare results computed with equal weights to bootstrap observed value.
  if(length(treatment))
    boot.observed <- (fit.func(equal.weights[treat1], data1,
			      statistic, args.stat) -
		     fit.func(equal.weights[!treat1], data2,
			      statistic, args.stat))
  else
    boot.observed <- fit.func(equal.weights, data, statistic, args.stat)
  if(all.equal(as.numeric(boot.observed), as.numeric(boot.obj$observed)) != T)
    warning("observed value in bootstrap object is not the same as the statistic evaluated with equal weights.  Your statistic (or modifiedStatistic) may not be functional.")
  
  # Getting parameter names and coercing matrix to vector.
  names.observed <- resampMakeNames(observed, substitute.stat)
  dim.obs <- dim(observed)
  if(!is.null(dim.obs))
    observed <- as.vector(observed)
  names(observed) <- names.observed


  means <- colMeans(boot.obj$replicates, na.rm = T, weights=weights)

  result <-
    c(list(call=match.call(),
	   observed=observed),
      boot.obj[c("replicates", "B", "n", "dim.obs")],
      list(data.weights=subjWeights, weights=weights,
	   estimate=data.frame(Mean = means,
	     Bias = means - observed,
	     SE = colStdevs(boot.obj$replicates, na.rm = T, weights=weights))))
  if(!is.function(statistic))
    result$modifiedStatistic <- statistic
  oldClass(result) <- c("reweight", "resamp")
  result
}
# add argument df.var.names to call to resampMakeFuncWeighted
# Remove argument substitute.stat from resampMakeFuncWeighted (not used)
# Add code for handling groups
# treatment (handle bootstrap2 object)
# args.stat, subset.statistic
# Use indexProducts( log = F ), else incorrect on windows.  Remove log arg.



tiltAfterBootstrap <-
  function(boot.obj,
           functional = "Quantiles",
           passObserved = F,
           probs = ppoints(10, .5),
           L,
           tau,
           column = 1,
           tilt = "exponential",
           subjectDivide = F,
           modifiedStatistic = NULL,
           ...,
           frame.eval = boot.obj$parent.frame)
{
  # Bootstrap tilting diagnostics -- how does the value of
  # a functional (some summary of the bootstrap distribution,
  # likes quantiles, mean, bias, or standard error).
  # functional is either:
  #   function(replicates, weights, ...)
  #   function(replicates, observed, weights, ...)
  # depending on value of passObserved.
  # probs = approximate quantiles of the bootstrap distribution to target
  # L = linear approximation used for the tilting
  # tau = tilting parameter; supply this or probs
  # column = which column to use for tilting

  message <- resampCheckIfOrderMatters(boot.obj$order.matters)
  if(!is.null(message))
    stop("Cannot do tiltAfterBootstrap--", message)
  func.call <- match.call()
  # Ensure that func.call includes the functional, even if not passed
  if(is.null(func.call$functional))
     func.call$functional <- "Quantiles"

  tilt <- ifelse1(is.null(tilt), "exponential",
                  match.arg(casefold(tilt), c("exponential", "ml", "mle")))
  if(tilt == "mle") tilt <- "ml"
  if(is.null(frame.eval))
    frame.eval <- sys.parent(1)
  group <- resampGetArgument(boot.obj, "group", frame.eval = frame.eval)
  treatment <- resampGetArgument(boot.obj, "treatment", frame.eval =frame.eval)
  subject <- resampGetArgument(boot.obj, "subject", frame.eval = frame.eval)

  # If both treatment and group, stratify on both;
  # if only treatment, then call it "group".  If only group, don't change!
  if(length(treatment)){
    if(length(group))
      group <- paste(match(group, unique(group)),
		    match(treatment, unique(treatment)))
    else
      group <- treatment
  }

  if(missing(L))
    L <- resampGetL(boot.obj, frame.eval=frame.eval)
  if(is.matrix(L) && ncol(L) > 1)
    L <- L[,column]
  if(length(L) != boot.obj$n)
    stop("L is the wrong length")
  if(missing(tau)){
    tau <- saddlepointPSolve(probs, L, group=group)
    if(tilt == "ml")
      tau <- tiltMeanSolve(tiltMean(tau, L, group=group)$q,
			  L, tilt = "ml", group=group)$tau.ml
  }
  else
    probs <- saddlepointP(tau, L, group=group)

  indices <- resampGetIndices(boot.obj, frame.eval=frame.eval)
  boot.observed <-boot.obj$observed
  obs.names <- names(boot.observed)
  bySubject <- length(subject)
  nObs <- if(bySubject) bySubject else boot.obj$n
  p <- length(boot.observed)
  k <- length(tau)
  probs.names <-
    as.character(round(probs, max(3, abs(-logb(range(probs, na.rm=T),10)))))

  # Handle built-in functionals
  functionalChoices <- c("Mean", "Bias", "SE", "Bias&SE", "Mean&SE",
                         "Quantiles", "Centered Quantiles",
                         "Standardized Quantiles")
  quantileBased <- F
  if(is.character(functional)){
    functional <- match.arg(casefold(functional),
                            casefold(functionalChoices))
    passObserved <- is.element(functional,
                               c("bias", "bias&se", "centered quantiles",
                                 "standardized quantiles"))
    quantileBased <- is.element(functional,
                                c("quantiles", "centered quantiles",
                                  "standardized quantiles"))
  }
  orig.functional <- functional
  if(is.character(functional))
    functional <- resampFunctionalList[[functional]]

  # Functional for whole data
  func.full <-
    ifelse1(passObserved,
            functional(boot.obj$replicates,
                       observed = boot.observed,
                       weights = boot.obj$weights, ...),
            functional(boot.obj$replicates,
                       weights = boot.obj$weights, ...))
  dim.Func <- dim(func.full)
  dimnames.Func <- dimnames(func.full)

  # Get functional corresponding to "Mean", "Bias", "SE", or "Bias&SE".
  if(is.character(orig.functional)){
    func.names <- switch(orig.functional,
                         mean = paste("Mean", obs.names, sep = "."),
                         bias = paste("Bias", obs.names, sep = "."),
                         se   = paste("SE", obs.names, sep = "."),
                         "bias&se" = c(rbind(
                           paste("Bias", obs.names, sep = "."),
                           paste("SE", obs.names, sep = "."))),
                         "mean&se" = c(rbind(
                           paste("Mean", obs.names, sep = "."),
                           paste("SE", obs.names, sep = "."))),
                         quantiles =,
                         "centered quantiles" =,
                         "standardized quantiles" = {
                           probs.f <- dimnames.Func[[1]]
                           ifelse1(length(boot.observed) == 1,
                                   probs.f,
                                   t(outer(obs.names, probs.f, paste)))
                         })
  }
  else
    func.names <- resampMakeNames(observed = func.full,
                                  statistic = substitute(functional),
                                  prepend.stat.name = F,
                                  default.stat.name = "Func")
  if(!is.null(dim.Func))
    func.full <- as.vector(func.full)
  names(func.full) <- func.names


  ###### Setup for evaluating the statistic with weights
  bootcall <- boot.obj$call
  substitute.stat <- bootcall$statistic

  # Get statistic: if modifiedStatistic is supplied, use that.
  # Otherwise, get it from the bootstrap call.
  computeStat <- T
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
  else { # modified statistic is not supplied; modify it here
    substitute.stat <- bootcall$statistic

    # If statistic isn't a function or name of function or expression,
    # store it as an expression to pass to fit.func.
    if(!is.element(mode(substitute.stat), c("name", "function")))
      statistic <- substitute.stat
    else
      statistic <- eval(substitute.stat, frame.eval)

    if(is.function(statistic)){
      if(!any(is.element(c("weights","..."), names(statistic)))){
        if(passObserved)
          stop("statistic does not have a weights argument")
        computeStat <- F
      }
    }
    else { # statistic is a call, not a function
      # create a modifiedStatistic, add a weights argument
      modifiedStatistic <- addWeightsInCall(statistic)
      if(identical(statistic, modifiedStatistic)){
        if(passObserved)
          stop("unable to find a function in your statistic with a `weights' argument")
        computeStat <- F
      }
      else
        statistic <- modifiedStatistic
    } # end of case where statistic is a call
  } # end of case where modified statistic is not called

  if(computeStat){
    # Do this block if we will need to compute the statistic
    # (if the functional depends on the bootstrap distn and the statistic)
    # Get args.stat
    args.stat <- eval(bootcall$args.stat, frame.eval)
    # Get name of data.
    data.name <- bootcall$data
    data <- eval(data.name, frame.eval)
    if(!is.name(data.name))
      data.name <- "data"
    # Assign data to frame 1 if needed
    assign.frame1 <- bootcall$assign.frame1
    if(is.null(assign.frame1))
      assign.frame1 <- F
    if(assign.frame1){
      on.exit(if(exists(data.name, frame = 1)) remove(data.name, frame = 1))
      assign(data.name, data, frame=1)
    }

    # Create a function to manipulate weights, if necessary, to handle
    # dividing or replicating subject weight to observations.
    # Call with first argument only; set defaults for others.
    if(bySubject){  
      if(length(group))
	stop("Combination of subject and group or treatment is not supported in tiltAfterBootstrap")
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
    else # do nothing in this case, when not bySubject
      fixWeights <- function(weight) weight
  
    # Get function to evaluate the statistic given data and weights
    is.df.data <- is.data.frame(data)
    fit.func <-
      resampMakeFuncWeighted(statistic = statistic,
                             data.name = data.name,
                             is.df.data = is.df.data,
                             df.var.names = if(is.df.data) names(data),
                             is.null.args.stat = is.null(args.stat),
                             assign.frame1 = assign.frame1)

    # Test for functionality: is observed value same as result using
    # equal weights?


    if(length(treatment)){  # dealing with a bootstrap2 objects
      treat1 <- treatment == treatment[1]
      n1 <- sum(treat1)
      n2 <- sum(!treat1)
      if(is.matrix(data)){
	data1 <- data[treat1,,drop=F]
	data2 <- data[!treat1,,drop=F]
      } else {
	data1 <- data[treat1]
	data2 <- data[!treat1]
      }
      observedW <- (fit.func(rep(1/n1,n1), data1, statistic, args.stat) -
		    fit.func(rep(1/n2,n2), data2, statistic, args.stat))
    }
    else
      observedW <- fit.func(rep(1/nObs, nObs), data, statistic, args.stat)
    if(all.equal(as.numeric(boot.observed), as.vector(observedW)) != T)
      warning("observed value in bootstrap object is not the same as the statistic evaluated with equal weights.  Your statistic (or modifiedStatistic) may not be functional.")

    tilt.observed <- matrix(NA, k, p, dimnames=list(probs.names, obs.names))
  }  # end of: if(computeStat){...}

  # Effective sample sizes
  effectiveB <- numeric(k)
  effectiveBfun <- function(w){
    w <- w/sum(w)
    1/sum(w^2)
  }
    

  # Evaluate functional on the bootstrap distribution for each set of weights
  func.replicates <- matrix(NA, k, length(func.full))
  for(j in 1:k){
    data.weights <- tiltWeights(tau=tau[j], L=L, tilt=tilt, group=group)
    if(computeStat)
      tilt.observed[j,] <-
        ifelse1(length(treatment),
                (fit.func(data.weights[treat1], data1, statistic, args.stat) -
                 fit.func(data.weights[!treat1], data2, statistic, args.stat)),
                fit.func(fixWeights(data.weights), data, statistic, args.stat))
    normalized.weights <-
      ifelse1(length(group),
              data.weights / groupMeans(data.weights, group=group, repeats=T),
              data.weights * boot.obj$n)
    weights <- indexProducts(normalized.weights, indices)
    if(!is.null(boot.obj$weights))
      weights <- weights * boot.obj$weights
    func.replicates[j,] <- 
      ifelse1(passObserved,
              functional(boot.obj$replicates,
                         tilt.observed[j,], weights=weights,...),
              functional(boot.obj$replicates, weights=weights, ...))
    effectiveB[j] <- effectiveBfun(weights)
  }

  dimnames(func.replicates) <- list(probs.names, func.names)

  result <- list(call = func.call,
                 Func = func.full, # observed value of the functional
                 Func.replicates = func.replicates, # replicates of the functional
                 observed = boot.obj$observed, # observed statistic
                 statistic = if(computeStat) {
                   tilt.observed}, # replicates of the statistic
                 tau=tau, probs=probs,
                 column = column,
                 effectiveB = effectiveB,
                 quantiles = quantileBased,
                 dim.Func = dim.Func,
                 dimnames.Func = dimnames.Func
                 )
  oldClass(result) <- "tiltAfterBootstrap"
  result
}
# add argument df.var.names to call to resampMakeFuncWeighted
# Remove argument substitute.stat from resampMakeFuncWeighted (no longer there)
# Make functional argument function same as in jackknifeAfterBootstrap
#  (remove argument pivot, add argument passObserved)
# Make cleaner names, for when tau supplied instead of probs.
# Add a class, and a print method.
# Allow ml tilting if probs supplied instead of tau.
# treatment (handle bootstrap2 object)
# stop if resampling residuals
# Major redesign, to parallel changes in jackknifeAfterBootstrap
#   favor quantile-based approach, add various arguments & return components

print.tiltAfterBootstrap <-
  function(x, digits = max(options()$digits - 3, 4), ...)
{
  old.digits <- options(digits = digits)
  on.exit(options(old.digits))
  cat("Call:\n")
  print(x$call)
  cat("\nTilting based on:",
      ifelse1(is.character(x$column), x$column, names(x$observed)[x$column]))
  cat("\n", sep="",
      if(!is.null(x$statistic)) "Weighted statistic, ",
      "probs, functional of Bootstrap Distribution,\nand effective bootstrap sample size:\n")
  print(data.frame(x$statistic, probs = x$probs,
                   x$Func.replicates,
                   "B" = trunc(x$effectiveB),
		   row.names=NULL, check.names = F), ...)
}

plot.tiltAfterBootstrap <-
  function(x, plots = NULL, ..., omit = T, minimumB = 100, superpose = NULL)
{
  # Plot tiltAfterBootstrap results
  # plots should be one or more of:
  #   "fs" = functional vs statistic (default, if statistic present)
  #   "fp" = functional vs probabilities
  #   "sp" = statistic vs probabilities
  #   "pairs" = everything
  # Note that:
  #   statistic = weighted statistic, calculated for a tilted sample
  #   probs = approximately what quantile that weighted statistic is,
  #     in terms of the original bootstrap distribution.  More rigorously,
  #     is is the approximate quantile of a bootstrap distribution for
  #     a linear approximation to the statistic.

  # Convert x$column to numeric, then subscript various quantities
  k <- ifelse1(is.numeric(x$column), x$column,
               is.character(x$column), match(x$column, names(x$observed)),
               which(x$column))
  statistic <- x$statistic[,k,drop=F]
  subset.functional <- T # all columns, but may override this next
  xd <- x$dim.Func
  if(length(xd) == 2 && xd[[2]] > 1){
    # subset based on k & original array dimensions
    subset.functional <- array(1:prod(xd), dim=xd)[,k]
  }
  Func.replicates <- x$Func.replicates[,subset.functional,drop=F]
  if(is.null(plots))
    plots <- ifelse1(is.null(statistic), "pair", "fs")
  probs <- x$probs
  if(is.null(statistic)){
    warning("The statistic component is missing, some plots are not possible")
    plots <- setdiff(plots, c("fs", "fp"))
  }
  smallB <- (x$effectiveB < minimumB)
  ok <- 1
  if(any(smallB)){
    warning(sum(smallB), " points have effective bootstrap sample size < ",
            minimumB)
    if(omit)
      ok <- c(1,NA)[smallB+1]
  }
  statName <- paste("Weighted", dimnames(statistic)[[2]])

  if(is.element(plots, c("fs", "fp"))){
    if(is.null(superpose))
      superpose <- any(x$quantiles)
    ylab <- names(x$Func)
    if(superpose){
      if(any(x$quantiles)){
        percs <- x$dimnames.Func[[1]]
        percs <- as.numeric(substring(percs, 1, nchar(percs)-1))
        ylab <- paste(paste(round(percs,2), collapse=", "), "% ", sep="",
                      x$call$functional,
                      " for ", x$dimnames.Func[[2]][k])
      } else { # less common case - superpose other statistics
        ylab <- paste(as.character(deparse(x$call$functional)),
                      "for", x$dimnames.Func[[2]][k])
      }
    }
    if(!superpose) ylab <- rep(ylab, length=ncol(Func.replicates))
    # In the following plots, if ... has xlim, ylim, xlab, ylab
    # the last copy (the one in ...) is used
    if(is.element("fs", plots)){
      if(superpose)
        matplot(statistic*ok, Func.replicates*ok, type = "l",
                xlim = range(statistic), ylim=range(Func.replicates),
                xlab = statName, ylab=ylab, lty = 1, ...)
      else
        for(j in 1:ncol(Func.replicates))
          plot(statistic*ok, Func.replicates[,j]*ok, type = "l",
                xlim = range(statistic), ylim=range(Func.replicates[,j]),
                xlab = statName, ylab=ylab[j], ...)
    }
    if(is.element("fp", plots)){
      if(superpose)
        matplot(probs*ok, Func.replicates*ok, type = "l",
                xlim = range(probs), ylim=range(Func.replicates),
                xlab = "probs", ylab=ylab, lty = 1, ...)
      else
        for(j in 1:ncol(Func.replicates))
          matplot(probs*ok, Func.replicates[,j]*ok, type = "l",
                  xlim = range(probs), ylim=range(Func.replicates[,j]),
                  xlab = "probs", ylab=ylab[j], lty = 1, ...)
    }
  }
  if(is.element("sp", plots)){ # statistic vs probs
    plot(statistic*ok, probs*ok, type="l",
         xlim = range(statistic), ylim = range(probs),
         ylab = statName, xlab = "probs", ...)
  }
  if(is.element("pairs", plots)){
    X <- cbind(probs = probs,
                x$Func.replicates[,subset.functional,drop=F],
                x$statistic)
    pairs(X * ok, ...)
  }
  invisible(NULL)
}
# Major redesign, to parallel changes in jackknifeAfterBootstrap
#   favor quantile-based approach, add various arguments & return components
