# bootstrap and other resampling functions

# Functions in this file:
#-----------------
# bootstrap
# bootstrap.default
# jackknife
# jackknife.default
# resampMakeFunc	(was resamp.get.fit.func)
# resampMakeFuncSimple	(was resamp.get.pfit.func) (new)
# resampMakeNames	(was resamp.get.dimnames)
# resampMakeArgument
#-----------------
# makeBootstrap
# jackstats
# resampGetIndices	(was resamp.get.indices)
# resampGetArgument	(new)
# limits.percentile
# limits.emp
# limits.bca
# plot.resamp
# print.resamp
# print.summary.bootstrap
# print.summary.resamp
# qqnorm.resamp
# summary.bootstrap
# summary.resamp
# pairs.resamp
# boxplot.resamp
# scatterPlot
# scatterPlot.default
# scatterPlot.resamp
# resampCheckIfOrderMatters

# Now in other files:
# jackknifeAfterBootstrap (print, plot)

# Now removed
# resamp.return.seed	(new, but in S+6.0)  (removed for Version 2.0)
# resamp.set.seed	(new)                (removed for Version 2.0)



# There are no functions defined in the main code tree bootstrap directory
# that are not represented here (possibly with different names)
# with the exception of
#	samp.jack		is not in the product.
#	addSamples.bootstrap	(was update.bootstrap, now in another file)


##########################################################
# bootstrap()
# now a generic function
##########################################################

# Generic function
# Methods are defined in file bootstrapMethods.q

bootstrap <-
function(data, statistic, B = 1000, args.stat = NULL,
	 group = NULL,
	 subject = NULL,
	 sampler = samp.bootstrap,
	 seed = .Random.seed,
	 sampler.prob = NULL,
	 sampler.args = NULL,
	 sampler.args.group = NULL,
	 resampleColumns = NULL,
	 label = NULL,
	 statisticNames = NULL,
	 block.size = min(100, B),
	 trace = resampleOptions()$trace,
	 assign.frame1 = F,
	 save.indices = F,
	 save.group = NULL,
	 save.subject = NULL,
	 statistic.is.random = NULL,
	 group.order.matters = T,
         order.matters = NULL,
	 seed.statistic = 500,
	 L = NULL,
	 model.mat = NULL,
	 argumentList = NULL,
         observed.indices = 1:n,
	 ...)
{
  if(missing(statistic))
    stop("Argument statistic is missing with no default")
  UseMethod("bootstrap")
}
# make argument list same as bootstrap.default (so gui code works)
# add argument resampleColumns

bootstrap.default <-
function(data, statistic, B = 1000, args.stat = NULL,
	 group = NULL,
	 subject = NULL,
	 sampler = samp.bootstrap,
	 seed = .Random.seed,
	 sampler.prob = NULL,
	 sampler.args = NULL,
	 sampler.args.group = NULL,
	 resampleColumns = NULL,
	 label = NULL,
	 statisticNames = NULL,
	 block.size = min(100, B),
	 trace = resampleOptions()$trace,
	 assign.frame1 = F,
	 save.indices = NULL,
	 save.group = NULL,
	 save.subject = NULL,
	 statistic.is.random = NULL,
	 group.order.matters = T,
         order.matters = NULL,
	 seed.statistic = 500,
	 L = NULL,
	 model.mat = NULL,
	 argumentList = NULL,
         observed.indices = 1:n,
	 ...)
{
  # Basic bootstrap function.
  #
  # Organizaton of this function:
  # Initialization:
  #   call, handle data and statistic, subject and group arguments
  #   create function to evaluate statistic
  #   sampler (can be a list of arguments to bootstrap)
  #   sampling by group (multiple samples, or stratified) and subject
  #   sampler probabilities (importance sampling)
  #   create call to sampler, & additional call in group case
  #   on.exit() call, in case interrupted during loops
  # Loop over values of B (if B is a vector)
  #   Loop over blocks (block.size replications at a time)
  #     Create matrix of indices  (may require looping over groups)
  #     Loop over individual observations
  #       calculate statistic
  #     calculate probRatioB (if importance sampling)
  #   end block loop
  # end B loop
  # Call makeBootstrap()

  if(missing(statistic))
    stop("Argument statistic is missing with no default")

  # Capture call
  func.call <- unpacked.call <- match.call()
  parent.frame <- sys.parent()

  # Check if `data' is the output of a modeling function; if so
  # take special action.  For some modeling functions (e.g. lm)
  # there are bootstrap methods that are used instead of this.
  if((is.call(fitCall <- data$call) ||
      is.call(fitCall <- attr(data,"call"))) &&
     !is.null(fitCall$data) && class(data) != "model.list")
    return(resampMakeFitObjResults(resampCall=func.call, fitCall=fitCall,
				   frame.eval=parent.frame, model.method=NULL))

  # If argumentList is supplied, unpack it
  # (supports arguments except: data, statistic, group, subject)
  for(i in names(argumentList)){
    assign(i, argumentList[[i]])
    unpacked.call[[i]] <- argumentList[[i]] # needed for resampGetL
  }

  # If statistic isn't a function or name of function or expression,
  # store it as an unevaluated expression to pass to fit.func
  substitute.stat <- substitute(statistic)
  if(!is.element(mode(substitute.stat), c("name", "function")))
    statistic <- substitute.stat

  # Get name of data
  data.name <- substitute(data)
  if(!is.name(data.name))
    data.name <- "data"
  is.df.data <- is.data.frame(data)

  # How many observations, or subjects if sampling by subject?
  n <- numRows(data)  # This will be changed if by subject

  # Save group or subject arguments?
  if(is.null(save.subject)) save.subject <- (n <= 10000)
  if(is.null(save.group))   save.group   <- (n <= 10000)

  # If data is a data frame, look for subject there first.
  subjectName <- func.call$subject
  if(!missing(subject) && is.df.data){
    subject <- resampMakeArgument(subjectName, data, parent.frame)
    subjectInDF <- is.name(subjectName) && is.element(subjectName, names(data))
  }
  else
    subjectInDF <- F
  bySubject <- length(subject)
  if(bySubject){
    if(bySubject != n)
      stop("length of subject must match number of observations in data")
    subjectIndices <- split(1:n, subject)
    subjectSizes <- sapply(subjectIndices, length)
    n <- length(subjectIndices)
  }
  else
    subjectIndices <- subjectSizes <- NULL

  # Create function to evaluate the statistic given data, indices, ...
  if(is.null(assign.frame1)) assign.frame1 <- F
  fit.func <- resampMakeFunc(statistic = statistic, data.name = data.name,
			     is.df.data = is.df.data,
			     is.null.args.stat = is.null(args.stat),
			     assign.frame1 = assign.frame1,
			     dimLength = length(dim(data)),
			     bySubject = bySubject,
                             subjectInDF = subjectInDF,
			     resampleColumns = resampleColumns)
  # is: function(inds, data, statistic, args.stat, subjectIndices,
  #              subjectSizes, nSubjects, subjectName)
  if(length(list(...)))
    warning("unrecognized arguments in function call")


  # Set seed in case statistic uses randomization
  seed.start <- seed
  if(is.null(statistic.is.random) || statistic.is.random) {
    set.seed(seed.statistic)
    prev.seed <- .Random.seed
  }

  # Calculate statistic for observed data
  if(assign.frame1)
    on.exit(if(exists(data.name, frame = 1)) remove(as.character(data.name),
				   frame = 1))
  observed <- fit.func(inds = observed.indices, data = data,
                       statistic = statistic, args.stat = args.stat,
                       subjectIndices = subjectIndices,
                       subjectSizes = subjectSizes, nSubjects = n,
                       subjectName = subjectName)
  if(length(observed) && !is.atomic(observed))
    stop("Statistic must return numerical (non-list) data.")
  if(mode(observed) == "missing"){
    observed <- NULL
    functionBody(fit.func) <- c(functionBody(fit.func), Quote({NULL}))
  }

  # Determine if statistic uses randomization; this may fail if
  #  a statistic sometimes use randomization.
  if(is.null(statistic.is.random))
    statistic.is.random <- any(.Random.seed != prev.seed)
  if(statistic.is.random)
    seed.statistic <- .Random.seed

  # Get statistic names and coerce to vector
  names.observed <- resampMakeNames(observed = observed,
				    statistic = substitute.stat,
				    statisticNames = statisticNames)
  dim.obs <- dim(observed)
  observed <- as.vector(observed)
  names(observed) <- names.observed

  set.seed(seed.start)

  # Modify group.order.matters, if necessary.  Related error checking.
  modifySize <- (any(is.element(c("size", "reduceSize"), names(args.stat))) ||
		any(is.element(c("size", "reduceSize"),
			       names(unlist(sampler.args.group, recursive=F)))))
  if(!is.null(resampleColumns)){
    if(modifySize)
      stop("Cannot resample with modified size when sampling only some columns")
    if(!is.matrix(data))
      stop("data must be a matrix or data frame to sample only some columns")
    if(!group.order.matters)
      stop("Need group.order.matters = TRUE to sample only some columns")
  }
  if(modifySize)
    group.order.matters <- F
  
		  

  # Sampling by group
  if(!missing(group) && is.df.data)
    group <- resampMakeArgument(func.call$group, data, parent.frame)
  byGroup <- length(group)
  if(byGroup){
    # Get indices
    if(byGroup != numRows(data))
      stop("length of group != number of observations in the data")
    if(bySubject){
      # get group for each subject (subject must be nested within group)
      subjectByGroup <- lapply(split(subject, group), unique)
      if(notNested(subjectByGroup))
	stop("some values of `subject' occur in multiple groups; subject must be nested within the `group' argument")
      if(group.order.matters){
	group.order.matters <- F
	if(!missing(group.order.matters))
	  warning("setting group.order.matters to F for sampling by subject and group")
      }
      group.inds <- lapply(subjectByGroup, match, table=names(subjectIndices))
      # contains indices into elements of subjectIndices
    }
    else
      group.inds <- split(1:n, group)
    groupNames <- names(group.inds)
    nGroups <- length(group.inds)
    groupSizes <- sapply(group.inds, length)
  }


  # Is the statistic a simple mean, or colMeans?
  # If so use indexMeans later, and convert a data frame to matrix now.
  # Also calculate L now.
  statistic.is.mean <- (is.name(substitute.stat) &&
		       is.element(substitute.stat,
				  c("mean", "colMeans")) &&
		       is.null(args.stat) &&
		       !bySubject &&
		       length(dim(data)) <= 2)
  if(is.null(trace))
    trace <- (!statistic.is.mean)
  if(statistic.is.mean && is.matrix(data)){
    if(ncol(data) > 1 && substitute.stat == "mean"){
      warning("You are computing the mean for a matrix with multiple columns; colMeans may be more appropriate.")
      statistic.is.mean <- F
    }
    else if(is.df.data)
      data <- numerical.matrix(data)
  }
  if(statistic.is.mean &&
     (is.null(func.call$L) || identical(func.call$L, "choose"))){
    # add L, if original call was NULL or "choose"
    L <- data
    if(byGroup) { # overall mean is a weighted sum of group means
      L <- subtractMeans(L, group)
      if(!is.matrix(L)) L <- as.matrix(L)
      for(i in 1:nGroups)
	L[group.inds[[i]], ] <- L[group.inds[[i]], ] * (groupSizes[i]/n)
    }
  }


  # If L is numerical, "jackknife", "influence", evaluate L now, in case
  # used by a sampler.
  # If "regression" or "ace", compute later.
  # If "choose", make same decision as in resampGetL (code is duplicated below)
  message <- resampCheckIfOrderMatters(order.matters)
  if(!is.null(message) && !is.null(L)){
    warning("Cannot calculate L--", message)
    L <- NULL
  }
  regressionL <- F
  if(is.numeric(L) || is.logical(L)){
    if(!is.matrix(L))
      L <- as.matrix(L)
    if(is.data.frame(L))
      L <- numerical.matrix(L)
    if(numRows(L) != n){
      warning("L has the wrong number of rows; will delete L")
      L <- NULL
    }
    if(bySubject){
      # In this case, require L to be in order of sorted groups,
      # and assume it is already normalized to mean 0 by group.
      if(length(dimnames(L)[[1]]) &&
	 !identical(dimnames(L)[[1]], names(subjectIndices)))
	stop("rows of L were not ordered by subject")
    }
    else {
      # Normalize columns of L to mean 0 (by group, if group present)
      L <- subtractMeans(L, group)
      dimnames(L) <- list(NULL, names.observed)
    }
  }
  else if(is.character(L)){
    modifiedStatistic <- NULL
    if(L == "choose"){
      if(!missing(model.mat))
        L <- "ace"
      else if(sum(B) > 2*n+100)
        L <- "ace"
      else if(missing(group))
        L <- "jackknife"  # no group
      else{# group
        if(is.function(statistic)){
          if(any(is.element(c("weights","..."), names(statistic))))
            L <- "influence"
          else
            L <- "jackknife"
        }
        else {
          # create a modifiedStatistic, add a weights argument
          modifiedStatistic <- addWeightsInCall(statistic)
          if(identical(statistic, modifiedStatistic))
            L <- "jackknife"
          else
            L <- "influence"
        }
      }
    }
    if(L == "jackknife"){
      if(trace)
	cat("Creating L using jackknife method\n")
      L <- resampGetL(unpacked.call, frame.eval = parent.frame,
		     method = L)
    }
    else if(L == "influence"){
      if(trace)
	cat("Creating L using influence method\n")
      L <- resampGetL(unpacked.call, frame.eval = parent.frame,
		     method = L, modifiedStatistic = modifiedStatistic)
    }
    else if(is.element(L, c("regression", "ace")))
      regressionL <- T
    else if(L == "none")
      L <- NULL
    else
      stop("unrecognized value for L")
  }
  else if(!is.null(L))
    stop("L is not numeric or character")



  # Check sampler probabilities
  nB <- length(B)
  sumB <- sum(B)
  if(is.null(sampler.prob)){
    anyProb <- F
    withProb <- rep(F, nB)
  }
  else {
    checkProb <- function(prob, n){
      m <- length(prob)
      if(m == 0)
	return(NULL)
      if(m != n)
	stop("length of sampling probabilities must match number of observations")
      if(any(prob < 0))
	stop("sampling probabilities may not be negative")
      NULL
    }
    if(is.list(sampler.prob)){
      if(length(sampler.prob) != nB)
	stop("sampler.prob must have length(B) components")
      lapply(sampler.prob, checkProb, n=n)
      withProb <- sapply(sampler.prob, length) > 0
    }
    else {
      checkProb(sampler.prob, n)
      withProb <- rep(T, nB)
      sampler.prob <- rep(list(sampler.prob), nB)
    }
    anyProb <- any(withProb)
    if(anyProb && !is.element("prob", names(sampler)))
      stop("sampling probabilities were provided but sampler does not accept probabilities")
  }

  # normalize probabilities (by group)
  # sampler.prob is list[nB], elements are normed to sum to 1 (by group)
  # probRatio     is list[nB], normed to average to 1 (by group)
  # probRatioB    is [sumB, nB] matrix of likelihood ratios
  if(anyProb){
    if(byGroup){
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
    probRatioB <- matrix(1, sumB, nB)
  }
  else
    probRatioB <- NULL

  # The sampler will be called in one of four ways, depending
  # on whether there are groups and probabilities:
  #    sampler(n, B2)
  #    sampler(n, B2, prob=probIB)
  #    sampler(groupSizes[iGroup], B2)
  #    sampler(groupSizes[iGroup], B2, prob=probIB[thisGroup])
  # Also, if sampler.args is present or e.g. sampler=samp.finite(N=50),
  # add those arguments to the call
  # Create call now, for use in the loop.
  # In addition, in the group case:
  #    use sampler.args.group[[iGroup]]
  #    create a second call to fill inds.mat, either
  #        inds.mat[thisGroup,] <-     (sampled indices)
  #        inds.mat <- rbind(inds.mat, (sampled indices))
  if(is.call(samplerCall <- substitute(sampler)) &&
     {
       newSampler <- get(samplerCall[[1]])
       all(names(newSampler)[1:2] == c("n","B"))
     }){
    # Have something like sampler = samp.finite(N=50)
    samplerCall[[1]] <- as.name("list")
    sampler.args <- c(sampler.args, eval(samplerCall, sys.parent()))
    sampler <- newSampler
  }
  samplerCall <- Quote(sampler(n=n, B=B2))
  if(length(sampler.args))
    samplerCall <- as.call(c(samplerCall, sampler.args))
  if(anyProb)
    samplerCall$prob <- Quote(probIB)
  # previous assignment will be overwritten if byGroup
  if(byGroup){
    samplerCall$n <- Quote(groupSizes[iGroup])
    if(anyProb)
      samplerCall$prob <- Quote(probIB[thisGroup])
    if(saglen <- length(sampler.args.group)){
      if(saglen != nGroups)
	stop("argument sampler.args.group should be of length equal to ",
	     "number of groups")
      if(length(sagNames <- names(sampler.args.group))){
	if(any(sort(sagNames) != sort(groupNames)))
	  stop("names of sampler.args.group must match group values")
        # rearrange to be in same order as sorted groups
	sampler.args.group <- sampler.args.group[groupNames]
      }
    }
    # Now construct an expression to get the sample indices for each stratum,
    # to be evaluated in the loop over groups
    samplerCall.no.sag <- samplerCall
    getGroupInds <-
      c(Quote({thisGroup <- group.inds[[iGroup]]}),
	if(saglen)
	Quote({samplerCall <- as.call(c(samplerCall.no.sag,
					sampler.args.group[[iGroup]]))}),
	if(group.order.matters)
	Quote({inds.mat[thisGroup, ] <- thisGroup[eval(samplerCall)]})
	else
	Quote({inds.mat <- rbind(inds.mat,
				 matrix(thisGroup[eval(samplerCall)],
					ncol = B2))}))
  }

  # Determine default value for save.indices
  if(is.null(save.indices)){
    save.indices <- {
      if(byGroup && group.order.matters) F
      else if(n * sum(B) <=  20000) T
      else if(n * sum(B) <= 500000) 2  # save a compressed version
      else F
    }
  }

  # Create objects for storing results
  p <- length(observed)
  reps <- matrix(NA, p, sumB)
  need.indices <- (save.indices || regressionL)
  if(need.indices)
    all.indices <- NULL
  Lstar <- NULL

  on.exit(
    if(totalB) {
      cat("\nDid ", totalB,
	  " replications, saving results in .bootstrap.partial.results, interrupt again to abort completely.\n")
      reps <- t(reps[, 1:totalB, drop = F])
      dimnames(reps) <- list(NULL, names.observed)
      B <- c(B[1:iB][-iB], doneB)
      func.call$B <- B
      anyProb <- any(withProb[1:iB])
      if(anyProb){
	probRatio <- probRatio[1:iB]
	probRatioB <- probRatioB[1:totalB, 1:iB, drop=F]
	prob <- 1/(probRatioB %*% (B/totalB))
      }
      assign(".bootstrap.partial.results", where = 1,
	     immediate = T,
	     makeBootstrap(replicates = reps,
		       observed = observed, n = n, call = func.call,
		       seed.start = seed.start,
		       seed.end = "Unknown, due to interrupt",
		       dim.obs = dim.obs,
		       group = if(save.group) group,
		       subject = if(save.subject) subject,
		       parent.frame = parent.frame,
		       indices = if(save.indices && save.indices<2)
		           all.indices[,1:totalB,drop=F],
		       compressedIndices = if(save.indices == 2)
		           compressIndices(all.indices[,1:totalB,drop=F], n),
		       weights = if(anyProb) as.vector(prob),
		       statistic.is.random = statistic.is.random,
		       label = label,
		       defaultLabel = resampMakeLabel("bootstrap",
			 func.call$data, func.call$statistic),
		       L = L,
		       Lstar=if(!is.null(Lstar)) Lstar[1:totalB,,drop=F]))
    }
	  , add = T)

  # Loop: iB=1:nB
  #         iBlock=1:nblocks
  #           iGroup=1:nGroups (for creating indices)
  #           lapply(1:B2, ...) (for calculating statistic)
  for(iB in seq(length = nB)){
    nblocks <- ceiling(B[iB]/block.size)
    B2 <- block.size
    previousB <- sum(B[1:iB])-B[iB]
    if(anyProb)
      probIB <- sampler.prob[[iB]]

  for(iBlock in 1:nblocks) {
    doneB <- (iBlock - 1) * block.size
    totalB <- previousB + doneB

    if(iBlock == nblocks && B[iB] %% block.size)
      B2 <- B[iB] %% block.size

    if(trace)
      cat("Forming replications ", totalB + 1,
	  " to ", totalB + B2, "\n")

    if(!byGroup)
      inds.mat <- eval(samplerCall)
    else {
      # byGroup case
      inds.mat <- {if(group.order.matters)
		     array(as.integer(0), c(n, B2))
		   else
		     NULL}
      for(iGroup in 1:nGroups)
	eval(getGroupInds)
    }
    if(need.indices)
      all.indices <- cbind(all.indices, inds.mat)
    if(is.numeric(L) && !bySubject)
      Lstar <- rbind(Lstar, indexMeans(L, inds.mat, group=group))

    if(statistic.is.random) {
      seed.sampler <- .Random.seed
      .Random.seed <<- seed.statistic
    }

    if(statistic.is.mean){
      reps[, totalB + 1:B2] <- t(indexMeans(data, inds.mat))
    }
    else {
      for(i in 1:B2){
	tempRep <- fit.func(inds = inds.mat[, i],
			   data = data, statistic = statistic,
			   args.stat = args.stat,
			   subjectIndices = subjectIndices,
			   subjectSizes = subjectSizes, nSubjects = n,
			   subjectName = subjectName)
	if(length(tempRep) != p)
	  stop("statistic returns result with varying length")
	reps[, totalB + i] <- tempRep
      }
    }

    if(anyProb)
      probRatioB[totalB + 1:B2, withProb] <-
	sapply(probRatio[withProb],
	       indexProducts, indices = inds.mat)
    if(statistic.is.random) {
      seed.statistic <- .Random.seed
      .Random.seed <<- seed.sampler
    }
  }
  }
  totalB <- sumB
  seed.end <- .Random.seed

  # Dimension and dimnames for replicates
  reps <- t(reps)
  dimnames(reps) <- list(NULL, names.observed)

  if(anyProb)
    prob <- 1/(probRatioB %*% (B/sumB))

  result <- makeBootstrap(replicates = reps, observed = observed, n = n,
	    call = func.call, seed.start = seed.start, seed.end = seed.end,
	    dim.obs = dim.obs,
	    group = if(save.group) group,
	    subject = if(save.subject) subject,
	    parent.frame = parent.frame,
	    indices = if(save.indices && save.indices<2) all.indices,
	    compressedIndices = if(save.indices == 2)
		compressIndices(all.indices, n),
	    weights = if(anyProb) as.vector(prob),
	    statistic.is.random = statistic.is.random,
            label = label,
	    defaultLabel = resampMakeLabel("bootstrap",
	      func.call$data, func.call$statistic),
	    L=L, Lstar=Lstar)

  on.exit()
  on.exit({cat("Saving results in .bootstrap.partial.results\n")
	   assign(".bootstrap.partial.results", result,
		  where = 1, immediate = T)
	 })
  result$order.matters <- order.matters

  # If L is "regression" or "ace", create L now.
  if(regressionL){
    if(trace)
      cat("Creating L using", L, "method\n")
    Lmethod <- L
    if(is(model.mat, "formula"))
      L <- linearApproxReg(reps, all.indices,
			   formula = model.mat, data = data,
			   n = n, group = group, subject = subject,
			   weights = if(anyProb) as.vector(prob),
			   transform = (L == "ace"))
    else
      L <- linearApproxReg(reps, all.indices,
			   model.mat = model.mat,
			   n = n, group = group, subject = subject,
			   weights = if(anyProb) as.vector(prob),
			   transform = (L == "ace"))
    attr(L, "method") <- Lmethod
    result$L <- L
    if(!bySubject)
      result$Lstar <- indexMeans(L, all.indices, group=group)
  }

  if(trace)
    cat("\n")
  on.exit()
  result
}

# See file "ToDo" for changes not listed here.
# Generic
# No coerce to matrix,
# new resamp.get.fit.func
# passes argument "DimLength" to resamp.get.fit.func (bug 13642)
# first argument to sampler is n, not 1:n
# drop=F so B=101 works
# Fix bug 13637 - seed.statistic is treated like an integer
# B can be vector (bug 13641)
# Importance sampling
# Allow expressions which are named, mode (, or mode { & length 1
# Changed how seed.start and resamp.set.seed are used, email 1/20/99
# Use is.element
# Use resamp.set.seed in place of set.seed
# Fix bug 15115; use is.name() rather than length==1 for determining
#   if substitute(data) is a name; avoid defining substitute.data
# Use numRows
# observed must be atomic, may be 3+array
# Allow 0-length data, useful for graphical bootstrapping (except summaries).
#   Old check did is.null, missed numeric(0).
# Do not include component group (Release Notes) (may change later)
# Do not define "temp"
# Removed warning about new argument for sampler
# Added argument "subject", for sampling by subject
#   (Release Notes - new argument order)
# Use ... in arguments to call.stat, make inds.mat the third argument
# Add parent.frame
# Eliminated must.swap (very rare that statistic is random & sampler is not)
# Define "stratified" to replace !missing(group) so another function
#   can call bootstrap(..., group=something) with something==NULL
# Handle case where both group and subject are supplied.
# Simplified handling of subject argument (do like group)
# Changed argument model.method in call to resampMakeFitObjResults
# Added save.subject and save.group arguments
# Changed all.indices/inds.mat assignments, to allow sample size != n
# Add arguments args.sampler.group, group.order.matters
# Removed dup.row.names code (now use resampSubDF in resampMakeFunc)
# Make default for statistic.is.random = NULL, don't use missing()
# Comments at beginning on the organization of the function
# Added argument L; may be:
#	matrix  (n * p)
#	"jackknife" or "influence": call resampGetL() immediately
#	"regression" or "ace":      call linearApproxReg() at end
# Added model.mat; if supplied is passed to linearApproxReg()
# Change default for save.indices to T if L supplied.
# Rearrange order of some parts near the end.
# Make resampled subjects have unique names, if subject a variable in
#   data frame (add some objects like subjectSizes, change resampMakeFunc).
# Make things work if L="regression" and save.indices=F.
# Let model.mat be a formula
# Add argument `argumentList'; may be a list containing values of arguments
# Use `argumentList' instead of making `sampler' a list.
# Add support for sampler=samp.finite(N=500)
# Set group.order.matters=F if sampling by group and subject,
# Add statistic.is.random to the output (for use in resampGetL and limits.tilt).
# Add options L="choose", "none"; L=NULL now means compute L if n is small
# When trace=T, print note when computing L.
# Revert to old behavior; L=NULL does not compute L.
# If L supplied as vector, convert to matrix.  Check for right # of rows.
# Sort both before comparing names of sampler.args.group against groupNames
# Add statistic.is.mean code, for faster handling of mean and colMeans.
# Stop if L is not numeric or character (e.g. a data frame)
# Add arguments label, statisticNames, components label and defaultLabel
# add option save.indices=2 to save compressed indices
# changed default to save.indices, now saves indices or compressed indices
# Return Lstar if L computed
# If statistic.is.mean and group is supplied, don't use the data as L.
# Lstar = indexMeans(..., group)
# If bySubject and numeric L, assume L is sorted and normalized.
# Add argument resampleColumns
# If sampling with different size, set group.order.matters = F
# Allow L to be logical (for mean of proportions)
# always do as.vector(observed), to handle coef(lm) that is vector w class
# default trace to F when statistic.is.mean
# Remove sampler.setup and sampler.wrapup arguments
#   (delete functions resamp.set.seed, resamp.return.seed)


# Want examples in help file for:
#   use subject
#   multi-way hierarchy
#   stratify on multiple variables




##########################################################
# Jackknife
# now a generic function
##########################################################

# My version:
# Generic function; see lm.S for an `lm' method

jackknife <-
function(data, statistic, args.stat = NULL,
	 group = NULL, subject = NULL,
	 label = NULL,
	 statisticNames = NULL,
	 seed = .Random.seed,
	 group.size = 1, assign.frame1 = F,
	 save.group = NULL, save.subject = NULL, ...)
{
  if(missing(statistic))
    stop("Argument statistic is missing with no default")
  UseMethod("jackknife")
}
# make argument list same as jackknife.default (so gui code works)

jackknife.default <-
function(data, statistic, args.stat = NULL,
	 group = NULL, subject = NULL,
	 label = NULL,
	 statisticNames = NULL,
	 seed = .Random.seed,
	 group.size = 1, assign.frame1 = F,
	 save.group = NULL, save.subject = NULL, ...)
{
  # Delete-one jackknifing (if subject supplied, delete-one-subject).
  # If a group.size is specified then random sets of size group.size
  # are omitted.

  if(missing(statistic))
    stop("Argument statistic is missing with no default")

  # Capture call
  func.call <- match.call()
  parent.frame <- sys.parent()

  # Check if `data' is the output of a modeling function; if so
  # take special action.  For some modeling functions (e.g. lm)
  # there are jackknife methods that are used instead of this.
  if((is.call(fitCall <- data$call) ||
      is.call(fitCall <- attr(data,"call"))) &&
     !is.null(fitCall$data) && class(data) != "model.list")
    return(resampMakeFitObjResults(resampCall=func.call, fitCall=fitCall,
				   frame.eval=parent.frame, model.method=NULL))

  # If statistic isn't a function or name of function or expression,
  # store it as an unevaluated expression to pass to fit.func
  substitute.stat <- substitute(statistic)
  if(!is.element(mode(substitute.stat), c("name", "function")))
    statistic <- substitute.stat

  # Get name of data
  data.name <- substitute(data)
  if(!is.name(data.name))
    data.name <- "data"

  # Faster subscripting for a data frame, allow duplicate row names
  is.df.data <- is.data.frame(data)
  if(is.df.data && is.null(attr(data, "dup.row.names")))
    attr(data, "dup.row.names") <- T

  # Set seed in case statistic uses randomization
  if(!missing(seed)) {
    orig.random.seed <- .Random.seed
    on.exit(.Random.seed <<- orig.random.seed)
    set.seed(seed)
  }

  # How many observations, or subjects if sampling by subject?
  n <- numRows(data)  # This may be changed to the number of subjects

  # Save group or subject arguments?
  if(is.null(save.subject)) save.subject <- (n <= 10000)
  if(is.null(save.group))   save.group   <- (n <= 10000)

  # If data is a data frame, look for subject there first.
  if(!missing(subject) && is.df.data)
    subject <- resampMakeArgument(func.call$subject, data, parent.frame)
  bySubject <- length(subject)
  if(bySubject){
    if(bySubject != n)
      stop("length of subject must match number of observations in data")
    subjectIndices <- split(1:n, subject)
    n <- length(subjectIndices)
  }
  else
    subjectIndices <- NULL

  # Create function to evaluate the statistic
  if(is.null(assign.frame1)) assign.frame1 <- F
  fit.func <- resampMakeFunc(statistic = statistic, data.name = data.name,
			     is.df.data = is.df.data,
			     is.null.args.stat = is.null(args.stat),
			     assign.frame1 = assign.frame1,
			     dimLength = length(dim(data)),
			     bySubject = bySubject)
  # is: function(inds, data, statistic, args.stat, subjectIndices)
  if(length(list(...)))
    warning("unrecognized arguments in function call")

  # Calculate statistic for observed data
  if(assign.frame1)
    on.exit(if(exists(data.name, frame = 1))
	    remove(as.character(data.name), frame = 1))
  observed <- fit.func(inds = 1:n, data = data, statistic = statistic,
		       args.stat = args.stat, subjectIndices = subjectIndices)
  if(!is.atomic(observed))
    stop("Statistic must return numerical (non-list) data.")

  # Get statistic names and coerce to vector
  names.observed <- resampMakeNames(observed = observed,
				    statistic = substitute.stat,
				    statisticNames = statisticNames)
  dim.obs <- dim(observed)
  if(!is.null(dim.obs))
    observed <- as.vector(observed)
  names(observed) <- names.observed


  # Stratified sampling ("group" argument)
  # Affects estimates, not replicates (except prevents group.size > 1).
  if(!missing(group)) {
    # If data is a data frame, look for group there first.
    if(is.df.data)
      group <- resampMakeArgument(func.call$group, data, parent.frame)
    orig.group <- group
    if(length(group)){
      if(group.size > 1){
	warning("resetting group.size = 1 (larger value not currently supported when group is supplied)")
	group.size <- 1
      }
      # Get indices
      if(bySubject){
	# get group for each subject (subject must be nested within group)
	subjectByGroup <- lapply(split(subject, group), unique)
	if(notNested(subjectByGroup))
	  stop("some values of `subject' occur in multiple groups; subject must be nested within the `group' argument")
	group <- rep(names(subjectByGroup), sapply(subjectByGroup, length))
      }
    }
  }
  else
    orig.group <- NULL


  # Delete group.size observations at a time.
  # In stratified (group) case, require group.size = 1, n=B
  # If group.size does not divide n, some points will never be left out.
  # The main motivation for different group.size is in approximating
  # acceleration without doing a full delete-1 jackknife.
  B <- floor(n/group.size)
  drop.inds <- if(group.size == 1) -(1:n)
              else split(sample(x = -(1:n), size = B*group.size),
			 rep(1:B, each = group.size))
  tempReps <- lapply(drop.inds, fit.func,
		 data = data, statistic = statistic,
		 args.stat = args.stat,
		 subjectIndices = subjectIndices)
  if(any(sapply(tempReps, length) != length(observed)))
    stop("statistic returns result with varying length")

  reps <- matrix(unlist(tempReps), ncol = length(observed),
                 nrow = length(tempReps),
		 dimnames = list(NULL, names.observed), byrow = T)

  if(anyMissing(reps))
    warning("NA's encountered in replicates")

  if(length(observed)){
    # Compute mean, bias, SE
    jack.mean <- colMeans(reps)
    jack.bias <- (B - 1) * (jack.mean - observed)
    # In some stratified sampling situations the bias should be NA,
    # but we cannot determine that here.
    if(!length(group)){  # usual case, no stratification
      jack.se <- sqrt((B - 1) * colVars(reps, unbiased = F))
    } else { # stratified sampling; currently require n == B
      f <- function(inds, reps)
        length(inds) * colVars(reps[inds,,drop=F])
      groupVars <- lapply(split(1:n, group), f, reps = reps)
      jack.se <- (n-1)/n * sqrt(rowSums(do.call("cbind", groupVars)))
    }
    est <- data.frame(Mean = jack.mean,
                      Bias = jack.bias,
                      SE = jack.se)
  }
  else
    est <- data.frame(list(Mean=numeric(0), Bias=numeric(0), SE=numeric(0)))

  result <- list(call = func.call,
		 observed = observed,
		 replicates = reps,
		 estimate = est,
		 B = B,
		 n = n,
		 dim.obs = dim.obs,
		 seed.start = seed,
		 defaultLabel = resampMakeLabel("jackknife",
		   func.call$data, func.call$statistic),
                 parent.frame = parent.frame)
  # add optional components, if they are non-null
  result <- c(result,
	     label = label,
	     group = if(save.group && !is.null(group)) list(orig.group),
	     subject = if(save.subject && !is.null(subject)) list(subject))
  oldClass(result) <- c("jackknife", "resamp")
  result
}
# No coerce to matrix, passes argument "DimLength" to resampMakeFunc
# allow data to be vector, matrix, or array.
# Allow expressions which are named, mode (, or mode { & length 1
# Use is.element
# Fix bug 15115; use is.name() rather than length==1 for determining
#   if substitute(data) is a name; avoid defining substitute.data
# Use numRows
# observed must be atomic, may be 3+array
# Allow 0-length data, useful for graphical bootstrapping (except summaries).
#   Old check did is.null, missed numeric(0).
# Add parent.frame
# Add argument subject.
# Change argument model.method in call to resampMakeFitObjResults.
# Add save.subject argument.
# Removed dup.row.names code (now use resampSubDF in resampMakeFunc)
# Do not call jackstats() -- need more info than we can easily pass.
# Add group argument, use when calculating estimates
# Change n.groups to B
# Add argument label and statisticNames, output label and defaultLabel
# Allow missing values in reps (may occur in only some statistics)
# allow zero-length result


###############################################################
# resampMakeFunc
#   Simpler, adds argument DimLength
#
#   For original version, comments, and the result of this function
#   under various circumstances see file
#	resamp.get.fit.func.S
###############################################################

resampMakeFunc <-
function(statistic, data.name, is.df.data, is.null.args.stat,
         assign.frame1 = F, dimLength, bySubject = F, subjectInDF = F,
	 resampleColumns = NULL)
{
  # Construct a function which takes arguments
  #  inds:      indices for subscripting
  #  data:      vector, matrix, array, or data frame
  #  statistic: function or expression
  #  args.stat: optional arguments or objects (may be NULL)
  #  subjectIndices: optional list containing observations by subject
  #  resampleColumns:  if not null, then subscript only these columns
  #
  #  Below are used only if subjectInDF = T (subject a variable in data frame),
  #  in which case the function replaces the subject variable in the data frame
  #  so that each resampled subject has a unique value.
  #
  #  subjectSizes: number of observations per subject
  #  nSubjects: number of subjects
  #  subjectName: name of the subject variable in the data frame

  ffs <- paste("function(inds, data, statistic, args.stat, subjectIndices,",
               "subjectSizes, nSubjects, subjectName){")

  # create expression to subscript data
  # Typically "data <- data[inds,,drop=F]", but may vary
  indsText <- "inds"
  if(bySubject)
    indsText <- "unlist(subjectIndices[inds])"
  
  if(is.null(resampleColumns)){
    ffs2 <- {
      if(is.df.data)		   
	paste("resampSubDF(data,", indsText, ")") # eg "resampSubDF(data, inds)"
      else
	paste("data[", indsText,
	      if(dimLength > 1)
	      paste(c(rep(",",dimLength),"drop=F"), collapse=""),
	      "]")
    }
    ffs <- paste(ffs, "data <-", ffs2, "\n")
  }
  else {  # subscript only some columns.  Only support matrix & data.frame
    # typical output here:  data[, cols] <- data[inds, cols]
    cols <- deparse(resampleColumns)
    if(is.df.data)
      ffs2 <- paste("data[ ", cols, "] <- resampSubDF(data[", cols, "], ",
		  indsText, ")")
    else
      ffs2 <- paste("data[,", cols, "] <- data[", indsText, ",", cols, "]")
    ffs <- paste(ffs, ffs2, "\n")
  }

  if(subjectInDF){
    # Replace subject variable with unique subject names (1, 2, ... nSubjects)
    ffs <- paste(ffs, "\"$\"(data, subjectName) <-",
                 "rep(1:nSubjects, subjectSizes[inds])\n")
  }

  if(assign.frame1)
    ffs <- paste(ffs, "assign('", data.name, "', data, frame=1)\n", sep="")

  # Handle function or expression
  if(is.function(statistic))
    ffs <- paste(ffs,
		 if(is.null.args.stat)
                   "statistic(data)}"
		 else
                   "do.call('statistic',c(list(data), args.stat))}")
  else
    ffs <- paste(ffs,
		 "eval(statistic,c(list(", data.name, "=data)",
		 if(is.df.data)
		   ", data",
		 if(!is.null.args.stat)
		   ", args.stat",
		 "))}")
  eval(parse(text = ffs))
}
# Let bootstrap take expressions which are named, mode (, or mode { & length 1
# Add argument bySubject, to support sampling by subject
# Removed argument substitute.stat (not used)
# Renamed from resamp.get.fit.func to resampMakeFunc
# Make assign.frame=NULL equivalent to F
# Add argument subjectInDF to resampMakeFunc, add arguments subjectSizes,
#   nSubjects, subjectName to returned function, to make subject names unique.
# add resampleColumns

# Could handle the case that statistic is a function by doing
# the subscripting in-line, e.g. statistic(data[inds]) rather than
# first data <- data[inds] then calling statistic.  Would that save
# a lot of memory (depending on whether statistic mucks with data)?


###############################################################
# resampMakeFuncSimple
#   Like resampMakeFunc, but without arguments
#	inds or dimLength
#   because result doesn't have arguments
#	indices, subjectIndices
###############################################################

resampMakeFuncSimple <- function(statistic, data.name, is.df.data,
					is.null.args.stat, assign.frame1 = F){
# Construct a function which takes arguments
#  data:      vector, matrix, or data frame
#  statistic: function or expression
#  args.stat: optional arguments or objects (may be NULL)

        ffs <- paste("function(data, statistic, args.stat){")

        if(assign.frame1) ffs <- paste(ffs, "assign('", data.name,
                        "', data, frame=1)\n", sep = "")

        # Handle function or expression
        if(is.function(statistic))
                ffs <- paste(ffs, if(is.null.args.stat) "statistic(data)}"
                         else "do.call('statistic',c(list(data), args.stat))}")
        else ffs <- paste(ffs, "eval(statistic, c(list(",
                        data.name, "=data)",
			if(is.df.data) ", data",
			if(!is.null.args.stat) ", args.stat",
			"))}")
        eval(parse(text = ffs))
	}
# Renamed from Steve Ellis's resamp.get.pfit.func
# removed as(...,"list")  (that SV4 bug is fixed)


##########################################################
# resampMakeNames
#	avoid creating names.observed
#	eliminate as.vector line
#	avoid call to outer()
#	use stat.name for matrix case
#	Work right if stat is length 1 & mode {
#	Support statistics which return 1-d arrays
##########################################################
resampMakeNames <- function(observed, statistic, default.stat.name = "Param",
                            prepend.stat.name = F,
			    statisticNames = NULL){
  # Create names for statistics, based on observed statistic, etc.
  dim.obs <- dim(observed)
  p <- length(observed)
  if(length(statisticNames) == p)
    return(statisticNames)
  
  stat.name <-
    if(is.name(statistic)) deparse(statistic) else default.stat.name

  if(length(dim.obs) == 2){
    # Observed is a matrix, use row and column names if present
    obs.dimnames <- dimnames(observed)
    f <- function(dn, n){
      if(length(dn))
        abbreviate(dn, 5)
      else seq(length = n)
    }
    paste(if(is.null(obs.dimnames))
	    paste(stat.name, seq(length = dim.obs[1]), sep = "")
          else{
            if(prepend.stat.name)
              paste(stat.name, f(obs.dimnames[[1]], dim.obs[1]), sep = ".")
            else
              f(obs.dimnames[[1]], dim.obs[1])
          },
	  rep(f(obs.dimnames[[2]], dim.obs[2]), each = dim.obs[1]),
	  sep = ".")
  }
  else {
    # Observed is a vector, 1-d array, or 3+dimensional array.
    if(is.null(names(observed)))
      paste(rep(stat.name, p), if(p > 1) 1:p, sep = "")
    else{
      if(prepend.stat.name)
        paste(stat.name, names(observed), sep = ".")
      else
        names(observed)
    }
  }
}
# Bad names if statistic name ends in a digit & returns a vector, but don't fix
# Fixed bugs 14263, 14252
# Support statistics which return 1-d arrays
# Support statistics which return 3+arrays
# Support statistics which are length 0 (useful for graphical bootstrapping).
# Arg name change from stat to statistic


##########################################################
# resampMakeLabel
##########################################################
resampMakeLabel <-
function(method, dataExpr, statisticExpr,
	 data2Expr = NULL,
	 treatmentNames = NULL, ratio = F){
  # Create default label for a resamp object, e.g.:
  #     bootstrap: x: mean
  #     bootstrap: x: mean: difference Female - Male
  # method is e.g. "bootstrap" or "permutations"
  # dataExpr, statisticExpr, data2Expr are expressions used in original call
  # treatmentNames, if present, are names for 2 treatments 
  #   (ignore it if data2Expr is present)
  # ratio is logical, indicating difference or ratio between 2 treatments

  # There are three main cases:  one sample, data2, treatment
  # Need text for data and statistic in any case

  # Create text for data
  if(is.name(dataExpr))
    dataText <- as.character(dataExpr)
  else {
    dataText <- deparse(dataExpr)
    if(length(dataText) > 1)
      dataText <- paste(substring(dataText[1], 1, 17), "...", sep="")
    else if(nchar(dataText) > 20)
      dataText <- paste(substring(dataText, 1, 17), "...", sep="")
  }

  # Create text for statistic
  if(is.name(statisticExpr))
    statisticText <- as.character(statisticExpr)
  else {
    statisticText <- deparse(statisticExpr)
    if(length(statisticText) > 1)
      statisticText <- paste(substring(statisticText[1], 1, 17), "...", sep="")
    else if(nchar(statisticText) > 20)
      statisticText <- paste(substring(statisticText, 1, 17), "...", sep="")
  }

  # Function to use for pasting together parts of label
  pasteC <- paste
  pasteC$sep <- " : "

  if(!is.null(data2Expr)){  # First of the three main cases
    # Create text for data2 (if present)
    if(is.name(data2Expr))
      data2Text <- as.character(data2Expr)
    else {
      data2Text <- deparse(data2Expr)
      if(length(data2Text) > 1)
	data2Text <- paste(substring(data2Text[1], 1, 17), "...", sep="")
      else if(nchar(data2Text) > 20)
	data2Text <- paste(substring(data2Text, 1, 17), "...", sep="")
    }
    result <- pasteC(statisticText,
		    paste(dataText, data2Text,
			  sep = if(ratio) " / " else " - "))
  }
  else { # Handle the second and third main cases together
    result <- pasteC(dataText, statisticText)
    # Now a special case; for e.g. bootstrap(x, mean(x, trim=.1))
    # make that just "mean(x, trim=.1")
    if(!is.name(statisticExpr) && is.name(dataExpr) &&
       !is.function(statisticExpr) &&
       is.element(dataExpr, all.names(statisticExpr)))
      result <- statisticText
    if(!is.null(treatmentNames)){ # Second of the three main cases
      result <- pasteC(result,
		      paste(as.character(treatmentNames[1:2]),
			    collapse = if(ratio) " / " else " - "))
    }
  }
  paste(method, result, sep=" : ")
}
# Add ratio argument



##########################################################
# resampMakeArgument
#	Evaluate an argument, with scope a data frame and parent frame
#	See ~timh/bootstrap/arguments.S for background
##########################################################
resampMakeArgument <- function(expr, data, parent.frame){
  # Evaluate the expression that defines an argument,
  # looking first in a data frame and then in another frame.
  # E.g. if f() calls bootstrap(data=df, ..., group=A)
  # then A could be a variable in df, or defined in f's frame.
  if(!any(is.element(all.names(expr), names(data)))){
    eval(expr, parent.frame)
  }
  else if(is.name(expr)){
    # The following expression might print as `data$expr', but
    # is really `"$"(data, expr)':
    "$"(data, expr)
  }
  else {
    eval(expr,
	 if(is.list(parent.frame)) c(parent.frame, data)
	 else if(parent.frame > 1) c(sys.frame(parent.frame), data)
	 else data)
  }
}





##########################################################
# makeBootstrap
##########################################################
makeBootstrap <-
function(replicates, observed, n, call, seed.start, seed.end,
	 dim.obs = NULL, weights = NULL,
	 ...)
{
  ## Calculate bootstrap statistics and create a "bootstrap" object.
  if(anyMissing(replicates)) {
    warning("NA's encountered in replicates.  Mean and standard error are calculated using remaining observations")
    B.missing <- sum(rowSums(is.na(replicates)) > 0)
  }
  else
    B.missing <- NULL
  if(length(observed)){
    means <- colMeans(replicates, na.rm = T, weights=weights)
    est <- data.frame(Mean = means,
		      Bias = means - observed,
		      SE = colStdevs(replicates, na.rm = T, weights=weights))
    colRange <- apply(replicates, 2, range, na.rm=T)
    if(anyMissing(colRange))
      warning("All values are missing for one or more statistics")
    if(is.logical(all.equal(colRange[1,], colRange[2,]))
       && !all(is.na(colRange)))
      warning("Bootstrap replicates are (nearly) identical; does your statistic depend on the resampled data?  In some cases bootstrap(..., assign.frame1=T) helps.")
  }
  else
    est <- data.frame(list(Mean=numeric(0), Bias=numeric(0), SE=numeric(0)))
  result <- list(call = call,
		 observed = observed,
		 replicates = replicates,
		 estimate = est,
		 B = dim(replicates)[1],
		 n = n,
		 dim.obs = dim.obs,
		 seed.start = seed.start,
		 seed.end = seed.end,
		 B.missing = B.missing,
		 weights=weights,
		 ...)
  # Components after seed.end may be NULL

  # Strip certain components if they are NULL
  for(i in c("B.missing", "group", "subject", "treatment",
	     "indices", "weights", "L", "label",
             "Lstar", "compressedIndices")){
    if(is.null(result[[i]]))
      result[[i]] <- NULL
  }

  oldClass(result) <- c("bootstrap", "resamp")
  result
}
# I avoided if(!is.null()) lines
# Move Bias to second column.
# Add ... to provide support for other functions.
# Do not include group   as argument or component (unless in ... list)
# Do not include indices as argument or component (unless in ... list)
# Strip some components if they are NULL, like group, subject, L
# Use apply(...,2,range) for checking constant columns instead of colStdevs.
# More informative message when all values missing.
# Catch approximately identical replicates.
# Add some additional components that are stripped if null
# add compressedIndices to list of components stripped if NULL


# I'm uneasy with calculating moments as a matter of course.
# That might be better done in print or summary.
# Moments are not appropriate for all problems (or at least, "bias"
#  is meaningless for permutation tests).


##########################################################
# jackstats
##########################################################

jackstats <- function(replicates, observed, n,
		      call, seed.start, dim.obs = NULL, B = n,
		      strata = NULL, ...)
{
  # This function is deprecated

  if(anyMissing(replicates)) {
    assign(".jackknife.replicates", replicates, where = 1, immediate = T)
    stop("NA's encountered in replicates.  Replicates stored as .jackknife.replicates.")
  }
  jack.mean <- colMeans(replicates)
  jack.bias <- (B - 1) * (jack.mean - observed)
  # In some stratified sampling situations the bias should be NA
  if(!length(strata)){  # usual case, no stratification
    jack.se <- sqrt((B - 1) * colVars(replicates, unbiased = F))
  }
  else { # stratified sampling;  currently require n == B
    f <- function(inds, reps)
      length(inds) * colVars(reps[inds,,drop=F])
    groupVars <- lapply(split(1:n, group), f, reps = reps)
    jack.se <- (n-1)/n * sqrt(rowSums(do.call("cbind", groupVars)))
  }
  est <- data.frame(Mean = jack.mean,
		    Bias = jack.bias,
		    SE = jack.se)

  result <- list(call = call,
		 observed = observed,
		 replicates = replicates,
		 estimate = est,
		 B = B,
		 n = n,
		 dim.obs = dim.obs,
		 seed.start = seed.start,
		 ...)
  oldClass(result) <- c("jackknife", "resamp")
  result
}
# Moved Bias to second column
# Make distinction between n and n.groups
# Quicker searching for missing values
# If stop save replicates as "jackknife.replicates"
# Add ... (use it for call.stat
# Add strata argument (group, or shortened version if subjects)
# Change n.groups to B
# Make this deprecated



##########################################################
# resampGetIndices
##########################################################
resampGetIndices <-
function(object, frame.eval = object$parent.frame)
{
  # Get the indices for a bootstrap object (create if necessary)
  if(!is.null(object$indices))
    return(object$indices)
  if(!is.null(object$compressedIndices))
    return(uncompressIndices(object$compressedIndices))
  # Indices were not stored; recreate them.
  boot.call <- object$call
  if(is.null(frame.eval))
    frame.eval <- sys.parent(1)
  boot.call$statistic <- function(...) NULL
  # fast function, assume actual values unneeded
  boot.call$assign.frame1 <- F
  boot.call$save.indices <- T
  boot.call$seed <- object$seed.start	#in case seed modified
  boot.call$L <- "none"
  boot.call$trace <- F
  eval(boot.call, frame.eval)$indices
}
# B can be vector (bug 13641), no update.B
# Change `x' to `boot.object', later to `object'.
# Make statistic return NULL, to avoid "identical values" warning
# Change default for frame.eval.boot
# Set assign.frame1 = F
# Changed frame.eval.boot to parent.frame
# look for compressed indices if indices not found


##########################################################
# resampGetArgument
##########################################################
resampGetArgument <- function(object, argument, checkObject = NULL,
			     frame.eval = object$parent.frame){
  # object is a resamp object (e.g. bootstrap)
  # argument is the quoted or unquoted name of an argument that
  #  may be included in the call that created object,
  #  such as "group" or "subject".
  # frame.eval is the frame in which to look for the definitions
  #  of the argument.
  # If the bootstrap object has a component with the same name, return that.
  # Otherwise, evaluate the argument as it would have been evaluated;
  #  if the data argument was a data frame, look there first.
  sarg <- substitute(argument)
  if(is.name(sarg))
    argument <- as.character(sarg)

  # If appropriate, look first inside object
  if(is.null(checkObject))
    checkObject <- (argument != "B")
  if(checkObject && !is.null(object[[argument]]))
    return(object[[argument]])

  # Need an expression to evaluate (maybe just the name of a variable)
  # Essentially object$call[[argument]], but avoid partial matching
  j <- match(argument, names(object$call))
  if(is.na(j))
    return(if(argument == "B") object$B else NULL)
  argCall <- object$call[[argument]]
  if(is.null(argCall))
    return(NULL)

  if(is.null(frame.eval))
    frame.eval <- sys.parent(1)

  # Evaluate the expression; look in a data frame first, then frame.eval
  data <- eval(object$call$data, frame.eval)
  if(is.data.frame(data))
    resampMakeArgument(argCall, data=data, parent.frame=frame.eval)
  else
    eval(argCall, frame.eval)
}
# add argument frame.eval (RT, 7/16/01)
# look in `object' first for component with same name as argument.
# add argument checkObject, special handling when argument = "B".
# avoid partial argument matching
# Fix handling of B, when B not in original call.


##########################################################
# limits.percentile
##########################################################
limits.percentile <- 
function(x, probs = c(25, 50, 950, 975)/1000,
	 subset.statistic=T, narrow=F){
  # Calculate percentiles of replicates in a "resamp" object.
  # When bootstrapping these are known as bootstrap percentile limits.
  y <- apply(x$replicates[,subset.statistic,drop=F], 2,
	    quantile,
	    probs = probs, na.rm = T, weights=x$weights,
	    alpha=narrow, rule=2)
  if(anyMissing(x$replicates[,subset.statistic,drop=F]))
    warning("Missing values omitted, percentiles may be biased")
  # return a matrix with length(probs) columns
  if(length(probs) == 1)
    matrix(y, ncol=1, dimnames=list(names(y),
			paste(100 * probs, "%", sep = "")))
  else {
    dimnames(y)[[1]] <- paste(100 * probs, "%", sep = "")
    t(y)
  }
}
# Use weights
# Create matrix in one step, rather than adding dimnames.  Avoid
#   unnecessary as.matrix() call.
# quantile adds names, so we don't need to.
# Use alpha=0 in computing quantiles by default
#  (will give wider quantiles than before).
# Add subset.statistic to allow computation for only part of statistic.
# Return matrix when length(probs)==1.
# Add names to quantile dimension of the result. (work around bug 24156)
# Change name to limits.percentile, initial comments.  (new limits.emp below)
# Add warning about missing values.


##########################################################
# limits.emp
##########################################################
limits.emp <- function(x, ...){
  # This function is deprecated.  Use limits.percentile() for
  # bootstrap percentile intervals, or for quantiles.

  limits.percentile(x, ...)
}


##########################################################
# limits.bca
##########################################################
limits.bca <-
function(boot.obj, probs = c(25, 50, 950, 975)/1000, details = F, z0 = NULL,
	 acceleration = NULL, group.size = NULL,
	 frame.eval = boot.obj$parent.frame,
	 subset.statistic=1:p, narrow=F, ...)
{
  # Calculate BCa confidence limits.
  # Or, specify acceleration = 0 to calculate BC limits (bias-corrected).
  # z0 is a bias correction.
  # If you do not specify the acceleration, jackknifing is used.
  # Use group.size to avoid doing all n jacknife replicates;
  # the default yields roughly 20 jackknife replicates, unless
  # stratified bootstrap sampling was used.
  if(!is(boot.obj, "bootstrap"))
    stop("obj must be a 'bootstrap' object.")
  if(is.null(frame.eval))
    frame.eval <- sys.parent(1)
  p <- length(boot.obj$observed)

  save.group.size <- F
  # Calculate acceleration if not specified.
  if(is.null(acceleration)) {
    # Calculate acceleration based on linear approximation L,
    # or on a partial jackknife.

    message <- resampCheckIfOrderMatters(boot.obj$order.matters)
    if(!is.null(message)){
      warning("Cannot calculate BCa limits--", message)
      return(matrix(NA, nrow=length(subset.statistic), ncol=length(probs)))
    }

    # Start with group and treatment vectors:
    group <- resampGetArgument(boot.obj, "group")
    treatment <- resampGetArgument(boot.obj, "treatment")
    # If L is present in boot.obj, use it for calculating acceleration
    L <- boot.obj$L
    # Otherwise calculate it
    if(is.null(L) && is(boot.obj, "bootstrap2")){
      # For a bootstrap2 object, do not use the jackknife methods below
      L <- resampGetL(boot.obj, frame.eval = frame.eval)
    }
    else if(is.null(L)){
      if(!is.null(boot.obj$treatment))
	stop("The bootstrap object has a treatment component.  That is not expected at this point.")
      save.group.size <- T
      # Calculate acceleration from jackknife values.
      # If stratified sampling, or group.size=1, then use resampGetL
      if(is.null(group.size))
	group.size <- max(1, floor(boot.obj$n/20))
      boot.call <- boot.obj$call
      if(group.size == 1 || !is.null(group)){
	# use group.size = 1
	if(!missing(group.size) && group.size > 1)
	  warning("group.size ignored, due to stratified bootstrap sampling")
	L <- resampGetL.jackknife(boot.call, frame.eval = frame.eval)
        group.size <- 1
      }
      else {
	# not stratified, leave-out-group.size jackkife
	jack.call <- modifyCall(boot.call, NEWFUN = "jackknife",
				"data", "statistic", "args.stat",
				"assign.frame1",
				seed = boot.obj$seed.end,
				group.size = group.size,
				subject = resampGetArgument(boot.obj, "subject",
				  frame.eval = frame.eval))
	jack.obj <- try(eval(jack.call, frame.eval))
	if(is(jack.obj, "Error")){
	  warning("Call to jackknife (to compute acceleration) failed,",
		  "\nso cannot compute BCa confidence limits.",
		  "\nThis is usually because only part of the relevant",
		  "\ndata is being resampled, so that lengths differ.",
		  "\nThe original error message is:\n",
		  unclass(jack.obj))
	  # return NULL, in case this is called from summary()
	  return(NULL)
	}
	# In this case L is not an influence function approximation
	L <- -subtractMeans(jack.obj$replicates)
	acceleration <- colSums(L^3)/(6 * colSums(L^2)^1.5)
      }
    }
    if(is.null(acceleration)){
      # Ensure that L has mean 0 (for each group and treatment)
      gg <- {if(is.null(group)) treatment
			    else if(is.null(treatment)) group
			    else paste(match(group, unique(group)),
				       match(treatment, unique(treatment)))}
      subject = resampGetArgument(boot.obj, "subject",
	frame.eval = frame.eval)
      if(!is.null(subject)){
	# Need one value of gg for each subject (to match L)
        subjectIndices <- split(1:length(subject), subject) 
	gg <- sapply(subjectIndices, function(i, gg) gg[i[[1]]], gg = gg)
      }
      L <- subtractMeans(L, gg)
      acceleration <- colSums(L^3)/(6 * colSums(L^2)^1.5)
    }
  }
  else
    acceleration <- rep(acceleration, length=p)

  # Get zprobs
  zprobs <- qnorm(probs)
  # Get z0
  if(is.null(z0))
    z0 <- qnorm(colMeans(boot.obj$replicates <
			rep(boot.obj$observed, each=boot.obj$B),
			na.rm=T, weights=boot.obj$weights) +
		colMeans(boot.obj$replicates ==
			rep(boot.obj$observed, each=boot.obj$B),
			na.rm=T, weights=boot.obj$weights) / 2)
  if(anyMissing(z0[subset.statistic]))
    warning("Missing values in z0, corresponding intervals will be NA")
  if(any(abs(z0[subset.statistic]) > 0.25, na.rm=T))
    warning("z0 is outside range (-.25, .25), this indicates extreme bias, assumptions underlying BCa interval may be violated; value(s) are:  ",
	    paste(format(z0[subset.statistic]), collapse=" "))
  names(z0) <- names(boot.obj$observed)

  # Check if acceleleration is undefined, take measures.
  if(anyMissing(acceleration[subset.statistic])){
    if(any(is.na(acceleration[subset.statistic]) &
	   !is.na(z0[subset.statistic])))
    warning("Default values of acceleration are undefined for ",
	    sum(is.na(acceleration[subset.statistic])), " variable(s).\n",
	    "Perhaps the jackknife replicates are identical?\n",
	    "BCa limits are undefined; switching to BC limits,\n",
	    "by setting acceleration = 0 for those variables.")
    acceleration[is.na(acceleration)] <- 0
  }

  # Check for missing values
  if(anyMissing(boot.obj$replicates[,
	subset.statistic[ !is.na((z0+acceleration)[subset.statistic])]]))
    warning(sum(rowSums(is.na(boot.obj$replicates[,subset.statistic,drop=F])) > 0),
	    " bootstrap samples have missing values; this may cause bias")

  # Get BCA percentiles
  emp.probs <- bca.percent <- matrix(nrow = p, ncol = length(probs))
  for(i in subset.statistic) {
    emp.probs[i, ] <- pnorm(z0[i] + (z0[i] + zprobs)/
			    pmax(0, (1 - acceleration[i] * (z0[i] + zprobs))))
    bca.percent[i, ] <- quantile(boot.obj$replicates[, i],
				 probs = emp.probs[i, ],
				 na.rm = T, weights=boot.obj$weights,
				 alpha=narrow, rule=2)
  }
  temp <- range(emp.probs[subset.statistic,], na.rm=T)
  if(any(temp[1], 1-temp[2]) * boot.obj$B < 5)
    warning("bootstrap sample size is inadequate for accurate BCA limits")
  dimnames(bca.percent) <- list(names(boot.obj$observed),
				paste(100 * probs, "%", sep = ""))
  if(!missing(subset.statistic)){
    bca.percent <- bca.percent[subset.statistic, , drop=F]
    emp.probs <- emp.probs[subset.statistic, , drop=F]
    z0 <- z0[subset.statistic]
    acceleration <- acceleration[subset.statistic]
  }
  if(!details)
    return(bca.percent)
  dimnames(emp.probs) <- dimnames(bca.percent)
  if(save.group.size)
    list(limits = bca.percent, emp.probs = emp.probs, z0 = z0,
         acceleration = acceleration, group.size = group.size)
  else
    list(limits = bca.percent, emp.probs = emp.probs, z0 = z0,
         acceleration = acceleration)
}
# Fix bug 12313 -- give warning if abs(z0) > .25
# Use full names for components, not $obs $rep $est
# Avoid loop in calculating z0.
# Warning about z0 applies also if user supplied z0.
# Minor change at end
# Don't bother eliminating the loop to computed bca.percent
# Slightly faster calculation of adj.rep, do "-" in next step.
# Make bca interval undefined (need to test this)
# Use alpha=0 in computing quantiles by default
#  (will give wider quantiles than before).
# Use rule=2 to prevent NA's.
# Do not return group.size (is not always calculated).
# Do not define accel; use acceleration instead
# Change default for frame.eval.jack
# Strip names from groupSizes to work around bug #19847; if that is
#   fixed then don't need to define groupSizes.
# Quicker calculations in stratified case when strata variable is sorted.
# Added frame.eval argument to call to resampGetArgument
# Add argument `frame.eval', replace internal use of `frame.eval.jack'.
# Look for L inside boot.obj
# Use resampGetL.jackknife in stratified case or when group.size=1
# Use modifyCall to create jack.call
# Avoid re-evaluating subject.
# Go ahead and return group.size after all, if it is used.
# Simplify handling of extreme limits (prevent negative denominator)
# Add ... argument.  Though extra args ignored, need for summary.bootstrap.
# Handle bootstrap2 objects
# normalize by group and treatment
# Use try() to avoid problem when calling jackknife
# Allow missing values in z0 (warn, not fail)
# Better warnings if NAs in some dimensions
# Remove frame.eval.jack
# Return NA if resampling residuals


# To do:
# make result undefined when emp.probs relationship is not increasing.
# Give warning if B*emp.probs < 5




##########################################################
# plot.resamp
##########################################################
plot.resamp <-
function(x, nrow = NULL, grid.layout = T, rugplot = F,
	 nclass.func = nclass.fd,
	 bandwidth.func = bandwidth.nrd,
	 subset.statistic = 1:p,
	 hist.col = if(getPlatform() == "WIN386") 6 else 2,
	 key = T, corner = c(1,1), key.args = list(between=2:0),
	 ..., xlim=NULL, ylim=NULL, inside=F, nclass = NULL,
	 main = NULL)
{
  # Histogram of resample values with mean and observed values marked.
  # Includes smooth density estimate and optionally a rug plot of the values.
  p <- length(x$observed)
  P <- length(subset.statistic)
  if(is.character(subset.statistic)){
    j <- match(subset.statistic, names(x$observed))
    if(anyMissing(j))
      stop("subset.statistic names not recognized:\n    ",
	   paste(subset.statistic[which.na(j)], collapse=", "))
    subset.statistic <- j
  }
  if(grid.layout && P != 1) {
    if(is.null(nrow)){
      if(length(x$dim.obs) == 2)
	old.par <- par(mfcol = x$dim.obs)
      else {
        nrow <- min(P, 2)
        old.par <- par(mfrow = c(nrow, ceiling(P/nrow)))
      }
    }
    else
      old.par <- par(mfrow = c(nrow, ceiling(P/nrow)))
    on.exit(par(old.par))
  }
  # Main title; if not supplied then use x$label or x$defaultLabel
  if(is.null(main))
    main <- x$label
  if(is.null(main))
    main <- x$defaultLabel

  # Record the input value (if any) of xlim and ylim, for reference when
  # plotting multiple statistics
  xlim.original <- xlim
  ylim.original <- ylim
  for(i in subset.statistic) {
    xi <- x$replicates[, i]
    if(length(wna <- c(which.na(xi), which.inf(xi)))){
      xi <- xi[-wna]
      weights <- x$weights[-wna]  # this works too if weights is NULL
    }
    else weights <- x$weights
    hist.vals <- hist(xi,
		      nclass = if(is.null(nclass)) nclass.func(xi) else nclass,
		      probability = T, plot = F,
		      include.lowest = T, weights=weights, ...)
    dens <- density(xi, width = bandwidth.func(xi),
		    from = hist.vals$breaks[1],
		    to = hist.vals$breaks[length(hist.vals$breaks)],
		    na.rm = T, weights=weights)

    ymax <- max(c(hist.vals$counts, dens$y))
    if(is.null(ylim.original))
      ylim <- c(0, ymax)
    if(is.null(xlim.original))
      xlim <- range( c(hist.vals$breaks), x$observed[i])
    if(inside) {
      barplot(hist.vals$counts, width = hist.vals$breaks, histo = T, ...,
	      xlim = xlim, ylim = ylim, xlab = names(x$observed)[i],
	      ylab = "Density",
	      inside = inside, col=hist.col,
	      main=if(is.null(main)) names(x$observed)[i] else main)
    }
    else {
      # Use polygon() instead of barplot(), because barplot gives
      # phantom white space between some bars.
      plot(xlim, ylim, type="n", ...,
	   xlim = xlim, ylim = ylim, xlab = names(x$observed)[i],
	   ylab = "Density",
	   main=if(is.null(main)) names(x$observed)[i] else main)
      dotArgs <- list(...)
      if(length(dotArgs))
	dotArgs <- dotArgs[is.element(names(dotArgs),
				     c("density","angle","border"))]
      # Set up polygon.  Avoid drawing height-zero bars.
      X <- hist.vals$breaks
      Y <- hist.vals$counts
      xx <- rep(X, each=3)
      yy <- c(0, rep(Y, each=3), 0,0)
      zero <- which(Y==0)
      if(length(zero))
	yy[3*zero] <- NA
      nd <- !duplicated(list(xx,yy))
      # Do not include col; use hist.col instead
      do.call("polygon", c(list(xx[nd], yy[nd],
				border = F, col = hist.col),
			   dotArgs))
    }
    # In the following line x$estimate$Mean may be NULL
    Mean <- x$estimate$Mean[i]
    # If mean was infinite, treat as NULL
    if(!length(Mean) || !is.finite(Mean))
      Mean <- NULL
    line.vals <- approx(dens, xout = c(x$observed[i], Mean))
    line.vals$y <- pmax(line.vals$y, max(dens$y)/4, na.rm=T)
    lines(rep(line.vals$x[1], 2), c(0, line.vals$y[1]))
    both <- length(line.vals$x) > 1
    if(both)
      lines(rep(line.vals$x[2], 2), c(0, line.vals$y[2]), lty = 2)
    lines(dens)
    if(rugplot)
      rug(xi)
    dots <- list(...)
    if(!is.element("axes", names(dots)) || dots$axes){
      if(length(dots)){
	do.call("box", dots[is.element(names(dots), names(par()))])
	# like box(...), but only pass par parameters like bty, col
      }
      else
	box()
    }
    if(key){
      points(x = x$observed[i], y = mean(0, par("usr")[1]), pch=1)
      if(both)
	points(x = Mean, y = mean(0, par("usr")[1]), pch=3)
      both2 <- c(T,both)
      if(is.null(key.args$between))
	key.args$between <- 2:0
      do.call("key", c(list(corner = corner,
			    text = c("Observed", "Mean")[both2],
			    lines = list(lty=(1:2)[both2]),
			    points = list(pch=c(1,3)[both2])),
		       key.args))
    }
  }
  invisible(NULL)
}
# Use full names for components, not $est
# Avoid using old.par <- par(new values) trick.
# Use {} to clarify if() if() else() else()
# Support dim.obs=scalar
# Change nstats to p
# Faster NA removal.
# Took out warning() if all missing values -- would have crashed anyway
#   due to missing paste()
# Use is.element
# Partial support for weights (density is wrong)
# Use nclass.fd instead of in-line definition (see Charlie email)
# Remove in-line definition of bandwidth.nrd - it is now built-in (Charlie).
# Don't fail if there is no x$estimate$Mean
# Oops, restore the old.par <- par(new values) trick (bug 24874)
# Add argument subset.statistic
# Added arguments key and corner, do key by default.
# Add argument inside; do not draw lines inside hist, by default.
# If axes=F, do not call box()
# Add xlim argument; by default include the observed value within xlim
# For the two vertical lines, make them go at least 1/4 of the way up.
# Use polygon() instead of barplot() if inside=F.
# Add key.args.  Use between=2:0 when calling key.
# Add argument col.hist; use it instead of col for the histograms;
#   this way you don't need to specify col, which affects axes too.
# Handle xlim and ylim when plotting multiple statistics
# Avoid drawing zero-height bars for histograms (matters for pdf.graph).
# Changes to support infinite values in replicates and Mean.
# Add argument nclass, if supplied is used instead of nclass.func
# Add title, change default xlab from "Value" to name of observed
# Use box(...) instead of box(), to support bty


##################################################
# nclass.fd
##################################################
# was here, now removed; corrected version is in S+



##########################################################
# print.resamp
##########################################################
print.resamp <-
function(x, digits = max(options()$digits - 3, 4), ...,
	 printCall = .Options.resample$printCall)
{
  # print a resample object
  # printCall: logical, whether to print the call.
  #            numeric, print if call has this many or fewer characters.
  #            Otherwise print the label or defaultLabel
  old.digits <- options(digits = digits)
  on.exit(options(old.digits))
  if(!is.null(x$label) && nchar(x$label))
    cat("Label:  ", x$label, "\n\n", sep="")
  if(is.null(printCall))
    printCall <- .Options.resample$printCall
  if(is.numeric(printCall))
    printCall <- (sum(nchar(deparse(x$call))) <= printCall)
  if(printCall){
    cat("Call:\n")
    print(x$call)
    if(length(x$actual.calls)){
      cat("\nOther calls used in creating the object:\n")
      for(i in seq(along=x$actual.calls))
	print(x$actual.calls[[i]])
    }
    cat("\n")
  }
  else if(is.null(x$label) && !is.null(x$defaultLabel))
    cat(x$defaultLabel, "\n\n")
  cat("Number of Replications:", x$B, "\n")
  if(!is.null(x$B.missing))
    cat("\nReplicates with missing values:", x$B.missing, "\n")
  cat("\nSummary Statistics:\n")
  if(length(x$observed))
    print(cbind(Observed = x$observed, x$estimate))
  else
    cat("(None - statistic had length 0.)\n")
  invisible(x)
}
# Use full component names, not $est
# Print a label, if present.
# Add option to not print the call
# If not print call or label, then print defaultLabel


##########################################################
# print.summary.bootstrap
##########################################################
print.summary.bootstrap <-
function(x, digits = max(options()$digits - 3, 4), ...)
{
  old.digits <- options(digits = digits)
  on.exit(options(old.digits))
  cat("Call:\n")
  print(x$call)
  cat("\nNumber of Replications:", x$B, "\n")
  if(!is.null(x$B.missing))
    cat("\nReplicates with missing values:", x$B.missing,
	"(this may cause substantial bias).\n")
  cat("\nSummary Statistics:\n")
  print(cbind(Observed = x$observed, x$estimate))
  f <- function(Name, value, ...)
    if(!is.null(value)){
      cat("\n", sep="", Name, ":\n")
      print(value, ...)
    }
  f("Percentiles", x$limits.emp, ...)
  f("BCa Confidence Intervals", x$limits.bca, ...)
  f("Bootstrap Tilting Confidence Intervals", x$limits.tilt, ...)
  f("t limits with bootstrap standard error", x$limits.t, ...)
  if(length(x$correlation) > 1) {
    cat("\nCorrelation of Replicates:\n")
    print(x$correlation)
  }
  invisible(x)
}
# I changed "BCa Percentiles" to "BCa Confidence Limits"
# I added warning about bias.
# Use full component names, not $est
# Change "Empirical Percentiles" to "Percentiles"
# Print limits and correlations if present.


##########################################################
# print.summary.resamp
##########################################################
print.summary.resamp <-
function(x, digits = max(options()$digits - 3, 4), ...)
{
  old.digits <- options(digits = digits)
  on.exit(options(old.digits))
  print.resamp(x, ...)

  # Also print empirical percentiles, and correlations
  cat("\nPercentiles:\n")
  print(x$limits.emp)
  if(length(x$correlation) > 1) {
    cat("\nCorrelation of Replicates:\n")
    print(x$correlation)
  }
  invisible(x)
}
# Use full component names, not $est
# Call print.resamp (it had some things this didn't -- will be easier
#   to maintain this way).
# Change "Empirical Percentiles" to "Percentiles"
# Print correlations if present.

##########################################################
# qqnorm.resamp
##########################################################
qqnorm.resamp <- function(x, nrow = NULL, grid.layout = T, lines = T,
			  subset.statistic = 1:p, ..., main=NULL){
  # qqnorm plots of resample replicates
  p <- ncol(x$replicates)
  P <- length(subset.statistic)
  if(is.character(subset.statistic)){
    j <- match(subset.statistic, names(x$observed))
    if(anyMissing(j))
      stop("subset.statistic names not recognized:\n    ",
	   paste(subset.statistic[which.na(j)], collapse=", "))
    subset.statistic <- j
  }
  if(grid.layout && P != 1) {
    if(is.null(nrow)){
      if(length(x$dim.obs) == 2)
	old.par <- par(mfcol = x$dim.obs)
      else {
        nrow <- min(P, 2)
        old.par <- par(mfrow = c(nrow, ceiling(P/nrow)))
      }
    }
    else
      old.par <- par(mfrow = c(nrow, ceiling(P/nrow)))
    on.exit(par(old.par))
  }
  # Main title; if not supplied then use x$label or x$defaultLabel
  if(is.null(main))
    main <- x$label
  if(is.null(main))
    main <- x$defaultLabel

  for(i in subset.statistic) {
    xi <- x$replicates[, i]
    qqnorm(xi,
	   ylab = names(x$observed)[i], weights = x$weights, ...,
	   main=if(is.null(main)) names(x$observed)[i] else main)
    if(lines)
      qqline(xi, weights = x$weights)
    dots <- list(...)
    if(!is.element("axes", names(dots)) || dots$axes){
      if(length(dots)){
	do.call("box", dots[is.element(names(dots), names(par()))])
	# like box(...), but only pass par parameters like bty, col
      }
      else
	box()
    }
  }
  invisible(NULL)
}
# Avoid using old.par <- par(new values) trick.
# Use {} to clarify if() if() else() else()
# Support dim.obs=scalar
# Change nstats to p
# Skip NA handling -- qqline is fixed, and call to warning would
#    have crashed anyway due to lack of paste()
# stop if x has weights
# Only add title if main not supplied.
# support weights
# Oops, restore the old.par <- par(new values) trick (bug 24874)
# Added argument subset.statistic
# If axes=F, do not call box()
# Add title, change ylab from "Quantiles of replicates" to name of observed
# Use box(...) instead of box(), to support bty

##########################################################
# summary.bootstrap
##########################################################
summary.bootstrap <-
function(object, probs = c(25, 50, 950, 975)/1000,
	 frame.eval = object$parent.frame, narrow=F, ...,
	 summaryBCa = resampleOptions()$summaryBCa,
	 summaryTilt = resampleOptions()$summaryTilt,
	 summaryT = resampleOptions()$summaryT,
	 correlation = resampleOptions()$summaryCorrelations)
{
  p <- length(object$observed)
  if(!p)
    stop("statistic had length 0")
  if(is.null(frame.eval))
    frame.eval <- sys.parent(1)
  if(is.null(summaryBCa)) summaryBCa <- T
  if(is.null(summaryTilt)) summaryTilt <- !is.null(object$L)
  if(is.null(summaryT)) summaryT <- F
  if(is.null(correlation)) correlation <- (p < 6)
  if(p == 1) correlation <- F

  if(!is.null(resampCheckIfOrderMatters(object$order.matters))){
    # Cannot calculate tilting or BCa limits; disable without warning
    summaryBCa <- F
    summaryTilt <- F
  }

  out <- list(call = object$call,
	      B = object$B,
	      observed = object$observed,
	      estimate = object$estimate,
	      limits.emp = limits.percentile(object, probs = probs,
		narrow = narrow))
  if(summaryBCa) {
    out$limits.bca <- try(limits.bca(object, probs = probs,
				     frame.eval = frame.eval,
				     narrow=narrow, ...))
    if(is(out$limits.bca, "Error")){
      warning("Cannot compute BCa limits, will omit")
      out$limits.bca <- NULL
    }
  }

  if(summaryTilt){
    out$limits.tilt <- try(limits.tilt(object, probs = probs,
				       frame.eval = frame.eval, ...))
    if(is(out$limits.tilt, "Error")){
      warning("Cannot compute tilting limits, will omit")
      out$limits.tilt <- NULL
    }
    # limits.tilt normally returns a character string if it cannot compute
    if(is.character(out$limits.tilt))
      out$limits.tilt <- NULL
  }
  if(summaryT)
    out$limits.t <- limits.t(object, probs = probs,
			     frame.eval = frame.eval, ...)
  if(correlation)
    out$correlation = cor(object$replicates, na.method = "available",
      weights = object$weights)
  if(!is.null(object$B.missing))
    out$B.missing <- sum(rowSums(is.na(object$replicates)) > 0)
  oldClass(out) <- "summary.bootstrap"
  out
}
# Use full component names, not $est
# Stop of statistic had length 0 (to prevent later noninformative stop).
# Change default for frame.eval.jack
# Add argument narrow=F
# Add argument `frame.eval', replace internal use of `frame.eval.jack'.
# Add arguments determining whether to do limits and correlation
# Remove frame.eval.jack
# Skip BCa and tilting if resampling residuals.
# Use try() to prevent stopping if limits.tilt or limits.bca fail.

# I've made a change to the Splus6 version on 1/23/2001 that should be
# backed out in favor of this version, once limits.bca handles stratification.

##########################################################
# summary.resamp
##########################################################
summary.resamp <-
function(object, probs = c(25, 50, 950, 975)/1000, narrow=F, ...,
	 correlation = resampleOptions()$summaryCorrelations)
	 
{
  # Summary method for "resamp" objects.
  out <- list(call = object$call,
	      B = object$B,
	      observed = object$observed,
	      estimate = object$estimate,
	      limits.emp = limits.percentile(object, probs = probs,
		narrow = narrow))
  p <- length(object$observed)
  if(is.null(correlation)) correlation <- (p < 6)
  if(p == 1) correlation <- F
  if(correlation)
    out$correlation <- cor(object$replicates, weights = object$weights)
  oldClass(out) <- "summary.resamp"
  out
}
# Use full component names, not $est
# Add narrow=F argument
# Added ... (avoid bug 20501)
# Add arguments determining whether to do correlation


##################################################
# pairs.resamp
##################################################
pairs.resamp <- function(x, ..., subset.statistic = T){
  # Do a pairs plot of x$replicates
  # By default include all columns.
  # May use subset.statistic to select only certain columns.
  x <- x$replicates[,subset.statistic]
  pairs(x, ...)
}

##################################################
# boxplot.resamp
##################################################
boxplot.resamp <- function(x, ...,
			  subset.statistic = 1:p,
			  nrow = NULL, grid.layout = T,
			  main = NULL){
  # Do a boxplot of  x$replicates
  # May use subset.statistic to select only certain columns.
  # grid.layout is F if each column should be shown in a separate plot
  p <- length(x$observed)
  P <- length(subset.statistic)
  if(is.character(subset.statistic)) {
    j <- match(subset.statistic, names(x$observed))
    if(anyMissing(j))
      stop("subset.statistic names not recognized:\n    ", paste(
        subset.statistic[which.na(j)], collapse = ", "))
    subset.statistic <- j
  }
  if(grid.layout && P != 1) {
    if(is.null(nrow))
      nrow <- floor(sqrt(P*2/3)) # 1 for P < 6
    old.par <- par(mfrow = c(nrow, ceiling(P/nrow)))
    on.exit(par(old.par))
  }
  # Main title; if not supplied then use x$label or x$defaultLabel
  if(is.null(main)) main <- x$label
  if(is.null(main))
    main <- x$defaultLabel
  for(i in subset.statistic)
    boxplot(x$replicates[, i], ..., main = main,
	    ylab = names(x$observed)[i])
  invisible(NULL)
}
# Note that nrow has a different default than for plot.resamp
# Here default is 1 row unless 6 or more columns.


resampCheckIfOrderMatters <- function(order.matters){
  # order.matters is NULL, logical, or character.
  # If NULL or FALSE then order doesn't matter, and return NULL
  # Else return character, a message that can be appended to
  # a stop or warning message
  if(is.null(order.matters) || identical(order.matters, FALSE))
    return(NULL)
  if(is.logical(order.matters))
    return("order of observations matters")
  order.matters
}
