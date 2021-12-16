# file contains one function:
# parametricBootstrap: used for parametric bootstrapping  (modified for Sv4)

# setOldClass(c("parametricBootstrap","resamp")) # now defined in meta/

parametricBootstrap <- function(data, statistic, rsampler, B = 1000, args.stat = NULL,
		       args.rsampler = NULL, dsampler = NULL,
		       seed = .Random.seed,
		       label = NULL,
		       statisticNames = NULL,
		       block.size = min(100, B),
		       trace = resampleOptions()$trace,
                       assign.frame1 = F, save.samples = F,
		       statistic.is.random, seed.statistic = 500){

  if(missing(statistic))
    stop("Argument statistic is missing with no default")

  if(missing(rsampler))
    stop("Argument rsampler is missing with no default")

  # Capture call.
  func.call <- match.call()

  # If statistic isn't function or name of function or expression,
  # store it as an expression to pass to fit.func.
  substitute.stat <- substitute(statistic)
  if(is.na(match(mode(substitute.stat), c("name", "function"))))
    statistic <- substitute.stat

  # Get name of data.
  data.name <- substitute(data)
  data.name <- if (!is.name(data.name)) "data"
               else as.character(data.name)

  # Handle assign.frame1; prevent overwriting object on frame 1
  if(assign.frame1){
    if(exists(data.name, frame = 1))
      stop("You specified assign.frame1, but the name of your data matches an object that already exists on frame 1; to prevent this object from being overwritten, assign it to a permanent database before calling parametricBootstrap, and remove it from frame 1")
    on.exit(if(exists(data.name, frame = 1))
	    remove(data.name, frame = 1))
  }

  # Create function to evaluate the statistic given data.
  is.df.data <- is.data.frame(data)
  fit.func <- resampMakeFuncSimple(statistic = statistic,
				   data.name = data.name,
				   is.df.data = is.df.data,
				   is.null.args.stat = is.null(args.stat),
				   assign.frame1 = assign.frame1)
  # is: function(data, statistic, args.stat)

  # Faster subscripting for a data frame, allow duplicate row names.
  if(is.df.data && is.null(attr(data, "dup.row.names")))
    attr(data, "dup.row.names") <- T

  # Set seed in case statistic uses randomization.
  seed <- eval(seed)
  if(missing(statistic.is.random) || statistic.is.random) {
    set.seed(seed.statistic)
    prev.seed <- .Random.seed
  }

  #  Get parameter values for observed data.
  n <- numRows(data)

  observed <- fit.func(data = data, statistic = statistic,
		       args.stat = args.stat)
  # Determine if statistic uses randomization; this may fail if
  # a statistic sometimes use randomization.
  if(missing(statistic.is.random))
    statistic.is.random <- any(.Random.seed != prev.seed)
  if(statistic.is.random) seed.statistic <- .Random.seed

  # Check that observed is vector or matrix.
  if(is.null(observed))
    stop("Statistic returned a NULL result on observed data.  It must return a vector or matrix.")
  if(!is.atomic(observed)) stop("Statistic must return a vector or matrix.")

  # Getting parameter names and coercing matrix to vector.
  names.observed <- resampMakeNames(observed = observed,
				    statistic = substitute.stat,
				    statisticNames = statisticNames)
  dim.obs <- dim(observed)
  if(!is.null(dim.obs)) observed <- as.vector(observed)
  names(observed) <- names.observed

  set.seed(seed)
  seed.start <- .Random.seed

  must.swap <- statistic.is.random & any(.Random.seed != seed.statistic)
  # Need to swap only if both the sampler and statistic use .Random.seed

  reps <- matrix(NA, nrow = length(observed), ncol = sum(B))
  if(save.samples) all.samples <- vector("list", sum(B))
  nB <- length(B)

  do.call.noeval <- function(what, args = list()){
    if(!is.list(args)) stop("args should be a list")
    this.call <- c(list(as.name(what)), args)
    mode(this.call) <- "call"
    this.call
  }
  samplerCall <- do.call.noeval("rsampler",
                                c(list(n), args.rsampler))
  # Below, assign .parametricBootstrap.partial.results, either with a
  # dummy list if no results were generated before the interrupt
  # (so that smoothedBootstrap can know that no results were generated),
  # or with the actual partial results.

  on.exit({
    if(!totalB){
      assign(".parametricBootstrap.partial.results", where = 1,
	     immediate = T, list(B = 0))
    }
    if(totalB) {
      cat("\nDid ", totalB,
" replications, saving results in .parametricBootstrap.partial.results, interrupt again to abort completely.\n")
      reps <- t(reps[, 1:totalB, drop = F])
      if(save.samples) all.samples <- all.samples[1:totalB]
      dimnames(reps) <- list(NULL, names.observed)
      func.call$B <- c(B[1:iB][ - iB], doneB)
      seed.end <- "Unknown, due to interrupt"
      assign(".parametricBootstrap.partial.results", where = 1,
	     immediate = T,
	     makeParametricBootstrap(replicates = reps,
			observed =  observed,
			n = n,
			call = func.call,
			seed.start = seed.start,
			seed.end = seed.end,
			dim.obs = dim.obs,
			parent.frame = sys.parent(),
			samples=if(save.samples) all.samples,
			rsampler = rsampler,
			args.rsampler = args.rsampler,
			dsampler = dsampler,
			defaultLabel = resampMakeLabel("parametric bootstrap",
			  func.call$data, func.call$statistic),
			label = label))
    }
  } , add = T)

  for(iB in seq(nB)) {
    nblocks <- ceiling(B[iB]/block.size)
    B2 <- block.size
    previousB <- sum(B[1:iB]) - B[iB]

    for(i in 1:nblocks) {
      doneB <- (i - 1) * block.size
      totalB <- previousB + doneB
      if(i == nblocks)
	if(B[iB] %% block.size)	B2 <- B[iB] %% block.size
      if(trace)
	cat("Forming replications ", totalB + 1, " to ",
	    totalB + B2, "\n")

      if(must.swap)
	for(ii in 1:B2){
	  tempData <- eval(samplerCall)
	  seed.sampler <- .Random.seed
	  .Random.seed <<- seed.statistic
	  reps[,totalB + ii] <-
	    fit.func(data = tempData, statistic = statistic,
		     args.stat = args.stat)
	  seed.statistic <- .Random.seed
	  .Random.seed <<- seed.sampler
	  if(save.samples==T)
	    all.samples[[totalB+ii]] <- tempData
	}
      else{
	if(save.samples == T)
	  for(ii in 1:B2){
	    tempData <- eval(samplerCall)
	    reps[,totalB + ii] <-
	      fit.func(data = tempData, statistic = statistic,
		       args.stat = args.stat)
	    all.samples[[totalB+ii]] <- tempData
	  }
	else
	  for(ii in 1:B2)
	    reps[,totalB + ii] <-
	      fit.func(data = eval(samplerCall), statistic = statistic,
		       args.stat = args.stat)
      }
    }
  }
  reps <- t(reps)
  dimnames(reps) <- list(NULL, names.observed)
  seed.end <- .Random.seed

  if(assign.frame1 && exists(data.name, frame = 1))
    remove(data.name, frame = 1)
  on.exit(add=F)

  if(trace)  cat("\n")

  makeParametricBootstrap(replicates = reps,
	     observed = observed,
	     n = n,
	     call = func.call,
	     seed.start = seed.start,
	     seed.end = seed.end,
	     dim.obs = dim.obs,
	     parent.frame = sys.parent(),
	     samples = if(save.samples) all.samples,
	     rsampler = rsampler,
	     args.rsampler = args.rsampler,
	     dsampler = dsampler,
	     defaultLabel = resampMakeLabel("influence",
	       func.call$data, func.call$statistic),
	     label = label)
}
# remove argument passStatistic (just use args.rsampler instead)
# add arguments label and statisticNames, output label and defaultLabel
