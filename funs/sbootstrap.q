# file contains one function:
# smoothedBootstrap: used for smoothed bootstrapping (modified for Sv4)

# setOldClass(c("smoothedBootstrap","resamp"))  # now defined in meta/

smoothedBootstrap <- function(data, statistic, B = 1000, args.stat = NULL,
		       sampler = samp.bootstrap, seed = .Random.seed,
		       smoother = rmvnorm,
		       args.smoother =
		       list(cov = as.matrix(var(data))/numRows(data),
			    d = numCols(data)),
		       label = NULL,
		       statisticNames = NULL,
		       block.size = min(100, B),
                       trace = resampleOptions()$trace, assign.frame1 = F,
		       save.samples = F, statistic.is.random,
		       seed.statistic = 500){

  if(missing(statistic))
    stop("Argument statistic is missing with no default")

  func.call <- match.call() 		# Capture call.

  # Create function to be passed as the rsampler argument to the
  # function parametricBootstrap().
  # (See resampGetFitFunc for analogous function for bootstrap.)
  dimLength <- length(dim(data))
  sss <- paste("function(n, data, smoother, args.smoother, sampler){",
	       "data[sampler(n,1)",
	       if(dimLength > 1) paste(c(rep(",", dimLength), "drop=F"),
				       collapse = ""),
	       "] + do.call('smoother', c(n = n, args.smoother))}")
  ssampler <- eval(parse(text = sss))

  # Build parametricBootstrap() function call.
  new.call <- func.call
  new.call$smoother <- new.call$sampler <- new.call$args.smoother <- NULL
  new.call$rsampler <- ssampler
  new.call$args.rsampler <- list(data = data, smoother = smoother,
				 args.smoother = args.smoother,
				 sampler = sampler)
  new.call[[1]] <- as.name("parametricBootstrap")

  # On early termination, convert any partial results from parametricBootstrap
  # into smoothedBootstrap format.
  on.exit({
    pbootB <- .parametricBootstrap.partial.results$B[1]
    if(pbootB){
      cat("\nTransferring results from \n.parametricBootstrap.partial.results to .smoothedBootstrap.partial.results\n")
      partRes <- .parametricBootstrap.partial.results
      partRes$rsampler <-  partRes$args.rsampler <-
	partRes$dsampler <- NULL
      partRes$call <- func.call
      oldClass(partRes) <- c("smoothedBootstrap", "resamp")
      remove(".parametricBootstrap.partial.results", where = 1)
      assign(".smoothedBootstrap.partial.results", where = 1,
	     immediate = T, partRes)
    }
  }, add = T)

  result <- eval(new.call, sys.parent())
  result$rsampler <- result$args.rsampler <- result$dsampler <- NULL
  result$call <- func.call
  result$defaultLabel <- resampMakeLabel("smoothed bootstrap",
					 func.call$data, func.call$statistic)
  oldClass(result) <- c("smoothedBootstrap", "resamp")
  on.exit()
  result
}



