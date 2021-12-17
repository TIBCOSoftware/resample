# Copyright 2021. TIBCO Software Inc.
# This file is subject to the license terms contained
# in the license file that is distributed with this file.

# file contains one function:
# parametricBootstrap: used for parametric bootstrapping (for S version 4)


# setOldClass(c("parametricBootstrapTest","resamp"))  # defined in meta/

print.parametricBootstrapTest <- function(x, digits = max(options()$digits - 3, 4), ...,
			    printCall = .Options.resample$printCall){

  old.digits <- options(digits = digits)
  on.exit(options(old.digits))

  # Print label, and maybe call or defaultLabel
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
  cat("\nIn the following table, each row contains the summary statistics\nfor the given test statistic:\n")
  est <- x$estimate
  if(is.numeric(x$estimate$null.value) &&
     !is.null(names(x$estimate$null.value)))
    est$null.value <- paste(names(x$estimate$null.value),
			   x$estimate$null.value, sep = " = ")
  # special case: named numeric null.value -- print as string
  print(cbind(Observed = x$observed, est))
  invisible(x)
}


parametricBootstrapTest <-
  function(data, statistic, rsampler, B=999, args.stat = NULL,
	   args.rsampler = NULL, null.value = NULL, alternative="two.sided",
	   label = NULL, ...){

  if(missing(statistic))
    stop("Argument statistic is missing with no default")

  if(missing(rsampler))
    stop("Argument rsampler is missing with no default")

  alternative <- sapply(alternative, match.arg,
			c("two.sided", "greater", "less"))

  func.call <- match.call()       # Capture call.
  new.call <- sys.call()
  if(any(is.element(names(new.call), "passStatistic")))
    stop("argument passStatistic is not supported in parametricBootstrapTest")
  new.call$alternative <- new.call$null.value <- NULL
  new.call$B <- B
  new.call[[1]] <- as.name("parametricBootstrap")

  postProcess <- function(obj,alternative, null.value){
    p <- length(obj$observed)
    if (!is.null(null.value) & length(null.value) != p)
      stop("null.value did not match statistic length")
    B <- numRows(obj$replicates)
    alternative <- rep(alternative, length = p)
    Pvalues <- rep(NA, p)

    for (i in 1:p){
      X <- obj$replicates[,i]
      Q <- obj$observed[i]
      counts <- switch(alternative[i],
		       less = sum(X <= Q) + 1,
		       greater = sum(X >= Q) + 1,
		       two.sided= (min(sum(X <= Q), sum(X >= Q))+1)*2)
      Pvalues[i] <- counts / (sum(B)+1)
    }
    rnames <- row.names(obj$estimate) #save names from makeParametricBootstrap()
    if (!is.null(null.value)){

      # use I() twice below to avoid conversion to factor
      if(!is.numeric(null.value))
	obj$estimate <- data.frame(null.value = I(null.value),
				   alternative= I(alternative),
				   p.value = Pvalues)
      else{
	# do the assigning in 2 steps to preserved
	# possible names for named numeric null.value
	# don't need to use I(null.value) in this case
	obj$estimate <- data.frame(null.value =
				   rep(0, length(null.value)),
				   alternative=I(alternative),
				   p.value = Pvalues)
	obj$estimate$null.value <- null.value
      }
    }
    else
      obj$estimate <- data.frame(alternative=I(alternative),
				 p.value = Pvalues)
    row.names(obj$estimate) <- rnames
    obj
  }

  on.exit({
    pbootB <- .parametricBootstrap.partial.results$B[1]
    if(pbootB){
      cat("\nTransferring results from \n.parametricBootstrap.partial.results to .parametricBootstrapTest.partial.results\n")
      partRes <- .parametricBootstrap.partial.results
      partRes$call <- func.call
      if(anyMissing(partRes$replicates))
	warning("missing values in replicates are causing missing p-values")

      partRes <- postProcess(obj = partRes, alternative = alternative,
			     null.value = null.value)
      oldClass(partRes) <- c("parametricBootstrapTest", "resamp")
      remove(".parametricBootstrap.partial.results",where = 1)
      assign(".parametricBootstrapTest.partial.results",where = 1,
	     immediate = T, partRes)
    }
  }
	  , add = T)

  obj <- eval(new.call, sys.parent())
  obj$call <- func.call
  if(anyMissing(obj$replicates))
    warning("missing values in replicates are causing missing p-values")
  obj <- postProcess(obj = obj, alternative = alternative,
		     null.value = null.value)
  on.exit()
  if(!is.null(label))
    obj$label <- label
  obj$defaultLabel <- resampMakeLabel("parametricBootstrapTest",
				      func.call$data, func.call$statistic)
  oldClass(obj) <- c("parametricBootstrapTest", "resamp")
  obj
}
