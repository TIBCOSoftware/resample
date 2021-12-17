# Copyright 2021. TIBCO Software Inc.
# This file is subject to the license terms contained
# in the license file that is distributed with this file.

# Here are functions:
#	permutationTest
#	print.permutationTest

# setOldClass(c("permutationTest", "resamp"))

permutationTest <- function(data, statistic, B = 999, ...,
			   alternative = "two.sided",
			   combine = NULL,
			   combinationFunction = combinePValues.Fisher
			   ){

  if(missing(statistic))
    stop("Argument statistic is missing with no default")

  alternative <- match.arg(alternative, c("two.sided", "greater", "less"))

  # Use `bootstrap' to do the sampling
  new.call <- match.call()
  if(any(w <- is.element(names(new.call),
			c("sampler", "sampler.prob"))))
    stop(paste("arguments",paste(names(new.call)[w], collapse=" "),
	       "are not supported in permutationTest"))
  new.call[[1]] <- as.name("bootstrap")
  new.call$sampler <- as.name("samp.permute")
  new.call$alternative <- NULL
  new.call$combine <- NULL
  new.call$combinationFunction <- NULL
  new.call$B <- B
  if(is.null(new.call$L))
    new.call$L <- "none"
  if(is.null(new.call$save.indices))
    new.call$save.indices <- F
  obj <- eval(new.call, sys.parent(1))
  obj$call <- match.call()
  obj$estimate$Bias <- NULL
  obj$defaultLabel <- resampMakeLabel("permutation",
				     new.call$data, new.call$statistic)
  if(anyMissing(obj$replicates)){
    oldClass(obj) <- NULL
    on.exit(assign(".permutationTest.partial.results", obj, immediate=T))
    stop("Missing values found in the replicates, partial results saved in .permutationTest.partial.results")
  }
  p <- length(obj$observed)
  alternative <- rep(alternative, length=p)

  if(length(combine) && p == 1){
    warning("No combinations are possible with only one statistic")
    combine <- NULL
  }
  if(length(combine)){
    Names <- names(obj$observed)
    if(!is.list(combine))
      combine <- list(combine)
    if(any(w <- sapply(combine, is.character)))
      combine[w] <- lapply(combine[w], match, table=Names)
    temp <- findPvalues(obj$observed, obj$replicates,
		       B, p, alternative, combine, combinationFunction)
    Pvalues <- temp$Pvalues
    names(temp$combP) <- sapply(combine, function(x, Names) {
	if(!is.character(x)) x <- Names[x]
	paste(x, collapse=", ")
    }, Names=Names)
    obj$"combined-p-value" <- temp$combP
  }
  else
    Pvalues <- findPvalues(obj$observed, obj$replicates,
			  B, p, alternative, NULL, NULL)

  obj$estimate <- cbind(obj$estimate,
		       alternative = alternative,
		       "p-value" = Pvalues)
  oldClass(obj) <- c("permutationTest", "resamp")
  obj
}
# 6/26/00: removed confidence intervals.
# Set L="none" if not specified.

# To do:
#   allow missing values?
#   support probabilities & optional samplers?
#   simplify the process of specifying a two-sample or multi-sample
#     statistic (i.e. a variable that is not subject to subscripting)

print.permutationTest <- function(x, ...){
  # Print a permutationTest object.
  # First use print.resamp (which prints p-values as part of the estimate
  # component).  Then if there is a combined p-value, print it too.
  NextMethod("print")
  if(length(x$"combined-p-value")){
    cat("\n")
    print(data.frame("Combined p-value:" = x$"combined-p-value",
		     check.names=F))
  }
  invisible(x)
}
