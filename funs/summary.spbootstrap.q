# Copyright 2021. TIBCO Software Inc.
# This file is subject to the license terms contained
# in the license file that is distributed with this file.

# This file contains
#	summary.parametricBootstrap
#	summary.smoothedBootstrap
#	print.summary.parametricBootstrap
#	print.summary.smoothedBootstrap (identical to previous)

# The functions below provide a summary of a parametricBootstrap or smoothedBootstrap
# object, respectively.


##########################################################
# summary.parametricBootstrap
##########################################################

summary.parametricBootstrap <- function(object, probs = c(25, 50, 950, 975)/1000, ...){
  if(!length(object$observed)) stop("statistic had length 0")
  out <- list(call = object$call,
              B = object$B,
              observed = object$observed,
              estimate = object$estimate,
              limits.emp = limits.percentile(object, probs),
              correlation = cor(object$replicates, na.method = "available"))
  if(!is.null(object$B.missing))
    out$B.missing <- sum(rowSums(is.na(object$replicates)) > 0)
  oldClass(out) <- "summary.parametricBootstrap"
  out
}

##########################################################
# summary.smoothedBootstrap
##########################################################

summary.smoothedBootstrap <- function(object, probs = c(25, 50, 950, 975)/1000, ...){
  if(!length(object$observed)) stop("statistic had length 0")
  out <- list(call = object$call,
              B = object$B,
              observed = object$observed,
              estimate = object$estimate,
              limits.emp = limits.percentile(object, probs),
              correlation = cor(object$replicates, na.method = "available"))
  if(!is.null(object$B.missing))
    out$B.missing <- sum(rowSums(is.na(object$replicates)) > 0)
  oldClass(out) <- "summary.smoothedBootstrap"
  out
}

# note: These functions will differ more once a BCa interval (or other
#       second-order-correct interval) is included in the procedure.



##########################################################
# print.summary.parametricBootstrap
##########################################################
print.summary.parametricBootstrap <- function(x, digits = max(options()$digits - 3, 4), 
				     ...){

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
  cat("\nPercentiles:\n")
  print(x$limits.emp)
  if(length(x$correlation) > 1) {
    cat("\nCorrelation of Replicates:\n")
    print(x$correlation)
  }
  invisible(x)
}

print.summary.smoothedBootstrap <- print.summary.parametricBootstrap

