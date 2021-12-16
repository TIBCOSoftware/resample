##########################################################
# controlVariates (S3 generic)
##########################################################

controlVariates <- function(x, ...)
  UseMethod("controlVariates")



##########################################################
# controlVariates.default
##########################################################

controlVariates.default <- function(x, mu, method="linear", positive=T,
			   weights=NULL){
  # x is a vector, or [n, p] matrix
  # Return weights such that
  #	colSums(x, weights) = mu
  # (unless positive=T, in which case weights are non-negative
  # and that equality may not hold exactly).
  # If x and y are correlated, then
  #     mean(y, weights=controlVariates(x,mu))
  # is the estimated mean of y, controlled for x.
  #
  # If input weights are supplied, the output weights are the
  # input weights times a function of x.
  #
  # For binary data, where x takes on values 0 or 1, all methods
  # yield equivalent results; the "binary" method is a shortcut.

  method <- match.arg(casefold(method), c("linear", "exponential", "ml", "binary"))
  if(is.data.frame(x))
    x <- numericalMatrix(x)
  p <- length(mu)
  if(p != numCols(x))
    stop("length of mu is not equal to number of columns of x")
  n <- numRows(x)
  if(length(weights)){
    if(length(weights) != n)
      stop("length of weights does not match number of rows of x")
    weights <- weights / sum(weights)
  }
  if(method == "binary"){
    if(p > 1)
      stop("x must be univariate for the binary method")
    if(is.logical(x))
      one <- x
    else
      one <- as.logical(x)
    if(mu <= 0 || mu >= 1 || (!is.logical(x) && any(x != one))){
      warning("conditions not met for binary method; switching to linear method")
      method <- "linear"
      positive <- T
    }
    else {
      # true binary data
      if(length(weights)){
	sum0 <- sum(weights[!one])
	sum1 <- sum(weights[one])
	weights <- ifelse(one, mu/sum1, (1-mu)/sum0) * weights
      }
      else {
	n1 <- sum(one)
	weights <- ifelse(one, mu/n1, (1-mu)/(n-n1))
      }
    }
  }
  if(is.element(method, c("exponential", "ml"))){
    if(p > 1)
      stop("only univariate exponential or ml tilting is currently available")
    if(method == "ml")
      tau <- tiltMeanSolve(q = mu, L = x, weights=weights, tilt="ml")$tau.ml
    else
      tau <- tiltMeanSolve(q = mu, L = x, weights=weights)$tau.exp
    if(is.na(tau)){
      warning("mu is outside the range of x; switching to linear method")
      method <- "linear"
      positive <- T
    }
    else
      weights <- tiltWeights(tau = tau, L = x,
			    tilt = method, weights = weights)
  }
  if(method == "linear"){
    meanx <- colMeans(x, weights=weights)
    fac <- try(solve(var(x, weights=weights, unbiased=F), mu - meanx))
    if(is(fac, "Error")){
      warning("controlVariates failed, with message: ", fac[[1]],
              "\nWill try a workaround.")
      fac <- lm.fit.svd(var(x, weights=weights, unbiased=F), mu-meanx)$coef
    }
    old <- (if(length(weights)) weights else 1/n)
    if(is.matrix(x))
      weights <- old * (1 + c((x - rep(meanx, each=n)) %*% fac))
    else
      weights <- old * (1 + fac * (x - meanx))
    if(positive && any(weights < 0)){
      warning("some weights were negative")
      weights <- pmax(weights, 0)
      weights <- weights/sum(weights)
    }
  }
  weights
}
# Workaround for singular x

# To do:
#	support multivariate cv's for ML and exp methods, once
#		tiltMeanSolve supports that.



##########################################################
# controlVariates.bootstrap
##########################################################

controlVariates.bootstrap <- 
function(x,
	 subset.covariates = 1:p,
         moments = 2,
         ...,
	 L = resampGetL(x, frame.eval = frame.eval),
	 group = resampGetArgument(x, "group",
           frame.eval = frame.eval),
	 treatment = resampGetArgument(x, "treatment",
	   frame.eval = frame.eval),
	 frame.eval = x$parent.frame){

  # x:                  bootstrap object
  # subset.covariates: can be used to specify that only some columns of
  #	L should be used for adjustments.
  # moments: scalar integer, 1 = use mean of L, 2 = mean & variance
  # ... = additional arguments to pass to controlVariates.default,
  #       in particular method and positive.  Should not include weights.
  #       Currently only method="linear" is supported with moments=2,
  #       or if subset.covariates has length > 1

  if(is(x, "concomitants"))
    stop("Cannot do concomitants adjustments on a bootstrap object that has already been adjusted")
  if(is.null(frame.eval))
    frame.eval <- sys.parent(1)
  n <- x$n
  B <- x$B
  p <- length(x$observed)

  # If both treatment and group, stratify on both;
  # if only treatment, then call it "group".
  if(length(treatment)){
    if(length(group))
      group <- paste(match(group, unique(group)),
		    match(treatment, unique(treatment)))
    else
      group <- treatment
  }
  # Note - for bootstrap2 objects, this function calculates weights
  # for the main object, not for each of the two contained bootstrap objects.


  if(missing(L) && !is.null(x$Lstar)){
    if(length(group) && moments >= 2)
      indices <- resampGetIndices(x, frame.eval = frame.eval)
    Lstar <- x$Lstar
  }
  else {
    indices <- resampGetIndices(x, frame.eval = frame.eval)
    Lstar <- indexMeans(L, indices, group = group)
  }
  if(length(subset.covariates) < p){
    L <- L[,subset.covariates, drop=F]
    Lstar <- Lstar[,subset.covariates, drop=F]
  }

  # Calculate expected values for covariates
  # In case of moments = 2, use both means and variances.
  X <- Lstar
  if(!length(group)){
    mu <- colMeans(L)
    if(moments >= 2){
      X <- cbind(Lstar, (Lstar - rep(mu, each=B))^2)
      mu <- c(mu, colVars(L, unbiased=F) / n)
    }
  }
  else { # Lstar is the sum of group means
    mu <- colSums(groupMeans(L, group))
    if(moments >= 2){
      X <- cbind(Lstar, (Lstar - rep(mu, each=B))^2)
      gp <- as.numeric(factor(group))
      nGroups <- max(gp)
      groupSizes <- tabulate(gp[indices[,1]], nGroups)
      # Used gp[indices[,1]] rather than gp, in case did sampling with
      # sizes different than the original group sizes.
      # We assume here that all columns of indices have the same groupSizes.
      gpVars <- groupVars(L, group=group, unbiased = F)
      mu <- c(mu, colSums( gpVars / groupSizes))
    }
  }
  weights <- controlVariates.default(x = X, mu = mu, weights = x$weights, ...)
  means <- colMeans(x$replicates, na.rm = T, weights = weights)
  x$original <- x[c("call", "estimate")]
  x$call <- match.call()
  x$estimate <- 
    data.frame(Mean = means,
	       Bias = means - x$observed,
	       SE = colStdevs(x$replicates, na.rm=T, weights=weights))
  x$weights <- weights
  oldClass(x) <- c("controlVariates.bootstrap", "bootstrap", "resamp")
  x
}
# Allow weights in the bootstrap object


##########################################################
# methods
##########################################################

print.controlVariates.bootstrap <-
function(x, ...)
{
  cat("Original bootstrap call:\n")
  print(x$original$call)
  cat("\n")
  NextMethod("print")
  invisible(x)
}

update.controlVariates.bootstrap <- function(object, ...){
  stop("A controlVariates.bootstrap object cannot be updated; update the original bootstrap object")
}
