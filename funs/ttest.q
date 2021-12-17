# Copyright 2021. TIBCO Software Inc.
# This file is subject to the license terms contained
# in the license file that is distributed with this file.

# Modified version of t-test -- add tilting intervals and p-value
# (but not for two-sample tests).  Exponential tilting only, implement with
# saddlepoint.

# This version adds an attribute to the htest object produced
# by the original t.test, with components "method", "p.value", and
# "conf.int", and has class htest.saddlepoint  There is a print
# method for this class.
saddlepoint.test <-
function(x, y = NULL, alternative = "two.sided", mu = 0, paired = FALSE,
	 var.equal = FALSE, conf.level = 0.95, treatment = NULL)
{
  alt.expanded <- if(!missing(alternative))
    char.expand(casefold(alternative),
                c("two.sided", "greater", "less"),
                stop("argument 'alternative' must match one of \"greater\", \"less\", \"two.sided\".")) else alternative
  if(!missing(mu))
    if((length(mu) != 1) || !is.finite(mu))
      stop("argument 'mu' must be a single finite numeric value.")
  if(!missing(conf.level))
    if((length(conf.level) != 1) || !is.finite(conf.level) || (conf.level <= 0)
       || (conf.level >= 1))
      stop("argument 'conf.level' must be a single number greater than zero and less than one.")
  alpha <- 1 - conf.level
  xname <- deparse(substitute(x))
  if(!is.null(treatment)){
    utreatment <- unique(treatment)
    if(length(utreatment) > 2)
      stop("Only two levels in treatment are allowed")
    y <- x[treatment == utreatment[2]]
    x <- x[treatment == utreatment[1]]
    data.name <- paste(xname, ":",
                       deparse(substitute(treatment)), "=",
                       as.character(utreatment[1]), "or",
                       as.character(utreatment[2]))
  } else {
    if(!is.null(y))
      data.name <- paste(xname, "and", deparse(substitute(y)))
  }
  if(is.null(y)) {
    # one-sample t-test.
    if(paired) stop("argument 'y' missing for paired test.")
    if(!missing(var.equal))
      warning("argument 'var.equal' ignored for one-sample test.")
    if((bad.obs <- sum(!(x.ok <- is.finite(x)))) > 0) {
      is.not.finite.warning(x)
      x <- x[x.ok]
      warning(paste(bad.obs, "observations with NA/NaN/Inf in 'x' removed."))
    }
    nx <- length(x)
    if(all(x == x[1]))
      warning("All values in 'x' are the same. Variance zero.")
    conf.int.xbar <- mean(x)
    conf.int.s <- sqrt(var(x)/nx)
    ret.val <- list(statistic = (conf.int.xbar - mu)/conf.int.s,
		    parameters = nx - 1,
		    estimate = conf.int.xbar,
		    null.value = mu,
		    alternative = alt.expanded,
		    method = "One-sample t-Test",
		    data.name = xname)
    names(ret.val$estimate) <- "mean of x"
    names(ret.val$null.value) <- "mean"
  }
  else if(paired) {
    # paired t-test
    if(!missing(var.equal))
      warning("argument 'var.equal' ignored for paired test.")
    if((nd <- length(x)) != length(y))
      stop("'x' and 'y' must have the same length when paired=TRUE.")
    d <- x - y
    if((bad.obs <- sum(!(both.ok <- is.finite(d)))) > 0) {
      if(!all(is.finite(x)))
        is.not.finite.warning(x)
      if(!all(is.finite(y)))
        is.not.finite.warning(y)
      d <- d[both.ok]
      warning(paste(bad.obs,
        "observations with NA/NaN/Inf in 'x' or 'y' removed."))
      nd <- length(d)
    }
    if(all(d == d[1]))
      warning("All values in 'x - y' are the same. Variance zero.")
    conf.int.xbar <- mean(d)
    conf.int.s <- sqrt(var(d)/nd)
    ret.val <- list(statistic = (conf.int.xbar - mu)/conf.int.s,
		    parameters = nd - 1,
		    estimate = conf.int.xbar,
		    null.value = mu,
		    alternative = alt.expanded,
		    method = "Paired t-Test",
		    data.name = data.name)
    names(ret.val$estimate) <- "mean of x - y"
    names(ret.val$null.value) <- "mean of differences"
  }
  else {
    # two-sample test
    if((bad.obs <- sum(!(x.ok <- is.finite(x)))) > 0) {
      is.not.finite.warning(x)
      x <- x[x.ok]
      warning(paste(bad.obs, "observations with NA/NaN/Inf in 'x' removed."))
    }
    nx <- length(x)
    if(all(x == x[1]) && nx > 1)
      warning("All values in 'x' are the same. Variance zero.")
    # Could omit that warning in the pooled case.  However, it
    # does indicate a violation of the normality assumption.
    if((bad.obs <- sum(!(y.ok <- is.finite(y)))) > 0) {
      is.not.finite.warning(y)
      y <- y[y.ok]
      warning(paste(bad.obs, "observations with NA/NaN/Inf in 'y' removed."))
    }
    ny <- length(y)
    if(all(y == y[1]) && ny > 1)
      warning("All values in 'y' are the same. Variance zero.")
    mean.x <- mean(x)
    mean.y <- mean(y)
    conf.int.xbar <- mean.x - mean.y
    if(!var.equal){
      var.x <- var(x)
      var.y <- var(y)
      conf.int.s <- sqrt((var.x/nx) + (var.y/ny))
    }
    else
      conf.int.s <- sqrt((1/nx + 1/ny) *
			 (sum((x-mean.x)^2) + sum((y-mean.y)^2)) /
                         (nx + ny - 2))
    ret.val <- c(if(var.equal)
                 list(method = "Pooled-Variance Two-Sample t-Test",
                      parameters = nx + ny - 2) else
                 list(method = "Welch Modified Two-Sample t-Test",
                      parameters = {
                        const <- 1/(1 + (nx * var.y)/(ny * var.x))
                        1/((const^2)/(nx - 1) + ((1 - const)^2)/(ny - 1))
                      }),
                 list(statistic = (conf.int.xbar - mu)/conf.int.s,
                      estimate = c(mean.x, mean.y),
                      null.value = mu,
                      alternative = alt.expanded,
                      data.name = data.name))
    names(ret.val$estimate) <- c("mean of x", "mean of y")
    names(ret.val$null.value) <- "difference in means"
  }
  ret.val <- c(ret.val,
	       switch(alt.expanded,
    two.sided = {
      conf.int.hw <- qt((1 - alpha/2), ret.val$parameters) * conf.int.s
      list(p.value = 2 * pt( - abs(ret.val$statistic), ret.val$parameters),
	   conf.int = c(conf.int.xbar - conf.int.hw,
	                conf.int.xbar + conf.int.hw))
    }
    ,
    greater = {
      list(p.value = 1 - pt(ret.val$statistic, ret.val$parameters),
	   conf.int = c(conf.int.xbar -
	     qt((1 - alpha), ret.val$parameters) * conf.int.s,
	     NA)
        )
    }
    ,
    less = {
      list(p.value = pt(ret.val$statistic, ret.val$parameters),
	   conf.int = c(NA,
	     conf.int.xbar +
	     qt((1 - alpha), ret.val$parameters) * conf.int.s))
    }
    ))
  names(ret.val$statistic) <- "t"
  names(ret.val$parameters) <- "df"
  attr(ret.val$conf.int, "conf.level") <- conf.level
  ret.val <- ret.val[c("statistic", "parameters", "p.value", "conf.int",
    "estimate", "null.value", "alternative", "method", "data.name")]
  oldClass(ret.val) <- "htest"
  # Saddlepoint inferences are not currently supported for two-sample problems
  if(!is.null(y) && !paired)
    return(ret.val)

  #### Saddlepoint Inferences
  if(!is.null(y))
    x <- d	# work with differences as the data
  if(mu <= min(x) || mu >= max(x)){
    # warning("mu is outside the range of the data")
    tau0 <- NA
    p <- ifelse(mu >= x[1],1,0)
  }
  else {
    tau0 <- tiltMeanSolve(mu, x)$tau.exp
    p <- revSaddlepointP(tau0, x)[,3]
  }
  n <- length(x)
  # Make one-sided p-value the larger of the calculated p-value
  # and (1/2)^n.  This lower bound is correct for a symmtric distribution,
  # otherwise is an approximation.  And calculate for three cases.
  switch(alt.expanded,
	 two.sided = {
	   tau.ci <- revSaddlepointPSolve(c(alpha/2, 1-alpha/2), x)[,2]
	   conf.int <- tiltMean(tau.ci, x)$q
	   p.value <- max(2 * min(p, 1-p), 2*.5^n)
	 },
	 greater = {
	   tau.ci <- revSaddlepointPSolve(alpha, x)[,2]
	   conf.int <- c(tiltMean(tau.ci, x)$q, NA)
	   p.value <- max(p, .5^n)
	 },
	 less = {
	   tau.ci <- revSaddlepointPSolve(1-alpha, x)[,2]
	   conf.int <- c(NA, tiltMean(tau.ci, x)$q)
	   p.value <- max(1-p, .5^n)
	 }
	 )
  attr(conf.int, "conf.level") <- conf.level
  ret.val2 <- list(tau0 = tau0, p.value = p.value,
		  tau.ci = tau.ci, conf.int = conf.int,
		  method = "saddlepoint")
  attr(ret.val, "saddlepoint") <- ret.val2
  oldClass(ret.val) <- c("htest.saddlepoint", "htest")
  ret.val
}
# Change from list to usual htest with attribute and class inheritance
# Change name from tilting to saddlepoint
# Change to default var.equal=F.  Change "ignored" warning for var.equal.
# No warn if mu outside range; minimum p-value is 1/n (or 2/n).
# Allow nx=1 or ny=1 if var.equal=T
# Fix warning about ignoring var.equal in paired case.
# Add a treatment argument
# Change F to FALSE, better line breaks
# Better data.name when treatment is supplied
# rename to saddlepoint.test, was t.test


print.htest.saddlepoint <-
function(x, ...)
{
  cat("Usual Students-t inferences (assume zero skewness):\n")
  # NextMethod("print") # This causes problems in some circumstances
  y <- x
  oldClass(y) <- oldClass(y)[-1]
  temp <- print(y, ...)
  
  cat("\nSaddlepoint inferences (allow non-zero skewness):\n")
  y <- attr(x, "saddlepoint")
  cat(paste("\t method = ", y$method, "\n", sep = ""))
  cat("p-value =", format(round(y$p.value, 4)), "\n")
  if(!is.null(y$conf.int)) {
    cat(format(100 * attr(y$conf.int, "conf.level")),
      "percent confidence interval:\n", format(c(y$conf.int[1], y$conf.int[2])),
      "\n")
  }
  invisible(x)
}
# work around bug 27786
# Different handling of arguments; pass ... to print.  This does mean that ...
# are not passed to format for the additional printout.
