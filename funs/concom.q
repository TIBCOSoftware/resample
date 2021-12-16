# concomitants function.

### In this file
# concomitants
# concomitants.default (returns vector)
# concomitants.bootstrap (returns resamp object)
# print.concomitants.bootstrap
# plot.concomitants.bootstrap


# concomitants.bootstrap could be built into bootstrap,
# built into aftermarket functions like print and summary,
# or called from those functions.


##########################################################
# concomitants (S3 generic)
##########################################################
concomitants <- function(x, ...)
  UseMethod("concomitants")

##########################################################
# concomitants.default
##########################################################
concomitants.default <-
function(x, y, qfun, args.qfun = NULL, qx = NULL,
         df = 3, weights = NULL){
  # Concomitants adjustment.  Observed values are x and y
  # from a joint distribution, for which
  # * the distribution of x is known (at least to better accuracy than
  #   given by the empirical distribution of x),
  # * y ~ f(x) + epsilon, where f is a smooth function 
  #   (approximated here by smooth.spline) and epsilon is small noise
  # Return adjusted values of y, based on descrepancy between observed
  # values of x ad its known quantiles (determined by qfun or qx).
  #
  # x, y = vectors of the same length "n"
  # qfun = function, qfun(p, ...) gives quantiles for the distribution of x
  # qfun.args = list, additional arguments to pass to qfun
  #  -or-  (pass either qfun or qx)
  # qx = vector of length n, giving quantiles for the distribution of x.
  #      e.g. qx = qfun(ppoints(n))
  #
  # If weights are present, presumably x was drawn by importance sampling.
  # Let F be the target distribution and G the design distribution for x;
  # i.e. mean(x <= a) ~= G(a) and mean(x <= a, weights=weights) ~= F(a)
  # In this case, qfun should be the inverse of F, and the weighted
  # distribution of qx should be approximately F
  # (the unweighted qx values correspond to G).
  # The output y are values from the weighted distribution for y.
  #
  # df = integer, degrees of freedom used for smoothing spline of y against x
  n <- length(x)
  if(length(y) != n)
    stop("length of x and y must be the same")
  nWeights <- length(weights)
  if(nWeights && nWeights != n)
    stop("length of x and weights must be the same")
  notSortedX <- notSorted(x)
  if(notSortedX)
    o <- order(x)
  if(is.null(qx)) {
    if(nWeights){               # need points corresponding to sort(x)
      qx <- ifelse1(is.null(args.qfun),
                    qfun(ppoints(n, weights = weights[o])),
                    do.call("qfun", c(list(ppoints(n, weights[o])),args.qfun)))
    } else {
      qx <- ifelse1(is.null(args.qfun),
                    qfun(ppoints(n)),
                    do.call("qfun", c(list(ppoints(n)), args.qfun)))
    }
    assign("internalQx", qx, where=1)
  } else {
    if(length(qx) != n)
      stop("length of x and qx must be the same")
    if(notSorted(qx))
      stop("qx should be sorted in increasing order")
  }
  if(notSorted(x))
    qx[o] <- qx
  if(df == 2){ # linear fit
    fit <- lsfit(x=x, y=y, weights=weights)
    return(y + fit$coef[2] * (qx - x))
  }
  # df > 2, use smoothing spline
  fit <- smooth.spline(x=x, y=y, df=df,
                       w = ifelse1(nWeights, weights, rep(1.0,n)))
  y + predict(fit, qx)$y - predict(fit, x)$y
}
# Add weights




##########################################################
# concomitants.bootstrap
##########################################################
concomitants.bootstrap <-
function(x,
	 subset.statistic = 1:p,
	 L = resampGetL(x, frame.eval = frame.eval),
         conv.factor = 0.1,
	 group = resampGetArgument(x, "group",
				    frame.eval = frame.eval),
	 treatment = resampGetArgument(x, "treatment",
	   frame.eval = frame.eval),
	 frame.eval = x$parent.frame)
{
  # x is a bootstrap object
  # subset.statistic can be used to specify that only some columns of
  #	x$replicates should be handled
  # conv.factor is passed to qDiscreteMean.  0.1 is usually enoug to
  #   improve the reliability of the estimates.  A somewhat larger value
  #   n/(n-1) corrects for downward bias of sample variance.
  if(is(x, "concomitants"))
    stop("Cannot do concomitants adjustments on a bootstrap object that has already been adjusted")
  if(is.null(frame.eval))
    frame.eval <- sys.parent(1)
  n <- x$n
  replicates <- x$replicates
  if(anyMissing(replicates))
    stop("missing values found in replicates, cannot do concomitants adjustment")
  p <- ncol(replicates)
  nWeights <- length(x$weights)

  # If both treatment and group, stratify on both;
  # if only treatment, then call it "group".
  if(length(treatment)){
    if(length(group))
      group <- paste(match(group, unique(group)),
		    match(treatment, unique(treatment)))
    else
      group <- treatment
  }
  # Note - for bootstrap2 objects, concomitants only adjusts the main
  # replicates, not the replicates for each sample.

  adjust <- function(y, L, Lstar, group, conv.factor, ...){
    fit <- smooth.spline(x = Lstar, y = y, df = 3)
    B <- length(y)
    n <- length(L)
    Ldagger <- qDiscreteMean(p=ppoints(B), values=L, group=group,
                             conv.factor = conv.factor)
    Ldagger[order(Lstar)] <- Ldagger
    y + predict(fit, Ldagger)$y - predict(fit, Lstar)$y
  }
  if(nWeights)
    adjust <- function(y, L, Lstar, group, conv.factor, weights){
      fit <- smooth.spline(x = Lstar, y = y, df = 3, w=weights)
      B <- length(y)
      n <- length(L)
      o <- order(Lstar)
      Ldagger <- qDiscreteMean(p=ppoints(B, weights = weights[o]),
                               values=L, group=group,
                               conv.factor = conv.factor)
      Ldagger[o] <- Ldagger
      y + predict(fit, Ldagger)$y - predict(fit, Lstar)$y
    }
  Lstar <- 
    ifelse1(missing(L) && !is.null(x$Lstar),
            x$Lstar,
            indexMeans(L, resampGetIndices(x, frame.eval = frame.eval),
                       group = group))

  for(j in subset.statistic)
    replicates[,j] <- adjust(replicates[,j], L[,j], Lstar[,j], group,
                             conv.factor = conv.factor,
                             weights = x$weights)

  means <- colMeans(replicates, na.rm = T, weights = x$weights)
  x$original <- x[c("call", "estimate", "replicates")]
  x$call <- match.call()
  x$replicates <- replicates
  x$L <- L
  x$Lstar <- Lstar
  x$estimate <- 
    data.frame(Mean = means,
	       Bias = means - x$observed[subset.statistic],
	       SE = colStdevs(replicates, na.rm=T, weights = x$weights))
  oldClass(x) <- c("concomitants.bootstrap", "bootstrap", "resamp")
  x
}
# Use group when calling indexMeans
# Add arguments group and treatment
# Add conv.factor argument
# allow weights

# Should allow subset.statistic to be logical or character.


##################################################
# methods
##################################################

print.concomitants.bootstrap <-
function(x, ...)
{
  cat("Original bootstrap call:\n")
  print(x$original$call)
  cat("\n")
  NextMethod("print")
  invisible(x)
}

update.concomitants.bootstrap <- function(object, ...){
  stop("A concomitants.bootstrap object cannot be updated; update the original bootstrap object")
}

# Deprected
print.concomitants <-
function(x, ...)
{
  warning('The concomitants class is deprecated, replaced by the concomitants.bootstrap class.  Please do: oldClass(x) <- c("concomitants.bootstrap", "bootstrap", "resamp") to convert this object.')
  cat("Original bootstrap call:\n")
  print(x$original$call)
  cat("\n")
  NextMethod("print")
  invisible(x)
}

# Deprected
update.concomitants <- function(object, ...){
  warning('The concomitants class is deprecated, replaced by the concomitants.bootstrap class.  Please do: oldClass(x) <- c("concomitants.bootstrap", "bootstrap", "resamp") to convert this object.')
  stop("A concomitants object cannot be updated; update the original bootstrap object")
}
