# bootstrap and other resampling functions

# Functions in this file:
#-----------------
# jackknifeAfterBootstrap
# plot.jackknifeAfterBootstrap
# print.jackknifeAfterBootstrap
# resampFunctionalList  (a list of functions)


##########################################################
# jackknifeAfterBootstrap
##########################################################

jackknifeAfterBootstrap <-
function(boot.obj,
         functional = NULL,
         graphical = NULL,
         passObserved = NULL,
         jack.obj = NULL,
         threshold = 2,
         subset.statistic = 1:p,
         control = "choose",
         moments = 2,
         crossCorr = FALSE,
         ...,
	 frame.eval = boot.obj$parent.frame)
{
  # This does one of two things:
  #   (1) graphical diagnostics for a functional
  #   (2) analytical estimates of bias and/or SE of a functional
  # By default, graphical = TRUE (unless functional is one of
  #   "Mean" "Bias" "SE" "Bias&SE" "Mean&SE"
  # Conversely, if graphical = TRUE the default functional is
  #    "Quantiles".
  #
  # The functional is a summary of a bootstrap distribution
  # (and may depend on the observed value).  Built-in options include:
  #    mean
  #    bias (mean - observed value)
  #    standard error   (standard deviation of the bootstrap distn)
  #    quantiles
  #    centered quantiles (quantiles - observed value)
  #    standardized quantiles (centered quantiles / sd)
  # It may also be a function with arguments:
  #    functional(x, weights, ...)
  #    functional(x, observed, weights, ...)
  # depending on whether passObserved is FALSE or TRUE
  #
  # Graphical diagnostics: by default, plot centered quantiles
  # of the leave-one-out bootstrap distributions against an
  # estimate of the influence function of observations.
  #
  # Analytical estimates: for example, one may use the bootstrap to
  # estimate the standard error of a statistic.  That standard error
  # is a functional of the bootstrap distribution (the standard deviation).
  # We may use the analytical JAB to estimate the bias or standard
  # error of that functional.
  #
  # Caution - there is a big problem with the analytical estimates:
  # jackknife formulas estimate bias or SE by multiplying small
  # quantities by n, and if those quantities are calculated with
  # random simulation error, the errors are magnified.  In this case,
  # the quantities are calculated using random bootstrap sampling.
  #
  # Terminology used here:
  #   jackknife value = leave out one observation, calculate statistic
  #   leave-one-out functional = leave out one observation, draw
  #      bootstrap distribution, calculate functional such as mean
  #      or standard error
  #   observed value = this has two meanings.  First, it is the
  #      statistic calculated for the original sample.  Second, when
  #      calculating the leave-one-out functionals, it is a jackknife
  #      value, i.e. the statistic calculated from a jackknife sample.
  #
  # Other arguments:
  # passObserved = should be true for functionals such as "bias" that
  #   require both the bootstrap replicates and the observed value.
  # crossCorr = calculate the correlation between leave-one-out
  #   functionals and jackknife values
  # jack.obj = jackknife object (used for jackknife statistics)
  #   set this to F to prevent jackknife from being called
  # threshold = flag observations whose influence on the functional
  #   exceeds this many standard deviations (for "non-graphical" case).
  # subset.statistic = if the observed value has length > 1, which
  #   dimensions to work with.  Only these dimensions are passed to
  #   the functional.
  # control = options are "choose", "concomitants", "controlVariates", "none"
  #   "choose" gives concomitants for quantiles (unless there are weights),
  #      and controlVariates for known moment-based statistics
  # moments = integer, 1 or 2; if controlVariates are used, how many
  #   moments to use.
  # ... = additional arguments to pass to the functional
  if(is.null(frame.eval))
    frame.eval <- sys.parent(1)
  if(!is(boot.obj, "bootstrap"))
    stop("boot.obj must be a 'bootstrap' object.")
  if(!is.null(jack.obj) && !(identical(jack.obj, FALSE) || is(jack.obj, "jackknife")))
    stop("jack.obj must be a 'jackknife' object, or FALSE.")
  if(any(boot.obj$resampleResiduals))
    stop("Cannot do jackknife after bootstrap when resampling residuals")
  if(!is.null(boot.obj$call$group))
    warning("jackknifeAfterBootstrap may give incorrect results when stratified sampling was used")
  func.call <- match.call()
  weights <- boot.obj$weights
  observed <- boot.obj$observed # This may be subsetted later
  obs.names <- names(observed)  # This may be subsetted later
  p <- length(observed)
  n <- boot.obj$n
  inds <- 1:n
  control <- match.arg(control, c("choose", "concomitants",
                                  "controlVariates", "none"))
  if(!is.null(boot.obj$call$group)){
    if(is.element(control, c("concomitants", "controlVariates"))){
      warning("control variates and concomitants are not supported when stratified sampling was used")
    }
    control <- "none"
  }

  # Pick default value for functional
  if(is.null(functional)){
    if(is.null(graphical))
      graphical <- TRUE
    functional <- ifelse1(graphical,
                          "Quantiles",
                          "Bias&SE")
    func.call$functional <- functional
  }

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
    if(is.null(graphical))
      graphical <- quantileBased
    if(control == "choose" && !quantileBased)
      control <- "controlVariates"  # later will set control to concomitants unless there are weights
  }
  else {
    if(!is.function(functional))
      stop("functional must be a function, or a recognized name of a functional")
    if(is.null(passObserved))
      passObserved <- is.element("observed", names(functional))
  }
  if(is.null(graphical))
    graphical <- TRUE
  if(passObserved)
    crossCor <- TRUE

  if(!is.null(weights)){
    if(control == "concomitants")
      warning("concomitants is not supported with weights, switching to controlVariates")
    if(is.element(control, c("choose", "concomitants")))
      control <- "controlVariates"
  }
  if(control == "choose" && quantileBased)
    control <- "concomitants"

  # subset.statistic
  if(is.character(subset.statistic)){
    j <- match(subset.statistic, names(x$observed))
    if(anyMissing(j))
      stop("subset.statistic names not recognized:\n    ",
	   paste(subset.statistic[which.na(j)], collapse=", "))
    subset.statistic <- j
  }
  observed <- observed[subset.statistic]
  obs.names <- obs.names[subset.statistic]
  
  orig.functional <- functional
  if(is.character(functional))
    functional <- resampFunctionalList[[functional]]

  # End of argument handling; now do calculations
  # Calculate functional for full sample
  # (May redo this later, with controlVariates or concomitants)
  func.full <-
    ifelse1(passObserved,
            functional(boot.obj$replicates[,subset.statistic,drop=F],
                       observed = observed[subset.statistic],
                       weights = weights, ...),
            functional(boot.obj$replicates[,subset.statistic,drop=F],
                       weights = weights, ...))
  dim.Func <- dim(func.full)
  dimnames.Func <- dimnames(func.full)

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
                           probs <- dimnames.Func[[1]]
                           ifelse1(length(boot.obj$observed) == 1,
                                   probs,
                                   t(outer(obs.names, probs, paste)))
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

  if(identical(jack.obj, FALSE)){
    if(graphical || passObserved || crossCorr)
      stop("A jackknife object is required for the arguments you supplied")
    if(is.element(control, c("concomitants", "controlVariates"))){
      warning("Cannot use ", control, " without a jackknife object")
      control <- "none"
    }
    L <- resampGetL(boot.obj)
  }
  else {
    if(is(boot.obj, "bootstrap2"))
      stop("cannot call jackknife for two-sample problems")
    # Get jackknife values, if they weren't passed in
    if(is.null(jack.obj)){
      boot.call <- boot.obj$call
      jack.call <- modifyCall(boot.call, NEWFUN="jackknife",
                              "data", "statistic", "args.stat", "group",
                              "subject", "label", "statisticNames",
                              "assign.frame1", "save.group", "save.subject",
                              seed = boot.obj$seed.end)
      jack.obj <- eval(jack.call, frame.eval)
    }
    else { # make sure dimensions of jackknife object match bootstrap results
      if(numRows(jack.obj$replicates) != n)
        stop("bootstrap/jackknife object mismatch; number of rows of jackknife replicates should be same as boot.obj$n")
      if(numCols(jack.obj$replicates) != p)
        stop("bootstrap and jackknife objects do not have the same number of columns")
    }
    L <- resampGetL(jack.obj)
  }

  # Get resampling indices, locate matches.
  inds.mat <- resampGetIndices(object = boot.obj, frame.eval = frame.eval)
  samplen <- numRows(inds.mat)
  has.match <- function(samp, inds, w)
    duplicated(c(samp, inds))[w]
  matches.mat <- apply(inds.mat, 2, has.match, inds, w=(samplen+1):(samplen+n))
  jabB <- ncol(matches.mat) - rowSums(matches.mat)
  if(any(jabB == 0))
    stop("At least one jackknife-after-bootstrap sample has no replications, cannot calculate influence.  Increase B or use a different sampler.")
  if(any(jabB < 5))
    warning("At least one jackknife-after-bootstrap sample is of size < 5.")

  # Calculate functional for jab samples
  # If control = "controlVariates", then modify weights
  # If control = "concomitants", modify values
  # When statistic is multi-dimensional, concomitants adjustments are
  # done independently but controlVariates reweighting is done for all.
  bootReps <- boot.obj$replicates[,subset.statistic,drop = F]
  L.original <- L
  if(!identical(subset.statistic, 1:p))
    L <- L[,subset.statistic,drop=F]
  if(is.element(control, c("concomitants", "controlVariates")))
    Lstar <- indexMeans(L, inds.mat)

  # Main loop - omit each observation in turn, and compute functional
  func.replicates <- matrix(nrow = n, ncol = length(func.full))
  concomitantsFail <- NULL
  for(i in 1:n){
    notInSamp <- !matches.mat[i,]
    w <- weights[notInSamp]
    y <- bootReps[notInSamp,,drop=F]
    if(control == "controlVariates"){
      mu <- colMeans(L[-i,])
      Lstari <- Lstar[notInSamp,,drop=F]
      Bi <- jabB[i]  # nrow(Lstari)
      if(moments == 2)
        w <- controlVariates(cbind(Lstari,
                                   (Lstari - rep(mu, each=Bi))^2),
                             mu = c(mu, colVars(L[-i,],unbiased=F)/(n-1)),
                             weights=weights[notInSamp])
      else
        w <- controlVariates(Lstari, mu, 
                             weights=weights[notInSamp])
    }
    else if(control == "concomitants"){
      Lstari <- Lstar[notInSamp,,drop=F]
      Bi <- jabB[i]  # nrow(Lstari)
      for(j in 1:ncol(y)){
        qx <- try(qDiscreteMean(p=ppoints(Bi), values=L[-i,j],
                                conv.factor = 1))
        # qDiscreteMean sometimes fails
        if(is(qx, "Error"))
          concomitantsFail <- c(concomitantsFail, i)
        else
          y[,j] <- concomitants(Lstari[,j], y[,j], qx = qx)
      }
    }
    func.replicates[i,] <-
      ifelse1(passObserved,
              functional(y,
                         observed = jack.obj$replicates[i,subset.statistic,drop = F],
                         weights = w, ...),
              functional(y, weights = w, ...))
  }

  if(all(is.na(func.replicates)))
    stop("All calculated values are missing; this could occur if there are missing values in the bootstrap replicates and your functional does not omit missing values")
  if(length(concomitantsFail))
    warning("concomitants failed in ",
            ifelse1(length(concomitantsFail) < 6,
                    paste("replications",
                          paste(concomitantsFail, collapse=" ")),
                    paste(length(concomitantsFail), "replications")),
            ", using unadjusted jab distributions")
  dimnames(func.replicates) <- list(inds, func.names)

  # Correct the functional for the full data, if necessary
  if(control == "controlVariates"){
    boot2 <- controlVariates(boot.obj, L = L.original,
                             subset.covariates=subset.statistic,
                             moments=moments)
    func.full[] <-
      ifelse1(passObserved,
              functional(boot.obj$replicates[,subset.statistic,drop=F],
                         observed = observed[subset.statistic],
                         weights = boot2$weights, ...),
              functional(boot.obj$replicates[,subset.statistic,drop=F],
                         weights = boot2$weights, ...))
  }
  else if(control == "concomitants"){
    boot2 <- try(concomitants(boot.obj, L = L.original, conv.factor = 1,
                              subset.statistic=subset.statistic))
    if(is(boot2, "Error"))
      warning("concomitants failed for original bootstrap distribution, using unadjusted bootstrap distribution")
    else
      func.full[] <-
        ifelse1(passObserved,
                functional(boot2$replicates[,subset.statistic,drop=F],
                           observed = observed[subset.statistic],
                           weights = NULL, ...),
                functional(boot2$replicates[,subset.statistic,drop=F],
                           weights = NULL, ...))
  }

  # Analytical calculations (do these even in the graphical case)
  # Calculate mean, bias, SE(s) of the functional
  func.mean <- colMeans(func.replicates, na.rm=T)
  func.bias <- (n - 1) * (func.mean - func.full)
  func.se <- sqrt( (n-1)/n * colVars(func.replicates, SumSquares=T, na.rm=T))

  # Calculate jackknife influence values
  rel.influence <- ( -(n-1)) * scale(func.replicates,
				     center = func.mean,
				     scale = sqrt(n) * func.se)
  # dimnames(rel.influence) <- dimnames(func.replicates)

  # Summary of relative influences
  lri.func <- function(x, rel.inf, thresh)
    rel.inf[abs(rel.inf[, x]) >= thresh, x, drop = F]
  large.rel.influence <- lapply(func.names, lri.func, rel.influence, threshold)
  names(large.rel.influence) <- func.names

  # Return results
  result <- list(call = func.call,
		 Func = func.full, # observed value of the functional
                 estimate = data.frame(
                   Mean.Func = func.mean,
                   Bias.Func = func.bias,
                   SE.Func = func.se),
		 Func.replicates = func.replicates, # replicates of functional
                 observed = boot.obj$observed, # observed statistic
                 jack.obj = jack.obj,          # replicates of the statistic
		 jabB = jabB,
                 graphical = graphical,
                 control = control,
		 rel.influence = rel.influence,
		 large.rel.influence = large.rel.influence,
		 threshold = threshold,
                 n = boot.obj$n,
                 B = boot.obj$B,
                 boot.call = boot.obj$call,
                 parent.frame = frame.eval,
                 L = L,
                 dim.Func = dim.Func,
                 dimnames.Func = dimnames.Func,
                 quantiles = quantileBased)
  if(crossCorr)
    result$cross.corr <- cor(func.replicates,
                             jack.obj$replicates, na.method = "omit")
  oldClass(result) <- "jackknifeAfterBootstrap"
  result
}
# rel.influence[n,p], values.functional[n,p]
#
# Use "switch"
# Use full names for components of bootstrap object
# B was unused
# Additional argument to has.match() for speed.
# Use na.rm for mean & variance
# Change n.param to p
# Change default for frame.eval.boot
# Add subject subject = boot.call$subject to jackknife call
# add frame.eval argument (to eventually supercede frame.eval.boot)
# modify def of matches.mat to handle case of sample-size != n.
# Remove frame.eval.boot

# Should change this to make a distinction between a functional and
# a pivot, and allow the user to supply either.
# Should allow functions to operate on multiple columns (e.g. for t-stat)
# stop if bootstrap2 object.  Use modifyCall.
# Major changes 1/12/07
#   Add graphical approach - quantile-based, like Canty's boot library.
#   New argment "graphical", 
#     whether to do the previous numerical calculations or a
#     quantile-based graphical approach instead. 
#   Add options "Quantiles" and "Centered Quantiles" for functional.   
#   Change default for functional to NULL (equivalent to "Centered Quantiles"
#     unless graphical=F).
#   Add argument subset.statistic
#   Add argument control
# Charlie (?) wrote (and had in the code):
#   The standard error estimates tend to be too large.  I'm interested in
#   finding a well-supported alternative, probably involving weighting.
#   1/12/07 - the latest "major changes" partially address the cause for the
#   standard errors being too large, and discuss in the help file.
# Major changes 4/13/07
#   Default to "Quantiles" in graphical case
#   Add return components
#   Handle built-in statistics differently - stored in resampFunctionalList
#   Add component dim.Func, remove dim.obs
#   If graphical=T, don't set subset.statistic to 1.
#   Allow jack.obj = F; otherwise compute it if not supplied




##########################################################
# plot.jackknifeAfterBootstrap
##########################################################
plot.jackknifeAfterBootstrap <-
function(x, nrow = NULL, grid.layout = T, id.outliers = T,
         threshold = x$threshold, ..., graphical = x$graphical,
         xaxis = ifelse1(graphical, "jackknife", "Observation"),
         type = ifelse1(graphical, "o", "h"),
         superpose = NULL, Func = T,
         subset.x = NULL, subset.plots = NULL, absolute = F)
{
  # Produces either
  # (1) if graphical=TRUE,
  #     plot of jab quantiles against L or jackknife values
  # (2) else, Cook's distance type of plot showing the influence of
  #     each observation upon the functional under consideration.
  # xaxis is one of
  #   "Observation" (default for (2))
  #   "L" - empirical influence function (default for (1))
  #   "jackknife" - jackknife statistics (leave-one-out statistics)
  #   "data" - original data  (must still exist to use this)
  # superpose = if TRUE, then functional should have returned a
  #   matrix; plot all rows of a column together
  # Func = if TRUE, add horizontal dotted lines at the value
  #   of the functional for the full data
  # subset.x = use to subscript a column of X
  #   (where X = original data, L, or jackknife statistics)
  # subset.plots = do only a subset of the possible plots
  # absolute = logical, if TRUE then (2) uses absolute values
  k <- length(x$Func)
  n <- x$n
  xaxis <- match.arg(xaxis, c("L", "jackknife", "data", "Observation"))
  if(xaxis == "jackknife" && is.null(x$jack.obj)){
    warning("jack.obj is not present, will use L instead of jackknife")
    xaxis <- "L"
  }
  if(xaxis == "L" && is.null(x$L)){
    warning("L is not present, will use Observation instead of L")
    xaxis <- "Observation"
  }
  # xlab and data for plotting on x axis
  xlab <- list(...)$xlab
  if(xaxis == "Observation"){
    # ignore subset.x and xlab in this case
    X <- 1:n
    xlab <- "Observation"
  }
  else if(xaxis == "data"){
    frame.eval <- x$parent.frame
    if(is.null(frame.eval))
      frame.eval <- sys.parent(1)
    X <- eval(x$boot.call$data, frame.eval)
    if(is.call(fitCall <- X$call) || is.call(fitCall <- attr(X,"call"))){
      if(is.null(fitCall$data))
        stop("Cannot find data")
      X <- eval(fitCall$data, frame.eval)
    }
    X <- ifelse1(is.data.frame(X), numerical.matrix(X), as.matrix(X))
    if(numCols(X) == 1 && is.null(colIds(X)))
      colIds(X) <- deparse(x$boot.call$data)
  }
  else  # xaxis = "L" or "jackknife"
    X <- ifelse1(xaxis == "L", x$L, x$jack.obj$replicates)
  if(xaxis != "Observation"){
    # It is possible to use L from one one column as the x value for
    # plotting against other columns.
    if(!is.null(subset.x))
      X <- X[ , subset.x, drop=F]
    # Now make xlab the same length as number of X columns.
    # This may differ from the number of plots - if so, will
    # use rep(each)
    if(is.null(xlab))
      xlab <- switch(xaxis,
                     L = paste("Empirical influence for", colIds(X)),
                     jackknife = paste("Jackknife statistics for", colIds(X)),
                     data = colIds(X))
    else
      xlab <- rep(xlab, length = ncol(X))
  }

  ylab <- list(...)$ylab
  if(graphical){  # y axis = functional values (usually quantiles)
    quantileBased <- x$quantiles
    if(is.null(quantileBased))
      quantileBased <- F
    if(is.null(superpose))
      superpose <- quantileBased
    if(superpose && is.null(x$dim.Func))
      warning("functional did not return a matrix, cannot superpose")
    nPlots <- ifelse1(superpose, x$dim.Func[2], k)
    Xcols <- ifelse1(xaxis == "Observation",
                     rep(1, nPlots),
                     rep(1:ncol(X), each = ceiling(nPlots / ncol(X)), length=nPlots))
    if(is.null(ylab)){
      if(superpose){
        if(quantileBased){  # collapse the percentiles, only one "%"
          percs <- x$dimnames.Func[[1]]
          percs <- as.numeric(substring(percs, 1, nchar(percs)-1))
          ylab <- paste(paste(round(percs,2), collapse=", "), "% ", sep="",
                        x$call$functional,
                        " for ", x$dimnames.Func[[2]])
        } else { # less common case - superpose other statistics
          ylab <- paste(as.character(deparse(x$call$functional)),
                        "for", x$dimnames.Func[[2]])
        }
      }
      else
        ylab <- names(x$Func)
    }
    ylab <- rep(ylab, each = nPlots / length(ylab))
    if(is.null(subset.plots))
      subset.plots <- T
    else if(any(subset.plots > nPlots))
      stop("requested nonexistent plot, only ", nPlots, " plots are available")
    for(j in (1:nPlots)[subset.plots]){
      Xj <- ifelse1(is.matrix(X), X[,Xcols[j]], X)
      o <- order(Xj)
      jj <- ifelse1(superpose, matrix(1:k, ncol=nPlots)[,j], j)
      if(superpose) # multiple columns per plot, usually quantiles
        matplot(Xj[o], x$Func.replicates[o, jj],
                type = type,
                xlab = xlab[Xcols[j]], ylab = ylab[j], lty = 1, ...)
      else
        plot(Xj[o], x$Func.replicates[o, j], type = type,
             xlab = xlab[Xcols[j]], ylab = ylab[j], ...)
      if(Func)
        abline(h = x$Func[jj], lty=2, ...)
    }
    return(invisible(NULL))
  }

  # Cook's distance plot
  rel.influence <- x$rel.influence
  if(absolute)
    rel.influence <- abs(rel.influence)
  n <- nrow(rel.influence)
  p <- ncol(rel.influence)
  if(is.null(subset.plots))
    subset.plots <- 1:p
  P <- length(subset.plots)
  if(grid.layout && P != 1) {
    if(is.null(nrow)){
      if(length(x$dim.Func) == 2 && p == P)
	old.par <- par(mfcol = x$dim.Func)
      else {
        nrow <- min(P, 2)
        old.par <- par(mfrow = c(nrow, ceiling(P/nrow)))
      }
    }
    else
      old.par <- par(mfrow = c(nrow, ceiling(P/nrow)))
    on.exit(par(old.par))
  }
  ylim <- c(ifelse1(absolute, 0, min(-3, rel.influence, na.rm=T)),
            max(3, rel.influence, na.rm = T))
  if(is.null(ylab))
    ylab <- ifelse1(absolute, "Absolute Standardized Influence",
                          "Standardized Influence")
  else
    ylab <- rep(ylab, length = k)
  for(j in (1:k)[subset.plots]) {
    # Cook's distance type of plot
    if(all(is.na(rel.influence[, j]))) {
      plot(X, rep(0, n), ..., type = "h", ylim = ylim,
           xlab = xlab, ylab = ylab,
	   main = names(x$Func)[j], sub = "All Values are NA")
      next
    }
    plot(X, rel.influence[, j], ..., type = "h", ylim = ylim,
         xlab = xlab, ylab = ylab,
         main = names(x$Func)[j])
    abline(h = ifelse1(absolute, 1, c(-1,1)) * threshold, lty = 2)
    if(!absolute) abline(h=0)
    if(id.outliers) {
      outliers <- (rel.influence[, j] >= threshold)
      text(X[outliers], rel.influence[outliers, j], which(outliers))
    }
  }
  invisible(NULL)
}
# Avoid using old.par <- par(new values) trick.
# Use {} to clarify if() if() else() else()
# Support dim.obs=scalar
# Change nstats to p
# Define n and p consistently
# Add ... argument
# Oops, restore the old.par <- par(new values) trick (bug 24874)
# Add support for graphical approach, plotting quantiles
# 4/13/07 - major redesign.  Choice of what to plot on x axis
#   More choice of y axis; in influence case, don't need absolute value.

# To consider/do:
# nrow is a bad name for the argument, gives warning
# Good default for nrow?  (What if high-dimensional result?)



##########################################################
# print.jackknifeAfterBootstrap
##########################################################
print.jackknifeAfterBootstrap <-
function(x, digits = max(options()$digits - 3, 4), ...,
         graphical = x$graphical)
{
  old.digits <- options(digits = 4)
  on.exit(options(old.digits))
  cat("Call:\n")
  print(x$call)
  if(is.null(graphical)){
    warning("x was created using an older version of jackknifeAfterBootstrap")
    graphical <- F
  }
  if(graphical){
    # Do not give analytical results
    if(any(x$graphical))
      cat("\njackknife.after.bootstrap was called for graphical diagnostics,\n",
          "so analytical results are not printed by default.\n", sep="")
    return(invisible(x))
  }
  # Remainder of this function is for analytical results
  cat("\nFunctional Under Consideration:\n")
  ifelse1(is.character(x$call$functional),
          cat(x$call$functional, "\n"),
          print(x$call$functional))
  cat("\nFunctional of Bootstrap Distribution of Parameters:\n")
  print(data.frame(Func = x$Func, x$estimate), digits = digits, ...)
  if(!is.null(x$cross.corr)){
    cat("\nCorrelation Between Jackknifed Parameter Values (columns)\n")
    cat("and Jackknifed Functional Values (rows):\n")
    print(x$cross.corr, digits = digits, ...)
  }
  cat("\nObservations with Large Influence on Functional:\n")
  print(x$large.rel.influence)
  if(any(x$jabB < 100)){
    cat("Number of jackknife-after-bootstrap samples for each observation:\n")
    print(matrix(x$jab, nrow=1, dimnames=list("", 1:length(x$jabB))), ...)
  }
  invisible(x)
}
# Add option graphical - if TRUE then print only the call


##################################################
# resampFunctionalList
#
# List of functionals, to be called by jackknifeAfterBootstrap or tiltAfterBootstrap
# These accept as input:
#       x = matrix[B,p] of resampled statistics
#       observed (optional) vector[p] of observed statistics
#       weights = vector[B] of weights, or NULL
#       ... = additional optional arguments
# Output = vector or matrix of summary statistics, for example
# quantiles of each column.  This need not have p columns; for
# example, there could be a summary that depends on two columns.
#
# These built-in versions do not necessarily return an object with names;
# they may be called thousands of times, and it is quicker to assign names
# once, inside the calling function.  However, a user-defined function
# should return a named object, because the names will be used.
resampFunctionalList <-
  list(
       quantiles = function(x, weights = NULL,
         probs = c(0.025, 0.16, 0.5, 0.84, 0.975), ...){
         apply(x, 2, quantile, probs = probs, weights = weights, na.rm = T, ...)
       },
       "centered quantiles" = function(x, observed, weights = NULL,
         probs = c(0.025, 0.16, 0.5, 0.84, 0.975), ...){
         apply(x, 2, quantile, probs = probs, weights = weights, na.rm = T, ...) -
           rep(observed, each = length(probs))
       },
       "standardized quantiles" = function(x, observed, weights = NULL,
         probs = c(0.025, 0.16, 0.5, 0.84, 0.975), ...){
         stdz <- function(x,m,s) (x-rep(m, each=nrow(x)))/rep(s,each=nrow(x))
         stdz(apply(x, 2, quantile, probs = probs, weights = weights,
                    na.rm = T, ...),
              observed,
              colStdevs(x, weights=weights, na.rm=T))
       },
       mean = function(x, weights = NULL, ...){
         mean = colMeans(x, ..., weights = weights, na.rm=T)
       },
       bias = function(x, observed, weights = NULL, ...){
         Bias = colMeans(x, ..., weights = weights, na.rm=T) - observed
       },
       se = function(x, weights = NULL, ...){
         SE = colStdevs(x, ..., weights = weights, na.rm=T)
       },
       "bias&se" = function(x, observed, weights = NULL, ...){
         rbind(Bias = colMeans(x, ..., weights = weights, na.rm=T) - observed,
               SE = colStdevs(x, ..., weights = weights, na.rm=T))
       },
       "mean&se" = function(x, weights = NULL, ...){
         rbind(Mean = colMeans(x, ..., weights = weights, na.rm=T),
               SE = colStdevs(x, ..., weights = weights, na.rm=T))
       }
)
