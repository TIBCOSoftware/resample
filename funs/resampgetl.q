# In this file:
#  resampGetL (S3 generic)
#  resampGetL.default
#  resampGetL.bootstrap
#  resampGetL.jackknife
#  resampGetL.influence
#  resampGetL.bootstrap2
#  linearApproxReg

##########################################################
# resampGetL
##########################################################

resampGetL <- function(x, ...) UseMethod("resampGetL")


##########################################################
# resampGetL.default
##########################################################

resampGetL.default <- function(x, ..., frame.eval = NULL){
  # x should be a resampling object with a call, or a call
  # The call should have argument data, statistic, ...
  # Use resampGetL.jackknife (which might call influence)
  Call <- ifelse1(is.call(x), x, x$call)
  if(is.null(Call))
    stop("x must be a resampling object with a call")
  if(is.null(frame.eval))
    frame.eval <- x$parent.frame
  if(is.null(frame.eval))
    frame.eval <- Call$frame.eval
  if(is.null(frame.eval))
    frame.eval <- sys.parent()
  resampGetL.jackknife(Call, frame.eval = frame.eval, ...)
}
# Change function completely -- now takes resamp object with a call.

##########################################################
# resampGetL.bootstrap
##########################################################

resampGetL.bootstrap <-
function(x, method, ...,
	 model.mat = NULL, formula = NULL, data = NULL,
	 frame.eval = x$parent.frame)
{
  # x is a bootstrap object
  # Calculate influence function approximation
  # For actual calculations call one of:
  #   resampGetL.jackknife   "jackknife" method
  #   influence(returnL=T)   "influence" method
  #   linearApproxReg       other methods (regression & model matrix)

  if(is.null(frame.eval))
    frame.eval <- sys.parent(1)
  n <- x$n
  modifiedStatistic <- NULL # potentially used for "influence" method

  message <- resampCheckIfOrderMatters(x$order.matters)
  if(!is.null(message))
    stop("Cannot calculate L--", message)

  if(missing(method)){
    if(!is.null(x$L))
      return(x$L)
    if(!missing(model.mat))
      method <- "ace"
    else if(x$B > 2*n+100)
      method <- "ace"
    else if(is.null(x$call$group))
      method <- "jackknife"  # no group
    else{# group
      statistic <- x$call$statistic
      if(is.name(statistic))
        statistic <- eval(statistic, x$parent.frame)
      if(is.function(statistic)){
        if(any(is.element(c("weights","..."), names(statistic))))
          method <- "influence"
        else
          method <- "jackknife"
      }
      else {
        # create a modifiedStatistic, add a weights argument
        modifiedStatistic <- addWeightsInCall(statistic)
        if(identical(statistic, modifiedStatistic))
          method <- "jackknife"
        else
          method <- "influence"
      }
    }
  }
  else
    method <- match.arg(method, c("ace", "regression",
				 "jackknife", "influence"))
  # Extract some arguments from x, to avoid re-evaluating them
  subject <- resampGetArgument(x, "subject", frame.eval=frame.eval)
  group   <- resampGetArgument(x, "group", frame.eval=frame.eval)

  if(method == "jackknife"){
    Call <- modifyCall(x$call,
		      "data", "statistic", "args.stat", "assign.frame1",
		      group = group,
		      subject = subject)
    L <- resampGetL.jackknife(Call, method = "jackknife",
                              frame.eval = frame.eval, ...)
  } else if(method == "influence") {
    Call <- modifyCall(x$call, NEWFUN = "influence",
		      "data", "statistic", "args.stat", "assign.frame1",
		      group = group,
		      subject = subject,
                      modifiedStatistic = modifiedStatistic,
		      returnL = T, ...)
    # ... are arguments to influence, like epsilon

    # If statistic is random, and epsilon not specified, make it larger
    if(any(x$statistic.is.random) # treat NULL as FALSE
       && is.null(Call$epsilon)){
      Call$epsilon <- max( 1/sqrt(n), .001)
      warning("Statistic is random; using larger value of epsilon = ",
	      Call$epsilon, " for computing influence function approximation")
    }
    L <- eval(Call, frame.eval)
  } else { # method not jackknife or influence
    # Calculate a regression or model matrix approximation
    if(!is.null(formula)){
      # Create the model matrix; pass it rather than formula and data
      Call <- modifyCall(match.call(), "formula", "data",
                        NEWFUN = "model.matrix")
      names(Call)[names(Call)=="formula"] <- ""
      # If original data was a data frame, use that as a default data
      if(is.null(data) &&
         is.data.frame(eval(x$call$data, frame.eval)))
        Call$data <- x$call$data
      model.mat <- eval(Call, frame.eval)
    }
    L <- linearApproxReg(x$replicates,
                         indices = resampGetIndices(object = x,
                           frame.eval = frame.eval),
                         model.mat = model.mat, n = n,
                         group = group, subject = subject,
                         weights = x$weights,
                         transform = (method == "ace"),
                         ...)
  }
  attr(L, "method") <- method
  L
}
# I renamed user.model.reg to model.matrix; want to avoid word "user"
# Use matrix capability of indexMeans
# Changed name from resamp.get.L to resampGetL
# Stop if weights detected, unless jackknife.
# change arg names: boot.obj to x, for consistency with other methods
# add sampling by subject capability for methods ace, regression, jackknife
# add stratified capability for methods model.matrix(.ace), jackknife
# NOTE: centered results in jackknife method, even in non-stratified case
# Call new function for regression and model.matrix calculations
# Call resampGetL.jackknife for jackknife
# Add method "influence" (call to influence()); is default in some cases
# Add ...; passed to other functions (e.g. df, epsilon)
# Removed methods "model.matrix*" -- just do them if model.mat is supplied.
# Move frame.eval to after ... in argument list, for consistency
# Add arguments formula and data; use to create model.mat
# If original data was a data frame, use that by default.
# If statistic is random, use larger epilon when calling influence()
# Stop if order matters (e.g. resampling residuals).

# In all regression and ace versions, could include the original sample
# as one observation.

# Should support weights.



##########################################################
# resampGetL.jackknife
##########################################################

resampGetL.jackknife <-
function(x, method = "jackknife", ..., frame.eval = x$parent.frame)
{
  # x is a jackknife object, or 
  #      a call with arguments (data, statistic, ...)
  # method "jackknife" performs the ordinary jackknife estimate
  #        "influence" calls influence with the same data. 

  if(!missing(method))
    method <- match.arg(method, c("jackknife", "influence"))
  if(is.null(frame.eval))
    frame.eval <- sys.parent()

  if(method == "influence"){
    Call <- (if(is.call(x)) x else x$call)
    Call$seed <- Call$group.size <- NULL
    Call[[1]] <- Quote(influence)
    if(length(l <- list(...)))
      for(j in names(l))
	Call[[j]] <- l[[j]]
    Call$returnL <- T
    L <- eval(Call, frame.eval)
    attr(L, "method") <- method
    return(L)
  }

  # Jackknife method
  if(is.call(x)){
    Call <- modifyCall(x, NEWFUN = "jackknife",
                       "data", "statistic", "args.stat",
                       "group", "subject", "assign.frame1",
                       save.group = 1, save.subject = 1,
                       group.size = 1)
    x <- eval(Call, frame.eval)
  }
  else if(!is(x, "jackknife"))
    stop("x is not a jackknife object or a call")

  # Have jackknife object, use jackknife method.
  n <- x$n
  if(n != x$B){
    # jackknife object originally created with group.size > 1
    Call <- x$call
    Call$group.size <- 1
    x <- eval(Call, frame.eval)
  }

  # Get group and subject; if both then get one group for each subject.
  group <- resampGetArgument(x, "group")
  subject <- resampGetArgument(x, "subject")
  if(length(subject)){
    # L is ordered by the sorted subject names (how split orders its results)
    # n = number of subjects; group and subject are longer
    if(length(group)){
      # need one value of group for each subject, sorted by subject
      subjectIndices <- split(1:length(subject), subject) 
      group <- sapply(subjectIndices, function(i, group) group[i[[1]]], 
                      group = group)
      subjectNames <- names(subjectIndices)
    }
    else
      subjectNames <- as.character(sort(unique(subject)))
  }

  # Core calculations
  if(length(group)){
    # multiply L for each group by (group size - 1)
    L <- -(x$replicates - rep(x$observed, each = n))
    groupIndices <- split(1:n, group)
    for(inds in groupIndices)
      L[inds,] <- L[inds,] * (length(inds) - 1)
  }
  else
    L <- (-(n - 1) * (x$replicates - rep(x$observed, each = n)))
  if(length(subject))
    dimnames(L)[[1]] <- subjectNames
  L <- subtractMeans(L, group = group)
  attr(L, "method") <- method
  L
}
# In influence case, returnL=T
# Allow x to be a call rather than a jackknife object
# Support group argument.
# Add ...
# Support subject argument.
# Move frame.eval to after ... in argument list, for consistency
# Trivial change, to save a level of indentation


##########################################################
# resampGetL.influence
##########################################################

resampGetL.influence <-
function(x, ...)
{
  # x is an influence object
  if(is.null(x$L))
    stop("x does not have a component L; should be an influence object")
  x$L
}


##########################################################
# resampGetL.bootstrap2
##########################################################
resampGetL.bootstrap2 <- function(x, ..., frame.eval = x$parent.frame){
  # x should be a bootstrap2 object
  # call resampGetL for each of the two bootstrap objects, take the difference
  if(is.null(frame.eval))
    frame.eval <- sys.parent(1)
  L1 <- resampGetL(x$bootstrap.objects[[1]], ..., frame.eval=frame.eval)
  L2 <- resampGetL(x$bootstrap.objects[[2]], ..., frame.eval=frame.eval)
  treatment <- resampGetArgument(x, "treatment", frame.eval=frame.eval)
  L <- matrix(nrow = x$n, ncol = length(x$observed))
  treat1 <- treatment == treatment[1]
  L[treat1,]  <-  L1
  L[!treat1,] <- -L2
  if(!is.null(x$ratio) && x$ratio){ # bootstrap ratio rather than difference
    L <- L / rep(x$bootstrap.objects[[2]]$observed, each=nrow(L))
    L[!treat1,] <- L[!treat1,] * rep(x$observed, each = sum(!treat1))
  }
  L
}
# Fix calculation of L in ratio case


##########################################################
# influence function approximation (regression or model matrix)
##########################################################
linearApproxReg <-
function(replicates, indices, n = max(indices),
	 model.mat = NULL, formula = NULL, data = NULL,
	 group = NULL, subject = NULL, weights = NULL,
	 transform = T, df = 3, details = T, ...)
{
  # Calculate regression approximation to influence function values,
  # using bootstrap (or other resampling) replicates and indices.
  # If model.mat is supplied then use model matrix variation.
  # In univariate case (replicates is a vector of univariate statistics),
  # return vector L such that
  #    replicates[i] ~= c + mean(L[indices[,i]])
  sqrtW <- sqrt(weights)			# may be NULL

  regress <- function(x, y, sqrtW){
    if(length(sqrtW))
      fit <- lm.fit.qr(x*sqrtW, y*sqrtW, singular.ok=T)
    else
      fit <- lm.fit.qr(x, y, singular.ok=T)
    bad <- which.na(fit$coefficients)
    fit$coefficients[bad] <- 0
    fit
  }

  if(length(group) && length(subject)){
    temp1 <- which(!duplicated(subject))
    group <- group[temp1][order(subject[temp1])]
  }

  if(!is.null(formula)){
    # Create the model matrix
    Call <- modifyCall(match.call(), "formula", "data", NEWFUN = "model.matrix")
    names(Call)[names(Call)=="formula"] <- ""
    model.mat <- eval(Call, sys.parent())
  }

  Mat <- function(x, colnames){
    # Convert a vector with names into a matrix with those row names
    # This is needed because starting with S+7.0, lm.fit.qr drops
    # a single-column matrix to a vector.
    if(is.matrix(x))
      return(x)
    matrix(x, dimnames = list(names(x), colnames))
  }
  Mat$colnames = colIds(replicates)

  if(is.null(model.mat)){
    #
    # regression or ace; regress replicates against frequencies
    #
    counts <- t(colTabulate(indices, nbins=n))

    # Want to regress against proportions rather than counts; but
    # that would require dividing this big counts matrix by n
    # (or by numbers in each group).  Instead multiply coefficients later.
    lmfit <- regress(counts, replicates, sqrtW)
    # regression with no intercept

    if(!transform)
      L <- Mat(lmfit$coefficients)
    else {
      # transform the replicates for linearity, and regress again
      # Make the derivative of the transformation equal 1 at the mean
      preds <- Mat(lmfit$fitted.values)
      replicates2 <- replicates * NA
      for(j in 1:ncol(replicates)){
	rj <- replicates[,j]
	if(anyMissing(rj))
	  next
	splineFit <- smooth.spline(x=rj, y=preds[,j], df=df)
	replicates2[,j] <-
	  predict(splineFit, rj)$y / predict(splineFit, mean(rj), deriv=1)$y
      }
      L <- Mat(regress(counts, replicates2, weights)$coefficients)
    }

    if(length(group)){
      # Fit is singular; coefficient for last observation in each group
      # should be NA (but regress() converts them to 0's).
      #
      # For each group, subtract group mean and multiply by length of group.
      L <- subtractMeans(L, group) * c(table(group)[as.character(group)])
    }
    else				#simple rescale
      L <- n * subtractMeans(L)
  } else {
    #
    # model matrix approximation
    #
    if(numRows(model.mat) != n)
      stop("Wrong number of rows in the model.mat, must match number of observations (or subjects)")

    # Add an intercept if there is no constant term
    if(any(model.mat[,1] != 1))
      model.mat <- cbind(1, model.mat)

    mm <- indexMeans(x = model.mat, indices = indices, group=group)
    lmfit <- regress(mm, replicates, sqrtW)

    if(!transform)
      L <- model.mat %*% Mat(lmfit$coefficients)
    else {
      # transform the replicates for linearity, and regress again
      # Make the derivative of the transformation equal 1 at the mean
      preds <- Mat(lmfit$fitted.values)
      replicates2 <- replicates * NA
      for(j in 1:ncol(replicates)){
	rj <- replicates[,j]
	if(anyMissing(rj))
	  next
	splineFit <- smooth.spline(x=rj, y=preds[,j], df=df)
	replicates2[,j] <-
	  predict(splineFit, rj)$y / predict(splineFit, mean(rj), deriv=1)$y
      }
      coef <- Mat(regress(mm, replicates2, sqrtW)$coefficients)
      L <- model.mat %*% coef
    }
    L <- subtractMeans(L, group = group)
  }
  if(details)
    attr(L, "correlation") <- cor(replicates, na.method = "available",
				 indexMeans(L, indices, group=group))
  if(length(subject))
    dimnames(L)[[1]] <- as.character(sort(unique(subject)))
  L
}
# Add support for multivariate statistics
# Add details argument; if T then add correlation as an attribute
# Fix handling of intercept in model.mat case (failed when supplied).
# Add arguments formula and data, use to create model.mat
# call indexMeans with the group argument
# Convert outputs from lm.fit.qr to matrices (needed for change in S+7.0)
# Missing values - give NA instead of stop when using spline.
#   (Could go further, calculate L from nonmissing observations.)
