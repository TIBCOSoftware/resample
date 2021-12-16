# file contains functions used for bootstrap prediction.

# Function names:
# bootstrapValidation
# bootstrapValidation.formula
# bootstrapValidation.default
# resampMakeFuncBootPred
# print.bootstrapValidation

bootstrapValidation <- function(x, ...) UseMethod("bootstrapValidation")

bootstrapValidation.formula <- function(x, data, modelFit, B,
			     group = NULL, subject = NULL,
			     args.modelFit = NULL,
			     predFun = function(object, newdata, ...)
			       predict(object, newdata = newdata, ...),
			     args.predFun = NULL, passOldData.predFun = F,
			     errFun = function(y, fitted){(y - fitted)^2},
			     args.errFun = NULL, seed = .Random.seed,
			     label = NULL,
			     trace = resampleOptions()$trace,
                             assign.frame1 = F, save.indices = F,
			     save.group = NULL, save.subject = NULL,
			     save.errors = FALSE){

  func.call <- match.call()		 # Capture call.
  parent.frame <- sys.parent()
  formula <- x
  # If modelFit isn't a function or name of function or expression,
  # store it as an unevaluated expression to pass to fit.func.
  substitute.modelFit <- substitute(modelFit)
  if(!is.element(mode(substitute.modelFit), c("name", "function")))
    modelFit <- substitute.modelFit

  # Get name of data.
  data.name <- substitute(data)
  if(!is.name(data.name)) data.name <- "data"

  # Set passOldData.predFun to TRUE for predict.gam
  # (which is the only known prediction function that refits the model).
  if(is.element("gam", all.names(func.call))) passOldData.predFun <- T

  # Get function to evaluate modelFit given data and indices.
  is.df.data <- is.data.frame(data)
  data.dimLength <- length(dim(data))
  fit.func <- resampMakeFuncErrF(modelFit = modelFit,
				 data.name = data.name,
				 is.df.data = is.df.data,
				 is.null.args.modelFit =
				 is.null(args.modelFit),
				 assign.frame1 = assign.frame1,
				 dimLength = data.dimLength)
  # is: function(inds, formula, data, modelFit, args.modelFit)
  pred.func <- resampMakeFuncBootPred(data.name = data.name,
				      is.null.args.predFun =
				      is.null(args.predFun),
				      assign.frame1 = assign.frame1,
				      dimLength = data.dimLength,
				      passOldData.predFun = passOldData.predFun)
  # is: function(inds, oldInds, obj, data, predFun, args.predFun)

  # Faster subscripting for a data frame, allow duplicate row names.
  if(is.df.data && is.null(attr(data, "dup.row.names")))
    attr(data, "dup.row.names") <- T

  if(assign.frame1)
    on.exit(if(exists(data.name, frame = 1))
	    remove(as.character(data.name), frame = 1))
  if(passOldData.predFun) on.exit(if(exists("oldData", frame = 1))
				  remove("oldData", frame = 1))
  set.seed(seed)
  seed.start <- .Random.seed

  n <- nObs <- numRows(data)  # n changed later to number of subjects
  B <- trunc(B)

  # To avoid scoping problems in certain prediction functions,
  # like predict.loess, assign modelFit to frame 1.
  assign("modelFit", modelFit, frame = 1)
  on.exit(if(exists("modelFit", frame = 1)) remove("modelFit", frame=1))

  fit.nullmodel <- fit.func(inds = 1:n, formula = formula,
			    data = data, modelFit = modelFit,
			    args.modelFit = args.modelFit)
  yhat.nullmodel <- pred.func(inds = 1:n, oldInds = 1:n,
			      obj = fit.nullmodel, data = data,
			      predFun = predFun,
			      args.predFun = args.predFun)

  # extract response
  y <- resampMakeArgument(formula[[2]], data, parent.frame)

  err.full <- ifelse1(is.null(args.errFun),
                      mean(errFun(y, yhat.nullmodel)),
                      mean(do.call("errFun", c(list(y, yhat.nullmodel),
                                               args.errFun))))
  if(is.null(args.errFun))
    gammahat <- sum(outer(y, yhat.nullmodel, errFun)) / n^2
  else {
    # check for clashes among argument lists
    # if so, use a wrapper function to avoid clash
    if(!any(is.element(names(outer), names(args.errFun))))
      gammahat <- sum(do.call("outer", c(list(X = y,
					      Y = yhat.nullmodel, FUN = errFun),
					 args.errFun))) / n^2
    else {
      errWrap <- function(x, y, errFun, args.errFun)
	do.call("errFun", c(list(x,y),args.errFun))
      gammahat <- sum(outer(y, yhat.nullmodel, errWrap,
			    errFun = errFun, args.errFun =args.errFun))/n^2
    }
  }

  err.optimistic <- matrix(0, nrow = n, ncol = B)
  err.bootsamp <- rep(0, B)

  # Save group or subject arguments?
  if(is.null(save.subject)) save.subject <- n <= 10000
  if(is.null(save.group))   save.group   <- n <= 10000

  #
  # Sampling by subject
  #
  # If data is a data frame, look for subject there first.
  if(!missing(subject) && is.df.data)
    subject <- resampMakeArgument(func.call$subject, data, parent.frame)
  bySubject <- length(subject)
  if(bySubject){
    if(bySubject != n)
      stop("length of subject must match number of observations in data")
    subjectIndices <- split(1:n, subject)
    isubject <- ifelse1(is.character(subject), match(subject, unique(subject)),
                        subject)
    n <- length(subjectIndices)
  }

  #
  # Stratified sampling ("group" argument)
  # The sampler will be called in one of two ways, depending
  # on whether there are groups:
  #    sampler(n, B)
  #    sampler(groupSizes[iGroup], B)
  #
  if(!missing(group)) {
    # If data is a data frame, look for group there first.
    if(is.df.data)
      group <- resampMakeArgument(func.call$group, data, parent.frame)
    stratified <- length(group)
    if(stratified){
      if(stratified != nObs)
	stop("length of group != number of observations in the data")
      if(bySubject){
	# get group for each subject (subject must be nested within group)
	subjectByGroup <- lapply(split(subject, group), unique)
	if(notNested(subjectByGroup))
	  stop("some values of `subject' occur in multiple groups; subject must be nested within the `group' argument")
        group.inds <- lapply(subjectByGroup, match, table=names(subjectIndices))
      }
      else
        group.inds <- split(1:n, group)
      nGroups <- length(group.inds)
      groupSizes <- sapply(group.inds, length)
      all.indices <- matrix(as.integer(0), n, B)
      for(iGroup in 1:nGroups){
	thisGroup <- group.inds[[iGroup]]
	all.indices[thisGroup, ] <-
	  thisGroup[samp.bootstrap(n=groupSizes[iGroup], B = B)]
      }
    }
  }
  else
    all.indices <- samp.bootstrap(n=n, B=B)
  
  for(j in 1:B){
    if(trace) cat("Fitting model ",j,"\n")
    currinds <- ifelse1(bySubject,
                        unlist(subjectIndices[all.indices[,j]]),
                        all.indices[,j])
    currfit <- fit.func(inds = currinds, formula = formula,
			data = data, modelFit = modelFit,
			args.modelFit = args.modelFit)
    yhat.fulldata <- pred.func(inds = 1:nObs, oldInds = currinds,
			       obj = currfit, data = data,
			       predFun = predFun,
			       args.predFun = args.predFun)
    yhat.bootsamp <- yhat.fulldata[currinds]
    if(is.null(args.errFun))
      err.optimistic[,j] <- errFun(y, yhat.fulldata)
    else
      err.optimistic[,j] <-
	do.call("errFun", c(list(y, yhat.fulldata), args.errFun))
    err.bootsamp[j] <- mean(err.optimistic[currinds, j], na.rm=T)
  }
  if(trace) cat("\n")
  replicatesMissing <- NULL
  observationsMissing <- NULL
  if(anyMissing(err.optimistic)){
    replicatesMissing <- which(colSums(is.na(err.optimistic)) > 0)
    observationsMissing <- which(rowSums(is.na(err.optimistic)) > 0)
    warning("Missing values in error estimates for ",
	    length(replicatesMissing), " bootstrap replications, and ",
            ifelse1(bySubject, " subjects", " observations."))
  }

  optimism <- mean(colMeans(err.optimistic, na.rm=T) - err.bootsamp)

  # Compute .632 estimate; use bootstrap samples without each observation
  has.match <- function(samp, inds, w)
    duplicated(c(samp, inds))[w]
  matches.mat <- apply(all.indices, 2, has.match, 1:n, w=(n+1):(n+n))
  # matches.mat[i,j] = 1 if observation/subject i is in bootstrap sample j
  # jabB[i] = number of bootstrap samples that do contain observation i
  jabB <- ncol(matches.mat) - rowSums(matches.mat)
  if(any(jabB == 0))
    warning("increase B for computation of the .632 and .632+ estimators")
  else if(min(jabB) < 5)
    warning("an observation is omitted from only ", min(jabB),
	    " bootstrap samples; .632 and .632+ estimates are inaccurate")
  # Compute means of rows, using columns without matches
  means.err.optimistic <-
    ifelse1(bySubject,
            rowSums(err.optimistic * !(matches.mat[isubject,,drop=F]),
                    na.rm=T) / jabB[isubject],
            rowSums(err.optimistic * !matches.mat, na.rm=T) / jabB)
  eps0 <- mean(means.err.optimistic)
  err632 <- .368 * err.full + .632 * eps0

  # compute .632+ estimator
  Rhat <- (eps0 - err.full) / (gammahat - err.full)
  omegahat <- .632 / (1 - .368 * Rhat)
  err632plus <- (1 - omegahat) * err.full + omegahat * eps0

  result <- list(call = func.call, B = B, apparent.error = err.full,
		 optimism = optimism, err632 = err632,
		 err632plus = err632plus,
		 seed.start = seed.start, seed.end = .Random.seed,
		 parent.frame = parent.frame,
		 label = label,
		 defaultLabel = resampMakeLabel("bootstrapValidation",
		   func.call$x, func.call$data),
		 group = if(save.group) group,
		 subject = if(save.subject) subject,
		 indices = if(save.indices) all.indices,
		 errors = if(save.errors) err.optimistic,
		 replicatesMissing = replicatesMissing,
		 observationsMissing = observationsMissing
                 )
  oldClass(result) <- "bootstrapValidation"
  result
}
# moved seed argument towards end of calling sequence
# added group/subject arguments
# added parent.frame component
# repair group+subject handling
# replace currinds eval with logical test within loop (less memory used)
# Add argument label, return label and defaultLabel
# If assign.frame1, remove data on.exit
# get yhat.bootsamp and err.bootsamp by subscripting rather than from scratch
# Add save.errors argument
# Allow missing values in errors
# Change order of outputs - optional components last
# Clean up code - ifelse1, subjectIndices, gammahat, ...
# 632 estimates with subjects: use all observations, not just 1:nSubjects

bootstrapValidation.default <- function(x, y, modelFit, B, group = NULL, subject = NULL,
			     args.modelFit = NULL,
			     predFun = function(object, newdata, ...)
			       predict(object, newdata = newdata, ...),
			     args.predFun = NULL, passOldData.predFun = F,
			     errFun = function(y, fitted){(y - fitted)^2},
			     args.errFun = NULL, seed = .Random.seed,
			     label = NULL,
			     trace = resampleOptions()$trace,
                             assign.frame1 = F, save.indices = F,
			     save.group = NULL, save.subject = NULL,
			     save.errors = FALSE){

  func.call <- match.call()		 # Capture call.
  parent.frame <- sys.parent()

  # If modelFit isn't a function or name of function or expression,
  # store it as an unevaluated expression to pass to fit.func.
  substitute.modelFit <- substitute(modelFit)
  if(!is.element(mode(substitute.modelFit), c("name", "function")))
    modelFit <- substitute.modelFit

  # Get name of x, y.
  x.name <- substitute(x)
  y.name <- substitute(y)
  if(!is.name(x.name)) x.name <- "x"
  if(!is.name(y.name)) y.name <- "y"

  # Set passOldData.predFun to TRUE for predict.gam
  # (which is the only known prediction function that refits the model).
  if(is.element("gam", all.names(func.call))) passOldData.predFun <- T

  # Get function to evaluate modelFit given x, y, and indices.
  is.df.x <- is.data.frame(x)
  is.df.y <- is.data.frame(y)
  x.dimLength <- length(dim(x))
  y.dimLength <- length(dim(y))
  fit.func <- resampMakeFuncErrXY(modelFit = modelFit,
				  x.name = x.name, y.name = y.name,
				  is.df.x = is.df.x, is.df.y = is.df.y,
				  is.null.args.modelFit =
				  is.null(args.modelFit),
				  assign.frame1 = assign.frame1,
				  x.dimLength = x.dimLength,
				  y.dimLength = y.dimLength)
  # is: function(inds, x, y, modelFit, args.modelFit)
  pred.func <- resampMakeFuncBootPred(data.name = x.name,
				      is.null.args.predFun =
				      is.null(args.predFun),
				      assign.frame1 = assign.frame1,
				      dimLength = x.dimLength,
				      passOldData.predFun =
				      passOldData.predFun)
  # is: function(inds, oldInds, obj, data, predFun, args.predFun)

  # Faster subscripting for a data frame, allow duplicate row names.
  if(is.df.x && is.null(attr(x, "dup.row.names")))
    attr(x, "dup.row.names") <- T
  if(is.df.y && is.null(attr(y, "dup.row.names")))
    attr(y, "dup.row.names") <- T
  if(assign.frame1)
    on.exit({
      if(exists(x.name, frame = 1)) remove(x.name, frame = 1)
      if(exists(y.name, frame = 1)) remove(y.name, frame = 1)
    })
  if(passOldData.predFun) on.exit(if(exists("oldData", frame = 1))
				  remove("oldData", frame = 1))
  set.seed(seed)
  seed.start <- .Random.seed

  n <- nObs <- numRows(x)  # n changed later to number of subjects
  B <- trunc(B)

  # To avoid scoping problems in certain prediction functions,
  # like predict.loess, assign modelFit to frame 1.
  assign("modelFit", modelFit, frame = 1)
  on.exit(if(exists("modelFit", frame = 1)) remove("modelFit", frame=1))

  fit.nullmodel <- fit.func(inds = 1:n, x = x, y = y,
			    modelFit = modelFit,
			    args.modelFit = args.modelFit)
  yhat.nullmodel <- pred.func(inds = 1:n, oldInds = 1:n,
			      obj = fit.nullmodel, data = x,
			      predFun = predFun,
			      args.predFun = args.predFun)

  err.full <- ifelse1(is.null(args.errFun),
                      mean(errFun(y, yhat.nullmodel)),
                      mean(do.call("errFun", c(list(y, yhat.nullmodel),
                                               args.errFun))))
  if(is.null(args.errFun))
    gammahat <- sum(outer(y, yhat.nullmodel, errFun)) / n^2
  else {
    # check for clashes among argument lists
    # if so, use a wrapper function to avoid clash
    if(!any(is.element(names(outer), names(args.errFun))))
      gammahat <- sum(do.call("outer", c(list(X = y,
					      Y = yhat.nullmodel, FUN = errFun),
					 args.errFun))) / n^2
    else {
      errWrap <- function(x, y, errFun, args.errFun)
	do.call("errFun", c(list(x,y),args.errFun))
      gammahat <- sum(outer(y, yhat.nullmodel, errWrap,
			    errFun = errFun, args.errFun =args.errFun))/n^2
    }
  }

  err.optimistic <- matrix(0, nrow = n, ncol = B)
  err.bootsamp <- rep(0, B)

  # Save group or subject arguments?
  if(is.null(save.subject)) save.subject <- n <= 10000
  if(is.null(save.group))   save.group   <- n <= 10000

  #
  # Sampling by subject
  #
  bySubject <- length(subject)
  if(bySubject){
    if(bySubject != n)
      stop("length of subject must match number of observations in data")
    subjectIndices <- split(1:n, subject)
    isubject <- ifelse1(is.character(subject), match(subject, unique(subject)),
                        subject)
    n <- length(subjectIndices)
  }
  #
  # Stratified sampling ("group" argument)
  # The sampler will be called in one of two ways, depending
  # on whether there are groups:
  #    sampler(n, B)
  #    sampler(groupSizes[iGroup], B)
  #
  if(stratified <- length(group)){
    if(stratified != nObs)
      stop("length of group != number of observations in the data")
    if(bySubject){
      # get group for each subject (subject must be nested within group)
      subjectByGroup <- lapply(split(subject, group), unique)
      if(notNested(subjectByGroup))
	stop("some values of `subject' occur in multiple groups; subject must be nested within the `group' argument")
      group.inds <- lapply(subjectByGroup, match, table=names(subjectIndices))
    }
    else
      group.inds <- split(1:n, group)
    nGroups <- length(group.inds)
    groupSizes <- sapply(group.inds, length)
    all.indices <- matrix(as.integer(0), n, B)
    for(iGroup in 1:nGroups){
      thisGroup <- group.inds[[iGroup]]
      all.indices[thisGroup, ] <-
	thisGroup[samp.bootstrap(n=groupSizes[iGroup], B = B)]
    }
  }
  else {
    stratified <- F
    all.indices <- samp.bootstrap(n=n, B=B)
  }

  for(j in 1:B){
    if(trace) cat("Fitting model ",j,"\n")
    currinds <-
      if(bySubject) unlist(subjectIndices[all.indices[,j]])
      else all.indices[,j]
    currfit <- fit.func(inds = currinds, x = x, y = y,
			modelFit = modelFit,
			args.modelFit = args.modelFit)
    yhat.fulldata <- pred.func(inds = 1:nObs, oldInds = currinds,
			       obj = currfit, data = x,
			       predFun = predFun,
			       args.predFun = args.predFun)
    yhat.bootsamp <- yhat.fulldata[currinds]
    if(is.null(args.errFun))
      err.optimistic[,j] <- errFun(y, yhat.fulldata)
    else
      err.optimistic[,j] <-
	do.call("errFun", c(list(y, yhat.fulldata), args.errFun))
    err.bootsamp[j] <- mean(err.optimistic[currinds, j], na.rm=T)
  }
  if(trace) cat("\n")
  replicatesMissing <- NULL
  observationsMissing <- NULL
  if(anyMissing(err.optimistic)){
    replicatesMissing <- which(colSums(is.na(err.optimistic)) > 0)
    observationsMissing <- which(rowSums(is.na(err.optimistic)) > 0)
    warning("Missing values in error estimates for ",
	    length(replicatesMissing), " bootstrap replications, and ",
	    length(observationsMissing),
            ifelse1(bySubject, " subjects", " observations."))
  }

  optimism <- mean(colMeans(err.optimistic, na.rm=T) - err.bootsamp)

  # Compute .632 estimate; use bootstrap samples without each observation
  has.match <- function(samp, inds, w)
    duplicated(c(samp, inds))[w]
  matches.mat <- apply(all.indices, 2, has.match, 1:n, w=(n+1):(n+n))
  # matches.mat[i,j] = 1 if observation/subject i is in bootstrap sample j
  # jabB[i] = number of bootstrap samples that do contain observation i
  jabB <- ncol(matches.mat) - rowSums(matches.mat)
  if(any(jabB == 0))
    warning("increase B for computation of the .632 and .632+ estimators")
  else if(min(jabB) < 5)
    warning("an observation is omitted from only ", min(jabB),
	    " bootstrap samples; .632 and .632+ estimates are inaccurate")
  # Compute means of rows, using columns without matches
  means.err.optimistic <-
    ifelse1(bySubject,
            rowSums(err.optimistic * !(matches.mat[isubject,,drop=F]),
                    na.rm=T) / jabB[isubject],
            rowSums(err.optimistic * !matches.mat, na.rm=T) / jabB)
  eps0 <- mean(means.err.optimistic)
  err632 <- .368 * err.full + .632 * eps0

  # compute .632+ estimator
  Rhat <- (eps0 - err.full) / (gammahat - err.full)
  omegahat <- .632 / (1 - .368 * Rhat)
  err632plus <- (1 - omegahat) * err.full + omegahat * eps0

  result <- list(call = func.call, B = B, apparent.error = err.full,
		 optimism = optimism, err632 = err632,
		 err632plus = err632plus,
		 seed.start = seed.start, seed.end = .Random.seed,
		 parent.frame = parent.frame,
		 label = label,
		 defaultLabel = resampMakeLabel("bootstrapValidation",
		   func.call$x, func.call$y),
		 group = if(save.group) group,
		 subject = if(save.subject) subject,
		 indices = if(save.indices) all.indices,
		 errors = if(save.errors) err.optimistic,
		 replicatesMissing = replicatesMissing,
		 observationsMissing = observationsMissing
                 )
  oldClass(result) <- "bootstrapValidation"
  result
}
# (see also list of changes for bootstrapValidation.formula)
# fixed calculation of err.bootsamp (from "err.bootsamp <- ..." to
#   "err.bootsamp[i] <- ...")



resampMakeFuncBootPred <- function(data.name, is.null.args.predFun,
				   assign.frame1=F, dimLength,
				   passOldData.predFun){

  # Construct a function which takes arguments
  # inds:         indices for subscripting--cases to do predictions on
  # oldInds:	indices of the data used to fit obj
  # obj:          the fitted model object
  # data:         data frame containing data
  # predFun:      prediction function
  # args.predFun: optional arguments

  pfs <- "function(inds, oldInds, obj, data, predFun, args.predFun){"

  if (passOldData.predFun)
    pfs <- paste(pfs, "oldData <- data[oldInds",
		 if(dimLength > 1)
		 paste(c(rep(",",dimLength),"drop=F"),collapse=""),
		 "]\n", "assign('oldData', oldData, frame = 1)\n",
		 "obj$call$data <- as.name('oldData')\n", sep = "")
  pfs <- paste(pfs, "data <- data[inds",
               if(dimLength > 1)
	       paste(c(rep(",",dimLength),"drop=F"),collapse=""),
	       "]\n")

  if (assign.frame1) pfs <- paste(pfs, "assign('", data.name,
				  "', data, frame=1)\n", sep = "")

  pfs <- paste(pfs,
	       if(is.null.args.predFun) "predFun(obj, data)}"
	       else "do.call('predFun',c(list(obj,data), args.predFun))}"
	       )
  eval(parse(text = pfs))
}

print.bootstrapValidation <- function(x, digits = max(options()$digits - 3, 4), ...,
			   printCall = .Options.resample$printCall){
  # Print a bootstrapValidation object.  Print call and/or label or defaultLabel,
  # Then error rates and optimism.
  old.digits <- options(digits = digits)
  on.exit(options(old.digits))
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
  if(!is.null(x$replicatesMissing)){
    if(length(x$replicatesMissing) > 10)
      cat("Number of Replications with Missing values:",
	  length(x$replicatesMissing),"\n")
    else
      cat("Missing values in replications:", x$replicatesMissing,"\n")
    if(length(x$observationsMissing) > 10)
      cat("Number of Observations (or subjects) with Missing values:",
	  length(x$observationsMissing),"\n")
    else
      cat("Missing values for observations:", x$observationsMissing,"\n")
  }
  cat("\nApparent Error Rate:", format(x$apparent.error,
				       digits = digits), "\n")
  cat("\noptimism:",format(x$optimism, digits = digits), "\n")
  cat("\n.632 Prediction Error:", format(x$err632, digits = digits),"\n")
  cat("\n.632+ Prediction Error:", format(x$err632plus, digits=digits),
      "\n")
  invisible(x)
}

# Add printCall (like print.resamp).  Print label, call or defaultLabel.
# Add support for missing values
