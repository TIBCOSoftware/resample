# Copyright 2021. TIBCO Software Inc.
# This file is subject to the license terms contained
# in the license file that is distributed with this file.

# This file contains functions for cross-validation.

# Functions:
# crossValidation
# crossValidation.formula
# crossValidation.default
# balancedSample  (replaced by version in sampmc.q, RT, 5/16/01)
# resampMakeFuncErrF 
# resampMakeFuncErrXY
# resampMakeFuncCrossVal 
# print.crossValidation

crossValidation <- function(x, ...) UseMethod("crossValidation")

crossValidation.formula <- function(x, data, modelFit, K = n, 
			     args.modelFit = NULL, 
			     predFun = function(object, newdata, ...) 
			       predict(object, newdata = newdata, ...),
			     args.predFun = NULL, passOldData.predFun = F, 
			     errFun = function(y, fitted) mean((y - fitted)^2),
			     args.errFun = NULL, seed = .Random.seed,  
			     label = NULL,
			     trace = resampleOptions()$trace,
                             assign.frame1 = F, save.indices = F){
	
  func.call <- match.call()		 # Capture call.
  parent.frame <- sys.parent(1)
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
  fit.func <- 
    resampMakeFuncErrF(modelFit = modelFit, data.name = data.name,
		       is.df.data = is.df.data,
		       is.null.args.modelFit = is.null(args.modelFit), 
		       assign.frame1 = assign.frame1, 
		       dimLength = data.dimLength)
  # is: function(inds, formula, data, modelFit, args.modelFit)
  pred.func <- 
    resampMakeFuncCrossVal(data.name = data.name, 
			   is.null.args.predFun = is.null(args.predFun), 
			   assign.frame1 = assign.frame1, 
			   dimLength = data.dimLength, 
			   passOldData.predFun = passOldData.predFun)
  # is: function(inds, obj, data, predFun, args.predFun)

  # Faster subscripting for a data frame, allow duplicate row names.
  if(is.df.data && is.null(attr(data, "dup.row.names")))
    attr(data, "dup.row.names") <- T
  
  if(assign.frame1) on.exit(if(exists(data.name, frame = 1)) 
			    remove(data.name, frame = 1)) 
  if(passOldData.predFun) on.exit(if(exists("oldData", frame = 1))
				  remove("oldData", frame = 1))
  set.seed(seed)
  seed.start <- .Random.seed

  n <- numRows(data)
  K <- trunc(K)
  if(K < 2) stop("K must be at least 2")
  if (K > n) stop("K must be no more than the number of observations")

  # No on.exit() code for early termination because the algorithm
  # is a balanced one--we need predictions for all the observations
  # so we are not going to try to make sense out of partial results.

  # Set up the groups:
  if(K == n) indices <- 1:n 
  else indices <- balancedSample(n = K, size = n)
  fit.cv <- rep(NA, n)
  # To avoid scoping problems in certain prediction functions,
  # like predict.loess, assign modelFit to frame 1.
  assign("modelFit", modelFit, frame = 1)

  on.exit(if(exists("modelFit", frame = 1)) remove("modelFit", frame=1))

  for(i in 1:K){
    if(trace) cat("Fitting group ",i,"\n")
    currinds <- (1:n)[indices != i]
    obj <- fit.func(inds = currinds, formula = formula, data = data, 
		    modelFit = modelFit, args.modelFit = args.modelFit)
    fit.cv[-currinds] <- pred.func(inds = currinds, obj = obj, data = data, 
				   predFun = predFun, 
				   args.predFun = args.predFun)
  }
  # Now extract the entire response (first removing truncated data set):

  if(exists(data.name, frame = 1)) remove(as.character(data.name), 
			 frame = 1)
  y <- resampMakeArgument(formula[[2]], data, parent.frame)

  if(is.null(args.errFun)) error <- errFun(y, fit.cv)
  else	error <- do.call("errFun", c(list(y, fit.cv), args.errFun))

  if(trace) cat("\n")

  result <- list(call = func.call, fitted = fit.cv, K = K, error = error,
		 seed.start = seed.start, seed.end = .Random.seed, 
		 indices = if(save.indices) indices,
		 label = label,
		 defaultLabel = resampMakeLabel("bootstrapValidation",
		   func.call$x, func.call$data))
  oldClass(result) <- "crossValidation"
  result
}
# moved seed argument towards end of calling sequence

crossValidation.default <- function(x, y, modelFit, K = n,
			     args.modelFit = NULL, 
			     predFun = function(object, newdata, ...) 
			       predict(object, newdata = newdata, ...),
			     args.predFun = NULL, passOldData.predFun = F,
			     errFun = function(y, fitted) mean((y - fitted)^2),
			     args.errFun = NULL, seed = .Random.seed,
			     label = NULL,
			     trace = resampleOptions()$trace,
                             assign.frame1 = F, save.indices = F){

  func.call <- match.call()		 # Capture call.
 
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
  fit.func <- 
    resampMakeFuncErrXY(modelFit = modelFit, x.name = x.name, y.name = y.name,
			is.df.x = is.df.x, is.df.y = is.df.y, 
			is.null.args.modelFit = is.null(args.modelFit), 
			assign.frame1 = assign.frame1, 
			x.dimLength = x.dimLength, y.dimLength = y.dimLength)
  # is: function(inds, x, y, modelFit, args.modelFit)
  pred.func <- 
    resampMakeFuncCrossVal(data.name = x.name, 
			   is.null.args.predFun = is.null(args.predFun), 
			   assign.frame1 = assign.frame1, 
			   dimLength = x.dimLength, 
			   passOldData.predFun = passOldData.predFun)
  # is: function(inds, obj, data, predFun, args.predFun)

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

  n <- numRows(x)

  K <- trunc(K)
  if(K < 2) stop("K must be at least 2")
  if (K > n) stop("K must be no more than the number of observations")

  # No on.exit() code for early termination because the algorithm
  # is a balanced one--we need predictions for all the observations
  # so we are not going to try to make sense out of partial results.
 
  # Set up the groups:
  if(K == n) indices <- 1:n
  else indices <- balancedSample(n = K, size = n)
  fit.cv <- rep(NA, n)
  # To avoid scoping problems in certain prediction functions,
  # like predict.loess, assign modelFit to frame 1.
  assign("modelFit", modelFit, frame = 1)
 
  on.exit(if(exists("modelFit", frame = 1)) remove("modelFit", frame=1))


  for(i in 1:K){
    if(trace) cat("Fitting group ",i,"\n")
    currinds <- (1:n)[indices != i]
    obj <- fit.func(inds = currinds, x = x, y = y, modelFit = modelFit, 
		    args.modelFit = args.modelFit)
    fit.cv[-currinds] <- pred.func(inds = currinds, obj = obj, data = x, 
				   predFun = predFun, 
				   args.predFun = args.predFun)
  }
  if(is.null(args.errFun)) error <- errFun(y, fit.cv)
  else error <- do.call("errFun", c(list(y, fit.cv), args.errFun))
  if(trace) cat("\n")
 
  result <- list(call = func.call, fitted = fit.cv, K = K, error = error,
		 seed.start = seed.start, seed.end = .Random.seed,
		 indices = if(save.indices) indices,
		 label = label,
		 defaultLabel = resampMakeLabel("bootstrapValidation",
		   func.call$x, func.call$y))
  oldClass(result) <- "crossValidation"
  result
}
# moved seed argument towards end of calling sequence


# balancedSample superceded by version in sampmc.q (RT, 5/16/01)
# 
# balancedSample <- function(n, size){
	# Sample of size `size' with treplacement from integers 1:n,
	# with no observation occuring more than once more than any other.

#	if(size <= n) sample(n, size)
#	else sample(c(rep(1:n, size %/% n), sample(n, size %% n)))
#	}
	
	
resampMakeFuncErrF <- function(modelFit, data.name, is.df.data, 
			       is.null.args.modelFit, assign.frame1 = F, 
			       dimLength){

  # Construct a function which takes arguments
  #  inds:          indices for subscripting
  #  formula:	    model formula
  #  data:          data frame containing data set to be modelled
  #  modelFit:      function or expression
  #  args.modelFit: optional arguments or objects (may be NULL)

  ffs <- "function(inds, formula, data, modelFit, args.modelFit){"

  # create expression to subscript data
  ffs2 <- {
    if(is.df.data) "resampSubDF(data, inds)\n"
    else paste("data[inds",
	       if(dimLength > 1)
	         paste(c(rep(",",dimLength),"drop=F"),collapse=""),
	       "]\n")
  }
  ffs <- paste(ffs, "data <-", ffs2)

  if(assign.frame1) ffs <- paste(ffs, "assign('", data.name, 
				 "', data, frame=1)\n", sep = "")
  # Handle function, or expression
  if(is.function(modelFit))
    ffs <- paste(ffs, 
		 if(is.null.args.modelFit) "modelFit(formula, data = data)}"
		 else 
		   "do.call('modelFit', c(list(formula,data = data), args.modelFit))}")

  else
    ffs <- paste(ffs, "eval(modelFit, as(c(list(formula,", 
		 data.name, "=data)",
		 if(is.df.data) ", data",
		 if(!is.null.args.modelFit) ", args.modelFit",
		 "),'list'))}")

  eval(parse (text = ffs))
}


resampMakeFuncErrXY <- function(modelFit, x.name, y.name, is.df.x, 
				is.df.y, is.null.args.modelFit, 
				assign.frame1 = F, x.dimLength, y.dimLength){
	
  # Construct a function which takes arguments	
  #  inds:          indices for subscripting
  #  x:		  independent variables of data set to be modelled
  #  y: 		  dependent variable of data set to be modelled
  #  modelFit:	  function or expression that does the modelling
  #  args.modelFit: optional arguments or objects (may be NULL)

  ffs <- "function(inds, x, y, modelFit, args.modelFit){"

  # create expressions to subscript x and y
  ffs2 <- {
    if(is.df.x) "resampSubDF(x, inds)\n"
    else paste("x[inds",
	       if(x.dimLength > 1)
	         paste(c(rep(",",x.dimLength), "drop=F"), collapse =""),
	       "]\n")
  }
  ffs3 <- {
    if(is.df.y) "resampSubDF(y, inds)\n"
    else paste("y[inds",
	       if(y.dimLength > 1)
	         paste(c(rep(",",y.dimLength), "drop=F"), collapse =""),
	       "]\n")
  }
  ffs <- paste(ffs, "x <-", ffs2, "y <-", ffs3)

  if(assign.frame1) 
    ffs <- paste(ffs, "assign('", x.name, "', x, frame=1)\n",
		 "assign('", y.name, "', y, frame=1)\n",
		 sep = "")

  # Handle function, or expression
  if(is.function(modelFit))
    ffs <- paste(ffs,
		 if(is.null.args.modelFit) "modelFit(x = x, y = y)}"
		 else
		 "do.call('modelFit', c(list(x=x,y=y),args.modelFit))}")
  else	
    ffs <- paste(ffs, "eval(modelFit, as(c(list(", 
		 x.name, "=x,", 	
		 y.name, "=y)",
		 if(is.df.x) ", x",
		 if(is.df.y) ", y",
		 if(!is.null.args.modelFit) ", args.modelFit", 
		 "),'list'))}")
  eval(parse (text = ffs))
}
			
resampMakeFuncCrossVal <- function(data.name, is.null.args.predFun, 
				   assign.frame1=F, dimLength, 
				   passOldData.predFun){

  # Construct a function which takes arguments
  # inds:         indices for subscripting--leave these cases out
  # obj:	        the fitted model object
  # data:         data frame containing data
  # predFun:      prediction function
  # args.predFun: optional arguments

  pfs <- "function(inds, obj, data, predFun, args.predFun){"

  if (passOldData.predFun)
    pfs <- paste(pfs,
		 "oldData <- data[inds",
		 if(dimLength > 1)
		   paste(c(rep(",",dimLength), "drop=F"), collapse=""),
		 "]\n",
		 "assign('oldData', oldData, frame = 1)\n", 
		 "obj$call$data <- as.name('oldData')\n", sep = "")

  pfs <- paste(pfs, "data <- data[-inds", 
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

print.crossValidation <- function(x, digits = max(options()$digits - 3, 4), ...,
			   printCall = .Options.resample$printCall){
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

  cat("Number of Groups:", x$K,  "\n")
  cat("\nPrediction Error:", format(x$error, digits=digits), "\n")
  invisible(x)
}
# Add printCall (like print.resamp).  Print label, call or defaultLabel.
