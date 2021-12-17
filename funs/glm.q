# Copyright 2021. TIBCO Software Inc.
# This file is subject to the license terms contained
# in the license file that is distributed with this file.

#
# New version of glm: (RT, 4/12/01)
# 
# glm
# glm.default
# glm.model.list
#
# Define as old-style generic.  Default method will be the current function in 
# Splus6, modified to handle `method=model.list'. New method glm.model.list.
# See help for glm.model.list and model.list.object. 
#
# Dispatch on the first argument, `formula'.  
# 
# Same as glm.IVc.q in ~rthurman/glm.
#
glm <- 
function(formula, family = gaussian, data, 
	 weights, subset, na.action, start, control = glm.control(...),
	 method = "glm.fit", model = F, x = F, y = T, contrasts = NULL, ...)
  UseMethod("glm")
# make argument list same as glm.default (so gui code works)

glm.default <-
function(formula, family = gaussian, data, 
	 weights, subset, na.action, start, control = glm.control(...),
	 method = "glm.fit", model = F, x = F, y = T, contrasts = NULL, ...)
{
  fit.call <- match.call()
  parent.frame <- sys.parent()
  m <- match.call(expand = F)
  m$family <- m$method <- m$model <- m$x <- m$y <- m$control <- 
    m$contrasts <- m$... <- NULL
  m$drop.unused.levels <- T
  m[[1]] <- as.name("model.frame")
  m <- eval(m, parent.frame)
  Terms <- attr(m, "terms")
  if(any(method == "model.frame"))
    return(m)
  xvars <- as.character(attr(Terms, "variables"))
  if(length(xvars) > 0) {
    xlevels <- lapply(m[xvars], levels)
    xlevels <- xlevels[!sapply(xlevels, is.null)]
    if(length(xlevels) == 0)
      xlevels <- NULL
  }
  else xlevels <- NULL

  Y <- model.extract(m, response)
  X <- model.matrix(Terms, m, contrasts)
  w <- model.extract(m, weights)
  if(!length(w))
    w <- rep(1, nrow(m))
  else if(any(w < 0))
    stop("negative weights not allowed")
  start       <- model.extract(m, start)
  offset      <- model.extract(m, offset)
  na.action   <- attr(m, "na.action")
  family      <- as.family(family)
  fit.class   <- c("glm","lm")
  fit.formula <- as.vector(attr(Terms, "formula"))

  modelList <- any(method == "model.list")

  # if method is a vector, then elements are "model.list" and fitting method
  if(length(method) > 1)
    method <- method[method != "model.list"]
  else
    if(missing(method) || modelList){
      method <- attr(family, "method")      
      if(is.null(method))
	method <- "glm.fit"
    }
  if(!existsFunction(method))
    stop(paste("unimplemented method:", method))

  #
  # glmMakeFitExp constructs an expression to perform the fit.
  # The expression will either be evaluated later in this routine, or if
  # method="model.list", returned as part of the model list for later
  # evaluation. 
  # 
  if(modelList){# return the model list

    fitexpr <- glmMakeFitExpr(fitOnly = F, method = method, offset = offset,
                             Terms = Terms, xlevels = xlevels, ...)
    # two calls will go in the model.list: 
    #   thiscall is the  call creating the model list, and 
    #   fit.call is the call that generates the fit.
    # they differ only in the method argument. 
    thiscall <- fit.call
    fit.call$method <- method  

    return(structure(list(data=list(X = X, Y = Y, weights = w), 
			  fit=list(fitexpr = fitexpr,
			    control = control, 
			    start = start, 
			    offset = offset, 
			    family = family,
			    method = method, 
			    class = fit.class,
			    call = fit.call,
			    Terms = Terms,
			    na.action = na.action, 
			    xlevels = xlevels,
                            dotargs = list(...)),
			  formula=fit.formula,
			  call=thiscall,
                          parent.frame = parent.frame),
		     class="model.list"))
  }
  fitexpr <- glmMakeFitExpr(fitOnly = T, method = method, offset = offset,
                            Terms = Terms, xlevels = xlevels, ...)
  #
  # Evaluate the fit expression
  #
  fit <- eval(fitexpr)

  #
  # Some final detailing not available to glm.model.list, so handled 
  # here rather than in fitexpr.
  #
  if(model)
    fit$model <- m
  if(x)
    fit$x <- X
  if(!y)
    fit$y <- NULL
  fit
}

#
# glm.model.list
#
# `formula' is an object of class `model.list'.
# If no other arguments are specified, evaluate the fit expression (fast). 
# Otherwise, evaluate the original call to lm with the new arguments
# appended (essentially starting from scratch).  A special case if only
# weights are to be added, for then we can simply adjust the model list and
# evaluate the fit expression.
#
glm.model.list <- 
function(formula = formula(data), family = gaussian, data = sys.parent(), 
	weights, subset, na.action, start = eta, control = glm.control(...),
	method = "glm.fit", model = F, x = F, y = T, contrasts = NULL, ...){
  thisCall <- match.call()
  if((callLen <- length(thisCall)) > 2){
    # extra arguments to add
    if(is.null(frame.eval <- formula$parent.frame))
      frame.eval <- sys.parent(1)
    if(callLen == 3 && names(thisCall)[[3]] == "weights"){
      # Special case when only weights are to be added.
      formula$fit$call$weights <- thisCall$weights
      formula$data$weights <- eval(weights, frame.eval) 
      names(formula$data$weights) <- dimnames(formula$data$X)[[1]]
      eval(formula$fit$fitexpr)
    }
    else{
      # evaluate the original call to lm, appending new arguments
      formula$call$"method" <- NULL
      eval(c(formula$call, thisCall[-(1:2)]), local = frame.eval)
    }
  }
  else
    eval(formula$fit$fitexpr)
}


glmMakeFitExpr <- function(fitOnly, method, offset, Terms, xlevels, ...){
  #
  # Make fit expression for glm
  #
  # Arguments:
  # fitOnly  logical; if T the returned expression is ready for evaluation
  #          by glm.default.  If F a data extraction expression is prepended
  #          to the fit expression; this version of the fit expression 
  #          typically becomes part of a model list to be evaluated from within
  #          glm.model.list.
  # method   
  # offset
  # Terms
  # xlevels
  # ...      same meaning as inside glm.default.

  fcall <- list(fcall=call(method, x = Quote(X), y = Quote(Y), w = Quote(w), 
		     start = Quote(start), offset = Quote(offset), 
		     family = Quote(family), maxit = Quote(control$maxit),
		     epsilon = Quote(control$epsilon),  
		     trace = Quote(control$trace), 
		     null.dev = T,...))

  fitexpr <- Quote({fit <- fcall})
  fitexpr <- substitute(fitexpr,fcall,evaluate=T)

  #
  # Second piece of the fit expression: null deviance. 
  #
  # If an offset and intercept are present, iterations are needed to
  # compute the null deviance; these are done by another call to the fitter
  # unless the model is NULL, in which case the computations have been done
  # by the original fit. 
  #
  if(any(offset) && attr(Terms, "intercept"))
    fitexpr <- 
      if(length(Terms)){
	null.dev.call  <- list(null.dev.call=call(method, 
				 x = Quote(X[, "(Intercept)",drop = F]), 
				 y = Quote(Y), w = Quote(w), 
				 offset = Quote(offset), 
				 family = Quote(family), 
				 maxit = Quote(control$maxit),
				 epsilon = Quote(control$epsilon),  
				 null.dev = NULL))
	fitexpr <- c(fitexpr,
		     Quote({fit$null.deviance <- null.dev.call$deviance}))
	substitute(fitexpr,null.dev.call,eval=T)
      }
      else
	c(fitexpr, Quote({fit$null.deviance <- fit$deviance}))
  #
  # Final piece: add descriptive information to the fitted object.
  #
  fitexpr <- c(fitexpr, Quote({oldClass(fit)        <- fit.class;
			       fit$terms            <- Terms;
			       fit$formula          <- fit.formula;
			       fit$call             <- fit.call;
			       fit$control          <- control;
			       fit$na.action        <- na.action}),
	       if(!is.null(xlevels)) Quote({attr(fit, "xlevels") <- xlevels}), 
	       Quote({fit}))

  if(fitOnly) return(fitexpr)

  # Prepend a data extraction expression for evaluation within
  # glm.model.list.  All of the objects needed to perform the fit are now
  # part of the model.list, and this part of the expression extracts those
  # objects in preparation for doing the fit. Argument `formula' is the first
  # argument to glm.model.list, and refers to the model.list. 
  c(Quote({X           <- formula$data$X;
           Y           <- formula$data$Y;
           w           <- formula$data$weights;
           control     <- formula$fit$control;
           start       <- formula$fit$start;
           offset      <- formula$fit$offset;
           family      <- formula$fit$family;
           method      <- formula$fit$method;
           fit.class   <- formula$fit$class;
           fit.call    <- formula$fit$call;
           Terms       <- formula$fit$Terms;
           na.action   <- formula$fit$na.action;
           xlevels     <- formula$fit$xlevels;
           fit.formula <- formula$formula}), 
    fitexpr)
}
