# Copyright 2021. TIBCO Software Inc.
# This file is subject to the license terms contained
# in the license file that is distributed with this file.

#
# New version of lm: (RT, 4/12/01)
#
# Define as old-style generic.  Default method will be the current function in
# Splus6, modified to handle `method=model.list'. New method lm.model.list.
# See help for lm.model.list and model.list.object.
#
# Originally adapted from ~timh/bootstrap/lm.S
#
# Dispatch on the first argument, `formula'.
#
# Changes, 2/01-3/01 (RT):
#
# * model.list structure changed: no attributes, all information in components.
#   A new data component with subcomponents X,Y, weights. numRows.q and
#   model.list.q changed to accomodate this.
# * new call component to contain call creating model list.  This is referred
#   to by lm.model.list to get arguments for fitting.  lm objects created from
#   a model list can now be reconstructed, using the call component.  Function
#   update now works on this type of lm object (tests in lm.t).
# * lm.model.frame deleted.  As something intermediate between old-style uses
#   of bootstrap and bootstrap.lm, it served no purpose.
#
# 4/01 (RT):
#
# * revamp model.list structure:  include expression to perform the fit and
#   fill out the lm object.
# * lm.model.list becomes the single line:   eval(formula$fit$fitexpr)
#
#
lm <-
function(formula, data, weights, subset, na.action, method = "qr", model = F,
	x = F, y = F, contrasts = NULL, ...)
  UseMethod("lm")
# make argument list same as lm.default (so gui code works)

lm.default <-
function(formula, data, weights, subset, na.action, method = "qr", model = F,
	x = F, y = F, contrasts = NULL, ...)
{
  fit.call <- match.call()
  parent.frame <- sys.parent()
  m <- match.call(expand = F)
  m$method <- m$model <- m$x <- m$y <- m$contrasts <- m$... <- NULL
  m$drop.unused.levels <- T
  m[[1]] <- as.name("model.frame")
  m <- eval(m, parent.frame)

  if(any(method=="model.frame")) # Note: slightly faster than is.element
    return(m)

  Terms <- attr(m, "terms")
  xvars <- as.character(attr(Terms, "variables"))
  if((yvar <- attr(Terms, "response")) > 0)
    xvars <- xvars[ - yvar]
  if(length(xvars) > 0) {
    xlevels <- lapply(m[xvars], levels)
    xlevels <- xlevels[!unlist(lapply(xlevels, is.null))]
    if(length(xlevels) == 0)
      xlevels <- NULL
  }
  else xlevels <- NULL

  weights     <- model.extract(m, weights)
  Y           <- model.extract(m, response)
  X           <- model.matrix(Terms, m, contrasts)
  na.message  <- attr(m, "na.message")
  na.action   <- attr(m, "na.action")

  #
  # lmMakeFitExp constructs an expression to perform the fit.
  # The expression will either be evaluated later in this routine, or if
  # method="model.list", returned as part of the model list for later
  # evaluation.
  #
  if(any(method == "model.list")){# return the model list

    fitexpr <- lmMakeFitExpr(fitOnly = F, x = x, weights = weights,
                           xlevels = xlevels, ...)
    # the fitting method
    method <- if(length(method) > 1) method[method != "model.list"]
             else args(lm.default)$method

    # two calls will go in the model.list:
    #   thiscall is the  call creating the model list, and
    #   fit.call is the call that performs the fit.
    # they differ only in the method argument.
    thiscall <- fit.call
    fit.call$method <- method

    return(structure(list(data = list(X = X, Y = Y, weights = weights),
			  fit = list(fitexpr = fitexpr,
			    method = method,
			    call = fit.call,
			    Terms = Terms,
			    na.message = na.message,
			    na.action = na.action,
			    xlevels = xlevels,
                            dotargs = list(...)),
			  formula = formula,
			  call = thiscall,
                          parent.frame = parent.frame),
		     class="model.list"))
  }

  fitexpr <- lmMakeFitExpr(fitOnly = T, x = x, weights = weights,
                         xlevels = xlevels, ...)
  #
  # Evaluate the fit expression
  #
  fit <- eval(fitexpr)

  #
  # Some final detailing not available to lm.model.list, so handled
  # here rather than in fitexpr.
  #
  if(!model)
    unset(m)
  if(model)
    fit$model <- m
  if(x)
    fit$x <- X
  if(y)
    fit$y <- Y
  fit
}

#
# lm.model.list
#
# `formula' is an object of class `model.list'.
# If no other arguments are specified, evaluate the fit expression (fast).
# Otherwise, evaluate the original call to lm with the new arguments
# appended (essentially starting from scratch).  A special case if only
# weights are to be added, for then we can simply adjust the model list and
# the fit expression, and evaluate the fit expression.
#
lm.model.list <-
function(formula, data, weights, subset, na.action, method = "qr", model = F,
         x = F, y = F, contrasts = NULL, ...){
  thisCall <- match.call()
  if((callLen <- length(thisCall)) > 2){
    # extra arguments to add
    if(is.null(frame.eval <- formula$parent.frame))
      frame.eval <- sys.parent(1)
    if(callLen == 3 && names(thisCall)[[3]] == "weights"){
      # Special case when only weights are to be added.
      formula$fit$call$weights <- thisCall$weights
      formula$data$weights <- weights <- eval(weights, frame.eval)
      names(formula$data$weights) <- dimnames(formula$data$X)[[1]]
      # Need to add ... from original call to below
      formula$fit$fitexpr <-
        do.call("lmMakeFitExpr",
                c(list(fitOnly = F, weights = weights,
                       xlevels = formula$fit$xlevels), formula$fit$dotargs))
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

#
# lmMakeFitExpr
#
lmMakeFitExpr <- function(fitOnly, x = F, weights = NULL, xlevels = NULL, ...){
  #
  # Make fit expression for lm
  #
  # Arguments:
  # fitOnly  logical; if T the returned expression is ready for evaluation
  #          by lm.default.  If F a data extraction expression is prepended
  #          to the fit expression; this version of the fit expression
  #          typically becomes part of a model list to be evaluated from within
  #          lm.model.list.
  # x
  # weights
  # xlevels
  # ...      same meaning as inside lm.default.

  fitexpr <- c(Quote({fit           <- fcall;
		     fit$terms     <- Terms;
		     fit$call      <- fit.call;
		     fit$na.action <- na.action;}),
	      if(!is.null(xlevels))
	      Quote({attr(fit, "xlevels") <- xlevels}),
	      Quote({attr(fit, "na.message") <- na.message;
		     fit}))
  fcall <- ifelse1(length(weights),
    list(fcall=call("lm.wfit", x=if(x) Quote(X) else Quote(unset(X)),
	   y=Quote(Y), w=Quote(weights), method=Quote(method), ...)),
    list(fcall=call("lm.fit", x=if(x) Quote(X) else Quote(unset(X)),
	   y=Quote(Y), method=Quote(method), ...)))
  fitexpr <- substitute(fitexpr,fcall,evaluate=T)

  if(fitOnly) return(fitexpr)

  # Prepend a data extraction expression for evaluation within
  # lm.model.list.  All of the objects needed to perform the fit are now
  # part of the model.list, and this part of the expression extracts those
  # objects in preparation for doing the fit. Argument `formula' is the first
  # argument to lm.model.list, and refers to the model.list.
  c(Quote({X          <- formula$data$X;
           Y          <- formula$data$Y;
           weights    <- formula$data$weights;
           method     <- formula$fit$method;
           fit.call   <- formula$fit$call;
           Terms      <- formula$fit$Terms;
           na.action  <- formula$fit$na.action;
           na.message <- formula$fit$na.message;
           xlevels    <- formula$fit$xlevels}),
    fitexpr)
}

