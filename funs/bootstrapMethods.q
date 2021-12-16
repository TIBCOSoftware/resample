#
# Methods and functions for bootstrapping fitted objects from
# modeling functions:
#
# bootstrap.lm
# bootstrap.glm
# bootstrap.censorReg
# jackknife.lm
# jackknife.glm
# jackknife.censorReg
#
# resampMakeFitObjResults
# resampMakeFitObjResults.censorReg
# resampMakeNewStat
#
#
# Most methods call resampMakeFitObjResults, which creates and evaluates a
# new call to the default method of bootstrap or jackknife.  One idea specific
# to bootstrap is the notion of bootstrapping residuals with lm, so that case
# is handled separately in bootstrap.lm.  In addition, censorReg can
# produce objects of class censorReg or censorRegList, so special
# consideration there is provided by resampMakeFitObjResults.censorReg.
#

#################################################################
# bootstrap.lm
#################################################################
bootstrap.lm <-
function(data, statistic, B = 1000, args.stat = NULL,
	 group = NULL,
	 subject = NULL,
	 sampler = samp.bootstrap,
	 seed = .Random.seed,
	 sampler.prob = NULL,
	 sampler.args = NULL,
	 sampler.args.group = NULL,
	 resampleColumns = NULL,
	 label = NULL,
	 statisticNames = NULL,
	 block.size = min(100, B),
	 trace = resampleOptions()$trace,
	 assign.frame1 = F,
	 save.indices = F,
	 save.group = NULL,
	 save.subject = NULL,
	 statistic.is.random = NULL,
	 group.order.matters = T,
         order.matters = NULL,
	 seed.statistic = 500,
	 L = NULL,
	 model.mat = NULL,
	 lmsampler = "observations", ...){

  this.call <- match.call()
  originalCall <- this.call
  lmsampler <- match.arg(casefold(lmsampler),
                         c("observations", "residuals", "wild", "wild-as"))

  # If statistic is coef, then change to coef.default
  if(identical(this.call$statistic, as.name("coef")))
    this.call$statistic <- as.name("coef.default")
  else if(is.call(this.call$statistic) &&
	  identical(this.call$statistic[[1]], as.name("coef")))
    this.call$statistic[[1]] <- as.name("coef.default")

  newcall <- this.call
  newcall$lmsampler <- NULL             # used only by bootstrap.lm

  if(lmsampler == "observations"){
    boot <- resampMakeFitObjResults(resampCall = newcall, fitCall = data$call,
                                    frame.eval = sys.parent(),
                                    model.method = if(class(data) == "lm") "lm")
    boot$actual.calls[[1]] <- originalCall
    boot$defaultLabel = resampMakeLabel("bootstrap",
      originalCall$data, originalCall$statistic)
    return(boot)
  }
  
  # Rest of this function is for "residuals" "wild" or "wild-as"
  #   For "residuals", new call will be
  #        bootstrap(R, stat(lm(resampChangeY.model.list(ml,preds+data))))
  #   where R & data are the residuals from the lm fit, preds are the
  #   predicted values, and ml is the model.list.
  #   For "wild",
  #        bootstrap(c(1,-1), ... preds+data*R), additional args)
  #   where data = c(1,-1)
  #         args.sampler: size=n
  #         observed.indices = rep(1,n), ...)
  #   For "wild-as" (asymmetric)
  #         data = c(-(sqrt(5)-1)/2, (sqrt(5)+1)/2)
  #         args.sampler: size=n, prob=c(p,1-p) with p=(sqrt(5)+1)/(2*sqrt(5))
  if(!is.null(newcall$L)){
    newcall$L <- "none"
    warning("L is not supported when resampling residuals")
  }

  # create the model list by modifying the call creating the lm object
  mlcall <- data$call        # gets call even if data is a name, since
  # data is evaluated (newcall$data is not)
  if(is.null(mlcall$method))
    mlcall$method <- "model.list"
  else                     # in case there is a fitting method present
    mlcall$method <- c(eval(mlcall$method,sys.parent()),"model.list")
  ml <- eval(mlcall,sys.parent())

  # compute predicted values
  preds <- unname(predict(data))

  # create new statistic: note that all data is passed directly in call
  # This is a list with one component (to be added to a call later)
  statarg <- list(call("lm", call("resampChangeY.model.list", ml,
        newY = call("+", preds, 
          ifelse1(lmsampler == "residuals",
                  Quote(data),
                  call("*", Quote(data), unname(resid(data))))))))
    
  names(statarg) <- ifelse1(is.name(newcall$data), newcall$data, "data")
  newcall$statistic <- resampMakeNewStat(statistic = newcall$statistic,
                                         newstatarg = statarg,
                                         args.stat = args.stat,
                                         frame.eval = sys.parent())
  # data component of new call
  if(lmsampler == "residuals") {
    newcall$data <- unname(resid(data))
  } else { # wild
    n <- numRows(preds)
    newcall$observed.indices <- rep(1, n)
    if(lmsampler == "wild"){
      newcall$data <- c(1,-1)
      newcall$sampler.args <- c(newcall$args.sampler, list(size = n))
    } else {
      newcall$data <- c(-(sqrt(5)-1)/2, (sqrt(5)+1)/2)
      p <- (sqrt(5)+1)/(2*sqrt(5))
      newcall$sampler.args <- c(newcall$args.sampler,
                                list(size = n, prob = c(p, 1-p)))
    }
  }
  # make sure call is to bootstrap, not bootstrap.lm
  newcall[[1]] <- Quote(bootstrap)
  # evaluate the call and set call component to be original call
  result <- eval(newcall,sys.parent())
  result$call <- this.call
  result$defaultLabel <- resampMakeLabel(paste("bootstrap", lmsampler),
    originalCall$data, originalCall$statistic)
  result$order.matters <- switch(lmsampler,
                                 residuals = "resampling residuals",
                                 "wild bootstrap")
  return(result)
}
# use match.arg instead of pmatch and switch.
# add "wild" and "wild-as" (asymmtric)


#################################################################
# bootstrap.glm
#################################################################
bootstrap.glm <-
function(data, statistic, B = 1000, args.stat = NULL,
	 group = NULL,
	 subject = NULL,
	 sampler = samp.bootstrap,
	 seed = .Random.seed,
	 sampler.prob = NULL,
	 sampler.args = NULL,
	 sampler.args.group = NULL,
	 resampleColumns = NULL,
	 label = NULL,
	 statisticNames = NULL,
	 block.size = min(100, B),
	 trace = resampleOptions()$trace,
	 assign.frame1 = F,
	 save.indices = F,
	 save.group = NULL,
	 save.subject = NULL,
	 statistic.is.random = NULL,
	 group.order.matters = T,
         order.matters = NULL,
	 seed.statistic = 500,
	 L = NULL,
	 model.mat = NULL,
	 ...){
  resampMakeFitObjResults(resampCall = match.call(), fitCall = data$call,
			  frame.eval = sys.parent(),
			  model.method = if(class(data) == "glm") "glm")
}


#################################################################
# bootstrap.censorReg
#################################################################
bootstrap.censorReg <-
function(data, statistic, B = 1000, args.stat = NULL,
	 group = T,
	 subject = NULL,
	 sampler = samp.bootstrap,
	 seed = .Random.seed,
	 sampler.prob = NULL,
	 sampler.args = NULL,
	 sampler.args.group = NULL,
	 resampleColumns = NULL,
	 label = NULL,
	 statisticNames = NULL,
	 block.size = min(100, B),
	 trace = resampleOptions()$trace,
	 assign.frame1 = F,
	 save.indices = F,
	 save.group = NULL,
	 save.subject = NULL,
	 statistic.is.random = NULL,
	 group.order.matters = T,
         order.matters = NULL,
	 seed.statistic = 500,
	 L = NULL,
	 model.mat = NULL,
	 ...){
  this.call <- match.call()
  # The default value of group must be forced into the call
  if(missing(group))
    this.call$group <- group

  if(is.element(class(data), c("censorReg", "censorRegList")))
    resampMakeFitObjResults.censorReg(resampCall = this.call,
				      frame.eval = sys.parent())
  else
    resampMakeFitObjResults(resampCall = this.call, fitCall = data$call,
			    frame.eval = sys.parent(), model.method = NULL)
}

#################################################################
# jackknife.lm
#################################################################
jackknife.lm <-
function(data, statistic, args.stat = NULL,
	 group = NULL, subject = NULL,
	 label = NULL,
	 statisticNames = NULL,
	 seed = .Random.seed,
	 group.size = 1, assign.frame1 = F,
	 save.group = NULL, save.subject = NULL, ...){
  resampMakeFitObjResults(resampCall = match.call(), fitCall = data$call,
			  frame.eval = sys.parent(),
			  model.method = if(class(data) == "lm") "lm")
}

#################################################################
# jackknife.glm
#################################################################
jackknife.glm <-
function(data, statistic, args.stat = NULL,
	 group = NULL, subject = NULL,
	 label = NULL,
	 statisticNames = NULL,
	 seed = .Random.seed,
	 group.size = 1, assign.frame1 = F,
	 save.group = NULL, save.subject = NULL, ...){
  resampMakeFitObjResults(resampCall = match.call(), fitCall = data$call,
			  frame.eval = sys.parent(),
			  model.method = if(class(data) == "glm") "glm")
}

#################################################################
# jackknife.censorReg
#################################################################
jackknife.censorReg <-
function(data, statistic, args.stat = NULL,
	 group = NULL, subject = NULL,
	 label = NULL,
	 statisticNames = NULL,
	 seed = .Random.seed,
	 group.size = 1, assign.frame1 = F,
	 save.group = NULL, save.subject = NULL, ...){
  this.call <- match.call()

  if(is.element(class(data), c("censorReg", "censorRegList")))
    resampMakeFitObjResults.censorReg(resampCall = this.call,
				      frame.eval = sys.parent())
  else
    resampMakeFitObjResults(resampCall = this.call, fitCall = data$call,
			    frame.eval = sys.parent(), model.method = NULL)
}


################################################################
# resampMakeFitObjResults
################################################################
resampMakeFitObjResults <-
function(resampCall, fitCall, frame.eval, model.method = NULL){
  # Various bootstrap and jackknife methods take as input calls like:
  #    bootstrap(fit, ...)
  #    jackknife(fit, ...)
  # then call this function to create a new call to bootstrap or jackknife.
  #
  # If model.method is not null, the new call is of the form:
  #    resampFUN(model.method(...,method="model.list"),
  #              stat(model.method(data), ...))
  # where resampFUN is the function called in resampCall.
  # I.e., data is replaced by a call to a modeling function that returns
  # a model.list, and the statistic acts on new argument `model.method(data)'.
  # For example, this
  #    fit <- lm(Mileage ~ Weight, data = fuel.frame)
  #    bootstrap(fit, coef(fit), B=500)
  # becomes
  #    bootstrap(lm(Mileage ~ Weight, data = fuel.frame, method="model.list"),
  #              coef(lm(data)), B=500)
  # (bootstrap.lm calls this function with model.method="lm").
  #
  # If model.method = NULL, the new call is of the form:
  #    resampFUN(fitCall$data, stat(fitCall), ...)
  # I.e., the data is taken from the data component of fitCall,
  # and the statistic acts on the fit.
  # For example, this
  #    fit <- ltsreg(ozone~wind+radiation+temperature, data=small.air)
  #    bootstrap(fit, coef)
  # becomes
  #    bootstrap(small.air, coef(ltsreg(ozone~wind+radiation+temperature,
  #                                     data=small.air))
  # (there is no bootstrap method for ltsreg; bootstrap.default
  # calls this function with model.method = F).

  origCall <- resampCall

  # Check that the resample function is recognized, and make the call
  # to the generic (needed in case bootstrap.lm is called directly, for
  # instance).
  resampFuns <- c("bootstrap", "jackknife")
  if(any(resampFunInd <- !is.na(pmatch(resampFuns, resampCall[[1]]))))
    resampCall[[1]] <- as.name(resampFuns[resampFunInd])
  else
    stop("unknown resample function")
  stat      <- resampCall$statistic
  args.stat <- as.list(resampCall$args.stat)[-1]

  # newData and statarg are different depending on whether or not there is
  # a model.method.
  if(is.null(model.method)){
    statarg <- list(fitCall)
    newData <- fitCall$data
  }
  else{
    statarg        <- list(call(model.method, Quote(data)))
    newData        <- fitCall
    newData[[1]]   <- as.name(model.method)
    newData$method <-
      if(is.null(newData$method)) "model.list"
      else c(eval(newData$method,frame.eval),"model.list")
  }

  # statarg will be substituted for the name of the data,  or for "data"
  # if the data is unnamed.
  names(statarg)  <- if(is.name(resampCall$data)) resampCall$data
                     else "data"
  resampCall$data <- newData

  # Create a new statistic acting on statarg.
  resampCall$statistic <- resampMakeNewStat(statistic = stat,
					   newstatarg = statarg,
					   args.stat = args.stat,
					   frame.eval = frame.eval)
  # Evaluate the call and assign calls to results
  result <- eval(resampCall, frame.eval)
  result$actual.calls <- list(origCall)
  result$call <- resampCall
  result
}
# change model.method to hold name of generic fitting routine.  This
# addresses the cases where match.call returns lm.default, for instance,
# even if lm is called (see Bug 21288).  I was using fitCall[[1]] as the
# value for statarg, which was causing errors when it was lm.default It
# needs to be lm or lm.model.list.  Now we use model.method. (RT, 7/17/01)


#################################################################
# resampMakeFitObjResults.censorReg
#################################################################
resampMakeFitObjResults.censorReg <-
function(resampCall, frame.eval){

  origCall        <- resampCall
  resampFuns      <- c("bootstrap", "jackknife")
  resampCall[[1]] <- as.name(resampFuns[!is.na(pmatch(resampFuns,
						      resampCall[[1]]))])
  stat            <- resampCall$statistic
  statarg         <- list(Quote(censorReg(data)))
  args.stat       <- as.list(resampCall$args.stat)[-1]

  # If the data is named, replace it with the actual call, and change the
  # method to "model.list" (or add "model.list" to the vector).  If data
  # is of class "censorRegList", the call comes from the "strata.call"
  # component of any element of the list (we use the first).
  # statarg will be substituted for the name of the data, or for "data"
  # if the data is unnamed.
  if(is.name(resampCall$data)){
    names(statarg) <- resampCall$data
    cr.obj <- eval(resampCall$data, frame.eval)
    resampCall$data <- if(class(cr.obj)=="censorRegList") cr.obj[[1]]$strata.call
                      else cr.obj$call
  }
  else
    names(statarg) <- "data"
  resampCall$data$method <-
    if(is.null(resampCall$data$method)) "model.list"
    else c(eval(resampCall$data$method, frame.eval),"model.list")

  # Create a new statistic acting on statarg.
  resampCall$statistic <- resampMakeNewStat(statistic = stat,
					   newstatarg = statarg,
					   args.stat = args.stat,
					   frame.eval = frame.eval)

  # Group argument (only applies to bootstrap)
  if(resampCall[[1]]=="bootstrap"){
    group <- eval(resampCall$group, frame.eval)
    if(is.logical(group) && length(group)==1){
      if(group){
        # Evaluate the model list: component data$istrats indicates
        # strata, and is appropriate for the group argument to bootstrap.
        # The original call needs to have this group argument, as opposed
        # to T or F, because it gets used in limits.bca.
	ml <- eval(resampCall$data, frame.eval)
	origCall$group <- resampCall$group <- if(length(ml$fit$cstrats) > 1)
	  ml$data$istrats else NULL
      }
      else
	origCall$group <- resampCall$group <- NULL
    }
  }
  # Evaluate the call and assign the calls to the results
  result <- eval(resampCall, frame.eval)
  result$actual.calls <- list(origCall)
  result$call <- resampCall
  result
}

#################################################################
# resampMakeNewStat
#################################################################
resampMakeNewStat <-
function(statistic, newstatarg, args.stat, modifyFunctions = T, frame.eval){
  # Create a new statistic with argument newstatarg.
  # Result will be an expression, unless modifyFunctions = F and the statistic
  # is (or evaluates to) a function, in which case just return the unevaluated,
  # un-modified statistic.
  #
  # If statistic is a function or name, then args.stat contains extra arguments,
  # so the new expression includes those arguments. Technically it would be
  # appropriate to set args.stat to NULL in the final call to bootstrap in
  # this case, but it doesn't hurt to leave args.stat as is.
  #
  # If statistic is an expression then args.stat is a list of objects to include
  # in the frame of evaluation, so nothing needs to be done: the call to
  # bootstrap will still have the same args.stat argument, and the new
  # expression will be evaluated using that list.

  # In case, for instance, statistic=temp, where temp=Quote(coef) or
  # Quote(coef(fit)), or some other expression, then make statistic
  # that expression.
  if(mode(statistic) == "name" &&
     mode(newstat <- eval(statistic,frame.eval)) != "function")
    statistic <- newstat

  # Different cases based on mode of the original statistic.
  switch(mode(statistic),

	 name =,     # e.g., coef

	 function = {# e.g., function(x, ...) predict(x, ...)
	   if(modifyFunctions)
             as.call(c(call(statistic, newstatarg[[1]]), args.stat))
           else
             statistic
	 },

	 call =,     # e.g., coef(glmobj), or coef(data)[1]-coef(data)[2])

	 "{" =,      # e.g., {result<-try(coef(fit));
	             #        if(is(result,"Error")) rep(NA,4) else result}

	 "(" = {     # e.g., ((coef(fit)))
	   substitute(statistic, newstatarg, evaluate=T)
	 },

	 stop(paste("Unsupported mode for statistic:", mode(statistic)))
	 )
}
