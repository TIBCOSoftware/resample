#
# Functions in this file:
#	censorReg
# 	censorReg.default
# 	censorReg.model.list
# 	censorReg.dofit
#

#
# censorReg
#
# Note: full calling sequence retained for the generic definition, rather 
#       than simply function(formula, ...).
# Reason: predict.censorReg() fails otherwise (match.call statement in 
#       qftdist.dist() fails to identify `data' argument). 
#
censorReg <-
function(formula = formula(data), data = sys.parent(), weights = rep(1, n), 
	 truncation, subset, na.action, distribution = "weibull", 
	 threshold = 0, initial, fixed = NULL, control = NULL, 
	 model = F, x = F, y = F, method = NULL, ...)
  UseMethod("censorReg")

#
# censorReg.default
#
censorReg.default <-
function(formula = formula(data), data = sys.parent(), weights = rep(1, n), 
	 truncation, subset, na.action, distribution = "weibull", 
	 threshold = 0, initial, fixed = NULL, control = NULL, 
	 model = F, x = F, y = F, method = NULL, ...)
{
	Date <- date()
	this.call <- Call <- match.call()
	m <- match.call(expand = F)
	m$distribution <- m$model <- m$x <- m$y <- m$method <- m$... <- NULL
	m$initial <- m$control <- m$fixed <- m$threshold <- NULL
	TERMS <- Terms <- if(missing(data)) terms(formula, "strata") else terms(
			formula, "strata", data = data)
	m$formula <- Terms
	m[[1]] <- as.name("model.frame")
	m$drop.unused.levels <- T
#       See bug 20994
#	if(missing(data))
		m <- eval(m, sys.parent())
#	else m <- eval(m, data)
	YT <- model.extract(m, "truncation")
	dd <- make.distribution(distribution)
	Distribution <- dd$Distribution
	distribution <- dd$family[1]
	link <- dd$family[2]
	fixed <- c(dd$fixed, fixed)
	Y <- model.extract(m, "response")
	if(!is(Y, "censor"))
		Y <- as.censor(Y)
	n <- nrow(Y)
	censor.codes <- Y[, ncol(Y)]
	outCodes <- attr(Y, "Codes")$outCodes
	type <- attr(Y, "type")	#
	na.action <- attr(m, "na.action")

#
#
#  case weights
#
	casewt <- model.extract(m, "weights")
	if(is.null(casewt) || length(casewt) == 0)
		casewt <- rep(1, n)
	if(is.logical(threshold))
		if(threshold == F)
			threshold <- 0
#
#       X and strata
#
	strats <- attr(Terms, "specials")$strata
	if(length(strats)) {
		strat <- untangle.specials(Terms, "strata", 1)
		X <- model.matrix(Terms[ - strat$terms], m)
		cstrats <- m[strat$vars]
		istrats <- as.numeric(data.matrix(cstrats))
		cstrats <- levels(cstrats[, 1])[sort(unique(istrats))]
	}
	else {
		istrats <- rep(1, length = dim(Y)[1])
		cstrats <- strata(istrats)
		cstrats <- levels(cstrats)
		X <- model.matrix(Terms, m)
	}
	nvar <- dim(X)[2]	

#
#       fixed parameters
#
	if(missing(initial))
		initial <- NULL
	if(length(X) == 0) {
		asgn <- NULL
		contrasts <- list()
	}
	else {
		asgn <- attr(X, "assign")
		contrasts <- attr(X, "contrasts")
	}
	if(length(X) > 0)
		parameter.fixed <- rep(F, dim(X)[2] + 1)
	else parameter.fixed <- F
	ifix <- NULL
	if(length(fixed) > 0) {
		if(length(X) > 0)
			nm <- c(dimnames(X)[[2]], "scale")
		else nm <- "scale"
		ifix <- match(names(fixed), nm, nomatch = 0)
		if(any(ifix == 0))
			stop(paste("parameter(s)", paste(nm[ifix == 0]), 
				"in the fixed list are not valid for this model"
				))
		names(ifix) <- nm[ifix]
		parameter.fixed[ifix] <- T
	}

#
#  extract offset
#
	offset <- attr(Terms, "offset")
	if(!is.null(offset)){
	  offset <- as.numeric(m[[offset]])
	  is.offset <- T
	}
	else{
	  offset <- rep(0, n)
	  is.offset <- F
	}

	term.labs <- dimnames(attr(Terms, "factors"))[[1]][-1]
#	offsets <- attr(Terms, "offset")	#not used
#
# modify formula and subset in the case of a stratified fit so that the 
#   call stored on the fit object can reproduce the fit for each stratum.
#
	formula <- attr(Terms, "formula")
	variables <- attr(Terms, "variables")
	rcall <- function(name, expr)
	{
		if(length(expr) <= 1) {
			return(expr[[1]])
		}
		else {
			Call <- call(name, expr[[1]], expr[[2]])
			expr <- expr[-1]
			expr[[1]] <- Call
			rcall(name, expr)
		}
	}
	assign("rcall", rcall, frame = 1)
	on.exit(if(exists("rcall", frame = 1))
	    remove(as.character("rcall"), frame = 1))
	if(length(cstrats) > 1) {
		strata.call <- Call
		if(length(term.labs) > 1) {
#
#  	#  	Stratified fit with more than one term 
#
			strat.term <- attr(Terms, "specials")$strata
			if(!attr(TERMS, "intercept"))
				formula <- call("~", formula[[2]], call("+", 
				  substitute(-1), rcall("+", variables[ - 
				  strat.term][-1])))
			else formula <- call("~", formula[[2]], rcall("+", 
				  variables[ - strat.term][-1]))
			Terms <- terms(formula)	#
#
#	#	Stratified with one term
#
		}
		else {
			if(!attr(TERMS, "intercept"))
				stop(paste("Can't remove intercept when", attr(
				  TERMS, "term.labels"), 
				  "is the only term in the model"))
			formula <- call("~", formula[[2]], substitute(1))
			Terms <- terms(formula)
		}
		Call[["formula"]] <- formula	#
#
#	#     Now subsets
#
		subsets <- lapply(cstrats, function(x, y)
		call("==", y, x), attr(TERMS, "variables")[attr(TERMS, 
			"specials")$strata][[1]])
	}
	else subsets <- strata.call <- NULL
#	linkfun <- glm.links["link", link][[1]] # not used

	sd <- censorReg.distributions[[distribution]]
	sd.init <- sd$init
	sd.deviance <- sd$deviance


        #
        # censorRegMakeFitExp constructs an expression to perform the fit.
        # The expression will either be evaluated later in this routine, or if
        # method="model.list", returned as part of the model list for later
        # evaluation. 
        # 
	if(!is.null(method)){
	  if(method == "model.list"){# return the model list
   	    fitexpr <- censorRegMakeFitExpr(fitOnly = F)
  	    return(structure(list(data=list(X=X, Y=Y, YT=YT, casewt=casewt, 
  				    istrats=istrats, offset=offset, 
  				    censor.codes=censor.codes), 
  				  fit=list(fitexpr=fitexpr,
  				    cstrats=cstrats,
  				    is.offset=is.offset,
  				    subsets=subsets,
  				    na.action=na.action,
  				    distribution=distribution,
  				    Distribution=Distribution,
  				    threshold=threshold,
  				    initial=initial,
  				    fixed=fixed,
  				    control=control, 
  				    asgn=asgn,
  				    contrasts=contrasts,
  				    link=link,
  				    outCodes=outCodes,
  				    Terms=Terms,
  				    strata.call=strata.call,
  				    Date=Date,
  				    parameter.fixed=parameter.fixed,
  				    type=type,
  				    sd.init=sd.init,
  				    sd.deviance=sd.deviance,
  				    nvar=nvar,
  				    ifix=ifix,
  				    Call=Call,
  				    model=model,
  				    x=x,
  				    y=y),
  			    formula=formula,
  			    call=this.call),
  		       class="model.list"))
	  }
	  else
	    stop("unrecognized value for method")
	}
        fitexpr <- censorRegMakeFitExpr(fitOnly = T)
	#
        # Evaluate the fit expression
        #
	eval(fitexpr)
      }


#
# censorReg.model.list
#
censorReg.model.list <-
function(formula = formula(data), data = sys.parent(), weights = rep(1, n), 
	 truncation, subset, na.action, distribution = "weibull", 
	 threshold = 0, initial, fixed = NULL, control = NULL, 
	 model = F, x = F, y = F, method=NULL, ...){
  thisCall <- match.call()
  if((callLen <- length(thisCall)) > 2){
    # extra arguments to add
    if(is.null(frame.eval <- formula$parent.frame))
      frame.eval <- sys.parent(1)
    if(callLen == 3 && names(thisCall)[[3]] == "weights"){
      # Special case when only weights are to be added.
      formula$fit$Call$weights <- thisCall$weights
      formula$data$casewt <- eval(weights, frame.eval) 
      names(formula$data$casewt) <- dimnames(formula$data$X)[[1]]
      eval(formula$fit$fitexpr)
    }
    else{
      # evaluate the original call, appending new arguments
      formula$call$"method" <- NULL
      eval(c(formula$call, thisCall[-(1:2)]), local = frame.eval)
    }
  }
  else
    eval(formula$fit$fitexpr)
}
  
censorRegMakeFitExpr <- function(fitOnly){
  #
  # Make fit expression for censorReg
  #
  # Arguments:
  # fitOnly  logical; if T the returned expression is ready for evaluation
  #          by censorReg.default.  If F a data extraction expression is 
  #          prepended to the fit expression; this version of the fit
  #          expression typically becomes part of a model list to be evaluated
  #          from within censorReg.model.list.
  # ...      same meaning as inside censorReg.default.
  fitexpr <- Quote({censorReg.dofit(formula, X, Y, YT, casewt, 
                                    istrats, cstrats, is.offset, offset,
                                    subsets, censor.codes, na.action, 
                                    distribution, Distribution, threshold,
                                    initial, fixed, control, asgn, contrasts,
                                    link, outCodes, Terms, strata.call, Date,
                                    parameter.fixed, type, sd.init, sd.deviance,
                                    nvar, ifix, Call, model, x, y)})
  if(fitOnly) return(fitexpr)

  # Prepend a data extraction expression for evaluation within
  # censorReg.model.list.  All of the objects needed to perform the fit are now
  # part of the model.list, and this part of the expression extracts those
  # objects in preparation for doing the fit. Argument `formula' is the first
  # argument to censorReg.model.list, and refers to the model.list. 
  c(Quote({X           <- formula$data$X;
           Y           <- formula$data$Y;
           YT          <- formula$data$YT;
           casewt      <- formula$data$casewt;
           istrats     <- formula$data$istrats;
           offset      <- formula$data$offset;
           censor.codes<- formula$data$censor.codes;
           cstrats     <- formula$fit$cstrats;
           is.offset <- formula$fit$is.offset;
           subsets <- formula$fit$subsets;
           na.action <- formula$fit$na.action;
           distribution <- formula$fit$distribution;
           Distribution <- formula$fit$Distribution;
           threshold <- formula$fit$threshold;
           initial <- formula$fit$initial;
           fixed <- formula$fit$fixed;
           control <- formula$fit$control;
           asgn <- formula$fit$asgn;
           contrasts <- formula$fit$contrasts;
           link <- formula$fit$link;
           outCodes <- formula$fit$outCodes;
           Terms <- formula$fit$Terms;
           strata.call <- formula$fit$strata.call;
           Date <- formula$fit$Date;
           parameter.fixed <- formula$fit$parameter.fixed;
           type <- formula$fit$type;
           sd.init <- formula$fit$sd.init;
           sd.deviance <- formula$fit$sd.deviance;
           nvar <- formula$fit$nvar;
           ifix <- formula$fit$ifix;
           Call <- formula$fit$Call;
           model <- formula$fit$model;
           x <- formula$fit$x;
           y <- formula$fit$y
         }), 
    fitexpr)
}

#
#  censorReg.dofit
#
censorReg.dofit <-
function(formula, X, Y, YT, casewt, istrats, cstrats, 
	 is.offset, offset, subsets, censor.codes, na.action, 
	 distribution, Distribution, threshold, initial, fixed, 
	 control, asgn, contrasts, link, outCodes, Terms, 
	 strata.call, Date, parameter.fixed, type, sd.init, 
	 sd.deviance, nvar, ifix, Call, 
	 model, x, y){
#
#  Set up Y
#
	Y.cens <- censorReg.make.Y(Y, Terms, threshold, casewt, link)
	Y <- Y.cens$Y
	Threshold <- Y.cens$threshold	#

#
#  truncation
#
	if(!is.null(YT)) {
		if(dim(YT)[2] == 2) {
			pVec <- YT[, 2] == outCodes["event"]
			YT[, 1] <- YT[, 1] - Threshold
			YT[pVec, 1] <- 0
			if(any(YT[, 1] < 0))
				stop("truncation point less than zero detected"
				  )
			YT.Y <- cbind(YT[, 1], YT[, 1])
			YT.codes <- YT[, 2]
		}
		else if(dim(YT)[2] == 3) {
			pVec <- YT[, 3] == outCodes["event"]
			YT[, 1:2] <- YT[, 1:2] - Threshold
			YT[pVec, 1:2] <- 0
			if(any(YT[, 1:2] < 0))
				stop("truncation point less than zero detected"
				  )
			YT.Y <- YT[, 1:2]
			YT.codes <- YT[, 3]
		}
		censorReg.check.code(censor.codes, Y, YT.codes, YT.Y, outCodes)
	}
	else {
		YT.Y <- NULL
		YT.codes <- NULL
	}

#
#  Check Y against offset (Needs to be done each resample since Y values
#  are adjusted by a potentially different threshold each resample.  offset
#  will not change, but gets resampled.)
#
	if (is.offset){
	  if(any(Y[, 1] < offset)) {
	    stop("offset value greater than a lower time bound")
	  }
	  if(!is.null(YT)) {
	    pVec <- YT.codes != outCodes["event"]
	    if(any(YT.Y[pVec, 1] < offset[pVec]))
	      stop("offset value greater than the truncation lower bound")
	  }
	}
  
#
#  loop over all strata
#    length of loop is one for non stratified fit
#
	fits <- list()
	ustrats <- sort(unique(istrats))
	for(j in seq(length(ustrats))) {
		i <- ustrats[j]
		Yi <- Y[istrats == i,  , drop = F]
		censor.codesi <- censor.codes[istrats == i]
		Calli <- Call
		if(length(unique(istrats)) > 1) {
			subseti <- subsets[[j]]
			if(!is.null(Calli$subset))
				Calli$subset <- call("&", Calli$subset, subseti
				  )
			else Calli$subset <- subseti
		}
		if(length(X) == 0)
			Xi <- NULL
		else Xi <- X[istrats == i,  , drop = F]
		casewti <- casewt[istrats == i]
		offseti <- offset[istrats == i]
		tfix <- sd.init(1.4, fixed, initial)
		fix2 <- tfix[, 2] == 1
		tfix[!fix2, 1] <- NA
		parms <- tfix[, 1]
		if(!censorReg.good.data(censor.codesi, casewti)) {
#
#  no fit - insufficient data
#
			coef <- as.numeric(rep(NA, ncol(Xi)))
			names(coef) <- dimnames(Xi)[[2]]
			warning(paste("Skipped", cstrats[i], 
				"because of too few failures"))
			fiti <- list(coefficients = coef, residuals = rep(
				as.numeric(NA), length(Yi)), fitted.values = 
				rep(as.numeric(NA), length(Yi)), effects = rep(
				as.numeric(NA), length(Yi)), stzd.residuals = 
				rep(as.numeric(NA), length(Yi)), censor.codes
				 = censor.codesi, R = as.numeric(NA), rank = 0, 
				assign = asgn, df.residual = as.numeric(NA), 
				contrasts = contrasts, weights = rep(as.numeric(
				NA), length(Yi)), distribution = Distribution, 
				family = c(name = distribution, link = link, ""
				), linear.predictors = rep(as.numeric(NA), 
				length(Yi)), iter = 0, parms = parms, var = 
				matrix(as.numeric(NA), ncol = length(coef) + 1, 
				nrow = length(coef) + 1, dimnames = list(c(
				names(coef), "Log(scale)"), c(names(coef), 
				"Log(scale)"))), fixed = fix2, dresiduals = rep(
				as.numeric(NA), length(Yi)), df.total = sum(
				casewti), case.weights = casewti, n.censored = 
				sum(casewti[censor.codesi != outCodes["event"]]
				), deviance = as.numeric(NA), null.deviance = 
				as.numeric(NA), loglik = as.numeric(c(NA, NA)), 
				terms = Terms, formula = formula, strata.call
				 = strata.call, call = Calli, correlation = 
				matrix(as.numeric(NA), ncol = 2, nrow = 2, 
				dimnames = list(c("(Intercept)", "Log(scale)"), 
				c("(Intercept)", "Log(scale)"))), threshold = 
				as.numeric(NA), offset = offseti, scale = 1, 
				truncation = YT)
			if(length(na.action))
				fiti$na.action <- na.action
			oldClass(fiti) <- c("censorReg", "glm", "lm")
			attr(fiti, "date") <- Date
			fits[[i]] <- fiti
			next
		}
		sfit <- censorReg.mlest(Yi, censor.codesi, Xi, casewti, 
			distribution, link, parameter.fixed = parameter.fixed, 
			fixed = fixed, control = control, offset = offseti, 
			YT.Y = YT.Y, YT.codes = YT.codes)
		thresholdi <- censorReg.make.Y(censor(Yi + Threshold, 
			censor.codesi, type = type, inCodes = "1-4"), Terms, 
			threshold, casewti, link)$threshold
		parms.fixed <- c(F, parameter.fixed[length(parameter.fixed)])
		fix2p <- fixed["scale"]
		nfit <- censorReg.mlest(Yi, censor.codesi, NULL, casewti, 
			distribution, link, parameter.fixed = parms.fixed, 
			fixed = fix2p, control = control, offset = offseti, 
			YT.Y = YT.Y, YT.codes = YT.codes)
		if(link == "log")
			yi <- cbind(log(Yi), censor.codesi)
		else yi <- cbind(Yi, censor.codesi)
		tmp1 <- as.list(sfit$coefficients)
		names(tmp1) <- names(sfit$coefficients)
		parms <- sd.init(yi, fixed, tmp1)[, 1]
		fixed2 <- sd.init(yi, fixed, tmp1)[, 2] == 1
		deriv <- sfit$deriv[, 2:4]
		scale <- sfit$coefficients["scale"]
		deriv[, 2:3] <- deriv[, 2:3]/scale
		deriv[, 3] <-  - deriv[, 3]/scale
		if(link == "log")
			deriv[censor.codesi == outCodes["event"], 1] <- deriv[
				censor.codesi == outCodes["event"], 1] + yi[
				censor.codesi == outCodes["event"], 1]
		if(is.null(Xi)) {
			ny <- 1
			etai <- rep(sfit$coefficients[1], length = dim(Yi)[1])
		}
		else {
			ny <- dim(Xi)[2]
			etai <- Xi %*% sfit$coefficients[1:nvar]
		}
		ii <- match("scale", names(sfit$coefficients), nomatch = 0)
		if(ii != 0)
			parameter.fixed2 <- parameter.fixed[ - ii]
		else parameter.fixed2 <- parameter.fixed
		if(sum(parameter.fixed2) != 0)
			etai <- etai - Xi[, parameter.fixed2, drop = F] %*% (
				sfit$coefficients[1:nvar][parameter.fixed2])
		fiti <- censorReg.wfit(Xi, etai + deriv[, 2]/abs(deriv[, 3]), 
			abs(deriv[, 3]), "qr", parameter.fixed = 
			parameter.fixed2)
		if(sum(parameter.fixed2) != 0) {
			tmp <- double(length(parameter.fixed2))
			tmp[!parameter.fixed2] <- fiti$coefficients
			tmp[parameter.fixed2] <- (sfit$coefficients[1:nvar][
				parameter.fixed2])
			names(tmp) <- names(sfit$coefficients)[1:nvar]
			fiti$coefficients <- tmp
		}
		fiti$weights <-  - deriv[, 3]
		fiti$residuals <-  - deriv[, 2]/deriv[, 3]
		fiti$stzd.residuals <- sfit$residuals
		fiti$censor.codes <- censor.codesi
		fiti$type.censoring <- type
		fiti$contrasts <- contrasts
		ifun <- glm.links["inverse", link][[1]]
		fiti$fitted.values <- ifun(fiti$fitted.values)
		if(is.null(casewti)) {
			fiti$df.residual <- dim(Yi)[1] - (nvar + 1) + sum(
				parameter.fixed)
			fiti$df.total <- dim(Yi)[1]
			fiti$n.censored <- dim(Yi[censor.codesi != outCodes[
				"event"],  ])[1]
		}
		else {
			fiti$df.residual <- sum(casewti) - (nvar + 1) + sum(
				parameter.fixed)
			fiti$df.total <- sum(casewti)
			fiti$n.censored <- sum(casewti[censor.codesi != 
				outCodes["event"]])
		}
		fiti$family <- c(name = distribution, link = link, "")
		fiti$distribution <- Distribution
		rnamesi <- dimnames(Yi)[[1]]
		if(!is.null(rnamesi) && length(rnamesi))
			names(etai) <- rnamesi
		fiti$linear.predictors <- etai
		fiti$iter <- sfit$iter
		fiti$iter <- NA
		fiti$parms <- parms
		var <- sfit$cov
		matnames <- dimnames(var)[[1]]
		itmp <- pmatch("scale", matnames, nomatch = 0)
		if(itmp > 0) {
			matnames[itmp] <- "Log(scale)"
			var[, "scale"] <- var[, "scale"]/scale
			var["scale",  ] <- var["scale",  ]/scale
			dimnames(var) <- list(matnames, matnames)
		}
		if(length(ifix) > 0)
			fiti$var <- var[ - ifix,  - ifix]
		else fiti$var <- var
		if(!is.null(fixed$scale))
			fixed$scale <- NULL
		if(length(fixed) > 0) {
			tmp <- rep(T, length(fixed))
			names(tmp) <- names(fixed)
			fixed2 <- c(fixed2, tmp)
		}
		fiti$fixed <- fixed2
		deviances <- deriv[, 1]
		dev <- sd.deviance(yi, parms, deviances)
		fiti$dresiduals <- sign(fiti$residuals) * sqrt(dev)
		fiti$deviance <- sum(dev)
		fiti$null.deviance <- fiti$deviance + 2 * (sfit$loglike - nfit$
			loglike)
		fiti$loglik <- c(nfit$loglike, sfit$loglike)
		fiti$first.derivative <- sfit$derivative
		fiti$case.weights <- casewti
		fiti$scale <- scale
		dimnames(sfit$correlation) <- list(matnames, matnames)
		if(length(ifix) > 0)
			fiti$correlation <- sfit$correlation[ - ifix,  - ifix]
		else fiti$correlation <- sfit$correlation
		if(length(na.action))
			fiti$na.action <- na.action
		fiti$terms <- Terms
		fiti$formula <- as.vector(attr(Terms, "formula"))
		fiti$strata.call <- strata.call
		fiti$call <- Calli
		if(model)
			fiti$model <- m[istrats == i,  , drop = F]
		if(x)
			fiti$x <- Xi
		if(y) {
			tmp <- attr(Terms, "variables")[[1]]
			if(length(tmp) == 1)
				yi <- yi[, -2, drop = F]
			else {
				if(mode(tmp) == "call" | mode(tmp) == 
				  "expression" | mode(tmp) == "name")
				  dimnames(yi)[[2]][dim(Yi)[2] + 1] <- deparse(
				    tmp[[dim(Yi)[2] + 2]])
				else dimnames(yi)[[2]][dim(Yi)[2] + 1] <- 
				    as.character(tmp[[dim(Yi)[2] + 2]])
			}
			fiti$y <- yi
		}
		fiti$threshold <- thresholdi
		fiti$offset <- offseti
		fiti$truncation <- YT
		oldClass(fiti) <- c("censorReg", "glm", "lm")
		attr(fiti, "date") <- Date
		fits[[i]] <- fiti
	}
	fits <- fits[!sapply(fits, is.null)]
	if(length(cstrats) == 1) {
		attr(fits[[1]], "date") <- Date
		return(fits[[1]])
	}
	else {
		names(fits) <- cstrats
		oldClass(fits) <- c("censorRegList", "censorReg")
		return(fits)
	}
      }

# See where most of the time is spent:  looks like calls to censorReg.mlest
# (sometimes the first, sometimes 2nd).
#

# a1
# fit <- censorReg(censor(days, event) ~ voltage, data = capacitor, 
#                 distribution = "loggaussian", threshold = "Lin")
#[1] 0.13 0.00
#[1] "Entering censorReg.dofit"
#[1] 0.19 0.00
#[1] "Before censorReg.mlest (sfit)"
#[1] 0.38 0.01
#[1] "Before censorReg.mlest (nfit)"
#[1] 0.17 0.00
#[1] "After censorReg.mlest (nfit)"
#[1] 0.16 0.00
#[1] "Done"

# a14
# fit <- censorReg(censor(days, event) ~ strata(voltage), data = capacitor2, 
#                  weights = weights, distribution = "weibull", 
#                  threshold = "Lin")
#[1] 0.25 0.00
#[1] "Entering censorReg.dofit"
#[1] 0.25 0.00
#[1] "Before censorReg.mlest (sfit)"
#[1] 0.39 0.00
#[1] "Before censorReg.mlest (nfit)"
#[1] 0.13 0.00
#[1] "After censorReg.mlest (nfit)"
#[1] 0.15 0.00
#[1] "Before censorReg.mlest (sfit)"
#[1] 0.37 0.00
#[1] "Before censorReg.mlest (nfit)"
#[1] 0.12 0.00
#[1] "After censorReg.mlest (nfit)"
#[1] 0.15 0.00
#[1] "Before censorReg.mlest (sfit)"
#[1] 0.39 0.00
#[1] "Before censorReg.mlest (nfit)"
#[1] 0.13 0.00
#[1] "After censorReg.mlest (nfit)"
#[1] 0.12 0.00
#[1] "Done"

# d0
# fit <- censorReg(censor(time, lower, status, type = "interval") ~ 1, 
#                  data = cwc, threshold=T)
#[1] 0.17 0.00
#[1] "Entering censorReg.dofit"
#[1] 0.05 0.00
#[1] "Before censorReg.mlest (sfit)"
#[1] 0.57 0.01
#[1] "Before censorReg.mlest (nfit)"
#[1] 1.49 0.00
#[1] "After censorReg.mlest (nfit)"
#[1] 0.18 0.00
#[1] "Done"

# trunc1
# fit <- censorReg(censor(time, tlower, status, type = "interval") ~ 1, 
#                  data = trunc.data, threshold = 0, truncation = 
#	           censor(lower, upper, tstatus, type = "interval"))
#[1] 0.19 0.00
#[1] "Entering censorReg.dofit"
#[1] 0.1 0.0
#[1] "Before censorReg.mlest (sfit)"
#[1] 0.19 0.00
#[1] "Before censorReg.mlest (nfit)"
#[1] 0.13 0.00
#[1] "After censorReg.mlest (nfit)"
#[1] 0.11 0.00
#[1] "Done"

# offset1
# fit <- censorReg(censor(time, tlower, status, type = "interval") ~ 1 + 
#                  offset(offseti), data = trunc.data, threshold = 0, 
#                  truncation = censor(lower, upper, tstatus, 
#	          type = "interval"))
#[1] 0.22 0.00
#[1] "Entering censorReg.dofit"
#[1] 0.11 0.00
#[1] "Before censorReg.mlest (sfit)"
#[1] 0.19 0.00
#[1] "Before censorReg.mlest (nfit)"
#[1] 0.12 0.00
#[1] "After censorReg.mlest (nfit)"
#[1] 0.13 0.00
#[1] "Done"


