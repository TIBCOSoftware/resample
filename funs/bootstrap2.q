#
# bootstrap2: bootstrap the difference or ratio
#             of statistic computed on two samples
#
bootstrap2 <-
  function(data, statistic, treatment, data2 = NULL, B = 1000, ratio = F,
	   group = NULL, subject = NULL, seed = .Random.seed,
	   trace = resampleOptions()$trace,
           save.group = NULL, save.subject = NULL,
	   save.treatment = NULL, L = NULL, twoSample.args=NULL, ...){
  # Call this with either of:
  #   bootstrap2(data, statistic, treatment, ...)
  #   bootstrap2(data, statistic, data2, ...)
  if(missing(statistic))
    stop("Argument statistic is missing with no default.")

  new.call <- this.call <- match.call()
  parent.frame <- sys.parent()

  # If data2 is used, construct a treatment vector to match, and set the
  # new data argument to concatenate data and data2.
  if(n2 <- numRows(data2)){
    # Note -- n2 is used as a flag below, is nonzero if data2 supplied.
    n1 <- numRows(data)
    n <- n1 + n2
    treatment <- rep(1:2, c(n1, n2))
    data.expr <-
      if((length(dim(data))))
        Quote(rbind(data, data2))
      else
        Quote(c(data, data2))
    data.expr <- substitute(data.expr, list(data = this.call$data,
                                            data2 = this.call$data2), eval = T)
    new.call$data2 <- NULL
  }
  else {
    n <- numRows(data)
    data.expr <- this.call$data
  }

  if(is.df <- is.data.frame(data)){
    # Look for treatment and group as variables in the data frame.
    # Also allow expressions like `region == "south"'.
    rectangular <- T
    if(is.element("data", names(data)))
      stop("The name \"data\" is reserved for internal use, and must not ",
           "be used as a name in the data.")
    if(!n2)
      treatment <- resampMakeArgument(this.call$treatment, data, parent.frame)
    if(!missing(group))
      group <- resampMakeArgument(this.call$group,
				  if(n2) rbind(data, data2) else data,
				  parent.frame)
    if(!missing(subject))
      subject <- resampMakeArgument(this.call$subject,
				    if(n2) rbind(data, data2) else data,
				    parent.frame)
  }
  else{
    rectangular <- is.matrix(data)
    if(rectangular && length(cnames <- dimnames(data)[[2]]) &&
       is.element("data", cnames))
         stop("The name \"data\" is reserved for internal use, and must not ",
              "be used as a name in the data.")
  }

  # New statistic should use "data" instead of the data name, since
  # the new data will be a subscripted expression.
  if(is.name(data.expr)){
    # Check that the original statistic does not use "data" when the
    # original data argument is a name.
    if(is.element("data", all.names(this.call$statistic)))
      stop("The name \"data\" can be used in statistic only when the ",
           "data argument is an expresssion or when data2 is used.")
    newstatarg <- list(Quote(data))
    names(newstatarg) <- data.expr
    new.call$statistic <- resampMakeNewStat(statistic = this.call$statistic,
                                            newstatarg = newstatarg,
                                            args.stat = NULL,
                                            modifyFunctions = F,
                                            frame.eval = sys.parent())
  }

  # Get the treatment subscript vectors
  u <- unique(treatment)
  if(length(u) != 2)
    stop("treatment must have exactly two unique values.")
  treat1 <- treatment == u[1]
  treatment.both <- list(treat1, !treat1)
  rle.treat1 <- rle(treat1)
  if(length(rle.treat1$values) * 2 + 10 < length(treat1)){
    treatment.both[[1]] <- substitute(rep(A,B),
                                      list(A = rle.treat1$values,
                                           B = rle.treat1$lengths))
    treatment.both[[2]] <- substitute(rep(A,B),
                                      list(A = !rle.treat1$values,
                                           B = rle.treat1$lengths))
  }

  # Set up group argument
  if(byGroup <- length(group)){
    if(byGroup != n)
      stop("group must have length equal to the total number of observations")
    # Factors are a problem when we pass a proper subset of them to
    # bootstrap, so just convert them to character.
    if(is.factor(group)) group <- as.character(group)
    treatment.group <- list(group[treat1], group[!treat1])
  }

  # Set up subject argument
  if(bySubject <- length(subject)){
    if(is.factor(subject)) subject <- as.character(subject)
    if(notNested(lapply(split(subject, treatment), unique)))
      stop("some values of `subject' occur in both treatments; ",
           "subject must be nested within the `treatment' argument")
    treatment.subject <- list(subject[treat1], subject[!treat1])
  }

  # Set up L argument, if supplied as a numerical object.
  # In this case we expect that L has the signs reversed for
  # the second groups, e.g. bootstrap2(x, mean, treatment, L=x)
  # would be incorrect, would need to be
  # L*ifelse(treatment==treatment[1], 1, -1)
  # This is ugly, but want to prevent contradiction between L passed
  # and L returned.
  # Do not encourage supplying a numerical object.
  # In the case of mean, the sub-bootstraps will compute L anyway.
  if(L.is.numeric <- is.numeric(L)){
    if(!is.matrix(L))
      L <- as.matrix(L)
    if(is.data.frame(L))
      L <- numerical.matrix(L)
    if(bySubject){
      # minimal support for this case; require that L have row names
      # that match subjects, and assume (without checking) that L
      # is normalized to mean 0 by group
      if(!identical(dimnames(L)[[1]], as.character(unique(subject))))
	stop("L must have row names matching the sorted subjects")
      treatment.L <- list(L[subject[treat1],,drop=F],
			 L[subject[!treat1],,drop=F])
    }
    else if(numRows(L) != n){
      warning("L has the wrong number of rows; will delete L")
      L <- NULL
      L.is.numeric <- F
    }
    else{
      L <- subtractMeans(L, treatment)
      treatment.L <- list(L[treat1,,drop=F], -L[!treat1,,drop=F])
    }
  }
  if(ratio && L.is.numeric){
    # Things get even more complicated when bootstrapping a ratio;
    # would need to know the two observed values to translate between L's
    # for individual samples and the overall L.  That information
    # is not available.  Our convention will be that the L that
    # is passed is correct for the two samples (aside from +-).
    # That means that L will be modified below to be correct for the ratio.
    warning("You supplied a numerical L when bootstrapping a ratio.",
	    "I'll assume that L is scaled correctly for the two",
	    "individual samples; I'll modify it below to be correct",
	    "for the ratio.")
  }

  # Each call to bootstrap will look something like
  #
  #   bootstrap(data[treatment,drop=F], statistic,
  #             group = group[treatment], L = (depends), ...)
  #
  # where treatment switches between the two groups between calls. The data
  # argument is a subscripted expression, not the actual subscripted
  # data.  This is to keep the call, which is essentially stored 3 times in
  # the returned object, manageable. The actual treatment and group vectors,
  # however, are used in the call (an attempt at using the expression
  # "[treat == unique(treat)[i]]" failed when treatment was a name of a
  # variable of data when data was a data frame -- this was getting hairy).
  new.call[[1]] <- as.name("bootstrap")
  new.call$treatment <- NULL
  new.call$data2 <- NULL
  new.call$ratio <- NULL
  new.call$twoSample.args <- NULL
  boot <- vector("list", 2)
  names(boot) <- u
  newdata <- substitute(if(rectangular) Quote((data)[treat,,drop=F])
                        else            Quote((data)[treat,drop=F]),
                        list(data = data.expr), eval = T)
  for(i in 1:2){
    i.call <- new.call
    i.call$data <- substitute(newdata, list(treat = treatment.both[[i]]),
                                eval = T)
    if(byGroup) i.call$group <- treatment.group[[i]]
    if(bySubject) i.call$subject <- treatment.subject[[i]]
    if(L.is.numeric) i.call$L <- treatment.L[[i]]
    if(trace) cat("Resampling data for treatment group", i, "of 2.\n")
    thisSample.args <- twoSample.args[[i]]
    if(length(thisSample.args))  # use loop to preserve class "call"
      for(j in names(thisSample.args))
        i.call[[j]] <- thisSample.args[[j]]
    boot[[i]] <- eval(i.call, sys.parent(1))
    new.call$seed <- NULL # so 2nd call picks up where 1st leaves off
  }

  # In case something goes wrong during combination:
  on.exit({
    cat("Problem occured while combining results for two samples.\n",
	"Saving results in .bootstrap2.partial.results\n")
    assign(".bootstrap2.partial.results", boot, where = 1, immediate = T)
  })

  # Save treatment, group, subject arguments?
  if(is.null(save.treatment)) save.treatment <- (n <= 10000) || n2
  if(is.null(save.group)) save.group <- (n <= 10000) || n2
  if(is.null(save.subject)) save.subject <- (n <= 10000) || n2

  # combine results into a single bootstrap object
  assign("%-%",  if(ratio) get("/") else get("-"))
  obj <- makeBootstrap(replicates = boot[[1]]$replicates %-% boot[[2]]$replicates,
                   observed = boot[[1]]$observed %-% boot[[2]]$observed,
                   n = n,
                   call = this.call,
                   seed.start = boot[[1]]$seed.start,
                   seed.end = boot[[2]]$seed.end,
                   dim.obs = boot[[1]]$dim.obs,
                   treatment = if(save.treatment) treatment,
                   group = if(save.group) group,
                   subject = if(save.subject) subject,
                   parent.frame = parent.frame,
		   label = boot[[1]]$label,
		   defaultLabel = resampMakeLabel("bootstrap",
		     this.call$data, this.call$statistic,
		     data2Expr = this.call$data2,
		     treatmentNames = u, ratio = ratio),
                   bootstrap.objects = boot,
		   ratio = if(ratio) T)
  # Add L, if supplied or found in the two objects
  if(L.is.numeric)
    obj$L <- L
  else if(!is.null(boot[[1]]$L)){
    L <- matrix(nrow = n, ncol = length(obj$observed))
    L[treat1,]  <-  boot[[1]]$L
    L[!treat1,] <- -boot[[2]]$L
    obj$L <- L
  }
  if(ratio && !is.null(obj$L)){
    # Need to modify L.  This includes the case where L is supplied;
    # see note (and warning) above.
    #  group 1:  Lratio =  Lindividual / observed2
    #  group 2:  Lratio = -Lindividual * ratio / observed2
    denominator <- boot[[2]]$observed
    L <- L / rep(denominator, each=nrow(L))
    L[!treat1,] <- L[!treat1,] * rep(obj$observed, each = sum(!treat1))
    obj$L <- L
  }
  # Add other components, if found in the two objects
  if(!is.null(boot[[1]]$Lstar)){
    if(ratio)
      obj$Lstar <- (1/denominator) *
	(boot[[1]]$Lstar - boot[[2]]$Lstar * rep(obj$observed, each=obj$B))
    else
      obj$Lstar <- boot[[1]]$Lstar - boot[[2]]$Lstar
  }
  if(!is.null(boot[[1]]$weights))
    obj$weights <- boot[[1]]$weights * boot[[2]]$weights

  # indices.  Check for both indices and compressedIndices;
  # if present in both objects, then combine.
  # If treatment is not sorted, then need some adjustment
  indices <- NULL
  if((!is.null(boot[[1]]$indices) ||
      !is.null(boot[[1]]$compressedIndices)) &&
     (!is.null(boot[[2]]$indices) ||
      !is.null(boot[[2]]$compressedIndices)))
    indices <- rbind(resampGetIndices(boot[[1]]),
		    resampGetIndices(boot[[2]]) + boot[[1]]$n)
  if(!is.null(indices)){
    if(bySubject){
      # reorder rows by sorted subjects, assuming rows currently
      # match sorted subjects within each treatment
      uniqueSubjects <- lapply(treatment.subject, function(x)sort(unique(x)))
      if(notSorted(unlist(uniqueSubjects)))
	indices <- indices[ order(unlist(uniqueSubjects)), ,drop=F]
    }
    else if(!n2 && notSorted(treat1, increasing=F)){
      ord <- order(!treat1)
      indices <- ord[indices]
      dim(indices) <- c(n, obj$B)
      indices[ord,] <- indices
    }
    if(!is.null(boot[[1]]$indices))
      obj$indices <- indices
    else if(!is.null(boot[[1]]$compressedIndices))
      obj$compressedIndices <- compressIndices(indices, n)
  }

  oldClass(obj) <- c("bootstrap2", "bootstrap", "resamp")
  on.exit()
  obj
}
# Added drop=F in some subscripting.
# Changed message when trace=T.
# Use resampMakeArgument for evaluating treatment, group, and subject
# Change hierarchical to bySubject, stratified to byGroup.
# New argument save.treatment
# By default save treatment, group, subject if data2 supplied of n <= 10K
# If L is supplied as a numerical object, expect it to be supplied
#   with the sign already reversed on the observations for the second
#   bootstrap call.  This is ugly, but prevents a discrepancy between
#   the L supplied and the L in the object.
# Do not copy all attributes of L.
# Do not set new.call$B (that was left over from permutationTestMeans)
# Add Lstar, indices, compressedIndices, weights to the output
# If data is a matrix with no column names, still consider it rectangular
# Do not assume that both objects have same type of indices (compressed or not)
# Use on.exit() in case something goes wrong.
# check for right length(group)
# Add ratio argument
# add twoSample.args - list of length 2, each a list of arguments
# Fix calculation of L in ratio case
# shorter calls in bootstrap.objects, if possible


print.bootstrap2 <-
function(x, digits = max(options()$digits - 3, 4), ...,
	 printCall = .Options.resample$printCall, printGroups = F)
{
  # print a bootstrap2 object, using the usual print.resamp, then
  # optionally print summary info for the two groups.
  NextMethod("print")
  if(printGroups && length(x$observed)){
    xx <- x$bootstrap.objects
    for(j in 1:2){
      cat("\nSummary Statistics for group ", names(xx)[j], ":\n")
      print(cbind(Observed = xx[[j]]$observed, xx[[j]]$estimate))
    }
  }
  invisible(x)
}


summary.bootstrap2 <-
function(object, probs = c(25, 50, 950, 975)/1000,
	 frame.eval = object$parent.frame, narrow=F, ...){
  # Summary of a bootstrap2 object.  
  # This calls summary.bootstrap for the parent bootstrap object, 
  # then adds short information for the two children.
  bo <- object$bootstrap.objects
  result <- list(combined = NextMethod("summary"),
		a = cbind(Observed = bo[[1]]$observed, bo[[1]]$estimate),
		b = cbind(Observed = bo[[2]]$observed, bo[[2]]$estimate))
  names(result)[2:3] <- names(bo)
  if(!is.null(object$ratio) && object$ratio)
    result$ratio <- T
  oldClass(result) <- "summary.bootstrap2"
  result
}

print.summary.bootstrap2 <-
function(x, digits = max(options()$digits - 3, 4), ...)
{
  # Print a summary.bootstrap2 object.
  # This prints the first component (a summary.bootstrap object),
  # then adds summary statistics for the two children.
  Names <- names(x)
  ratio <- (!is.null(x$ratio) && x$ratio)
  cat(if(ratio) "Results for the ratio for samples "
      else "Results for the difference for samples ",
      Names[2], " and ",
      Names[3], "\n", sep="")
  print(x[[1]], digits=digits, ...)

  cat("\nSummary Statistics for sample ", Names[2], ":\n", sep="")
  print(x[[2]], digits=digits, ...)
  cat("\nSummary Statistics for sample ", Names[3], ":\n", sep="")
  print(x[[3]], digits=digits, ...)
  invisible(x)
}

# boot2 <- bootstrap2(Verizon, mean(Time), treatment = Group)
# temp1 <- summary(boot2)
# temp1

