# Define an influence() function, similar to jackknife(),
# and a print method.


##########################################################
# influence
##########################################################
influence <- 
function(data, statistic, args.stat = NULL,
	 group = NULL, subject = NULL,
	 label = NULL,
	 statisticNames = NULL,
	 assign.frame1 = F,
	 weights = NULL,
	 epsilon = min(.001, 1/n), unbiased = F, returnL = F,
	 save.group = NULL, save.subject = NULL,
	 subjectDivide = F, modifiedStatistic = NULL, ...)
{
  # Influence function, estimated using a finite value of epsilon.
  # epsilon = -1/(n-1)   -> usual (negative) jackknife (if no group)
  # epsilon =  1/(n+1)   -> positive jackknife         (if no group)
  # epsilon =  1/sqrt(n) -> butcher knife              (if no group)
  # weights may be a vector of unequal weights
  # If returnL=T, then only return L values, not everything else

  # For subject data, the model used here is that weights are
  # weights on the subject (not observations; if weights are supplied
  # they must be identical for all observations within a subject).
  # if(!subjectDivide), then
  #   * when calling function, subject weight is replicated to every
  #     observation in the subject
  #   * default weights (1/nSubjects on each subject) give total weight
  #     (added up across subjects) greater than 1, and larger total
  #     weight on larger subjects
  # if(subjectDivide), then
  #   * when calling function, subject weight is divided among the
  #     observations for that subject
  #   * default weights (1/nSubjects on each subject) effectively put
  #     smaller weight on observations in larger subjects
  #   * if weights are supplied as a vector of length nObservations,
  #     they are added up for each subject to give subject weights
  #
  # The ... argument is intended to handle extraneous arguments.  This is
  # needed because in bootstrap the call to influence (via resampGetL) is
  # constructed from the bootstrap call, which contains arguments influence
  # does not take.  There may be a better way of handling this.

  if(missing(statistic))
    stop("Argument statistic is missing with no default")

  # Capture call.
  func.call <- match.call()
  parent.frame <- sys.parent(1)

  # If statistic isn't a function or name of function or expression,
  # store it as an expression to pass to fit.func.
  substitute.stat <- substitute(statistic)
  if(!is.element(mode(substitute.stat), c("name", "function")))
    statistic <- substitute.stat

  # Get name of data.
  data.name <- substitute(data)
  if(!is.name(data.name))
    data.name <- "data"
  is.df.data <- is.data.frame(data)

  # How many observations, or subjects if sampling by subject?
  n <- nObservations <- numRows(data)  # n will change later if by subject

  # Save group or subject arguments?
  if(is.null(save.subject)) save.subject <- (n <= 10000)
  if(is.null(save.group))   save.group   <- (n <= 10000)

  # Assign data to frame 1 if needed
  if(is.null(assign.frame1)) assign.frame1 <- F
  if(assign.frame1)
    on.exit(if(exists(data.name, frame = 1))
	    remove(as.character(data.name), frame = 1))

  # Does statistic have a weights argument?  If expression, add weights.
  if(is.function(statistic)){
    if(!any(is.element(c("weights","..."), names(statistic))))
      stop("statistic does not have a weights argument")
  }
  else {
    if(missing(modifiedStatistic)){
      # create a modifiedStatistic, add a weights argument
      modifiedStatistic <- addWeightsInCall(statistic)
      if(identical(statistic, modifiedStatistic))
	stop("unable to find a function in your statistic with a `weights' argument")
    }
    else {
      # modifiedStatistic was supplied; either an expression like
      # mean(data, weights = Splus.resamp.weights)
      # in which case we use substitute(), or name of object containing
      # a stored expression.
      substitute.modstat <- substitute(modifiedStatistic)
      if(mode(substitute.modstat) != "name")
	modifiedStatistic <- substitute.modstat
      if(!is.element("Splus.resamp.weights",
		     all.names(modifiedStatistic)))
	stop('Your modifiedStatistic does not use a weights vector with name "Splus.resamp.weights"')
    }
    statistic <- modifiedStatistic
  }

  # group, subject, and weights arguments
  # Define byGroup, bySubject, nWeights (use below as length or logical)
  if(!missing(group) && is.df.data)
    group <- resampMakeArgument(func.call$group, data, parent.frame)
  byGroup <- length(group)
  if(!missing(subject) && is.df.data)
    subject <- resampMakeArgument(func.call$subject, data, parent.frame)
  bySubject <- length(subject)
  if(!missing(weights) && is.df.data)
    weights <- resampMakeArgument(func.call$weights, data, parent.frame)
  nWeights <- length(weights) # length may match number of subjects, or n

  # Set up sampling by `subject'; and do some weights calculations
  if(bySubject){
    # May redefine:  n = nSubjects, weights to have length nSubjects
    if(bySubject != nObservations)
      stop("length of subject must match number of observations in data")
    subjectIndices <- split(1:nObservations, subject)
    subjectSizes <- sapply(subjectIndices, length)
    subjectNames <- names(subjectIndices)
    subjectI <- match(subject, subjectNames)
    nSubjects <- length(subjectIndices)
    if(nWeights){
      if(nWeights == nSubjects) {
	# weights has length matching number of subjects;
	# weights should be ordered to match sort(subjects)
	if(length(namesWeights <- names(weights))){
	  if(notSorted(namesWeights))
	    weights <- weights[order(namesWeights)]
	  if(any(names(weights) != subjectNames))
	    stop("names of weights do not match names of subjects")
	}
      }
      else if(nWeights == nObservations){
	# weights should be identical for each observation for a subject
	if(notNested(outer = weights, inner = subject))
	  stop("Weights must be identical for all observations for each subject")
	if(subjectDivide)
	  weights <- groupSums(weights, subject)
	else # one weight for each subject, order to match sorted subjects
	  weights <- weights[ sapply(subjectIndices, function(x) x[1]) ]
      }
      else
	stop("length of weights must be number of observations or number of subjects")
    }
    else{  # weights not provided
      weights <- rep(1/nSubjects, nSubjects)
    }
    n <- nSubjects
  }
  else {
    # not by subject
    if(!nWeights)
      weights <- rep(1/n,n)
    else if(nWeights != n)
      stop("length of weights must match number of observations")
    subjectIndices <- NULL
  }

  # Create a function to manipulate weights, if necessary, to handle
  # dividing or replicating subject weight to observations.
  # Call with first argument only; set defaults for others.
  if(!bySubject)  # do nothing in this case
    fixWeights <- function(weight) weight
  else if(subjectDivide){
    # by subject, and divide subject weight among observations
    fixWeights <- function(weight, Order, div){(weight/div)[Order]}
    fixWeights$Order <- subjectI
    fixWeights$div <- subjectSizes
  }
  else {
    # by subject, and replicate subject weight to observations
    fixWeights <- function(weight, Order) {weight[Order]}
    fixWeights$Order <- subjectI
  }

  # Set up sampling by group (multiple-sample, or stratified)
  if(byGroup) {
    if(byGroup != nObservations)
      stop("length of group != number of observations in the data")
    if(bySubject){
      # get group for each subject (subject must be nested within group)
      subjectByGroup <- lapply(split(subject, group), unique)
      if(notNested(subjectByGroup))
	stop("some values of `subject' occur in multiple groups; subject must be nested within the `group' argument")
      group.inds <- lapply(subjectByGroup, match, table=subjectNames)
      # contains indices into elements of subjectIndices
      Sgroup <- sapply(subjectIndices, function(i, group) group[i[[1]]],
		      group = group)
      # shortened version of groups with length nSubjects
    }
    else {
      group.inds <- split(1:n, group)
      Sgroup <- group
    }
    nGroups <- length(group.inds)
    groupSizes <- sapply(group.inds, length)
  }
  else
    Sgroup <- NULL

  # Sum of weights (by group)
  # groupSums returns a single sum if group is NULL
  if(nWeights)
    sumW <- groupSums(weights, group)
  else
    sumW <- (if(byGroup) groupSizes/n else 1)

  # Set up: epsilon.
  if(epsilon == 0)
    stop("epsilon must be non-zero")
  if(epsilon > 1)
    stop("epsilon must be <= 1 (should usually be much smaller)")
  if(epsilon < 0){
    if(byGroup){
      # Don't have a nice way to check for epsilon too small in this case
    }
    else {
      if(nWeights){
	wmin <- min(weights)
	if(wmin <= 0 && epsilon < 0)
	  stop("Weights must be greater than zero when epsilon is negative")
      }
      else
	wmin <- 1/n
      minResultingW <- (1-epsilon) * wmin + epsilon * sumW
      if(minResultingW < 0){
	# epsilon too negative; change it
	epsilon <- max(epsilon, -wmin / (sumW-wmin))
	if(minResultingW < -1e-5)
	  warning("epsilon was adjusted to", epsilon,
		  "to prevent negative weights")

      }
    }
  }

  # Get function to evaluate the statistic given data and weights
  fit.func <- resampMakeFuncWeighted(statistic = statistic,
				     data.name = data.name,
				     is.df.data = is.df.data,
				     df.var.names = if(is.df.data) names(data),
				     is.null.args.stat = is.null(args.stat),
				     assign.frame1 = assign.frame1)
  # is: function(weights, data, statistic, args.stat)

  #
  # Ready to compute:  first, the observed values.
  observed <- fit.func(weights = fixWeights(weights),
		       data = data, statistic = statistic,
		       args.stat = args.stat)
  if(!is.atomic(observed))
    stop("Statistic must return numerical (non-list) data.")

  # Get parameter names and coerce matrix to vector.
  names.observed <- resampMakeNames(observed = observed,
				    statistic = substitute.stat,
				    statisticNames = statisticNames)
  dim.obs <- dim(observed)
  if(!is.null(dim.obs))
    observed <- as.vector(observed)
  names(observed) <- names.observed
  p <- length(observed)

  # Basic formula for weights for i'th replication is
  # weights * (1-epsilon) + sum(weights) * epsilon * (i'th observation)
  incrementW <- function(w, i, add){
    w[i] <- w[i] + add
    w
  }

  # Calculate replicates
  replicates <- matrix(NA, nrow=n, ncol=p,
		      dimnames = list(if(bySubject) subjectNames,
			              names.observed))

  if(!byGroup){  # no group.
    w <- (1 - epsilon) * weights
    epsilon2 <- epsilon * sumW
    for(i in 1:n)
      replicates[i,] <-
	fit.func(weights = fixWeights(incrementW(w, i, epsilon2)),
		 data = data, statistic = statistic, args.stat = args.stat)
  }
  else {
    for(g in 1:nGroups){
      w <- weights
      gi <- group.inds[[g]]
      w[gi] <- w[gi] * (1-epsilon)
      epsilon2 <- epsilon * sumW[g]
      for(i in gi)
	replicates[i,] <-
	  fit.func(weights = fixWeights(incrementW(w, i, epsilon2)),
		   data = data, statistic = statistic, args.stat = args.stat)
    }
  }

  if(anyMissing(replicates))
    warning("NA's encountered in replicates.")

  # The empirical influence function is defined as the limit of the next
  # expression, as epsilon -> 0.  For nonzero epsilon, subtract
  # means to get a better estimate of L, and use the means to estimate bias.
  L.eps <- (replicates - rep(observed, each=n))/epsilon
  Lmeans <- groupMeans(L.eps, group = Sgroup, weights = weights, repeats = 2)

  L <- L.eps - Lmeans$repeat.means
  attr(L, "epsilon") <- epsilon
  attr(L, "method")  <- "influence"
  if(returnL)
    return(L)

  n2 <- (if(byGroup) groupSizes else n) - unbiased
  if(!byGroup)
    est <- data.frame(Mean = colMeans(replicates),
		     Bias = t(Lmeans$means) / (n2 * epsilon),
		     SE = sqrt(colVars(L, weights=weights) / n2))
  else
    est <- data.frame(Mean = colMeans(replicates),
		     Bias = colSums(Lmeans$means / n2) / epsilon,
		     SE = sqrt(colSums(groupVars(L, weights=weights) / n2)))

  result <- list(call = func.call,
		 observed = observed,
		 replicates = replicates,
		 estimate = est,
		 B = n,
		 n = n,
		 dim.obs = dim.obs,
		 L = L,
		 epsilon = epsilon,
		 defaultLabel = resampMakeLabel("influence",
		   func.call$data, func.call$statistic))
  # add optional components, if they are non-null
  result <- c(result,
	     label = label,
	     groupSizes = if(byGroup) list(groupSizes),
	     group = if(save.group && !is.null(group)) list(group),
	     subject = if(save.subject && !is.null(subject)) list(subject))
  if(!is.function(statistic))
    result$modifiedStatistic <- statistic

  if(!nWeights){
    # Calculate acceleration and z0
    replicates2 <- matrix(NA, nrow = 2, ncol = p,
			 dimnames=list(NULL, names.observed))

    if(byGroup)
      L <- L * (n / groupSizes[as.character(Sgroup)])
    # see Davison & Hinkley p 217, for why we do that
    Lnorms <- sqrt(colSums(L^2))
    eps <- min(abs(epsilon) / Lnorms,
	      1/(1+n*max(abs(range(L))))) # avoid negative weights

    for(j in 1:p){
      # Use appropriate column of L for each dimension of statistic
      replicates2[1,j] <- fit.func(weights = fixWeights(weights + eps * L[,j]),
				  data = data, statistic = statistic,
				  args.stat = args.stat)[j]
      replicates2[2,j] <- fit.func(weights = fixWeights(weights - eps * L[,j]),
				  data = data, statistic = statistic,
				  args.stat = args.stat)[j]
    }
    z01 <- colSums(L^3)/(6*Lnorms^3)
    z02a <- ((replicates2[1,] + replicates2[2,] - 2*observed)/(Lnorms*eps)^2
	    ) / (2*n*Lnorms)		# curvature in least favorable direction
    z02b <- -(n-unbiased) * est$Bias / Lnorms # bias term
    z02 <- z02a + z02b

    result$estimate <- cbind(est, acceleration = z01, z0 = z01 + z02, cq = z02a)
    result$replicates2 <- replicates2
    result$epsilon2 <- eps
  }

  oldClass(result) <- c("influence", "resamp")
  result
}
# Use weighted version of colMeans later, instead of colSums; see "####" above
# Should subtract (weighted) mean of L before doing some calculations.
# Added subject/group arguments (RT, 8/10/01)
# Fix acceleration/z0 calculations in sampling by subject case
# Change default value of unbiased to F
# Stop if data frame uses variable name conflicting with resampMakeFuncWeighted
# Move above check inside resampMakeFuncWeighted, change its calling sequence
# Add returnL, option to only return the L values (rather than list with class)
# Check for weights=0 when epsilon<0.
# Avoid warning about epsilon when just numerical tolerance.
# Change output component from accel to acceleration
# Make assign.frame=NULL equivalent to F
# Remove commented-out error check for resampMakeFuncWeightsName
# Change default for epsilon
# Replace resampMakeFuncWeightsName with "Splus.resamp.weights"
# Remove argument substitute.stat from resampMakeFuncWeighted (not used)
# Change order of parts, to be more consistent with bootstrap, and so
#   fit.func is created after bySubject is known (may use this later).
# Change how subject is handled; instead of oWeights, use fixWeights,
#   add argument subjectDivide
# add arguments label and statisticNames, output label and defaultLabel

##########################################################
# influence
##########################################################
print.influence <- function(x, ...){
  if(!is.null(x$statistic))
    cat("Modified statistic:\n", deparse(x$statistic), "\n\n")
  NextMethod("print")
}
