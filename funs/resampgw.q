#
# Defined in this file:
# functions: resampMakeFuncWeighted
#            addWeightsInCall
#            resampSubjectWeights
#


#
# resampMakeFuncWeighted
#
resampMakeFuncWeighted <-
function(statistic, data.name, is.df.data, df.var.names, 
	 is.null.args.stat, assign.frame1 = F)
{
  # Construct a function which takes arguments
  #  weights:   weights on arguments
  #  data:      vector, matrix, array, or data frame
  #  statistic: function or expression  (expression should use "weights")
  #  args.stat: optional arguments or objects (may be NULL)

  # Make sure that the name used internally for weights is not a variable 
  # in the data frame.
  if(is.df.data && 
     is.element("Splus.resamp.weights", df.var.names))
    stop("Variable Splus.resamp.weights is found in your data, ",
	 "but resampMakeFuncWeighted reserves that name for its own use.")

  ffs <- paste("function(weights, data, statistic, args.stat){\n")
  if(assign.frame1)
    ffs <- paste(ffs, 
                 "assign('", data.name, 
                 "', data, frame=1)\n", 
                 sep = "")
  # Handle function, or expression
  if(is.function(statistic))
    ffs <- 
      paste(ffs, 
	    if(is.null.args.stat) 
	      "statistic(data, weights=weights)}" 
	    else
	      "do.call('statistic',c(list(data, weights=weights), args.stat))}")
  else 
    ffs <- 
      paste(ffs, 
	    "eval(statistic, c(list(", data.name, " = data, ", 
	    "Splus.resamp.weights = weights)", 
	    if(is.df.data) ", data", 
	    if(!is.null.args.stat) ", args.stat",
	    "))}")

  eval(parse(text = ffs))
}
# add resampMakeFuncWeightsName (RT, 9/17/01)
# add argument df.var.names, to enable check against resampMakeFuncWeightsName 
# (RT, 9/18/01)
# replace resampMakeFuncWeightsName with "Splus.resamp.weights"
# removed argument substitute.stat (not used)


#
# addWeightsInCall
#
addWeightsInCall <- function(expr){
  # expr may include one or more function calls;
  # for every call to a function that has a "weights" argument,
  # modify the call to add the weights argument, in the form:
  #   weights = Splus.resamp.weights
  #
  # Note -- this does not add a weights argument to functions that accept
  # weights under a different name or as part of a ... list (this is
  # intentional).
  #
  # The list resampFunctionsWithWeights is a list of some functions
  # with weights.  You may add other functions to that list, that
  # should be viewed as having weights.
  moex <- mode(expr)
  switch(mode(expr),
	 "call" = {
	   if(is.element(expr[[1]], resampFunctionsWithWeights) ||
	      is.element("weights", names(get(expr[[1]])))){
	     # The function has an argument "weights"
	     if(is.element("weights", names(expr)))
	       stop("argument weights is already present")
	     else
	       expr$weights <- as.name("Splus.resamp.weights")
	   }
	   expr[-1] <- lapply(expr[-1], addWeightsInCall)
	   mode(expr) <- moex
	 },
	 "{" =, "if" =, "while" =, "repeat" =, "return" =,
	 "comment.expression" =, "expression" = {
 	   expr[] <- lapply(expr, addWeightsInCall)
	   mode(expr) <- moex
	 },
	 "(" =, "<-" =, "<<-" =, "for" = {
	   expr[-1] <- lapply(expr[-1], addWeightsInCall)
	   mode(expr) <- moex
	 })
  expr
}
# An example showing why you don't want to add weights to all functions with
# ... is the expression
#
#   c(mean(x), stdev(x))
#
# Note that stdev() accepts weights only through "...".  But c()
# also accepts ..., so if we add weights to all functions with weights or
# ..., we get the unintended result
#
#   c(mean(x, weights = Splus.resamp.weights),
#     stdev(x, weights = Splus.resamp.weights),
#     weights = Splus.resamp.weights)


# Here is a list of some functions with weights, that are likely to
# be called as a statistic for resampling functions.
# Exclude:  graphics functions, resamp utilities, sum & colSums, tabulate.
# Most common functions first.

resampFunctionsWithWeights <- c("mean",
			       "median",
			       "colMeans",
			       "colStdevs",
			       "colVars",
			       "lm",
			       "cor",
			       "var",
			       "glm",
			       "quantile",

			       "cdf",
			       "censorReg",
			       "coxph",
			       "discrim",
			       "gam",
			       "kaplanMeier",
			       "kurtosis",
			       "lmRobMM",
			       "lmsreg",
			       "loess",
			       "lsfit",
			       "ltsreg",
			       "skewness",
			       "stdev",
			       "survReg",
			       "survexp",
			       "survfit",
			       "tree",
			       "varcomp")



#
# resampSubjectWeights
#
resampSubjectWeights <-
  function(subject, weights = NULL, weights.out = "both", subjectDivide = F)
{
  #
  # Given subject vector and weights vector, of length = nSubjects
  # (number of unique values in subject) or nObservations (length(subject)),
  # return
  #
  #    weights vector of length = nSubjects (weights.out = "subject")
  #    weights vector of length = # of observations
  #      (weights.out = "observation")
  #    list of both, with names "subject" and "observation"
  #      (weights.out = "both")
  #
  # if(!subjectDivide), then
  #   * subject weight is replicated to every observation in the subject
  #   * default weights (1/nSubjects on each subject) give total weight
  #     (added up across subjects) greater than 1, and larger total
  #     weight on larger subjects
  #   * if weights are supplied as a vector of length 
  # if(subjectDivide), then
  #   * subject weight is divided among the observations for that subject
  #   * default weights (1/nSubjects on each subject) effectively put
  #     smaller weight on observations in larger subjects
  #   * if weights are supplied as a vector of length nObservations,
  #     they are added up for each subject to give subject weights
  # 
  # If original length is nSubjects, and named, check that names match
  # subject values, and reorder if necessary to match sorted values of subject.
  #
  # If original length is nObservations, check that weights are identical for
  # all observations in a subject.
  n <- length(subject)
  subjectIndices <- split(1:n, subject)
  subjectSizes <- sapply(subjectIndices, length)
  subjectNames <- names(subjectIndices)
  subjectI <- match(subject, subjectNames)
  nSubjects <- length(subjectIndices)
  if(nWeights <- length(weights)){
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
    else if(nWeights == n){
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

  if(pmatch(weights.out, "subject", nomatch = 0)) return(weights)

  obsWeights <-
    if(subjectDivide) (weights / subjectSizes)[subjectI]
    else weights[subjectI]

  if(pmatch(weights.out, "observation", nomatch = 0)) return(obsWeights)
  if(pmatch(weights.out, "both", nomatch = 0))
     return(list(subject = weights, observation = obsWeights))
  stop("uncrecognized value of weights.out")
}
