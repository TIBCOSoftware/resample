# Functions
#
#   permutationTestMeans
#   print.permutationTestMeans
#
# Thurman's revision (12/5/2001) of Hesterberg's revision of Salmaso's code,
# 6/26/2000.
#
# Original code and additional investigations are in
# ~timh/bootstrap/permutation

permutationTestMeans <- function(data, treatment, data2 = NULL,
				B = 999,
				alternative = "two.sided",
				ratio = F,
				paired = F,
				group = NULL,
				combine = NULL,
				combineGroup = F,
				combinationFunction = combinePValues.Fisher,
				seed = .Random.seed,
				diffMeans = T,
				label = NULL,
				statisticNames = NULL)
{
  CALL <- match.call()
  parent.frame <- sys.parent()
  alternative <- match.arg(casefold(alternative), c("two.sided", "greater", "less"))

  # If data2 is used, construct a treatment vector to match, and set the
  # new data argument to concatenate data and data2. 
  if(n2 <- numRows(data2)){
    n1 <- numRows(data)
    n <- n1 + n2
    # Create vector treatment; use names of data and data2, if given
    Name1 <- (if(is.name(CALL$data)) as.character(CALL$data)
	     else "1")
    Name2 <- (if(is.name(CALL$data2)) as.character(CALL$data2)
	     else "2")
    treatment <- rep(c(Name1, Name2), c(n1, n2))
    if((length(dim(data)))){
      data <- rbind(data, data2) 
    }
    else{
      data <- c(data, data2)
    }
  }
  else
    n <- numRows(data)

  if(is.data.frame(data)){
    # Evaluate treatment and group as variables in the data frame,
    # then remove them before calculating means.
    # Also allow expressions like `region == "south"' (but no removal)
    # Then convert the data frame to an ordinary matrix
    #
    # If `combine' contains numerical values, the values should be
    # column numbers after treatment & group are removed.
    dfnames <- names(data)

    expr <- CALL$treatment
    if(any(is.element(all.names(expr), dfnames))){
      treatment <- resampMakeArgument(expr, data, parent.frame)
      if(is.name(expr) && is.element(expr, dfnames))
	data[[as.character(expr)]] <- NULL
    }
    expr <- CALL$group
    if(!missing(group)
       && any(is.element(all.names(expr), dfnames))){
      group <- resampMakeArgument(expr, data, parent.frame)
      if(is.name(expr) && is.element(expr, dfnames))
	data[[as.character(expr)]] <- NULL
    }
    data <- numerical.matrix(data, remove=T)
  }
  else if(!is.matrix(data))
    dim(data) <- c(n, length(data)/n)
  p <- ncol(data)

  # Get names for dimensions.  Remove dimnames, for speed.
  if(!is.null(statisticNames) && length(statisticNames) == p)
    Names <- statisticNames
  else {
    Names <- dimnames(data)[[2]]
    if(is.null(Names)){
      Names <- if(is.name(CALL$data)) as.character(CALL$data) else "Var"
      if(p>1)
	Names <- paste(Names, 1:p, sep="")
    }
  }
  dimnames(data) <- NULL


  # Argument `treatment'
  if(length(treatment) != n)
    stop("length(treatment) must equal number of observations in data")
  u <- unique(treatment)
  if(length(u) != 2)
    stop("permutationTestMeans compares two means, requires that treatment have two unique values")
  set1 <- (treatment == u[1])

  # Argument `paired'
  if(paired){
    if(length(group))
      warning("Argument group is ignored when paired=T")
    if(n2){  # data2 was supplied
      if(n1 != n2)
	stop("A paired permutation test requires n1 == n2")
      group <- rep(1:n1, 2)
    }
    else {  # argument treatment was supplied
      n1 <- sum(treatment == u[1])
      if(2*n1 != n)
	stop("A paired permutation test requires n1 == n2")
      # Create a group vector - unique value for each pair
      # Assume that there are individuals 1:n1, same order for each treatment
      group <- rep(0, n)
      group[treatment == u[1]] <- 1:n1
      group[treatment == u[2]] <- 1:n1
    }
  }


  stratified <- !is.null(group)
  if(stratified){
    # naming convention for objects:
    #   "s.*" = list, one component for each group
    #   "S.*" = array for which one dimension corresponds to groups
    s.inds <- split(1:n, group)
    s.set1 <- split(set1, group)
    S.n <- sapply(s.inds, length)
    s.data <- lapply(s.inds, function(data, i) data[i,,drop=F], data=data)
    K <- length(S.n)
    if(!paired){
      S.observed <- array(NA, c(K,p))
      S.replicates <- array(as.double(NA), c(B, K, p))
    }
  }
  else K <- 0  # may later create arrays with dimension K

  seed.start <- seed
  if(!missing(seed)) set.seed(seed)

  B1 <- B + 1
  sums <- array(0, dim = c(B1, p)) # first row for observed

  if(stratified){ # Calculate and save sums for each stratum
    ssums <- array(dim = c(B1, p)) # first row for observed
    for(k in 1:K){
      ssums[] <- .C("S_pTestSums",
                    as.double(s.data[[k]]),
                    as.integer(S.n[[k]]),
                    as.integer(p),
                    as.integer(s.set1[[k]]),
                    as.integer(B),
                    sums = double(B1*p),
                    NAOK=T, specialsok=T)$sums
      sums <- sums + ssums  # add all group sums, even for observed.
      if(!paired){
	S.observed[k,] <- ssums[1,]
	S.replicates[, k, ] <- ssums[-1,]
      }
    }
  }
  else
    sums[] <- .C("S_pTestSums",
                 as.double(data),
                 as.integer(n),
                 as.integer(p),
                 as.integer(set1),
                 as.integer(B),
                 sums = double(B1*p),
                 NAOK=T, specialsok=T)$sums

  observed <- sums[1,]
  replicates <- sums[-1,,drop = F]  # dim = c(B,p)

  # Convert to ratio, if bootstrapping the ratio
  if(ratio){
    # Converting group 1 sums to ratio of means
    # Convert observed and replicates, not groupObj$observed
    # xbar1-xbar2 = sum1/n1 - sum2/n2 = sum1(1/n1 + 1/n2) - sum/n2
    # xbar1/xbar2 = sum1/n1 / sum2/n2 = n2/n1 * sum1/(sum-sum1)
    n1 <- sum(set1)
    n2 <- n-n1
    sum0 <- colSums(data)
    observed <- n2/n1 * observed/(sum0-observed)
    replicates <- n2/n1 * replicates/(rep(sum0, each=B)-replicates)
  }


  #############################################################
  # Find individual p-values, and possibly combined p-values
  alternative <- rep(alternative, length=p)
  if(length(combine) && p == 1){
    warning("No combinations are possible with only one statistic")
    combine <- NULL
  }
  if(length(combine)){
    if(!is.list(combine))
      combine <- list(combine)
    if(any(w <- sapply(combine, is.character)))
      combine[w] <- lapply(combine[w], match, table=Names)
    temp <- findPvalues(observed, replicates,
		       B, p, alternative, combine, combinationFunction)
    Pvalues <- temp$Pvalues
    combP <- temp$combP
    names(combP) <- sapply(combine, function(x, Names) {
	if(!is.character(x)) x <- Names[x]
	paste(x, collapse=", ")
    }, Names=Names)
  }
  else
    Pvalues <- findPvalues(observed, replicates,
			  B, p, alternative, NULL, NULL)

  # Attach names
  names(observed) <- Names
  dimnames(replicates) <- list(NULL, Names)
  names(Pvalues) <- Names

  # Find P values for each stratum
  if(stratified && !paired && !ratio){
    # Get p-values for each stratum (and perhaps combine across variables)
    S.Pvalues <- matrix(NA, K, p, dimnames = list(names(S.n), Names))
    if(length(combine)){
      S.combP <- matrix(NA, K, length(combP),
		       dimnames = list(names(S.n), names(combP)))
      for(k in 1:K){
	temp <- findPvalues(S.observed[k,], as.matrix(S.replicates[,k,]),
			   B, p, alternative, combine, combinationFunction)

	S.Pvalues[k,] <- temp$Pvalues
	S.combP[k,] <- temp$combP
      }
    }
    else
      for(k in 1:K)
	S.Pvalues[k,] <-
	  findPvalues(S.observed[k,], as.matrix(S.replicates[,k,]),
		      B, p, alternative, NULL, NULL)

    # Also get p-values by combining p-values for each stratum.
    # This has low power for any dimensions with 2-sided tests
    if(combineGroup){
      T.combP <- rep(NA, p)
      names(T.combP) <- Names

      for(j in 1:p)
        T.combP[j] <-
          findPvalues(S.observed[,j], S.replicates[,,j],
                      B, K, rep(alternative[j], K),
                      combine=list(1:K), combinePValues.Fisher)$combP

    }
    groupObj <- list(observed = S.observed,
                    "p-value" = S.Pvalues,
                    "combined-p-value" = if(length(combine)) S.combP,
                    "combineGroup-p-values" = if(combineGroup) T.combP)
  }

  if(diffMeans && !ratio){
    # Converting group 1 sums to difference in means
    # Convert observed and replicates, not groupObj$observed
    # xbar1-xbar2 = sum1/n1 - sum2/n2 = sum1(1/n1 + 1/n2) - sum/n2
    n1 <- sum(set1)
    n2 <- n-n1
    beta0 <- -colSums(data) / n2
    beta1 <- (1/n1 + 1/n2)
    observed <- beta0 + beta1 * observed
    replicates <- (rep(beta0, each=B) +
		  rep(beta1, each=B) * replicates)
  }

  obj <- c(list(call = CALL,
	       observed = observed,
	       replicates = replicates,
	       estimate = data.frame(Mean = colMeans(replicates, na.rm=T),
		 SE = colStdevs(replicates, na.rm=T),
		 alternative = alternative,
		 "p-value" = Pvalues),
	       B = B,
	       n = n,
	       dim.obs = NULL,
	       "p-value" = Pvalues,
	       parent.frame = parent.frame,
	       defaultLabel = resampMakeLabel("permutation",
		 CALL$data, as.name("mean"),
		 data2Expr = CALL$data2,
		 treatmentNames = u, ratio = ratio),
	       seed.start = seed.start,
	       seed.end = .Random.seed),
	  if(!is.null(label))
	    list(label = label),
	  if(length(combine))
	    list("combined-p-value" = combP),
	  if(stratified && !paired)
	    list(group = groupObj))

  oldClass(obj) <- c("permutationTestMeans", "resamp")
  obj
}
# Argument treatment before B
# Use match.arg() for argument alternative
# Don't define set2
# If data is a data frame, look for treatment and group there.  (later changed)
# Add argument group.
# Add names to the combined-p-value
# Changed (again) how treatment and group are extracted from a df.
# Changed name of argument from group to strata.
# Replace computation of sums by C code
# Change strata back to group!
# Add argument diffMeans; if true (default) then convert observed
#  and replicates to difference in group means.
# Add arguments paired and data2
# In paired case, avoid many calculations on each stratum.
# Change order of arguments for consistency with permutationTestMeans
# Add a label component (used when printing)
# Add arguments label and statisticNames.
# Use resampMakeArgument for evaluating treatment and group
# Add parent.frame and defaultLabel to returned object.
# Add components parent.frame and defaultLabel
# Add argument ratio

# For help file, note:
# Elements of combine may be numeric, logical, or character.
#   For numeric and logical, number refers to data after
#   omitting treatment, group, and non-numeric variables.
# To do tests for only some variables in a matrix or data frame,
#   pass only those columns to this function.




print.permutationTestMeans <- function(x, ...){
  # Print a permutationTest object.
  # First use print.resamp (which prints p-values as part of the estimate
  # component).
  # Then if there is a combined p-value, print it too.
  # Finally, if there are p-values for individual groups, print those.
  NextMethod("print")
  if(length(x$"combined-p-value")){
    cat("\n")
    print(data.frame("Combined p-value:" = x$"combined-p-value",
		     check.names=F))
  }
  if(length(y <- x$group)){
    cat("\nP-values for each variable and stratum:\n")
    print(data.frame(alternative = x$estimate$alternative,
		     t(y$"p-value"), check.names=F))
    if(length(x$"combined-p-value")){
      cat("\nCombined p-values for each stratum:\n")
      print(t(y$"combined-p-value"))
    }
    if(length(y$"combineGroup-p-values")){
      cat(paste("\nNonparametric combination of individual stratum results",
                "\n(This is an alternative to using the overall sum as a test statistic.",
                if(any(x$estimate$alternative == "two.sided"))
                "\nFor two sided tests, this has good power for detecting interactions\nbetween the stratum variable and treatment effect,\nbut poor power for detecting differences between the two treatments;\nyou may want to switch to using one-sided alternatives.",
                "):\n", sep=""))
      print(y$"combineGroup-p-values")
    }
  }
  invisible(x)
}


