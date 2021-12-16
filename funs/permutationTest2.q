#
# permutationTest2
#
permutationTest2 <- function(data, statistic, treatment, data2 = NULL,
			     B = 999, 
                             alternative = "two.sided",
			     ratio = F,
			     paired = F,
			     group = NULL,
			     combine = NULL,
                             combinationFunction = combinePValues.Fisher,
                             seed = .Random.seed,
			     trace = resampleOptions()$trace,
			     save.group = NULL,
			     save.treatment = NULL,
			     ..., L = "none")
{
  # Call this with either of:
  #   permutationTest2(data, statistic, treatment, ...)
  #   permutationTest2(data, statistic, data2, ...)
  alternative <- match.arg(alternative, c("two.sided", "greater", "less"))
  new.call <- this.call <- match.call()
  parent.frame <- sys.parent()

  # If data2 is used, construct a treatment vector to match, and set the
  # new data argument to concatenate data and data2. 
  if(n2 <- numRows(data2)){
    # Note -- n2 is used as a flag below, is nonzero if data2 supplied.
    n1 <- numRows(data)
    n <- n1 + n2
    # Create vector treatment; use names of data and data2, if given
    Name1 <- (if(is.name(this.call$data)) as.character(this.call$data)
	     else "1")
    Name2 <- (if(is.name(this.call$data2)) as.character(this.call$data2)
	     else "2")
    treatment <- rep(c(Name1, Name2), c(n1, n2))
    if((length(dim(data)))){
      data <- rbind(data, data2) 
      data.expr <- Quote(rbind(data, data2))
    }
    else{
      data <- c(data, data2)
      data.expr <- Quote(c(data, data2))
    }
    new.call$data <-
      substitute(data.expr, list(data = new.call$data, data2 = new.call$data2),
                 eval = T)
    new.call$data2 <- NULL
  }
  else
    n <- numRows(data)

  if(is.df <- is.data.frame(data)){
    # Look for `treatment' and `group' as variables in the data frame.
    # Also allow expressions like `region == "south"'.
    if(!n2)
      treatment <- resampMakeArgument(new.call$treatment, data, parent.frame)
    if(!missing(group))
      group <- resampMakeArgument(new.call$group, data, parent.frame)
  }

  # Argument `treatment'
  if(length(treatment) != n)
    stop("length(treatment) must equal number of observations in data")
  u <- unique(treatment)
  if(length(u) != 2)
    stop("permutation tests require that treatment have two unique values")
  
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

  # Check for unsupported arguments (related to sampling)
  if(any(w <- is.element(names(new.call),
			c("sampler", "sampler.prob", "sampler.args.group"))))
    stop(paste("arguments",paste(names(new.call)[w], collapse=" "),
	       "are not supported"))

  # Construct calls to `bootstrap' to do the sampling; one call for
  # each treatment value. Use `sampler = samp.permute'.  Although we are
  # doing separate calls for each treatment, the `full.partition' argument
  # ensures that the sampled indices work together to form essentially
  # full permutations.
  new.call[[1]] <- as.name("bootstrap")
  new.call$alternative <- NULL
  new.call$combine <- NULL
  new.call$combinationFunction <- NULL
  new.call$treatment <- NULL
  new.call$data2 <- NULL
  new.call$paired <- NULL
  new.call$ratio <- NULL
  new.call$resampleColumns <- NULL # in case this was passed by GUI
  if(is.null(new.call$L))
    new.call$L <- "none"
  if(is.null(new.call$save.indices))
    new.call$save.indices <- F
  if(paired)
    new.call$group <- group
  # Following supplies correct default values. In the case of the random
  # seed, the seed needs to be the same for both calls to bootstrap,
  # even if no seed is provided, in order to get the index sampling right.
  new.call$B <- B
  new.call$seed <- seed
  # This is also important to get index sampling correct.
  if(byGroup <- length(group)){
    if(byGroup != n)
      stop("group must have length equal to the total number of observations")
    new.call$group.order.matters <- F
    Table <- table(group, factor(treatment, levels=u))
  }
  part <- c("first", "last")
  boot <- vector("list", 2)
  names(boot) <- u

  for(i in 1:2){
    inds <- which(treatment == u[i])
    new.call$observed.indices <- {
      if(!n2){
        if(identical(inds, min(inds):max(inds)))
          call(":", min(inds), max(inds))
        else
          inds
      }
      else
        ifelse1(i == 1,
                call(":", 1, n1),
                call(":", n1+1, n))
    }
    new.call$sampler <- call("samp.permute", full.partition = part[i])
    # The size argument to sampler goes in the call to the sampler for
    # the no-group case, and in sampler.args.group for sampling by group.
    if(byGroup){
      # Need the number of times this treatment occurs in each group
      new.call$sampler.args.group <-
	as.list(lapply(Table[,i], function(n) list(size=n)))
    }
    else
      new.call$sampler$size <- length(inds)

    if(trace) cat("Resampling data for treatment group", i, "of 2.\n")
    if(trace && i==2) cat("(Samples are the complement of those for treatment group 1.)\n")

    boot[[i]] <- eval(new.call, sys.parent(1))
  }

  # Combine the two results into obj.
  assign("%-%",  if(ratio) get("/") else get("-"))
  obj <- list(call = this.call, 
              observed = boot[[1]]$observed %-% boot[[2]]$observed,
              replicates = boot[[1]]$replicates %-% boot[[2]]$replicates,
              B = B,
              n = n,
              dim.obs = boot[[1]]$dim.obs,
	      parent.frame = parent.frame,
	      label = boot[[1]]$label,
	      defaultLabel = resampMakeLabel("permutation",
		this.call$data, this.call$statistic,
		data2Expr = this.call$data2,
		treatmentNames = u, ratio = ratio),
              seed.start = boot[[1]]$seed.start,
              seed.end = boot[[1]]$seed.end,
              bootstrap.objects = boot,
	      ratio = if(ratio) T)
  if(is.null(boot[[1]]$label))
    obj$label <- NULL


  # Save treatment and group arguments?
  if(is.null(save.treatment)) save.treatment <- (n <= 10000) || n2
  if(is.null(save.group)) save.group <- (n <= 10000) || n2
  if(save.group)
    obj$group <- group
  if(save.treatment)
    obj$treatment <- treatment

  if(anyMissing(obj$replicates)){
    oldClass(obj) <- NULL
    on.exit(assign(".permutationTest2.partial.results", obj, immediate=T))
    stop("Missing values found in the replicates, partial results saved in .permutationTest2.partial.results")
  }

  p <- length(obj$observed)
  alternative <- rep(alternative, length=p)

  if(length(combine) && p == 1){
    warning("No combinations are possible with only one statistic")
    combine <- NULL
  }
  if(length(combine)){
    Names <- names(obj$observed)
    if(!is.list(combine))
      combine <- list(combine)
    if(any(w <- sapply(combine, is.character)))
      combine[w] <- lapply(combine[w], match, table=Names)
    temp <- findPvalues(obj$observed, obj$replicates,
		       B, p, alternative, combine, combinationFunction)
    Pvalues <- temp$Pvalues
    names(temp$combP) <- sapply(combine, function(x, Names) {
	if(!is.character(x)) x <- Names[x]
	paste(x, collapse=", ")
    }, Names=Names)
    obj$"combined-p-value" <- temp$combP
  }
  else
    Pvalues <- findPvalues(obj$observed, obj$replicates,
			  B, p, alternative, NULL, NULL)

  obj$estimate <- data.frame(Mean = colMeans(obj$replicates, na.rm = T),
                             SE   = colStdevs(obj$replicates, na.rm = T),
                             alternative = alternative,
                             "p-value" = Pvalues)
  oldClass(obj) <- c("permutationTest2", "permutationTest", "resamp")
  obj
}
# Changed message when trace=T.
# Add argument "paired"
# Use names of data and data2.
# Prevent failure when not all values of group are found within one treatment
# Change order of arguments for consistency with permutationTest2
# Add a label component (used when printing)
# Fix a problem when `data' & `data2' are data frames and `group'
#    is a variable in the data frames.
# Use resampMakeArgument for evaluating treatment, group.
# New arguments save.group, save.treatment; 
#   by default save if data2 supplied or n <= 10K
# Add parent.frame and defaultLabel to returned object.
# Add argument ratio.
# Set L="none" if not specified.
# set resampleColumns to NULL (in case argument passed by GUI)
# shorter calls in bootstrap.objects, if possible
