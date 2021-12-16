
# Simplified version of  ~timh/tilting/tilt.S


############################################################
# tiltWeights - find weights for given tau
############################################################

tiltWeights <- function(tau, L, tilt="exponential", weights = NULL,
                       group = NULL, lambda = NULL, details = F,
                       tol = 1e-6, maxiter = 50){
  # Compute weights using exponential or ML tilting
  if(is.matrix(L) && dim(L)[[2]] > 1)
    stop("L must be a vector or single-column matrix\n",
	 "See help(tiltDetails) for how to handle multivariate L.")

  Ltau <- outer(as.vector(subtractMeans(L, group = group, weights = weights)),
	       tau)
  N <- nrow(Ltau)
  if(length(weights) && length(weights) != N)
    stop("length of weights must match number of observations")

  # Get group information
  if(byGroup <- length(group)){
    if(byGroup != N)
      stop("Length of group must match number of observations in L")
    group.inds  <- split(1:N, group)
    nGroups <- length(group.inds)
    if(nGroups > 1){	# (If only one group, then like no group.)
      group.names <- names(group.inds)
      group.char <- as.character(group)
      groupSizes <- sapply(group.inds, length)
      # Divide L_gi by g_i = n_g/n, n_g = size of group g
      Ltau <- Ltau/(groupSizes[group.char] / N)
    }
  }
  else
    nGroups <- 1

  # Exponential tilting
  if(pmatch(casefold(tilt), "exponential", nomatch = 0) || is.null(tilt)){
    w <- exp(Ltau)
    if(length(weights))
      w <- weights * w
    return(w / groupSums(w, group = group, repeats = T))
  }
  
  # ML tilting (also recognize "mle" for "ml", but don't encourage this)
  else if(pmatch(casefold(tilt), "mle", nomatch = 0)){
    K <- length(tau)
    if(nGroups <= 1){			# ML tilting, no groups 
      w <- 1 / (1 - Ltau)
      if(max(Ltau) >= 1.0){
        warning("tau too large, weights set to NA")
	w[ Ltau >= 1.0] <- NA
      }
      if(length(weights))
        w <- weights * w
      return(w / rep(colSums(w), each=N))
    }
    # Have groups, use the group parameterization
    # normalize prior weights to sum to 1 (in each group)
    if(length(weights))
      weights <-  weights / groupSums(weights, group = group, repeats = T)
    else
      weights <- 1/as.vector(groupSizes[group.char])

    # If the normalization constants lambda are provided, use them.
    # Otherwise solve for them.
    if(have.lambda <- length(lambda)){
      if(!is.matrix(lambda))
	lambda <- as.matrix(lambda)
      if(length(lambda.names <- dimnames(lambda)[[1]])){
	if(!all(sort(lambda.names) == group.names))
	  stop("\n names(lambda) must match unique values of group")
      }
      else
	dimnames(lambda) <- list(group.names, NULL)
    }
    else{				# solve for lambda
      lambda <- matrix(NA, nGroups, K, dimnames=list(group.names,NULL))
      ord <- as.vector(unlist(group.inds))
      for(k in 1:K)
	lambda[,k] <- .C("S_tiltFindLambdas",
			 as.double(Ltau[ord, k]),
			 as.double(weights[ord]),
			 as.integer(nGroups),
			 as.integer(groupSizes),
			 as.double(tol),
			 as.integer(maxiter),
			 lambda = double(nGroups),
			 COPY = c(rep(F, 6), lambda = T))$lambda
    }
    w <- weights/(lambda[group.char,,drop=F] - Ltau)
    dimnames(w) <- list(NULL, names(tau))
    wsums <- groupSums(w, group = group)
    if(have.lambda && (any(abs(wsums - 1) > 1e-4)))
      warning("lambda gives weights for group which sum to ", wsums,
	      " not 1, before normalization.")
    # extra normalization just to be sure 
    w <- w / wsums[group.char,,drop=F]
    return(if(details) list(tau = tau, w=w, lambda = lambda) else w)
  }
  else
    stop("tilting method ", tilt, " is not recognized.")
}
# Do not support multivariate case (L %*% tau) (instead see help(tiltDetails)).
# Vectorize tau.
# In group case, do not add row names to output.
# Treat one group like no groups.
