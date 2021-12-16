# Functions in this file
#  tiltMean
#  tiltMeanSolve


##################################################
# tiltMean
##################################################
tiltMean <- function(tau, L, tilt = "exponential", weights=NULL, group=NULL,
                    lambda = NULL, tol = 1e-6, maxiter = 50,
		    debug=F)
{
  # Return vector of tilted means, one for each tau
  # In multiple group case, return sum of tilted means.
  n <- length(L)
  nw <- length(weights)
  if(nw){
    if(nw != n)
      stop("length of weights must be same as length of L")
    if(min(weights) < 0)
      stop("weights must be nonnegative")
  }

  if(byGroup <- length(group)){
    if(byGroup != n)
      stop("length of group must be same as length of L")
    group.inds <- split(1:n, group)
    nGroups <- length(group.inds)
    groupSizes <- unlist(lapply(group.inds, length))
    if(nGroups == 1){
      # Treat like no-group case
      byGroup <- F
    }
    else {
      group.names <- names(group.inds)
      if(!nw){
	weights <- rep(1,n)
	nw <- n
      }
      # normalize weights to sum to 1 in each group
      weights <- weights / groupSums(weights, group = group, repeats = T)
      if(notSorted(ord <- unlist(group.inds))){
	# Reorder so that groups are contiguous
	L <- L[ord]
	weights <- weights[ord]
      }
    }
  }
  else{  # No group
    nGroups <- 1
    groupSizes <- n
  }

  k <- length(tau)

  tilt <- casefold(tilt)
  if(pmatch(tilt, "exponential", nomatch = 0)){
    result <- .C("S_tiltMean_exp",
		 as.integer(n),
		 as.integer(k),
		 as.integer(nw),
		 as.double(L),
		 as.double(weights),
		 as.integer(nGroups),
		 as.integer(groupSizes),
		 as.integer(debug),
		 q = as.double(tau),
		 specialsok = T,  # for tau = +- infinity
		 COPY = c( rep(F,8), tau=T))$q
    return(list(tau=tau, q=result))
  }

  if(!pmatch(tilt, "mle", nomatch = 0)){  # officially "ml", but allow "mle"
    stop("\n Supplied tilting method ", tilt, " is not recognized.")
  }

  # rest of file is for ml case
  if(byGroup){
    have.lambda <- numRows(lambda)
    if(have.lambda){
	lambda <- as.matrix(lambda)
	if(nrow(lambda) != nGroups)
	  stop("lambda must have as many rows as there are unique group values")
	if(ncol(lambda) != k)
	  stop("lambda must have one colun for each value of tau")
	if(length(lambda.names <- dimnames(lambda)[[1]])){
	  # re-order lambda to be same as sorted group values
	  lambda <- lambda[order(lambda.names), , drop=F]
	  if(!all(dimnames(lambda)[[1]] == group.names))
	    stop("Row names of lambda must match unique values of group.")
	}
    }
    else
      lambda <- vector("numeric", k * nGroups)
  }
  else
    have.lambda <- F

  # Subtract mean, or group means
  if(byGroup){
    Lbars <- groupMeans(L, group = rep(1:nGroups, groupSizes),
		       weights = weights)
    L <- L - rep(Lbars, groupSizes)
  }
  else {
    Lbar <- mean(L, weights=weights)
    L <- L - Lbar
  }

  result <- .C("S_tiltMean_ml",
               as.integer(n),
               as.integer(k),
               as.integer(nw),
               as.double(L),
               as.double(weights),
               as.integer(nGroups),
               as.integer(groupSizes),
               as.integer(have.lambda),
               as.double(tol), # solving for lambda
               as.integer(maxiter), # solving for lambda
	       as.integer(debug),
               q = as.double(tau),
               lambda = as.double(lambda),
               specialsok = T,  # for tau = +- infinity
               COPY = c( rep(F,3), L = T, rep(F,5),
		   maxiter = T, F, q = T, lambda = T))
  if(byGroup){
    lambda <- matrix(result$lambda, nrow = nGroups, ncol = k,
		       dimnames = list(names(group.inds), NULL))
    result <- list(tau = tau, q = result$q + sum(Lbars), lambda = lambda)
  }
  else
    result <- list(tau = tau, q = result$q + Lbar)

  result
}
# If only one group, then treat like no-group case.
# lambda is a matrix (allow vector if length(tau)=1), return a matrix 
# Do not pass byGroup to C code; instead infer from number of groups.


##################################################
# tiltMeanSolve
##################################################
tiltMeanSolve <-
  function(q, L, tilt = "exponential", weights = NULL, group = NULL, initial,
	   tol = 1E-6, tol.tau = tol, maxiter = 50, debug=F)
{
  # Find tilting parameters tau, for which tilted mean equals q
  # If multiple group scae, solve q = sum of tilted means.
  tilt <- casefold(tilt)
  if(pmatch(tilt, "mle", nomatch = 0)){
    codeTilt <- 1
  }
  else if(pmatch(tilt, "exponential", nomatch = 0))
    codeTilt <- 0
  else
    stop("\n Supplied tilting method ", tilt, " is not recognized.")


  # find tau such that sum(w * L) = Q, where
  #    w = c weights exp(tau * (L-Lbar))	for exponential tilting
  #    w = c weights / (1 - tau * (L-Lbar))	for ML tilting, nGroups = 1
  #    w =   weights/ (lambda_i - tau * (L - Lbar))  for ML tilting, groups
  n <- length(L)
  nw <- length(weights)
  if(nw){
    if(nw != n)
      stop("length of weights must be same as length of L")
    if(min(weights) < 0)
      stop("weights must be nonnegative")
  }

  if(byGroup <- length(group)){
    if(byGroup != n)
      stop("length of group must be same as length of L")
    group.inds <- split(1:n, group)
    nGroups <- length(group.inds)
    groupSizes <- unlist(lapply(group.inds, length))
    if(nGroups == 1){
      # Treat like no-group case
      byGroup <- F
    }
    else {
      if(!nw){
	weights <- rep(1,n)
	nw <- n
      }
      # normalize weights to sum to 1 in each group
      weights <- weights / groupSums(weights, group = group, repeats = T)
      if(notSorted(ord <- unlist(group.inds))){
	# Reorder so that groups are contiguous
	L <- L[ord]
	weights <- weights[ord]
      }
    }
  }
  else{  # No group
    nGroups <- 1
    groupSizes <- n
  }

  k <- length(q)

  # Subtract the mean (or group means) for numerical stability,
  # and solve for adjusted q
  if(nGroups == 1){
    Lbars <- mean(L, weights=weights)
    L <- L - Lbars
  } else {
    Lbars <- groupMeans(L, rep(1:nGroups, groupSizes), weights=weights)
    L <- L - rep(Lbars, groupSizes)
  }
  qnew <- q - sum(Lbars)

  if(missing(initial) || length(initial) != k)
    initial <- qnew / sum(groupVars(L, rep(1:nGroups, groupSizes),
				   weights=weights))

  # The two "tau" components are used for both input and output
  resultList <- .C("S_tiltMean_solve",
                  as.integer(n),
                  as.integer(k),
                  as.integer(nw),
                  as.double(L),
                  as.double(weights),
                  as.integer(nGroups),
                  as.integer(groupSizes),
                  as.double(tol.tau),
                  as.double(tol),
                  as.integer(maxiter),
                  as.integer(codeTilt), # 1 => return both, 0 => exp only
		  as.integer(debug),
                  tau.exp = as.double(initial),
                  tau.ml = as.double(qnew),
                  lambda = double(k * nGroups),
                  # maxiter also needs COPY = T (was hard to track down)
                  COPY = c(rep(F, 9), maxiter = T, F,
		           tau.exp = T, tau.ml = T, lambda = T)
                  )
  if(codeTilt){
    if(byGroup){
        resultList$lambda <- matrix(resultList$lambda, nrow = nGroups, ncol = k,
                                    dimnames = list(names(group.inds), NULL))
      result <- resultList[c("tau.exp", "tau.ml", "lambda")]
    }
    else
      result <- resultList[c("tau.exp", "tau.ml")]
  }
  else
    result <- resultList["tau.exp"]

  c(list(q=q), result)
}
# If only one group, then treat like no-group case.
# returned lambda is a matrix
# Add q to the returned list.
# Do not pass byGroup to C code; instead infer from number of groups.
