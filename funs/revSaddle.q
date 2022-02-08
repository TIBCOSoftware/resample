# Copyright 2021. TIBCO Software Inc.
# This file is subject to the license terms contained
# in the license file that is distributed with this file.

# In this file:
#	revSaddlepointP
#	revSaddlepointPSolve


##################################################
revSaddlepointP <- function(tau, L, tilt = "exponential", weights = NULL,
			   group = NULL, Q = NULL, tau2 = NULL,
                           tol = 1e-6, tol.tau = 1e-6, debug=F)
{
  # Given "tilting parameters" {tau, tau2}, find P(Y > Q), where
  # Y is the mean of n = length(L) observations chosen from L_tau, the
  # tilted distribution with tilting parameter tau, and where Q is the
  # tilted mean of L_tau with tilting parameter tau2,
  # Q = tiltMean(L_tau, tau2).  If Q is provided, then find
  # tau2 = tiltMeanSolve(L_tau, Q).  The default Q is mean(L, weights).

  if(pmatch(casefold(tilt), "mle", nomatch=0))
    tilt <- "ml"
  else if(pmatch(casefold(tilt), "exponential", nomatch=0))
    tilt <- "exp"
  else
    stop("argument `tilt' should be either \"exponential\" or \"ml\".")

  n <- length(L)
  nw <- length(weights)
  if(nw && nw != n)
    stop("length of weights must be same as length of L")

  k <- length(tau)
  probs <- double(k)

  if(is.null(tau2)) {
    tau2 <- double(k)
    if(is.null(Q)) 		# Q then to be default, i.e. tau + tau2 = 0
      for(i in 1:k) {
	w <- tiltWeights(tau[i], L, tilt = tilt, weights = weights,
                         group = group, tol=tol)
	tau2[i] <- -tau[i]
	probs[i] <- 1 - saddlepointP(tau2[i], L, weights = w, group = group)
      }
    else
      for(i in 1:k) {
	w <- tiltWeights(tau[i], L, tilt = tilt, weights = weights,
                         group = group, tol=tol)
	tau2[i] <- tiltMeanSolve(Q, L, weights = w, group = group,
                                initial = -tau[i],
                                tol = tol, tol.tau = tol.tau,
				debug = debug)$tau.exp
	probs[i] <- 1 - saddlepointP(tau2[i], L, weights = w, group = group)
      }
  }
  else
    for(i in 1:k) {
      w <- tiltWeights(tau[i], L, tilt = tilt, weights = weights, group = group,
		       tol = tol)
      probs[i] <- 1 - saddlepointP(tau2[i], L, weights = w, group = group)
    }

  ans <- cbind(tau, tau2, probs)
  dimnames(ans) <- list(NULL, c("tau", "tau2", "probs"))

  ans
}
# Add argument debug.
# Allow groups in ml case.
# Pass tol to tiltWeights.


## May later add argument "trace = F", if "True", then print last 10
## (or user supplied) iteration results:
## {tau, tau2, f1 (diff to final mean), f2 (diff to final saddle, prob).



##################################################
revSaddlepointPSolve <-
  function(probs, L, tilt = "exponential", weights = NULL, group = NULL, 
	   Q = NULL, initial = NULL, initial.tau2 = NULL,
           useExpForInit = T,
	   tol = 1E-6, tol.tau = tol, maxiter = 100, debug=F)
{
  # Given {probs, Q}, find {tau, tau2} such that P(Y > Q) = probs, where
  # tau2 = tiltMeanSolve(L_tau, Q) for the tilted distribution L_tau with
  # tilting parameter tau, and Y is the mean of n = length(L) observations
  # chosen from L_tau.
  #
  # The default value for `Q' is the original weighted/sample mean:
  # Q = mean(L, weights=weights)

  if(pmatch(casefold(tilt), "mle", nomatch=0))
    tilt <- "ml"
  else if(pmatch(casefold(tilt), "exponential", nomatch=0))
    tilt <- "exp"
  else
    stop("Argument \"tilt\" should be either \"exponential\" or \"ml\".")

  n <- length(L)
  if(nw <- length(weights))
    if(nw != n)
      stop("length of weights must be same as length of L")

  good <- (!is.na(probs) & (probs > 0) & (probs < 1))
  if(!all(good))
    probs <- probs[good]
  k <- length(probs)

  if(lg <- length(group)){
    if(lg != n)
      stop("length of group must be same as length of L")
    if(tilt == "ml")
      stop("maximum likelihood tilting not supported with groups")
    group.inds <- split(1:n, group)
    nGroups <- length(group.inds)
    groupSizes <- unlist(lapply(group.inds, length))
    L <- L[unlist(group.inds)]
    weights <- weights[unlist(group.inds)]
    group <- group[unlist(group.inds)]
  }
  else{
    nGroups <- 1
    groupSizes <- n
  }

  ## Case of exponential tilting:
  if(tilt == "exp" || (tilt == "ml" && useExpForInit)) {

    ## define s as a mark for Q, where s = tau + tau2
    if(!length(Q))
      s <- 0 	## default, equivalent to: Q <- mean(L, weights=weights)
    else
      s <- tiltMeanSolve(Q, L, weights = weights, group = group,
                        tol = tol, tol.tau = tol.tau)$tau.exp

    if(nGroups > 1)
      L <- L * n / rep(groupSizes, groupSizes)

    if(is.null(initial))
      initial <- qnorm(probs) / stdev(as.vector(L), SumSquares=T)

    out <- .C("S_revSaddlepointPSolve_exp",
	     as.integer(k),
	     as.integer(nw),
	     as.double(L),
	     as.double(weights),
	     as.double(probs),
	     as.double(s),
             as.integer(nGroups),
             as.integer(groupSizes),
	     as.double(initial),
	     as.double(tol.tau),
	     as.double(tol),
	     as.integer(maxiter),
	     as.integer(debug),
	     numEval = integer(1),
	     tau = double(k),
	     COPY=c(rep(F,13), numEval=T, tau=T))

    tau <- out$tau
    tau2 <- s - tau
  }

  ## Case of "ml" tilting:
  if(tilt == "ml") {

    if(is.null(Q))
      Q <- mean(L, weights= weights)

    # set up initial
    if(useExpForInit) {
      initial  <- tau
      initial.tau2 <- tau2
    }
    if(is.null(initial))
      initial <- qnorm(probs) / sqrt(var(as.vector(L), SumSquares=T))
    if(is.null(initial.tau2))
      initial.tau2 <- -initial

    out <- .C("S_revSaddlepointPSolve_ml",
	     as.integer(n),
	     as.integer(k),
	     as.integer(nw),
	     as.double(L),
	     as.double(weights),
	     as.double(probs),
	     as.double(Q),
	     as.double(initial),
	     as.double(initial.tau2),
	     as.double(tol.tau),
	     as.double(tol),
	     as.integer(maxiter),
	     as.integer(debug),
	     numEval = integer(1),
	     tau = double(k),
	     tau2 = double(k),
	     COPY=c(rep(F,13), numEval = T, tau = T, tau2 = T))

    tau  <- out$tau
    tau2 <- out$tau2
  }

  ans <- cbind(probs, tau, tau2)
  if(!all(good)){
    ans2 <- ans
    ans <- matrix(NA, length(good), 3)
    ans[good,] <- ans2
  }
  dimnames(ans) <- list(NULL, c("probs", "tau", "tau2"))

  ans
}
# add argument debug
# support groups in exponential case (divide L instead of group tau's)
# return NA for probs <= 0 or >= 1
