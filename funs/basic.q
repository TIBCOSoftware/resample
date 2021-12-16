#
# Basic bootstrap functions,
#
#   indexSums
#   indexMeans
#   indexProducts
#   indexVars
#
# Input variables `x' and `indices', both matrices.  The columns of 
# `x' are variables, and rows are observations.  Columns of `indices' 
# are bootstrap sample indices. 
#
# Given x[n,p], indices[m, b], return y[b,p], 
# where y[i,j] = fun( x[indices[, i],j]), and `fun' is `sum', `mean', 
# `prod', or `var'. 
#  
# By Tim H, modified to use matrices by RT, 12/27/00.
# Modified 2/16/01 (RT):  indices may have different number of rows than x.
#


indexSums <- function(x, indices){
  # Calculate answer[b, p] = sum( x[ indices[,b], p] )
  if(length(dim(x)) < 2){
    n <- length(x)
    p <- 1
  }
  else{
    n <- nrow(x)
    p <- ncol(x)
  }
  if(length(dim(indices)) < 2){
    m <- length(indices)
    b <- 1
  }
  else{
    m <- nrow(indices)
    b <- ncol(indices)
  }

  if(is.character(indices))
    indices <- match(indices, rowIds(x))
  else if(is.logical(indices)){
    if(m != n) stop("logical indices are wrong length")
    # convert to integer; convert F to zero
    # do not use apply(indices, 2, which) because number T varies
    indices <- indices * seq(length = m)
  }

  matrix(.C("S_bootstrapSums",
	    as.integer(n),
	    as.integer(p),
	    as.integer(m),
	    as.integer(b),
	    as.double(x),
	    as.integer(indices),
	    answer = double(p*b),
	    specialsok = T, COPY=c(rep(F,6), T))$answer, nrow=b)
}
# Allow specials, use copy=F
# support character and logical indices

indexMeans <- function(x, indices, group=NULL){
  # Calculate answer[b, p] = sum( x[ indices[,b], p] )
  # However, if there are groups, then calculate the sum of group means
  if(length(dim(x)) < 2){
    n <- length(x)
    p <- 1
  }
  else{
    n <- nrow(x)
    p <- ncol(x)
  }
  if(length(dim(indices)) < 2){
    m <- length(indices)
    b <- 1
  }
  else{
    m <- nrow(indices)
    b <- ncol(indices)
  }

  if(is.character(indices))
    indices <- match(indices, rowIds(x))
  else if(is.logical(indices)){
    if(m != n) stop("logical indices are wrong length")
    # convert to integer; convert F to zero
    # do not use apply(indices, 2, which) because number T varies
    indices <- indices * seq(length = m)
  }

  if(!length(group))
    return(matrix(.C("S_bootstrapMeans",
		     as.integer(n),
		     as.integer(p),
		     as.integer(m),
		     as.integer(b),
		     as.double(x),
		     as.integer(indices),
		     answer = double(p*b),
		     specialsok = T, COPY=c(rep(F,6), T))$answer, nrow=b))
  # Have groups.  Need to calculate sum of group means.
  # Divide each x value by the corresponding group size, then indexSums.
  if(length(group) != n)
    stop("length of group must match number of rows of x")
  gp <- as.numeric(factor(group))
  nGroups <- max(gp)
  groupSizes <- tabulate(gp[indices[1:m]], nGroups)
  # Used gp[indices[1:m]] rather than gp, in case did sampling with
  # sizes different than the original group sizes.
  # We assume here that all columns of indices have the same groupSizes.
  for(i in 1:nGroups)
    x[gp == i] <- x[gp == i] / groupSizes[i]
  # If x is a matrix that does the same as x[gp==i,] = ...
  matrix(.C("S_bootstrapSums",
	    as.integer(n),
	    as.integer(p),
	    as.integer(m),
	    as.integer(b),
	    as.double(x),
	    as.integer(indices),
	    answer = double(p*b),
	    specialsok = T, COPY=c(rep(F,6), T))$answer, nrow=b)
}
# Allow specials, use copy=F
# Add group argument
# support character and logical indices

indexProducts <- function(x, indices){
  # Calculate answer[b, p] = prod( x[ indices[,b], p] )
  if(length(dim(x)) < 2){
    n <- length(x)
    p <- 1
  }
  else{
    n <- nrow(x)
    p <- ncol(x)
  }
  if(length(dim(indices)) < 2){
    m <- length(indices)
    b <- 1
  }
  else{
    m <- nrow(indices)
    b <- ncol(indices)
  }

  if(is.character(indices))
    indices <- match(indices, rowIds(x))
  else if(is.logical(indices)){
    if(m != n) stop("logical indices are wrong length")
    # convert to integer; convert F to zero
    # do not use apply(indices, 2, which) because number T varies
    indices <- indices * seq(length = m)
  }

  matrix(.C("S_bootstrapProducts" ,
	    as.integer(n),
	    as.integer(p),
	    as.integer(m),
	    as.integer(b),
	    as.double(x),
	    as.integer(indices),
	    answer = double(p*b),
	    specialsok = T, COPY=c(rep(F,6), T))$answer, nrow=b)

}
# Allow specials, use copy=F
# support character and logical indices
# Remove log argument


indexVars <- function(x, indices){
  # Calculate answer[b, p] = var( x[ indices[,b], p] )
  # use 1/(m-1) weighting
  if(length(dim(x)) < 2){
    n <- length(x)
    p <- 1
  }
  else{
    n <- nrow(x)
    p <- ncol(x)
  }
  if(length(dim(indices)) < 2){
    m <- length(indices)
    b <- 1
  }
  else{
    m <- nrow(indices)
    b <- ncol(indices)
  }

  if(is.character(indices))
    indices <- match(indices, rowIds(x))
  else if(is.logical(indices)){
    if(m != n) stop("logical indices are wrong length")
    # convert to integer; convert F to zero
    # do not use apply(indices, 2, which) because number T varies
    indices <- indices * seq(length = m)
  }

  matrix(.C("S_bootstrapVars",
	    as.integer(n),
	    as.integer(p),
	    as.integer(m),
	    as.integer(b),
	    as.double(x),
	    as.integer(indices),
	    answer = double(p*b),
	    specialsok = T, COPY=c(rep(F,6), T))$answer, nrow=b)
}
# Allow specials, use copy=F
# support character and logical indices
