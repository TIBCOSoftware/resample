# Copyright 2021. TIBCO Software Inc.
# This file is subject to the license terms contained
# in the license file that is distributed with this file.

# Miscellaneous functions


resampCor <- function(x, y, ..., resampleColumns = NULL){
  # Like cor(x), but instead of returning a matrix, return a vector
  #   containing the elements below the diagonal.
  # If y supplied, then like cor(x, y) except return vector w/ nice names.
  # If resampleColumns supplied then treat those columns like x, others y,
  if(missing(y)){
    M <- cor(x, ...)
    Names <- dimnames(M)[[1]]
    if(is.null(resampleColumns))
      lower <- col(M) < row(M)
    else {
      # make lower T if row is in resampleColumns and column is not
      lower <- M == (M+1)  # matrix of F's of right size & dimnames
      lower[, resampleColumns] <- T
      lower[resampleColumns, ] <- F
    }
    result <- M[lower]
    if(is.null(Names))
      Names <- paste("x", seq(along=M), sep="")
    names(result) <- paste("cor(",
			   Names[col(M)][lower], ",",
			   Names[row(M)][lower], ")", sep="")
    # order of names is switched, nicer in usual sub-diagonal case
  }
  else {
    M <- cor(x, y, ...)
    Names <- dimnames(M)
    if(!length(Names[[1]]))
      Names[[1]] <- paste("x", seq(nrow(M)), sep="")
    if(!length(Names[[2]]))
      Names[[2]] <- paste("y", seq(nrow(M)), sep="")
    result <- as.vector(M)
    names(result) <- paste("cor(",
			   Names[[1]], ",",
			   Names[[2]], ")", sep="")
  }
  result
}
# support argument y
# add argument resampleColumns


resampVar <- function(x, y, ..., unbiased = F, resampleColumns = NULL){
  # Like var(x), but instead of returning a matrix, return a vector
  #   containing the elements at and below the diagonal,
  #   and use "unbiased = F" (to give a functional statistic for resampling)
  # If y supplied, then like var(x, y, ...) except return vector w/ nice names.
  # If resampleColumns supplied then treat those columns like x, others y,
  if(missing(y)){
    M <- var(x, ...)
    Names <- dimnames(M)[[1]]
    if(is.null(resampleColumns))
      lower <- col(M) < row(M)
    else {
      # make lower T if row is in resampleColumns and column is not
      lower <- M == (M+1)  # matrix of F's of right size & dimnames
      lower[, resampleColumns] <- T
      lower[resampleColumns, ] <- F
    }
    result <- M[lower]
    if(is.null(Names))
      Names <- paste("x", seq(along=M), sep="")
    names(result) <- paste("cov(",
			   Names[col(M)][lower], ",",
			   Names[row(M)][lower], ")", sep="")
    # order of names is switched, nicer in usual sub-diagonal case
    names(result)[ col(M)[lower] == row(M)[lower] ] <-
      paste("var(", Names, ")", sep="")
  }
  else {
    M <- var(x, y, ...)
    Names <- dimnames(M)
    if(!length(Names[[1]]))
      Names[[1]] <- paste("x", seq(nrow(M)), sep="")
    if(!length(Names[[2]]))
      Names[[2]] <- paste("y", seq(nrow(M)), sep="")
    result <- as.vector(M)
    names(result) <- paste("cov(",
			   Names[[1]], ",",
			   Names[[2]], ")", sep="")
  }
  result
}
# support argument y
# add argument resampleColumns

print.resampHtest <- function(x, digits = 4, ...){
  # a resampHtest object is intended to be used in addition to a regular
  # hTest object, so not all information is here (e.g. alternative, data.name)
  if(!is.null(x$method))
    cat(paste("\n\t", x$method, "\n\n", sep = ""))
  if(!is.null(x$data.name))
    cat("data: ", x$data.name, "\n")
  for(i in setdiff(names(x), c("method", "data.name")))
    cat(i, ": ", sep="",
	paste(format(x[[i]], digits=digits, ...), collapse=" "), "\n")
}
