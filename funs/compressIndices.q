# Copyright 2021. TIBCO Software Inc.
# This file is subject to the license terms contained
# in the license file that is distributed with this file.

# In this file:
#	compressIndices
#	uncompressIndices
#	print.compressIndices
#	"[.compressIndices"
#	cbind.compressIndices


compressIndices <- function(indices, n = nrow(indices)){
  # Convert an [m,B] matrix of indices into compressed form
  # n should be the maximum value possible in indices
  # (usually n = nrow(indices))
  originalDim <- dim(indices)
  m <- originalDim[1]
  B <- originalDim[2]
  k <- ceiling((n+m)/32)  # number of columns of the result
  result <- .C("S_compressIndices",
	      as.integer(originalDim),
	      as.integer(indices),
	      as.integer(n),
	      as.integer(k),
	      result = integer(k * B),
	      COPY=c(F,F,F,F,T))$result
  attributes(result) <- list(dim=c(k,B),
			    originalDim = originalDim,
			    n = n,
			    class = "compressIndices")
  result
}

uncompressIndices <- function(x){
  # Convert an object compressed by compressIndices into
  # an indices matrix
  ax <- attributes(x)
  Dim <- ax$originalDim
  if(length(Dim) != 2)
    stop("this is not a legal compressIndices object")
  result <- .C("S_uncompressIndices",
	      as.integer(Dim),
	      as.integer(x),
	      as.integer(ax$n),
	      as.integer(nrow(x)),
	      result = integer(prod(Dim)),
	      COPY=c(F,F,F,F,T))$result
  dim(result) <- Dim
  result
}


print.compressIndices <- function(x, ...){
  ax <- attributes(x)
  cat("Compressed indices, original dimensions ",
      ax$originalDim[1], " by ", 
      ax$originalDim[2], " with values up to ",
      ax$n,"\n", sep="")
  invisible(x)
}


"[.compressIndices" <- function(x, ..., drop=T){
  if(missing(..1) && (nargs() - !missing(drop)) == 3.) {
    # subscripting columns only, return another compressIndices object
    ax <- attributes(x)
    drop <- F
    x <- NextMethod("[")
    # restore attributes:  n, class, part of originalDim
    attributes(x) <- list(dim=dim(x),
			 originalDim = c(ax$originalDim[1], ncol(x)),
			 n = ax$n,
			 class = ax$class)
  }
  else x <- NextMethod("[")
  x
}


cbind.compressIndices <- function(..., deparse.level = 1){
  L <- list(...)
  if(!all(sapply(L, is, "compressIndices")))
    stop("Not all input values are compressIndices objects, cannot cbind.")
  originalDims <- sapply(L, function(x) attr(x, "originalDim"))
  originalN <-    sapply(L, function(x) attr(x, "n"))
  if(length(unique(originalDims[1,])) > 1 ||
     length(unique(originalN)) > 1)
    stop("Incompatible dimensions or n for compressed indices")
  x <- do.call("cbind", lapply(L, oldUnclass))
  if(ncol(x) != sum(originalDims[2,]))
    stop("Problem, actual dimensions do not match originalDims attribute.")
  attributes(x) <- list(dim=dim(x),
		       originalDim = c(originalDims[1,1], dim(x)[2]),
		       n = originalN[1],
		       class = "compressIndices")
  x
}

