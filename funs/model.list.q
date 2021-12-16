#
# Functions in this file:
#   [.model.list
#   dim.model.list
#   resampChangeY.model.list

#
# [.model.list
#
#   subscripting method for model.list
#
# Basically, only two types supported:
#
# x[inds] or x[inds,] 
#
# In the first case, a list of the specified components are returned, 
# In the second case, all components are returned, with only the 
# specified rows of each component of x$data returned in $data.
# Column indices are not supported.
#
# Adapted from [.data.frame (see by doing get("[.data.frame"))
#
"[.model.list" <- function(x,...,drop=T){
  Nargs <- nargs() - !missing(drop)

# If at most one index argument, return list of appropriate components.

  if(Nargs < 3) { 
    x <- as.list(x)[...,drop=drop] # regular list subscripting
    return(x)
  }

  if(!missing(..2) && mode(..2) != "missing") { 
    stop("column indices not supported for class model.list")
  }

# Get rows from each component

  if(!missing(..1) && mode(..1) != "missing") { # row indices
    i <- ..1
    if(is.character(i)){
      rows <- names(x$data$Y)
      i <- pmatch(i, rows, duplicates.ok = T)
    }
    for(j in seq(along = x$data)) {
      if(!is.null(xj <- x$data[[j]])){
	x$data[[j]] <- (if(length(dim(xj)) == 2) xj[i,  , drop = F]
	else xj[i])
      }
    }
  }
  x
}
#)

#
# model.list method for dim (old-style generic? the def. of dim does not 
# have UseMethod, so not sure why this works, but it does)
#
dim.model.list <- function(x) c(numRows(x),length(x))


#
# update Y component of a model.list
#
# ml is a model.list
# newY is the new Y component
#

resampChangeY.model.list <- function(obj, newY){
  # change the $data$Y component of a model.list
  obj$data$Y <- newY
  obj
}
# old name was updatemlY
  
