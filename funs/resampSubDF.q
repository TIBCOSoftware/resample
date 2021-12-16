resampSubDF <- function(x, i){
  # subscript rows of a data frame.
  # Intended for use in bootstrapping or other resampling.
  # i should be numeric or logical, not character
  ax <- attributes(x)
  f <- function(xj, i) if(length(dim(xj)) == 2) xj[i,,drop=F] else xj[i]
  x <- lapply(x, f, i=i)
  ax$row.names <- ax$row.names[i]
  attributes(x) <- ax
  x
}
