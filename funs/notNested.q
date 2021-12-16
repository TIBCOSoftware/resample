notNested <- function(outer, inner){
  # Check if inner is nested within outer; return T if not.
  # In normal use, outer and inner are two vectors.
  # However split(inner, outer) has already been computed, then
  # pass that result as `outer' (and inner is ignored)
  if(!is.list(outer))
    outer <- split(inner, outer)
  any(duplicated(unlist(lapply(outer, unique))))
}
