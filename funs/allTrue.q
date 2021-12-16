# A couple of utility routines useful in loop tests:
#
#   allTrue (now removed, is in S+8)
#   all.equal.excluding
#   all.but.call

all.equal.excluding <- function(x, y, ..., excluding=NULL, attrs=NULL){
  # Like all.equal, but exclude components in `excluding',
  #   and excluding attributes named in `attrs'.
  #
  # `excluding' and `attrs' should be character, names of components 
  #   and attributes.
  # 
  # For example:
  #   all.equal.excluding(obj1, obj2, excluding = c("call", "residuals"))
  for(i in intersect(names(x), excluding)) x[[i]] <- NULL
  for(i in intersect(names(y), excluding)) y[[i]] <- NULL
  for(i in intersect(names(attributes(x)), attrs)) attr(x,i) <- NULL
  for(i in intersect(names(attributes(y)), attrs)) attr(y,i) <- NULL
  all.equal(x,y, ...)
}
# Changed algorithm to catch components listed in excluding that
#   were present but NULL; previously had loop if(!is.null(x[[i]])) ...
#   Changed algorithm for attributes similarly.  However, there
#   may be attributes that show up using attr but not attributes; if
#   so then add this code back in (after current loops):
#    for(i in attrs){
#      if(!is.null(attr(x,i))) attr(x,i) <- NULL
#      if(!is.null(attr(y,i))) attr(y,i) <- NULL
#    }



all.but.call <- function(x, y, ...){
  # Compare objects, except for call, actual.calls, and defaultLabel components
  # This is for the resample library, not intended for other use
  x$call <- NULL
  y$call <- NULL
  x$actual.calls <- NULL
  y$actual.calls <- NULL
  x$defaultLabel <- NULL
  y$defaultLabel <- NULL
  all.equal.excluding(x, y, ...)
}
# Did calculations directly, rather than calling all.equal.excluding
# Call all.equal.excluding, in case there are additional components to ignore.
