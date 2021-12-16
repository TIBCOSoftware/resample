callTo <- function(f, fname){
  # Create a function with the same arguments as f, which just
  # evaluates the functionBody of f.
  # Supply either f or fname.
  if(missing(f))
    f <- getFunction(fname)
  else {
    fname <- substitute(f)
    if(!is.name(fname))
      stop("f must be the name of a function, not an expression which creates a function")
  }
  if(!is.function(f))
    stop("f (or get(fname)) is not a function")
  # two cases, hard case is where fname is something like "[.foo"
  if(make.names(fname) == fname)
    functionBody(f) <- substitute(eval(functionBody(replaceThis)),
				 list(replaceThis=as.name(fname)))
  else
    functionBody(f) <- substitute(eval(functionBody(get(replaceThis))),
				 list(replaceThis=fname))
  f
}
# Bug 17781 -- setting the functionBody removes any self-doc
# Note that if my fix there is accepted then I might need to modify this
# function to replace the proper part of the functionBody.

# Usage for `callTo' when defining a generic function and methods:
#       f.default <- function(args){body}
#       f <- f.default
#	setGeneric("f")
#	setDefaultMethod("f", callTo(f.default))
#	setMethod("f", "matrix", callTo(f.matrix))

# The method is of the form:
#	function(same arguments as f.matrix when the method was set)
#	  eval(functionBody(f.matrix))
