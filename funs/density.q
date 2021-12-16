#
# density.q, function density
#
# Previously here, now in S+:
#	density (S3 generic)
#	density.default
#	bandwidth.sj.old


# Modified from S+6.0 version:
#   1) add weights argument
#   2) address Bug 22150 in two ways --
#      a) switch to double precision
#      b) increase gaussian cutoff window to 4s.d's. But allow old 3 s.d's
#         behavior by adding an extra window option, "3gaussian".  Note that
#         since there is a switch to double precision, you still can't
#         duplicate current S+6.0 results exactly.
#   3) address Bug 22126 by calling bandwidth.sj.old (which calls our own
#      versions of the corrupted C functions) if version == 6.0.


density.resamp <- function(x, ...){
  if(ncol(x$replicates) > 1)
    warning("There are multiple columns in replicates, density may not make sense")
  density(x$replicates, weights = x$weights, ...)
}


