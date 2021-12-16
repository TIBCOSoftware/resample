hist.resamp <- function(x, ...) {
  invisible(plot.resamp(x, ...))
}

# That isn't quite right--returns invisibly even when plot=F
