# Copyright 2021. TIBCO Software Inc.
# This file is subject to the license terms contained
# in the license file that is distributed with this file.

hist.resamp <- function(x, ...) {
  invisible(plot.resamp(x, ...))
}

# That isn't quite right--returns invisibly even when plot=F
