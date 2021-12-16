#
#
# Subscripting method for model.frames: same as data.frames 
# (for new-style generics, it doesn't inherit the method, even though 
# model.frame inherits from data.frame as a class).
#
# setMethod("[","model.frame",getMethod("[","data.frame"))
"[.model.frame" <- function(x, ..., drop = T)
  "[.data.frame"(x, ..., drop=drop)

