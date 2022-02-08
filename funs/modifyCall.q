# Copyright 2021. TIBCO Software Inc.
# This file is subject to the license terms contained
# in the license file that is distributed with this file.

##########################################################
# modifyCall
##########################################################
modifyCall <- function(Call, ..., NEWFUN){
  # Create a new function call by modifying an existing call
  #
  # Call is a function call to be modified, say mean(x=Y, trim=.2)
  # It should be created by match.call(), so arguments are given by
  # name, not position; e.g. mean(Y, .2) will not work.
  #
  # ... include
  # (1) Character strings, interpreted as argument names to be kept; e.g.
  #        modifyCall(Call, "x")         # mean(x=Y)
  #     Argument names are the intersection of those given and those
  #     present in Call.  If you want to have arguments in Call kept
  #     by default then use update().
  # (2) Arguments given in name=value form, e.g.
  #       modifyCall(Call, "x", na.rm=T) # mean(x=Y, na.rm=T)
  #     These are added to those created in step (1).
  #     These may be given quoted, e.g. if Z=1:99, then compare
  #        modifyCall(Call, x=Quote(Z))  # mean(x=Z)
  #        modifyCall(Call, x=Z)         # long, contains actual data
  #
  # NEWFUN, if present, is a new function name, e.g.
  #        modifyCall(Call, "x", NEWFUN="median") # median(x=Y)
  l <- list(...)
  Names <- names(l)
  if(is.null(Names)) Names <- rep("", length(l))
  set1 <- (nchar(Names) == 0)
  newCall <- Call[ c(1, match(unlist(l[set1]), names(Call), nomatch=0))]
  for(i in Names[!set1])
    newCall[[i]] <- l[[i]]
  if(!missing(NEWFUN))
    newCall[[1]] <- {if(is.character(NEWFUN)) as.name(NEWFUN)
                    else if(is.name(substitute(NEWFUN))) substitute(NEWFUN)
		    else NEWFUN
		  }
  newCall
}
