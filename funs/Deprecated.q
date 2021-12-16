deprecateFun <- function(Old, New){
  # Old and New are character strings.
  # Create a function (to be assigned to the old name) that
  # gives a warning once per session, then calls the new function
  Old <- Old # force evaluation
  New <- New  
  substitute(function(...)
             {
               if(!exists(paste("warning.", OLD, sep=""), frame = 0)){
                 warning(OLD, " is deprecated, use:  ", NEW)
                 assign(paste("warning.", OLD, sep=""), "gaveWarning", frame = 0)
               }
               Call <- sys.call()
               Call[[1]] <- as.name(NEW)
               eval(Call, sys.parent())
             },
             list(OLD = Old, NEW = New))
}

limits.emp <- deprecateFun("limits.emp", "limits.percentile")
samp.boot.mc <- deprecateFun("samp.boot.mc", "samp.bootstrap")

rm(deprecateFun)
