# This is modified from options()
resampleOptions <- 
function(..., TEMPORARY = FALSE)
{
  current <- .Options.resample
  if(nargs() == 0)
    return(current)
  temp <- list(...)
  if(length(temp) == 1 && is.null(names(temp))) {
    arg <- temp[[1]]
    switch(mode(arg),
      list = temp <- arg,
      character = return(.Options.resample[arg]),
      stop(paste("invalid argument:", arg)))
  }
  if(length(temp) == 0)
    return()
  n <- names(temp)
  if(is.null(n))
    stop("options must be given by name")
  changed <- current[n]
  names(changed) <- n
  current[n] <- temp
  assign(".Options.resample", current, frame = if(TEMPORARY) 1 else 0)
  invisible(changed)
}



.Options.resample <-
  list(printCall = 800,  # if call exceeds this many characters is not printed
       guiPrintCorrelations = T, # should tabSummary.bootstrap print cor?
       summaryCorrelations = T,  # should summary methods do correlations?
       # The following affect summary.bootstrap()
       summaryBCa = T,
       summaryTilt = NULL,  # decide based on whether L is present
       summaryT = F,
       trace = F
       )
# Add summary* components
# add trace=F
