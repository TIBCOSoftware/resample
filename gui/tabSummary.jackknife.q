tabSummary.jackknife <-
function(obj, probs = c(25, 50, 950, 975)/1000, print.short.p = T, 
	 percentiles = F, ABC = F, tLimits = F,
	 correlations = F)
{
  # This is used for both jackknife and influence function objects
  label <- (if(is(obj, "jackknife")) "Jackknife"
	   else if(is(obj, "influence")) "Influence"
	   else class(obj)
	   )

  cat("\n\t*** ", label, " Results ***\n")
  if(print.short.p)
    print(obj)

  if(!length(obj$observed)){
    if(!print.short.p && (percentiles || ABC || tLimits))
      warning("Statistic was length 0, cannot calculate summaries")
    return(invisible(NULL))
  }

  if(is.character(probs))
    probs <- eval(parse(text = paste("c(", probs, ")"))[[1]])    

  result <- list(object = obj)
  
  if(percentiles) {
    cat("\nPercentiles:\n")
    result$percentiles <- limits.percentile(obj, probs)
    print(result$percentiles)
  }

  if(ABC){
    x <- limits.abc(obj, probs=probs)
    result$ABC <- x
    cat("\nABC confidence limits:\n")
    print(x$abc.limits)
    if(!is.null(x$exp.limits)) {
      cat("\nABC/exponential confidence limits:\n")
      print(x$exp.limits)
    }
    if(!is.null(x$ml.limits)) {
      cat("ABC/ML confidence limits:\n")
      print(x$ml.limits)
    }
  }

  if(tLimits){
    cat( "\nT Confidence Intervals using ", label, " Standard Errors:\n" )
    result$tLimits <- limits.t(obj, probs=probs)
    print(result$tLimits)
  }

  if(correlations) {
    obj.cor <- cor(obj$replicates, na.method = "available")
    if(length(obj.cor) > 1) {
      cat("\nCorrelation of Replicates:\n")
      print(obj.cor)
      result$correlations <- obj.cor
    }
  }
  invisible(result)
}


"Done with tabSummary.jackknife.q" # for script file input
