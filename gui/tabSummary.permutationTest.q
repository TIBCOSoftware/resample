tabSummary.permutationTest <-
function(obj, 
	 probs = c(25, 50, 950, 975)/1000, 
	 percentiles = F, 
	 correlations = F, 
	 printSummary = T,
	 printCall = NULL,
	 ...)
{
  cat("\n\t*** Permutation Test Results ***\n")
  if(printSummary)
    print(obj , printCall = printCall, ...)

  if(!length(obj$observed)){
    if(!printSummary && percentiles)
      warning("Statistic was length 0, cannot calculate summaries")
    return(invisible(NULL))
  }

  if(is.character(probs))
    probs <- eval(parse(text = paste("c(", probs, ")"))[[1]])    

  result <- list()

  if( percentiles ) {
    result$percentiles <- limits.percentile(obj, probs)
    cat("\nPercentiles:\n")
    print( result$percentiles )
  }

  if(correlations && length(obj$observed) >= 2) {
    result$cor <- cor(obj$replicates, na.method = "available")    
    cat("\nCorrelation of Replicates:\n")
    print( result$cor )
  }

  invisible( result )
}
# Rename arguments, and reorder.
# Add option to not print call.  
# Move cor calculations later.
# Do not print correlations if length < 2.  Only print correlations if asked.


"Done with tabSummary.permutationTest.q" # for script file input
