# In this file:
#   tabSummary.bootstrap  (menu function in contextBootstrap.ssc)
# 

# Arguments through correlations are consistent with old tabSummary.bootstrap,
# albeit with different argument names.


tabSummary.bootstrap <-
function(obj,
	 probs = c(25, 50, 950, 975)/1000,
	 printSummary = T,
	 percentiles = T,
	 bca = T,
	 correlations = F,
	 tilting = F,
	 tLimits = F,
	 save.obj = "",
	 viewReplicates = F,
	 replicatesName = "",  # save using this name
	 ...,
	 printCall = NULL)
{
  # Do a summary of a bootstrap object.
  # This does some printing, and returns the result invisibly.
  # printSummary = print the original object
  # percentiles, bca, tilting, tLimits = for confidence limits
  # correlations = correlations between columns of replicates
  # save.obj = save results under this name
  # viewReplicates = for use when right-click on bootstrap object in Obj Expl
  # replicateName = name to save replicates, so can view them
  # printCall = passed if do printSummary
  # ... = passed if do printSummary

  cat("\n\t*** Bootstrap Results ***\n")
  if(printSummary)
    print(obj , printCall = printCall, ...)

  if(!length(obj$observed)){
    if(!printSummary && (percentiles || bca || tilting || tLimits))
      warning("Statistic was length 0, cannot calculate summaries")
    return(invisible(NULL))
  }

  if(is.character(probs))
    probs <- eval(parse(text = paste("c(", probs, ")"))[[1]])

  result <- list()

  if(percentiles) {
    result$percentiles <- limits.percentile(obj, probs)
    cat("\nPercentiles:\n")
    print( result$percentiles )
  }

  if(bca) {
    result$bca <- limits.bca( obj, probs )
    cat("\nBCa Confidence Intervals:\n")
    print( result$bca )
  }

  if(tilting){
    result$tilt <- limits.tilt( obj, probs )
    cat( "\nTilting Confidence Intervals:\n" )
    print( result$tilt )
  }

  if(tLimits){
    result$tLimits <- limits.t( obj, probs )
    cat( "\nT Confidence Intervals using Bootstrap Standard Errors:\n" )
    print( result$tLimits )
  }

  if(correlations && length(obj$observed) >= 2) {
    result$cor <- cor(obj$replicates, na.method = "available")
    cat("\nCorrelation of Replicates:\n")
    print( result$cor )
  }

  if(replicatesName != ""){
    assign(replicatesName, obj$replicates, where=1)
  }

  if(viewReplicates){
    if(replicatesName == ""){
      warning("Must save replicates to view them;",
	      "I'll save them as .resample.replicates.")
      replicatesName <- ".resample.replicates"
    }
    guiOpenView("matrix", Name = replicatesName)
  }

  invisible( result )
}
# Rename arguments; keep in same order as original (may regret that).
# Add option to not print call.
# Move cor calculations later.
# Do not print correlations if length < 2.  Only print correlations if asked.
# Add arguments for tilting and t with boot SE.
# Add options to save and view replicates.


"Done with tabSummary.bootstrap.q" # for script file input
