limits.t <-
function(x, probs = c(25, 50, 950, 975)/1000,
	 df = "choose", adjust = T, z = F,
	 subset.statistic=1:p,
	 frame.eval = x$parent.frame, ...)
{
  # x should be a resamp object (e.g. bootstrap, or jackknife)
  # Compute confidence limits of the form
  #    estimate +- t * bootstrapStandardError
  #
  # This is not what is commonly known as a bootstrap t interval;
  # use bootstrapT for that.  Instead, this is a t interval using
  # a bootstrap standard error.
  #
  # df (degrees of freedom) may be one of:
  #   a number  (Inf gives a normal approximation)
  #   "smaller" (smaller of the sample or group sizes, minus 1);
  #             (this is n-1 in the simplest case)
  #   "normal"  (usual formula for non-pooled variance, assumes
  #             a normal sampling distribution; n-1 for simplest case;
  #   "pooled"  (n1+n2-2, for two-sample problems with no group)
  #   "choose"  normal if not sampling by group, smaller if by group.
  # In the group case, pooled assumes that the variance of the
  # contribution of each group to the overall statistic is proportional
  # to 1/(groupsize);
  # "normal" assumes it is proportional to c/(groupsize), where
  # c may differ between the two samples.
  #   
  # if(adjust=T), then adjust the degrees of freedom for number of
  #             bootstrap samples
  # if(z=T), then df and adjust are ignored; df is set to inf
  if(!is(x, "resamp"))
    stop("x must be a 'resamp' object.")
  if(is.null(frame.eval))
    frame.eval <- sys.parent(1)
  p <- length(x$observed)
  if(z){
    if(!missing(df) || !missing(adjust))
      warning("Arguments df and adjust are ignored when z=TRUE")
    df <- Inf
    adjust <- F
  }
  else if(!is.numeric(df)){
    # Calculate degrees of freedom.
    dfChoices <- c("choose", "smaller", "pooled", "normal")
    if(!is.element(df, dfChoices))
      stop("dfChoices must be one of:", paste(dfChoices, collapse=", "))
    if(!is(x, "bootstrap2")){
      # 1-sample problem (bootstrap, jackknife, influence)
      g <- resampGetArgument(x, argument="group", frame.eval=frame.eval)
      if(is.null(g))		# not sampling by group
	df <- x$n - 1
      else  			# sampling by group
	df <- switch(df,
		     choose =,
		     smaller = min(table(g)) - 1,
		     normal =,
		     pooled = x$n - length(unique(g)))
    }
    else { # bootstrap2 object
      # This is complicated if there is also sampling by group.
      nn <- c(x$bootstrap.objects[[1]]$n, x$bootstrap.objects[[2]]$n)
      gg <- lapply(x$bootstrap.objects, resampGetArgument,
		   argument="group", frame.eval=frame.eval)
      if(is.null(unlist(gg))){	# not sampling by group
	df <- switch(df,
		     smaller = min(nn)-1,
		     pooled = sum(nn)-2,
		     choose =, # same as normal
		     normal = {
		       v1 <- x$bootstrap.objects[[1]]$estimate$SE^2/nn[1]
		       v2 <- x$bootstrap.objects[[2]]$estimate$SE^2/nn[2]
		       (v1 + v2)^2/(v1^2/(nn[1]-1) + v2^2/(nn[2]-1))
		     })
      }
      else {			# sampling by group, in addtion to two samples
	counts <- lapply(gg, table)
	df <- switch(df,
		     choose =, # same as smaller
		     smaller = min(unlist(counts)) - 1,
		     pooled = sum(nn) - length(unlist(counts)),
		     normal = {
		       v1 <- x$bootstrap.objects[[1]]$estimate$SE^2/nn[1]
		       v2 <- x$bootstrap.objects[[2]]$estimate$SE^2/nn[2]
		       (v1 + v2)^2/(v1^2/(nn[1]-length(counts[[1]])) +
				    v2^2/(nn[2]-length(counts[[2]])))
		     })
      }
    }
  }
  # adjust the degrees of freedom to reflect limited bootstrap sampling
  if(adjust)
    df <- 1/(1/df + 1/(x$B-1))	# do not simplify; df may be Inf
  result <- x$observed + outer(x$estimate$SE, qt(probs, df))
  dimnames(result) <- list(names(x$observed),
			   paste(100 * probs, "%", sep = ""))
  result[subset.statistic,,drop=F]
}
# Add arguments "df" and "adjust".
# Handle df for bootstrap2 objects correctly.
# Add option "normal".  Revamp handling of df.
