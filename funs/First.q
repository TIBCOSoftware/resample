# Keep this synchronized with .First.lib in gui/first.q
# The code here should match the end of .First.lib in gui/first.q
.First.lib <- function(lib.loc, section)
{
  # Initialization for the resample package
  #

  # Warn about conflicting versions.
  cat("The resample package contains newer versions of some functions.\n",
      "You may get warnings about masked functions.\n")

  # Make sure that resample is loaded at top of search list.
  searchList <- search()
  if(which(searchList == section)[1] > which(searchList == "splus")){
    attach(as.character(section), 2)
    searchList <- search()
    if(which(searchList == section)[1] > which(searchList == "splus")){
      # Tried to move, but failed
      cat("\n\nYou must load the resample package early in the search list\n",
	  "for it to work correctly.  Please unload using:\n",
	  '    detach("resample")\n',
	  "then reload using\n",
	  "    library(resample, first=T)\n")
    }
  }

  # set RESAMPTEST to location of loop test code (for internal use)
  setenv("RESAMPTEST", paste(lib.loc, section, "loop", sep="/"))

  invisible()
}
