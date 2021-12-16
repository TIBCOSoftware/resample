# Keep this synchronized with .First.lib in funs/First.q
.First.lib <- function(lib.loc, section)
{
  # Initialization for the resample package
  #

  # Warn about conflicting versions.
  cat("Welcome to the S+Resample package.\n",
      "This package contains newer versions of some functions in S-PLUS.\n",
      "You may get warnings about masked functions.\n")

  if(is.sgui.app()){
    # Windows version, add graphical interface
    prop.file <- "guiResample.prp"
    info.file <- "guiResample.fni"
    path <- paste(lib.loc, section, sep = "/")
    prop.file <- paste(path, ".Prefs", prop.file, sep = "/")
    info.file <- paste(path, ".Prefs", info.file, sep = "/")

    # Load the gui objects into current session.
    guiLoadDefaultObjects("Property", FileName = prop.file)
    guiLoadDefaultObjects("FunctionInfo", FileName = info.file)

    # Add menu items, and better names and defaults
    if(interactive()) {
      addResampleMenus()

      # Change "Dependent" and "Independent" to "Response" and "Explanatory"
      guiModify("Property",
		Name = "SPropDependent",
		DialogPrompt = "&Response:")
      guiModify("Property",
		Name = "SPropIndependent",
		DialogPrompt = "&Explanatory:")
    }

    # Set option so that Help button commands not echoed to command line
    guiSetOption(option.name = "EchoExecuteStringInCmdWnd",
		 value.string = "F")
  }

  # Make sure that resample is loaded at top of search list.
  searchList <- search()
  if(which(searchList == section)[1] > which(searchList == "splus")){
    attach(as.character(section), 2)
    searchList <- search()
    if(which(searchList == section)[1] > which(searchList == "splus")){
      # Tried to move, but failed
      cat("\n\nYou must load the S+Resample package early in the search list\n",
	  "for it to work correctly.  Please unload using:\n",
	  '    detach("resample")\n',
	  "then reload using\n",
	  "    library(resample, first=T)\n")
      if(is.sgui.app() && interactive())
	cat("or check the 'Attach at top of search list' box in the\n",
	    "File:Load Library dialog.\n")
    }
  }

  # set RESAMPTEST to location of loop test code (for internal use)
  setenv("RESAMPTEST", paste(lib.loc, section, "loop", sep="/"))

  invisible()
}


.Last.lib <- function(...){
  if(is.sgui.app() && interactive()) {
    removeResampleMenus()
    # Remove our copy of three properties saved in 
  }
}
# I tried doing this, but then usual dialogs can't find these
#    guiRemove("Property", Name = "SPropDependent")
#    guiRemove("Property", Name = "SPropIndependent")




"Done with first.q" # for script file input
