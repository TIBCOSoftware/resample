# This file contains:
#  moveOldResampleMenus
#  restoreOldResampleMenus
#  addResampleMenus
#  removeResampleMenus
#  resampAddArgumentsToCall

moveOldResampleMenus <- function(){
  # This function hides the original jackknife and bootstrap menus
  resample.name <- "SPlusMenuBar$Statistics$Resample"
  guiModify( "MenuItem",
            Name = paste(resample.name,"$Bootstrap", sep=""),
            Hide = T)
  guiModify( "MenuItem",
            Name = paste(resample.name,"$Jackknife", sep=""),
            Hide = T)
}
# Hide, rather than change the MenuItemText

restoreOldResampleMenus <- function(){
  # This function unhides original jackknife and bootstrap menus
  resample.name <- "SPlusMenuBar$Statistics$Resample"
  guiModify( "MenuItem",
            Name = paste(resample.name,"$Bootstrap", sep=""),
	    Hide = F)
  guiModify( "MenuItem",
            Name = paste(resample.name,"$Jackknife", sep=""),
	    Hide = F)
}


addResampleMenus <- function(){
  # This function is called when running under a GUI to add resample menus.

  # Move old menus (call them e.g. "Original bootstrap..."
  moveOldResampleMenus()

  # Under Statistics/Resample:
  createMenuBootstrap()
  createMenuJackknife()
  createMenuPerm()
  g <- function(Name, MenuItemText, Command, Index,
	       location = "SPlusMenuBar$Statistics$Resample"){
    guiCreate( "MenuItem",
	      Name = paste(location, Name, sep="$"),
	      Type = "MenuItem",
	      Action = "Function",
	      MenuItemText = MenuItemText,
	      Command = Command,
	      Overwrite = F,
	      Index = Index)
  }
  g("SummaryStatResample", "Summary Statistics", "menuDescribeResample", 6)
  g("CorResample", "Correlations", "menuCorResample",7)
  g("t-testResample", "One-sample t", "menuTTest1resample", 8)
  g("T-testResample", "Two-sample t", "menuTTest2resample", 9)
  g("PropResample", "Proportions", "menuPropResample", 10)
  g("LinearResample", "Linear Regression", "menuLmResample", 11)


  # Under other menus:
  createMenuCorResample()
  createMenuDescribeResample()
  createMenuLinearResample()
  createMenuPropResample()
  createMenuTTest1resample()
  createMenuTTest2resample()

  # Add manual to menus
  guiCreate( "MenuItem",
	    Name = "SPlusMenuBar$Help$On_Line_Manuals$ResampleSeparator",
	    Type = "Separator")

  libraryName <- {
    if(is.element("resample", search()))
      "resample"
    else
      find("bootstrap.default")[1]
  }
  guiCreate( "MenuItem",
	    Name = "SPlusMenuBar$Help$On_Line_Manuals$Resample",
	    Type = "MenuItem",
	    DocumentType = "Any Documents",
	    Action = "Open",
	    MenuItemText = "&Resample Library User's Manual",
	    Command = paste(database.path(libraryName),
	      "/../resample.pdf", sep=""))
  guiCreate( "MenuItem",
	    Name = "SPlusMenuBar$Help$On_Line_Manuals$ResampleTutorial",
	    Type = "MenuItem",
	    DocumentType = "Any Documents",
	    Action = "Open",
	    MenuItemText = "&Resample Library Tutorial",
	    Command = paste(database.path(libraryName),
	      "/../doc/tutorial.ssc", sep=""))
  guiCreate( "MenuItem",
	    Name = "SPlusMenuBar$Help$On_Line_Manuals$ResampleExamples",
	    Type = "MenuItem",
	    DocumentType = "Any Documents",
	    Action = "Open",
	    MenuItemText = "&Resample Manual Examples",
	    Command = paste(database.path(libraryName),
	      "/../doc/resampleExamples.ssc", sep=""))
}



removeResampleMenus <- function(){
  # This function is called when the resample library exits
  # in order to remove menu items.

  # Under Statistics/Resample:
  removeMenuBootstrap()
  removeMenuJackknife()
  removeMenuPerm()
  g <- function(Name, ...,
	       location = "SPlusMenuBar$Statistics$Resample"){
    guiRemove( "MenuItem", Name = paste(location, Name, sep="$"))
  }
  g("SummaryStatResample", "Summary Statistics", "menuDescribeResample", 6)
  g("CorResample", "Correlations", "menuCorResample",7)
  g("t-testResample", "One-sample t Test", "menuTTest1resample", 8)
  g("T-testResample", "Two-sample t Test", "menuTTest2resample", 9)
  g("PropResample", "Proportions", "menuPropResample", 10)
  g("LinearResample", "Linear Regression", "menuLmResample", 11)

  # Under other menus:
  removeMenuCorResample()
  removeMenuDescribeResample()
  removeMenuLinearResample()
  removeMenuPropResample()
  removeMenuTTest1resample()
  removeMenuTTest2resample()

  # Restore old menus:
  restoreOldResampleMenus()

  menu.item.name <- "SPlusMenuBar$Help$On_Line_Manuals$ResampleSeparator"
  if(is.element(menu.item.name, guiGetObjectNames("MenuItem")))
    guiRemove( "MenuItem", Name = menu.item.name)
  menu.item.name <- "SPlusMenuBar$Help$On_Line_Manuals$Resample"
  if(is.element(menu.item.name, guiGetObjectNames("MenuItem")))
    guiRemove( "MenuItem", Name = menu.item.name)
  menu.item.name <- "SPlusMenuBar$Help$On_Line_Manuals$ResampleTutorial"
  if(is.element(menu.item.name, guiGetObjectNames("MenuItem")))
    guiRemove( "MenuItem", Name = menu.item.name)
  menu.item.name <- "SPlusMenuBar$Help$On_Line_Manuals$ResampleExamples"
  if(is.element(menu.item.name, guiGetObjectNames("MenuItem")))
    guiRemove( "MenuItem", Name = menu.item.name)
}



resampAddArgumentsToCall <- function(Call, additional.args){
  # Call should be a call to a function, like bootstrap(x, mean)
  # additional.args should be a string, e.g. "label='My Label'"
  # Additional arguments must be named.
  # Add the arguments to the call using
  # deparse, substring, paste, parse

  if(additional.args == "")
    return(Call)

  # Check that all additional arguments are given by name
  tmp <- names(parse(text=paste("f(", additional.args, ")"))[[1]])
  if(is.null(tmp) || sum(tmp=="") > 1)
    stop("Additional arguments must be given by name, e.g.",
	 "'block = 50' rather than '50'")

  tmp <- deparse( Call )
  tmp[length(tmp)] <- substring( tmp[length(tmp)], 1,(nchar(tmp[length(tmp)])-1))
  tmp[length(tmp)] <- paste( tmp[length(tmp)], ",", additional.args, ")", sep="" )
  parse(text=tmp)[[1]]
}



"Done with menufun.q" # for script file input
