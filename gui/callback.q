# This file contains
# (1) functions to simplify creating and modifying properties
#  newProperty
#  newFunctionInfo
#  copyProperty
#  newClassInfo
#  newContextMenu
#  newContextMenuItem

# (2) functions that support callback functions
#  cbIsNewOrChange
#  cbIsNew
#  cbIsChange
#  cbIsDone
#  cbSetEnableFlagL
#  callbackVarsFromDataFrame
#  cbGetColumnNames
#  callbackData2Resample

# The cbIs functions are based on building blocks:
#  cbIsInitDialogMessage
#  cbIsRollbackMessage
#  cbIsUpdateMessage
#  cbIsApplyMessage
#  cbIsOkMessage
#  cbIsCancelMessage


##################################################
# Define some functions to to simplify creating and modifying properties
# and menu items.  This uses a trick suggested by Chris Disdero
# to define properties cleanly - define, remove, then define

# The functions below use the following two names:
#  prop.file _ "guiResample.prp"
#  info.file _ "guiResample.fni"
##################################################

newProperty <- function(Name, SavePathName = "guiResample.prp", ...){
  # Call guiCreate for a Property
  # Don't know if it already existed, so first define then remove it.
  guiCreate("Property", Name = Name, SavePathName = SavePathName, ...)
  guiRemove("Property", Name = Name)
  guiCreate("Property", Name = Name, SavePathName = SavePathName, ...)
}

newFunctionInfo <- function(Name, SavePathName = "guiResample.fni", ...){
  # Call guiCreate for FunctionInfo
  # Don't know if it already existed, so first define then remove it.
  guiCreate("FunctionInfo", Name = Name, SavePathName = SavePathName, ...)
  guiRemove("FunctionInfo", Name = Name)
  guiCreate("FunctionInfo", Name = Name, SavePathName = SavePathName, ...)
}

copyProperty <- function(Name, OldName, SavePathName = "guiResample.prp", ...){
  # Call guiCopy for a Property
  # Name refers to the object being defined (like guiCreate, not guiCopy

  # Don't know if it already existed, so first define then remove it.
  guiCreate("Property", Name = Name, SavePathName = SavePathName, ...)
  guiRemove("Property", Name = Name)
  # Now copy to make the new object
  guiCopy("Property", Name = OldName, NewName = Name,
	  SavePathName = SavePathName, ...)
}


# Need to define ClassInfo, even for inheriting classes
newClassInfo <- function(class, action="print", ...){
  guiCreate( "ClassInfo", Name = class,
            ContextMenu = class,
            SavePathName = "",
            ImageFileName = "boot.bmp",
            DoubleClickAction = action, ...)
}
# use boot.bmp; hope to edit that file to be an actual boot
# current icon has no meaning to me.

newContextMenu <- function(class, ...){
  guiCreate( "MenuItem", Name = class,
            Type = "Menu",
            DocumentType = class, ...)
}

newContextMenuItem <- function(class, fun, text = fun){
  # e.g. "bootstrap", "limits.bca", "BCa confidence interval"
  guiCreate("MenuItem", Name = paste(class, fun, sep="$"),
            Type = "MenuItem",
            DocumentType = class,
            Action = "Function",
            Command = fun,
            MenuItemText = text,
            ShowDialogOnRun = F)
}


##################################################
# Define functions to simplify callback work
##################################################

cbIsNewOrChange <- function(dfProps, properties = ""){
  # True if called on dialog initialization, rollback, or
  # when updating any property listed in properties.
  cbIsNew(dfProps) || cbIsChange(dfProps, properties)
}

cbIsNew <- function(dfProps){
  # True if called on dialog initialization and rollback
  (cbIsInitDialogMessage(dfProps) ||
   cbIsRollbackMessage(dfProps))
}

cbIsChange <- function(dfProps, properties = ""){
  # True if update and one of properties is active (i.e. just changed)
  (cbIsUpdateMessage(dfProps) &&
   is.element(cbGetActiveProp(dfProps), properties))
}

cbIsDone <- function(dfProps){
  # True if OK or Apply was clicked
  # when updating any property listed in properties.
  (cbIsApplyMessage(dfProps) ||
   cbIsOkMessage(dfProps))
}


cbCopyCurrValue <- function(dfProps, to = "", from){
  # Copy the value from property `from' to properties `to'
  # like cbSetCurrValue, except last argument is the name of a property
  # Returns a modified value of dfProps
  cbSetEnableFlag(dfProps, to, cbGetCurrValue(dfProps, from))
}
cbCopyEnableFlag <- function(dfProps, to = "", from){
  # Copy the enable flag from property `from' to properties `to'
  # Returns a modified value of dfProps
  cbSetEnableFlag(dfProps, to, cbGetEnableFlag(dfProps, from))
}
cbAssignEnableFlag <- function(dfProps, to = "", from){
  # Set enable flag for properties `to' according to the
  # value (not the enable flag) for `from'
  # Returns a modified value of dfProps
  cbSetEnableFlag(dfProps, to, cbGetCurrValue(dfProps, from))
}

callbackVarsFromDataFrame <- function(dfp, dataFrameProperty, properties){
  # If the property `dataFrameProperty' points to a data frame,
  # set the optionList for all properties given in `properties'
  # to the variable names.
  # Then set the properties blank unless this is a rollback.
  dataName <- cbGetCurrValue(dfp, dataFrameProperty)
  x.names <- ""
  if(exists(dataName)) {
    data <- get(dataName)	     # This assumes the name is quoted
    if(is.data.frame(data))
      x.names <- paste(names(data), collapse = ",")
  }

  for(i in properties)
    dfp <- cbSetOptionList(dfp, i, x.names)
  if(!cbIsRollbackMessage(dfp)) {
    for(i in properties)
      dfp <- cbSetCurrValue(dfp, i, "")
  }
  dfp
}

cbGetColumnNames <- function(data.name, append, prepend, local = F){
  # Like cbGetColumnNamesString, except return a vector rather than
  # string with commas.
  # Also, replace parent argument with local (parent ignored by eval?)
  if(exists(data.name)) {
    name.vec <- names(eval(as.name(data.name), local = local))
    if(!missing(append))
      name.vec <- c(name.vec, append)
    if(!missing(prepend))
      name.vec <- c(prepend, name.vec)
  }
  else name.vec <- ""
  name.vec
}

callbackData2Resample <- function(dfp){
  # This is called by callBackPerm and callBackBootstrap
  # to handle second dataset fields (data or treatment)

  # If either of resamp_treatment or resamp_Data2 is not blank, disable other
  if(cbIsChange(dfp, "resamp_treatment") || cbIsRollbackMessage(dfp))
    dfp <- cbSetEnableFlag(dfp, "resamp_Data2",
			  cbGetCurrValue(dfp, "resamp_treatment") == "")
  if(cbIsChange(dfp, "resamp_Data2") || cbIsRollbackMessage(dfp))
    dfp <- cbSetEnableFlag(dfp, "resamp_treatment",
			  cbGetCurrValue(dfp, "resamp_Data2") == "")
  dfp
}



"Done with callback.q" # for script file input
