# Copyright 2021. TIBCO Software Inc.
# This file is subject to the license terms contained
# in the license file that is distributed with this file.

#
#  Kluge fix to Bug 24751 (was 23860) (help files are not searched
#  correctly when there are duplicate functions in different
#  libraries), to work with our library.
#
#  If resample is loaded, always return that help file.  This assumes
#  resample is loaded in first position, which it should be in order
#  to work.
#
#  This fix is for 6.1 for Windows.  6.0 is also broken, but always returns
#  the non-splus help file, so we'll take it.  Define find.doc.resample,
#  then rely on .First.lib to define find.doc as find.doc.resample if 6.1.
#
#  Problem is still there in S+6.2, S+7.0.
#  Starting with Resample beta 8, which only supports S+7 and later,
#  call this find.doc rather than find.doc.resample (previously .First copied
#  this to find.doc for S+6.1 and S+6.2).
"find.doc"<-
  function(what, defaultname = paste(getenv("SHOME"), "\\cmd\\splus.chm",
                   sep = ""), find.all = F, suppress.display = T, 
           hidden.keyword.search = F)
{
  if(getPlatform() != "WIN386") {
    paths <- search("paths")
    deflibs <- c("splus", "stat", "data", "trellis", "nlme3",
                 "main")
    defpath <- paste(getenv("SHOME"), "library", deflibs, sep = 
                     "/")
    defpath <- match(paths, defpath, 0)
    paths <- paths[!defpath]
    paths <- paths[paths != ""]
    paths <- paste(sep = "/", unlist(lapply(paths, database.path)), "__Jhelp")
    return(unix(paste(c("$SHOME/cmd/find.doc",
                        paste("'", as.character(what), "'", sep = ""),
                        paths), collapse = " ")))
  }

  # The rest of this function is for Windows.
  if(!is.name(what) && !is.character(what) || length(what) != 1)
    stop("what argument must be a character vector of length 1\n")
  sname <- what
  found.doc <- character(0)
  index.help <- 1
  
  if(length(defaultname) > 0) {
    index.help <- grep("splus.chm", defaultname)
  }
  
  if(!is.null(index.help) && length(index.help) > 0 && index.help[1] !=
     0) {
    found.doc <- c(oldFindDoc(sname), findHtml(sname))
    find.list <- character(0)
    if(find.all || length(found.doc) == 0) {
      for(search.list in c(T, F)) {
        sys.list <- sysfile.listing(pattern = 
                                    "\\.chm$", search.list = search.list)
        if(search.list != T) {
          index.help <- grep("splus.chm", 
                             sys.list)
          sys.list <- c(sys.list[index.help],
                        sys.list[1:length(sys.list) !=
                                 index.help])
        }
        if(access(defaultname, 0) == 0) {
          if(!is.na(i.default <- match.path(
                                            defaultname, sys.list)))
            sys.list <- c(defaultname,
                          sys.list[ - i.default]
                          )
        }
        find.list <- c(sys.list, find.list)
      }
      # Make sure that resample.chm is before splus.chm
      if(length(index.resample <- grep("resample.chm", find.list)))
	find.list <- c(find.list[index.resample], find.list[-index.resample])
    }
    if(getPlatform() == "WIN386") {
      if(length(find.list) <= 0) {
        warning("unable to locate any help (*.chm) files"
                )
      }
    }
    if(find.all && !suppress.display) {
      suppress.display <- T
      warning("'suppress.display = F' not allowed when 'find.all = T'"
              )
    }
    
    for(helpfile in find.list) {
      if(.C("S_query_for_help",
            helpfile,
            sname,
            as.logical(suppress.display),
            ret = logical(1),
            hidden.keyword.search)$ret) {
        found.doc <- c(found.doc, helpfile)
        if(!find.all)
          break
      }
    }
  }
  else {
    if(defaultname == "") {
      defaultname = paste(getenv("SHOME"), "\\cmd\\splus.chm", sep = "")
    }
    find.all <- F
    if(.C("S_query_for_help",
          defaultname,
          sname,
          as.logical(suppress.display),
          ret = logical(1),
          hidden.keyword.search)$ret) {
      found.doc <- defaultname
    }
  }
  
  #if(access(defaultname, 0) == 0 && (length(found.doc) == 0 || (find.all && is.na(match.path(defaultname, find.list)))))
  #  found.doc <- c(found.doc, defaultname)
  found.doc
}

