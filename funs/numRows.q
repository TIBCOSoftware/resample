# Copyright 2021. TIBCO Software Inc.
# This file is subject to the license terms contained
# in the license file that is distributed with this file.

#
# model.list method for numRows 
#
# numRows is already a new-style generic in SV4.  
#
# Define this as new-style method in meta/numRowsMethods.q
#
# 3/7/01 (RT)
#
numRows.model.list <- function(x,columns) numRows(x$data$X)

