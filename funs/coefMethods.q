# Copyright 2021. TIBCO Software Inc.
# This file is subject to the license terms contained
# in the license file that is distributed with this file.

#
# censorRegList method for coef
#
# Without this, coef.lm is the method for censorReg objects, which 
# fails on clas censorRegList.
#

coef.censorRegList <-
function(object)
  unlist(lapply(object, "coef"))

