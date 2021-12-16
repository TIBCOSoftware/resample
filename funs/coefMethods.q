#
# censorRegList method for coef
#
# Without this, coef.lm is the method for censorReg objects, which 
# fails on clas censorRegList.
#

coef.censorRegList <-
function(object)
  unlist(lapply(object, "coef"))

