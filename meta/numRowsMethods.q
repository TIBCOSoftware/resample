#
# Set model.list new-style method for numRows.
#
# numRows already a new-style generic in SV4.
#
setMethod("numRows","model.list",numRows.model.list)




