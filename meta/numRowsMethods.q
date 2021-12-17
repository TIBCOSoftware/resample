# Copyright 2021. TIBCO Software Inc.
# This file is subject to the license terms contained
# in the license file that is distributed with this file.

#
# Set model.list new-style method for numRows.
#
# numRows already a new-style generic in SV4.
#
setMethod("numRows","model.list",numRows.model.list)




