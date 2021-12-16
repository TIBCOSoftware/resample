#include <cdefs.h>
POUND @(#) Copyright (c), 2001, 2005 Insightful, Inc.  All rights reserved.

DEST=..

QFILES=  *.q

all : $(QFILES)
install : install.funs

install.funs : $(QFILES) $(DEST)/$(DATA_DIR)
	QINSTALL $(DEST) $(QFILES)




