#include <cdefs.h>

POUND Copyright (c), 2001, 2005 Insightful, Inc.  All rights reserved.

DEST=..

C_SOURCES= basic.c \
           compressIndices.c \
           pTestMeans.c revSaddleSolveExp.c \
           revSaddleSolveMl.c saddle.c  \
           tilt.c tiltBoot.c tiltMean.c 

  /* Some .c files depend on tilt.h */
  /* That dependency and rule is not encoded here, but should be. */

OBJS=$(C_SOURCES:.c=.o)

/* TODO: Routine names should be listed here...  (is this important?)*/
C_OBJ_NAMES=

all : $(OBJS)

install : update
update : $(OBJS)
	$(C)/LIBRARY $(DEST)/$(S_SO) $(OBJS) $(LIBS)

virgin : clean virgin.std

