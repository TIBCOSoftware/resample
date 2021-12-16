#include <cdefs.h>
POUND @(#) $RCSfile: Smake.cpp,v $: $Revision: #7 $, $Date: 2010/10/21 $
POUND Copyright (c), 2001 Insightful, Inc.  All rights reserved.

DEST=.
LIBRARY=resample

COMMON_TARGS=
COMMON_DIRS=funs src meta data sgml

POUND S7 CHANGE For Windows 64-bit we don't need to build gui directory

#if defined(WIN386)
#if defined(WIN64) || defined(SPLUS_CONTINUOUS_BUILD)
DIRS=$(COMMON_DIRS) 
#else
DIRS=$(COMMON_DIRS) gui
#endif
#endif
#if defined(WIN386)/*(*/
EXTRA_TARGS=$(COMMON_TARGS) Uglist
VERSION_PATH=$(CMDDIR)
FUNS_DIRS=funs meta
#else /*)(*/
DIRS=$(COMMON_DIRS)
EXTRA_TARGS=$(COMMON_TARGS)
VERSION_PATH=$(SHOME)/adm/cmd
FUNS_DIRS=funs meta
#endif /*)*/

DATA_DIR=.Data

default: all install

install: install.subs install.all $(EXTRA_TARGS)
install.funs: install.subs $(EXTRA_TARGS)

all :
	@for dir in $(DIRS); do \
	  (cd $$dir; echo ===== Making $@ for $(LIBRARY) in $$dir; $(S_MAKE) $@) || exit 1 ; \
	done

install.subs:
	@for dir in $(FUNS_DIRS); do \
	  (cd $$dir; echo ===== Making $* for $(LIBRARY) in `pwd`; $(S_MAKE) $*) || exit 1 ; \
	done

install.all:
	@for dir in $(DIRS); do \
	  (cd $$dir; echo ===== Making install for $(LIBRARY) in `pwd`; $(S_MAKE) install) || exit 1 ; \
	done

virgin:
POUND leave virgin blank; for sysgen-update only.
