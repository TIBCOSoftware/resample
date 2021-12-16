#include <cdefs.h>

POUND #
POUND # Smakefile for building resample GUI
POUND # Created by copying and modifying SeqTrial version
POUND #
DEST=..
BUILDLOG=build.log
WINQFILES=callback.q first.q menufun.q tabSummary.bootstrap.q tabSummary.jackknife.q tabSummary.permutationTest.q tabplot.q 
GUIFILES=context.ssc contextBootstrap.ssc \
	Cor.ssc Describe.ssc Linear.ssc Prop.ssc TTest1.ssc TTest2.ssc \
	guiBootstrap.ssc guiJackknife.ssc guiPermutations.ssc guiResamp.ssc \
	guiTabs.ssc \
	remove.ssc
# remove.ssc should be last
ALLQS=allfiles.q
ALLGUI=allguifiles.q
ERRLOG=error.log

POUND #
POUND # For gui target, batch in the ssc files.  The result should be
POUND # gui.fni and gui.pp files in the _Prefs dir.
POUND #
POUND # all: is run first; it concatenates the .q files
POUND #	and the .ssc files
POUND # install: is run second; it calls install.funs (uses the .q files)
POUND #	then runs the .ssc files
all : $(ALLQS) $(ALLGUI)
install : gui

gui: install.funs $(PREFS_DIR)
	$(SBATCH) $(MODULEBLD)/gui/$(ALLGUI) $(DEST)/gui/$(ERRLOG)
	@if `grep -iq error $(ERRLOG)`; 	\
	  then echo "Errors occurred while building gui. See $(ERRLOG)";  \
	else touch gui; cd $(DEST)/$(PREFS_DIR) &&  find . \( -type f -a \! -name 'gui*' \) | xargs rm -f; \
	fi

install.funs : $(PREFS_DIR)
	$(SBATCH) $(MODULEBLD)/gui/$(ALLQS) $(DEST)/gui/$(ERRLOG)
	@if `grep -iq error $(ERRLOG)`; 	\
	  then echo "Errors occurred while building gui. See $(ERRLOG)";  \
	else touch install.funs; \
	fi

$(ALLQS): $(WINQFILES)
	cat $(WINQFILES) > $@

$(ALLGUI): $(GUIFILES)
	cat $(GUIFILES) > $@

$(PREFS_DIR):
	test -d $(DEST)/$(PREFS_DIR) || mkdir $(DEST)/$(PREFS_DIR)

virgin:
	rm -f $(ALLQS) install.funs gui
	rm -rf $(DEST)/$(PREFS_DIR)
