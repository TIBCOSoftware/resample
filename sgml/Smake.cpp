#include <cdefs.h>

DEST=..
LIBRARY=resample

HFILEDEFAULT=$(MODULEBLD)/gui/help/Resample

HFILES= Defunct.sgml Deprecated.sgml addSamples.sgml \
	balancedSample.sgml bootPred.sgml bootstats.sgml \
	bootstrap.args.sgml bootstrap.censorReg.sgml \
	bootstrap.sgml bootstrap.lm.sgml \
	bootstrap2.sgml bootstrapMeans.sgml bootstrapT.sgml \
	cdf.sgml cdf.resamp.sgml censorReg.dofit.sgml \
	censrreg.sgml coef.censorRegList.sgml colTabulate.sgml \
	combinePValues.sgml \
	compressIndices.sgml \
	concomitants.sgml concomitants.bootstrap.sgml \
	controlVariates.sgml controlVariates.bootstrap.sgml \
	crossVal.sgml \
	frame.eval.sgml glm.sgml \
	influence.sgml  \
	jack.after.bootstrap.sgml jackknife.sgml jackknife.lm.sgml \
	limits.abc.sgml limits.bca.sgml \
	limits.percentile.sgml limits.t.sgml limits.tilt.sgml \
	linearApproxReg.sgml lm.sgml \
	lm.model.list.sgml model.list.obj.sgml \
	pbootTest.sgml pbootstrap.sgml \
	permutationTest.sgml permutationTest2.sgml permutationTestMeans.sgml \
	plot.jack.after.bootstrap.sgml plot.resamp.sgml \
	plot.tilt.after.bootstrap.sgml \
	plotCDF.sgml plotCDF.resamp.sgml \
	print.jack.after.bootstrap.sgml print.resamp.sgml \
	print.summary.bootstrap.sgml print.summary.pbootstrap.sgml \
	print.summary.resamp.sgml print.summary.sbootstrap.sgml \
	qqnorm.resamp.sgml \
	randomSample.sgml resamp.object.sgml resamp.problems.sgml \
	resampFunctionalList.sgml \
	resampGet.sgml resampGetL.sgml \
	resampMake.sgml resampMakeFitObjResults.sgml \
	resampMakeNewStat.sgml resampPivot.sgml resampleOptions.sgml \
	resampleOverview.sgml \
	revSaddlepointP.sgml \
	reweight.sgml saddlepointP.sgml samp.sgml sbootstrap.sgml \
	summary.bootstrap.sgml summary.pbootstrap.sgml \
	summary.resamp.sgml summary.sbootstrap.sgml t.test.sgml \
	tilt.after.bootstrap.sgml \
	tiltBootProbs.sgml tiltDetails.sgml tiltMean.sgml tiltWeights.sgml \
	Verizon.sgml

all : $(HFILES)

install : $(DEST)/$(DATA_DIR) help

#if defined(WIN386)
help : $(HFILES)
	@for file in $(HFILES); do \
		( ( test -d bldtmp || mkdir bldtmp ) && cp $$file bldtmp/$$file ) || exit 1 ; \
	done
	( chmod -R u+w bldtmp )
	$(CHMTOOLS)/build_helpfiles -a $(LIBRARY) $(MODULEBLD)/sgml/bldtmp $(MODULEBLD)/sgml/bldtmp $(CHMTOOLS) $(CYG_HTML_HELPDIR) $(HFILEDEFAULT).html $(MODULEBLD) $(CYG_CYGWIN_BIN) "$(LIBRARY) Reference" F notused $(MODULEBLD)/gui/help F
#else
help : $(HFILES)
	HINSTALL $(DEST)/.Data $(HFILES)
	cd $(DEST); BUILD_JHELP -nogifs $(LIBRARY)
#endif

update :

virgin : clean virgin.std
