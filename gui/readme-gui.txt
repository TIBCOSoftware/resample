Resample menus occur in the following places under the Statistics menu:
	Data Summaries
		Summary Statistics with Bootstrap	
		Correlations with Bootstrap		
	Compare Samples
		One Sample
			t Test with Bootstrap		
		Two Samples
			t Test with Bootstrap		
		Counts and Proportions
			Proportions Parameters with Bootstrap

	Regression
		Linear with Bootstrap			
	Resample
		Bootstrap
		Jackknife 
		Permutation tests
		Summary Statistics
		Correlations
		One-sample t
		Two-sample t
		Proportions
		Linear Regression

The first three under the resample menu are "generic menus" (can use
them for resampling with any statistic).

The final six are "special-purpose with supplemental tabs for
resampling".  These are identical to the six menus located under
Data Summaries, Compare Samples, and Regression (albeit with different names).

--------------------------------------------------
The menus are defined by the following files:

*** Define functions, not specific to one menu
tabplot.q
tabSummary.bootstrap.q
tabSummary.jackknife.q
callback.q		routines called by callback functions
menufun.q		from gui/, addResampleMenus(), removeResampleMenus()


*** Define properties (and a few functions)
--- Code to be reused by the 9 menus
guiResamp.ssc		properties & groups used in generic menus & supp. tabs
guiTabs.ssc		properties & groups used in supplemental tabs
-- Following six "special-purpose" use both guiResamp.ssc and guiTabs.ssc --
Cor.ssc			
Describe.ssc		
Linear.ssc		
Prop.ssc		
TTest1.ssc		
TTest2.ssc		
-- Following three "generic menus" use guiResamp.ssc --
guiBootstrap.ssc
guiJackknife.ssc
guiPermutations.ssc

--------------------------------------------------
Context menus (when in the Object Explorer and you right-click an object)

contextBootstrap.ssc
(still need to do jackknife and permutation test versions)
