# resamp objects
setOldClass(c("bootstrap", "resamp"), where="DESTDIR")
setOldClass(c("concomitants", "bootstrap", "resamp"), where="DESTDIR") # deprecated
setOldClass(c("concomitants.bootstrap", "bootstrap", "resamp"),
            where="DESTDIR")
setOldClass(c("controlVariates.bootstrap", "bootstrap", "resamp"),
            where="DESTDIR")
setOldClass(c("influence", "resamp"), where="DESTDIR")
setOldClass(c("jackknife", "resamp"), where="DESTDIR")
setOldClass(c("parametricBootstrapTest", "resamp"), where="DESTDIR")
setOldClass(c("parametricBootstrap", "resamp"), where="DESTDIR")
setOldClass(c("permutationTest", "resamp"), where="DESTDIR")
setOldClass(c("smoothedBootstrap", "resamp"), where="DESTDIR")
setOldClass(c("reweight", "resamp"), where="DESTDIR")
setOldClass(c("permutationTestMeans", "resamp"), where="DESTDIR")
setOldClass(c("permutationTest2", "permutationTest", "resamp"), where="DESTDIR")
setOldClass(c("bootstrap2", "bootstrap", "resamp"), where="DESTDIR")
setOldClass(c("limits.abc", "influence", "resamp"), where="DESTDIR")

# other objects
setOldClass(c("compressIndices", "matrix"))
setOldClass(c("htest.saddlepoint", "htest"), where="DESTDIR")
