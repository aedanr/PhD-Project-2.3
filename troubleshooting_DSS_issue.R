library(DSS)
library(edgeR)
rawdata <- readRDS('raw.counts.DEDD50.3.rds')
group <- factor(c(rep(1,50), rep(2,50)))
libsizes <- colSums(rawdata)
nf.tmm <- calcNormFactors(rawdata, method="TMM")
els.tmm <- nf.tmm * libsizes
sf.tmm <- els.tmm / exp(mean(log(libsizes)))

dat.DSS.tmm <- newSeqCountSet(counts=matrix(rawdata, ncol=length(group)), 
                              designs=as.numeric(group), 
                              normalizationFactor=sf.tmm)
dat.DSS.tmm <- estDispersion(dat.DSS.tmm)
res.DSS.tmm <- waldTest(dat.DSS.tmm, 1, 2)
# Error in solve.default(G0) : 
#   Lapack routine dgesv: system is computationally singular: reciprocal condition number = 4.47602e-18
# In addition: Warning messages:
#   1: In locfdr(normstat, plot = 0) :
#   f(z) misfit = 114.9.  Rerun with increased df
# 2: In locfdr(normstat, plot = 0) :
#   CM estimation failed, middle of histogram non-normal

dat.DSS.default <- newSeqCountSet(counts=matrix(rawdata, ncol=length(group)), 
                                  designs=as.numeric(group))
dat.DSS.default <- estNormFactors(dat.DSS.default)
dat.DSS.default <- estDispersion(dat.DSS.default)
res.DSS.default <- waldTest(dat.DSS.default, 1, 2)
# Error in solve.default(G0) : 
#   system is computationally singular: reciprocal condition number = 1.49125e-18
# In addition: Warning messages:
#   1: In locfdr(normstat, plot = 0) :
#   f(z) misfit = 114.8.  Rerun with increased df
# 2: In locfdr(normstat, plot = 0) :
#   CM estimation failed, middle of histogram non-normal

