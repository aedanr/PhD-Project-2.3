library(DSS)
library(edgeR)
rawdata <- readRDS('raw.counts.DE50.1.rds')
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
#   Lapack routine dgesv: system is exactly singular: U[3,3] = 0
# In addition: Warning messages:
# 1: glm.fit: fitted rates numerically 0 occurred 
# 2: In locfdr(normstat, plot = 0) :
#   f(z) misfit = 73.6.  Rerun with increased df
# 3: In locfdr(normstat, plot = 0) :
#   CM estimation failed, middle of histogram non-normal

dat.DSS.default <- newSeqCountSet(counts=matrix(rawdata, ncol=length(group)), 
                                  designs=as.numeric(group))
dat.DSS.default <- estNormFactors(dat.DSS.default)
dat.DSS.default <- estDispersion(dat.DSS.default)
res.DSS.default <- waldTest(dat.DSS.default, 1, 2)
# Warning messages:
# 1: glm.fit: fitted rates numerically 0 occurred 
# 2: In locfdr(normstat, plot = 0) :
#   f(z) misfit = 74.9.  Rerun with increased df

