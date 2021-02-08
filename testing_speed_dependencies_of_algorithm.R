library(compcodeR)
counts <- generateSyntheticData(dataset='data', n.vars=1000, samples.per.cond=10, n.diffexp=100)
library(DESeq2)
group <- factor(c(rep(1,10), rep(2,10)))
dat.DESeq <- DESeqDataSetFromMatrix(countData=counts@count.matrix, colData=data.frame(group), design=~group)
dat.DESeq <- estimateSizeFactors(dat.DESeq)
counts <- t(counts@count.matrix) / dat.DESeq$sizeFactor
rm(dat.DESeq)
system.time(ln_hmm_adapt_3_chains(counts=counts, groups=group, chain.length=1))
# 60.07
system.time(ln_hmm_adapt_3_chains(counts=counts, groups=group, chain.length=10))
# 60.59
system.time(ln_hmm_adapt_3_chains(counts=counts, groups=group, chain.length=100))
# 82.08
system.time(ln_hmm_adapt_3_chains(counts=counts, groups=group, chain.length=1000))
# 160.08
system.time(ln_hmm_adapt_3_chains(counts=counts, groups=group, chain.length=10000))
# 1006.00
g1000_i20 <- c(60.07, 60.59, 82.08, 160.08, 1006.00)


counts <- generateSyntheticData(dataset='data', n.vars=100, samples.per.cond=10, n.diffexp=100)
group <- factor(c(rep(1,10), rep(2,10)))
dat.DESeq <- DESeqDataSetFromMatrix(countData=counts@count.matrix, colData=data.frame(group), design=~group)
dat.DESeq <- estimateSizeFactors(dat.DESeq)
counts <- t(counts@count.matrix) / dat.DESeq$sizeFactor
rm(dat.DESeq)
system.time(ln_hmm_adapt_3_chains(counts=counts, groups=group, chain.length=1))
# 9.34
system.time(ln_hmm_adapt_3_chains(counts=counts, groups=group, chain.length=10))
# 9.06
system.time(ln_hmm_adapt_3_chains(counts=counts, groups=group, chain.length=100))
# 10.58
system.time(ln_hmm_adapt_3_chains(counts=counts, groups=group, chain.length=1000))
# 22.25
system.time(ln_hmm_adapt_3_chains(counts=counts, groups=group, chain.length=10000))
# 131.41
g100_i20 <- c(9.34, 9.06, 10.58, 22.25, 131.41)


counts <- generateSyntheticData(dataset='data', n.vars=10, samples.per.cond=10, n.diffexp=10)
group <- factor(c(rep(1,10), rep(2,10)))
dat.DESeq <- DESeqDataSetFromMatrix(countData=counts@count.matrix, colData=data.frame(group), design=~group)
dat.DESeq <- estimateSizeFactors(dat.DESeq)
counts <- t(counts@count.matrix) / dat.DESeq$sizeFactor
rm(dat.DESeq)
system.time(ln_hmm_adapt_3_chains(counts=counts, groups=group, chain.length=1))
# 3.07
system.time(ln_hmm_adapt_3_chains(counts=counts, groups=group, chain.length=10))
# 3.17
system.time(ln_hmm_adapt_3_chains(counts=counts, groups=group, chain.length=100))
# 3.16
system.time(ln_hmm_adapt_3_chains(counts=counts, groups=group, chain.length=1000))
# 6.75
system.time(ln_hmm_adapt_3_chains(counts=counts, groups=group, chain.length=10000))
# 43.35
g10_i20 <- c(3.07, 3.17, 3.16, 6.75, 43.35)


# 10 genes to 100 genes increases ~ 3x
# 100 genes to 1000 genes increases ~7x
# 100 to 1000 iterations increases ~2x
# 1000 to 10,000 iterations increases ~6x
# After accounting for baseline, increases ~linearly with number of iterations



counts <- generateSyntheticData(dataset='data', n.vars=100, samples.per.cond=2, n.diffexp=0)
group <- factor(c(rep(1,2), rep(2,2)))
dat.DESeq <- DESeqDataSetFromMatrix(countData=counts@count.matrix, colData=data.frame(group), design=~group)
dat.DESeq <- estimateSizeFactors(dat.DESeq)
counts <- t(counts@count.matrix) / dat.DESeq$sizeFactor
rm(dat.DESeq)
sample.means0 <- pmax(colMeans(counts), 0.01)
sample.vars0 <- apply(counts, 2, var)
sample.disps0 <- pmax(sample.vars0 - sample.means0, 0.01) / pmax(sample.means0^2, 0.1)
inits=list("means0" = sample.means0, "means1" = sample.means0, "means2" = sample.means0, 
           "disps0" = sample.disps0, "disps1" = sample.disps0, "disps2" = sample.disps0, 
           "mean.prior.location" = 1, "disp.prior.location" = 1, "mean.prior.scale" = 1, "disp.prior.scale" = 1)
system.time(ln_hmm_1_chain(counts=counts, groups=group, inits=inits, chain.length=100))
# 0.65
system.time(ln_hmm_1_chain(counts=counts, groups=group, inits=inits, chain.length=1000))
# 5.69
system.time(ln_hmm_1_chain(counts=counts, groups=group, inits=inits, chain.length=10000))
# 49.14

counts <- generateSyntheticData(dataset='data', n.vars=100, samples.per.cond=10, n.diffexp=0)
group <- factor(c(rep(1,10), rep(2,10)))
dat.DESeq <- DESeqDataSetFromMatrix(countData=counts@count.matrix, colData=data.frame(group), design=~group)
dat.DESeq <- estimateSizeFactors(dat.DESeq)
counts <- t(counts@count.matrix) / dat.DESeq$sizeFactor
rm(dat.DESeq)
sample.means0 <- pmax(colMeans(counts), 0.01)
sample.vars0 <- apply(counts, 2, var)
sample.disps0 <- pmax(sample.vars0 - sample.means0, 0.01) / pmax(sample.means0^2, 0.1)
inits=list("means0" = sample.means0, "means1" = sample.means0, "means2" = sample.means0, 
           "disps0" = sample.disps0, "disps1" = sample.disps0, "disps2" = sample.disps0, 
           "mean.prior.location" = 1, "disp.prior.location" = 1, "mean.prior.scale" = 1, "disp.prior.scale" = 1)
system.time(ln_hmm_1_chain(counts=counts, groups=group, inits=inits, chain.length=100))
# 0.81
system.time(ln_hmm_1_chain(counts=counts, groups=group, inits=inits, chain.length=1000))
# 5.93
system.time(ln_hmm_1_chain(counts=counts, groups=group, inits=inits, chain.length=10000))
# 56.20

counts <- generateSyntheticData(dataset='data', n.vars=100, samples.per.cond=50, n.diffexp=0)
group <- factor(c(rep(1,50), rep(2,50)))
dat.DESeq <- DESeqDataSetFromMatrix(countData=counts@count.matrix, colData=data.frame(group), design=~group)
dat.DESeq <- estimateSizeFactors(dat.DESeq)
counts <- t(counts@count.matrix) / dat.DESeq$sizeFactor
rm(dat.DESeq)
sample.means0 <- pmax(colMeans(counts), 0.01)
sample.vars0 <- apply(counts, 2, var)
sample.disps0 <- pmax(sample.vars0 - sample.means0, 0.01) / pmax(sample.means0^2, 0.1)
inits=list("means0" = sample.means0, "means1" = sample.means0, "means2" = sample.means0, 
           "disps0" = sample.disps0, "disps1" = sample.disps0, "disps2" = sample.disps0, 
           "mean.prior.location" = 1, "disp.prior.location" = 1, "mean.prior.scale" = 1, "disp.prior.scale" = 1)
system.time(ln_hmm_1_chain(counts=counts, groups=group, inits=inits, chain.length=100))
# 1.08
system.time(ln_hmm_1_chain(counts=counts, groups=group, inits=inits, chain.length=1000))
# 9.66
system.time(ln_hmm_1_chain(counts=counts, groups=group, inits=inits, chain.length=10000))
# 96.10

counts <- generateSyntheticData(dataset='data', n.vars=100, samples.per.cond=250, n.diffexp=0)
group <- factor(c(rep(1,250), rep(2,250)))
dat.DESeq <- DESeqDataSetFromMatrix(countData=counts@count.matrix, colData=data.frame(group), design=~group)
dat.DESeq <- estimateSizeFactors(dat.DESeq)
counts <- t(counts@count.matrix) / dat.DESeq$sizeFactor
rm(dat.DESeq)
sample.means0 <- pmax(colMeans(counts), 0.01)
sample.vars0 <- apply(counts, 2, var)
sample.disps0 <- pmax(sample.vars0 - sample.means0, 0.01) / pmax(sample.means0^2, 0.1)
inits=list("means0" = sample.means0, "means1" = sample.means0, "means2" = sample.means0, 
           "disps0" = sample.disps0, "disps1" = sample.disps0, "disps2" = sample.disps0, 
           "mean.prior.location" = 1, "disp.prior.location" = 1, "mean.prior.scale" = 1, "disp.prior.scale" = 1)
system.time(ln_hmm_1_chain(counts=counts, groups=group, inits=inits, chain.length=100))
# 3.08
system.time(ln_hmm_1_chain(counts=counts, groups=group, inits=inits, chain.length=1000))
# 28.75
system.time(ln_hmm_1_chain(counts=counts, groups=group, inits=inits, chain.length=10000))
# 288.75


counts <- generateSyntheticData(dataset='data', n.vars=1000, samples.per.cond=2, n.diffexp=0)
group <- factor(c(rep(1,2), rep(2,2)))
dat.DESeq <- DESeqDataSetFromMatrix(countData=counts@count.matrix, colData=data.frame(group), design=~group)
dat.DESeq <- estimateSizeFactors(dat.DESeq)
counts <- t(counts@count.matrix) / dat.DESeq$sizeFactor
rm(dat.DESeq)
sample.means0 <- pmax(colMeans(counts), 0.01)
sample.vars0 <- apply(counts, 2, var)
sample.disps0 <- pmax(sample.vars0 - sample.means0, 0.01) / pmax(sample.means0^2, 0.1)
inits=list("means0" = sample.means0, "means1" = sample.means0, "means2" = sample.means0, 
           "disps0" = sample.disps0, "disps1" = sample.disps0, "disps2" = sample.disps0, 
           "mean.prior.location" = 1, "disp.prior.location" = 1, "mean.prior.scale" = 1, "disp.prior.scale" = 1)
system.time(ln_hmm_1_chain(counts=counts, groups=group, inits=inits, chain.length=100))
# 4.29
system.time(ln_hmm_1_chain(counts=counts, groups=group, inits=inits, chain.length=1000))
# 40.43
system.time(ln_hmm_1_chain(counts=counts, groups=group, inits=inits, chain.length=10000))
# 396.96

counts <- generateSyntheticData(dataset='data', n.vars=1000, samples.per.cond=10, n.diffexp=0)
group <- factor(c(rep(1,10), rep(2,10)))
dat.DESeq <- DESeqDataSetFromMatrix(countData=counts@count.matrix, colData=data.frame(group), design=~group)
dat.DESeq <- estimateSizeFactors(dat.DESeq)
counts <- t(counts@count.matrix) / dat.DESeq$sizeFactor
rm(dat.DESeq)
sample.means0 <- pmax(colMeans(counts), 0.01)
sample.vars0 <- apply(counts, 2, var)
sample.disps0 <- pmax(sample.vars0 - sample.means0, 0.01) / pmax(sample.means0^2, 0.1)
inits=list("means0" = sample.means0, "means1" = sample.means0, "means2" = sample.means0, 
           "disps0" = sample.disps0, "disps1" = sample.disps0, "disps2" = sample.disps0, 
           "mean.prior.location" = 1, "disp.prior.location" = 1, "mean.prior.scale" = 1, "disp.prior.scale" = 1)
system.time(ln_hmm_1_chain(counts=counts, groups=group, inits=inits, chain.length=100))
# 5.26
system.time(ln_hmm_1_chain(counts=counts, groups=group, inits=inits, chain.length=1000))
# 49.90
system.time(ln_hmm_1_chain(counts=counts, groups=group, inits=inits, chain.length=10000))
# 571.45

counts <- generateSyntheticData(dataset='data', n.vars=1000, samples.per.cond=50, n.diffexp=0)
group <- factor(c(rep(1,50), rep(2,50)))
dat.DESeq <- DESeqDataSetFromMatrix(countData=counts@count.matrix, colData=data.frame(group), design=~group)
dat.DESeq <- estimateSizeFactors(dat.DESeq)
counts <- t(counts@count.matrix) / dat.DESeq$sizeFactor
rm(dat.DESeq)
sample.means0 <- pmax(colMeans(counts), 0.01)
sample.vars0 <- apply(counts, 2, var)
sample.disps0 <- pmax(sample.vars0 - sample.means0, 0.01) / pmax(sample.means0^2, 0.1)
inits=list("means0" = sample.means0, "means1" = sample.means0, "means2" = sample.means0, 
           "disps0" = sample.disps0, "disps1" = sample.disps0, "disps2" = sample.disps0, 
           "mean.prior.location" = 1, "disp.prior.location" = 1, "mean.prior.scale" = 1, "disp.prior.scale" = 1)
system.time(ln_hmm_1_chain(counts=counts, groups=group, inits=inits, chain.length=100))
# 9.44
system.time(ln_hmm_1_chain(counts=counts, groups=group, inits=inits, chain.length=1000))
# 91.51
system.time(ln_hmm_1_chain(counts=counts, groups=group, inits=inits, chain.length=10000))
# 894.67

counts <- generateSyntheticData(dataset='data', n.vars=1000, samples.per.cond=250, n.diffexp=0)
group <- factor(c(rep(1,250), rep(2,250)))
dat.DESeq <- DESeqDataSetFromMatrix(countData=counts@count.matrix, colData=data.frame(group), design=~group)
dat.DESeq <- estimateSizeFactors(dat.DESeq)
counts <- t(counts@count.matrix) / dat.DESeq$sizeFactor
rm(dat.DESeq)
sample.means0 <- pmax(colMeans(counts), 0.01)
sample.vars0 <- apply(counts, 2, var)
sample.disps0 <- pmax(sample.vars0 - sample.means0, 0.01) / pmax(sample.means0^2, 0.1)
inits=list("means0" = sample.means0, "means1" = sample.means0, "means2" = sample.means0, 
           "disps0" = sample.disps0, "disps1" = sample.disps0, "disps2" = sample.disps0, 
           "mean.prior.location" = 1, "disp.prior.location" = 1, "mean.prior.scale" = 1, "disp.prior.scale" = 1)
system.time(ln_hmm_1_chain(counts=counts, groups=group, inits=inits, chain.length=100))
# 26.97
system.time(ln_hmm_1_chain(counts=counts, groups=group, inits=inits, chain.length=1000))
# 279.12
system.time(ln_hmm_1_chain(counts=counts, groups=group, inits=inits, chain.length=10000))
# 2773.05

counts <- generateSyntheticData(dataset='data', n.vars=10000, samples.per.cond=2, n.diffexp=0)
group <- factor(c(rep(1,2), rep(2,2)))
dat.DESeq <- DESeqDataSetFromMatrix(countData=counts@count.matrix, colData=data.frame(group), design=~group)
dat.DESeq <- estimateSizeFactors(dat.DESeq)
counts <- t(counts@count.matrix) / dat.DESeq$sizeFactor
rm(dat.DESeq)
sample.means0 <- pmax(colMeans(counts), 0.01)
sample.vars0 <- apply(counts, 2, var)
sample.disps0 <- pmax(sample.vars0 - sample.means0, 0.01) / pmax(sample.means0^2, 0.1)
inits=list("means0" = sample.means0, "means1" = sample.means0, "means2" = sample.means0, 
           "disps0" = sample.disps0, "disps1" = sample.disps0, "disps2" = sample.disps0, 
           "mean.prior.location" = 1, "disp.prior.location" = 1, "mean.prior.scale" = 1, "disp.prior.scale" = 1)
system.time(ln_hmm_1_chain(counts=counts, groups=group, inits=inits, chain.length=100))
# 35.08
system.time(ln_hmm_1_chain(counts=counts, groups=group, inits=inits, chain.length=1000))
# 356.83
system.time(ln_hmm_1_chain(counts=counts, groups=group, inits=inits, chain.length=10000))
# 3718.22

counts <- generateSyntheticData(dataset='data', n.vars=10000, samples.per.cond=10, n.diffexp=0)
group <- factor(c(rep(1,10), rep(2,10)))
dat.DESeq <- DESeqDataSetFromMatrix(countData=counts@count.matrix, colData=data.frame(group), design=~group)
dat.DESeq <- estimateSizeFactors(dat.DESeq)
counts <- t(counts@count.matrix) / dat.DESeq$sizeFactor
rm(dat.DESeq)
sample.means0 <- pmax(colMeans(counts), 0.01)
sample.vars0 <- apply(counts, 2, var)
sample.disps0 <- pmax(sample.vars0 - sample.means0, 0.01) / pmax(sample.means0^2, 0.1)
inits=list("means0" = sample.means0, "means1" = sample.means0, "means2" = sample.means0, 
           "disps0" = sample.disps0, "disps1" = sample.disps0, "disps2" = sample.disps0, 
           "mean.prior.location" = 1, "disp.prior.location" = 1, "mean.prior.scale" = 1, "disp.prior.scale" = 1)
system.time(ln_hmm_1_chain(counts=counts, groups=group, inits=inits, chain.length=100))
# 44.83
system.time(ln_hmm_1_chain(counts=counts, groups=group, inits=inits, chain.length=1000))
# 448.45
system.time(ln_hmm_1_chain(counts=counts, groups=group, inits=inits, chain.length=10000))
# 4298.47

counts <- generateSyntheticData(dataset='data', n.vars=10000, samples.per.cond=50, n.diffexp=0)
group <- factor(c(rep(1,50), rep(2,50)))
dat.DESeq <- DESeqDataSetFromMatrix(countData=counts@count.matrix, colData=data.frame(group), design=~group)
dat.DESeq <- estimateSizeFactors(dat.DESeq)
counts <- t(counts@count.matrix) / dat.DESeq$sizeFactor
rm(dat.DESeq)
sample.means0 <- pmax(colMeans(counts), 0.01)
sample.vars0 <- apply(counts, 2, var)
sample.disps0 <- pmax(sample.vars0 - sample.means0, 0.01) / pmax(sample.means0^2, 0.1)
inits=list("means0" = sample.means0, "means1" = sample.means0, "means2" = sample.means0, 
           "disps0" = sample.disps0, "disps1" = sample.disps0, "disps2" = sample.disps0, 
           "mean.prior.location" = 1, "disp.prior.location" = 1, "mean.prior.scale" = 1, "disp.prior.scale" = 1)
system.time(ln_hmm_1_chain(counts=counts, groups=group, inits=inits, chain.length=100))
# 80.03
system.time(ln_hmm_1_chain(counts=counts, groups=group, inits=inits, chain.length=1000))
# 804.06
system.time(ln_hmm_1_chain(counts=counts, groups=group, inits=inits, chain.length=10000))
# 7921.67

counts <- generateSyntheticData(dataset='data', n.vars=10000, samples.per.cond=250, n.diffexp=0)
group <- factor(c(rep(1,250), rep(2,250)))
dat.DESeq <- DESeqDataSetFromMatrix(countData=counts@count.matrix, colData=data.frame(group), design=~group)
dat.DESeq <- estimateSizeFactors(dat.DESeq)
counts <- t(counts@count.matrix) / dat.DESeq$sizeFactor
rm(dat.DESeq)
sample.means0 <- pmax(colMeans(counts), 0.01)
sample.vars0 <- apply(counts, 2, var)
sample.disps0 <- pmax(sample.vars0 - sample.means0, 0.01) / pmax(sample.means0^2, 0.1)
inits=list("means0" = sample.means0, "means1" = sample.means0, "means2" = sample.means0, 
           "disps0" = sample.disps0, "disps1" = sample.disps0, "disps2" = sample.disps0, 
           "mean.prior.location" = 1, "disp.prior.location" = 1, "mean.prior.scale" = 1, "disp.prior.scale" = 1)
system.time(ln_hmm_1_chain(counts=counts, groups=group, inits=inits, chain.length=100))
# 268.82
system.time(ln_hmm_1_chain(counts=counts, groups=group, inits=inits, chain.length=1000))
# 2684.75
system.time(ln_hmm_1_chain(counts=counts, groups=group, inits=inits, chain.length=10000))
# 26743.11





