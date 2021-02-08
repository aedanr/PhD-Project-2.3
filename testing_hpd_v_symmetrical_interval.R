hpd.pval <- function(x, m=0) {
  require(HDInterval)
  y1 <- numeric(9)
  for (t in 1:9) {
    y1[t] <- sign(hdi(x, credMass = 1 - t/10)[1] - m) == sign(hdi(x, credMass = 1 - t/10)[2] + m)
    if (y1[t] == 1) {break}
  }
  if (sum(y1) == 0) {x1 <- 9} else {x1 <- min(which(y1 == 1)) - 1}
  iter <- x1/10
  y2 <- numeric(9)
  for (t in 1:9) {
    y2[t] <- sign(hdi(x, credMass = 1 - (iter + t/100))[1] - m) == sign(hdi(x, credMass = 1 - (iter + t/100))[2] + m)
    if (y2[t] == 1) {break}
  }
  if (sum(y2) == 0) {x2 <- 9} else {x2 <- min(which(y2 == 1)) - 1}
  iter <- iter + x2/100
  y3 <- numeric(9)
  for (t in 1:9) {
    y3[t] <- sign(hdi(x, credMass = 1 - (iter + t/1000))[1] - m) == sign(hdi(x, credMass = 1 - (iter + t/1000))[2] + m)
    if (y3[t] == 1) {break}
  }
  if (sum(y3) == 0) {x3 <- 9} else {x3 <- min(which(y3 == 1)) - 1}
  iter <- iter + x3/1000
  y4 <- numeric(9)
  for (t in 1:9) {
    y4[t] <- sign(hdi(x, credMass = 1 - (iter + t/10000))[1] - m) == sign(hdi(x, credMass = 1 - (iter + t/10000))[2] + m)
    if (y4[t] == 1) {break}
  }
  if (sum(y4) == 0) {x4 <- 10} else {x4 <- min(which(y4 == 1))}
  iter <- iter + x4/10000
  return(iter)
}

# hdi() computes highest density interval for a probability distribution for a given probability mass.
# Aim is to find density of smallest interval that doesn't contain zero.
# If used a symmetric interval, could just count proportion of samples above and below zero and 
# use 2*min as a two-sided p-value, but if don't want to assume symmetry, need to find boundaries of 
# highest density interval

library(here)
library(HDInterval)
library(compcodeR)
library(ROCR)
source(here('scripts','2019-03-27_lognormal_hmm_adaptive_proposals_three_chains_function.R'))
counts <- generateSyntheticData(dataset='DEonly', n.vars=2400, samples.per.cond=5, 
                                n.diffexp=120, fraction.upregulated=0.5,
                                filter.threshold.mediancpm=0.5)
DE <- counts@variable.annotations$differential.expression
lnHM <- ln_hmm_adapt_3_chains(counts=t(counts@count.matrix), groups=c(rep(1,5),rep(2,5)))
samples <- log(as.matrix(lnHM$means1)) - log(as.matrix(lnHM$means2))
plot(density(samples[,1]))
plot(density(samples[,2]))
hpd.pval(samples[,1]) # 0.0001
hpd.pval(samples[,2]) # 0.937
hpd.pval(samples[,14]) # 0.0012
test <- samples[,14] - 0.1
hpd.pval(test) # 0.0012 # 0.0002
system.time(apply(samples, 2, hpd.pval)) #  46 s for 2400 genes; > 5 min for 20,000
2 * min(mean(samples[,1] < 0), mean(samples[,1] > 0)) # 0
2 * min(mean(samples[,2] < 0), mean(samples[,2] > 0)) # 0.813
2 * min(mean(samples[,14] < 0), mean(samples[,14] > 0)) # 0.0017
2 * min(mean(test < 0), mean(test > 0)) # 0.0003333

p.hpd <- apply(samples, 2, hpd.pval)
p.sym <- apply(samples, 2, function(x) {2 * min(mean(x < 0), mean(x > 0))})
performance(prediction(1 - p.hpd, DE), measure='auc')@y.values[[1]]
performance(prediction(1 - p.sym, DE), measure='auc')@y.values[[1]]


auc.hpd <- numeric(10)
auc.sym <- numeric(10)
fdr.hpd <- numeric(10)
fdr.sym <- numeric(10)
for (i in 1:10) {
  counts <- generateSyntheticData(dataset='DEonly', n.vars=10000, samples.per.cond=5, 
                                  n.diffexp=200, fraction.upregulated=0.5,
                                  filter.threshold.mediancpm=0.5)
  DE <- counts@variable.annotations$differential.expression
  lnHM <- ln_hmm_adapt_3_chains(counts=t(counts@count.matrix), groups=c(rep(1,5),rep(2,5)))
  samples <- log(as.matrix(lnHM$means1)) - log(as.matrix(lnHM$means2))
  p.hpd <- apply(samples, 2, hpd.pval)
  p.sym <- apply(samples, 2, function(x) {2 * min(mean(x < 0), mean(x > 0))})
  auc.hpd[i] <- performance(prediction(1 - p.hpd, DE), measure='auc')@y.values[[1]]
  auc.sym[i] <- performance(prediction(1 - p.sym, DE), measure='auc')@y.values[[1]]
  fdr.hpd[i] <- sum(p.adjust(p.hpd, method="BH") < 0.05 & DE == 0) / sum(p.adjust(p.hpd, method="BH") < 0.05)
  fdr.sym[i] <- sum(p.adjust(p.sym, method="BH") < 0.05 & DE == 0) / sum(p.adjust(p.sym, method="BH") < 0.05)
}
boxplot(auc.hpd, auc.sym)
plot(auc.hpd, auc.sym); lines(c(0.5,1),c(0.5,1))
plot(auc.hpd-auc.sym)
mean(auc.hpd-auc.sym) # 0.0007; -2e07; 7e-5; 0.0005
mean(auc.hpd>auc.sym) # 0.72; 0.5; 0.3; 0.6
boxplot(fdr.hpd, fdr.sym)
plot(abs(fdr.hpd - 0.05), abs(fdr.sym - 0.05)); lines(c(-0.05,1),c(-0.05,1))
mean(fdr.hpd-fdr.sym) # 0.0025; 0.016; 0.006; 0.0009
c(mean(fdr.hpd > 0.05), mean(fdr.sym > 0.05)) # (0.17, 0.17), (0.7, 0.6), (0.6, 0.5), (0.9,0.9)
c(mean(fdr.hpd), mean(fdr.sym)) # (0.035, 0.032); (0.13, 0.11); (0.059, 0.053); (0.182,0.181)
# above for 2k genes, 1k DE; 2k genes, 120 DE; 10k genes, 1k DE; 10k genes, 200 DE



hdi.me <- function(x, credMass=0.95) {
  n <- length(x)
  exclude <- n - floor(n * credMass)   # Number of values to exclude
  low.poss <- x[1:exclude]             # Possible lower limits...
  upp.poss <- x[(n - exclude + 1):n]   # ... and corresponding upper limits
  best <- which.min(upp.poss - low.poss)  # Combination giving the narrowest interval
  result <- c(low.poss[best], upp.poss[best])
  return(result)
}
# Uses code from hdi() from HDInterval, but removed sorting part so can sort
# data once in hpd.pval() instead of sorting for every iteration

hpd.pval <- function(x, m=0) {
  x <- sort.int(x, method='quick')
  y1 <- numeric(9)
  for (t in 1:9) {
    y1[t] <- sign(hdi.me(x, credMass=1-t/10)[1]-m)==sign(hdi.me(x, credMass=1-t/10)[2]+m)
    if (y1[t]==1) {break}
  }
  if (sum(y1)==0) {x1 <- 9} else {x1 <- min(which(y1==1))-1}
  iter <- x1/10
  y2 <- numeric(9)
  for (t in 1:9) {
    y2[t] <- sign(hdi.me(x, credMass=1-(iter + t/100))[1]-m)==sign(hdi.me(x, credMass=1-(iter + t/100))[2]+m)
    if (y2[t]==1) {break}
  }
  if (sum(y2)==0) {x2 <- 9} else {x2 <- min(which(y2==1))-1}
  iter <- iter + x2/100
  y3 <- numeric(9)
  for (t in 1:9) {
    y3[t] <- sign(hdi.me(x, credMass=1-(iter + t/1000))[1]-m)==sign(hdi.me(x, credMass=1-(iter + t/1000))[2]+m)
    if (y3[t]==1) {break}
  }
  if (sum(y3)==0) {x3 <- 9} else {x3 <- min(which(y3==1))-1}
  iter <- iter + x3/1000
  y4 <- numeric(9)
  for (t in 1:9) {
    y4[t] <- sign(hdi.me(x, credMass=1-(iter + t/10000))[1]-m)==sign(hdi.me(x, credMass=1-(iter + t/10000))[2]+m)
    if (y4[t]==1) {break}
  }
  if (sum(y4)==0) {x4 <- 10} else {x4 <- min(which(y4==1))}
  iter <- iter + x4/10000
  return(iter)
}


library(here)
library(HDInterval)
library(compcodeR)
library(ROCR)
library(edgeR)
library(DESeq2)
source(here('scripts','2019-03-27_lognormal_hmm_adaptive_proposals_three_chains_function.R'))

for (i in 1:50) {
  counts <- generateSyntheticData(dataset='DEonly', n.vars=20000, samples.per.cond=5, 
                                  n.diffexp=1000, fraction.upregulated=0.5,
                                  filter.threshold.mediancpm=0.5)
  DE <- counts@variable.annotations$differential.expression
  group <- factor(c(rep(1,5), rep(2,5)))
  nf.TMM <- calcNormFactors(counts@count.matrix)
  norm.TMM <- t(t(counts@count.matrix) / nf.TMM)
  dat.DESeq <- DESeqDataSetFromMatrix(countData=counts@count.matrix, colData=data.frame(group), 
                                      design=~group)
  dat.DESeq <- estimateSizeFactors(dat.DESeq)
  nf.DESeq <- dat.DESeq$sizeFactor
  rm(dat.DESeq)
  norm.DESeq <- t(t(counts@count.matrix) / nf.DESeq)
  lnHM.raw <- ln_hmm_adapt_3_chains(counts=t(counts@count.matrix), groups=c(rep(1,5),rep(2,5)))
  rawpost <- log(as.matrix(lnHM.raw$means1)) - log(as.matrix(lnHM.raw$means2))
  lnHM.TMM <- ln_hmm_adapt_3_chains(counts=t(norm.TMM), groups=c(rep(1,5),rep(2,5)))
  TMMpost <- log(as.matrix(lnHM.TMM$means1)) - log(as.matrix(lnHM.TMM$means2))
  lnHM.DESeq <- ln_hmm_adapt_3_chains(counts=t(norm.DESeq), groups=c(rep(1,5),rep(2,5)))
  DESeqpost <- log(as.matrix(lnHM.DESeq$means1)) - log(as.matrix(lnHM.DESeq$means2))
  saveRDS(rawpost, paste0('rawpost',i,'.rds'))
  saveRDS(TMMpost, paste0('TMMpost',i,'.rds'))
  saveRDS(DESeqpost, paste0('DESeqpost',i,'.rds'))
  saveRDS(DE, paste0('DE',i,'.rds'))
}

for (i in 1:50) {
  rawpost <- readRDS(paste0('rawpost',i,'.rds'))
  TMMpost <- readRDS(paste0('TMMpost',i,'.rds'))
  DESeqpost <- readRDS(paste0('DESeqpost',i,'.rds'))
  DE <- readRDS(paste0('DE',i,'.rds'))
  p.raw.hpd <- apply(rawpost, 2, hpd.pval)
  p.TMM.hpd <- apply(TMMpost, 2, hpd.pval)
  p.DESeq.hpd <- apply(DESeqpost, 2, hpd.pval)
  p.raw.sym <- apply(rawpost, 2, function(x) {2 * min(mean(x < 0), mean(x > 0))})
  p.TMM.sym <- apply(TMMpost, 2, function(x) {2 * min(mean(x < 0), mean(x > 0))})
  p.DESeq.sym <- apply(DESeqpost, 2, function(x) {2 * min(mean(x < 0), mean(x > 0))})
  auc.raw.hpd[i] <- performance(prediction(1 - p.raw.hpd, DE), measure='auc')@y.values[[1]]
  auc.raw.sym[i] <- performance(prediction(1 - p.raw.sym, DE), measure='auc')@y.values[[1]]
  fdr.raw.hpd[i] <- sum(p.adjust(p.raw.hpd, method="BH") < 0.05 & DE == 0) / sum(p.adjust(p.raw.hpd, method="BH")< 0.05)
  fdr.raw.sym[i] <- sum(p.adjust(p.raw.sym, method="BH") < 0.05 & DE == 0) / sum(p.adjust(p.raw.sym, method="BH") < 0.05)
  auc.TMM.hpd[i] <- performance(prediction(1 - p.TMM.hpd, DE), measure='auc')@y.values[[1]]
  auc.TMM.sym[i] <- performance(prediction(1 - p.TMM.sym, DE), measure='auc')@y.values[[1]]
  fdr.TMM.hpd[i] <- sum(p.adjust(p.TMM.hpd, method="BH") < 0.05 & DE == 0) / sum(p.adjust(p.TMM.hpd, method="BH")< 0.05)
  fdr.TMM.sym[i] <- sum(p.adjust(p.TMM.sym, method="BH") < 0.05 & DE == 0) / sum(p.adjust(p.TMM.sym, method="BH")< 0.05)
  auc.DESeq.hpd[i] <- performance(prediction(1 - p.DESeq.hpd, DE), measure='auc')@y.values[[1]]
  auc.DESeq.sym[i] <- performance(prediction(1 - p.DESeq.sym, DE), measure='auc')@y.values[[1]]
  fdr.DESeq.hpd[i] <- sum(p.adjust(p.DESeq.hpd, method="BH") < 0.05 & DE == 0) / sum(p.adjust(p.DESeq.hpd, method="BH")< 0.05)
  fdr.DESeq.sym[i] <- sum(p.adjust(p.DESeq.sym, method="BH") < 0.05 & DE == 0) / sum(p.adjust(p.DESeq.sym, method="BH")< 0.05)
}

boxplot(auc.raw.hpd, auc.raw.sym, auc.TMM.hpd, auc.TMM.sym, auc.DESeq.hpd, auc.DESeq.sym)
c(mean(auc.raw.hpd), mean(auc.raw.sym), mean(auc.TMM.hpd), mean(auc.TMM.sym), mean(auc.DESeq.hpd), mean(auc.DESeq.sym))
# [1] 0.8865605 0.8862184 0.8869993 0.8866074 0.9068692 0.9063386
c(round(100 * ((mean(auc.raw.hpd) - mean(auc.raw.sym)) / mean(auc.raw.sym)), 2), 
  round(100 * ((mean(auc.TMM.hpd) - mean(auc.raw.sym)) / mean(auc.TMM.sym)), 2), 
  round(100 * ((mean(auc.DESeq.hpd) - mean(auc.DESeq.sym)) / mean(auc.DESeq.sym)), 2))
# hpd better in all cases by very small amount: 0.04%, 0.09%, 0.06% for raw, TMM, DESeq
# (difference relative to sym, not actual difference in AUC)
boxplot(fdr.raw.hpd, fdr.raw.sym, fdr.TMM.hpd, fdr.TMM.sym, fdr.DESeq.hpd, fdr.DESeq.sym)
abline(h = seq(0.03, 0.07, 0.02), col='lightgrey')
c(mean(fdr.raw.hpd), mean(fdr.raw.sym), mean(fdr.TMM.hpd), mean(fdr.TMM.sym), mean(fdr.DESeq.hpd), mean(fdr.DESeq.sym))
# [1] 0.09070200 0.08616247 0.09122712 0.08882409 0.06217058 0.05454797
# sym closer to 0.05 on average for all normalisations.
c(range(fdr.DESeq.hpd), range(fdr.DESeq.sym))
# Maximum similar for hpd and sym, but min for sym much lower.
c(mean(fdr.DESeq.hpd > 0.05), mean(fdr.DESeq.sym > 0.05))
# hpd above 0.05 in 38/50 runs, sym in 31/50.
c(mean(abs(fdr.DESeq.hpd - 0.05) < 0.01), mean(abs(fdr.DESeq.sym - 0.05) < 0.01))
# hpd in (0.04, 0.06) in 18/50 runs, sym in 23/50.
c(mean(abs(fdr.DESeq.hpd - 0.05) < 0.02), mean(abs(fdr.DESeq.sym - 0.05) < 0.02))
# hpd in (0.03, 0.07) in 34/50 runs, sym in 38/50.
c(mean(abs(fdr.DESeq.hpd - 0.05) < 0.03), mean(abs(fdr.DESeq.sym - 0.05) < 0.03))
# hpd in (0.02, 0.08) in 41/50 runs, sym in 44/50.

hpd_v_sym_auc_fdr <- list(auc.raw.hpd = auc.raw.hpd, 
                          auc.raw.sym = auc.raw.sym, 
                          auc.TMM.hpd = auc.TMM.hpd, 
                          auc.TMM.sym = auc.TMM.sym, 
                          auc.DESeq.hpd = auc.DESeq.hpd, 
                          auc.DESeq.sym = auc.DESeq.sym, 
                          fdr.raw.hpd = fdr.raw.hpd, 
                          fdr.raw.sym = fdr.raw.sym, 
                          fdr.TMM.hpd = fdr.TMM.hpd, 
                          fdr.TMM.sym = fdr.TMM.sym, 
                          fdr.DESeq.hpd = fdr.DESeq.hpd, 
                          fdr.DESeq.sym = fdr.DESeq.sym)
saveRDS(hpd_v_sym_auc_fdr, file='2019-09-24_hpd_v_sym_auc_fdr')

for (i in 1:50) {
  rawpost <- readRDS(paste0('rawpost',i,'.rds'))
  saveRDS(rawpost, file=paste0('log.diff.mean.DE5.',i,'.nonorm'))
  TMMpost <- readRDS(paste0('TMMpost',i,'.rds'))
  saveRDS(TMMpost, file=paste0('log.diff.mean.DE5.',i,'.TMM'))
  DESeqpost <- readRDS(paste0('DESeqpost',i,'.rds'))
  saveRDS(DESeqpost, file=paste0('log.diff.mean.DE5.',i,'.DESeqnorm'))
  DE <- readRDS(paste0('DE',i,'.rds'))
  saveRDS(DE, file=paste0('DE.DE5.',i,'.rds'))
}

# This was for log-transformed samples. Untransformed samples generally look to give slightly 
# lower AUCs and much better FDRs with hpd. Should see what looks like with sym, although 
# would expect symmetric intervals to not work as well with untransformed samples.
# Therefore probably expect lower AUCs generally for untransformed samples, and particularly 
# lower with sym, and better FDRs for untransformed samples with hpd, but probably worse with 
# sym. So still likely to have a conflict - best AUCs for transformed with hpd, and best FDRs 
# either for transformed with sym or untransformed with hpd.

group <- factor(c(rep(1,5), rep(2,5)))
for (i in 1:50) {
  counts <- generateSyntheticData(dataset='DEonly', n.vars=20000, samples.per.cond=5, 
                                  n.diffexp=1000, fraction.upregulated=0.5,
                                  filter.threshold.mediancpm=0.5)
  DE <- counts@variable.annotations$differential.expression
  dat.DESeq <- DESeqDataSetFromMatrix(countData=counts@count.matrix, colData=data.frame(group), 
                                      design=~group)
  dat.DESeq <- estimateSizeFactors(dat.DESeq)
  nf.DESeq <- dat.DESeq$sizeFactor
  rm(dat.DESeq)
  norm.DESeq <- t(t(counts@count.matrix) / nf.DESeq)
  lnHM.DESeq <- ln_hmm_adapt_3_chains(counts=t(norm.DESeq), groups=c(rep(1,5),rep(2,5)))
  DESeqpost <- as.matrix(lnHM.DESeq$means1) - as.matrix(lnHM.DESeq$means2)
  saveRDS(DESeqpost, file=paste0('diff.mean.DE5.',i,'.DESeqnorm.rds'))
  saveRDS(DE, file=paste0('DE.DE5.',i,'.rds'))
}

auc.DESeq.hpd <- numeric(50); auc.DESeq.sym <- numeric(50)
fdr.DESeq.hpd <- numeric(50); fdr.DESeq.sym <- numeric(50)
for (i in 1:50) {
  DESeqpost <- readRDS(paste0('diff.mean.DE5.',i,'.DESeqnorm.rds'))
  DE <- readRDS(paste0('DE.DE5.',i,'.rds'))
  p.DESeq.hpd <- apply(DESeqpost, 2, hpd.pval)
  p.DESeq.sym <- apply(DESeqpost, 2, function(x) {2 * min(mean(x < 0), mean(x > 0))})
  auc.DESeq.hpd[i] <- performance(prediction(1 - p.DESeq.hpd, DE), measure='auc')@y.values[[1]]
  auc.DESeq.sym[i] <- performance(prediction(1 - p.DESeq.sym, DE), measure='auc')@y.values[[1]]
  fdr.DESeq.hpd[i] <- sum(p.adjust(p.DESeq.hpd, method="BH") < 0.05 & DE == 0) / sum(p.adjust(p.DESeq.hpd, method="BH")< 0.05)
  fdr.DESeq.sym[i] <- sum(p.adjust(p.DESeq.sym, method="BH") < 0.05 & DE == 0) / sum(p.adjust(p.DESeq.sym, method="BH")< 0.05)
}
boxplot(auc.DESeq.hpd, auc.DESeq.sym)
c(mean(auc.DESeq.hpd), mean(auc.DESeq.sym))
# sym higher, as expected for untransformed samples
round(100 * ((mean(auc.DESeq.sym) - mean(auc.DESeq.hpd)) / mean(auc.DESeq.hpd)), 2)
# sym better by 0.28% of hpd value (i.e. difference relative to hpd, not actual difference in AUC)
boxplot(fdr.DESeq.hpd, fdr.DESeq.sym)
abline(h = seq(0.03, 0.07, 0.02), col='lightgrey')
c(mean(fdr.DESeq.hpd), mean(fdr.DESeq.sym))
# [1] 0.04948032 0.05596006
# hpd very slightly closer to 0.05 on average, and below where sym is above
c(range(fdr.DESeq.hpd), range(fdr.DESeq.sym))
# max further from 0.05 for sym, min further from 0.05 for hpd
c(mean(fdr.DESeq.hpd > 0.05), mean(fdr.DESeq.sym > 0.05))
# hpd above 0.05 in 27/50 runs, sym in 30/50.
c(mean(abs(fdr.DESeq.hpd - 0.05) < 0.01), mean(abs(fdr.DESeq.sym - 0.05) < 0.01))
# hpd in (0.04, 0.06) in 22/50 runs, sym in 17/50.
c(mean(abs(fdr.DESeq.hpd - 0.05) < 0.02), mean(abs(fdr.DESeq.sym - 0.05) < 0.02))
# hpd in (0.03, 0.07) in 36/50 runs, sym in 32/50.
c(mean(abs(fdr.DESeq.hpd - 0.05) < 0.03), mean(abs(fdr.DESeq.sym - 0.05) < 0.03))
# hpd in (0.02, 0.08) in 47/50 runs, sym in 45/50.


## For DESeq norm, directly compare hpd and sym intervals for untransformed and log posterior samples ####
# These are for different data (same characteristics and same for hpd and sym within each posterior type, but 
# untransformed and log transformed are from different datasets; don't expect to affect qualitative results over 
# 50 simulations)
auc.hpd <- numeric(50); auc.sym <- numeric(50); auc.hpd.log <- numeric(50); auc.sym.log <- numeric(50)
fdr.hpd <- numeric(50); fdr.sym <- numeric(50); fdr.hpd.log <- numeric(50); fdr.sym.log <- numeric(50)
for (i in 1:50) {
  DESeqpost <- readRDS(here('Results/HPD v symmetric posterior intervals Sept 2019',
                            paste0('diff.mean.DE5.',i,'.DESeqnorm.rds')))
  DE <- readRDS(here('Results/HPD v symmetric posterior intervals Sept 2019',
                     paste0('DE.DE5.',i,'.rds')))
  p.hpd <- apply(DESeqpost, 2, hpd.pval)
  p.sym <- apply(DESeqpost, 2, function(x) {2 * min(mean(x < 0), mean(x > 0))})
  auc.hpd[i] <- performance(prediction(1 - p.hpd, DE), measure='auc')@y.values[[1]]
  auc.sym[i] <- performance(prediction(1 - p.sym, DE), measure='auc')@y.values[[1]]
  fdr.hpd[i] <- sum(p.adjust(p.hpd, method="BH") < 0.05 & DE == 0) / sum(p.adjust(p.hpd, method="BH")< 0.05)
  fdr.sym[i] <- sum(p.adjust(p.sym, method="BH") < 0.05 & DE == 0) / sum(p.adjust(p.sym, method="BH")< 0.05)
  DESeqpost.log <- readRDS(here('Results/HPD v symmetric posterior intervals Sept 2019', 
                                paste0('log.diff.mean.DE5.',i,'.DESeqnorm.rds')))
  DE.log <- readRDS(here('Results/HPD v symmetric posterior intervals Sept 2019', 
                         paste0('DE.DE5.log.',i,'.rds')))
  p.hpd.log <- apply(DESeqpost.log, 2, hpd.pval)
  p.sym.log <- apply(DESeqpost.log, 2, function(x) {2 * min(mean(x < 0), mean(x > 0))})
  auc.hpd.log[i] <- performance(prediction(1 - p.hpd.log, DE.log), measure='auc')@y.values[[1]]
  auc.sym.log[i] <- performance(prediction(1 - p.sym.log, DE.log), measure='auc')@y.values[[1]]
  fdr.hpd.log[i] <- sum(p.adjust(p.hpd.log, method="BH") < 0.05 & DE.log == 0) / sum(p.adjust(p.hpd.log, method="BH")< 0.05)
  fdr.sym.log[i] <- sum(p.adjust(p.sym.log, method="BH") < 0.05 & DE.log == 0) / sum(p.adjust(p.sym.log, method="BH")< 0.05)
}
boxplot(auc.hpd, auc.sym, auc.hpd.log, auc.sym.log)
c(mean(auc.hpd), mean(auc.sym), mean(auc.hpd.log), mean(auc.sym.log))
# [1] 0.9037177 0.9062207 0.9068692 0.9063386
c(median(auc.hpd), median(auc.sym), median(auc.hpd.log), median(auc.sym.log))
# [1] 0.9036432 0.9061415 0.9063813 0.9064621
# hpd untransformed clearly worst, surprisingly. Others similar but smaller variance for sym untransformed.
# Highest mean for hpd log by a very small margin, highest median for sym log by an even smaller margin.
c(round(100 * ((mean(auc.hpd.log) - mean(auc.hpd)) / mean(auc.hpd)), 2), 
  round(100 * ((mean(auc.sym.log) - mean(auc.hpd)) / mean(auc.hpd)), 2), 
  round(100 * ((mean(auc.sym) - mean(auc.hpd)) / mean(auc.hpd)), 2))
# Differences very small - 0.35%, 0.29%, 0.28% better than worst.
boxplot(fdr.hpd, fdr.sym, fdr.hpd.log, fdr.sym.log)
abline(h = seq(0.03, 0.07, 0.02), col='lightgrey')
c(mean(fdr.hpd), mean(fdr.sym), mean(fdr.hpd.log), mean(fdr.sym.log))
# [1] 0.04948032 0.05596006 0.06217058 0.05454797
c(abs(0.05 - mean(fdr.hpd)), abs(0.05 - mean(fdr.sym)), abs(0.05 - mean(fdr.hpd.log)), abs(0.05 - mean(fdr.sym.log)))
# [1] 0.0005196766 0.0059600586 0.0121705770 0.0045479727
c(median(fdr.hpd), median(fdr.sym), median(fdr.hpd.log), median(fdr.sym.log))
# [1] 0.05147138 0.05495169 0.06068977 0.05399015
# Only hpd untransformed has mean below 0.05, which has clearly the worst AUC.
# Also closest to 0.05, followed by sym log.
# Also closest median to 0.05 although slighly above, followed by sym log.
c(range(fdr.hpd), range(fdr.sym), range(fdr.hpd.log), range(fdr.sym.log))
# Both untransformed have extremes closer to 0.05 than both transformed.
c(mean(fdr.hpd > 0.05), mean(fdr.sym > 0.05), mean(fdr.hpd.log > 0.05), mean(fdr.sym.log > 0.05))
# [1] 0.54 0.60 0.76 0.62
# hpd untransformed best, hpd log worst
c(mean(abs(fdr.hpd - 0.05) < 0.01), mean(abs(fdr.sym - 0.05) < 0.01), 
  mean(abs(fdr.hpd.log - 0.05) < 0.01), mean(abs(fdr.sym.log - 0.05) < 0.01))
# [1] 0.44 0.34 0.36 0.46
# sym log best, followed by hpd untransformed
c(mean(abs(fdr.hpd - 0.05) < 0.02), mean(abs(fdr.sym - 0.05) < 0.02), 
  mean(abs(fdr.hpd.log - 0.05) < 0.02), mean(abs(fdr.sym.log - 0.05) < 0.02))
# sym log best, followed by hpd untransformed
c(mean(abs(fdr.hpd - 0.05) < 0.03), mean(abs(fdr.sym - 0.05) < 0.03), 
  mean(abs(fdr.hpd.log - 0.05) < 0.03), mean(abs(fdr.sym.log - 0.05) < 0.03))
# hpd untransformed best, followed by sym untransformed

# Partial AUCs may give a better idea - not really interested in very low ranked genes.
# What should be max FPR? The one paper I've seen use partial AUCs (Chu 2015, deGPS) uses 
# 0.05, but even with such a low FPR, the FDR is likely to be very high, so it seems that 
# a lower cutoff would be a better reflection of how a method is going to be used in 
# practice. A very brief exploration shows that a threshold p-value that gives FPR around 
# 0.05 gives FDR around 0.5-0.6, and the thresholds are around 0.06-0.09, which would give 
# far too many positives in practice. Without knowing the actual number of positives, it's 
# difficult to identify a good threshold, but it surely should be less than 0.05 since 
# is the equivalent of not using any FDR adjustment. However, an FPR cutoff of 0.05 is 
# definitely a better reflection than full AUC, so maybe I should start with that anyway.
pauc.hpd <- numeric(50); pauc.sym <- numeric(50); 
pauc.hpd.log <- numeric(50); pauc.sym.log <- numeric(50)
for (i in 1:50) {
  DESeqpost <- readRDS(here('Results/HPD v symmetric posterior intervals Sept 2019',
                            paste0('diff.mean.DE5.',i,'.DESeqnorm.rds')))
  DE <- readRDS(here('Results/HPD v symmetric posterior intervals Sept 2019',
                     paste0('DE.DE5.',i,'.rds')))
  p.hpd <- apply(DESeqpost, 2, hpd.pval)
  p.sym <- apply(DESeqpost, 2, function(x) {2 * min(mean(x < 0), mean(x > 0))})
  pauc.hpd[i] <- performance(prediction(1 - p.hpd, DE), measure='auc', fpr.stop=0.05)@y.values[[1]]
  pauc.sym[i] <- performance(prediction(1 - p.sym, DE), measure='auc', fpr.stop=0.05)@y.values[[1]]
  DESeqpost.log <- readRDS(here('Results/HPD v symmetric posterior intervals Sept 2019', 
                                paste0('log.diff.mean.DE5.',i,'.DESeqnorm.rds')))
  DE.log <- readRDS(here('Results/HPD v symmetric posterior intervals Sept 2019', 
                         paste0('DE.DE5.log.',i,'.rds')))
  p.hpd.log <- apply(DESeqpost.log, 2, hpd.pval)
  p.sym.log <- apply(DESeqpost.log, 2, function(x) {2 * min(mean(x < 0), mean(x > 0))})
  pauc.hpd.log[i] <- performance(prediction(1 - p.hpd.log, DE.log), measure='auc', fpr.stop=0.05)@y.values[[1]]
  pauc.sym.log[i] <- performance(prediction(1 - p.sym.log, DE.log), measure='auc', fpr.stop=0.05)@y.values[[1]]
}
boxplot(pauc.hpd, pauc.sym, pauc.hpd.log, pauc.sym.log)
c(mean(pauc.hpd), mean(pauc.sym), mean(pauc.hpd.log), mean(pauc.sym.log))
# [1] 0.02951159 0.02910296 0.02938042 0.02922486
c(median(pauc.hpd), median(pauc.sym), median(pauc.hpd.log), median(pauc.sym.log))
# [1] 0.02951261 0.02909361 0.02946230 0.02933079
# hpd > hpd log > sym log > sym
# Untransformed hpd worst for overall AUC but best for AUC up to FPR 0.05, as well as best FDR.
c(round(100 * ((mean(pauc.hpd) - mean(pauc.sym)) / mean(pauc.sym)), 2), 
  round(100 * ((mean(pauc.hpd.log) - mean(pauc.sym)) / mean(pauc.sym)), 2), 
  round(100 * ((mean(pauc.sym.log) - mean(pauc.sym)) / mean(pauc.sym)), 2))
# Differences bigger than for overall AUC - 1.4%, 0.95%, 0.42%.


## This isn't a really comprehensive evaluation, but it looks like enough to say that:
## 1. There aren't big differences between symmetric and HPD intervals or with and without transforming.
## 2. Untransformed HPD intervals are probably best overall - best FDR performance and pAUC.

