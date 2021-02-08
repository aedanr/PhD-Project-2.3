library(here)
library(edgeR)
library(DESeq2)
library(limma)
library(DSS)
library(baySeq)
library(MDSeq)
library(ROCR)
library(caret)
source(here('scripts','2019-04-03_exponential_hmm_adaptive_proposals_three_chains_function.R'))
source(here('scripts','2019-03-27_lognormal_hmm_adaptive_proposals_three_chains_function.R'))
source(here('scripts','2019-06-26_hpd_tail_prob_function.R'))

i <- 5
samples.per.cond <- 10
filename <- paste0('DE', samples.per.cond, '.', i)
data <- readRDS(here('Simulated data', paste0(filename,'.rds')))

counts <- data@count.matrix
DE <- data@variable.annotations$differential.expression
group <- factor(c(rep(1,samples.per.cond), rep(2,samples.per.cond)))
design <- model.matrix(~group)

dge <- DGEList(counts=counts, group=group)
dds <- DESeqDataSetFromMatrix(countData=counts, colData=data.frame(group), design=~group)
dds.rle <- estimateSizeFactors(dds)

nf.tmm <- calcNormFactors(counts, method="TMM")
dge.tmm <- DGEList(counts=counts, norm.factors=nf.tmm, group=group)
els.tmm <- nf.tmm * colSums(counts)
norm.counts.nf.tmm <- t(t(counts) * exp(mean(log(els.tmm))) / els.tmm)

nf.rle <- calcNormFactors(counts, method="RLE")
dge.rle <- DGEList(counts=counts, norm.factors=nf.rle, group=group)
els.rle <- nf.rle * colSums(counts)
norm.counts.nf.rle <- t(t(counts) * exp(mean(log(els.rle))) / els.rle)

sf.rle <- sizeFactors(dds.rle)
sf.rle.adj <- sf.rle / exp(mean(log(sf.rle)))
norm.counts.sf.rle <- t(t(counts) / sf.rle)
norm.counts.sf.rle.adj <- t(t(counts) / sf.rle.adj)

sf.tmm <- els.tmm / exp(mean(log(els.tmm)))
norm.counts.sf.tmm <- t(t(counts) / sf.tmm)

mean(norm.counts.nf.rle - norm.counts.sf.rle.adj < 1e-6) # 1
mean(norm.counts.nf.tmm - norm.counts.sf.tmm < 1e-6) # 1


# To get from size factors to normalised counts, divide counts by size factors.
# To get from normalisation factors to effective library sizes, multiply by library sizes.
# To get from effective library sizes to size factors, divide by geometric mean of effective library sizes.
# To get from effective library sizes to normalised counts, divide counts by (effective library sizes divided by geometric mean of effective library sizes).
# To get from effective library sizes to normalised counts, multiply counts by geometric mean of effective library sizes and divide by effective library sizes.
# To make size factors match between edgeR and DESeq2, divide factors from DESeq2 by their geometric mean.


dge.default.tmm <- DGEList(counts=counts, group=group, norm.factors=nf.tmm)
dge.default.rle <- DGEList(counts=counts, group=group, norm.factors=nf.rle)
dge.tmm <- DGEList(counts=norm.counts.nf.tmm, group=group)
dge.rle <- DGEList(counts=norm.counts.sf.rle, group=group)
dge.rle.adj <- DGEList(counts=norm.counts.sf.rle.adj, group=group)

dds.default.tmm <- DESeqDataSetFromMatrix(countData=counts, colData=data.frame(group), design=~group)
sizeFactors(dds.default.tmm) <- sf.tmm
dds.default.rle <- DESeqDataSetFromMatrix(countData=counts, colData=data.frame(group), design=~group)
sizeFactors(dds.default.rle) <- sf.rle
dds.tmm <- DESeqDataSetFromMatrix(countData=round(norm.counts.nf.tmm), colData=data.frame(group), design=~group)
sizeFactors(dds.tmm) <- rep(1, ncol(counts))
dds.rle <- DESeqDataSetFromMatrix(countData=round(norm.counts.sf.rle), colData=data.frame(group), design=~group)
sizeFactors(dds.rle) <- rep(1, ncol(counts))
dds.rle.adj <- DESeqDataSetFromMatrix(countData=round(norm.counts.sf.rle.adj), colData=data.frame(group), design=~group)
sizeFactors(dds.rle.adj) <- rep(1, ncol(counts))
# Might not be really comparable as DESeq2 needs integer input

for (i in c('default.tmm', 'default.rle', 'tmm', 'rle', 'rle.adj')) {
  assign(paste0('dge.', i), estimateDisp(get(paste0('dge.', i)), design))
  assign(paste0('fit.', i, '.edgeR'), glmFit(get(paste0('dge.', i)), design))
  assign(paste0('test.', i, '.edgeR'), glmLRT(get(paste0('fit.', i, '.edgeR'))))
  assign(paste0('p.', i, '.edgeR'), get(paste0('test.', i, '.edgeR'))$table$PValue)
  assign(paste0('q.', i, '.edgeR'), topTags(get(paste0('test.', i, '.edgeR')),
                                            n=nrow(counts), sort='none')$table$FDR)
  assign(paste0('dat.', i, '.DESeq2'), DESeq(get(paste0('dds.', i)), minReplicatesForReplace=Inf))
  assign(paste0('res.', i, '.DESeq2'), results(get(paste0('dat.', i, '.DESeq2')), cooksCutoff=F))
  assign(paste0('p.', i, '.DESeq2'), get(paste0('res.', i, '.DESeq2'))$pvalue)
  assign(paste0('q.', i, '.DESeq2'), get(paste0('res.', i, '.DESeq2'))$padj)
}

for (i in c('default.tmm', 'default.rle', 'tmm', 'rle', 'rle.adj')) {
  for (j in c('.edgeR', '.DESeq2')) {
    assign(paste0('fdr.', i, j), 1 - precision(factor(get(paste0('q.', i, j)) < 0.05, levels=c("TRUE", "FALSE")), 
                                               factor(DE == 1, levels=c("TRUE", "FALSE"))))
    assign(paste0('tpr.', i, j), sensitivity(factor(get(paste0('q.', i, j)) < 0.05, levels=c("TRUE", "FALSE")), 
                                             factor(DE == 1, levels=c("TRUE", "FALSE"))))
    assign(paste0('pred.', i, j), prediction(1 - get(paste0('p.', i, j)), DE))
    assign(paste0('auc.', i, j), performance(get(paste0('pred.', i, j)), measure='auc')@y.values[[1]])
    assign(paste0('pauc.', i, j), performance(get(paste0('pred.', i, j)), 
                                              measure='auc', fpr.stop=0.02)@y.values[[1]])
  }
}
fdr.list <- c(fdr.default.tmm.edgeR, fdr.default.rle.edgeR, fdr.tmm.edgeR, fdr.rle.edgeR, fdr.rle.adj.edgeR,
              fdr.default.tmm.DESeq2, fdr.default.rle.DESeq2, fdr.tmm.DESeq2, fdr.rle.DESeq2, fdr.rle.adj.DESeq2)
tpr.list <- c(tpr.default.tmm.edgeR, tpr.default.rle.edgeR, tpr.tmm.edgeR, tpr.rle.edgeR, tpr.rle.adj.edgeR,
              tpr.default.tmm.DESeq2, tpr.default.rle.DESeq2, tpr.tmm.DESeq2, tpr.rle.DESeq2, tpr.rle.adj.DESeq2)
auc.list <- c(auc.default.tmm.edgeR, auc.default.rle.edgeR, auc.tmm.edgeR, auc.rle.edgeR, auc.rle.adj.edgeR,
              auc.default.tmm.DESeq2, auc.default.rle.DESeq2, auc.tmm.DESeq2, auc.rle.DESeq2, auc.rle.adj.DESeq2)
results <- data.frame(fdr = fdr.list, tpr = tpr.list, auc = auc.list)
rownames(results) <- c('default.tmm.edgeR', 'default.rle.edgeR', 'tmm.edgeR', 'rle.edgeR', 'rle.adj.edgeR',
                       'default.tmm.DESeq2', 'default.rle.DESeq2', 'tmm.DESeq2', 'rle.DESeq2', 'rle.adj.DESeq2')
#                            fdr       tpr       auc
# default.tmm.edgeR   0.04099379 0.9566295 0.9931840
# default.rle.edgeR   0.04337051 0.9566295 0.9931537
# tmm.edgeR           0.08033573 0.9504337 0.9924105
# rle.edgeR           0.08033573 0.9504337 0.9924205
# rle.adj.edgeR       0.08033573 0.9504337 0.9924107
# default.tmm.DESeq2  0.04926108 0.9566295 0.9930765
# default.rle.DESeq2  0.04926108 0.9566295 0.9930486
# tmm.DESeq2          0.04938272 0.9541512 0.9930723
# rle.DESeq2          0.04932182 0.9553903 0.9930227
# rle.adj.DESeq2      0.04938272 0.9541512 0.9930907

results.DE10.rle <- readRDS(here('Results/DE compcodeR data results DESeq norm Aug-Sept 2019', 
                                 'DE.results.DE10.DESeqnorm.rds'))
results.DE10.tmm <- readRDS(here('Results/DE compcodeR data results July-Aug 2019', 
                                 'DE.results.DE10.TMM.rds'))
rbind(c(results.DE10.tmm$fdr$fdr.lr.edgeR[5], results.DE10.tmm$tpr$fdr.lr.edgeR[5], results.DE10.tmm$auc$lr.edgeR[5]), 
      c(results.DE10.rle$fdr$fdr.lr.edgeR[5], results.DE10.rle$tpr$fdr.lr.edgeR[5], results.DE10.rle$auc$lr.edgeR[5]), 
      c(results.DE10.tmm$fdr$fdr.if.DESeq[5], results.DE10.tmm$tpr$fdr.if.DESeq[5], results.DE10.tmm$auc$if.DESeq[5]), 
      c(results.DE10.rle$fdr$fdr.if.DESeq[5], results.DE10.rle$tpr$fdr.if.DESeq[5], results.DE10.rle$auc$if.DESeq[5]))
#            [,1]      [,2]      [,3]
# [1,] 0.04099379 0.9566295 0.9931840 matches default.tmm.edgeR
# [2,] 0.08564536 0.9392813 0.9898013 doesn't match anything
# [3,] 0.06341463 0.9516729 0.9932447 doesn't match anything
# [4,] 0.04926108 0.9566295 0.9930486 matches default.rle.DESeq2

# Results with native normalisation methods are correct, so final results I've used for edgeR and DESeq2 will be correct 
# since I chose those normalisation methods.
# Expect that results for voom will be correct because I've used DGEList objects.
# Expect DSS results to be correct as well because I've used size factors from DESeq2 and DSS treats its normalisation
# factors the same way DESeq2 uses size factors.
# Not sure about baySeq - set libsizes using nf, which, if that means relative library size factors, would be wrong since I 
# used factors from edgeR.
# MDSeq results should be correct because used my normalised matrices created by dividing counts by DESeq2 size factors.
# HM results should also be correct for the same reason.
# But for all cases, if the results I've used are correct, the results for the other normalisation method will be wrong.
# 
# A separate issue seems to be that using pre-normalised counts and specifying normalisation factors as 1 doesn't work for 
# edgeR. 
# It's not clear if it's valid for DESeq2 because results are similar but don't expect them to be identical since I had to 
# round normalised counts for it to run.
# Not really an issue though since I can specify normalisation or size factors as appropriate. Can also do this for DSS and 
# baySeq. Not sure about MDSeq. For HM I know that specifying pre-normalised counts is fine.

dat.baySeq <- new('countData', data=counts, replicates=group, groups=list(NDE=rep(1,length(group)), DE=as.numeric(group)))
getLibsizes(dat.baySeq, estimationType='edgeR')
# getLibsizes() with estimationType='edgeR' returns same as effective library sizes from edgeR, i.e. nf*libsizes.
# Looking at code, it calls calcNormFactors and multiplies norm.factors by library sizes.
# I've specified nf, which might mean that I've done it really really wrongly, but that seems unlikely because 
# the results were ok. Will need to test, but probably need to run all again. At least expect to need to run RLE 
# again because it seems that, if specifying nf works, it probably treats them the same way as edgeR.

cnf <- calcNormFactors(counts, method="TMM") # identical to nf.tmm
libsize <- colSums(counts)
rellibsize <- libsize / exp(mean(log(libsize)))
nf <- cnf * rellibsize # identical to sf.tmm to 6 sf
# To go from norm factors to size factors, just multiply by library size divided by geometric mean of library size.
# This is what is used by MDSeq to include normalisation factors as offsets. To use normalised counts, it uses its 
# own function normalize.counts():
norm.counts.MDSeq.tmm <- normalize.counts(counts, group=group, method="TMM")
# Integers and doesn't look at all like any of the other normalised counts.
# Looking at code, calls calcNormFactors, divides counts by nf and rounds. So, apart from rounding, looks to be same 
# approach as I used initially, but which is wrong.
norm.counts.MDSeq.rle <- normalize.counts(counts, group=group, method="RLE")
# It seems that MDSeq does normalisation correctly if using normalisation factors as offsets, but wrongly if using 
# normalised counts. But it could also be that I've misunderstood how normalisation factors should be used to 
# normalise counts (although seems unlikely since my normalised counts agree with DESeq2's normalised counts). It 
# could also possibly be that MDSeq uses the counts in a different way such that their normalised counts end up 
# giving the same results as the offsets.
contrasts <- get.model.matrix(group)
fit.MDSeq.norm.counts.MDSeq <- MDSeq(norm.counts.MDSeq.tmm, contrast=contrasts, mc.cores=detectCores()-1)
res.MDSeq.norm.counts.MDSeq <- extract.ZIMD(fit.MDSeq.norm.counts.MDSeq, compare=list(A="1",B="2"))
p.mean.MDSeq.norm.counts.MDSeq <- res.MDSeq.norm.counts.MDSeq$Pvalue.mean
q.mean.MDSeq.norm.counts.MDSeq <- res.MDSeq.norm.counts.MDSeq$FDR.mean
fit.MDSeq.norm.counts.nf.tmm <- MDSeq(norm.counts.nf.tmm, contrast=contrasts, mc.cores=detectCores()-1)
res.MDSeq.norm.counts.nf.tmm <- extract.ZIMD(fit.MDSeq.norm.counts.nf.tmm, compare=list(A="1",B="2"))
p.mean.MDSeq.norm.counts.nf.tmm <- res.MDSeq.norm.counts.nf.tmm$Pvalue.mean
q.mean.MDSeq.norm.counts.nf.tmm <- res.MDSeq.norm.counts.nf.tmm$FDR.mean
fit.MDSeq.offsets <- MDSeq(counts, contrast=contrasts, mc.cores=detectCores()-1, offsets=nf)
res.MDSeq.offsets <- extract.ZIMD(fit.MDSeq.offsets, compare=list(A="1",B="2"))
p.mean.MDSeq.offsets <- res.MDSeq.offsets$Pvalue.mean
q.mean.MDSeq.offsets <- res.MDSeq.offsets$FDR.mean
# p-values completely different between offsets and MDSeq's normalised counts, and highly correlated between 
# offsets and my normalised counts, so this seems to confirm that the way I'm doing it now is right and 
# MDSeq uses wrongly normalised counts.
for (j in c('norm.counts.MDSeq', 'norm.counts.nf.tmm', 'offsets')) {
  assign(paste0('fdr.', j), 1 - precision(factor(get(paste0('q.mean.MDSeq.', j)) < 0.05, levels=c("TRUE", "FALSE")), 
                                          factor(DE == 1, levels=c("TRUE", "FALSE"))))
  assign(paste0('tpr.', j), sensitivity(factor(get(paste0('q.mean.MDSeq.', j)) < 0.05, levels=c("TRUE", "FALSE")), 
                                        factor(DE == 1, levels=c("TRUE", "FALSE"))))
  assign(paste0('pred.', j), prediction(1 - get(paste0('p.mean.MDSeq.', j)), DE))
  assign(paste0('auc.', j), performance(get(paste0('pred.', j)), measure='auc')@y.values[[1]])
  assign(paste0('pauc.', j), performance(get(paste0('pred.', j)), 
                                         measure='auc', fpr.stop=0.02)@y.values[[1]])
}
rbind(c(fdr.norm.counts.MDSeq, tpr.norm.counts.MDSeq, auc.norm.counts.MDSeq), 
      c(fdr.norm.counts.nf.tmm, tpr.norm.counts.nf.tmm, auc.norm.counts.nf.tmm), 
      c(fdr.offsets, tpr.offsets, auc.offsets))
# Results for my normalised counts computed from TMM norm factors close to but still different from results 
# using TMM norm factors as offsets, and both quite different from - and better than - results using MDSeq's 
# normalised counts.

i <- 20
filename <- paste0("blood_", i, "_DE")
data <- readRDS(here("recount data/GTEx", paste0(filename, ".rds")))
counts <- data$counts
DE <- data$DEindex
group <- factor(c(rep(1,i/2), rep(2,i/2)))

i <- 5
samples.per.cond <- 10
filename <- paste0('DE', samples.per.cond, '.', i)
data <- readRDS(here('Simulated data', paste0(filename,'.rds')))
counts <- data@count.matrix
DE <- data@variable.annotations$differential.expression
group <- factor(c(rep(1,samples.per.cond), rep(2,samples.per.cond)))

contrasts <- get.model.matrix(group)
norm.counts.MDSeq.tmm <- normalize.counts(counts, group=group, method="TMM")
norm.counts.MDSeq.rle <- normalize.counts(counts, group=group, method="RLE")
nf.tmm <- calcNormFactors(counts, method="TMM")
els.tmm <- nf.tmm * colSums(counts)
sf.tmm <- els.tmm / exp(mean(log(els.tmm)))
norm.counts.nf.tmm <- t(t(counts) * exp(mean(log(els.tmm))) / els.tmm)
nf.rle <- calcNormFactors(counts, method="RLE")
els.rle <- nf.rle * colSums(counts)
sf.rle <- els.rle / exp(mean(log(els.rle)))
norm.counts.nf.rle <- t(t(counts) * exp(mean(log(els.rle))) / els.rle)

fit.MDSeq.norm.counts.MDSeq.tmm <- MDSeq(norm.counts.MDSeq.tmm, contrast=contrasts, mc.cores=3)
res.MDSeq.norm.counts.MDSeq.tmm <- extract.ZIMD(fit.MDSeq.norm.counts.MDSeq.tmm, compare=list(A="1",B="2"))
p.MDSeq.norm.counts.MDSeq.tmm <- res.MDSeq.norm.counts.MDSeq.tmm$Pvalue.mean
q.MDSeq.norm.counts.MDSeq.tmm <- res.MDSeq.norm.counts.MDSeq.tmm$FDR.mean
fit.MDSeq.norm.counts.MDSeq.rle <- MDSeq(norm.counts.MDSeq.rle, contrast=contrasts, mc.cores=3)
res.MDSeq.norm.counts.MDSeq.rle <- extract.ZIMD(fit.MDSeq.norm.counts.MDSeq.rle, compare=list(A="1",B="2"))
p.MDSeq.norm.counts.MDSeq.rle <- res.MDSeq.norm.counts.MDSeq.rle$Pvalue.mean
q.MDSeq.norm.counts.MDSeq.rle <- res.MDSeq.norm.counts.MDSeq.rle$FDR.mean
fit.MDSeq.norm.counts.nf.tmm <- MDSeq(norm.counts.nf.tmm, contrast=contrasts, mc.cores=3)
res.MDSeq.norm.counts.nf.tmm <- extract.ZIMD(fit.MDSeq.norm.counts.nf.tmm, compare=list(A="1",B="2"))
p.MDSeq.norm.counts.nf.tmm <- res.MDSeq.norm.counts.nf.tmm$Pvalue.mean
q.MDSeq.norm.counts.nf.tmm <- res.MDSeq.norm.counts.nf.tmm$FDR.mean
fit.MDSeq.norm.counts.nf.rle <- MDSeq(norm.counts.nf.rle, contrast=contrasts, mc.cores=3)
res.MDSeq.norm.counts.nf.rle <- extract.ZIMD(fit.MDSeq.norm.counts.nf.rle, compare=list(A="1",B="2"))
p.MDSeq.norm.counts.nf.rle <- res.MDSeq.norm.counts.nf.rle$Pvalue.mean
q.MDSeq.norm.counts.nf.rle <- res.MDSeq.norm.counts.nf.rle$FDR.mean
fit.MDSeq.offsets.tmm <- MDSeq(counts, contrast=contrasts, mc.cores=3, offsets=sf.tmm)
res.MDSeq.offsets.tmm <- extract.ZIMD(fit.MDSeq.offsets.tmm, compare=list(A="1",B="2"))
p.MDSeq.offsets.tmm <- res.MDSeq.offsets.tmm$Pvalue.mean
q.MDSeq.offsets.tmm <- res.MDSeq.offsets.tmm$FDR.mean
fit.MDSeq.offsets.rle <- MDSeq(counts, contrast=contrasts, mc.cores=3, offsets=sf.rle)
res.MDSeq.offsets.rle <- extract.ZIMD(fit.MDSeq.offsets.rle, compare=list(A="1",B="2"))
p.MDSeq.offsets.rle <- res.MDSeq.offsets.rle$Pvalue.mean
q.MDSeq.offsets.rle <- res.MDSeq.offsets.rle$FDR.mean

for (j in c('norm.counts.MDSeq', 'norm.counts.nf', 'offsets')) {
  for (i in c('.tmm', '.rle')) {
    assign(paste0('fdr.', j, i), 1 - precision(factor(get(paste0('q.MDSeq.', j, i)) < 0.05, levels=c("TRUE", "FALSE")), 
                                               factor(DE == 1, levels=c("TRUE", "FALSE"))))
    assign(paste0('tpr.', j, i), sensitivity(factor(get(paste0('q.MDSeq.', j, i)) < 0.05, levels=c("TRUE", "FALSE")), 
                                             factor(DE == 1, levels=c("TRUE", "FALSE"))))
    assign(paste0('pred.', j, i), prediction(1 - get(paste0('p.MDSeq.', j, i)), DE))
    assign(paste0('auc.', j, i), performance(get(paste0('pred.', j, i)), measure='auc')@y.values[[1]])
    assign(paste0('pauc.', j, i), performance(get(paste0('pred.', j, i)), 
                                              measure='auc', fpr.stop=0.02)@y.values[[1]])
  }
}
rbind(c(fdr.norm.counts.MDSeq.tmm, tpr.norm.counts.MDSeq.tmm, auc.norm.counts.MDSeq.tmm), 
      c(fdr.norm.counts.nf.tmm, tpr.norm.counts.nf.tmm, auc.norm.counts.nf.tmm), 
      c(fdr.offsets.tmm, tpr.offsets.tmm, auc.offsets.tmm), 
      c(fdr.norm.counts.MDSeq.rle, tpr.norm.counts.MDSeq.rle, auc.norm.counts.MDSeq.rle), 
      c(fdr.norm.counts.nf.rle, tpr.norm.counts.nf.rle, auc.norm.counts.nf.rle), 
      c(fdr.offsets.rle, tpr.offsets.rle, auc.offsets.rle))


