# put results together, i.e. from ShrinkBayes and others (can't apply yet)
# how to organise results to get metrics of interest in most useful format
# 50 values for each metric for each method for each dataset
# but first need to get FDR values that aren't included in each method's own analysis

library(here)
library(compcodeR)
library(ROCR)
library(caret)
library(qvalue)
source(here('scripts','2019-05-03_bfdr_function.R'))
res <- readRDS(here('Results/DEDD compcodeR data results July 2019','results.DEDD10.1.rds'))
names(res)

# Data:
# data DE DD lfcm lfcm lfcd lfcd

# edgeR:
#            p.ql.lfc1.edgeR p.ql.lfc2.edgeR p.lr.edgeR p.lr.lfc1.edgeR p.lr.lfc2.edgeR p.et.edgeR
# q.ql.edgeR q.ql.lfc1.edgeR q.ql.lfc2.edgeR q.lr.edgeR q.lr.lfc1.edgeR q.lr.lfc2.edgeR q.et.edgeR
# edgeR user guide says it uses BH.
# mean(res$q.lr.edgeR == p.adjust(res$p.lr.edgeR, method='BH')) is 1, so confirms.
### Need to re-run DEDD, DE to save p.ql.edgeR ###

# DESeq
# p.noif.DESeq p.if.DESeq p.lfc1.DESeq p.lfc2.DESeq
# q.noif.DESeq q.if.DESeq q.lfc1.DESeq q.lfc2.DESeq
# mean(res$q.noif.DESeq == p.adjust(res$p.noif.DESeq, method='BH')) is also 1, and so is 
# mean(res$q.if.DESeq == p.adjust(res$p.if.DESeq, method='BH')).

# voom
# p.voom p.lfc1.voom p.lfc1.voom
# q.voom q.lfc2.voom q.lfc2.voom
# mean(res$q.voom == p.adjust(res$p.voom, method='BH')) also 1.

# DSS
# p.notrend.DSS    p.trend.DSS
# q.notrend.DSS    q.trend.DSS
# lfdr.notrend.DSS lfdr.trend.DSS
# mean(res$q.notrend.DSS == p.adjust(res$p.notrend.DSS, method='BH')) is 0. Not clear what 
# fdr (q here) and lfdr from DSS are, but they don't seem to be increasing with p-value. Should also 
# compute BH q-values from p-values.

# baySeq
# prob.baySeq
# q.baySeq
### Need to re-run DEDD, DE with $Likelihood instead of $likes for prob.baySeq ###

# MDSeq
# p.mean.zi.MDSeq   p.mean.zi.lfc1.MDSeq   p.mean.zi.lfc2.MDSeq
# p.mean.nozi.MDSeq p.mean.nozi.lfc1.MDSeq p.mean.nozi.lfc2.MDSeq
# p.disp.zi.MDSeq   p.disp.zi.lfc1.MDSeq   p.disp.zi.lfc2.MDSeq
# p.disp.nozi.MDSeq p.disp.nozi.lfc1.MDSeq p.disp.nozi.lfc2.MDSeq
# Haven't included FDR-corrected values but should. Better to use output from method directly where possible.
# MDSeq uses p.adjust() with method ='BY' by default, but still best to use output from MDSeq.
### Re-run DEDD, DE to include FDR for MDSeq ###

# expHM
# prob.expHM prop.expHM
# p.mean.expHM p.lmean.expHM p.mean.lfc1.expHM p.mean.lfc2.expHM
# p.disp.expHM p.ldisp.expHM p.disp.lfc1.expHM p.disp.lfc2.expHM
# Only have p-values so need to compute q-values. May use BY instead of BH as seems to be less conservative.
# Also need to use proportions to define thresholds for posterior probability.

# lnHM
# prob.lnHM prop.lnHM
# p.mean.lnHM p.lmean.lnHM p.mean.lfc1.lnHM p.mean.lfc2.lnHM
# p.disp.lnHM p.ldisp.lnHM p.disp.lfc1.lnHM p.disp.lfc2.lnHM
# Only have p-values so need to compute q-values. May use BY instead of BH as seems to be less conservative.
# Also need to use proportions to define thresholds for posterior probability.




raw.mean.expHM <- res$p.mean.expHM < 0.05
raw.lmean.expHM <- res$p.lmean.expHM < 0.05
raw.disp.expHM <- res$p.disp.expHM < 0.05
raw.ldisp.expHM <- res$p.ldisp.expHM < 0.05
raw.mean.lnHM <- res$p.mean.lnHM < 0.05
raw.lmean.lnHM <- res$p.lmean.lnHM < 0.05
raw.disp.lnHM <- res$p.disp.lnHM < 0.05
raw.ldisp.lnHM <- res$p.ldisp.lnHM < 0.05
bh.mean.expHM <- p.adjust(res$p.mean.expHM, method='BH') < 0.05
bh.lmean.expHM <- p.adjust(res$p.lmean.expHM, method='BH') < 0.05
bh.disp.expHM <- p.adjust(res$p.disp.expHM, method='BH') < 0.05
bh.ldisp.expHM <- p.adjust(res$p.ldisp.expHM, method='BH') < 0.05
bh.mean.lnHM <- p.adjust(res$p.mean.lnHM, method='BH') < 0.05
bh.lmean.lnHM <- p.adjust(res$p.lmean.lnHM, method='BH') < 0.05
bh.disp.lnHM <- p.adjust(res$p.disp.lnHM, method='BH') < 0.05
bh.ldisp.lnHM <- p.adjust(res$p.ldisp.lnHM, method='BH') < 0.05
by.mean.expHM <- p.adjust(res$p.mean.expHM, method='BY') < 0.05
by.lmean.expHM <- p.adjust(res$p.lmean.expHM, method='BY') < 0.05
by.disp.expHM <- p.adjust(res$p.disp.expHM, method='BY') < 0.05
by.ldisp.expHM <- p.adjust(res$p.ldisp.expHM, method='BY') < 0.05
by.mean.lnHM <- p.adjust(res$p.mean.lnHM, method='BY') < 0.05
by.lmean.lnHM <- p.adjust(res$p.lmean.lnHM, method='BY') < 0.05
by.disp.lnHM <- p.adjust(res$p.disp.lnHM, method='BY') < 0.05
by.ldisp.lnHM <- p.adjust(res$p.ldisp.lnHM, method='BY') < 0.05
q.mean.expHM <- qvalue(res$p.mean.expHM)$qval < 0.05
q.lmean.expHM <- qvalue(res$p.lmean.expHM)$qval < 0.05
q.disp.expHM <- qvalue(res$p.disp.expHM)$qval < 0.05
q.ldisp.expHM <- qvalue(res$p.ldisp.expHM)$qval < 0.05
q.mean.lnHM <- qvalue(res$p.mean.lnHM)$qval < 0.05
q.lmean.lnHM <- qvalue(res$p.lmean.lnHM)$qval < 0.05
q.disp.lnHM <- qvalue(res$p.disp.lnHM)$qval < 0.05
q.ldisp.lnHM <- qvalue(res$p.ldisp.lnHM)$qval < 0.05

tp.mean.expHM <- c(sum(raw.mean.expHM==1 & res$DE==1), sum(bh.mean.expHM==1 & res$DE==1), 
                   sum(by.mean.expHM==1 & res$DE==1), sum(q.mean.expHM==1 & res$DE==1))
tp.lmean.expHM <- c(sum(raw.lmean.expHM==1 & res$DE==1), sum(bh.lmean.expHM==1 & res$DE==1), 
                    sum(by.lmean.expHM==1 & res$DE==1), sum(q.lmean.expHM==1 & res$DE==1))
tp.disp.expHM <- c(sum(raw.disp.expHM==1 & res$DD==1), sum(bh.disp.expHM==1 & res$DD==1), 
                   sum(by.disp.expHM==1 & res$DD==1), sum(q.disp.expHM==1 & res$DD==1))
tp.ldisp.expHM <- c(sum(raw.ldisp.expHM==1 & res$DD==1), sum(bh.ldisp.expHM==1 & res$DD==1), 
                    sum(by.ldisp.expHM==1 & res$DD==1), sum(q.ldisp.expHM==1 & res$DD==1))
tp.mean.lnHM <- c(sum(raw.mean.lnHM==1 & res$DE==1), sum(bh.mean.lnHM==1 & res$DE==1), 
                  sum(by.mean.lnHM==1 & res$DE==1), sum(q.mean.lnHM==1 & res$DE==1))
tp.lmean.lnHM <- c(sum(raw.lmean.lnHM==1 & res$DE==1), sum(bh.lmean.lnHM==1 & res$DE==1), 
                   sum(by.lmean.lnHM==1 & res$DE==1), sum(q.lmean.lnHM==1 & res$DE==1))
tp.disp.lnHM <- c(sum(raw.disp.lnHM==1 & res$DD==1), sum(bh.disp.lnHM==1 & res$DD==1), 
                  sum(by.disp.lnHM==1 & res$DD==1), sum(q.disp.lnHM==1 & res$DD==1))
tp.ldisp.lnHM <- c(sum(raw.ldisp.lnHM==1 & res$DD==1), sum(bh.ldisp.lnHM==1 & res$DD==1), 
                   sum(by.ldisp.lnHM==1 & res$DD==1), sum(q.ldisp.lnHM==1 & res$DD==1))
fp.mean.expHM <- c(sum(raw.mean.expHM==1 & res$DE==0), sum(bh.mean.expHM==1 & res$DE==0), 
                   sum(by.mean.expHM==1 & res$DE==0), sum(q.mean.expHM==1 & res$DE==0))
fp.lmean.expHM <- c(sum(raw.lmean.expHM==1 & res$DE==0), sum(bh.lmean.expHM==1 & res$DE==0), 
                    sum(by.lmean.expHM==1 & res$DE==0), sum(q.lmean.expHM==1 & res$DE==0))
fp.disp.expHM <- c(sum(raw.disp.expHM==1 & res$DD==0), sum(bh.disp.expHM==1 & res$DD==0), 
                   sum(by.disp.expHM==1 & res$DD==0), sum(q.disp.expHM==1 & res$DD==0))
fp.ldisp.expHM <- c(sum(raw.ldisp.expHM==1 & res$DD==0), sum(bh.ldisp.expHM==1 & res$DD==0), 
                    sum(by.ldisp.expHM==1 & res$DD==0), sum(q.ldisp.expHM==1 & res$DD==0))
fp.mean.lnHM <- c(sum(raw.mean.lnHM==1 & res$DE==0), sum(bh.mean.lnHM==1 & res$DE==0), 
                  sum(by.mean.lnHM==1 & res$DE==0), sum(q.mean.lnHM==1 & res$DE==0))
fp.lmean.lnHM <- c(sum(raw.lmean.lnHM==1 & res$DE==0), sum(bh.lmean.lnHM==1 & res$DE==0), 
                   sum(by.lmean.lnHM==1 & res$DE==0), sum(q.lmean.lnHM==1 & res$DE==0))
fp.disp.lnHM <- c(sum(raw.disp.lnHM==1 & res$DD==0), sum(bh.disp.lnHM==1 & res$DD==0), 
                  sum(by.disp.lnHM==1 & res$DD==0), sum(q.disp.lnHM==1 & res$DD==0))
fp.ldisp.lnHM <- c(sum(raw.ldisp.lnHM==1 & res$DD==0), sum(bh.ldisp.lnHM==1 & res$DD==0), 
                   sum(by.ldisp.lnHM==1 & res$DD==0), sum(q.ldisp.lnHM==1 & res$DD==0))
fn.mean.expHM <- c(sum(raw.mean.expHM==0 & res$DE==1), sum(bh.mean.expHM==0 & res$DE==1), 
                   sum(by.mean.expHM==0 & res$DE==1), sum(q.mean.expHM==0 & res$DE==1))
fn.lmean.expHM <- c(sum(raw.lmean.expHM==0 & res$DE==1), sum(bh.lmean.expHM==0 & res$DE==1), 
                    sum(by.lmean.expHM==0 & res$DE==1), sum(q.lmean.expHM==0 & res$DE==1))
fn.disp.expHM <- c(sum(raw.disp.expHM==0 & res$DD==1), sum(bh.disp.expHM==0 & res$DD==1), 
                   sum(by.disp.expHM==0 & res$DD==1), sum(q.disp.expHM==0 & res$DD==1))
fn.ldisp.expHM <- c(sum(raw.ldisp.expHM==0 & res$DD==1), sum(bh.ldisp.expHM==0 & res$DD==1), 
                    sum(by.ldisp.expHM==0 & res$DD==1), sum(q.ldisp.expHM==0 & res$DD==1))
fn.mean.lnHM <- c(sum(raw.mean.lnHM==0 & res$DE==1), sum(bh.mean.lnHM==0 & res$DE==1), 
                  sum(by.mean.lnHM==0 & res$DE==1), sum(q.mean.lnHM==0 & res$DE==1))
fn.lmean.lnHM <- c(sum(raw.lmean.lnHM==0 & res$DE==1), sum(bh.lmean.lnHM==0 & res$DE==1), 
                   sum(by.lmean.lnHM==0 & res$DE==1), sum(q.lmean.lnHM==0 & res$DE==1))
fn.disp.lnHM <- c(sum(raw.disp.lnHM==0 & res$DD==1), sum(bh.disp.lnHM==0 & res$DD==1), 
                  sum(by.disp.lnHM==0 & res$DD==1), sum(q.disp.lnHM==0 & res$DD==1))
fn.ldisp.lnHM <- c(sum(raw.ldisp.lnHM==0 & res$DD==1), sum(bh.ldisp.lnHM==0 & res$DD==1), 
                   sum(by.ldisp.lnHM==0 & res$DD==1), sum(q.ldisp.lnHM==0 & res$DD==1))
tn.mean.expHM <- c(sum(raw.mean.expHM==0 & res$DE==0), sum(bh.mean.expHM==0 & res$DE==0), 
                   sum(by.mean.expHM==0 & res$DE==0), sum(q.mean.expHM==0 & res$DE==0))
tn.lmean.expHM <- c(sum(raw.lmean.expHM==0 & res$DE==0), sum(bh.lmean.expHM==0 & res$DE==0), 
                    sum(by.lmean.expHM==0 & res$DE==0), sum(q.lmean.expHM==0 & res$DE==0))
tn.disp.expHM <- c(sum(raw.disp.expHM==0 & res$DD==0), sum(bh.disp.expHM==0 & res$DD==0), 
                   sum(by.disp.expHM==0 & res$DD==0), sum(q.disp.expHM==0 & res$DD==0))
tn.ldisp.expHM <- c(sum(raw.ldisp.expHM==0 & res$DD==0), sum(bh.ldisp.expHM==0 & res$DD==0), 
                    sum(by.ldisp.expHM==0 & res$DD==0), sum(q.ldisp.expHM==0 & res$DD==0))
tn.mean.lnHM <- c(sum(raw.mean.lnHM==0 & res$DE==0), sum(bh.mean.lnHM==0 & res$DE==0), 
                  sum(by.mean.lnHM==0 & res$DE==0), sum(q.mean.lnHM==0 & res$DE==0))
tn.lmean.lnHM <- c(sum(raw.lmean.lnHM==0 & res$DE==0), sum(bh.lmean.lnHM==0 & res$DE==0), 
                   sum(by.lmean.lnHM==0 & res$DE==0), sum(q.lmean.lnHM==0 & res$DE==0))
tn.disp.lnHM <- c(sum(raw.disp.lnHM==0 & res$DD==0), sum(bh.disp.lnHM==0 & res$DD==0), 
                  sum(by.disp.lnHM==0 & res$DD==0), sum(q.disp.lnHM==0 & res$DD==0))
tn.ldisp.lnHM <- c(sum(raw.ldisp.lnHM==0 & res$DD==0), sum(bh.ldisp.lnHM==0 & res$DD==0), 
                   sum(by.ldisp.lnHM==0 & res$DD==0), sum(q.ldisp.lnHM==0 & res$DD==0))

fdr.mean.expHM <- fp.mean.expHM/(fp.mean.expHM+tp.mean.expHM)
fdr.lmean.expHM <- fp.lmean.expHM/(fp.lmean.expHM+tp.lmean.expHM)
fdr.disp.expHM <- fp.disp.expHM/(fp.disp.expHM+tp.disp.expHM)
fdr.ldisp.expHM <- fp.ldisp.expHM/(fp.ldisp.expHM+tp.ldisp.expHM)
fdr.mean.lnHM <- fp.mean.lnHM/(fp.mean.lnHM+tp.mean.lnHM)
fdr.lmean.lnHM <- fp.lmean.lnHM/(fp.lmean.lnHM+tp.lmean.lnHM)
fdr.disp.lnHM <- fp.disp.lnHM/(fp.disp.lnHM+tp.disp.lnHM)
fdr.ldisp.lnHM <- fp.ldisp.lnHM/(fp.ldisp.lnHM+tp.ldisp.lnHM)
tpr.mean.expHM <- tp.mean.expHM/(tp.mean.expHM+fn.mean.expHM)
tpr.lmean.expHM <- tp.lmean.expHM/(tp.lmean.expHM+fn.lmean.expHM)
tpr.disp.expHM <- tp.disp.expHM/(tp.disp.expHM+fn.disp.expHM)
tpr.ldisp.expHM <- tp.ldisp.expHM/(tp.ldisp.expHM+fn.ldisp.expHM)
tpr.mean.lnHM <- tp.mean.lnHM/(tp.mean.lnHM+fn.mean.lnHM)
tpr.lmean.lnHM <- tp.lmean.lnHM/(tp.lmean.lnHM+fn.lmean.lnHM)
tpr.disp.lnHM <- tp.disp.lnHM/(tp.disp.lnHM+fn.disp.lnHM)
tpr.ldisp.lnHM <- tp.ldisp.lnHM/(tp.ldisp.lnHM+fn.ldisp.lnHM)

data.frame(method=c('none','BY','BH','q'), 
           fdr.mean.expHM, tpr.mean.expHM, fdr.lmean.expHM, tpr.lmean.expHM, 
           fdr.mean.lnHM, tpr.mean.lnHM, fdr.lmean.lnHM, tpr.lmean.lnHM)
# q best FDR control for mean but slightly liberal for lmean; BY best for lmean and always better than BH.
# TPR only slightly better for q than BY, both much better than BH.
data.frame(method=c('none','BY','BH','q'), 
           fdr.disp.expHM, tpr.disp.expHM, fdr.ldisp.expHM, tpr.ldisp.expHM, 
           fdr.disp.lnHM, tpr.disp.lnHM, fdr.ldisp.lnHM, tpr.ldisp.lnHM)
# No positives for any method.

pr.mean.expHM <- prediction(1-res$p.mean.expHM, res$DE); pr.lmean.expHM <- prediction(1-res$p.lmean.expHM, res$DE)
pr.disp.expHM <- prediction(1-res$p.disp.expHM, res$DD); pr.ldisp.expHM <- prediction(1-res$p.ldisp.expHM, res$DD)
pr.mean.lnHM <- prediction(1-res$p.mean.lnHM, res$DE); pr.lmean.lnHM <- prediction(1-res$p.lmean.lnHM, res$DE)
pr.disp.lnHM <- prediction(1-res$p.disp.lnHM, res$DD); pr.ldisp.lnHM <- prediction(1-res$p.ldisp.lnHM, res$DD)
auc.mean.expHM <- performance(pr.mean.expHM, measure='auc')@y.values[[1]]
auc.lmean.expHM <- performance(pr.lmean.expHM, measure='auc')@y.values[[1]]
auc.disp.expHM <- performance(pr.disp.expHM, measure='auc')@y.values[[1]]
auc.ldisp.expHM <- performance(pr.ldisp.expHM, measure='auc')@y.values[[1]]
auc.mean.lnHM <- performance(pr.mean.lnHM, measure='auc')@y.values[[1]]
auc.lmean.lnHM <- performance(pr.lmean.lnHM, measure='auc')@y.values[[1]]
auc.disp.lnHM <- performance(pr.disp.lnHM, measure='auc')@y.values[[1]]
auc.ldisp.lnHM <- performance(pr.ldisp.lnHM, measure='auc')@y.values[[1]]
c(auc.mean.expHM, auc.lmean.expHM, auc.mean.lnHM, auc.lmean.lnHM)
c(auc.disp.expHM, auc.ldisp.expHM, auc.disp.lnHM, auc.ldisp.lnHM)





names.edgeR <- c('ql.edgeR', 'ql.lfc1.edgeR', 'ql.lfc2.edgeR', 'lr.edgeR', 'lr.lfc1.edgeR', 'lr.lfc2.edgeR', 'et.edgeR')
names.DESeq2 <- c('noif.DESeq', 'if.DESeq', 'lfc1.DESeq', 'lfc2.DESeq')
names.voom <- c('voom', 'lfc1.voom', 'lfc2.voom')
names.DSS <- c('notrend.DSS', 'trend.DSS')
names.baySeq <- 'baySeq'
names.ShrinkBayes <- c('mix.ShrinkBayes', 'np.ShrinkBayes', 'np.lfc1.ShrinkBayes', 'np.lfc2.ShrinkBayes')
names.MDSeq <- c('mean.zi.MDSeq', 'mean.nozi.MDSeq', 'mean.zi.lfc1.MDSeq', 'mean.nozi.lfc1.MDSeq', 
                 'mean.zi.lfc2.MDSeq', 'mean.nozi.lfc2.MDSeq', 
                 'disp.zi.MDSeq', 'disp.nozi.MDSeq', 'disp.zi.lfc1.MDSeq', 'disp.nozi.lfc1.MDSeq', 
                 'disp.zi.lfc2.MDSeq', 'disp.nozi.lfc2.MDSeq')
names.expHM <- c('hmm.expHM', 'mean.expHM', 'lmean.expHM', 'mean.lfc1.expHM', 'mean.lfc2.expHM', 
                 'disp.expHM', 'ldisp.expHM', 'disp.lfc1.expHM', 'disp.lfc2.expHM')
names.lnHM <- c('hmm.lnHM', 'mean.lnHM', 'lmean.lnHM', 'mean.lfc1.lnHM', 'mean.lfc2.lnHM', 
                'disp.lnHM', 'ldisp.lnHM', 'disp.lfc1.lnHM', 'disp.lfc2.lnHM')
results.list <- c(names.edgeR, names.DESeq2, names.voom, names.DSS, names.baySeq, names.ShrinkBayes, 
                  names.MDSeq, names.expHM, names.lnHM)
for (i in c(names.edgeR, names.DESeq2, names.voom, names.DSS, names.baySeq, names.MDSeq)) {
  for (j in c('fpr.raw.','fdr.raw.','tpr.raw.','fpr.fdr.','fdr.fdr.','tpr.fdr.')) {
    assign(paste0(j,i), numeric(50))
  }
}
for (i in names.DSS) {
  for (j in c('fpr.lfdr.','fdr.lfdr.','tpr.lfdr.','fpr.bh.','fdr.bh.','tpr.bh.')) {
    assign(paste0(j,i), numeric(50))
  }
}
for (i in names.ShrinkBayes) {
  for (j in c('fpr.lfdr.','fdr.lfdr.','tpr.lfdr.','fpr.bfdr.','fdr.bfdr.','tpr.bfdr.')) {
    assign(paste0(j,i), numeric(50))
  }
}
for (i in c(names.expHM[1],names.lnHM[1])) {
  for (j in c('fpr.5.','fdr.5.','tpr.5.','fpr.thr.','fdr.thr.','tpr.thr.','fpr.bfdr.','fdr.bfdr.','tpr.bfdr.')) {
    assign(paste0(j,i), numeric(50))
  }
}
for (i in c(names.expHM[-1],names.lnHM[-1])) {
  for (j in c('fpr.raw.','fdr.raw.','tpr.raw.','fpr.bh.','fdr.bh.','tpr.bh.', 
              'fpr.by.','fdr.by.','tpr.by.','fpr.q.','fdr.q.','tpr.q.')) {
    assign(paste0(j,i), numeric(50))
  }
}
for (i in results.list) {
  assign(paste0('auc.',i), numeric(50))
  assign(paste0('discoveries.',i), matrix(nrow=50, ncol=20000))
  assign(paste0('false.discoveries.',i), matrix(nrow=50, ncol=20000))
}
for (i in 1:50) {
  import <- paste0('results.DEDD2.',i,'.rds')
  res <- readRDS(here('Results/DEDD compcodeR data results July 2019',import))
  DE <- factor(res$DE, levels=c('1','0')); DD <- factor(res$DD, levels=c('1','0'))
  lfcm1 <- factor(as.numeric(res$lfcm1), levels=c('1','0'))
  lfcm2 <- factor(as.numeric(res$lfcm2), levels=c('1','0'))
  lfcd1 <- factor(as.numeric(res$lfcd1), levels=c('1','0'))
  lfcd2 <- factor(as.numeric(res$lfcd2), levels=c('1','0'))
  for (j in names.edgeR[-1]) {
    assign(paste0('call.raw.',j), factor(as.numeric(get(paste0('p.',j), res) < 0.05), levels=c('1','0')))
    assign(paste0('call.fdr.',j), factor(as.numeric(get(paste0('q.',j), res) < 0.05), levels=c('1','0')))
    if (grepl('lfc1',j)) {
      assign(paste0('pred.',j), prediction(1-get(paste0('p.',j), res), lfcm1))
      assign(paste0('fpr.raw.',j), `[<-`(get(paste0('fpr.raw.',j)), i, 
                                         value=1-specificity(get(paste0('call.raw.',j)), lfcm1)))
      assign(paste0('fpr.fdr.',j), `[<-`(get(paste0('fpr.fdr.',j)), i, 
                                         value=1-specificity(get(paste0('call.fdr.',j)), lfcm1)))
      assign(paste0('fdr.raw.',j), `[<-`(get(paste0('fdr.raw.',j)), i, 
                                         value=1-precision(get(paste0('call.raw.',j)), lfcm1)))
      assign(paste0('fdr.fdr.',j), `[<-`(get(paste0('fdr.fdr.',j)), i, 
                                         value=1-precision(get(paste0('call.fdr.',j)), lfcm1)))
      assign(paste0('tpr.raw.',j), `[<-`(get(paste0('tpr.raw.',j)), i, 
                                         value=sensitivity(get(paste0('call.raw.',j)), lfcm1)))
      assign(paste0('tpr.fdr.',j), `[<-`(get(paste0('tpr.fdr.',j)), i, 
                                         value=sensitivity(get(paste0('call.fdr.',j)), lfcm1)))
    }
    else if (grepl('lfc2',j)) {
      assign(paste0('pred.',j), prediction(1-get(paste0('p.',j), res), lfcm2))
      assign(paste0('fpr.raw.',j), `[<-`(get(paste0('fpr.raw.',j)), i, 
                                         value=1-specificity(get(paste0('call.raw.',j)), lfcm2)))
      assign(paste0('fpr.fdr.',j), `[<-`(get(paste0('fpr.fdr.',j)), i, 
                                         value=1-specificity(get(paste0('call.fdr.',j)), lfcm2)))
      assign(paste0('fdr.raw.',j), `[<-`(get(paste0('fdr.raw.',j)), i, 
                                         value=1-precision(get(paste0('call.raw.',j)), lfcm2)))
      assign(paste0('fdr.fdr.',j), `[<-`(get(paste0('fdr.fdr.',j)), i, 
                                         value=1-precision(get(paste0('call.fdr.',j)), lfcm2)))
      assign(paste0('tpr.raw.',j), `[<-`(get(paste0('tpr.raw.',j)), i, 
                                         value=sensitivity(get(paste0('call.raw.',j)), lfcm2)))
      assign(paste0('tpr.fdr.',j), `[<-`(get(paste0('tpr.fdr.',j)), i, 
                                         value=sensitivity(get(paste0('call.fdr.',j)), lfcm2)))
    }
    else {
      assign(paste0('pred.',j), prediction(1-get(paste0('p.',j), res), DE))
      assign(paste0('fpr.raw.',j), `[<-`(get(paste0('fpr.raw.',j)), i, 
                                         value=1-specificity(get(paste0('call.raw.',j)), DE)))
      assign(paste0('fpr.fdr.',j), `[<-`(get(paste0('fpr.fdr.',j)), i, 
                                         value=1-specificity(get(paste0('call.fdr.',j)), DE)))
      assign(paste0('fdr.raw.',j), `[<-`(get(paste0('fdr.raw.',j)), i, 
                                         value=1-precision(get(paste0('call.raw.',j)), DE)))
      assign(paste0('fdr.fdr.',j), `[<-`(get(paste0('fdr.fdr.',j)), i, 
                                         value=1-precision(get(paste0('call.fdr.',j)), DE)))
      assign(paste0('tpr.raw.',j), `[<-`(get(paste0('tpr.raw.',j)), i, 
                                         value=sensitivity(get(paste0('call.raw.',j)), DE)))
      assign(paste0('tpr.fdr.',j), `[<-`(get(paste0('tpr.fdr.',j)), i, 
                                         value=sensitivity(get(paste0('call.fdr.',j)), DE)))
    }
    assign(paste0('auc.',j), `[<-`(get(paste0('auc.',j)), i, 
                                   value=performance(get(paste0('pred.',j)), measure='auc')@y.values[[1]]))
    assign(paste0('false.discoveries.',j), `[<-`(get(paste0('false.discoveries.',j)), i,, 
                                                  value=c(get(paste0('pred.',j))@fp[[1]], 
                                                          rep(NA, 20000-length(get(paste0('pred.',j))@fp[[1]])))))
    assign(paste0('discoveries.',j), `[<-`(get(paste0('discoveries.',j)), i,, 
                                                  value=c(get(paste0('pred.',j))@n.pos.pred[[1]], 
                                                          rep(NA, 20000-length(get(paste0('pred.',j))@n.pos.pred[[1]])))))
  }
}
#mean.discoveries.ql.edgeR <- colMeans(discoveries.ql.edgeR)
mean.discoveries.ql.lfc1.edgeR <- colMeans(discoveries.ql.lfc1.edgeR)
mean.discoveries.ql.lfc2.edgeR <- colMeans(discoveries.ql.lfc2.edgeR)
mean.discoveries.lr.edgeR <- colMeans(discoveries.lr.edgeR)
mean.discoveries.lr.lfc1.edgeR <- colMeans(discoveries.lr.lfc1.edgeR)
mean.discoveries.lr.lfc2.edgeR <- colMeans(discoveries.lr.lfc2.edgeR)
mean.discoveries.et.edgeR <- colMeans(discoveries.et.edgeR)
#mean.fdr.ql.edgeR <- colMeans(false.discoveries.ql.edgeR) / colMeans(discoveries.ql.edgeR)
mean.fdr.ql.lfc1.edgeR <- colMeans(false.discoveries.ql.lfc1.edgeR) / colMeans(discoveries.ql.lfc1.edgeR)
mean.fdr.ql.lfc2.edgeR <- colMeans(false.discoveries.ql.lfc2.edgeR) / colMeans(discoveries.ql.lfc2.edgeR)
mean.fdr.lr.edgeR <- colMeans(false.discoveries.lr.edgeR) / colMeans(discoveries.lr.edgeR)
mean.fdr.lr.lfc1.edgeR <- colMeans(false.discoveries.lr.lfc1.edgeR) / colMeans(discoveries.lr.lfc1.edgeR)
mean.fdr.lr.lfc2.edgeR <- colMeans(false.discoveries.lr.lfc2.edgeR) / colMeans(discoveries.lr.lfc2.edgeR)
mean.fdr.et.edgeR <- colMeans(false.discoveries.et.edgeR) / colMeans(discoveries.et.edgeR)

par(mfrow=c(4,3))
#plot(mean.discoveries.ql.edgeR, mean.fdr.ql.edgeR, type='l'); abline(h=0.05)
plot(mean.discoveries.ql.lfc1.edgeR, mean.fdr.ql.lfc1.edgeR, type='l'); abline(h=0.05)
plot(mean.discoveries.ql.lfc2.edgeR, mean.fdr.ql.lfc2.edgeR, type='l'); abline(h=0.05)
plot(mean.discoveries.lr.edgeR, mean.fdr.lr.edgeR, type='l'); abline(h=0.05)
plot(mean.discoveries.lr.lfc1.edgeR, mean.fdr.lr.lfc1.edgeR, type='l'); abline(h=0.05)
plot(mean.discoveries.lr.lfc2.edgeR, mean.fdr.lr.lfc2.edgeR, type='l'); abline(h=0.05)
plot(mean.discoveries.et.edgeR, mean.fdr.et.edgeR, type='l'); abline(h=0.05)

par(mfcol=c(4,2), mar=c(3,3,2,2))
boxplot(fpr.raw.ql.edgeR, fpr.raw.lr.edgeR, fpr.raw.et.edgeR, fpr.raw.ql.lfc1.edgeR, 
        fpr.raw.lr.lfc1.edgeR, fpr.raw.ql.lfc2.edgeR, fpr.raw.lr.lfc2.edgeR, 
        fpr.fdr.ql.edgeR, fpr.fdr.lr.edgeR, fpr.fdr.et.edgeR, fpr.fdr.ql.lfc1.edgeR, 
        fpr.fdr.lr.lfc1.edgeR, fpr.fdr.ql.lfc2.edgeR, fpr.fdr.lr.lfc2.edgeR)
boxplot(fdr.raw.ql.edgeR, fdr.raw.lr.edgeR, fdr.raw.et.edgeR, fdr.raw.ql.lfc1.edgeR, 
        fdr.raw.lr.lfc1.edgeR, fdr.raw.ql.lfc2.edgeR, fdr.raw.lr.lfc2.edgeR, 
        fdr.fdr.ql.edgeR, fdr.fdr.lr.edgeR, fdr.fdr.et.edgeR, fdr.fdr.ql.lfc1.edgeR, 
        fdr.fdr.lr.lfc1.edgeR, fdr.fdr.ql.lfc2.edgeR, fdr.fdr.lr.lfc2.edgeR)
boxplot(tpr.raw.ql.edgeR, tpr.raw.lr.edgeR, tpr.raw.et.edgeR, tpr.raw.ql.lfc1.edgeR, 
        tpr.raw.lr.lfc1.edgeR, tpr.raw.ql.lfc2.edgeR, tpr.raw.lr.lfc2.edgeR, 
        tpr.fdr.ql.edgeR, tpr.fdr.lr.edgeR, tpr.fdr.et.edgeR, tpr.fdr.ql.lfc1.edgeR, 
        tpr.fdr.lr.lfc1.edgeR, tpr.fdr.ql.lfc2.edgeR, tpr.fdr.lr.lfc2.edgeR)
boxplot(auc.ql.edgeR, auc.lr.edgeR, auc.et.edgeR, auc.ql.lfc1.edgeR, 
        auc.lr.lfc1.edgeR, auc.ql.lfc2.edgeR, auc.lr.lfc2.edgeR)


