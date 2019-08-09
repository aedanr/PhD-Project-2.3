# put results together, i.e. from ShrinkBayes and others (can't apply yet)
# how to organise results to get metrics of interest in most useful format
# 50 values for each metric for each method for each dataset
# but first need to get FDR values that aren't included in each method's own analysis

library(here)
library(compcodeR)
for (i in 1:50) {
#  res.ShrinkBayes <- readRDS(here('Results/DEDD compcodeR data results July 2019', 
#                                  paste0('results.DEDD2.',i,'.ShrinkBayes.rds')))
#  res.ShrinkBayes <- res.ShrinkBayes[-(1:7)]
  res.edgeR_baySeq_MDSeq <- readRDS(here('Results/DE compcodeR data results July 2019', 
                                  paste0('results.DE2.',i,'.edgeR_baySeq_MDSeq.rds')))
  res.lfc1.lnHM_lfc2.lnHM <- readRDS(here('Results/DE compcodeR data results July 2019', 
                                          paste0('results.DE2.',i,'.lfc1.lnHM_lfc2.rds')))
  res.lfc1.lnHM_lfc2.lnHM <- res.lfc1.lnHM_lfc2.lnHM[-(1:4)]
  res.rest <- readRDS(here('Results/DE compcodeR data results July 2019', 
                                  paste0('results.DE2.',i,'.rds')))
  overlap <- which(names(res.rest) %in% names(res.edgeR_baySeq_MDSeq))
  res.rest <- res.rest[-overlap]
  res.all <- c(res.edgeR_baySeq_MDSeq, res.rest, res.lfc1.lnHM_lfc2.lnHM)
  saveRDS(res.all, file=here('Results/DE compcodeR data results July 2019', 
                             paste0('results.DE2.',i,'.all.rds')))
}
rm(list=c('i', 'res.ShrinkBayes', 'res.edgeR_baySeq_MDSeq', 'res.rest', 'overlap', 'res.all'))


library(ROCR)
library(caret)
library(qvalue)
source(here('scripts','2019-05-03_bfdr_function.R'))

# Data:
# data DE DD lfcm1 lfcm2 lfcd1 lfcd2

# edgeR:
# p.ql.edgeR p.ql.lfc1.edgeR p.ql.lfc2.edgeR p.lr.edgeR p.lr.lfc1.edgeR p.lr.lfc2.edgeR p.et.edgeR
# q.ql.edgeR q.ql.lfc1.edgeR q.ql.lfc2.edgeR q.lr.edgeR q.lr.lfc1.edgeR q.lr.lfc2.edgeR q.et.edgeR
# edgeR user guide says it uses BH.
# mean(res$q.lr.edgeR == p.adjust(res$p.lr.edgeR, method='BH')) is 1, so confirms.

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
# trend in DEDD2, not DEDD5, DEDD10; in DE2, DE5, not DE10

# baySeq
# prob.baySeq
# q.baySeq

# ShrinkBayes (only DEDD2)
# lfdr.mix.ShrinkBayes lfdr.np.ShrinkBayes lfdr.np.lfc1.ShrinkBayes lfdr.np.lfc2.ShrinkBayes
# bfdr.mix.ShrinkBayes bfdr.np.ShrinkBayes bfdr.np.lfc1.ShrinkBayes bfdr.np.lfc2.ShrinkBayes

# MDSeq
# p.mean.zi.MDSeq   p.mean.zi.lfc1.MDSeq   p.mean.zi.lfc2.MDSeq
# p.mean.nozi.MDSeq p.mean.nozi.lfc1.MDSeq p.mean.nozi.lfc2.MDSeq
# p.disp.zi.MDSeq   p.disp.zi.lfc1.MDSeq   p.disp.zi.lfc2.MDSeq
# p.disp.nozi.MDSeq p.disp.nozi.lfc1.MDSeq p.disp.nozi.lfc2.MDSeq
# q.mean.zi.MDSeq   q.mean.zi.lfc1.MDSeq   q.mean.zi.lfc2.MDSeq
# q.mean.nozi.MDSeq q.mean.nozi.lfc1.MDSeq q.mean.nozi.lfc2.MDSeq
# q.disp.zi.MDSeq   q.disp.zi.lfc1.MDSeq   q.disp.zi.lfc2.MDSeq
# q.disp.nozi.MDSeq q.disp.nozi.lfc1.MDSeq q.disp.nozi.lfc2.MDSeq
# Haven't included FDR-corrected values but should. Better to use output from method directly where possible.
# MDSeq uses p.adjust() with method ='BY' by default, but still best to use output from MDSeq.

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
names.expHM <- c('expHM', 'mean.expHM', 'lmean.expHM', 'mean.lfc1.expHM', 'mean.lfc2.expHM', 
                 'disp.expHM', 'ldisp.expHM', 'disp.lfc1.expHM', 'disp.lfc2.expHM')
names.lnHM <- c('lnHM', 'mean.lnHM', 'lmean.lnHM', 'mean.lfc1.lnHM', 'mean.lfc2.lnHM', 
                'disp.lnHM', 'ldisp.lnHM', 'disp.lfc1.lnHM', 'disp.lfc2.lnHM')
results.list <- c(names.edgeR, names.DESeq2, names.voom, names.DSS, names.baySeq, names.ShrinkBayes, 
                  names.MDSeq, names.expHM, names.lnHM)
for (i in c(names.edgeR, names.DESeq2, names.voom, names.DSS, names.MDSeq)) {
  for (j in c('fpr.raw.','fdr.raw.','tpr.raw.','fpr.fdr.','fdr.fdr.','tpr.fdr.')) {
    assign(paste0(j,i), numeric(50))
  }
}
for (i in names.DSS) {
  for (j in c('fpr.lfdr.','fdr.lfdr.','tpr.lfdr.','fpr.bh.','fdr.bh.','tpr.bh.')) {
    assign(paste0(j,i), numeric(50))
  }
}
for (i in names.baySeq) {
  for (j in c('fpr.5.','fdr.5.','tpr.5.','fpr.fdr.','fdr.fdr.','tpr.fdr.')) {
    assign(paste0(j,i), numeric(50))
  }
}
for (i in names.ShrinkBayes) { 
  for (j in c('fpr.lfdr.','fdr.lfdr.','tpr.lfdr.','fpr.bfdr.','fdr.bfdr.','tpr.bfdr.')) {
    assign(paste0(j,i), numeric(50))
  } # Only for DEDD2
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
  import <- paste0('results.DEDD2.',i,'.all.rds')
  res <- readRDS(here('Results/DEDD compcodeR data results July 2019',import))
  DE <- factor(res$DE, levels=c('1','0')); DD <- factor(res$DD, levels=c('1','0'))
  DEDD <- factor(as.numeric(res$DE==1 | res$DD==1), levels=c('1','0'))
  lfcm1 <- factor(as.numeric(res$lfcm1), levels=c('1','0'))
  lfcm2 <- factor(as.numeric(res$lfcm2), levels=c('1','0'))
  lfcd1 <- factor(as.numeric(res$lfcd1), levels=c('1','0'))
  lfcd2 <- factor(as.numeric(res$lfcd2), levels=c('1','0'))
  bh.notrend.DSS <- p.adjust(res$p.notrend.DSS, method='BH')
  bh.trend.DSS <- p.adjust(res$p.trend.DSS, method='BH')
  q.mean.expHM <- qvalue(res$p.mean.expHM)$qval
  for (j in c(names.expHM[-1], names.lnHM[-1])) {
    assign(paste0('bh.',j), p.adjust(get(paste0('p.',j), res), method='BH'))
    assign(paste0('by.',j), p.adjust(get(paste0('p.',j), res), method='BY'))
    assign(paste0('q.',j), qvalue(get(paste0('p.',j), res))$qval)
  }
  thr.expHM <- sort(res$prob.expHM, decreasing=T)[round(nrow(res$data@count.matrix)*res$prop.expHM)]
  thr.lnHM <- sort(res$prob.lnHM, decreasing=T)[round(nrow(res$data@count.matrix)*res$prop.lnHM)]
  for (j in c(names.edgeR, names.DESeq2, names.voom, names.DSS, names.MDSeq)) {
    assign(paste0('call.raw.',j), factor(as.numeric(get(paste0('p.',j), res) < 0.05), levels=c('1','0')))
    assign(paste0('call.fdr.',j), factor(as.numeric(get(paste0('q.',j), res) < 0.05), levels=c('1','0')))
    if (grepl('disp',j)) {
      if (grepl('lfc1',j)) {
        assign(paste0('pred.',j), prediction(1-get(paste0('p.',j), res), lfcd1))
        assign(paste0('fpr.raw.',j), `[<-`(get(paste0('fpr.raw.',j)), i, 
                                           value=1-specificity(get(paste0('call.raw.',j)), lfcd1)))
        assign(paste0('fpr.fdr.',j), `[<-`(get(paste0('fpr.fdr.',j)), i, 
                                           value=1-specificity(get(paste0('call.fdr.',j)), lfcd1)))
        assign(paste0('fdr.raw.',j), `[<-`(get(paste0('fdr.raw.',j)), i, 
                                           value=1-precision(get(paste0('call.raw.',j)), lfcd1)))
        assign(paste0('fdr.fdr.',j), `[<-`(get(paste0('fdr.fdr.',j)), i, 
                                           value=1-precision(get(paste0('call.fdr.',j)), lfcd1)))
        assign(paste0('tpr.raw.',j), `[<-`(get(paste0('tpr.raw.',j)), i, 
                                           value=sensitivity(get(paste0('call.raw.',j)), lfcd1)))
        assign(paste0('tpr.fdr.',j), `[<-`(get(paste0('tpr.fdr.',j)), i, 
                                           value=sensitivity(get(paste0('call.fdr.',j)), lfcd1)))
      }
      else if (grepl('lfc2',j)) {
        assign(paste0('pred.',j), prediction(1-get(paste0('p.',j), res), lfcd2))
        assign(paste0('fpr.raw.',j), `[<-`(get(paste0('fpr.raw.',j)), i, 
                                           value=1-specificity(get(paste0('call.raw.',j)), lfcd2)))
        assign(paste0('fpr.fdr.',j), `[<-`(get(paste0('fpr.fdr.',j)), i, 
                                           value=1-specificity(get(paste0('call.fdr.',j)), lfcd2)))
        assign(paste0('fdr.raw.',j), `[<-`(get(paste0('fdr.raw.',j)), i, 
                                           value=1-precision(get(paste0('call.raw.',j)), lfcd2)))
        assign(paste0('fdr.fdr.',j), `[<-`(get(paste0('fdr.fdr.',j)), i, 
                                           value=1-precision(get(paste0('call.fdr.',j)), lfcd2)))
        assign(paste0('tpr.raw.',j), `[<-`(get(paste0('tpr.raw.',j)), i, 
                                           value=sensitivity(get(paste0('call.raw.',j)), lfcd2)))
        assign(paste0('tpr.fdr.',j), `[<-`(get(paste0('tpr.fdr.',j)), i, 
                                           value=sensitivity(get(paste0('call.fdr.',j)), lfcd2)))
      }
      else {
        assign(paste0('pred.',j), prediction(1-get(paste0('p.',j), res), DD))
        assign(paste0('fpr.raw.',j), `[<-`(get(paste0('fpr.raw.',j)), i, 
                                           value=1-specificity(get(paste0('call.raw.',j)), DD)))
        assign(paste0('fpr.fdr.',j), `[<-`(get(paste0('fpr.fdr.',j)), i, 
                                           value=1-specificity(get(paste0('call.fdr.',j)), DD)))
        assign(paste0('fdr.raw.',j), `[<-`(get(paste0('fdr.raw.',j)), i, 
                                           value=1-precision(get(paste0('call.raw.',j)), DD)))
        assign(paste0('fdr.fdr.',j), `[<-`(get(paste0('fdr.fdr.',j)), i, 
                                           value=1-precision(get(paste0('call.fdr.',j)), DD)))
        assign(paste0('tpr.raw.',j), `[<-`(get(paste0('tpr.raw.',j)), i, 
                                           value=sensitivity(get(paste0('call.raw.',j)), DD)))
        assign(paste0('tpr.fdr.',j), `[<-`(get(paste0('tpr.fdr.',j)), i, 
                                           value=sensitivity(get(paste0('call.fdr.',j)), DD)))
      }
    }
    else if (grepl('lfc1',j)) {
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
  for (j in names.DSS) {
    assign(paste0('call.bh.',j), factor(as.numeric(get(paste0('bh.',j)) < 0.05), levels=c('1','0')))
    assign(paste0('call.lfdr.',j), factor(as.numeric(get(paste0('lfdr.',j), res) < 0.05), levels=c('1','0')))
    assign(paste0('fpr.bh.',j), `[<-`(get(paste0('fpr.bh.',j)), i, 
                                      value=1-specificity(get(paste0('call.bh.',j)), DE)))
    assign(paste0('fpr.lfdr.',j), `[<-`(get(paste0('fpr.lfdr.',j)), i, 
                                        value=1-specificity(get(paste0('call.lfdr.',j)), DE)))
    assign(paste0('fdr.bh.',j), `[<-`(get(paste0('fdr.bh.',j)), i, 
                                      value=1-precision(get(paste0('call.bh.',j)), DE)))
    assign(paste0('fdr.lfdr.',j), `[<-`(get(paste0('fdr.lfdr.',j)), i, 
                                        value=1-precision(get(paste0('call.lfdr.',j)), DE)))
    assign(paste0('tpr.bh.',j), `[<-`(get(paste0('tpr.bh.',j)), i, 
                                      value=sensitivity(get(paste0('call.bh.',j)), DE)))
    assign(paste0('tpr.lfdr.',j), `[<-`(get(paste0('tpr.lfdr.',j)), i, 
                                        value=sensitivity(get(paste0('call.lfdr.',j)), DE)))
  }
  for (j in names.baySeq) {
    assign(paste0('call.5.',j), factor(as.numeric(get(paste0('prob.',j), res) > 0.5), levels=c('1','0')))
    assign(paste0('call.fdr.',j), factor(as.numeric(get(paste0('q.',j), res) < 0.05), levels=c('1','0')))
    assign(paste0('pred.',j), prediction(get(paste0('prob.',j), res), DE))
    assign(paste0('fpr.5.',j), `[<-`(get(paste0('fpr.5.',j)), i, 
                                     value=1-specificity(get(paste0('call.5.',j)), DE)))
    assign(paste0('fpr.fdr.',j), `[<-`(get(paste0('fpr.fdr.',j)), i, 
                                       value=1-specificity(get(paste0('call.fdr.',j)), DE)))
    assign(paste0('fdr.5.',j), `[<-`(get(paste0('fdr.5.',j)), i, 
                                     value=1-precision(get(paste0('call.5.',j)), DE)))
    assign(paste0('fdr.fdr.',j), `[<-`(get(paste0('fdr.fdr.',j)), i, 
                                       value=1-precision(get(paste0('call.fdr.',j)), DE)))
    assign(paste0('tpr.5.',j), `[<-`(get(paste0('tpr.5.',j)), i, 
                                     value=sensitivity(get(paste0('call.5.',j)), DE)))
    assign(paste0('tpr.fdr.',j), `[<-`(get(paste0('tpr.fdr.',j)), i, 
                                       value=sensitivity(get(paste0('call.fdr.',j)), DE)))
    assign(paste0('auc.',j), `[<-`(get(paste0('auc.',j)), i, 
                                   value=performance(get(paste0('pred.',j)), measure='auc')@y.values[[1]]))
    assign(paste0('false.discoveries.',j), `[<-`(get(paste0('false.discoveries.',j)), i,, 
                                                 value=c(get(paste0('pred.',j))@fp[[1]], 
                                                         rep(NA, 20000-length(get(paste0('pred.',j))@fp[[1]])))))
    assign(paste0('discoveries.',j), `[<-`(get(paste0('discoveries.',j)), i,, 
                                           value=c(get(paste0('pred.',j))@n.pos.pred[[1]], 
                                                   rep(NA, 20000-length(get(paste0('pred.',j))@n.pos.pred[[1]])))))
  }
  for (j in names.ShrinkBayes) {
    assign(paste0('call.lfdr.',j), factor(as.numeric(get(paste0('lfdr.',j), res) < 0.05), levels=c('1','0')))
    assign(paste0('call.bfdr.',j), factor(as.numeric(get(paste0('bfdr.',j), res) < 0.05), levels=c('1','0')))
    if (grepl('lfc1',j)) {
      assign(paste0('pred.',j), prediction(1-get(paste0('lfdr.',j), res), lfcm1))
      assign(paste0('fpr.lfdr.',j), `[<-`(get(paste0('fpr.lfdr.',j)), i, 
                                          value=1-specificity(get(paste0('call.lfdr.',j)), lfcm1)))
      assign(paste0('fpr.bfdr.',j), `[<-`(get(paste0('fpr.bfdr.',j)), i, 
                                          value=1-specificity(get(paste0('call.bfdr.',j)), lfcm1)))
      assign(paste0('fdr.lfdr.',j), `[<-`(get(paste0('fdr.lfdr.',j)), i, 
                                          value=1-precision(get(paste0('call.lfdr.',j)), lfcm1)))
      assign(paste0('fdr.bfdr.',j), `[<-`(get(paste0('fdr.bfdr.',j)), i, 
                                          value=1-precision(get(paste0('call.bfdr.',j)), lfcm1)))
      assign(paste0('tpr.lfdr.',j), `[<-`(get(paste0('tpr.lfdr.',j)), i, 
                                          value=sensitivity(get(paste0('call.lfdr.',j)), lfcm1)))
      assign(paste0('tpr.bfdr.',j), `[<-`(get(paste0('tpr.bfdr.',j)), i, 
                                          value=sensitivity(get(paste0('call.bfdr.',j)), lfcm1)))
    }
    else if (grepl('lfc2',j)) {
      assign(paste0('pred.',j), prediction(1-get(paste0('lfdr.',j), res), lfcm2))
      assign(paste0('fpr.lfdr.',j), `[<-`(get(paste0('fpr.lfdr.',j)), i, 
                                          value=1-specificity(get(paste0('call.lfdr.',j)), lfcm2)))
      assign(paste0('fpr.bfdr.',j), `[<-`(get(paste0('fpr.bfdr.',j)), i, 
                                          value=1-specificity(get(paste0('call.bfdr.',j)), lfcm2)))
      assign(paste0('fdr.lfdr.',j), `[<-`(get(paste0('fdr.lfdr.',j)), i, 
                                          value=1-precision(get(paste0('call.lfdr.',j)), lfcm2)))
      assign(paste0('fdr.bfdr.',j), `[<-`(get(paste0('fdr.bfdr.',j)), i, 
                                          value=1-precision(get(paste0('call.bfdr.',j)), lfcm2)))
      assign(paste0('tpr.lfdr.',j), `[<-`(get(paste0('tpr.lfdr.',j)), i, 
                                          value=sensitivity(get(paste0('call.lfdr.',j)), lfcm2)))
      assign(paste0('tpr.bfdr.',j), `[<-`(get(paste0('tpr.bfdr.',j)), i, 
                                          value=sensitivity(get(paste0('call.bfdr.',j)), lfcm2)))
    }
    else {
      assign(paste0('pred.',j), prediction(1-get(paste0('lfdr.',j), res), DE))
      assign(paste0('fpr.lfdr.',j), `[<-`(get(paste0('fpr.lfdr.',j)), i, 
                                          value=1-specificity(get(paste0('call.lfdr.',j)), DE)))
      assign(paste0('fpr.bfdr.',j), `[<-`(get(paste0('fpr.bfdr.',j)), i, 
                                          value=1-specificity(get(paste0('call.bfdr.',j)), DE)))
      assign(paste0('fdr.lfdr.',j), `[<-`(get(paste0('fdr.lfdr.',j)), i, 
                                          value=1-precision(get(paste0('call.lfdr.',j)), DE)))
      assign(paste0('fdr.bfdr.',j), `[<-`(get(paste0('fdr.bfdr.',j)), i, 
                                          value=1-precision(get(paste0('call.bfdr.',j)), DE)))
      assign(paste0('tpr.lfdr.',j), `[<-`(get(paste0('tpr.lfdr.',j)), i, 
                                          value=sensitivity(get(paste0('call.lfdr.',j)), DE)))
      assign(paste0('tpr.bfdr.',j), `[<-`(get(paste0('tpr.bfdr.',j)), i, 
                                          value=sensitivity(get(paste0('call.bfdr.',j)), DE)))
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
  for (j in c(names.expHM[-1],names.lnHM[-1])) {
    assign(paste0('call.raw.',j), factor(as.numeric(get(paste0('p.',j), res) < 0.05), levels=c('1','0')))
    assign(paste0('call.q.',j), factor(as.numeric(get(paste0('q.',j)) < 0.05), levels=c('1','0')))
    assign(paste0('call.by.',j), factor(as.numeric(get(paste0('by.',j)) < 0.05), levels=c('1','0')))
    assign(paste0('call.bh.',j), factor(as.numeric(get(paste0('bh.',j)) < 0.05), levels=c('1','0')))
    if (grepl('disp',j)) {
      if (grepl('lfc1',j)) {
        assign(paste0('pred.',j), prediction(1-get(paste0('p.',j), res), lfcd1))
        assign(paste0('fpr.raw.',j), `[<-`(get(paste0('fpr.raw.',j)), i, 
                                           value=1-specificity(get(paste0('call.raw.',j)), lfcd1)))
        assign(paste0('fpr.q.',j), `[<-`(get(paste0('fpr.q.',j)), i, 
                                           value=1-specificity(get(paste0('call.q.',j)), lfcd1)))
        assign(paste0('fpr.by.',j), `[<-`(get(paste0('fpr.by.',j)), i, 
                                          value=1-specificity(get(paste0('call.by.',j)), lfcd1)))
        assign(paste0('fpr.bh.',j), `[<-`(get(paste0('fpr.bh.',j)), i, 
                                          value=1-specificity(get(paste0('call.bh.',j)), lfcd1)))
        assign(paste0('fdr.raw.',j), `[<-`(get(paste0('fdr.raw.',j)), i, 
                                           value=1-precision(get(paste0('call.raw.',j)), lfcd1)))
        assign(paste0('fdr.q.',j), `[<-`(get(paste0('fdr.q.',j)), i, 
                                           value=1-precision(get(paste0('call.q.',j)), lfcd1)))
        assign(paste0('fdr.by.',j), `[<-`(get(paste0('fdr.by.',j)), i, 
                                          value=1-precision(get(paste0('call.by.',j)), lfcd1)))
        assign(paste0('fdr.bh.',j), `[<-`(get(paste0('fdr.bh.',j)), i, 
                                          value=1-precision(get(paste0('call.bh.',j)), lfcd1)))
        assign(paste0('tpr.raw.',j), `[<-`(get(paste0('tpr.raw.',j)), i, 
                                           value=sensitivity(get(paste0('call.raw.',j)), lfcd1)))
        assign(paste0('tpr.q.',j), `[<-`(get(paste0('tpr.q.',j)), i, 
                                           value=sensitivity(get(paste0('call.q.',j)), lfcd1)))
        assign(paste0('tpr.by.',j), `[<-`(get(paste0('tpr.by.',j)), i, 
                                          value=sensitivity(get(paste0('call.by.',j)), lfcd1)))
        assign(paste0('tpr.bh.',j), `[<-`(get(paste0('tpr.bh.',j)), i, 
                                          value=sensitivity(get(paste0('call.bh.',j)), lfcd1)))
      }
      else if (grepl('lfc2',j)) {
        assign(paste0('pred.',j), prediction(1-get(paste0('p.',j), res), lfcd2))
        assign(paste0('fpr.raw.',j), `[<-`(get(paste0('fpr.raw.',j)), i, 
                                           value=1-specificity(get(paste0('call.raw.',j)), lfcd2)))
        assign(paste0('fpr.q.',j), `[<-`(get(paste0('fpr.q.',j)), i, 
                                           value=1-specificity(get(paste0('call.q.',j)), lfcd2)))
        assign(paste0('fpr.by.',j), `[<-`(get(paste0('fpr.by.',j)), i, 
                                          value=1-specificity(get(paste0('call.by.',j)), lfcd2)))
        assign(paste0('fpr.bh.',j), `[<-`(get(paste0('fpr.bh.',j)), i, 
                                          value=1-specificity(get(paste0('call.bh.',j)), lfcd2)))
        assign(paste0('fdr.raw.',j), `[<-`(get(paste0('fdr.raw.',j)), i, 
                                           value=1-precision(get(paste0('call.raw.',j)), lfcd2)))
        assign(paste0('fdr.q.',j), `[<-`(get(paste0('fdr.q.',j)), i, 
                                           value=1-precision(get(paste0('call.q.',j)), lfcd2)))
        assign(paste0('fdr.by.',j), `[<-`(get(paste0('fdr.by.',j)), i, 
                                          value=1-precision(get(paste0('call.by.',j)), lfcd2)))
        assign(paste0('fdr.bh.',j), `[<-`(get(paste0('fdr.bh.',j)), i, 
                                          value=1-precision(get(paste0('call.bh.',j)), lfcd2)))
        assign(paste0('tpr.raw.',j), `[<-`(get(paste0('tpr.raw.',j)), i, 
                                           value=sensitivity(get(paste0('call.raw.',j)), lfcd2)))
        assign(paste0('tpr.q.',j), `[<-`(get(paste0('tpr.q.',j)), i, 
                                           value=sensitivity(get(paste0('call.q.',j)), lfcd2)))
        assign(paste0('tpr.by.',j), `[<-`(get(paste0('tpr.by.',j)), i, 
                                          value=sensitivity(get(paste0('call.by.',j)), lfcd2)))
        assign(paste0('tpr.bh.',j), `[<-`(get(paste0('tpr.bh.',j)), i, 
                                          value=sensitivity(get(paste0('call.bh.',j)), lfcd2)))
      }
      else {
        assign(paste0('pred.',j), prediction(1-get(paste0('p.',j), res), DD))
        assign(paste0('fpr.raw.',j), `[<-`(get(paste0('fpr.raw.',j)), i, 
                                           value=1-specificity(get(paste0('call.raw.',j)), DD)))
        assign(paste0('fpr.q.',j), `[<-`(get(paste0('fpr.q.',j)), i, 
                                           value=1-specificity(get(paste0('call.q.',j)), DD)))
        assign(paste0('fpr.by.',j), `[<-`(get(paste0('fpr.by.',j)), i, 
                                          value=1-specificity(get(paste0('call.by.',j)), DD)))
        assign(paste0('fpr.bh.',j), `[<-`(get(paste0('fpr.bh.',j)), i, 
                                          value=1-specificity(get(paste0('call.bh.',j)), DD)))
        assign(paste0('fdr.raw.',j), `[<-`(get(paste0('fdr.raw.',j)), i, 
                                           value=1-precision(get(paste0('call.raw.',j)), DD)))
        assign(paste0('fdr.q.',j), `[<-`(get(paste0('fdr.q.',j)), i, 
                                           value=1-precision(get(paste0('call.q.',j)), DD)))
        assign(paste0('fdr.by.',j), `[<-`(get(paste0('fdr.by.',j)), i, 
                                          value=1-precision(get(paste0('call.by.',j)), DD)))
        assign(paste0('fdr.bh.',j), `[<-`(get(paste0('fdr.bh.',j)), i, 
                                          value=1-precision(get(paste0('call.bh.',j)), DD)))
        assign(paste0('tpr.raw.',j), `[<-`(get(paste0('tpr.raw.',j)), i, 
                                           value=sensitivity(get(paste0('call.raw.',j)), DD)))
        assign(paste0('tpr.q.',j), `[<-`(get(paste0('tpr.q.',j)), i, 
                                           value=sensitivity(get(paste0('call.q.',j)), DD)))
        assign(paste0('tpr.by.',j), `[<-`(get(paste0('tpr.by.',j)), i, 
                                          value=sensitivity(get(paste0('call.by.',j)), DD)))
        assign(paste0('tpr.bh.',j), `[<-`(get(paste0('tpr.bh.',j)), i, 
                                          value=sensitivity(get(paste0('call.bh.',j)), DD)))
      }
    }
    else if (grepl('lfc1',j)) {
      assign(paste0('pred.',j), prediction(1-get(paste0('p.',j), res), lfcm1))
      assign(paste0('fpr.raw.',j), `[<-`(get(paste0('fpr.raw.',j)), i, 
                                         value=1-specificity(get(paste0('call.raw.',j)), lfcm1)))
      assign(paste0('fpr.q.',j), `[<-`(get(paste0('fpr.q.',j)), i, 
                                         value=1-specificity(get(paste0('call.q.',j)), lfcm1)))
      assign(paste0('fpr.by.',j), `[<-`(get(paste0('fpr.by.',j)), i, 
                                        value=1-specificity(get(paste0('call.by.',j)), lfcm1)))
      assign(paste0('fpr.bh.',j), `[<-`(get(paste0('fpr.bh.',j)), i, 
                                        value=1-specificity(get(paste0('call.bh.',j)), lfcm1)))
      assign(paste0('fdr.raw.',j), `[<-`(get(paste0('fdr.raw.',j)), i, 
                                         value=1-precision(get(paste0('call.raw.',j)), lfcm1)))
      assign(paste0('fdr.q.',j), `[<-`(get(paste0('fdr.q.',j)), i, 
                                         value=1-precision(get(paste0('call.q.',j)), lfcm1)))
      assign(paste0('fdr.by.',j), `[<-`(get(paste0('fdr.by.',j)), i, 
                                        value=1-precision(get(paste0('call.by.',j)), lfcm1)))
      assign(paste0('fdr.bh.',j), `[<-`(get(paste0('fdr.bh.',j)), i, 
                                        value=1-precision(get(paste0('call.bh.',j)), lfcm1)))
      assign(paste0('tpr.raw.',j), `[<-`(get(paste0('tpr.raw.',j)), i, 
                                         value=sensitivity(get(paste0('call.raw.',j)), lfcm1)))
      assign(paste0('tpr.q.',j), `[<-`(get(paste0('tpr.q.',j)), i, 
                                         value=sensitivity(get(paste0('call.q.',j)), lfcm1)))
      assign(paste0('tpr.by.',j), `[<-`(get(paste0('tpr.by.',j)), i, 
                                        value=sensitivity(get(paste0('call.by.',j)), lfcm1)))
      assign(paste0('tpr.bh.',j), `[<-`(get(paste0('tpr.bh.',j)), i, 
                                        value=sensitivity(get(paste0('call.bh.',j)), lfcm1)))
    }
    else if (grepl('lfc2',j)) {
      assign(paste0('pred.',j), prediction(1-get(paste0('p.',j), res), lfcm2))
      assign(paste0('fpr.raw.',j), `[<-`(get(paste0('fpr.raw.',j)), i, 
                                         value=1-specificity(get(paste0('call.raw.',j)), lfcm2)))
      assign(paste0('fpr.q.',j), `[<-`(get(paste0('fpr.q.',j)), i, 
                                         value=1-specificity(get(paste0('call.q.',j)), lfcm2)))
      assign(paste0('fpr.by.',j), `[<-`(get(paste0('fpr.by.',j)), i, 
                                        value=1-specificity(get(paste0('call.by.',j)), lfcm2)))
      assign(paste0('fpr.bh.',j), `[<-`(get(paste0('fpr.bh.',j)), i, 
                                        value=1-specificity(get(paste0('call.bh.',j)), lfcm2)))
      assign(paste0('fdr.raw.',j), `[<-`(get(paste0('fdr.raw.',j)), i, 
                                         value=1-precision(get(paste0('call.raw.',j)), lfcm2)))
      assign(paste0('fdr.q.',j), `[<-`(get(paste0('fdr.q.',j)), i, 
                                         value=1-precision(get(paste0('call.q.',j)), lfcm2)))
      assign(paste0('fdr.by.',j), `[<-`(get(paste0('fdr.by.',j)), i, 
                                        value=1-precision(get(paste0('call.by.',j)), lfcm2)))
      assign(paste0('fdr.bh.',j), `[<-`(get(paste0('fdr.bh.',j)), i, 
                                        value=1-precision(get(paste0('call.bh.',j)), lfcm2)))
      assign(paste0('tpr.raw.',j), `[<-`(get(paste0('tpr.raw.',j)), i, 
                                         value=sensitivity(get(paste0('call.raw.',j)), lfcm2)))
      assign(paste0('tpr.q.',j), `[<-`(get(paste0('tpr.q.',j)), i, 
                                         value=sensitivity(get(paste0('call.q.',j)), lfcm2)))
      assign(paste0('tpr.by.',j), `[<-`(get(paste0('tpr.by.',j)), i, 
                                        value=sensitivity(get(paste0('call.by.',j)), lfcm2)))
      assign(paste0('tpr.bh.',j), `[<-`(get(paste0('tpr.bh.',j)), i, 
                                        value=sensitivity(get(paste0('call.bh.',j)), lfcm2)))
    }
    else {
      assign(paste0('pred.',j), prediction(1-get(paste0('p.',j), res), DE))
      assign(paste0('fpr.raw.',j), `[<-`(get(paste0('fpr.raw.',j)), i, 
                                         value=1-specificity(get(paste0('call.raw.',j)), DE)))
      assign(paste0('fpr.q.',j), `[<-`(get(paste0('fpr.q.',j)), i, 
                                         value=1-specificity(get(paste0('call.q.',j)), DE)))
      assign(paste0('fpr.by.',j), `[<-`(get(paste0('fpr.by.',j)), i, 
                                        value=1-specificity(get(paste0('call.by.',j)), DE)))
      assign(paste0('fpr.bh.',j), `[<-`(get(paste0('fpr.bh.',j)), i, 
                                        value=1-specificity(get(paste0('call.bh.',j)), DE)))
      assign(paste0('fdr.raw.',j), `[<-`(get(paste0('fdr.raw.',j)), i, 
                                         value=1-precision(get(paste0('call.raw.',j)), DE)))
      assign(paste0('fdr.q.',j), `[<-`(get(paste0('fdr.q.',j)), i, 
                                         value=1-precision(get(paste0('call.q.',j)), DE)))
      assign(paste0('fdr.by.',j), `[<-`(get(paste0('fdr.by.',j)), i, 
                                        value=1-precision(get(paste0('call.by.',j)), DE)))
      assign(paste0('fdr.bh.',j), `[<-`(get(paste0('fdr.bh.',j)), i, 
                                        value=1-precision(get(paste0('call.bh.',j)), DE)))
      assign(paste0('tpr.raw.',j), `[<-`(get(paste0('tpr.raw.',j)), i, 
                                         value=sensitivity(get(paste0('call.raw.',j)), DE)))
      assign(paste0('tpr.q.',j), `[<-`(get(paste0('tpr.q.',j)), i, 
                                         value=sensitivity(get(paste0('call.q.',j)), DE)))
      assign(paste0('tpr.by.',j), `[<-`(get(paste0('tpr.by.',j)), i, 
                                        value=sensitivity(get(paste0('call.by.',j)), DE)))
      assign(paste0('tpr.bh.',j), `[<-`(get(paste0('tpr.bh.',j)), i, 
                                        value=sensitivity(get(paste0('call.bh.',j)), DE)))
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
  for (j in c(names.expHM[1],names.lnHM[1])) {
    assign(paste0('call.5.',j), factor(as.numeric(get(paste0('prob.',j), res) > 0.5), levels=c('1','0')))
    assign(paste0('call.thr.',j), factor(as.numeric(get(paste0('prob.',j), res) > get(paste0('thr.',j))), levels=c('1','0')))
    assign(paste0('pred.',j), prediction(get(paste0('prob.',j), res), DEDD))
    assign(paste0('fpr.5.',j), `[<-`(get(paste0('fpr.5.',j)), i, 
                                     value=1-specificity(get(paste0('call.5.',j)), DEDD)))
    assign(paste0('fpr.thr.',j), `[<-`(get(paste0('fpr.thr.',j)), i, 
                                       value=1-specificity(get(paste0('call.thr.',j)), DEDD)))
    assign(paste0('fdr.5.',j), `[<-`(get(paste0('fdr.5.',j)), i, 
                                     value=1-precision(get(paste0('call.5.',j)), DEDD)))
    assign(paste0('fdr.thr.',j), `[<-`(get(paste0('fdr.thr.',j)), i, 
                                       value=1-precision(get(paste0('call.thr.',j)), DEDD)))
    assign(paste0('tpr.5.',j), `[<-`(get(paste0('tpr.5.',j)), i, 
                                     value=sensitivity(get(paste0('call.5.',j)), DEDD)))
    assign(paste0('tpr.thr.',j), `[<-`(get(paste0('tpr.thr.',j)), i, 
                                       value=sensitivity(get(paste0('call.thr.',j)), DEDD)))
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

for (i in results.list) {
  assign(paste0('mean.discoveries.',i), colMeans(get(paste0('discoveries.',i))))
  assign(paste0('mean.fdr.',i), colMeans(get(paste0('false.discoveries.',i)) / get(paste0('discoveries.',i))))
}

rm(list=ls()[grep('^pred',ls())])
rm(list=ls()[grep('^bh',ls())])
rm(list=ls()[grep('^by',ls())])
rm(list=ls()[grep('^q',ls())])
rm(list=ls()[grep('^thr',ls())])
rm(list=ls()[grep('^call',ls())])
rm(list=ls()[grep('^discoveries',ls())])
rm(list=ls()[grep('^false.discoveries',ls())])
rm(list=c('i','j','import','res','DE','DD','DEDD','lfcd1','lfcd2','lfcm1','lfcm2'))

names.DE <- c('ql.edgeR', 'lr.edgeR', 'et.edgeR', 'noif.DESeq', 'if.DESeq', 'voom', 'notrend.DSS', 'trend.DSS', 
                'baySeq', 'mean.zi.MDSeq', 'mean.nozi.MDSeq', 'mean.expHM', 'lmean.expHM', 'mean.lnHM', 'lmean.lnHM', 
                'mix.ShrinkBayes', 'np.ShrinkBayes')
names.DE.lfc1 <- c('ql.lfc1.edgeR', 'lr.lfc1.edgeR', 'lfc1.DESeq', 'lfc1.voom', 'np.lfc1.ShrinkBayes', 
                     'mean.zi.lfc1.MDSeq', 'mean.nozi.lfc1.MDSeq', 'mean.lfc1.expHM', 'mean.lfc1.lnHM')
names.DE.lfc2 <- c('ql.lfc2.edgeR', 'lr.lfc2.edgeR', 'lfc2.DESeq', 'lfc2.voom', 'np.lfc2.ShrinkBayes', 
                     'mean.zi.lfc2.MDSeq', 'mean.nozi.lfc2.MDSeq', 'mean.lfc2.expHM', 'mean.lfc2.lnHM')
names.DD <- c('disp.zi.MDSeq', 'disp.nozi.MDSeq', 'disp.expHM', 'ldisp.expHM', 'disp.lnHM', 'ldisp.lnHM')
names.DD.lfc1 <- c('disp.zi.lfc1.MDSeq', 'disp.nozi.lfc1.MDSeq', 'disp.lfc1.expHM', 'disp.lfc1.lnHM')
names.DD.lfc2 <- c('disp.zi.lfc2.MDSeq', 'disp.nozi.lfc2.MDSeq', 'disp.lfc2.expHM', 'disp.lfc2.lnHM')
names.DEDD <- c('expHM', 'lnHM')

length(grep('^auc',ls())) # 51 vectors length 50
length(grep('^fpr',ls())) # 140 vectors length 50
length(grep('^fdr',ls())) # 140 vectors length 50
length(grep('^tpr',ls())) # 140 vectors length 50
length(grep('^mean.fdr',ls())) # 51 vectors length 20,000
length(grep('^mean.discoveries',ls())) # 51 vectors length 20,000

# 7 results files, each a list of 6 dataframes
# DE.DEDD2
# DE.lfc1.DEDD2
# DE.lfc2.DEDD2
# DD.DEDD2
# DD.lfc1.DEDD2
# DD.lfc2.DEDD2
# DEDD.DEDD2

for (i in c('DE', 'DE.lfc1', 'DE.lfc2', 'DD', 'DD.lfc1', 'DD.lfc2', 'DEDD')) {
    assign(paste0('results.',i), list())
  for (j in c('auc', 'fpr', 'fdr', 'tpr')) {
    assign(paste0('results.',i), `[[<-`(get(paste0('results.',i)), j, data.frame(matrix(nrow=50, ncol=0))))
  }
  for (j in c('mean.fdr', 'mean.discoveries')) {
    assign(paste0('results.',i), `[[<-`(get(paste0('results.',i)), j, data.frame(matrix(nrow=20000, ncol=0))))
  }
}

names.fpr.fdr.tpr <- vector()
for (i in c(names.edgeR, names.DESeq2, names.voom, names.DSS, names.MDSeq)) {
  for (j in c('raw.','fdr.')) {
    names.fpr.fdr.tpr <- c(names.fpr.fdr.tpr,paste0(j,i))
  }
}
for (i in names.DSS) {
  for (j in c('lfdr.','bh.')) {
    names.fpr.fdr.tpr <- c(names.fpr.fdr.tpr,paste0(j,i))
  }
}
for (i in names.baySeq) {
  for (j in c('5.','fdr.')) {
    names.fpr.fdr.tpr <- c(names.fpr.fdr.tpr,paste0(j,i))
  }
}
for (i in names.ShrinkBayes) { 
  for (j in c('lfdr.','bfdr.')) {
    names.fpr.fdr.tpr <- c(names.fpr.fdr.tpr,paste0(j,i))
  } # Only for DEDD2
}
for (i in c(names.expHM[1],names.lnHM[1])) {
  for (j in c('5.','thr.','bfdr.')) {
    names.fpr.fdr.tpr <- c(names.fpr.fdr.tpr,paste0(j,i))
  }
}
for (i in c(names.expHM[-1],names.lnHM[-1])) {
  for (j in c('raw.','bh.', 'by.','q.')) {
    names.fpr.fdr.tpr <- c(names.fpr.fdr.tpr,paste0(j,i))
  }
}

for (i in names.DE) {
  results.DE$auc[[i]] <- get(paste0('auc.',i))
  results.DE$mean.fdr[[i]] <- get(paste0('mean.fdr.',i))
  results.DE$mean.discoveries[[i]] <- get(paste0('mean.discoveries.',i))
  for (j in grep(i, names.fpr.fdr.tpr, value=T)) {
    results.DE$fpr[[j]] <- get(paste0('fpr.',j))
    results.DE$fdr[[j]] <- get(paste0('fdr.',j))
    results.DE$tpr[[j]] <- get(paste0('tpr.',j))
  }
}
for (i in names.DE.lfc1) {
  results.DE.lfc1$auc[[i]] <- get(paste0('auc.',i))
  results.DE.lfc1$mean.fdr[[i]] <- get(paste0('mean.fdr.',i))
  results.DE.lfc1$mean.discoveries[[i]] <- get(paste0('mean.discoveries.',i))
  for (j in grep(i, names.fpr.fdr.tpr, value=T)) {
    results.DE.lfc1$fpr[[j]] <- get(paste0('fpr.',j))
    results.DE.lfc1$fdr[[j]] <- get(paste0('fdr.',j))
    results.DE.lfc1$tpr[[j]] <- get(paste0('tpr.',j))
  }
}
for (i in names.DE.lfc2) {
  results.DE.lfc2$auc[[i]] <- get(paste0('auc.',i))
  results.DE.lfc2$mean.fdr[[i]] <- get(paste0('mean.fdr.',i))
  results.DE.lfc2$mean.discoveries[[i]] <- get(paste0('mean.discoveries.',i))
  for (j in grep(i, names.fpr.fdr.tpr, value=T)) {
    results.DE.lfc2$fpr[[j]] <- get(paste0('fpr.',j))
    results.DE.lfc2$fdr[[j]] <- get(paste0('fdr.',j))
    results.DE.lfc2$tpr[[j]] <- get(paste0('tpr.',j))
  }
}
for (i in names.DD) {
  results.DD$auc[[i]] <- get(paste0('auc.',i))
  results.DD$mean.fdr[[i]] <- get(paste0('mean.fdr.',i))
  results.DD$mean.discoveries[[i]] <- get(paste0('mean.discoveries.',i))
  for (j in grep(i, names.fpr.fdr.tpr, value=T)) {
    results.DD$fpr[[j]] <- get(paste0('fpr.',j))
    results.DD$fdr[[j]] <- get(paste0('fdr.',j))
    results.DD$tpr[[j]] <- get(paste0('tpr.',j))
  }
}
for (i in names.DD.lfc1) {
  results.DD.lfc1$auc[[i]] <- get(paste0('auc.',i))
  results.DD.lfc1$mean.fdr[[i]] <- get(paste0('mean.fdr.',i))
  results.DD.lfc1$mean.discoveries[[i]] <- get(paste0('mean.discoveries.',i))
  for (j in grep(i, names.fpr.fdr.tpr, value=T)) {
    results.DD.lfc1$fpr[[j]] <- get(paste0('fpr.',j))
    results.DD.lfc1$fdr[[j]] <- get(paste0('fdr.',j))
    results.DD.lfc1$tpr[[j]] <- get(paste0('tpr.',j))
  }
}
for (i in names.DD.lfc2) {
  results.DD.lfc2$auc[[i]] <- get(paste0('auc.',i))
  results.DD.lfc2$mean.fdr[[i]] <- get(paste0('mean.fdr.',i))
  results.DD.lfc2$mean.discoveries[[i]] <- get(paste0('mean.discoveries.',i))
  for (j in grep(i, names.fpr.fdr.tpr, value=T)) {
    results.DD.lfc2$fpr[[j]] <- get(paste0('fpr.',j))
    results.DD.lfc2$fdr[[j]] <- get(paste0('fdr.',j))
    results.DD.lfc2$tpr[[j]] <- get(paste0('tpr.',j))
  }
}
for (i in names.DEDD) {
  results.DEDD$auc[[i]] <- get(paste0('auc.',i))
  results.DEDD$mean.fdr[[i]] <- get(paste0('mean.fdr.',i))
  results.DEDD$mean.discoveries[[i]] <- get(paste0('mean.discoveries.',i))
  for (j in grep(i, names.fpr.fdr.tpr, value=T)) {
    results.DEDD$fpr[[j]] <- get(paste0('fpr.',j))
    results.DEDD$fdr[[j]] <- get(paste0('fdr.',j))
    results.DEDD$tpr[[j]] <- get(paste0('tpr.',j))
  }
}







par(mfrow=c(4,3))
plot(mean.discoveries.ql.edgeR, mean.fdr.ql.edgeR, type='l'); abline(h=0.05)
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


