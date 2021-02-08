library(here)
library(compcodeR)
library(ROCR)
library(caret)
library(qvalue)
source(here('scripts','2019-05-03_bfdr_function.R'))

# Data:
# data DE DD lfcm1 lfcm2 lfcd1 lfcd2

# edgeR:
# p.ql.edgeR  p.ql.lfc1.edgeR  p.ql.lfc2.edgeR  p.lr.edgeR  p.lr.lfc1.edgeR  p.lr.lfc2.edgeR  p.et.edgeR
# q.ql.edgeR  q.ql.lfc1.edgeR  q.ql.lfc2.edgeR  q.lr.edgeR  q.lr.lfc1.edgeR  q.lr.lfc2.edgeR  q.et.edgeR

# DESeq
# p.noif.DESeq  p.if.DESeq  p.lfc1.DESeq  p.lfc2.DESeq
# q.noif.DESeq  q.if.DESeq  q.lfc1.DESeq  q.lfc2.DESeq

# voom
# p.voom  p.lfc1.voom  p.lfc1.voom
# q.voom  q.lfc2.voom  q.lfc2.voom

# DSS
# p.notrend.DSS
# q.notrend.DSS
# lfdr.notrend.DSS
# mean(res$q.notrend.DSS == p.adjust(res$p.notrend.DSS, method='BH')) is 0. Not clear what 
# fdr (q here) and lfdr from DSS are, but they don't seem to be increasing with p-value. 
# Should also compute BH q-values from p-values.
# trend in DEDD2, not DEDD2, DEDD10; in DE2, DE5, not DE10

# baySeq
# prob.baySeq
# q.baySeq

# MDSeq
# p.mean.zi.MDSeq   p.mean.zi.lfc1.MDSeq   p.mean.zi.lfc2.MDSeq
# p.mean.nozi.MDSeq p.mean.nozi.lfc1.MDSeq p.mean.nozi.lfc2.MDSeq
# p.disp.zi.MDSeq   p.disp.zi.lfc1.MDSeq   p.disp.zi.lfc2.MDSeq
# p.disp.nozi.MDSeq p.disp.nozi.lfc1.MDSeq p.disp.nozi.lfc2.MDSeq
# q.mean.zi.MDSeq   q.mean.zi.lfc1.MDSeq   q.mean.zi.lfc2.MDSeq
# q.mean.nozi.MDSeq q.mean.nozi.lfc1.MDSeq q.mean.nozi.lfc2.MDSeq
# q.disp.zi.MDSeq   q.disp.zi.lfc1.MDSeq   q.disp.zi.lfc2.MDSeq
# q.disp.nozi.MDSeq q.disp.nozi.lfc1.MDSeq q.disp.nozi.lfc2.MDSeq

# expHM
# prob.expHM    prop.expHM
# p.mean.expHM  p.lmean.expHM  p.mean.lfc1.expHM  p.mean.lfc2.expHM
# p.disp.expHM  p.ldisp.expHM  p.disp.lfc1.expHM  p.disp.lfc2.expHM
# Only have p-values so need to compute q-values. 
# Compute q, BH and BY.

# lnHM
# prob.lnHM    prop.lnHM
# p.mean.lnHM  p.lmean.lnHM  p.mean.lfc1.lnHM  p.mean.lfc2.lnHM
# p.disp.lnHM  p.ldisp.lnHM  p.disp.lfc1.lnHM  p.disp.lfc2.lnHM
# Only have p-values so need to compute q-values. 
# Compute q, BH and BY.


# Compute AUC, FPR, FDR, TPR for each run, and mean discoveries and FDRs ####
# Use existing FDR-corrected values where exist (i.e. where were included in method), 
# and compute others, i.e. BH for DSS and q, BH and BY for HMs, and results using 
# 50% posterior probability threshold for baySeq and HMMs, and posterior proportion 
# threshold for HMMs.

# First create variables to store results
{
# names.edgeR <- c('ql.edgeR', 'ql.lfc1.edgeR', 'ql.lfc2.edgeR', 'lr.edgeR', 'lr.lfc1.edgeR', 'lr.lfc2.edgeR', 'et.edgeR')
# names.DESeq2 <- c('noif.DESeq', 'if.DESeq', 'lfc1.DESeq', 'lfc2.DESeq')
# names.voom <- c('voom', 'lfc1.voom', 'lfc2.voom')
# names.DSS <- 'notrend.DSS'
# names.baySeq <- 'baySeq'
names.MDSeq <- c(
  # 'mean.zi.MDSeq', 'mean.nozi.MDSeq', 'mean.zi.lfc1.MDSeq', 'mean.nozi.lfc1.MDSeq',
  # 'mean.zi.lfc2.MDSeq', 'mean.nozi.lfc2.MDSeq',
  'disp.zi.MDSeq', 'disp.nozi.MDSeq', 'disp.zi.lfc1.MDSeq', 'disp.nozi.lfc1.MDSeq',
  'disp.zi.lfc2.MDSeq', 'disp.nozi.lfc2.MDSeq'
  )
names.expHM <- c('expHM',
                 # 'mean.expHM', 'lmean.expHM', 'mean.lfc1.expHM', 'mean.lfc2.expHM',
                 'disp.expHM', 'ldisp.expHM', 'disp.lfc1.expHM', 'disp.lfc2.expHM'
                 )
names.lnHM <- c('lnHM',
                # 'mean.lnHM', 'lmean.lnHM', 'mean.lfc1.lnHM', 'mean.lfc2.lnHM',
                'disp.lnHM', 'ldisp.lnHM', 'disp.lfc1.lnHM', 'disp.lfc2.lnHM'
                )
results.list <- c(
  # names.edgeR, names.DESeq2, names.voom, names.DSS, names.baySeq, 
  names.MDSeq, names.expHM, names.lnHM)
}

{
for (i in c(
  # names.edgeR, names.DESeq2, names.voom, names.DSS, 
            names.MDSeq)) {
  for (j in c('fpr.raw.','fdr.raw.','tpr.raw.','fpr.fdr.','fdr.fdr.','tpr.fdr.')) {
    assign(paste0(j,i), numeric(50))
  }
}
# for (i in names.DSS) {
#  for (j in c('fpr.lfdr.','fdr.lfdr.','tpr.lfdr.','fpr.bh.','fdr.bh.','tpr.bh.')) {
#    assign(paste0(j,i), numeric(50))
#  }
# }
# for (i in names.baySeq) {
#  for (j in c('fpr.5.','fdr.5.','tpr.5.','fpr.fdr.','fdr.fdr.','tpr.fdr.')) {
#    assign(paste0(j,i), numeric(50))
#  }
# }
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
}

# Then import each run in turn and save results
for (i in 1:50) {
  import <- paste0('results.DD2.',i,'.DESeqnorm.rds')
  res <- readRDS(here('Results/DD compcodeR data results DESeq norm Sept 2019', import))
  # DE <- factor(res$DE, levels=c('1','0'))
  DD <- factor(res$DD, levels=c('1','0'))
  # DEDD <- factor(as.numeric(res$DE==1 | res$DD==1), levels=c('1','0'))
  DEDD <- DD
  # lfcm1 <- factor(as.numeric(res$lfcm1), levels=c('1','0'))
  # lfcm2 <- factor(as.numeric(res$lfcm2), levels=c('1','0'))
  lfcd1 <- factor(as.numeric(res$lfcd1), levels=c('1','0'))
  lfcd2 <- factor(as.numeric(res$lfcd2), levels=c('1','0'))
  # bh.notrend.DSS <- p.adjust(res$p.notrend.DSS, method='BH')
  for (j in c(names.expHM[-1], names.lnHM[-1])) {
    assign(paste0('bh.',j), p.adjust(get(paste0('p.',j), res), method='BH'))
    assign(paste0('by.',j), p.adjust(get(paste0('p.',j), res), method='BY'))
    assign(paste0('q.',j), qvalue(get(paste0('p.',j), res))$qval)
  }
  thr.expHM <- sort(res$prob.expHM, decreasing=T)[round(nrow(res$data@count.matrix)*res$prop.expHM)]
  thr.lnHM <- sort(res$prob.lnHM, decreasing=T)[round(nrow(res$data@count.matrix)*res$prop.lnHM)]
  for (j in c(
    # names.edgeR, names.DESeq2, names.voom, names.DSS, 
    names.MDSeq)) {
    assign(paste0('call.raw.',j), factor(as.numeric(get(paste0('p.',j), res) < 0.05), levels=c('1','0')))
    assign(paste0('call.fdr.',j), factor(as.numeric(get(paste0('q.',j), res) < 0.05), levels=c('1','0')))
    # if (grepl('disp',j)) {
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
   #  }
   #  else
   #    if (grepl('lfc1',j)) {
   #   assign(paste0('pred.',j), prediction(1-get(paste0('p.',j), res), lfcm1))
   #   assign(paste0('fpr.raw.',j), `[<-`(get(paste0('fpr.raw.',j)), i,
   #                                      value=1-specificity(get(paste0('call.raw.',j)), lfcm1)))
   #   assign(paste0('fpr.fdr.',j), `[<-`(get(paste0('fpr.fdr.',j)), i,
   #                                      value=1-specificity(get(paste0('call.fdr.',j)), lfcm1)))
   #   assign(paste0('fdr.raw.',j), `[<-`(get(paste0('fdr.raw.',j)), i,
   #                                      value=1-precision(get(paste0('call.raw.',j)), lfcm1)))
   #   assign(paste0('fdr.fdr.',j), `[<-`(get(paste0('fdr.fdr.',j)), i,
   #                                      value=1-precision(get(paste0('call.fdr.',j)), lfcm1)))
   #   assign(paste0('tpr.raw.',j), `[<-`(get(paste0('tpr.raw.',j)), i,
   #                                      value=sensitivity(get(paste0('call.raw.',j)), lfcm1)))
   #   assign(paste0('tpr.fdr.',j), `[<-`(get(paste0('tpr.fdr.',j)), i,
   #                                      value=sensitivity(get(paste0('call.fdr.',j)), lfcm1)))
   # }
   # else if (grepl('lfc2',j)) {
   #   assign(paste0('pred.',j), prediction(1-get(paste0('p.',j), res), lfcm2))
   #   assign(paste0('fpr.raw.',j), `[<-`(get(paste0('fpr.raw.',j)), i,
   #                                      value=1-specificity(get(paste0('call.raw.',j)), lfcm2)))
   #   assign(paste0('fpr.fdr.',j), `[<-`(get(paste0('fpr.fdr.',j)), i,
   #                                      value=1-specificity(get(paste0('call.fdr.',j)), lfcm2)))
   #   assign(paste0('fdr.raw.',j), `[<-`(get(paste0('fdr.raw.',j)), i,
   #                                      value=1-precision(get(paste0('call.raw.',j)), lfcm2)))
   #   assign(paste0('fdr.fdr.',j), `[<-`(get(paste0('fdr.fdr.',j)), i,
   #                                      value=1-precision(get(paste0('call.fdr.',j)), lfcm2)))
   #   assign(paste0('tpr.raw.',j), `[<-`(get(paste0('tpr.raw.',j)), i,
   #                                      value=sensitivity(get(paste0('call.raw.',j)), lfcm2)))
   #   assign(paste0('tpr.fdr.',j), `[<-`(get(paste0('tpr.fdr.',j)), i,
   #                                      value=sensitivity(get(paste0('call.fdr.',j)), lfcm2)))
   # }
   # else {
   #   assign(paste0('pred.',j), prediction(1-get(paste0('p.',j), res), DE))
   #   assign(paste0('fpr.raw.',j), `[<-`(get(paste0('fpr.raw.',j)), i,
   #                                      value=1-specificity(get(paste0('call.raw.',j)), DE)))
   #   assign(paste0('fpr.fdr.',j), `[<-`(get(paste0('fpr.fdr.',j)), i,
   #                                      value=1-specificity(get(paste0('call.fdr.',j)), DE)))
   #   assign(paste0('fdr.raw.',j), `[<-`(get(paste0('fdr.raw.',j)), i,
   #                                      value=1-precision(get(paste0('call.raw.',j)), DE)))
   #   assign(paste0('fdr.fdr.',j), `[<-`(get(paste0('fdr.fdr.',j)), i,
   #                                      value=1-precision(get(paste0('call.fdr.',j)), DE)))
   #   assign(paste0('tpr.raw.',j), `[<-`(get(paste0('tpr.raw.',j)), i,
   #                                      value=sensitivity(get(paste0('call.raw.',j)), DE)))
   #   assign(paste0('tpr.fdr.',j), `[<-`(get(paste0('tpr.fdr.',j)), i,
   #                                      value=sensitivity(get(paste0('call.fdr.',j)), DE)))
   # }
    assign(paste0('auc.',j), `[<-`(get(paste0('auc.',j)), i, 
                                   value=performance(get(paste0('pred.',j)), measure='auc')@y.values[[1]]))
    assign(paste0('false.discoveries.',j), `[<-`(get(paste0('false.discoveries.',j)), i,, 
                                                 value=c(get(paste0('pred.',j))@fp[[1]], 
                                                         rep(NA, 20000-length(get(paste0('pred.',j))@fp[[1]])))))
    assign(paste0('discoveries.',j), `[<-`(get(paste0('discoveries.',j)), i,, 
                                           value=c(get(paste0('pred.',j))@n.pos.pred[[1]], 
                                                   rep(NA, 20000-length(get(paste0('pred.',j))@n.pos.pred[[1]])))))
  }
 # for (j in names.DSS) {
 #   assign(paste0('call.bh.',j), factor(as.numeric(get(paste0('bh.',j)) < 0.05), levels=c('1','0')))
 #   assign(paste0('call.lfdr.',j), factor(as.numeric(get(paste0('lfdr.',j), res) < 0.05), levels=c('1','0')))
 #   assign(paste0('fpr.bh.',j), `[<-`(get(paste0('fpr.bh.',j)), i,
 #                                     value=1-specificity(get(paste0('call.bh.',j)), DE)))
 #   assign(paste0('fpr.lfdr.',j), `[<-`(get(paste0('fpr.lfdr.',j)), i,
 #                                       value=1-specificity(get(paste0('call.lfdr.',j)), DE)))
 #   assign(paste0('fdr.bh.',j), `[<-`(get(paste0('fdr.bh.',j)), i,
 #                                     value=1-precision(get(paste0('call.bh.',j)), DE)))
 #   assign(paste0('fdr.lfdr.',j), `[<-`(get(paste0('fdr.lfdr.',j)), i,
 #                                       value=1-precision(get(paste0('call.lfdr.',j)), DE)))
 #   assign(paste0('tpr.bh.',j), `[<-`(get(paste0('tpr.bh.',j)), i,
 #                                     value=sensitivity(get(paste0('call.bh.',j)), DE)))
 #   assign(paste0('tpr.lfdr.',j), `[<-`(get(paste0('tpr.lfdr.',j)), i,
 #                                       value=sensitivity(get(paste0('call.lfdr.',j)), DE)))
 # }
 # for (j in names.baySeq) {
 #   assign(paste0('call.5.',j), factor(as.numeric(get(paste0('prob.',j), res) > 0.5), levels=c('1','0')))
 #   assign(paste0('call.fdr.',j), factor(as.numeric(get(paste0('q.',j), res) < 0.05), levels=c('1','0')))
 #   assign(paste0('pred.',j), prediction(get(paste0('prob.',j), res), DE))
 #   assign(paste0('fpr.5.',j), `[<-`(get(paste0('fpr.5.',j)), i,
 #                                    value=1-specificity(get(paste0('call.5.',j)), DE)))
 #   assign(paste0('fpr.fdr.',j), `[<-`(get(paste0('fpr.fdr.',j)), i,
 #                                      value=1-specificity(get(paste0('call.fdr.',j)), DE)))
 #   assign(paste0('fdr.5.',j), `[<-`(get(paste0('fdr.5.',j)), i,
 #                                    value=1-precision(get(paste0('call.5.',j)), DE)))
 #   assign(paste0('fdr.fdr.',j), `[<-`(get(paste0('fdr.fdr.',j)), i,
 #                                      value=1-precision(get(paste0('call.fdr.',j)), DE)))
 #   assign(paste0('tpr.5.',j), `[<-`(get(paste0('tpr.5.',j)), i,
 #                                    value=sensitivity(get(paste0('call.5.',j)), DE)))
 #   assign(paste0('tpr.fdr.',j), `[<-`(get(paste0('tpr.fdr.',j)), i,
 #                                      value=sensitivity(get(paste0('call.fdr.',j)), DE)))
 #   assign(paste0('auc.',j), `[<-`(get(paste0('auc.',j)), i,
 #                                  value=performance(get(paste0('pred.',j)), measure='auc')@y.values[[1]]))
 #   assign(paste0('false.discoveries.',j), `[<-`(get(paste0('false.discoveries.',j)), i,,
 #                                                value=c(get(paste0('pred.',j))@fp[[1]],
 #                                                        rep(NA, 20000-length(get(paste0('pred.',j))@fp[[1]])))))
 #   assign(paste0('discoveries.',j), `[<-`(get(paste0('discoveries.',j)), i,,
 #                                          value=c(get(paste0('pred.',j))@n.pos.pred[[1]],
 #                                                  rep(NA, 20000-length(get(paste0('pred.',j))@n.pos.pred[[1]])))))
 # }
  for (j in c(names.expHM[-1],names.lnHM[-1])) {
    assign(paste0('call.raw.',j), factor(as.numeric(get(paste0('p.',j), res) < 0.05), levels=c('1','0')))
    assign(paste0('call.q.',j), factor(as.numeric(get(paste0('q.',j)) < 0.05), levels=c('1','0')))
    assign(paste0('call.by.',j), factor(as.numeric(get(paste0('by.',j)) < 0.05), levels=c('1','0')))
    assign(paste0('call.bh.',j), factor(as.numeric(get(paste0('bh.',j)) < 0.05), levels=c('1','0')))
    # if (grepl('disp',j)) {
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
   #  }
   #  else
   #    if (grepl('lfc1',j)) {
   #   assign(paste0('pred.',j), prediction(1-get(paste0('p.',j), res), lfcm1))
   #   assign(paste0('fpr.raw.',j), `[<-`(get(paste0('fpr.raw.',j)), i,
   #                                      value=1-specificity(get(paste0('call.raw.',j)), lfcm1)))
   #   assign(paste0('fpr.q.',j), `[<-`(get(paste0('fpr.q.',j)), i,
   #                                      value=1-specificity(get(paste0('call.q.',j)), lfcm1)))
   #   assign(paste0('fpr.by.',j), `[<-`(get(paste0('fpr.by.',j)), i,
   #                                     value=1-specificity(get(paste0('call.by.',j)), lfcm1)))
   #   assign(paste0('fpr.bh.',j), `[<-`(get(paste0('fpr.bh.',j)), i,
   #                                     value=1-specificity(get(paste0('call.bh.',j)), lfcm1)))
   #   assign(paste0('fdr.raw.',j), `[<-`(get(paste0('fdr.raw.',j)), i,
   #                                      value=1-precision(get(paste0('call.raw.',j)), lfcm1)))
   #   assign(paste0('fdr.q.',j), `[<-`(get(paste0('fdr.q.',j)), i,
   #                                      value=1-precision(get(paste0('call.q.',j)), lfcm1)))
   #   assign(paste0('fdr.by.',j), `[<-`(get(paste0('fdr.by.',j)), i,
   #                                     value=1-precision(get(paste0('call.by.',j)), lfcm1)))
   #   assign(paste0('fdr.bh.',j), `[<-`(get(paste0('fdr.bh.',j)), i,
   #                                     value=1-precision(get(paste0('call.bh.',j)), lfcm1)))
   #   assign(paste0('tpr.raw.',j), `[<-`(get(paste0('tpr.raw.',j)), i,
   #                                      value=sensitivity(get(paste0('call.raw.',j)), lfcm1)))
   #   assign(paste0('tpr.q.',j), `[<-`(get(paste0('tpr.q.',j)), i,
   #                                      value=sensitivity(get(paste0('call.q.',j)), lfcm1)))
   #   assign(paste0('tpr.by.',j), `[<-`(get(paste0('tpr.by.',j)), i,
   #                                     value=sensitivity(get(paste0('call.by.',j)), lfcm1)))
   #   assign(paste0('tpr.bh.',j), `[<-`(get(paste0('tpr.bh.',j)), i,
   #                                     value=sensitivity(get(paste0('call.bh.',j)), lfcm1)))
   # }
   # else if (grepl('lfc2',j)) {
   #   assign(paste0('pred.',j), prediction(1-get(paste0('p.',j), res), lfcm2))
   #   assign(paste0('fpr.raw.',j), `[<-`(get(paste0('fpr.raw.',j)), i,
   #                                      value=1-specificity(get(paste0('call.raw.',j)), lfcm2)))
   #   assign(paste0('fpr.q.',j), `[<-`(get(paste0('fpr.q.',j)), i,
   #                                      value=1-specificity(get(paste0('call.q.',j)), lfcm2)))
   #   assign(paste0('fpr.by.',j), `[<-`(get(paste0('fpr.by.',j)), i,
   #                                     value=1-specificity(get(paste0('call.by.',j)), lfcm2)))
   #   assign(paste0('fpr.bh.',j), `[<-`(get(paste0('fpr.bh.',j)), i,
   #                                     value=1-specificity(get(paste0('call.bh.',j)), lfcm2)))
   #   assign(paste0('fdr.raw.',j), `[<-`(get(paste0('fdr.raw.',j)), i,
   #                                      value=1-precision(get(paste0('call.raw.',j)), lfcm2)))
   #   assign(paste0('fdr.q.',j), `[<-`(get(paste0('fdr.q.',j)), i,
   #                                      value=1-precision(get(paste0('call.q.',j)), lfcm2)))
   #   assign(paste0('fdr.by.',j), `[<-`(get(paste0('fdr.by.',j)), i,
   #                                     value=1-precision(get(paste0('call.by.',j)), lfcm2)))
   #   assign(paste0('fdr.bh.',j), `[<-`(get(paste0('fdr.bh.',j)), i,
   #                                     value=1-precision(get(paste0('call.bh.',j)), lfcm2)))
   #   assign(paste0('tpr.raw.',j), `[<-`(get(paste0('tpr.raw.',j)), i,
   #                                      value=sensitivity(get(paste0('call.raw.',j)), lfcm2)))
   #   assign(paste0('tpr.q.',j), `[<-`(get(paste0('tpr.q.',j)), i,
   #                                      value=sensitivity(get(paste0('call.q.',j)), lfcm2)))
   #   assign(paste0('tpr.by.',j), `[<-`(get(paste0('tpr.by.',j)), i,
   #                                     value=sensitivity(get(paste0('call.by.',j)), lfcm2)))
   #   assign(paste0('tpr.bh.',j), `[<-`(get(paste0('tpr.bh.',j)), i,
   #                                     value=sensitivity(get(paste0('call.bh.',j)), lfcm2)))
   # }
   # else {
   #   assign(paste0('pred.',j), prediction(1-get(paste0('p.',j), res), DE))
   #   assign(paste0('fpr.raw.',j), `[<-`(get(paste0('fpr.raw.',j)), i,
   #                                      value=1-specificity(get(paste0('call.raw.',j)), DE)))
   #   assign(paste0('fpr.q.',j), `[<-`(get(paste0('fpr.q.',j)), i,
   #                                      value=1-specificity(get(paste0('call.q.',j)), DE)))
   #   assign(paste0('fpr.by.',j), `[<-`(get(paste0('fpr.by.',j)), i,
   #                                     value=1-specificity(get(paste0('call.by.',j)), DE)))
   #   assign(paste0('fpr.bh.',j), `[<-`(get(paste0('fpr.bh.',j)), i,
   #                                     value=1-specificity(get(paste0('call.bh.',j)), DE)))
   #   assign(paste0('fdr.raw.',j), `[<-`(get(paste0('fdr.raw.',j)), i,
   #                                      value=1-precision(get(paste0('call.raw.',j)), DE)))
   #   assign(paste0('fdr.q.',j), `[<-`(get(paste0('fdr.q.',j)), i,
   #                                      value=1-precision(get(paste0('call.q.',j)), DE)))
   #   assign(paste0('fdr.by.',j), `[<-`(get(paste0('fdr.by.',j)), i,
   #                                     value=1-precision(get(paste0('call.by.',j)), DE)))
   #   assign(paste0('fdr.bh.',j), `[<-`(get(paste0('fdr.bh.',j)), i,
   #                                     value=1-precision(get(paste0('call.bh.',j)), DE)))
   #   assign(paste0('tpr.raw.',j), `[<-`(get(paste0('tpr.raw.',j)), i,
   #                                      value=sensitivity(get(paste0('call.raw.',j)), DE)))
   #   assign(paste0('tpr.q.',j), `[<-`(get(paste0('tpr.q.',j)), i,
   #                                      value=sensitivity(get(paste0('call.q.',j)), DE)))
   #   assign(paste0('tpr.by.',j), `[<-`(get(paste0('tpr.by.',j)), i,
   #                                     value=sensitivity(get(paste0('call.by.',j)), DE)))
   #   assign(paste0('tpr.bh.',j), `[<-`(get(paste0('tpr.bh.',j)), i,
   #                                     value=sensitivity(get(paste0('call.bh.',j)), DE)))
   # }
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

# Clean up temporary variables no longer needed
{rm(list=ls()[grep('^pred',ls())])
rm(list=ls()[grep('^bh',ls())])
rm(list=ls()[grep('^by',ls())])
rm(list=ls()[grep('^q',ls())])
rm(list=ls()[grep('^thr',ls())])
rm(list=ls()[grep('^call',ls())])
rm(list=ls()[grep('^discoveries',ls())])
rm(list=ls()[grep('^false.discoveries',ls())])
rm(list=c('i','j','import','res','DE','DD','DEDD','lfcd1','lfcd2','lfcm1','lfcm2'))
}

# Re-format results to save by experiment type ####

# Set up variables
{
# names.DE <- c('ql.edgeR', 'lr.edgeR', 'et.edgeR', 'noif.DESeq', 'if.DESeq', 'voom', 'notrend.DSS', 
#                'baySeq', 'mean.zi.MDSeq', 'mean.nozi.MDSeq', 'mean.expHM', 'lmean.expHM', 'mean.lnHM', 'lmean.lnHM')
# names.DE.lfc1 <- c('ql.lfc1.edgeR', 'lr.lfc1.edgeR', 'lfc1.DESeq', 'lfc1.voom', 
#                     'mean.zi.lfc1.MDSeq', 'mean.nozi.lfc1.MDSeq', 'mean.lfc1.expHM', 'mean.lfc1.lnHM')
# names.DE.lfc2 <- c('ql.lfc2.edgeR', 'lr.lfc2.edgeR', 'lfc2.DESeq', 'lfc2.voom', 
#                     'mean.zi.lfc2.MDSeq', 'mean.nozi.lfc2.MDSeq', 'mean.lfc2.expHM', 'mean.lfc2.lnHM')
names.DD <- c('disp.zi.MDSeq', 'disp.nozi.MDSeq', 'disp.expHM', 'ldisp.expHM', 'disp.lnHM', 'ldisp.lnHM')
names.DD.lfc1 <- c('disp.zi.lfc1.MDSeq', 'disp.nozi.lfc1.MDSeq', 'disp.lfc1.expHM', 'disp.lfc1.lnHM')
names.DD.lfc2 <- c('disp.zi.lfc2.MDSeq', 'disp.nozi.lfc2.MDSeq', 'disp.lfc2.expHM', 'disp.lfc2.lnHM')
names.DEDD <- c('expHM', 'lnHM')
}

for (i in c(
  # 'DE', 'DE.lfc1', 'DE.lfc2', 
  'DD', 'DD.lfc1', 'DD.lfc2',
  'DEDD')) {
  assign(paste0('results.',i), list())
  for (j in c('auc', 'fpr', 'fdr', 'tpr')) {
    assign(paste0('results.',i), `[[<-`(get(paste0('results.',i)), j, data.frame(matrix(nrow=50, ncol=0))))
  }
  for (j in c('mean.fdr', 'mean.discoveries')) {
    assign(paste0('results.',i), `[[<-`(get(paste0('results.',i)), j, data.frame(matrix(nrow=20000, ncol=0))))
  }
}

{names.fpr.fdr.tpr <- vector()
for (i in c(
  # names.edgeR, names.DESeq2, names.voom, names.DSS, 
  names.MDSeq)) {
  for (j in c('raw.','fdr.')) {
    names.fpr.fdr.tpr <- c(names.fpr.fdr.tpr,paste0(j,i))
  }
}
# for (i in names.DSS) {
#  for (j in c('lfdr.','bh.')) {
#    names.fpr.fdr.tpr <- c(names.fpr.fdr.tpr,paste0(j,i))
#  }
# }
# for (i in names.baySeq) {
#  for (j in c('5.','fdr.')) {
#    names.fpr.fdr.tpr <- c(names.fpr.fdr.tpr,paste0(j,i))
#  }
# }
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
}

# Put results into variables
{
# for (i in names.DE) {
#  results.DE$auc[[i]] <- get(paste0('auc.',i))
#  results.DE$mean.fdr[[i]] <- get(paste0('mean.fdr.',i))
#  results.DE$mean.discoveries[[i]] <- get(paste0('mean.discoveries.',i))
#  for (j in grep(i, names.fpr.fdr.tpr, value=T)) {
#    results.DE$fpr[[j]] <- get(paste0('fpr.',j))
#    results.DE$fdr[[j]] <- get(paste0('fdr.',j))
#    results.DE$tpr[[j]] <- get(paste0('tpr.',j))
#  }
# }
# for (i in names.DE.lfc1) {
#  results.DE.lfc1$auc[[i]] <- get(paste0('auc.',i))
#  results.DE.lfc1$mean.fdr[[i]] <- get(paste0('mean.fdr.',i))
#  results.DE.lfc1$mean.discoveries[[i]] <- get(paste0('mean.discoveries.',i))
#  for (j in grep(i, names.fpr.fdr.tpr, value=T)) {
#    results.DE.lfc1$fpr[[j]] <- get(paste0('fpr.',j))
#    results.DE.lfc1$fdr[[j]] <- get(paste0('fdr.',j))
#    results.DE.lfc1$tpr[[j]] <- get(paste0('tpr.',j))
#  }
# }
# for (i in names.DE.lfc2) {
#  results.DE.lfc2$auc[[i]] <- get(paste0('auc.',i))
#  results.DE.lfc2$mean.fdr[[i]] <- get(paste0('mean.fdr.',i))
#  results.DE.lfc2$mean.discoveries[[i]] <- get(paste0('mean.discoveries.',i))
#  for (j in grep(i, names.fpr.fdr.tpr, value=T)) {
#    results.DE.lfc2$fpr[[j]] <- get(paste0('fpr.',j))
#    results.DE.lfc2$fdr[[j]] <- get(paste0('fdr.',j))
#    results.DE.lfc2$tpr[[j]] <- get(paste0('tpr.',j))
#  }
# }
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
}

# Clean up temporary variables no longer needed and save results
{rm(list=ls()[grep('^auc',ls())])
rm(list=ls()[grep('^fdr',ls())])
rm(list=ls()[grep('^fpr',ls())])
rm(list=ls()[grep('^mean',ls())])
rm(list=ls()[grep('^names',ls())])
rm(list=ls()[grep('^tpr',ls())])
rm(list=c('i','j','results.list'))}

# Fix results to ensure no extraneous results in each file ####
# Extra results in some files due to imprecise naming meaning that matching to place 
# results in the right files sometimes caught more data than it should have, for fpr, 
# fdr and tpr, e.g. "voom.lfc1" contains "voom", so lfc results for voom were included 
# in DE results, and all DE/DD results for HMs were included in DEDD results because 
# all contain "expHM" or "lnHM". There was also some doubling up, e.g. "noif.DESeq" 
# contains "if.DESeq" as well as "noif.DESeq", but because results are filled by name, 
# these duplications aren't in the results files as they're overwritten by the next 
# match in the loop.
{
  results.DE$fpr <- results.DE$fpr[,-grep('lfc', names(results.DE$fpr))]
  results.DE$fdr <- results.DE$fdr[,-grep('lfc', names(results.DE$fdr))]
  results.DE$tpr <- results.DE$tpr[,-grep('lfc', names(results.DE$tpr))]
  results.DEDD$fpr <- results.DEDD$fpr[,-c(grep('mean', names(results.DEDD$fpr)), grep('disp', names(results.DEDD$fpr)))]
  results.DEDD$fdr <- results.DEDD$fdr[,-c(grep('mean', names(results.DEDD$fdr)), grep('disp', names(results.DEDD$fdr)))]
  results.DEDD$tpr <- results.DEDD$tpr[,-c(grep('mean', names(results.DEDD$tpr)), grep('disp', names(results.DEDD$tpr)))]
}

saveRDS(results.DE, file=here('Results/DD compcodeR data results DESeq norm Sept 2019','DE.results.DD2.DESeqnorm.rds'))
saveRDS(results.DE.lfc1, file=here('Results/DD compcodeR data results DESeq norm Sept 2019','DE.lfc1.results.DD2.DESeqnorm.rds'))
saveRDS(results.DE.lfc2, file=here('Results/DD compcodeR data results DESeq norm Sept 2019','DE.lfc2.results.DD2.DESeqnorm.rds'))
saveRDS(results.DD, file=here('Results/DD compcodeR data results DESeq norm Sept 2019','DD.results.DD2.DESeqnorm.rds'))
saveRDS(results.DD.lfc1, file=here('Results/DD compcodeR data results DESeq norm Sept 2019','DD.lfc1.results.DD2.DESeqnorm.rds'))
saveRDS(results.DD.lfc2, file=here('Results/DD compcodeR data results DESeq norm Sept 2019','DD.lfc2.results.DD2.DESeqnorm.rds'))
saveRDS(results.DEDD, file=here('Results/DD compcodeR data results DESeq norm Sept 2019','DEDD.results.DD2.DESeqnorm.rds'))

