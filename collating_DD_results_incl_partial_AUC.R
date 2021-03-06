library(here)
library(compcodeR)
library(ROCR)
library(caret)
library(qvalue)
source(here('scripts','2019-05-03_bfdr_function.R'))

# First create variables to store results
{
  names.MDSeq <- c('disp.zi.MDSeq', 'disp.nozi.MDSeq', 'disp.zi.lfc1.MDSeq', 'disp.nozi.lfc1.MDSeq',
    'disp.zi.lfc2.MDSeq', 'disp.nozi.lfc2.MDSeq')
  names.expHM <- c('expHM', 'disp.expHM', 'ldisp.expHM', 'disp.lfc1.expHM', 'disp.lfc2.expHM')
  names.lnHM <- c('lnHM', 'disp.lnHM', 'ldisp.lnHM', 'disp.lfc1.lnHM', 'disp.lfc2.lnHM')
  results.list <- c(names.MDSeq, names.expHM, names.lnHM)
}

{
  for (i in names.MDSeq) {
    for (j in c('fpr.raw.','fdr.raw.','tpr.raw.','fpr.fdr.','fdr.fdr.','tpr.fdr.')) {
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
    assign(paste0('pauc.',i), numeric(50))
    assign(paste0('discoveries.',i), matrix(nrow=50, ncol=20000))
    assign(paste0('false.discoveries.',i), matrix(nrow=50, ncol=20000))
  }
}

# Then import each run in turn and save results
for (i in 1:50) {
  # import <- paste0('results.DD2.',i,'.rds')
  import <- paste0('results.DD2.',i,'.DESeqnorm.rds')
  # res <- readRDS(here('Results/DD compcodeR data results July-Aug 2019',import))
  res <- readRDS(here('Results/DD compcodeR data results DESeq norm Sept 2019',import))
  DD <- factor(res$DD, levels=c('1','0'))
  DEDD <- DD
  lfcd1 <- factor(as.numeric(res$lfcd1), levels=c('1','0'))
  lfcd2 <- factor(as.numeric(res$lfcd2), levels=c('1','0'))
  for (j in c(names.expHM[-1], names.lnHM[-1])) {
    assign(paste0('bh.',j), p.adjust(get(paste0('p.',j), res), method='BH'))
    assign(paste0('by.',j), p.adjust(get(paste0('p.',j), res), method='BY'))
    assign(paste0('q.',j), qvalue(get(paste0('p.',j), res))$qval)
  }
  thr.expHM <- sort(res$prob.expHM, decreasing=T)[round(nrow(res$data@count.matrix)*res$prop.expHM)]
  thr.lnHM <- sort(res$prob.lnHM, decreasing=T)[round(nrow(res$data@count.matrix)*res$prop.lnHM)]
  for (j in names.MDSeq) {
    assign(paste0('call.raw.',j), factor(as.numeric(get(paste0('p.',j), res) < 0.05), levels=c('1','0')))
    assign(paste0('call.fdr.',j), factor(as.numeric(get(paste0('q.',j), res) < 0.05), levels=c('1','0')))
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
    assign(paste0('auc.',j), `[<-`(get(paste0('auc.',j)), i, 
                                   value=performance(get(paste0('pred.',j)), measure='auc')@y.values[[1]]))
    assign(paste0('pauc.',j), `[<-`(get(paste0('pauc.',j)), i, 
                                   value=performance(get(paste0('pred.',j)), measure='auc', 
                                                     fpr.stop=0.05)@y.values[[1]]))
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
    assign(paste0('auc.',j), `[<-`(get(paste0('auc.',j)), i, 
                                   value=performance(get(paste0('pred.',j)), measure='auc')@y.values[[1]]))
    assign(paste0('pauc.',j), `[<-`(get(paste0('pauc.',j)), i,
                                    value=performance(get(paste0('pred.',j)), measure='auc', 
                                                      fpr.stop=0.05)@y.values[[1]]))
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
    assign(paste0('pauc.',j), `[<-`(get(paste0('pauc.',j)), i,
                                    value=performance(get(paste0('pred.',j)), measure='auc', 
                                                      fpr.stop=0.05)@y.values[[1]]))
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
rm(list=c('i','j','import','res','DD','DEDD','lfcd1','lfcd2'))
}

# Re-format results to save by experiment type ####

# Set up variables
{
  names.DD <- c('disp.zi.MDSeq', 'disp.nozi.MDSeq', 'disp.expHM', 'ldisp.expHM', 'disp.lnHM', 'ldisp.lnHM')
  names.DD.lfc1 <- c('disp.zi.lfc1.MDSeq', 'disp.nozi.lfc1.MDSeq', 'disp.lfc1.expHM', 'disp.lfc1.lnHM')
  names.DD.lfc2 <- c('disp.zi.lfc2.MDSeq', 'disp.nozi.lfc2.MDSeq', 'disp.lfc2.expHM', 'disp.lfc2.lnHM')
  names.DEDD <- c('expHM', 'lnHM')
}

for (i in c(
  'DD', 'DD.lfc1', 'DD.lfc2',
  'DEDD')) {
  assign(paste0('results.',i), list())
  for (j in c('pauc','auc', 'fpr', 'fdr', 'tpr')) {
    assign(paste0('results.',i), `[[<-`(get(paste0('results.',i)), j, data.frame(matrix(nrow=50, ncol=0))))
  }
  for (j in c('mean.fdr', 'mean.discoveries')) {
    assign(paste0('results.',i), `[[<-`(get(paste0('results.',i)), j, data.frame(matrix(nrow=20000, ncol=0))))
  }
}

{
  names.fpr.fdr.tpr <- vector()
  for (i in names.MDSeq) {
    for (j in c('raw.','fdr.')) {
      names.fpr.fdr.tpr <- c(names.fpr.fdr.tpr,paste0(j,i))
    }
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
}

# Put results into variables
{
  for (i in names.DD) {
    results.DD$auc[[i]] <- get(paste0('auc.',i))
    results.DD$pauc[[i]] <- get(paste0('pauc.',i))
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
    results.DD.lfc1$pauc[[i]] <- get(paste0('pauc.',i))
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
    results.DD.lfc2$pauc[[i]] <- get(paste0('pauc.',i))
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
    results.DEDD$pauc[[i]] <- get(paste0('pauc.',i))
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
{
  rm(list=ls()[grep('^auc',ls())])
  rm(list=ls()[grep('^pauc',ls())])
  rm(list=ls()[grep('^fdr',ls())])
  rm(list=ls()[grep('^fpr',ls())])
  rm(list=ls()[grep('^mean',ls())])
  rm(list=ls()[grep('^names',ls())])
  rm(list=ls()[grep('^tpr',ls())])
  rm(list=c('i','j','results.list'))
}

# Fix and re-save results to ensure no extraneous results in each file ####
# Extra results in some files due to imprecise naming meaning that matching to place 
# results in the right files sometimes caught more data than it should have, for fpr, 
# fdr and tpr, e.g. "voom.lfc1" contains "voom", so lfc results for voom were included 
# in DE results, and all DE/DD results for HMs were included in DEDD results because 
# all contain "expHM" or "lnHM". There was also some doubling up, e.g. "noif.DESeq" 
# contains "if.DESeq" as well as "noif.DESeq", but because results are filled by name, 
# these duplications aren't in the results files as they're overwritten by the next 
# match in the loop.
{
  results.DEDD$fpr <- results.DEDD$fpr[,-c(grep('mean', names(results.DEDD$fpr)), grep('disp', names(results.DEDD$fpr)))]
  results.DEDD$fdr <- results.DEDD$fdr[,-c(grep('mean', names(results.DEDD$fdr)), grep('disp', names(results.DEDD$fdr)))]
  results.DEDD$tpr <- results.DEDD$tpr[,-c(grep('mean', names(results.DEDD$tpr)), grep('disp', names(results.DEDD$tpr)))]
}

saveRDS(results.DD, file=here('Results/DD compcodeR data results DESeq norm Sept 2019','DD.results.DD2.DESeqnorm.rds'))
saveRDS(results.DD.lfc1, file=here('Results/DD compcodeR data results DESeq norm Sept 2019','DD.lfc1.results.DD2.DESeqnorm.rds'))
saveRDS(results.DD.lfc2, file=here('Results/DD compcodeR data results DESeq norm Sept 2019','DD.lfc2.results.DD2.DESeqnorm.rds'))
saveRDS(results.DEDD, file=here('Results/DD compcodeR data results DESeq norm Sept 2019','DEDD.results.DD2.DESeqnorm.rds'))

