library(here)
library(RColorBrewer)

for (i in c('DE', 'DD', 'DEDD')) {
  for (j in c('DEDD2', 'DEDD5', 'DEDD10', 'DEDD20', 'DEDD50')) {
    assign(paste0(i,'.results.TMM.',j), 
           readRDS(here('Results/DEDD compcodeR data results July-Aug 2019', 
                        paste0(i,'.results.',j,'.TMM.rds'))))
  }
}

for (i in c('DE', 'DEDD')) {
  for (j in c('DE2', 'DE5', 'DE10', 'DE20', 'DE50')) {
    assign(paste0(i,'.results.TMM.',j), 
           readRDS(here('Results/DE compcodeR data results July-Aug 2019', 
                        paste0(i,'.results.',j,'.TMM.rds'))))
  }
}

for (i in c('DD', 'DEDD')) {
  for (j in c('DD2', 'DD5', 'DD10', 'DD20', 'DD50')) {
    assign(paste0(i,'.results.TMM.',j), 
           readRDS(here('Results/DD compcodeR data results July-Aug 2019', 
                        paste0(i,'.results.',j,'.TMM.rds'))))
  }
}

for (i in c('DE', 'DD', 'DEDD')) {
  for (j in c('DEDD2', 'DEDD5', 'DEDD10', 'DEDD20', 'DEDD50')) {
    assign(paste0(i,'.results.DESeq.',j), 
           readRDS(here('Results/DEDD compcodeR data results DESeq norm Sept 2019', 
                        paste0(i,'.results.',j,'.DESeqnorm.rds'))))
  }
}

for (i in c('DE', 'DEDD')) {
  for (j in c('DE2', 'DE5', 'DE10', 'DE20', 'DE50')) {
    assign(paste0(i,'.results.DESeq.',j), 
           readRDS(here('Results/DE compcodeR data results DESeq norm Aug-Sept 2019', 
                        paste0(i,'.results.',j,'.DESeqnorm.rds'))))
  }
}

for (i in c('DD', 'DEDD')) {
  for (j in c('DD2', 'DD5', 'DD10', 'DD20', 'DD50')) {
    assign(paste0(i,'.results.DESeq.',j), 
           readRDS(here('Results/DD compcodeR data results DESeq norm Sept 2019', 
                        paste0(i,'.results.',j,'.DESeqnorm.rds'))))
  }
}

for (j in c('DEDD2', 'DEDD5', 'DEDD10', 'DEDD20', 'DEDD50')) {
  assign(paste0('combined_DEDD.results.',j), 
         readRDS(here('Results/Diff dist combining DE, DD predictions Nov 2019', 
                      paste0('diff_dist_results_',j,'.rds'))))
}
rm(i,j)

## DE results for best version of each method (as decided previously - 
## QL for edgeR, untransformed ln for HM, defaults for DESeq2 and MDSeq because doesn't 
## make a difference; hadn't previously decided on DSS but using FDR as has generally ok 
## FDRs where LFDR is very conservative for bigger samples, and FDR has much higher power).
keep.TMM <- c('ql.edgeR', 'voom')
keep.fdr.TMM <- c('fdr.ql.edgeR', 'fdr.voom')
keep.DESeq <- c('if.DESeq', 'notrend.DSS', 'baySeq', 'mean.zi.MDSeq', 'mean.lnHM')
keep.DESeq50 <- c('if.DESeq', 'baySeq', 'mean.zi.MDSeq', 'mean.lnHM')
keep.fdr.DESeq <- c('fdr.if.DESeq', 'fdr.notrend.DSS', 'fdr.baySeq', 'fdr.mean.zi.MDSeq', 'bh.mean.lnHM')
keep.fdr.DESeq50 <- c('fdr.if.DESeq', 'fdr.baySeq', 'fdr.mean.zi.MDSeq', 'bh.mean.lnHM')
DE.results.DE2 <- DE.results.TMM.DE2
DE.results.DE5 <- DE.results.TMM.DE5
DE.results.DE10 <- DE.results.TMM.DE10
DE.results.DE20 <- DE.results.TMM.DE20
DE.results.DE50 <- DE.results.TMM.DE50
DE.results.DEDD2 <- DE.results.TMM.DE2
DE.results.DEDD5 <- DE.results.TMM.DE5
DE.results.DEDD10 <- DE.results.TMM.DE10
DE.results.DEDD20 <- DE.results.TMM.DE20
DE.results.DEDD50 <- DE.results.TMM.DE50
DE.results.DE2$auc <- cbind(DE.results.TMM.DE2$auc[, keep.TMM], DE.results.DESeq.DE2$auc[, keep.DESeq])
DE.results.DE2$pauc <- cbind(DE.results.TMM.DE2$pauc[, keep.TMM], DE.results.DESeq.DE2$pauc[, keep.DESeq])
DE.results.DE2$mean.fdr <- cbind(DE.results.TMM.DE2$mean.fdr[, keep.TMM], DE.results.DESeq.DE2$mean.fdr[, keep.DESeq])
DE.results.DE2$mean.discoveries <- cbind(DE.results.TMM.DE2$mean.discoveries[, keep.TMM], DE.results.DESeq.DE2$mean.discoveries[, keep.DESeq])
DE.results.DE2$fpr <- cbind(DE.results.TMM.DE2$fpr[, keep.fdr.TMM], DE.results.DESeq.DE2$fpr[, keep.fdr.DESeq])
DE.results.DE2$fdr <- cbind(DE.results.TMM.DE2$fdr[, keep.fdr.TMM], DE.results.DESeq.DE2$fdr[, keep.fdr.DESeq])
DE.results.DE2$tpr <- cbind(DE.results.TMM.DE2$tpr[, keep.fdr.TMM], DE.results.DESeq.DE2$tpr[, keep.fdr.DESeq])
DE.results.DE5$auc <- cbind(DE.results.TMM.DE5$auc[, keep.TMM], DE.results.DESeq.DE5$auc[, keep.DESeq])
DE.results.DE5$pauc <- cbind(DE.results.TMM.DE5$pauc[, keep.TMM], DE.results.DESeq.DE5$pauc[, keep.DESeq])
DE.results.DE5$mean.fdr <- cbind(DE.results.TMM.DE5$mean.fdr[, keep.TMM], DE.results.DESeq.DE5$mean.fdr[, keep.DESeq])
DE.results.DE5$mean.discoveries <- cbind(DE.results.TMM.DE5$mean.discoveries[, keep.TMM], DE.results.DESeq.DE5$mean.discoveries[, keep.DESeq])
DE.results.DE5$fpr <- cbind(DE.results.TMM.DE5$fpr[, keep.fdr.TMM], DE.results.DESeq.DE5$fpr[, keep.fdr.DESeq])
DE.results.DE5$fdr <- cbind(DE.results.TMM.DE5$fdr[, keep.fdr.TMM], DE.results.DESeq.DE5$fdr[, keep.fdr.DESeq])
DE.results.DE5$tpr <- cbind(DE.results.TMM.DE5$tpr[, keep.fdr.TMM], DE.results.DESeq.DE5$tpr[, keep.fdr.DESeq])
DE.results.DE10$auc <- cbind(DE.results.TMM.DE10$auc[, keep.TMM], DE.results.DESeq.DE10$auc[, keep.DESeq])
DE.results.DE10$pauc <- cbind(DE.results.TMM.DE10$pauc[, keep.TMM], DE.results.DESeq.DE10$pauc[, keep.DESeq])
DE.results.DE10$mean.fdr <- cbind(DE.results.TMM.DE10$mean.fdr[, keep.TMM], DE.results.DESeq.DE10$mean.fdr[, keep.DESeq])
DE.results.DE10$mean.discoveries <- cbind(DE.results.TMM.DE10$mean.discoveries[, keep.TMM], DE.results.DESeq.DE10$mean.discoveries[, keep.DESeq])
DE.results.DE10$fpr <- cbind(DE.results.TMM.DE10$fpr[, keep.fdr.TMM], DE.results.DESeq.DE10$fpr[, keep.fdr.DESeq])
DE.results.DE10$fdr <- cbind(DE.results.TMM.DE10$fdr[, keep.fdr.TMM], DE.results.DESeq.DE10$fdr[, keep.fdr.DESeq])
DE.results.DE10$tpr <- cbind(DE.results.TMM.DE10$tpr[, keep.fdr.TMM], DE.results.DESeq.DE10$tpr[, keep.fdr.DESeq])
DE.results.DE20$auc <- cbind(DE.results.TMM.DE20$auc[, keep.TMM], DE.results.DESeq.DE20$auc[, keep.DESeq])
DE.results.DE20$pauc <- cbind(DE.results.TMM.DE20$pauc[, keep.TMM], DE.results.DESeq.DE20$pauc[, keep.DESeq])
DE.results.DE20$mean.fdr <- cbind(DE.results.TMM.DE20$mean.fdr[, keep.TMM], DE.results.DESeq.DE20$mean.fdr[, keep.DESeq])
DE.results.DE20$mean.discoveries <- cbind(DE.results.TMM.DE20$mean.discoveries[, keep.TMM], DE.results.DESeq.DE20$mean.discoveries[, keep.DESeq])
DE.results.DE20$fpr <- cbind(DE.results.TMM.DE20$fpr[, keep.fdr.TMM], DE.results.DESeq.DE20$fpr[, keep.fdr.DESeq])
DE.results.DE20$fdr <- cbind(DE.results.TMM.DE20$fdr[, keep.fdr.TMM], DE.results.DESeq.DE20$fdr[, keep.fdr.DESeq])
DE.results.DE20$tpr <- cbind(DE.results.TMM.DE20$tpr[, keep.fdr.TMM], DE.results.DESeq.DE20$tpr[, keep.fdr.DESeq])
DE.results.DE50$auc <- cbind(DE.results.TMM.DE50$auc[, keep.TMM], DE.results.DESeq.DE50$auc[, keep.DESeq50])
DE.results.DE50$pauc <- cbind(DE.results.TMM.DE50$pauc[, keep.TMM], DE.results.DESeq.DE50$pauc[, keep.DESeq50])
DE.results.DE50$mean.fdr <- cbind(DE.results.TMM.DE50$mean.fdr[, keep.TMM], DE.results.DESeq.DE50$mean.fdr[, keep.DESeq50])
DE.results.DE50$mean.discoveries <- cbind(DE.results.TMM.DE50$mean.discoveries[, keep.TMM], DE.results.DESeq.DE50$mean.discoveries[, keep.DESeq50])
DE.results.DE50$fpr <- cbind(DE.results.TMM.DE50$fpr[, keep.fdr.TMM], DE.results.DESeq.DE50$fpr[, keep.fdr.DESeq50])
DE.results.DE50$fdr <- cbind(DE.results.TMM.DE50$fdr[, keep.fdr.TMM], DE.results.DESeq.DE50$fdr[, keep.fdr.DESeq50])
DE.results.DE50$tpr <- cbind(DE.results.TMM.DE50$tpr[, keep.fdr.TMM], DE.results.DESeq.DE50$tpr[, keep.fdr.DESeq50])
DE.results.DEDD2$auc <- cbind(DE.results.TMM.DEDD2$auc[, keep.TMM], DE.results.DESeq.DEDD2$auc[, keep.DESeq])
DE.results.DEDD2$pauc <- cbind(DE.results.TMM.DEDD2$pauc[, keep.TMM], DE.results.DESeq.DEDD2$pauc[, keep.DESeq])
DE.results.DEDD2$mean.fdr <- cbind(DE.results.TMM.DEDD2$mean.fdr[, keep.TMM], DE.results.DESeq.DEDD2$mean.fdr[, keep.DESeq])
DE.results.DEDD2$mean.discoveries <- cbind(DE.results.TMM.DEDD2$mean.discoveries[, keep.TMM], DE.results.DESeq.DEDD2$mean.discoveries[, keep.DESeq])
DE.results.DEDD2$fpr <- cbind(DE.results.TMM.DEDD2$fpr[, keep.fdr.TMM], DE.results.DESeq.DEDD2$fpr[, keep.fdr.DESeq])
DE.results.DEDD2$fdr <- cbind(DE.results.TMM.DEDD2$fdr[, keep.fdr.TMM], DE.results.DESeq.DEDD2$fdr[, keep.fdr.DESeq])
DE.results.DEDD2$tpr <- cbind(DE.results.TMM.DEDD2$tpr[, keep.fdr.TMM], DE.results.DESeq.DEDD2$tpr[, keep.fdr.DESeq])
DE.results.DEDD5$auc <- cbind(DE.results.TMM.DEDD5$auc[, keep.TMM], DE.results.DESeq.DEDD5$auc[, keep.DESeq])
DE.results.DEDD5$pauc <- cbind(DE.results.TMM.DEDD5$pauc[, keep.TMM], DE.results.DESeq.DEDD5$pauc[, keep.DESeq])
DE.results.DEDD5$mean.fdr <- cbind(DE.results.TMM.DEDD5$mean.fdr[, keep.TMM], DE.results.DESeq.DEDD5$mean.fdr[, keep.DESeq])
DE.results.DEDD5$mean.discoveries <- cbind(DE.results.TMM.DEDD5$mean.discoveries[, keep.TMM], DE.results.DESeq.DEDD5$mean.discoveries[, keep.DESeq])
DE.results.DEDD5$fpr <- cbind(DE.results.TMM.DEDD5$fpr[, keep.fdr.TMM], DE.results.DESeq.DEDD5$fpr[, keep.fdr.DESeq])
DE.results.DEDD5$fdr <- cbind(DE.results.TMM.DEDD5$fdr[, keep.fdr.TMM], DE.results.DESeq.DEDD5$fdr[, keep.fdr.DESeq])
DE.results.DEDD5$tpr <- cbind(DE.results.TMM.DEDD5$tpr[, keep.fdr.TMM], DE.results.DESeq.DEDD5$tpr[, keep.fdr.DESeq])
DE.results.DEDD10$auc <- cbind(DE.results.TMM.DEDD10$auc[, keep.TMM], DE.results.DESeq.DEDD10$auc[, keep.DESeq])
DE.results.DEDD10$pauc <- cbind(DE.results.TMM.DEDD10$pauc[, keep.TMM], DE.results.DESeq.DEDD10$pauc[, keep.DESeq])
DE.results.DEDD10$mean.fdr <- cbind(DE.results.TMM.DEDD10$mean.fdr[, keep.TMM], DE.results.DESeq.DEDD10$mean.fdr[, keep.DESeq])
DE.results.DEDD10$mean.discoveries <- cbind(DE.results.TMM.DEDD10$mean.discoveries[, keep.TMM], DE.results.DESeq.DEDD10$mean.discoveries[, keep.DESeq])
DE.results.DEDD10$fpr <- cbind(DE.results.TMM.DEDD10$fpr[, keep.fdr.TMM], DE.results.DESeq.DEDD10$fpr[, keep.fdr.DESeq])
DE.results.DEDD10$fdr <- cbind(DE.results.TMM.DEDD10$fdr[, keep.fdr.TMM], DE.results.DESeq.DEDD10$fdr[, keep.fdr.DESeq])
DE.results.DEDD10$tpr <- cbind(DE.results.TMM.DEDD10$tpr[, keep.fdr.TMM], DE.results.DESeq.DEDD10$tpr[, keep.fdr.DESeq])
DE.results.DEDD20$auc <- cbind(DE.results.TMM.DEDD20$auc[, keep.TMM], DE.results.DESeq.DEDD20$auc[, keep.DESeq])
DE.results.DEDD20$pauc <- cbind(DE.results.TMM.DEDD20$pauc[, keep.TMM], DE.results.DESeq.DEDD20$pauc[, keep.DESeq])
DE.results.DEDD20$mean.fdr <- cbind(DE.results.TMM.DEDD20$mean.fdr[, keep.TMM], DE.results.DESeq.DEDD20$mean.fdr[, keep.DESeq])
DE.results.DEDD20$mean.discoveries <- cbind(DE.results.TMM.DEDD20$mean.discoveries[, keep.TMM], DE.results.DESeq.DEDD20$mean.discoveries[, keep.DESeq])
DE.results.DEDD20$fpr <- cbind(DE.results.TMM.DEDD20$fpr[, keep.fdr.TMM], DE.results.DESeq.DEDD20$fpr[, keep.fdr.DESeq])
DE.results.DEDD20$fdr <- cbind(DE.results.TMM.DEDD20$fdr[, keep.fdr.TMM], DE.results.DESeq.DEDD20$fdr[, keep.fdr.DESeq])
DE.results.DEDD20$tpr <- cbind(DE.results.TMM.DEDD20$tpr[, keep.fdr.TMM], DE.results.DESeq.DEDD20$tpr[, keep.fdr.DESeq])
DE.results.DEDD50$auc <- cbind(DE.results.TMM.DEDD50$auc[, keep.TMM], DE.results.DESeq.DEDD50$auc[, keep.DESeq50])
DE.results.DEDD50$pauc <- cbind(DE.results.TMM.DEDD50$pauc[, keep.TMM], DE.results.DESeq.DEDD50$pauc[, keep.DESeq50])
DE.results.DEDD50$mean.fdr <- cbind(DE.results.TMM.DEDD50$mean.fdr[, keep.TMM], DE.results.DESeq.DEDD50$mean.fdr[, keep.DESeq50])
DE.results.DEDD50$mean.discoveries <- cbind(DE.results.TMM.DEDD50$mean.discoveries[, keep.TMM], DE.results.DESeq.DEDD50$mean.discoveries[, keep.DESeq50])
DE.results.DEDD50$fpr <- cbind(DE.results.TMM.DEDD50$fpr[, keep.fdr.TMM], DE.results.DESeq.DEDD50$fpr[, keep.fdr.DESeq50])
DE.results.DEDD50$fdr <- cbind(DE.results.TMM.DEDD50$fdr[, keep.fdr.TMM], DE.results.DESeq.DEDD50$fdr[, keep.fdr.DESeq50])
DE.results.DEDD50$tpr <- cbind(DE.results.TMM.DEDD50$tpr[, keep.fdr.TMM], DE.results.DESeq.DEDD50$tpr[, keep.fdr.DESeq50])

## DD results
keep.DESeq <- c('disp.zi.MDSeq', 'disp.lnHM')
keep.fdr.DESeq <- c('fdr.disp.zi.MDSeq', 'bh.disp.lnHM')
DD.results.DD2 <- DD.results.DESeq.DD2
DD.results.DD5 <- DD.results.DESeq.DD5
DD.results.DD10 <- DD.results.DESeq.DD10
DD.results.DD20 <- DD.results.DESeq.DD20
DD.results.DD50 <- DD.results.DESeq.DD50
DD.results.DEDD2 <- DD.results.DESeq.DD2
DD.results.DEDD5 <- DD.results.DESeq.DD5
DD.results.DEDD10 <- DD.results.DESeq.DD10
DD.results.DEDD20 <- DD.results.DESeq.DD20
DD.results.DEDD50 <- DD.results.DESeq.DD50
DD.results.DD2$auc <- DD.results.DESeq.DD2$auc[, keep.DESeq]
DD.results.DD2$pauc <- DD.results.DESeq.DD2$pauc[, keep.DESeq]
DD.results.DD2$mean.fdr <- DD.results.DESeq.DD2$mean.fdr[, keep.DESeq]
DD.results.DD2$mean.discoveries <- DD.results.DESeq.DD2$mean.discoveries[, keep.DESeq]
DD.results.DD2$fpr <- DD.results.DESeq.DD2$fpr[, keep.fdr.DESeq]
DD.results.DD2$fdr <- DD.results.DESeq.DD2$fdr[, keep.fdr.DESeq]
DD.results.DD2$tpr <- DD.results.DESeq.DD2$tpr[, keep.fdr.DESeq]
DD.results.DD5$auc <- DD.results.DESeq.DD5$auc[, keep.DESeq]
DD.results.DD5$pauc <- DD.results.DESeq.DD5$pauc[, keep.DESeq]
DD.results.DD5$mean.fdr <- DD.results.DESeq.DD5$mean.fdr[, keep.DESeq]
DD.results.DD5$mean.discoveries <- DD.results.DESeq.DD5$mean.discoveries[, keep.DESeq]
DD.results.DD5$fpr <- DD.results.DESeq.DD5$fpr[, keep.fdr.DESeq]
DD.results.DD5$fdr <- DD.results.DESeq.DD5$fdr[, keep.fdr.DESeq]
DD.results.DD5$tpr <- DD.results.DESeq.DD5$tpr[, keep.fdr.DESeq]
DD.results.DD10$auc <- DD.results.DESeq.DD10$auc[, keep.DESeq]
DD.results.DD10$pauc <- DD.results.DESeq.DD10$pauc[, keep.DESeq]
DD.results.DD10$mean.fdr <- DD.results.DESeq.DD10$mean.fdr[, keep.DESeq]
DD.results.DD10$mean.discoveries <- DD.results.DESeq.DD10$mean.discoveries[, keep.DESeq]
DD.results.DD10$fpr <- DD.results.DESeq.DD10$fpr[, keep.fdr.DESeq]
DD.results.DD10$fdr <- DD.results.DESeq.DD10$fdr[, keep.fdr.DESeq]
DD.results.DD10$tpr <- DD.results.DESeq.DD10$tpr[, keep.fdr.DESeq]
DD.results.DD20$auc <- DD.results.DESeq.DD20$auc[, keep.DESeq]
DD.results.DD20$pauc <- DD.results.DESeq.DD20$pauc[, keep.DESeq]
DD.results.DD20$mean.fdr <- DD.results.DESeq.DD20$mean.fdr[, keep.DESeq]
DD.results.DD20$mean.discoveries <- DD.results.DESeq.DD20$mean.discoveries[, keep.DESeq]
DD.results.DD20$fpr <- DD.results.DESeq.DD20$fpr[, keep.fdr.DESeq]
DD.results.DD20$fdr <- DD.results.DESeq.DD20$fdr[, keep.fdr.DESeq]
DD.results.DD20$tpr <- DD.results.DESeq.DD20$tpr[, keep.fdr.DESeq]
DD.results.DD50$auc <- DD.results.DESeq.DD50$auc[, keep.DESeq]
DD.results.DD50$pauc <- DD.results.DESeq.DD50$pauc[, keep.DESeq]
DD.results.DD50$mean.fdr <- DD.results.DESeq.DD50$mean.fdr[, keep.DESeq]
DD.results.DD50$mean.discoveries <- DD.results.DESeq.DD50$mean.discoveries[, keep.DESeq]
DD.results.DD50$fpr <- DD.results.DESeq.DD50$fpr[, keep.fdr.DESeq]
DD.results.DD50$fdr <- DD.results.DESeq.DD50$fdr[, keep.fdr.DESeq]
DD.results.DD50$tpr <- DD.results.DESeq.DD50$tpr[, keep.fdr.DESeq]
DD.results.DEDD2$auc <- DD.results.DESeq.DEDD2$auc[, keep.DESeq]
DD.results.DEDD2$pauc <- DD.results.DESeq.DEDD2$pauc[, keep.DESeq]
DD.results.DEDD2$mean.fdr <- DD.results.DESeq.DEDD2$mean.fdr[, keep.DESeq]
DD.results.DEDD2$mean.discoveries <- DD.results.DESeq.DEDD2$mean.discoveries[, keep.DESeq]
DD.results.DEDD2$fpr <- DD.results.DESeq.DEDD2$fpr[, keep.fdr.DESeq]
DD.results.DEDD2$fdr <- DD.results.DESeq.DEDD2$fdr[, keep.fdr.DESeq]
DD.results.DEDD2$tpr <- DD.results.DESeq.DEDD2$tpr[, keep.fdr.DESeq]
DD.results.DEDD5$auc <- DD.results.DESeq.DEDD5$auc[, keep.DESeq]
DD.results.DEDD5$pauc <- DD.results.DESeq.DEDD5$pauc[, keep.DESeq]
DD.results.DEDD5$mean.fdr <- DD.results.DESeq.DEDD5$mean.fdr[, keep.DESeq]
DD.results.DEDD5$mean.discoveries <- DD.results.DESeq.DEDD5$mean.discoveries[, keep.DESeq]
DD.results.DEDD5$fpr <- DD.results.DESeq.DEDD5$fpr[, keep.fdr.DESeq]
DD.results.DEDD5$fdr <- DD.results.DESeq.DEDD5$fdr[, keep.fdr.DESeq]
DD.results.DEDD5$tpr <- DD.results.DESeq.DEDD5$tpr[, keep.fdr.DESeq]
DD.results.DEDD10$auc <- DD.results.DESeq.DEDD10$auc[, keep.DESeq]
DD.results.DEDD10$pauc <- DD.results.DESeq.DEDD10$pauc[, keep.DESeq]
DD.results.DEDD10$mean.fdr <- DD.results.DESeq.DEDD10$mean.fdr[, keep.DESeq]
DD.results.DEDD10$mean.discoveries <- DD.results.DESeq.DEDD10$mean.discoveries[, keep.DESeq]
DD.results.DEDD10$fpr <- DD.results.DESeq.DEDD10$fpr[, keep.fdr.DESeq]
DD.results.DEDD10$fdr <- DD.results.DESeq.DEDD10$fdr[, keep.fdr.DESeq]
DD.results.DEDD10$tpr <- DD.results.DESeq.DEDD10$tpr[, keep.fdr.DESeq]
DD.results.DEDD20$auc <- DD.results.DESeq.DEDD20$auc[, keep.DESeq]
DD.results.DEDD20$pauc <- DD.results.DESeq.DEDD20$pauc[, keep.DESeq]
DD.results.DEDD20$mean.fdr <- DD.results.DESeq.DEDD20$mean.fdr[, keep.DESeq]
DD.results.DEDD20$mean.discoveries <- DD.results.DESeq.DEDD20$mean.discoveries[, keep.DESeq]
DD.results.DEDD20$fpr <- DD.results.DESeq.DEDD20$fpr[, keep.fdr.DESeq]
DD.results.DEDD20$fdr <- DD.results.DESeq.DEDD20$fdr[, keep.fdr.DESeq]
DD.results.DEDD20$tpr <- DD.results.DESeq.DEDD20$tpr[, keep.fdr.DESeq]
DD.results.DEDD50$auc <- DD.results.DESeq.DEDD50$auc[, keep.DESeq]
DD.results.DEDD50$pauc <- DD.results.DESeq.DEDD50$pauc[, keep.DESeq]
DD.results.DEDD50$mean.fdr <- DD.results.DESeq.DEDD50$mean.fdr[, keep.DESeq]
DD.results.DEDD50$mean.discoveries <- DD.results.DESeq.DEDD50$mean.discoveries[, keep.DESeq]
DD.results.DEDD50$fpr <- DD.results.DESeq.DEDD50$fpr[, keep.fdr.DESeq]
DD.results.DEDD50$fdr <- DD.results.DESeq.DEDD50$fdr[, keep.fdr.DESeq]
DD.results.DEDD50$tpr <- DD.results.DESeq.DEDD50$tpr[, keep.fdr.DESeq]


n <- 23
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual' & 
                                  brewer.pal.info$colorblind == T,]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# pie(rep(1,n), col=col_vector)
col_vector <- col_vector[-c(7,10,11,12,20)]

#################################
#### Differential dispersion ####
#################################
DD.results.2 <- list()
DD.results.5 <- list()
DD.results.10 <- list()
DD.results.20 <- list()
DD.results.50 <- list()

###########
### AUC ###
DD.results.2$auc <- cbind(DD.results.DD2$auc, DD.results.DEDD2$auc)
names(DD.results.2$auc) <- c('MDSeq.DD', 'HM.DD', 'MDSeq.DEDD', 'HM.DEDD')
DD.results.2$auc <- DD.results.2$auc[, c(1,3,2,4)]
DD.results.5$auc <- cbind(DD.results.DD5$auc, DD.results.DEDD5$auc)
names(DD.results.5$auc) <- c('MDSeq.DD', 'HM.DD', 'MDSeq.DEDD', 'HM.DEDD')
DD.results.5$auc <- DD.results.5$auc[, c(1,3,2,4)]
DD.results.10$auc <- cbind(DD.results.DD10$auc, DD.results.DEDD10$auc)
names(DD.results.10$auc) <- c('MDSeq.DD', 'HM.DD', 'MDSeq.DEDD', 'HM.DEDD')
DD.results.10$auc <- DD.results.10$auc[, c(1,3,2,4)]
DD.results.20$auc <- cbind(DD.results.DD20$auc, DD.results.DEDD20$auc)
names(DD.results.20$auc) <- c('MDSeq.DD', 'HM.DD', 'MDSeq.DEDD', 'HM.DEDD')
DD.results.20$auc <- DD.results.20$auc[, c(1,3,2,4)]
DD.results.50$auc <- cbind(DD.results.DD50$auc, DD.results.DEDD50$auc)
names(DD.results.50$auc) <- c('MDSeq.DD', 'HM.DD', 'MDSeq.DEDD', 'HM.DEDD')
DD.results.50$auc <- DD.results.50$auc[, c(1,3,2,4)]

par(mfrow=c(1,2), mar=c(1,2.5,4,1), mgp=c(3,0.7,0))
boxplot(cbind(DD.results.2$auc[c(1,3)], 
        DD.results.5$auc[c(1,3)], 
        DD.results.10$auc[c(1,3)], 
        DD.results.20$auc[c(1,3)], 
        DD.results.50$auc[c(1,3)]), 
        at=1:10 + c(0,0,0.5,0.5,1.0,1.0,1.5,1.5,2.0,2.0), 
        names=NA, col=col_vector[1:2], xaxt='n', ylim=c(0.5,0.9), 
        pch=20, cex.axis=2, 
        main=paste0("Differences in dispersion only\n", 
                            "2, 5, 10, 20, 50 samples per group"), 
        cex.main=2)
legend("topleft", fill=col_vector[1:2], bty='n', 
       legend=c("MDSeq", "HM"), cex=2)
boxplot(cbind(DD.results.2$auc[c(2,4)], 
              DD.results.5$auc[c(2,4)], 
              DD.results.10$auc[c(2,4)], 
              DD.results.20$auc[c(2,4)], 
              DD.results.50$auc[c(2,4)]), 
        at=1:10 + c(0,0,0.5,0.5,1.0,1.0,1.5,1.5,2.0,2.0), 
        names=NA, col=col_vector[1:2], xaxt='n', ylim=c(0.5,0.9), 
        pch=20, cex.axis=2, 
        main=paste0("Differences in mean and dispersion\n", 
                            "2, 5, 10, 20, 50 samples per group"), 
        cex.main=2)
legend("topleft", fill=col_vector[1:2], bty='n', 
       legend=c("MDSeq", "HM"), cex=2)


###########
### pAUC ###
DD.results.2$pauc <- cbind(DD.results.DD2$pauc, DD.results.DEDD2$pauc)
names(DD.results.2$pauc) <- c('MDSeq.DD', 'HM.DD', 'MDSeq.DEDD', 'HM.DEDD')
DD.results.2$pauc <- DD.results.2$pauc[, c(1,3,2,4)]
DD.results.5$pauc <- cbind(DD.results.DD5$pauc, DD.results.DEDD5$pauc)
names(DD.results.5$pauc) <- c('MDSeq.DD', 'HM.DD', 'MDSeq.DEDD', 'HM.DEDD')
DD.results.5$pauc <- DD.results.5$pauc[, c(1,3,2,4)]
DD.results.10$pauc <- cbind(DD.results.DD10$pauc, DD.results.DEDD10$pauc)
names(DD.results.10$pauc) <- c('MDSeq.DD', 'HM.DD', 'MDSeq.DEDD', 'HM.DEDD')
DD.results.10$pauc <- DD.results.10$pauc[, c(1,3,2,4)]
DD.results.20$pauc <- cbind(DD.results.DD20$pauc, DD.results.DEDD20$pauc)
names(DD.results.20$pauc) <- c('MDSeq.DD', 'HM.DD', 'MDSeq.DEDD', 'HM.DEDD')
DD.results.20$pauc <- DD.results.20$pauc[, c(1,3,2,4)]
DD.results.50$pauc <- cbind(DD.results.DD50$pauc, DD.results.DEDD50$pauc)
names(DD.results.50$pauc) <- c('MDSeq.DD', 'HM.DD', 'MDSeq.DEDD', 'HM.DEDD')
DD.results.50$pauc <- DD.results.50$pauc[, c(1,3,2,4)]

par(mfrow=c(1,2), mar=c(1,2,3,1), mgp=c(3,0.7,0))
boxplot(cbind(DD.results.2$pauc[c(1,3)], 
              DD.results.5$pauc[c(1,3)], 
              DD.results.10$pauc[c(1,3)], 
              DD.results.20$pauc[c(1,3)], 
              DD.results.50$pauc[c(1,3)]), 
        at=1:10 + c(0,0,0.5,0.5,1.0,1.0,1.5,1.5,2.0,2.0), 
        names=NA, col=col_vector[1:2], xaxt='n', ylim=c(0,0.03), 
        pch=20, main=paste0("Parial AUC, differences in dispersion only\n", 
                            "2, 5, 10, 20, 50 samples per group"))
legend("topleft", fill=col_vector[1:2], bty='n', 
       legend=c("MDSeq", "HM"), cex=1.2)
boxplot(cbind(DD.results.2$pauc[c(2,4)], 
              DD.results.5$pauc[c(2,4)], 
              DD.results.10$pauc[c(2,4)], 
              DD.results.20$pauc[c(2,4)], 
              DD.results.50$pauc[c(2,4)]), 
        at=1:10 + c(0,0,0.5,0.5,1.0,1.0,1.5,1.5,2.0,2.0), 
        names=NA, col=col_vector[1:2], xaxt='n', ylim=c(0,0.03), 
        pch=20, main=paste0("Parial AUC, differences in mean and dispersion\n", 
                            "2, 5, 10, 20, 50 samples per group"))
legend("topleft", fill=col_vector[1:2], bty='n', 
       legend=c("MDSeq", "HM"), cex=1.2)


###########
### FDR ###
DD.results.2$fdr <- cbind(DD.results.DD2$fdr, DD.results.DEDD2$fdr)
names(DD.results.2$fdr) <- c('MDSeq.DD', 'HM.DD', 'MDSeq.DEDD', 'HM.DEDD')
DD.results.2$fdr <- DD.results.2$fdr[, c(1,3,2,4)]
DD.results.5$fdr <- cbind(DD.results.DD5$fdr, DD.results.DEDD5$fdr)
names(DD.results.5$fdr) <- c('MDSeq.DD', 'HM.DD', 'MDSeq.DEDD', 'HM.DEDD')
DD.results.5$fdr <- DD.results.5$fdr[, c(1,3,2,4)]
DD.results.10$fdr <- cbind(DD.results.DD10$fdr, DD.results.DEDD10$fdr)
names(DD.results.10$fdr) <- c('MDSeq.DD', 'HM.DD', 'MDSeq.DEDD', 'HM.DEDD')
DD.results.10$fdr <- DD.results.10$fdr[, c(1,3,2,4)]
DD.results.20$fdr <- cbind(DD.results.DD20$fdr, DD.results.DEDD20$fdr)
names(DD.results.20$fdr) <- c('MDSeq.DD', 'HM.DD', 'MDSeq.DEDD', 'HM.DEDD')
DD.results.20$fdr <- DD.results.20$fdr[, c(1,3,2,4)]
DD.results.50$fdr <- cbind(DD.results.DD50$fdr, DD.results.DEDD50$fdr)
names(DD.results.50$fdr) <- c('MDSeq.DD', 'HM.DD', 'MDSeq.DEDD', 'HM.DEDD')
DD.results.50$fdr <- DD.results.50$fdr[, c(1,3,2,4)]

par(mfrow=c(1,2), mar=c(1,2.5,4,1), mgp=c(3,0.7,0))
boxplot(cbind(DD.results.2$fdr[c(1,3)], 
              DD.results.5$fdr[c(1,3)], 
              DD.results.10$fdr[c(1,3)], 
              DD.results.20$fdr[c(1,3)], 
              DD.results.50$fdr[c(1,3)]), 
        at=1:10 + c(0,0,0.5,0.5,1.0,1.0,1.5,1.5,2.0,2.0),
        names=NA, col=col_vector[1:2], xaxt='n', ylim=c(0,1), 
        pch=20, cex.axis=2, 
        main=paste0("Differences in dispersion only\n", 
                    "2, 5, 10, 20, 50 samples per group"), 
        cex.main=2)
abline(h=0.05, col='grey')
legend("topright", fill=col_vector[1:2], bty='n', cex=2, 
       legend=c("MDSeq", "HM"))
boxplot(cbind(DD.results.2$fdr[c(2,4)], 
              DD.results.5$fdr[c(2,4)], 
              DD.results.10$fdr[c(2,4)], 
              DD.results.20$fdr[c(2,4)], 
              DD.results.50$fdr[c(2,4)]), 
        at=1:10 + c(0,0,0.5,0.5,1.0,1.0,1.5,1.5,2.0,2.0),
        names=NA, col=col_vector[1:2], xaxt='n', ylim=c(0,1), 
        pch=20, cex.axis=2, 
        main=paste0("Differences in mean and dispersion\n", 
                    "2, 5, 10, 20, 50 samples per group"), 
        cex.main=2)
abline(h=0.05, col='grey')
legend("topright", fill=col_vector[1:2], bty='n', cex=2, 
       legend=c("MDSeq", "HM"))


###########
### TPR ###
DD.results.2$tpr <- cbind(DD.results.DD2$tpr, DD.results.DEDD2$tpr)
names(DD.results.2$tpr) <- c('MDSeq.DD', 'HM.DD', 'MDSeq.DEDD', 'HM.DEDD')
DD.results.2$tpr <- DD.results.2$tpr[, c(1,3,2,4)]
DD.results.5$tpr <- cbind(DD.results.DD5$tpr, DD.results.DEDD5$tpr)
names(DD.results.5$tpr) <- c('MDSeq.DD', 'HM.DD', 'MDSeq.DEDD', 'HM.DEDD')
DD.results.5$tpr <- DD.results.5$tpr[, c(1,3,2,4)]
DD.results.10$tpr <- cbind(DD.results.DD10$tpr, DD.results.DEDD10$tpr)
names(DD.results.10$tpr) <- c('MDSeq.DD', 'HM.DD', 'MDSeq.DEDD', 'HM.DEDD')
DD.results.10$tpr <- DD.results.10$tpr[, c(1,3,2,4)]
DD.results.20$tpr <- cbind(DD.results.DD20$tpr, DD.results.DEDD20$tpr)
names(DD.results.20$tpr) <- c('MDSeq.DD', 'HM.DD', 'MDSeq.DEDD', 'HM.DEDD')
DD.results.20$tpr <- DD.results.20$tpr[, c(1,3,2,4)]
DD.results.50$tpr <- cbind(DD.results.DD50$tpr, DD.results.DEDD50$tpr)
names(DD.results.50$tpr) <- c('MDSeq.DD', 'HM.DD', 'MDSeq.DEDD', 'HM.DEDD')
DD.results.50$tpr <- DD.results.50$tpr[, c(1,3,2,4)]

par(mfrow=c(1,2), mar=c(1,2,3,1), mgp=c(3,0.7,0))
boxplot(cbind(DD.results.2$tpr[c(1,3)], 
              DD.results.5$tpr[c(1,3)], 
              DD.results.10$tpr[c(1,3)], 
              DD.results.20$tpr[c(1,3)], 
              DD.results.50$tpr[c(1,3)]), 
        at=1:10 + c(0,0,0.5,0.5,1.0,1.0,1.5,1.5,2.0,2.0),
        names=NA, col=col_vector[1:2], xaxt='n', ylim=c(0,0.4), 
        pch=20, main=paste0("TPR, differences in dispersion only\n", 
                            "2, 5, 10, 20, 50 samples per group"))
legend("topleft", fill=col_vector[1:2], bty='n', cex=1.2, 
       legend=c("MDSeq", "HM"))
boxplot(cbind(DD.results.2$tpr[c(2,4)], 
              DD.results.5$tpr[c(2,4)], 
              DD.results.10$tpr[c(2,4)], 
              DD.results.20$tpr[c(2,4)], 
              DD.results.50$tpr[c(2,4)]), 
        at=1:10 + c(0,0,0.5,0.5,1.0,1.0,1.5,1.5,2.0,2.0),
        names=NA, col=col_vector[1:2], xaxt='n', ylim=c(0,0.4), 
        pch=20, main=paste0("TPR, differences in mean and dispersion\n", 
                            "2, 5, 10, 20, 50 samples per group"))
legend("topleft", fill=col_vector[1:2], bty='n', cex=1.2, 
       legend=c("MDSeq", "HM"))


#############################
### False discovery plots ###
DD.results.2$mean.discoveries <- cbind(DD.results.DD2$mean.discoveries, DD.results.DEDD2$mean.discoveries)
names(DD.results.2$mean.discoveries) <- c('MDSeq.DD', 'HM.DD', 'MDSeq.DEDD', 'HM.DEDD')
DD.results.2$mean.discoveries <- DD.results.2$mean.discoveries[, c(1,3,2,4)]
DD.results.5$mean.discoveries <- cbind(DD.results.DD5$mean.discoveries, DD.results.DEDD5$mean.discoveries)
names(DD.results.5$mean.discoveries) <- c('MDSeq.DD', 'HM.DD', 'MDSeq.DEDD', 'HM.DEDD')
DD.results.5$mean.discoveries <- DD.results.5$mean.discoveries[, c(1,3,2,4)]
DD.results.10$mean.discoveries <- cbind(DD.results.DD10$mean.discoveries, DD.results.DEDD10$mean.discoveries)
names(DD.results.10$mean.discoveries) <- c('MDSeq.DD', 'HM.DD', 'MDSeq.DEDD', 'HM.DEDD')
DD.results.10$mean.discoveries <- DD.results.10$mean.discoveries[, c(1,3,2,4)]
DD.results.20$mean.discoveries <- cbind(DD.results.DD20$mean.discoveries, DD.results.DEDD20$mean.discoveries)
names(DD.results.20$mean.discoveries) <- c('MDSeq.DD', 'HM.DD', 'MDSeq.DEDD', 'HM.DEDD')
DD.results.20$mean.discoveries <- DD.results.20$mean.discoveries[, c(1,3,2,4)]
DD.results.50$mean.discoveries <- cbind(DD.results.DD50$mean.discoveries, DD.results.DEDD50$mean.discoveries)
names(DD.results.50$mean.discoveries) <- c('MDSeq.DD', 'HM.DD', 'MDSeq.DEDD', 'HM.DEDD')
DD.results.50$mean.discoveries <- DD.results.50$mean.discoveries[, c(1,3,2,4)]
DD.results.2$mean.fdr <- cbind(DD.results.DD2$mean.fdr, DD.results.DEDD2$mean.fdr)
names(DD.results.2$mean.fdr) <- c('MDSeq.DD', 'HM.DD', 'MDSeq.DEDD', 'HM.DEDD')
DD.results.2$mean.fdr <- DD.results.2$mean.fdr[, c(1,3,2,4)]
DD.results.5$mean.fdr <- cbind(DD.results.DD5$mean.fdr, DD.results.DEDD5$mean.fdr)
names(DD.results.5$mean.fdr) <- c('MDSeq.DD', 'HM.DD', 'MDSeq.DEDD', 'HM.DEDD')
DD.results.5$mean.fdr <- DD.results.5$mean.fdr[, c(1,3,2,4)]
DD.results.10$mean.fdr <- cbind(DD.results.DD10$mean.fdr, DD.results.DEDD10$mean.fdr)
names(DD.results.10$mean.fdr) <- c('MDSeq.DD', 'HM.DD', 'MDSeq.DEDD', 'HM.DEDD')
DD.results.10$mean.fdr <- DD.results.10$mean.fdr[, c(1,3,2,4)]
DD.results.20$mean.fdr <- cbind(DD.results.DD20$mean.fdr, DD.results.DEDD20$mean.fdr)
names(DD.results.20$mean.fdr) <- c('MDSeq.DD', 'HM.DD', 'MDSeq.DEDD', 'HM.DEDD')
DD.results.20$mean.fdr <- DD.results.20$mean.fdr[, c(1,3,2,4)]
DD.results.50$mean.fdr <- cbind(DD.results.DD50$mean.fdr, DD.results.DEDD50$mean.fdr)
names(DD.results.50$mean.fdr) <- c('MDSeq.DD', 'HM.DD', 'MDSeq.DEDD', 'HM.DEDD')
DD.results.50$mean.fdr <- DD.results.50$mean.fdr[, c(1,3,2,4)]

par(mfrow=c(2,5), mar=c(2,2,3,1), mgp=c(3,0.7,0))
plot(DD.results.2$mean.discoveries$MDSeq.DD, 
     DD.results.2$mean.fdr$MDSeq.DD, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("2 samples per group\n", "Differences in dispersion only"))
lines(DD.results.2$mean.discoveries$HM.DD, 
      DD.results.2$mean.fdr$HM.DD, 
      type='l', col=col_vector[2])
legend("topright", bty='n', legend=c('MDSeq','HM'), col=col_vector[1:2], 
       lty=1, ncol=2)
abline(h=0.05, col='lightgrey')
plot(DD.results.5$mean.discoveries$MDSeq.DD, 
     DD.results.5$mean.fdr$MDSeq.DD, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("5 samples per group\n", "Differences in dispersion only"))
lines(DD.results.5$mean.discoveries$HM.DD, 
      DD.results.5$mean.fdr$HM.DD, 
      type='l', col=col_vector[2])
legend("topright", bty='n', legend=c('MDSeq','HM'), col=col_vector[1:2], 
       lty=1, ncol=2)
abline(h=0.05, col='lightgrey')
plot(DD.results.10$mean.discoveries$MDSeq.DD, 
     DD.results.10$mean.fdr$MDSeq.DD, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("10 samples per group\n", "Differences in dispersion only"))
lines(DD.results.10$mean.discoveries$HM.DD, 
      DD.results.10$mean.fdr$HM.DD, 
      type='l', col=col_vector[2])
legend("topright", bty='n', legend=c('MDSeq','HM'), col=col_vector[1:2], 
       lty=1, ncol=2)
abline(h=0.05, col='lightgrey')
plot(DD.results.20$mean.discoveries$MDSeq.DD, 
     DD.results.20$mean.fdr$MDSeq.DD, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("20 samples per group\n", "Differences in dispersion only"))
lines(DD.results.20$mean.discoveries$HM.DD, 
      DD.results.20$mean.fdr$HM.DD, 
      type='l', col=col_vector[2])
legend("topright", bty='n', legend=c('MDSeq','HM'), col=col_vector[1:2], 
       lty=1, ncol=2)
abline(h=0.05, col='lightgrey')
plot(DD.results.50$mean.discoveries$MDSeq.DD, 
     DD.results.50$mean.fdr$MDSeq.DD, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("50 samples per group\n", "Differences in dispersion only"))
lines(DD.results.50$mean.discoveries$HM.DD, 
      DD.results.50$mean.fdr$HM.DD, 
      type='l', col=col_vector[2])
legend("topright", bty='n', legend=c('MDSeq','HM'), col=col_vector[1:2], 
       lty=1, ncol=2)
abline(h=0.05, col='lightgrey')
plot(DD.results.2$mean.discoveries$MDSeq.DEDD, 
     DD.results.2$mean.fdr$MDSeq.DEDD, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("2 samples per group\n", "Differences in mean and dispersion"))
lines(DD.results.2$mean.discoveries$HM.DEDD, 
      DD.results.2$mean.fdr$HM.DEDD, 
      type='l', col=col_vector[2])
legend("topright", bty='n', legend=c('MDSeq','HM'), col=col_vector[1:2], 
       lty=1, ncol=2)
abline(h=0.05, col='lightgrey')
plot(DD.results.5$mean.discoveries$MDSeq.DEDD, 
     DD.results.5$mean.fdr$MDSeq.DEDD, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("5 samples per group\n", "Differences in mean and dispersion"))
lines(DD.results.5$mean.discoveries$HM.DEDD, 
      DD.results.5$mean.fdr$HM.DEDD, 
      type='l', col=col_vector[2])
legend("topright", bty='n', legend=c('MDSeq','HM'), col=col_vector[1:2], 
       lty=1, ncol=2)
abline(h=0.05, col='lightgrey')
plot(DD.results.10$mean.discoveries$MDSeq.DEDD, 
     DD.results.10$mean.fdr$MDSeq.DEDD, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("10 samples per group\n", "Differences in mean and dispersion"))
lines(DD.results.10$mean.discoveries$HM.DEDD, 
      DD.results.10$mean.fdr$HM.DEDD, 
      type='l', col=col_vector[2])
legend("topright", bty='n', legend=c('MDSeq','HM'), col=col_vector[1:2], 
       lty=1, ncol=2)
abline(h=0.05, col='lightgrey')
plot(DD.results.20$mean.discoveries$MDSeq.DEDD, 
     DD.results.20$mean.fdr$MDSeq.DEDD, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("20 samples per group\n", "Differences in mean and dispersion"))
lines(DD.results.20$mean.discoveries$HM.DEDD, 
      DD.results.20$mean.fdr$HM.DEDD, 
      type='l', col=col_vector[2])
legend("topright", bty='n', legend=c('MDSeq','HM'), col=col_vector[1:2], 
       lty=1, ncol=2)
abline(h=0.05, col='lightgrey')
plot(DD.results.50$mean.discoveries$MDSeq.DEDD, 
     DD.results.50$mean.fdr$MDSeq.DEDD, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("50 samples per group\n", "Differences in mean and dispersion"))
lines(DD.results.50$mean.discoveries$HM.DEDD, 
      DD.results.50$mean.fdr$HM.DEDD, 
      type='l', col=col_vector[2])
legend("topright", bty='n', legend=c('MDSeq','HM'), col=col_vector[1:2], 
       lty=1, ncol=2)
abline(h=0.05, col='lightgrey')


#################################
#### Differential expression ####
#################################
DE.results.2 <- list()
DE.results.5 <- list()
DE.results.10 <- list()
DE.results.20 <- list()
DE.results.50 <- list()

###########
### AUC ###
DE.results.2$auc <- cbind(DE.results.DE2$auc, DE.results.DEDD2$auc)
names(DE.results.2$auc) <- c('edgeR.DE', 'voom.DE', 'DESeq2.DE', 'DSS.DE', 
                             'baySeq.DE', 'MDSeq.DE', 'HM.DE', 
                             'edgeR.DEDD', 'voom.DEDD', 'DESeq2.DEDD', 'DSS.DEDD', 
                             'baySeq.DEDD', 'MDSeq.DEDD', 'HM.DEDD')
DE.results.2$auc <- DE.results.2$auc[, c(1,8,2,9,3,10,4,11,5,12,6,13,7,14)]
DE.results.5$auc <- cbind(DE.results.DE5$auc, DE.results.DEDD5$auc)
names(DE.results.5$auc) <- c('edgeR.DE', 'voom.DE', 'DESeq2.DE', 'DSS.DE', 
                             'baySeq.DE', 'HMDSeq.DE', 'HM.DE', 
                             'edgeR.DEDD', 'voom.DEDD', 'DESeq2.DEDD', 'DSS.DEDD', 
                             'baySeq.DEDD', 'HMDSeq.DEDD', 'HM.DEDD')
DE.results.5$auc <- DE.results.5$auc[, c(1,8,2,9,3,10,4,11,5,12,6,13,7,14)]
DE.results.10$auc <- cbind(DE.results.DE10$auc, DE.results.DEDD10$auc)
names(DE.results.10$auc) <- c('edgeR.DE', 'voom.DE', 'DESeq2.DE', 'DSS.DE', 
                              'baySeq.DE', 'HMDSeq.DE', 'HM.DE', 
                              'edgeR.DEDD', 'voom.DEDD', 'DESeq2.DEDD', 'DSS.DEDD', 
                              'baySeq.DEDD', 'HMDSeq.DEDD', 'HM.DEDD')
DE.results.10$auc <- DE.results.10$auc[, c(1,8,2,9,3,10,4,11,5,12,6,13,7,14)]
DE.results.20$auc <- cbind(DE.results.DE20$auc, DE.results.DEDD20$auc)
names(DE.results.20$auc) <- c('edgeR.DE', 'voom.DE', 'DESeq2.DE', 'DSS.DE', 
                              'baySeq.DE', 'HMDSeq.DE', 'HM.DE', 
                              'edgeR.DEDD', 'voom.DEDD', 'DESeq2.DEDD', 'DSS.DEDD', 
                              'baySeq.DEDD', 'HMDSeq.DEDD', 'HM.DEDD')
DE.results.20$auc <- DE.results.20$auc[, c(1,8,2,9,3,10,4,11,5,12,6,13,7,14)]
DE.results.50$auc <- cbind(DE.results.DE50$auc, DE.results.DEDD50$auc)
names(DE.results.50$auc) <- c('edgeR.DE', 'voom.DE', 'DESeq2.DE', 
                              'baySeq.DE', 'HMDSeq.DE', 'HM.DE', 
                              'edgeR.DEDD', 'voom.DEDD', 'DESeq2.DEDD', 
                              'baySeq.DEDD', 'HMDSeq.DEDD', 'HM.DEDD')
DE.results.50$auc <- DE.results.50$auc[, c(1,7,2,8,3,9,4,10,5,11,6,12)]

par(mfrow=c(1,2), mar=c(1,2.5,4,1), mgp=c(3,0.7,0))
boxplot(DE.results.2$auc[c(1,3,5,7,9,11,13)], 
        names=NA, col=col_vector[1:7], 
        xaxt='n', ylim=c(0.75,1), cex.axis=2, 
        pch=20, main=paste0("2 samples per group"), 
        cex.main=2)
# legend("topleft", fill=col_vector[1:7], bty='n', cex=1.5, ncol=4, 
#        legend=c("edgeR", "voom", "DESeq2", "DSS", "baySeq", "MDSeq", "HM"))
# boxplot(DE.results.5$auc[c(1,3,5,7,9,11,13)], 
#         names=NA, col=col_vector[1:7], 
#         xaxt='n', ylim=c(0.75,1), 
#         pch=20, main=paste0("AUC, 5 samples per group\n", 
#                             "Differences in mean only"), 
#         cex.main=1.5)
# legend("topleft", fill=col_vector[1:7],bty='n', cex=1.5, ncol=4, 
#        legend=c("edgeR", "voom", "DESeq2", "DSS", "baySeq", "MDSeq", "HM"))
boxplot(DE.results.10$auc[c(1,3,5,7,9,11,13)], 
        names=NA, col=col_vector[1:7], 
        xaxt='n', ylim=c(0.75,1), cex.axis=2, 
        pch=20, main=paste0("10 samples per group"), 
        cex.main=2)
# legend("bottomleft", fill=col_vector[1:7], bty='n', cex=1.5, ncol=4, 
#        legend=c("edgeR", "voom", "DESeq2", "DSS", "baySeq", "MDSeq", "HM"))
# boxplot(DE.results.2$auc[c(2,4,6,8,10,12,14)], 
#         names=NA, col=col_vector[1:7], 
#         xaxt='n', ylim=c(0.75,1), 
#         pch=20, main=paste0("AUC, 2 samples per group\n", 
#                             "Differences in mean and dispersion"), 
#         cex.main=1.5)
# legend("topleft", fill=col_vector[1:7], bty='n', cex=1.5, ncol=4, 
#        legend=c("edgeR", "voom", "DESeq2", "DSS", "baySeq", "MDSeq", "HM"))
# boxplot(DE.results.5$auc[c(2,4,6,8,10,12,14)], 
#         names=NA, col=col_vector[1:7], 
#         xaxt='n', ylim=c(0.75,1), 
#         pch=20, main=paste0("AUC, 5 samples per group\n", 
#                             "Differences in mean and dispersion"), 
#         cex.main=1.5)
# legend("topleft", fill=col_vector[1:7],bty='n', cex=1.5, ncol=4, 
#        legend=c("edgeR", "voom", "DESeq2", "DSS", "baySeq", "MDSeq", "HM"))
# boxplot(DE.results.10$auc[c(2,4,6,8,10,12,14)], 
#         names=NA, col=col_vector[1:7], 
#         xaxt='n', ylim=c(0.75,1), 
#         pch=20, main=paste0("AUC, 10 samples per group\n", 
#                             "Differences in mean and dispersion"), 
#         cex.main=1.5)
# legend("bottomleft", fill=col_vector[1:7], bty='n', cex=1.5, ncol=4, 
#        legend=c("edgeR", "voom", "DESeq2", "DSS", "baySeq", "MDSeq", "HM"))
# 
# 
###########
### pAUC ###
DE.results.2$pauc <- cbind(DE.results.DE2$pauc, DE.results.DEDD2$pauc)
names(DE.results.2$pauc) <- c('edgeR.DE', 'voom.DE', 'DESeq2.DE', 'DSS.DE', 
                             'baySeq.DE', 'MDSeq.DE', 'HM.DE', 
                             'edgeR.DEDD', 'voom.DEDD', 'DESeq2.DEDD', 'DSS.DEDD', 
                             'baySeq.DEDD', 'MDSeq.DEDD', 'HM.DEDD')
DE.results.2$pauc <- DE.results.2$pauc[, c(1,8,2,9,3,10,4,11,5,12,6,13,7,14)]
DE.results.5$pauc <- cbind(DE.results.DE5$pauc, DE.results.DEDD5$pauc)
names(DE.results.5$pauc) <- c('edgeR.DE', 'voom.DE', 'DESeq2.DE', 'DSS.DE', 
                             'baySeq.DE', 'HMDSeq.DE', 'HM.DE', 
                             'edgeR.DEDD', 'voom.DEDD', 'DESeq2.DEDD', 'DSS.DEDD', 
                             'baySeq.DEDD', 'HMDSeq.DEDD', 'HM.DEDD')
DE.results.5$pauc <- DE.results.5$pauc[, c(1,8,2,9,3,10,4,11,5,12,6,13,7,14)]
DE.results.10$pauc <- cbind(DE.results.DE10$pauc, DE.results.DEDD10$pauc)
names(DE.results.10$pauc) <- c('edgeR.DE', 'voom.DE', 'DESeq2.DE', 'DSS.DE', 
                              'baySeq.DE', 'HMDSeq.DE', 'HM.DE', 
                              'edgeR.DEDD', 'voom.DEDD', 'DESeq2.DEDD', 'DSS.DEDD', 
                              'baySeq.DEDD', 'HMDSeq.DEDD', 'HM.DEDD')
DE.results.10$pauc <- DE.results.10$pauc[, c(1,8,2,9,3,10,4,11,5,12,6,13,7,14)]
DE.results.20$pauc <- cbind(DE.results.DE20$pauc, DE.results.DEDD20$pauc)
names(DE.results.20$pauc) <- c('edgeR.DE', 'voom.DE', 'DESeq2.DE', 'DSS.DE', 
                              'baySeq.DE', 'HMDSeq.DE', 'HM.DE', 
                              'edgeR.DEDD', 'voom.DEDD', 'DESeq2.DEDD', 'DSS.DEDD', 
                              'baySeq.DEDD', 'HMDSeq.DEDD', 'HM.DEDD')
DE.results.20$pauc <- DE.results.20$pauc[, c(1,8,2,9,3,10,4,11,5,12,6,13,7,14)]
DE.results.50$pauc <- cbind(DE.results.DE50$pauc, DE.results.DEDD50$pauc)
names(DE.results.50$pauc) <- c('edgeR.DE', 'voom.DE', 'DESeq2.DE', 
                               'baySeq.DE', 'HMDSeq.DE', 'HM.DE', 
                               'edgeR.DEDD', 'voom.DEDD', 'DESeq2.DEDD', 
                               'baySeq.DEDD', 'HMDSeq.DEDD', 'HM.DEDD')
DE.results.50$pauc <- DE.results.50$pauc[, c(1,7,2,8,3,9,4,10,5,11,6,12)]

par(mfrow=c(2,3), mar=c(1,2,3,1), mgp=c(3,0.7,0))
boxplot(DE.results.2$pauc[c(1,3,5,7,9,11,13)], 
        names=NA, col=col_vector[1:7], 
        xaxt='n', ylim=c(0,0.04), 
        pch=20, main=paste0("Partial AUC, 2 samples per group\n", 
                            "Differences in mean only"), 
        cex.main=1.5)
legend("topleft", fill=col_vector[1:7], bty='n', cex=1.5, ncol=4, 
       legend=c("edgeR", "voom", "DESeq2", "DSS", "baySeq", "MDSeq", "HM"))
boxplot(DE.results.5$pauc[c(1,3,5,7,9,11,13)], 
        names=NA, col=col_vector[1:7], 
        xaxt='n', ylim=c(0,0.04), 
        pch=20, main=paste0("Partial AUC, 5 samples per group\n", 
                            "Differences in mean only"), 
        cex.main=1.5)
legend("topleft", fill=col_vector[1:7],bty='n', cex=1.5, ncol=4, 
       legend=c("edgeR", "voom", "DESeq2", "DSS", "baySeq", "MDSeq", "HM"))
boxplot(DE.results.10$pauc[c(1,3,5,7,9,11,13)], 
        names=NA, col=col_vector[1:7], 
        xaxt='n', ylim=c(0,0.04), 
        pch=20, main=paste0("Partial AUC, 10 samples per group\n", 
                            "Differences in mean only"), 
        cex.main=1.5)
legend("bottomleft", fill=col_vector[1:7], bty='n', cex=1.5, ncol=4, 
       legend=c("edgeR", "voom", "DESeq2", "DSS", "baySeq", "MDSeq", "HM"))
boxplot(DE.results.2$pauc[c(2,4,6,8,10,12,14)], 
        names=NA, col=col_vector[1:7], 
        xaxt='n', ylim=c(0,0.04), 
        pch=20, main=paste0("Partial AUC, 2 samples per group\n", 
                            "Differences in mean and dispersion"), 
        cex.main=1.5)
legend("topleft", fill=col_vector[1:7], bty='n', cex=1.5, ncol=4, 
       legend=c("edgeR", "voom", "DESeq2", "DSS", "baySeq", "MDSeq", "HM"))
boxplot(DE.results.5$pauc[c(2,4,6,8,10,12,14)], 
        names=NA, col=col_vector[1:7], 
        xaxt='n', ylim=c(0,0.04), 
        pch=20, main=paste0("Partial AUC, 5 samples per group\n", 
                            "Differences in mean and dispersion"), 
        cex.main=1.5)
legend("topleft", fill=col_vector[1:7],bty='n', cex=1.5, ncol=4, 
       legend=c("edgeR", "voom", "DESeq2", "DSS", "baySeq", "MDSeq", "HM"))
boxplot(DE.results.10$pauc[c(2,4,6,8,10,12,14)], 
        names=NA, col=col_vector[1:7], 
        xaxt='n', ylim=c(0,0.04), 
        pch=20, main=paste0("Partial AUC, 10 samples per group\n", 
                            "Differences in mean and dispersion"), 
        cex.main=1.5)
legend("bottomleft", fill=col_vector[1:7], bty='n', cex=1.5, ncol=4, 
       legend=c("edgeR", "voom", "DESeq2", "DSS", "baySeq", "MDSeq", "HM"))

###########
### FDR ###
DE.results.2$fdr <- cbind(DE.results.DE2$fdr, DE.results.DEDD2$fdr)
names(DE.results.2$fdr) <- c('edgeR.DE', 'voom.DE', 'DESeq2.DE', 'DSS.DE', 
                              'baySeq.DE', 'MDSeq.DE', 'HM.DE', 
                              'edgeR.DEDD', 'voom.DEDD', 'DESeq2.DEDD', 'DSS.DEDD', 
                              'baySeq.DEDD', 'MDSeq.DEDD', 'HM.DEDD')
DE.results.2$fdr <- DE.results.2$fdr[, c(1,8,2,9,3,10,4,11,5,12,6,13,7,14)]
DE.results.5$fdr <- cbind(DE.results.DE5$fdr, DE.results.DEDD5$fdr)
names(DE.results.5$fdr) <- c('edgeR.DE', 'voom.DE', 'DESeq2.DE', 'DSS.DE', 
                              'baySeq.DE', 'HMDSeq.DE', 'HM.DE', 
                              'edgeR.DEDD', 'voom.DEDD', 'DESeq2.DEDD', 'DSS.DEDD', 
                              'baySeq.DEDD', 'HMDSeq.DEDD', 'HM.DEDD')
DE.results.5$fdr <- DE.results.5$fdr[, c(1,8,2,9,3,10,4,11,5,12,6,13,7,14)]
DE.results.10$fdr <- cbind(DE.results.DE10$fdr, DE.results.DEDD10$fdr)
names(DE.results.10$fdr) <- c('edgeR.DE', 'voom.DE', 'DESeq2.DE', 'DSS.DE', 
                               'baySeq.DE', 'HMDSeq.DE', 'HM.DE', 
                               'edgeR.DEDD', 'voom.DEDD', 'DESeq2.DEDD', 'DSS.DEDD', 
                               'baySeq.DEDD', 'HMDSeq.DEDD', 'HM.DEDD')
DE.results.10$fdr <- DE.results.10$fdr[, c(1,8,2,9,3,10,4,11,5,12,6,13,7,14)]
DE.results.20$fdr <- cbind(DE.results.DE20$fdr, DE.results.DEDD20$fdr)
names(DE.results.20$fdr) <- c('edgeR.DE', 'voom.DE', 'DESeq2.DE', 'DSS.DE', 
                               'baySeq.DE', 'HMDSeq.DE', 'HM.DE', 
                               'edgeR.DEDD', 'voom.DEDD', 'DESeq2.DEDD', 'DSS.DEDD', 
                               'baySeq.DEDD', 'HMDSeq.DEDD', 'HM.DEDD')
DE.results.20$fdr <- DE.results.20$fdr[, c(1,8,2,9,3,10,4,11,5,12,6,13,7,14)]
DE.results.50$fdr <- cbind(DE.results.DE50$fdr, DE.results.DEDD50$fdr)
names(DE.results.50$fdr) <- c('edgeR.DE', 'voom.DE', 'DESeq2.DE', 
                              'baySeq.DE', 'HMDSeq.DE', 'HM.DE', 
                              'edgeR.DEDD', 'voom.DEDD', 'DESeq2.DEDD', 
                              'baySeq.DEDD', 'HMDSeq.DEDD', 'HM.DEDD')
DE.results.50$fdr <- DE.results.50$fdr[, c(1,7,2,8,3,9,4,10,5,11,6,12)]

par(mfrow=c(1,2), mar=c(1,2.5,4,1), mgp=c(3,0.7,0))
boxplot(DE.results.2$fdr[c(1,3,5,7,9,11,13)], 
        names=NA, col=col_vector[1:7], 
        xaxt='n', ylim=c(0,0.8), cex.axis=2, 
        pch=20, main=paste0("2 samples per group"), 
        cex.main=2)
abline(h=0.05, col='lightgrey')
# legend("topleft", fill=col_vector[1:7], bty='n', cex=1.5, ncol=4, 
#        legend=c("edgeR", "voom", "DESeq2", "DSS", "baySeq", "MDSeq", "HM"))
# boxplot(DE.results.5$fdr[c(1,3,5,7,9,11,13)], 
#         names=NA, col=col_vector[1:7], 
#         xaxt='n', ylim=c(0,1), 
#         pch=20, main=paste0("FDR, 5 samples per group\n", 
#                             "Differences in mean only"), 
#         cex.main=1.5)
# abline(h=0.05, col='lightgrey')
# legend("topleft", fill=col_vector[1:7],bty='n', cex=1.5, ncol=4, 
#        legend=c("edgeR", "voom", "DESeq2", "DSS", "baySeq", "MDSeq", "HM"))
boxplot(DE.results.10$fdr[c(1,3,5,7,9,11,13)], 
        names=NA, col=col_vector[1:7], 
        xaxt='n', ylim=c(0,0.8), cex.axis=2, 
        pch=20, main=paste0("10 samples per group"), 
        cex.main=2)
abline(h=0.05, col='grey')
# legend("topleft", fill=col_vector[1:7], bty='n', cex=1.5, ncol=4, 
#        legend=c("edgeR", "voom", "DESeq2", "DSS", "baySeq", "MDSeq", "HM"))
# boxplot(DE.results.2$fdr[c(2,4,6,8,10,12,14)], 
#         names=NA, col=col_vector[1:7], 
#         xaxt='n', ylim=c(0,1), 
#         pch=20, main=paste0("FDR, 2 samples per group\n", 
#                             "Differences in mean and dispersion"), 
#         cex.main=1.5)
# abline(h=0.05, col='lightgrey')
# legend("topleft", fill=col_vector[1:7], bty='n', cex=1.5, ncol=4, 
#        legend=c("edgeR", "voom", "DESeq2", "DSS", "baySeq", "MDSeq", "HM"))
# boxplot(DE.results.5$fdr[c(2,4,6,8,10,12,14)], 
#         names=NA, col=col_vector[1:7], 
#         xaxt='n', ylim=c(0,1), 
#         pch=20, main=paste0("FDR, 5 samples per group\n", 
#                             "Differences in mean and dispersion"), 
#         cex.main=1.5)
# abline(h=0.05, col='lightgrey')
# legend("topleft", fill=col_vector[1:7],bty='n', cex=1.5, ncol=4, 
#        legend=c("edgeR", "voom", "DESeq2", "DSS", "baySeq", "MDSeq", "HM"))
# boxplot(DE.results.10$fdr[c(2,4,6,8,10,12,14)], 
#         names=NA, col=col_vector[1:7], 
#         xaxt='n', ylim=c(0,1), 
#         pch=20, main=paste0("FDR, 10 samples per group\n", 
#                             "Differences in mean and dispersion"), 
#         cex.main=1.5)
# abline(h=0.05, col='lightgrey')
# legend("topleft", fill=col_vector[1:7], bty='n', cex=1.5, ncol=4, 
#        legend=c("edgeR", "voom", "DESeq2", "DSS", "baySeq", "MDSeq", "HM"))

###########
### TPR ###
DE.results.2$tpr <- cbind(DE.results.DE2$tpr, DE.results.DEDD2$tpr)
names(DE.results.2$tpr) <- c('edgeR.DE', 'voom.DE', 'DESeq2.DE', 'DSS.DE', 
                             'baySeq.DE', 'DSeq.DE', 'HM.DE', 
                             'edgeR.DEDD', 'voom.DEDD', 'DESeq2.DEDD', 'DSS.DEDD', 
                             'baySeq.DEDD', 'DSeq.DEDD', 'HM.DEDD')
DE.results.2$tpr <- DE.results.2$tpr[, c(1,8,2,9,3,10,4,11,5,12,6,13,7,14)]
DE.results.5$tpr <- cbind(DE.results.DE5$tpr, DE.results.DEDD5$tpr)
names(DE.results.5$tpr) <- c('edgeR.DE', 'voom.DE', 'DESeq2.DE', 'DSS.DE', 
                             'baySeq.DE', 'MDSeq.DE', 'HM.DE', 
                             'edgeR.DEDD', 'voom.DEDD', 'DESeq2.DEDD', 'DSS.DEDD', 
                             'baySeq.DEDD', 'MDSeq.DEDD', 'HM.DEDD')
DE.results.5$tpr <- DE.results.5$tpr[, c(1,8,2,9,3,10,4,11,5,12,6,13,7,14)]
DE.results.10$tpr <- cbind(DE.results.DE10$tpr, DE.results.DEDD10$tpr)
names(DE.results.10$tpr) <- c('edgeR.DE', 'voom.DE', 'DESeq2.DE', 'DSS.DE', 
                              'baySeq.DE', 'MDSeq.DE', 'HM.DE', 
                              'edgeR.DEDD', 'voom.DEDD', 'DESeq2.DEDD', 'DSS.DEDD', 
                              'baySeq.DEDD', 'MDSeq.DEDD', 'HM.DEDD')
DE.results.10$tpr <- DE.results.10$tpr[, c(1,8,2,9,3,10,4,11,5,12,6,13,7,14)]
DE.results.20$tpr <- cbind(DE.results.DE20$tpr, DE.results.DEDD20$tpr)
names(DE.results.20$tpr) <- c('edgeR.DE', 'voom.DE', 'DESeq2.DE', 'DSS.DE', 
                              'baySeq.DE', 'MDSeq.DE', 'HM.DE', 
                              'edgeR.DEDD', 'voom.DEDD', 'DESeq2.DEDD', 'DSS.DEDD', 
                              'baySeq.DEDD', 'MDSeq.DEDD', 'HM.DEDD')
DE.results.20$tpr <- DE.results.20$tpr[, c(1,8,2,9,3,10,4,11,5,12,6,13,7,14)]
DE.results.50$tpr <- cbind(DE.results.DE50$tpr, DE.results.DEDD50$tpr)
names(DE.results.50$tpr) <- c('edgeR.DE', 'voom.DE', 'DESeq2.DE', 
                              'baySeq.DE', 'MDSeq.DE', 'HM.DE', 
                              'edgeR.DEDD', 'voom.DEDD', 'DESeq2.DEDD', 
                              'baySeq.DEDD', 'MDSeq.DEDD', 'HM.DEDD')
DE.results.50$tpr <- DE.results.50$tpr[, c(1,7,2,8,3,9,4,10,5,11,6,12)]

par(mfrow=c(2,3), mar=c(1,2,3,1), mgp=c(3,0.7,0))
boxplot(DE.results.2$tpr[c(1,3,5,7,9,11,13)], 
        names=NA, col=col_vector[1:7], 
        xaxt='n', ylim=c(0,0.75), 
        pch=20, main=paste0("TPR, 2 samples per group\n", 
                            "Differences in mean only"), 
        cex.main=1.5)
legend("topleft", fill=col_vector[1:7], bty='n', cex=1.5, ncol=4, 
       legend=c("edgeR", "voom", "DESeq2", "DSS", "baySeq", "MDSeq", "HM"))
boxplot(DE.results.5$tpr[c(1,3,5,7,9,11,13)], 
        names=NA, col=col_vector[1:7], 
        xaxt='n', ylim=c(0,0.75), 
        pch=20, main=paste0("TPR, 5 samples per group\n", 
                            "Differences in mean only"), 
        cex.main=1.5)
legend("topleft", fill=col_vector[1:7],bty='n', cex=1.5, ncol=4, 
       legend=c("edgeR", "voom", "DESeq2", "DSS", "baySeq", "MDSeq", "HM"))
boxplot(DE.results.10$tpr[c(1,3,5,7,9,11,13)], 
        names=NA, col=col_vector[1:7], 
        xaxt='n', ylim=c(0,0.75), 
        pch=20, main=paste0("TPR, 10 samples per group\n", 
                            "Differences in mean only"), 
        cex.main=1.5)
legend("bottomright", fill=col_vector[1:7], bty='n', cex=1.5, ncol=4, 
       legend=c("edgeR", "voom", "DESeq2", "DSS", "baySeq", "MDSeq", "HM"))
boxplot(DE.results.2$tpr[c(2,4,6,8,10,12,14)], 
        names=NA, col=col_vector[1:7], 
        xaxt='n', ylim=c(0,0.75), 
        pch=20, main=paste0("TPR, 2 samples per group\n", 
                            "Differences in mean and dispersion"), 
        cex.main=1.5)
legend("topleft", fill=col_vector[1:7], bty='n', cex=1.5, ncol=4, 
       legend=c("edgeR", "voom", "DESeq2", "DSS", "baySeq", "MDSeq", "HM"))
boxplot(DE.results.5$tpr[c(2,4,6,8,10,12,14)], 
        names=NA, col=col_vector[1:7], 
        xaxt='n', ylim=c(0,0.75), 
        pch=20, main=paste0("TPR, 5 samples per group\n", 
                            "Differences in mean and dispersion"), 
        cex.main=1.5)
legend("topleft", fill=col_vector[1:7],bty='n', cex=1.5, ncol=4, 
       legend=c("edgeR", "voom", "DESeq2", "DSS", "baySeq", "MDSeq", "HM"))
boxplot(DE.results.10$tpr[c(2,4,6,8,10,12,14)], 
        names=NA, col=col_vector[1:7], 
        xaxt='n', ylim=c(0,0.75), 
        pch=20, main=paste0("TPR, 10 samples per group\n", 
                            "Differences in mean and dispersion"), 
        cex.main=1.5)
legend("bottomright", fill=col_vector[1:7], bty='n', cex=1.5, ncol=4, 
       legend=c("edgeR", "voom", "DESeq2", "DSS", "baySeq", "MDSeq", "HM"))


#############################
### False discovery plots ###
DE.results.2$mean.discoveries <- cbind(DE.results.DE2$mean.discoveries, DE.results.DEDD2$mean.discoveries)
names(DE.results.2$mean.discoveries) <- c('edgeR.DE', 'voom.DE', 'DESeq2.DE', 'DSS.DE', 
                                          'baySeq.DE', 'MDSeq.DE', 'HM.DE', 
                                          'edgeR.DEDD', 'voom.DEDD', 'DESeq2.DEDD', 'DSS.DEDD', 
                                          'baySeq.DEDD', 'MDSeq.DEDD', 'HM.DEDD')
DE.results.2$mean.discoveries <- DE.results.2$mean.discoveries[, c(1,8,2,9,3,10,4,11,5,12,6,13,7,14)]
DE.results.5$mean.discoveries <- cbind(DE.results.DE5$mean.discoveries, DE.results.DEDD5$mean.discoveries)
names(DE.results.5$mean.discoveries) <- c('edgeR.DE', 'voom.DE', 'DESeq2.DE', 'DSS.DE', 
                                          'baySeq.DE', 'MDSeq.DE', 'HM.DE', 
                                          'edgeR.DEDD', 'voom.DEDD', 'DESeq2.DEDD', 'DSS.DEDD', 
                                          'baySeq.DEDD', 'MDSeq.DEDD', 'HM.DEDD')
DE.results.5$mean.discoveries <- DE.results.5$mean.discoveries[, c(1,8,2,9,3,10,4,11,5,12,6,13,7,14)]
DE.results.10$mean.discoveries <- cbind(DE.results.DE10$mean.discoveries, DE.results.DEDD10$mean.discoveries)
names(DE.results.10$mean.discoveries) <- c('edgeR.DE', 'voom.DE', 'DESeq2.DE', 'DSS.DE', 
                                           'baySeq.DE', 'MDSeq.DE', 'HM.DE', 
                                           'edgeR.DEDD', 'voom.DEDD', 'DESeq2.DEDD', 'DSS.DEDD', 
                                           'baySeq.DEDD', 'MDSeq.DEDD', 'HM.DEDD')
DE.results.10$mean.discoveries <- DE.results.10$mean.discoveries[, c(1,8,2,9,3,10,4,11,5,12,6,13,7,14)]
DE.results.20$mean.discoveries <- cbind(DE.results.DE20$mean.discoveries, DE.results.DEDD20$mean.discoveries)
names(DE.results.20$mean.discoveries) <- c('edgeR.DE', 'voom.DE', 'DESeq2.DE', 'DSS.DE', 
                                           'baySeq.DE', 'MDSeq.DE', 'HM.DE', 
                                           'edgeR.DEDD', 'voom.DEDD', 'DESeq2.DEDD', 'DSS.DEDD', 
                                           'baySeq.DEDD', 'MDSeq.DEDD', 'HM.DEDD')
DE.results.20$mean.discoveries <- DE.results.20$mean.discoveries[, c(1,8,2,9,3,10,4,11,5,12,6,13,7,14)]
DE.results.50$mean.discoveries <- cbind(DE.results.DE50$mean.discoveries, DE.results.DEDD50$mean.discoveries)
names(DE.results.50$mean.discoveries) <- c('edgeR.DE', 'voom.DE', 'DESeq2.DE', 
                                           'baySeq.DE', 'MDSeq.DE', 'HM.DE', 
                                           'edgeR.DEDD', 'voom.DEDD', 'DESeq2.DEDD', 
                                           'baySeq.DEDD', 'MDSeq.DEDD', 'HM.DEDD')
DE.results.50$mean.discoveries <- DE.results.50$mean.discoveries[, c(1,7,2,8,3,9,4,10,5,11,6,12)]
DE.results.2$mean.fdr <- cbind(DE.results.DE2$mean.fdr, DE.results.DEDD2$mean.fdr)
names(DE.results.2$mean.fdr) <- c('edgeR.DE', 'voom.DE', 'DESeq2.DE', 'DSS.DE', 
                                  'baySeq.DE', 'MDSeq.DE', 'HM.DE', 
                                  'edgeR.DEDD', 'voom.DEDD', 'DESeq2.DEDD', 'DSS.DEDD', 
                                  'baySeq.DEDD', 'MDSeq.DEDD', 'HM.DEDD')
DE.results.2$mean.fdr <- DE.results.2$mean.fdr[, c(1,8,2,9,3,10,4,11,5,12,6,13,7,14)]
DE.results.5$mean.fdr <- cbind(DE.results.DE5$mean.fdr, DE.results.DEDD5$mean.fdr)
names(DE.results.5$mean.fdr) <- c('edgeR.DE', 'voom.DE', 'DESeq2.DE', 'DSS.DE', 
                                  'baySeq.DE', 'MDSeq.DE', 'HM.DE', 
                                  'edgeR.DEDD', 'voom.DEDD', 'DESeq2.DEDD', 'DSS.DEDD', 
                                  'baySeq.DEDD', 'MDSeq.DEDD', 'HM.DEDD')
DE.results.5$mean.fdr <- DE.results.5$mean.fdr[, c(1,8,2,9,3,10,4,11,5,12,6,13,7,14)]
DE.results.10$mean.fdr <- cbind(DE.results.DE10$mean.fdr, DE.results.DEDD10$mean.fdr)
names(DE.results.10$mean.fdr) <- c('edgeR.DE', 'voom.DE', 'DESeq2.DE', 'DSS.DE', 
                                   'baySeq.DE', 'MDSeq.DE', 'HM.DE', 
                                   'edgeR.DEDD', 'voom.DEDD', 'DESeq2.DEDD', 'DSS.DEDD', 
                                   'baySeq.DEDD', 'MDSeq.DEDD', 'HM.DEDD')
DE.results.10$mean.fdr <- DE.results.10$mean.fdr[, c(1,8,2,9,3,10,4,11,5,12,6,13,7,14)]
DE.results.20$mean.fdr <- cbind(DE.results.DE20$mean.fdr, DE.results.DEDD20$mean.fdr)
names(DE.results.20$mean.fdr) <- c('edgeR.DE', 'voom.DE', 'DESeq2.DE', 'DSS.DE', 
                                   'baySeq.DE', 'MDSeq.DE', 'HM.DE', 
                                   'edgeR.DEDD', 'voom.DEDD', 'DESeq2.DEDD', 'DSS.DEDD', 
                                   'baySeq.DEDD', 'MDSeq.DEDD', 'HM.DEDD')
DE.results.20$mean.fdr <- DE.results.20$mean.fdr[, c(1,8,2,9,3,10,4,11,5,12,6,13,7,14)]
DE.results.50$mean.fdr <- cbind(DE.results.DE50$mean.fdr, DE.results.DEDD50$mean.fdr)
names(DE.results.50$mean.fdr) <- c('edgeR.DE', 'voom.DE', 'DESeq2.DE', 
                                   'baySeq.DE', 'MDSeq.DE', 'HM.DE', 
                                   'edgeR.DEDD', 'voom.DEDD', 'DESeq2.DEDD', 
                                   'baySeq.DEDD', 'MDSeq.DEDD', 'HM.DEDD')
DE.results.50$mean.fdr <- DE.results.50$mean.fdr[, c(1,7,2,8,3,9,4,10,5,11,6,12)]

par(mfrow=c(2,5), mar=c(2,2,3,1), mgp=c(3,0.7,0))
plot(DE.results.2$mean.discoveries$edgeR.DE, 
     DE.results.2$mean.fdr$edgeR.DE, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("2 samples per group\n", "Differences in mean only"))
lines(DE.results.2$mean.discoveries$voom.DE, DE.results.2$mean.fdr$voom.DE, 
      type='l', col=col_vector[2])
lines(DE.results.2$mean.discoveries$DESeq2.DE, DE.results.2$mean.fdr$DESeq2.DE, 
      type='l', col=col_vector[3])
lines(DE.results.2$mean.discoveries$DSS.DE, DE.results.2$mean.fdr$DSS.DE, 
      type='l', col=col_vector[4])
lines(DE.results.2$mean.discoveries$baySeq.DE, DE.results.2$mean.fdr$baySeq.DE, 
      type='l', col=col_vector[5])
lines(DE.results.2$mean.discoveries$MDSeq.DE, DE.results.2$mean.fdr$MDSeq.DE, 
      type='l', col=col_vector[6])
lines(DE.results.2$mean.discoveries$HM.DE, DE.results.2$mean.fdr$HM.DE, 
      type='l', col=col_vector[7])
legend("topright", bty='n', col=col_vector[1:7], 
       legend=c('edgeR', 'voom', 'DESeq2', 'DSS', 'baySeq', 'MDSeq','HM'), 
       lty=1, ncol=2)
abline(h=0.05, col='lightgrey')
plot(DE.results.5$mean.discoveries$edgeR.DE, 
     DE.results.5$mean.fdr$edgeR.DE, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("5 samples per group\n", "Differences in mean only"))
lines(DE.results.5$mean.discoveries$voom.DE, DE.results.5$mean.fdr$voom.DE, 
      type='l', col=col_vector[2])
lines(DE.results.5$mean.discoveries$DESeq2.DE, DE.results.5$mean.fdr$DESeq2.DE, 
      type='l', col=col_vector[3])
lines(DE.results.5$mean.discoveries$DSS.DE, DE.results.5$mean.fdr$DSS.DE, 
      type='l', col=col_vector[4])
lines(DE.results.5$mean.discoveries$baySeq.DE, DE.results.5$mean.fdr$baySeq.DE, 
      type='l', col=col_vector[5])
lines(DE.results.5$mean.discoveries$MDSeq.DE, DE.results.5$mean.fdr$MDSeq.DE, 
      type='l', col=col_vector[6])
lines(DE.results.5$mean.discoveries$HM.DE, DE.results.5$mean.fdr$HM.DE, 
      type='l', col=col_vector[7])
legend("topright", bty='n', col=col_vector[1:7], 
       legend=c('edgeR', 'voom', 'DESeq2', 'DSS', 'baySeq', 'MDSeq','HM'), 
       lty=1, ncol=2)
abline(h=0.05, col='lightgrey')
plot(DE.results.10$mean.discoveries$edgeR.DE, 
     DE.results.10$mean.fdr$edgeR.DE, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("10 samples per group\n", "Differences in mean only"))
lines(DE.results.10$mean.discoveries$voom.DE, DE.results.10$mean.fdr$voom.DE, 
      type='l', col=col_vector[2])
lines(DE.results.10$mean.discoveries$DESeq2.DE, DE.results.10$mean.fdr$DESeq2.DE, 
      type='l', col=col_vector[3])
lines(DE.results.10$mean.discoveries$DSS.DE, DE.results.10$mean.fdr$DSS.DE, 
      type='l', col=col_vector[4])
lines(DE.results.10$mean.discoveries$baySeq.DE, DE.results.10$mean.fdr$baySeq.DE, 
      type='l', col=col_vector[5])
lines(DE.results.10$mean.discoveries$MDSeq.DE, DE.results.10$mean.fdr$MDSeq.DE, 
      type='l', col=col_vector[6])
lines(DE.results.10$mean.discoveries$HM.DE, DE.results.10$mean.fdr$HM.DE, 
      type='l', col=col_vector[7])
legend("topright", bty='n', col=col_vector[1:7], 
       legend=c('edgeR', 'voom', 'DESeq2', 'DSS', 'baySeq', 'MDSeq','HM'), 
       lty=1, ncol=2)
abline(h=0.05, col='lightgrey')
plot(DE.results.20$mean.discoveries$edgeR.DE, 
     DE.results.20$mean.fdr$edgeR.DE, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("20 samples per group\n", "Differences in mean only"))
lines(DE.results.20$mean.discoveries$voom.DE, DE.results.20$mean.fdr$voom.DE, 
      type='l', col=col_vector[2])
lines(DE.results.20$mean.discoveries$DESeq2.DE, DE.results.20$mean.fdr$DESeq2.DE, 
      type='l', col=col_vector[3])
lines(DE.results.20$mean.discoveries$DSS.DE, DE.results.20$mean.fdr$DSS.DE, 
      type='l', col=col_vector[4])
lines(DE.results.20$mean.discoveries$baySeq.DE, DE.results.20$mean.fdr$baySeq.DE, 
      type='l', col=col_vector[5])
lines(DE.results.20$mean.discoveries$MDSeq.DE, DE.results.20$mean.fdr$MDSeq.DE, 
      type='l', col=col_vector[6])
lines(DE.results.20$mean.discoveries$HM.DE, DE.results.20$mean.fdr$HM.DE, 
      type='l', col=col_vector[7])
legend("topright", bty='n', col=col_vector[1:7], 
       legend=c('edgeR', 'voom', 'DESeq2', 'DSS', 'baySeq', 'MDSeq','HM'), 
       lty=1, ncol=2)
abline(h=0.05, col='lightgrey')
plot(DE.results.50$mean.discoveries$edgeR.DE, 
     DE.results.50$mean.fdr$edgeR.DE, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("20 samples per group\n", "Differences in mean only"))
lines(DE.results.50$mean.discoveries$voom.DE, DE.results.50$mean.fdr$voom.DE, 
      type='l', col=col_vector[2])
lines(DE.results.50$mean.discoveries$DESeq2.DE, DE.results.50$mean.fdr$DESeq2.DE, 
      type='l', col=col_vector[3])
lines(DE.results.50$mean.discoveries$baySeq.DE, DE.results.50$mean.fdr$baySeq.DE, 
      type='l', col=col_vector[5])
lines(DE.results.50$mean.discoveries$MDSeq.DE, DE.results.50$mean.fdr$MDSeq.DE, 
      type='l', col=col_vector[6])
lines(DE.results.50$mean.discoveries$HM.DE, DE.results.50$mean.fdr$HM.DE, 
      type='l', col=col_vector[7])
legend("topright", bty='n', col=col_vector[c(1:3,5:7)], 
       legend=c('edgeR', 'voom', 'DESeq2', 'baySeq', 'MDSeq','HM'), 
       lty=1, ncol=2)
abline(h=0.05, col='lightgrey')
plot(DE.results.2$mean.discoveries$edgeR.DEDD, 
     DE.results.2$mean.fdr$edgeR.DEDD, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("2 samples per group\n", "Differences in mean and dispersion"))
lines(DE.results.2$mean.discoveries$voom.DEDD, DE.results.2$mean.fdr$voom.DEDD, 
      type='l', col=col_vector[2])
lines(DE.results.2$mean.discoveries$DESeq2.DEDD, DE.results.2$mean.fdr$DESeq2.DEDD, 
      type='l', col=col_vector[3])
lines(DE.results.2$mean.discoveries$DSS.DEDD, DE.results.2$mean.fdr$DSS.DEDD, 
      type='l', col=col_vector[4])
lines(DE.results.2$mean.discoveries$baySeq.DEDD, DE.results.2$mean.fdr$baySeq.DEDD, 
      type='l', col=col_vector[5])
lines(DE.results.2$mean.discoveries$MDSeq.DEDD, DE.results.2$mean.fdr$MDSeq.DEDD, 
      type='l', col=col_vector[6])
lines(DE.results.2$mean.discoveries$HM.DEDD, DE.results.2$mean.fdr$HM.DEDD, 
      type='l', col=col_vector[7])
legend("topright", bty='n', col=col_vector[1:7], 
       legend=c('edgeR', 'voom', 'DESeq2', 'DSS', 'baySeq', 'MDSeq','HM'), 
       lty=1, ncol=2)
abline(h=0.05, col='lightgrey')
plot(DE.results.5$mean.discoveries$edgeR.DEDD, 
     DE.results.5$mean.fdr$edgeR.DEDD, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("5 samples per group\n", "Differences in mean and dispersion"))
lines(DE.results.5$mean.discoveries$voom.DEDD, DE.results.5$mean.fdr$voom.DEDD, 
      type='l', col=col_vector[2])
lines(DE.results.5$mean.discoveries$DESeq2.DEDD, DE.results.5$mean.fdr$DESeq2.DEDD, 
      type='l', col=col_vector[3])
lines(DE.results.5$mean.discoveries$DSS.DEDD, DE.results.5$mean.fdr$DSS.DEDD, 
      type='l', col=col_vector[4])
lines(DE.results.5$mean.discoveries$baySeq.DEDD, DE.results.5$mean.fdr$baySeq.DEDD, 
      type='l', col=col_vector[5])
lines(DE.results.5$mean.discoveries$MDSeq.DEDD, DE.results.5$mean.fdr$MDSeq.DEDD, 
      type='l', col=col_vector[6])
lines(DE.results.5$mean.discoveries$HM.DEDD, DE.results.5$mean.fdr$HM.DEDD, 
      type='l', col=col_vector[7])
legend("topright", bty='n', col=col_vector[1:7], 
       legend=c('edgeR', 'voom', 'DESeq2', 'DSS', 'baySeq', 'MDSeq','HM'), 
       lty=1, ncol=2)
abline(h=0.05, col='lightgrey')
plot(DE.results.10$mean.discoveries$edgeR.DEDD, 
     DE.results.10$mean.fdr$edgeR.DEDD, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("10 samples per group\n", "Differences in mean and dispersion"))
lines(DE.results.10$mean.discoveries$voom.DEDD, DE.results.10$mean.fdr$voom.DEDD, 
      type='l', col=col_vector[2])
lines(DE.results.10$mean.discoveries$DESeq2.DEDD, DE.results.10$mean.fdr$DESeq2.DEDD, 
      type='l', col=col_vector[3])
lines(DE.results.10$mean.discoveries$DSS.DEDD, DE.results.10$mean.fdr$DSS.DEDD, 
      type='l', col=col_vector[4])
lines(DE.results.10$mean.discoveries$baySeq.DEDD, DE.results.10$mean.fdr$baySeq.DEDD, 
      type='l', col=col_vector[5])
lines(DE.results.10$mean.discoveries$MDSeq.DEDD, DE.results.10$mean.fdr$MDSeq.DEDD, 
      type='l', col=col_vector[6])
lines(DE.results.10$mean.discoveries$HM.DEDD, DE.results.10$mean.fdr$HM.DEDD, 
      type='l', col=col_vector[7])
legend("topright", bty='n', col=col_vector[1:7], 
       legend=c('edgeR', 'voom', 'DESeq2', 'DSS', 'baySeq', 'MDSeq','HM'), 
       lty=1, ncol=2)
abline(h=0.05, col='lightgrey')
plot(DE.results.20$mean.discoveries$edgeR.DEDD, 
     DE.results.20$mean.fdr$edgeR.DEDD, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("20 samples per group\n", "Differences in mean and dispersion"))
lines(DE.results.20$mean.discoveries$voom.DEDD, DE.results.20$mean.fdr$voom.DEDD, 
      type='l', col=col_vector[2])
lines(DE.results.20$mean.discoveries$DESeq2.DEDD, DE.results.20$mean.fdr$DESeq2.DEDD, 
      type='l', col=col_vector[3])
lines(DE.results.20$mean.discoveries$DSS.DEDD, DE.results.20$mean.fdr$DSS.DEDD, 
      type='l', col=col_vector[4])
lines(DE.results.20$mean.discoveries$baySeq.DEDD, DE.results.20$mean.fdr$baySeq.DEDD, 
      type='l', col=col_vector[5])
lines(DE.results.20$mean.discoveries$MDSeq.DEDD, DE.results.20$mean.fdr$MDSeq.DEDD, 
      type='l', col=col_vector[6])
lines(DE.results.20$mean.discoveries$HM.DEDD, DE.results.20$mean.fdr$HM.DEDD, 
      type='l', col=col_vector[7])
legend("topright", bty='n', col=col_vector[1:7], 
       legend=c('edgeR', 'voom', 'DESeq2', 'DSS', 'baySeq', 'MDSeq','HM'), 
       lty=1, ncol=2)
abline(h=0.05, col='lightgrey')
plot(DE.results.50$mean.discoveries$edgeR.DEDD, 
     DE.results.50$mean.fdr$edgeR.DEDD, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("50 samples per group\n", "Differences in mean and dispersion"))
lines(DE.results.50$mean.discoveries$voom.DEDD, DE.results.50$mean.fdr$voom.DEDD, 
      type='l', col=col_vector[2])
lines(DE.results.50$mean.discoveries$DESeq2.DEDD, DE.results.50$mean.fdr$DESeq2.DEDD, 
      type='l', col=col_vector[3])
lines(DE.results.50$mean.discoveries$baySeq.DEDD, DE.results.50$mean.fdr$baySeq.DEDD, 
      type='l', col=col_vector[5])
lines(DE.results.50$mean.discoveries$MDSeq.DEDD, DE.results.50$mean.fdr$MDSeq.DEDD, 
      type='l', col=col_vector[6])
lines(DE.results.50$mean.discoveries$HM.DEDD, DE.results.50$mean.fdr$HM.DEDD, 
      type='l', col=col_vector[7])
legend("topright", bty='n', col=col_vector[c(1:3,5:7)], 
       legend=c('edgeR', 'voom', 'DESeq2', 'baySeq', 'MDSeq','HM'), 
       lty=1, ncol=2)
abline(h=0.05, col='lightgrey')

par(mfrow=c(2,4), mar=c(2,2,3,1), mgp=c(3,0.7,0))
plot(DE.results.5$mean.discoveries$edgeR.DE, 
     DE.results.5$mean.fdr$edgeR.DE, 
     type='l', ylim=c(0,0.1), xlim=c(0,300), col=col_vector[1], 
     main=paste0("5 samples per group\n", "Differences in mean only"))
lines(DE.results.5$mean.discoveries$voom.DE, DE.results.5$mean.fdr$voom.DE, 
      type='l', col=col_vector[2])
lines(DE.results.5$mean.discoveries$DESeq2.DE, DE.results.5$mean.fdr$DESeq2.DE, 
      type='l', col=col_vector[3])
lines(DE.results.5$mean.discoveries$DSS.DE, DE.results.5$mean.fdr$DSS.DE, 
      type='l', col=col_vector[4])
lines(DE.results.5$mean.discoveries$baySeq.DE, DE.results.5$mean.fdr$baySeq.DE, 
      type='l', col=col_vector[5])
lines(DE.results.5$mean.discoveries$MDSeq.DE, DE.results.5$mean.fdr$MDSeq.DE, 
      type='l', col=col_vector[6])
lines(DE.results.5$mean.discoveries$HM.DE, DE.results.5$mean.fdr$HM.DE, 
      type='l', col=col_vector[7])
legend("topleft", bty='n', col=col_vector[1:7], 
       legend=c('edgeR', 'voom', 'DESeq2', 'DSS', 'baySeq', 'MDSeq','HM'), 
       lty=1, ncol=2)
abline(h=0.05, col='lightgrey')
plot(DE.results.10$mean.discoveries$edgeR.DE, DE.results.10$mean.fdr$edgeR.DE, 
     type='l', ylim=c(0,0.1), xlim=c(300,600), col=col_vector[1], 
     main=paste0("10 samples per group\n", "Differences in mean only"))
lines(DE.results.10$mean.discoveries$voom.DE, DE.results.10$mean.fdr$voom.DE, 
      type='l', col=col_vector[2])
lines(DE.results.10$mean.discoveries$DESeq2.DE, DE.results.10$mean.fdr$DESeq2.DE, 
      type='l', col=col_vector[3])
lines(DE.results.10$mean.discoveries$DSS.DE, DE.results.10$mean.fdr$DSS.DE, 
      type='l', col=col_vector[4])
lines(DE.results.10$mean.discoveries$baySeq.DE, DE.results.10$mean.fdr$baySeq.DE, 
      type='l', col=col_vector[5])
lines(DE.results.10$mean.discoveries$MDSeq.DE, DE.results.10$mean.fdr$MDSeq.DE, 
      type='l', col=col_vector[6])
lines(DE.results.10$mean.discoveries$HM.DE, DE.results.10$mean.fdr$HM.DE, 
      type='l', col=col_vector[7])
legend("topleft", bty='n', col=col_vector[1:7], 
       legend=c('edgeR', 'voom', 'DESeq2', 'DSS', 'baySeq', 'MDSeq','HM'), 
       lty=1, ncol=2)
abline(h=0.05, col='lightgrey')
plot(DE.results.20$mean.discoveries$edgeR.DE, 
     DE.results.20$mean.fdr$edgeR.DE, 
     type='l', ylim=c(0,0.1), xlim=c(500,750), col=col_vector[1], 
     main=paste0("20 samples per group\n", "Differences in mean only"))
lines(DE.results.20$mean.discoveries$voom.DE, DE.results.20$mean.fdr$voom.DE, 
      type='l', col=col_vector[2])
lines(DE.results.20$mean.discoveries$DESeq2.DE, DE.results.20$mean.fdr$DESeq2.DE, 
      type='l', col=col_vector[3])
lines(DE.results.20$mean.discoveries$DSS.DE, DE.results.20$mean.fdr$DSS.DE, 
      type='l', col=col_vector[4])
lines(DE.results.20$mean.discoveries$baySeq.DE, DE.results.20$mean.fdr$baySeq.DE, 
      type='l', col=col_vector[5])
lines(DE.results.20$mean.discoveries$MDSeq.DE, DE.results.20$mean.fdr$MDSeq.DE, 
      type='l', col=col_vector[6])
lines(DE.results.20$mean.discoveries$HM.DE, DE.results.20$mean.fdr$HM.DE, 
      type='l', col=col_vector[7])
legend("topleft", bty='n', col=col_vector[1:7], 
       legend=c('edgeR', 'voom', 'DESeq2', 'DSS', 'baySeq', 'MDSeq','HM'), 
       lty=1, ncol=2)
abline(h=0.05, col='lightgrey')
plot(DE.results.50$mean.discoveries$edgeR.DE, 
     DE.results.50$mean.fdr$edgeR.DE, 
     type='l', ylim=c(0,0.1), xlim=c(650,850), col=col_vector[1], 
     main=paste0("50 samples per group\n", "Differences in mean only"))
lines(DE.results.50$mean.discoveries$voom.DE, DE.results.50$mean.fdr$voom.DE, 
      type='l', col=col_vector[2])
lines(DE.results.50$mean.discoveries$DESeq2.DE, DE.results.50$mean.fdr$DESeq2.DE, 
      type='l', col=col_vector[3])
lines(DE.results.50$mean.discoveries$baySeq.DE, DE.results.50$mean.fdr$baySeq.DE, 
      type='l', col=col_vector[5])
lines(DE.results.50$mean.discoveries$MDSeq.DE, DE.results.50$mean.fdr$MDSeq.DE, 
      type='l', col=col_vector[6])
lines(DE.results.50$mean.discoveries$HM.DE, DE.results.50$mean.fdr$HM.DE, 
      type='l', col=col_vector[7])
legend("topleft", bty='n', col=col_vector[c(1:3,5:7)], 
       legend=c('edgeR', 'voom', 'DESeq2', 'baySeq', 'MDSeq','HM'), 
       lty=1, ncol=2)
abline(h=0.05, col='lightgrey')
plot(DE.results.5$mean.discoveries$edgeR.DEDD, 
     DE.results.5$mean.fdr$edgeR.DEDD, 
     type='l', ylim=c(0,0.1), xlim=c(0,800), col=col_vector[1], 
     main=paste0("5 samples per group\n", "Differences in mean and dispersion"))
lines(DE.results.5$mean.discoveries$voom.DEDD, DE.results.5$mean.fdr$voom.DEDD, 
      type='l', col=col_vector[2])
lines(DE.results.5$mean.discoveries$DESeq2.DEDD, DE.results.5$mean.fdr$DESeq2.DEDD, 
      type='l', col=col_vector[3])
lines(DE.results.5$mean.discoveries$DSS.DEDD, DE.results.5$mean.fdr$DSS.DEDD, 
      type='l', col=col_vector[4])
lines(DE.results.5$mean.discoveries$baySeq.DEDD, DE.results.5$mean.fdr$baySeq.DEDD, 
      type='l', col=col_vector[5])
lines(DE.results.5$mean.discoveries$MDSeq.DEDD, DE.results.5$mean.fdr$MDSeq.DEDD, 
      type='l', col=col_vector[6])
lines(DE.results.5$mean.discoveries$HM.DEDD, DE.results.5$mean.fdr$HM.DEDD, 
      type='l', col=col_vector[7])
legend("topleft", bty='n', col=col_vector[1:7], 
       legend=c('edgeR', 'voom', 'DESeq2', 'DSS', 'baySeq', 'MDSeq','HM'), 
       lty=1, ncol=2)
abline(h=0.05, col='lightgrey')
plot(DE.results.10$mean.discoveries$edgeR.DEDD, DE.results.10$mean.fdr$edgeR.DEDD, 
     type='l', ylim=c(0,0.1), xlim=c(650,1250), col=col_vector[1], 
     main=paste0("10 samples per group\n", "Differences in mean and dispersion"))
lines(DE.results.10$mean.discoveries$voom.DEDD, DE.results.10$mean.fdr$voom.DEDD, 
      type='l', col=col_vector[2])
lines(DE.results.10$mean.discoveries$DESeq2.DEDD, DE.results.10$mean.fdr$DESeq2.DEDD, 
      type='l', col=col_vector[3])
lines(DE.results.10$mean.discoveries$DSS.DEDD, DE.results.10$mean.fdr$DSS.DEDD, 
      type='l', col=col_vector[4])
lines(DE.results.10$mean.discoveries$baySeq.DEDD, DE.results.10$mean.fdr$baySeq.DEDD, 
      type='l', col=col_vector[5])
lines(DE.results.10$mean.discoveries$MDSeq.DEDD, DE.results.10$mean.fdr$MDSeq.DEDD, 
      type='l', col=col_vector[6])
lines(DE.results.10$mean.discoveries$HM.DEDD, DE.results.10$mean.fdr$HM.DEDD, 
      type='l', col=col_vector[7])
legend("topleft", bty='n', col=col_vector[1:7], 
       legend=c('edgeR', 'voom', 'DESeq2', 'DSS', 'baySeq', 'MDSeq','HM'), 
       lty=1, ncol=2)
abline(h=0.05, col='lightgrey')
plot(DE.results.20$mean.discoveries$edgeR.DEDD, 
     DE.results.20$mean.fdr$edgeR.DEDD, 
     type='l', ylim=c(0,0.1), xlim=c(800,1550), col=col_vector[1], 
     main=paste0("20 samples per group\n", "Differences in mean and dispersion"))
lines(DE.results.20$mean.discoveries$voom.DEDD, DE.results.20$mean.fdr$voom.DEDD, 
      type='l', col=col_vector[2])
lines(DE.results.20$mean.discoveries$DESeq2.DEDD, DE.results.20$mean.fdr$DESeq2.DEDD, 
      type='l', col=col_vector[3])
lines(DE.results.20$mean.discoveries$DSS.DEDD, DE.results.20$mean.fdr$DSS.DEDD, 
      type='l', col=col_vector[3])
lines(DE.results.20$mean.discoveries$baySeq.DEDD, DE.results.20$mean.fdr$baySeq.DEDD, 
      type='l', col=col_vector[5])
lines(DE.results.20$mean.discoveries$MDSeq.DEDD, DE.results.20$mean.fdr$MDSeq.DEDD, 
      type='l', col=col_vector[6])
lines(DE.results.20$mean.discoveries$HM.DEDD, DE.results.20$mean.fdr$HM.DEDD, 
      type='l', col=col_vector[7])
legend("topleft", bty='n', col=col_vector[1:7], 
       legend=c('edgeR', 'voom', 'DESeq2', 'DSS', 'baySeq', 'MDSeq','HM'), 
       lty=1, ncol=2)
abline(h=0.05, col='lightgrey')
plot(DE.results.50$mean.discoveries$edgeR.DEDD, 
     DE.results.50$mean.fdr$edgeR.DEDD, 
     type='l', ylim=c(0,0.1), xlim=c(800,1700), col=col_vector[1], 
     main=paste0("50 samples per group\n", "Differences in mean and dispersion"))
lines(DE.results.50$mean.discoveries$voom.DEDD, DE.results.50$mean.fdr$voom.DEDD, 
      type='l', col=col_vector[2])
lines(DE.results.50$mean.discoveries$DESeq2.DEDD, DE.results.50$mean.fdr$DESeq2.DEDD, 
      type='l', col=col_vector[3])
lines(DE.results.50$mean.discoveries$baySeq.DEDD, DE.results.50$mean.fdr$baySeq.DEDD, 
      type='l', col=col_vector[5])
lines(DE.results.50$mean.discoveries$MDSeq.DEDD, DE.results.50$mean.fdr$MDSeq.DEDD, 
      type='l', col=col_vector[6])
lines(DE.results.50$mean.discoveries$HM.DEDD, DE.results.50$mean.fdr$HM.DEDD, 
      type='l', col=col_vector[7])
legend("topleft", bty='n', col=col_vector[c(1:3,5:7)], 
       legend=c('edgeR', 'voom', 'DESeq2', 'baySeq', 'MDSeq','HM'), 
       lty=1, ncol=2)
abline(h=0.05, col='lightgrey')


###################################
#### Differential distribution ####
###################################
DEDD.results.2 <- list()
DEDD.results.5 <- list()
DEDD.results.10 <- list()
DEDD.results.20 <- list()
DEDD.results.50 <- list()

DEDD.results.2$auc <- cbind(DEDD.results.DESeq.DD2$auc, DEDD.results.DESeq.DE2$auc, 
                            DEDD.results.DESeq.DEDD2$auc, combined_DEDD.results.DEDD2$auc)
DEDD.results.2$auc <- DEDD.results.2$auc[, c(2,4,6,9,11,12)]
names(DEDD.results.2$auc) <- c('HMM.DD', 'HMM.DE', 'HMM.DEDD', 'HM-HM.DEDD', 'MDSeq-MDSeq.DEDD', 'edgeR-HM.DEDD')
DEDD.results.2$auc <- cbind(DEDD.results.2$auc, DD.results.2$auc, DE.results.2$auc)
DEDD.results.5$auc <- cbind(DEDD.results.DESeq.DD5$auc, DEDD.results.DESeq.DE5$auc, 
                            DEDD.results.DESeq.DEDD5$auc, combined_DEDD.results.DEDD5$auc)
DEDD.results.5$auc <- DEDD.results.5$auc[, c(2,4,6,9,11,12)]
names(DEDD.results.5$auc) <- c('HMM.DD', 'HMM.DE', 'HMM.DEDD', 'HM-HM.DEDD', 'MDSeq-MDSeq.DEDD', 'edgeR-HM.DEDD')
DEDD.results.5$auc <- cbind(DEDD.results.5$auc, DD.results.5$auc, DE.results.5$auc)
DEDD.results.10$auc <- cbind(DEDD.results.DESeq.DD10$auc, DEDD.results.DESeq.DE10$auc, 
                             DEDD.results.DESeq.DEDD10$auc, combined_DEDD.results.DEDD10$auc)
DEDD.results.10$auc <- DEDD.results.10$auc[, c(2,4,6,9,11,12)]
names(DEDD.results.10$auc) <- c('HMM.DD', 'HMM.DE', 'HMM.DEDD', 'HM-HM.DEDD', 'MDSeq-MDSeq.DEDD', 'edgeR-HM.DEDD')
DEDD.results.10$auc <- cbind(DEDD.results.10$auc, DD.results.10$auc, DE.results.10$auc)
DEDD.results.20$auc <- cbind(DEDD.results.DESeq.DD20$auc, DEDD.results.DESeq.DE20$auc, 
                             DEDD.results.DESeq.DEDD20$auc, combined_DEDD.results.DEDD20$auc)
DEDD.results.20$auc <- DEDD.results.20$auc[, c(2,4,6,9,11,12)]
names(DEDD.results.20$auc) <- c('HMM.DD', 'HMM.DE', 'HMM.DEDD', 'HM-HM.DEDD', 'MDSeq-MDSeq.DEDD', 'edgeR-HM.DEDD')
DEDD.results.20$auc <- cbind(DEDD.results.20$auc, DD.results.20$auc, DE.results.20$auc)
DEDD.results.50$auc <- cbind(DEDD.results.DESeq.DD50$auc, DEDD.results.DESeq.DE50$auc, 
                             DEDD.results.DESeq.DEDD50$auc, combined_DEDD.results.DEDD50$auc)
DEDD.results.50$auc <- DEDD.results.50$auc[, c(2,4,6,9,11,12)]
names(DEDD.results.50$auc) <- c('HMM.DD', 'HMM.DE', 'HMM.DEDD', 'HM-HM.DEDD', 'MDSeq-MDSeq.DEDD', 'edgeR-HM.DEDD')
DEDD.results.50$auc <- cbind(DEDD.results.50$auc, DD.results.50$auc, DE.results.50$auc)

DEDD.results.2$pauc <- cbind(DEDD.results.DESeq.DD2$pauc, DEDD.results.DESeq.DE2$pauc, 
                            DEDD.results.DESeq.DEDD2$pauc, combined_DEDD.results.DEDD2$pauc)
DEDD.results.2$pauc <- DEDD.results.2$pauc[, c(2,4,6,9,11,12)]
names(DEDD.results.2$pauc) <- c('HMM.DD', 'HMM.DE', 'HMM.DEDD', 'HM-HM.DEDD', 'MDSeq-MDSeq.DEDD', 'edgeR-HM.DEDD')
DEDD.results.2$pauc <- cbind(DEDD.results.2$pauc, DD.results.2$pauc, DE.results.2$pauc)
DEDD.results.5$pauc <- cbind(DEDD.results.DESeq.DD5$pauc, DEDD.results.DESeq.DE5$pauc, 
                            DEDD.results.DESeq.DEDD5$pauc, combined_DEDD.results.DEDD5$pauc)
DEDD.results.5$pauc <- DEDD.results.5$pauc[, c(2,4,6,9,11,12)]
names(DEDD.results.5$pauc) <- c('HMM.DD', 'HMM.DE', 'HMM.DEDD', 'HM-HM.DEDD', 'MDSeq-MDSeq.DEDD', 'edgeR-HM.DEDD')
DEDD.results.5$pauc <- cbind(DEDD.results.5$pauc, DD.results.5$pauc, DE.results.5$pauc)
DEDD.results.10$pauc <- cbind(DEDD.results.DESeq.DD10$pauc, DEDD.results.DESeq.DE10$pauc, 
                             DEDD.results.DESeq.DEDD10$pauc, combined_DEDD.results.DEDD10$pauc)
DEDD.results.10$pauc <- DEDD.results.10$pauc[, c(2,4,6,9,11,12)]
names(DEDD.results.10$pauc) <- c('HMM.DD', 'HMM.DE', 'HMM.DEDD', 'HM-HM.DEDD', 'MDSeq-MDSeq.DEDD', 'edgeR-HM.DEDD')
DEDD.results.10$pauc <- cbind(DEDD.results.10$pauc, DD.results.10$pauc, DE.results.10$pauc)
DEDD.results.20$pauc <- cbind(DEDD.results.DESeq.DD20$pauc, DEDD.results.DESeq.DE20$pauc, 
                             DEDD.results.DESeq.DEDD20$pauc, combined_DEDD.results.DEDD20$pauc)
DEDD.results.20$pauc <- DEDD.results.20$pauc[, c(2,4,6,9,11,12)]
names(DEDD.results.20$pauc) <- c('HMM.DD', 'HMM.DE', 'HMM.DEDD', 'HM-HM.DEDD', 'MDSeq-MDSeq.DEDD', 'edgeR-HM.DEDD')
DEDD.results.20$pauc <- cbind(DEDD.results.20$pauc, DD.results.20$pauc, DE.results.20$pauc)
DEDD.results.50$pauc <- cbind(DEDD.results.DESeq.DD50$pauc, DEDD.results.DESeq.DE50$pauc, 
                              DEDD.results.DESeq.DEDD50$pauc, combined_DEDD.results.DEDD50$pauc)
DEDD.results.50$pauc <- DEDD.results.50$pauc[, c(2,4,6,9,11,12)]
names(DEDD.results.50$pauc) <- c('HMM.DD', 'HMM.DE', 'HMM.DEDD', 'HM-HM.DEDD', 'MDSeq-MDSeq.DEDD', 'edgeR-HM.DEDD')
DEDD.results.50$pauc <- cbind(DEDD.results.50$pauc, DD.results.50$pauc, DE.results.50$pauc)

DEDD.results.2$fdr <- cbind(DEDD.results.DESeq.DD2$fdr, DEDD.results.DESeq.DE2$fdr, 
                            DEDD.results.DESeq.DEDD2$fdr, combined_DEDD.results.DEDD2$fdr)
DEDD.results.2$fdr <- DEDD.results.2$fdr[, c(5,11,17,21,23,24)]
names(DEDD.results.2$fdr) <- c('HMM.DD', 'HMM.DE', 'HMM.DEDD', 'HM-HM.DEDD', 'MDSeq-MDSeq.DEDD', 'edgeR-HM.DEDD')
DEDD.results.2$fdr <- cbind(DEDD.results.2$fdr, DD.results.2$fdr, DE.results.2$fdr)
DEDD.results.5$fdr <- cbind(DEDD.results.DESeq.DD5$fdr, DEDD.results.DESeq.DE5$fdr, 
                            DEDD.results.DESeq.DEDD5$fdr, combined_DEDD.results.DEDD5$fdr)
DEDD.results.5$fdr <- DEDD.results.5$fdr[, c(5,11,17,21,23,24)]
names(DEDD.results.5$fdr) <- c('HMM.DD', 'HMM.DE', 'HMM.DEDD', 'HM-HM.DEDD', 'MDSeq-MDSeq.DEDD', 'edgeR-HM.DEDD')
DEDD.results.5$fdr <- cbind(DEDD.results.5$fdr, DD.results.5$fdr, DE.results.5$fdr)
DEDD.results.10$fdr <- cbind(DEDD.results.DESeq.DD10$fdr, DEDD.results.DESeq.DE10$fdr, 
                             DEDD.results.DESeq.DEDD10$fdr, combined_DEDD.results.DEDD10$fdr)
DEDD.results.10$fdr <- DEDD.results.10$fdr[, c(5,11,17,21,23,24)]
names(DEDD.results.10$fdr) <- c('HMM.DD', 'HMM.DE', 'HMM.DEDD', 'HM-HM.DEDD', 'MDSeq-MDSeq.DEDD', 'edgeR-HM.DEDD')
DEDD.results.10$fdr <- cbind(DEDD.results.10$fdr, DD.results.10$fdr, DE.results.10$fdr)
DEDD.results.20$fdr <- cbind(DEDD.results.DESeq.DD20$fdr, DEDD.results.DESeq.DE20$fdr, 
                             DEDD.results.DESeq.DEDD20$fdr, combined_DEDD.results.DEDD20$fdr)
DEDD.results.20$fdr <- DEDD.results.20$fdr[, c(5,11,17,21,23,24)]
names(DEDD.results.20$fdr) <- c('HMM.DD', 'HMM.DE', 'HMM.DEDD', 'HM-HM.DEDD', 'MDSeq-MDSeq.DEDD', 'edgeR-HM.DEDD')
DEDD.results.20$fdr <- cbind(DEDD.results.20$fdr, DD.results.20$fdr, DE.results.20$fdr)
DEDD.results.50$fdr <- cbind(DEDD.results.DESeq.DD50$fdr, DEDD.results.DESeq.DE50$fdr, 
                             DEDD.results.DESeq.DEDD50$fdr, combined_DEDD.results.DEDD50$fdr)
DEDD.results.50$fdr <- DEDD.results.50$fdr[, c(5,11,17,21,23,24)]
names(DEDD.results.50$fdr) <- c('HMM.DD', 'HMM.DE', 'HMM.DEDD', 'HM-HM.DEDD', 'MDSeq-MDSeq.DEDD', 'edgeR-HM.DEDD')
DEDD.results.50$fdr <- cbind(DEDD.results.50$fdr, DD.results.50$fdr, DE.results.50$fdr)

DEDD.results.2$tpr <- cbind(DEDD.results.DESeq.DD2$tpr, DEDD.results.DESeq.DE2$tpr, 
                            DEDD.results.DESeq.DEDD2$tpr, combined_DEDD.results.DEDD2$tpr)
DEDD.results.2$tpr <- DEDD.results.2$tpr[, c(5,11,17,21,23,24)]
names(DEDD.results.2$tpr) <- c('HMM.DD', 'HMM.DE', 'HMM.DEDD', 'HM-HM.DEDD', 'MDSeq-MDSeq.DEDD', 'edgeR-HM.DEDD')
DEDD.results.2$tpr <- cbind(DEDD.results.2$tpr, DD.results.2$tpr, DE.results.2$tpr)
DEDD.results.5$tpr <- cbind(DEDD.results.DESeq.DD5$tpr, DEDD.results.DESeq.DE5$tpr, 
                            DEDD.results.DESeq.DEDD5$tpr, combined_DEDD.results.DEDD5$tpr)
DEDD.results.5$tpr <- DEDD.results.5$tpr[, c(5,11,17,21,23,24)]
names(DEDD.results.5$tpr) <- c('HMM.DD', 'HMM.DE', 'HMM.DEDD', 'HM-HM.DEDD', 'MDSeq-MDSeq.DEDD', 'edgeR-HM.DEDD')
DEDD.results.5$tpr <- cbind(DEDD.results.5$tpr, DD.results.5$tpr, DE.results.5$tpr)
DEDD.results.10$tpr <- cbind(DEDD.results.DESeq.DD10$tpr, DEDD.results.DESeq.DE10$tpr, 
                             DEDD.results.DESeq.DEDD10$tpr, combined_DEDD.results.DEDD10$tpr)
DEDD.results.10$tpr <- DEDD.results.10$tpr[, c(5,11,17,21,23,24)]
names(DEDD.results.10$tpr) <- c('HMM.DD', 'HMM.DE', 'HMM.DEDD', 'HM-HM.DEDD', 'MDSeq-MDSeq.DEDD', 'edgeR-HM.DEDD')
DEDD.results.10$tpr <- cbind(DEDD.results.10$tpr, DD.results.10$tpr, DE.results.10$tpr)
DEDD.results.20$tpr <- cbind(DEDD.results.DESeq.DD20$tpr, DEDD.results.DESeq.DE20$tpr, 
                             DEDD.results.DESeq.DEDD20$tpr, combined_DEDD.results.DEDD20$tpr)
DEDD.results.20$tpr <- DEDD.results.20$tpr[, c(5,11,17,21,23,24)]
names(DEDD.results.20$tpr) <- c('HMM.DD', 'HMM.DE', 'HMM.DEDD', 'HM-HM.DEDD', 'MDSeq-MDSeq.DEDD', 'edgeR-HM.DEDD')
DEDD.results.20$tpr <- cbind(DEDD.results.20$tpr, DD.results.20$tpr, DE.results.20$tpr)
DEDD.results.50$tpr <- cbind(DEDD.results.DESeq.DD50$tpr, DEDD.results.DESeq.DE50$tpr, 
                             DEDD.results.DESeq.DEDD50$tpr, combined_DEDD.results.DEDD50$tpr)
DEDD.results.50$tpr <- DEDD.results.50$tpr[, c(5,11,17,21,23,24)]
names(DEDD.results.50$tpr) <- c('HMM.DD', 'HMM.DE', 'HMM.DEDD', 'HM-HM.DEDD', 'MDSeq-MDSeq.DEDD', 'edgeR-HM.DEDD')
DEDD.results.50$tpr <- cbind(DEDD.results.50$tpr, DD.results.50$tpr, DE.results.50$tpr)

DEDD.results.2$mean.discoveries <- cbind(DEDD.results.DESeq.DD2$mean.discoveries, DEDD.results.DESeq.DE2$mean.discoveries, 
                            DEDD.results.DESeq.DEDD2$mean.discoveries, combined_DEDD.results.DEDD2$mean.discoveries)
DEDD.results.2$mean.discoveries <- DEDD.results.2$mean.discoveries[, c(2,4,6,9,11,12)]
names(DEDD.results.2$mean.discoveries) <- c('HMM.DD', 'HMM.DE', 'HMM.DEDD', 'HM-HM.DEDD', 'MDSeq-MDSeq.DEDD', 'edgeR-HM.DEDD')
DEDD.results.2$mean.discoveries <- cbind(DEDD.results.2$mean.discoveries, DD.results.2$mean.discoveries, DE.results.2$mean.discoveries)
DEDD.results.5$mean.discoveries <- cbind(DEDD.results.DESeq.DD5$mean.discoveries, DEDD.results.DESeq.DE5$mean.discoveries, 
                            DEDD.results.DESeq.DEDD5$mean.discoveries, combined_DEDD.results.DEDD5$mean.discoveries)
DEDD.results.5$mean.discoveries <- DEDD.results.5$mean.discoveries[, c(2,4,6,9,11,12)]
names(DEDD.results.5$mean.discoveries) <- c('HMM.DD', 'HMM.DE', 'HMM.DEDD', 'HM-HM.DEDD', 'MDSeq-MDSeq.DEDD', 'edgeR-HM.DEDD')
DEDD.results.5$mean.discoveries <- cbind(DEDD.results.5$mean.discoveries, DD.results.5$mean.discoveries, DE.results.5$mean.discoveries)
DEDD.results.10$mean.discoveries <- cbind(DEDD.results.DESeq.DD10$mean.discoveries, DEDD.results.DESeq.DE10$mean.discoveries, 
                             DEDD.results.DESeq.DEDD10$mean.discoveries, combined_DEDD.results.DEDD10$mean.discoveries)
DEDD.results.10$mean.discoveries <- DEDD.results.10$mean.discoveries[, c(2,4,6,9,11,12)]
names(DEDD.results.10$mean.discoveries) <- c('HMM.DD', 'HMM.DE', 'HMM.DEDD', 'HM-HM.DEDD', 'MDSeq-MDSeq.DEDD', 'edgeR-HM.DEDD')
DEDD.results.10$mean.discoveries <- cbind(DEDD.results.10$mean.discoveries, DD.results.10$mean.discoveries, DE.results.10$mean.discoveries)
DEDD.results.20$mean.discoveries <- cbind(DEDD.results.DESeq.DD20$mean.discoveries, DEDD.results.DESeq.DE20$mean.discoveries, 
                             DEDD.results.DESeq.DEDD20$mean.discoveries, combined_DEDD.results.DEDD20$mean.discoveries)
DEDD.results.20$mean.discoveries <- DEDD.results.20$mean.discoveries[, c(2,4,6,9,11,12)]
names(DEDD.results.20$mean.discoveries) <- c('HMM.DD', 'HMM.DE', 'HMM.DEDD', 'HM-HM.DEDD', 'MDSeq-MDSeq.DEDD', 'edgeR-HM.DEDD')
DEDD.results.20$mean.discoveries <- cbind(DEDD.results.20$mean.discoveries, DD.results.20$mean.discoveries, DE.results.20$mean.discoveries)
DEDD.results.50$mean.discoveries <- cbind(DEDD.results.DESeq.DD50$mean.discoveries, DEDD.results.DESeq.DE50$mean.discoveries, 
                                          DEDD.results.DESeq.DEDD50$mean.discoveries, combined_DEDD.results.DEDD50$mean.discoveries)
DEDD.results.50$mean.discoveries <- DEDD.results.50$mean.discoveries[, c(2,4,6,9,11,12)]
names(DEDD.results.50$mean.discoveries) <- c('HMM.DD', 'HMM.DE', 'HMM.DEDD', 'HM-HM.DEDD', 'MDSeq-MDSeq.DEDD', 'edgeR-HM.DEDD')
DEDD.results.50$mean.discoveries <- cbind(DEDD.results.50$mean.discoveries, DD.results.50$mean.discoveries, DE.results.50$mean.discoveries)

DEDD.results.2$mean.fdr <- cbind(DEDD.results.DESeq.DD2$mean.fdr, DEDD.results.DESeq.DE2$mean.fdr, 
                            DEDD.results.DESeq.DEDD2$mean.fdr, combined_DEDD.results.DEDD2$mean.fdr)
DEDD.results.2$mean.fdr <- DEDD.results.2$mean.fdr[, c(2,4,6,9,11,12)]
names(DEDD.results.2$mean.fdr) <- c('HMM.DD', 'HMM.DE', 'HMM.DEDD', 'HM-HM.DEDD', 'MDSeq-MDSeq.DEDD', 'edgeR-HM.DEDD')
DEDD.results.2$mean.fdr <- cbind(DEDD.results.2$mean.fdr, DD.results.2$mean.fdr, DE.results.2$mean.fdr)
DEDD.results.5$mean.fdr <- cbind(DEDD.results.DESeq.DD5$mean.fdr, DEDD.results.DESeq.DE5$mean.fdr, 
                            DEDD.results.DESeq.DEDD5$mean.fdr, combined_DEDD.results.DEDD5$mean.fdr)
DEDD.results.5$mean.fdr <- DEDD.results.5$mean.fdr[, c(2,4,6,9,11,12)]
names(DEDD.results.5$mean.fdr) <- c('HMM.DD', 'HMM.DE', 'HMM.DEDD', 'HM-HM.DEDD', 'MDSeq-MDSeq.DEDD', 'edgeR-HM.DEDD')
DEDD.results.5$mean.fdr <- cbind(DEDD.results.5$mean.fdr, DD.results.5$mean.fdr, DE.results.5$mean.fdr)
DEDD.results.10$mean.fdr <- cbind(DEDD.results.DESeq.DD10$mean.fdr, DEDD.results.DESeq.DE10$mean.fdr, 
                             DEDD.results.DESeq.DEDD10$mean.fdr, combined_DEDD.results.DEDD10$mean.fdr)
DEDD.results.10$mean.fdr <- DEDD.results.10$mean.fdr[, c(2,4,6,9,11,12)]
names(DEDD.results.10$mean.fdr) <- c('HMM.DD', 'HMM.DE', 'HMM.DEDD', 'HM-HM.DEDD', 'MDSeq-MDSeq.DEDD', 'edgeR-HM.DEDD')
DEDD.results.10$mean.fdr <- cbind(DEDD.results.10$mean.fdr, DD.results.10$mean.fdr, DE.results.10$mean.fdr)
DEDD.results.20$mean.fdr <- cbind(DEDD.results.DESeq.DD20$mean.fdr, DEDD.results.DESeq.DE20$mean.fdr, 
                             DEDD.results.DESeq.DEDD20$mean.fdr, combined_DEDD.results.DEDD20$mean.fdr)
DEDD.results.20$mean.fdr <- DEDD.results.20$mean.fdr[, c(2,4,6,9,11,12)]
names(DEDD.results.20$mean.fdr) <- c('HMM.DD', 'HMM.DE', 'HMM.DEDD', 'HM-HM.DEDD', 'MDSeq-MDSeq.DEDD', 'edgeR-HM.DEDD')
DEDD.results.20$mean.fdr <- cbind(DEDD.results.20$mean.fdr, DD.results.20$mean.fdr, DE.results.20$mean.fdr)
DEDD.results.50$mean.fdr <- cbind(DEDD.results.DESeq.DD50$mean.fdr, DEDD.results.DESeq.DE50$mean.fdr, 
                                  DEDD.results.DESeq.DEDD50$mean.fdr, combined_DEDD.results.DEDD50$mean.fdr)
DEDD.results.50$mean.fdr <- DEDD.results.50$mean.fdr[, c(2,4,6,9,11,12)]
names(DEDD.results.50$mean.fdr) <- c('HMM.DD', 'HMM.DE', 'HMM.DEDD', 'HM-HM.DEDD', 'MDSeq-MDSeq.DEDD', 'edgeR-HM.DEDD')
DEDD.results.50$mean.fdr <- cbind(DEDD.results.50$mean.fdr, DD.results.50$mean.fdr, DE.results.50$mean.fdr)

##############################################################
### AUC, false discovery plots, FDR, TPR for HMM with DEDD ###
par(mfrow=c(2,2), mar=c(2,2,3,1), mgp=c(3,0.7,0))
boxplot(DEDD.results.2$auc[,3], DEDD.results.5$auc[,3], 
        DEDD.results.10$auc[,3], DEDD.results.20$auc[,3], 
        DEDD.results.50$auc[,3], 
        names=NA, col=col_vector[1:5], xaxt='n', ylim=c(0.6,0.95), 
        pch=20, main=paste0("AUC"))
legend("topleft", fill=col_vector[1:5], bty='n', 
       legend=c("2 samples per group", "5 samples per group", 
                "10 samples per group", "20 samples per group", 
                "50 samples per group"))
plot(DE.results.2$mean.discoveries[,3], DE.results.2$mean.fdr[,3], 
     type='l', ylim=c(0,0.8), xlim=c(0,1000), col=col_vector[1], 
     main=paste0("False discovery curves"))
lines(DE.results.5$mean.discoveries[,3], DE.results.5$mean.fdr[,3], 
      type='l', col=col_vector[2])
lines(DE.results.10$mean.discoveries[,3], DE.results.10$mean.fdr[,3], 
      type='l', col=col_vector[3])
lines(DE.results.20$mean.discoveries[,3], DE.results.20$mean.fdr[,3], 
      type='l', col=col_vector[4])
lines(DE.results.50$mean.discoveries[,3], DE.results.50$mean.fdr[,3], 
      type='l', col=col_vector[5])
legend("topleft", bty='n', col=col_vector[1:5], lty=1, 
       legend=c("2 samples per group", "5 samples per group", 
                "10 samples per group", "20 samples per group", 
                "50 samples per group"))
abline(h=0.05, col='lightgrey')
boxplot(DEDD.results.2$fdr[,3], DEDD.results.5$fdr[,3], 
        DEDD.results.10$fdr[,3], DEDD.results.20$fdr[,3], 
        DEDD.results.50$fdr[,3], 
        names=NA, col=col_vector[1:5], xaxt='n', ylim=c(0,0.75), 
        pch=20, main=paste0("FDR"))
legend("topright", fill=col_vector[1:5], bty='n', 
       legend=c("2 samples per group", "5 samples per group", 
                "10 samples per group", "20 samples per group", 
                "50 samples per group"))
abline(h=0.05, col='lightgrey')
boxplot(DEDD.results.2$tpr[,3], DEDD.results.5$tpr[,3], 
        DEDD.results.10$tpr[,3], DEDD.results.20$tpr[,3], 
        DEDD.results.50$tpr[,3], 
        names=NA, col=col_vector[1:5], xaxt='n', ylim=c(0,0.7), 
        pch=20, main=paste0("TPR"))
legend("topleft", fill=col_vector[1:5], bty='n', 
       legend=c("2 samples per group", "5 samples per group", 
                "10 samples per group", "20 samples per group", 
                "50 samples per group"))


###############################################
### AUC for HMM v combined methods for DEDD ###
par(mfrow=c(2,2), mar=c(1,2.5,3,1), mgp=c(3,0.7,0))
# boxplot(DEDD.results.2$auc[,3:6], 
#         names=NA, col=col_vector[1:5], xaxt='n', ylim=c(0.64,0.96), 
#         pch=20, main=paste0("AUC, 2 samples per group"))
# legend("topleft", fill=col_vector[1:5], bty='n', cex=1.4, 
#        legend=c("Mixture model", "HM/HM", "MDSeq/MDSeq", "edgeR/HM"))
boxplot(DEDD.results.5$auc[,3:6], 
        names=NA, col=col_vector[1:5], xaxt='n', ylim=c(0.73,0.96), 
        pch=20, cex.axis=2, 
        main=paste0("AUC, 5 samples per group"), cex.main=2)
# legend("topleft", fill=col_vector[1:5], bty='n', cex=1.4, 
#        legend=c("Mixture model", "HM/HM", "MDSeq/MDSeq", "edgeR/HM"))
# boxplot(DEDD.results.10$auc[,3:6], 
#         names=NA, col=col_vector[1:5], xaxt='n', ylim=c(0.73,0.96), 
#         pch=20, main=paste0("AUC, 10 samples per group"), cex.main=1.5)
# legend("bottomright", fill=col_vector[1:5], bty='n', cex=1.4, 
#        legend=c("Mixture model", "HM/HM", "MDSeq/MDSeq", "edgeR/HM"))
# boxplot(DEDD.results.20$auc[,3:6], 
#         names=NA, col=col_vector[1:5], xaxt='n', ylim=c(0.73,0.96), 
#         pch=20, main=paste0("AUC, 20 samples per group"), cex.main=1.5)
# legend("bottomright", fill=col_vector[1:5], bty='n', cex=1.4, 
#        legend=c("Mixture model", "HM/HM", "MDSeq/MDSeq", "edgeR/HM"))
boxplot(DEDD.results.50$auc[,3:6], 
        names=NA, col=col_vector[1:5], xaxt='n', ylim=c(0.73,0.96), 
        pch=20, cex.axis=2, 
        main=paste0("AUC, 50 samples per group"), cex.main=2)
# legend("bottomright", fill=col_vector[1:5], bty='n', cex=1.4, 
#        legend=c("Mixture model", "HM/HM", "MDSeq/MDSeq", "edgeR/HM"))


################################################
### pAUC for HMM v combined methods for DEDD ###
par(mfrow=c(1,4), mar=c(1,2,2,1), mgp=c(3,0.7,0))
# boxplot(DEDD.results.2$pauc[,3:6], 
#         names=NA, col=col_vector[1:4], xaxt='n', ylim=c(0.006,0.041), 
#         pch=20, main=paste0("Partial AUC, 2 samples per group"))
# legend("topleft", fill=col_vector[1:4], bty='n', cex=1.4, 
#        legend=c("Mixture model", "HM/HM", "MDSeq/MDSeq", "edgeR/HM"))
boxplot(DEDD.results.5$pauc[,3:6], 
        names=NA, col=col_vector[1:4], xaxt='n', ylim=c(0.006,0.041), 
        pch=20, main=paste0("Partial AUC, 5 samples per group"))
legend("topleft", fill=col_vector[1:4], bty='n', cex=1.4, 
       legend=c("Mixture model", "HM/HM", "MDSeq/MDSeq", "edgeR/HM"))
boxplot(DEDD.results.10$pauc[,3:6], 
        names=NA, col=col_vector[1:4], xaxt='n', ylim=c(0.006,0.041), 
        pch=20, main=paste0("Partial AUC, 10 samples per group"))
legend("bottomright", fill=col_vector[1:4], bty='n', cex=1.4, 
       legend=c("Mixture model", "HM/HM", "MDSeq/MDSeq", "edgeR/HM"))
boxplot(DEDD.results.20$pauc[,3:6], 
        names=NA, col=col_vector[1:4], xaxt='n', ylim=c(0.006,0.041), 
        pch=20, main=paste0("Partial AUC, 20 samples per group"))
legend("bottomright", fill=col_vector[1:4], bty='n', cex=1.4, 
       legend=c("Mixture model", "HM/HM", "MDSeq/MDSeq", "edgeR/HM"))
boxplot(DEDD.results.50$pauc[,3:6], 
        names=NA, col=col_vector[1:4], xaxt='n', ylim=c(0.006,0.041), 
        pch=20, main=paste0("Partial AUC, 50 samples per group"))
legend("bottomright", fill=col_vector[1:4], bty='n', cex=1.4, 
       legend=c("Mixture model", "HM/HM", "MDSeq/MDSeq", "edgeR/HM"))


###############################################
### FDR for HMM v combined methods for DEDD ###
par(mfrow=c(1,4), mar=c(1,2,2,1), mgp=c(3,0.7,0))
# boxplot(DEDD.results.2$fdr[,3:6], 
#         names=NA, col=col_vector[1:4], xaxt='n', ylim=c(0,0.76), 
#         pch=20, main=paste0("FDR, 2 samples per group"), cex.main=1.5)
# legend("bottomright", fill=col_vector[1:4], bty='n', cex=1.4, 
#        legend=c("Mixture model", "HM/HM", "MDSeq/MDSeq", "edgeR/HM"))
# abline(h=0.05, col='lightgrey')
boxplot(DEDD.results.5$fdr[,3:6], 
        names=NA, col=col_vector[1:4], xaxt='n', ylim=c(0,0.5), 
        pch=20, cex.axis=2, 
        main=paste0("FDR, 5 samples per group"), cex.main=2)
# legend("topleft", fill=col_vector[1:4], bty='n', cex=1.4, 
#        legend=c("Mixture model", "HM/HM", "MDSeq/MDSeq", "edgeR/HM"))
abline(h=0.05, col='grey')
# boxplot(DEDD.results.10$fdr[,3:6], 
#         names=NA, col=col_vector[1:4], xaxt='n', ylim=c(0,0.76), 
#         pch=20, main=paste0("FDR, 10 samples per group"), cex.main=1.5)
# legend("topleft", fill=col_vector[1:4], bty='n', cex=1.4, 
#        legend=c("Mixture model", "HM/HM", "MDSeq/MDSeq", "edgeR/HM"))
# abline(h=0.05, col='lightgrey')
# boxplot(DEDD.results.20$fdr[,3:6], 
#         names=NA, col=col_vector[1:4], xaxt='n', ylim=c(0,0.76), 
#         pch=20, main=paste0("FDR, 20 samples per group"), cex.main=1.5)
# legend("topleft", fill=col_vector[1:4], bty='n', cex=1.4, 
#        legend=c("Mixture model", "HM/HM", "MDSeq/MDSeq", "edgeR/HM"))
# abline(h=0.05, col='lightgrey')
boxplot(DEDD.results.50$fdr[,3:6], 
        names=NA, col=col_vector[1:4], xaxt='n', ylim=c(0,0.5), 
        pch=20, cex.axis=2, 
        main=paste0("FDR, 50 samples per group"), cex.main=2)
# legend("topleft", fill=col_vector[1:4], bty='n', cex=1.4, 
#        legend=c("Mixture model", "HM/HM", "MDSeq/MDSeq", "edgeR/HM"))
abline(h=0.05, col='grey')


###############################################
### TPR for HMM v combined methods for DEDD ###
par(mfrow=c(1,4), mar=c(1,2,2,1), mgp=c(3,0.7,0))
# boxplot(DEDD.results.2$tpr[,3:6], 
#         names=NA, col=col_vector[1:4], xaxt='n', ylim=c(0,0.81), 
#         pch=20, main=paste0("TPR, 2 samples per group"), cex.main=1.5)
# legend("topleft", fill=col_vector[1:4], bty='n', cex=1.4, 
#        legend=c("Mixture model", "HM/HM", "MDSeq/MDSeq", "edgeR/HM"))
boxplot(DEDD.results.5$tpr[,3:6], 
        names=NA, col=col_vector[1:4], xaxt='n', ylim=c(0,0.81), 
        pch=20, main=paste0("TPR, 5 samples per group"), cex.main=1.5)
legend("topleft", fill=col_vector[1:4], bty='n', cex=1.4, 
       legend=c("Mixture model", "HM/HM", "MDSeq/MDSeq", "edgeR/HM"))
boxplot(DEDD.results.10$tpr[,3:6], 
        names=NA, col=col_vector[1:4], xaxt='n', ylim=c(0,0.81), 
        pch=20, main=paste0("TPR, 10 samples per group"), cex.main=1.5)
legend("bottomright", fill=col_vector[1:4], bty='n', cex=1.4, 
       legend=c("Mixture model", "HM/HM", "MDSeq/MDSeq", "edgeR/HM"))
boxplot(DEDD.results.20$tpr[,3:6], 
        names=NA, col=col_vector[1:4], xaxt='n', ylim=c(0,0.81), 
        pch=20, main=paste0("TPR, 20 samples per group"), cex.main=1.5)
legend("bottomright", fill=col_vector[1:4], bty='n', cex=1.4, 
       legend=c("Mixture model", "HM/HM", "MDSeq/MDSeq", "edgeR/HM"))
boxplot(DEDD.results.50$tpr[,3:6], 
        names=NA, col=col_vector[1:4], xaxt='n', ylim=c(0,0.81), 
        pch=20, main=paste0("TPR, 50 samples per group"), cex.main=1.5)
legend("bottomright", fill=col_vector[1:4], bty='n', cex=1.4, 
       legend=c("Mixture model", "HM/HM", "MDSeq/MDSeq", "edgeR/HM"))


#####################################################
### False discovery plots for HMM v combined methods for DEDD ###
par(mfrow=c(1,4), mar=c(2,2,2,1), mgp=c(3,0.7,0))
plot(DEDD.results.5$mean.discoveries$HMM.DEDD, 
     DEDD.results.5$mean.fdr$HMM.DEDD, 
     type='l', ylim=c(0,0.1), xlim=c(0,800), col=col_vector[1], 
     main=paste0("5 samples per group"))
lines(DEDD.results.5$mean.discoveries$'HM-HM.DEDD', 
      DEDD.results.5$mean.fdr$'HM-HM.DEDD', 
      type='l', col=col_vector[2])
lines(DEDD.results.5$mean.discoveries$'MDSeq-MDSeq.DEDD', 
      DEDD.results.5$mean.fdr$'MDSeq-MDSeq.DEDD', 
      type='l', col=col_vector[3])
lines(DEDD.results.5$mean.discoveries$'edgeR-HM.DEDD', 
      DEDD.results.5$mean.fdr$'edgeR-HM.DEDD', 
      type='l', col=col_vector[4])
legend("topleft", bty='n', col=col_vector[1:4], lwd=1, cex=1.2, 
       legend=c("Mixture model", "HM/HM", "MDSeq/MDSeq", "edgeR/HM"))
abline(h=0.05, col='lightgrey')
plot(DEDD.results.10$mean.discoveries$HMM.DEDD, 
     DEDD.results.10$mean.fdr$HMM.DEDD, 
     type='l', ylim=c(0,0.1), xlim=c(400,1300), col=col_vector[1], 
     main=paste0("10 samples per group"))
lines(DEDD.results.10$mean.discoveries$'HM-HM.DEDD', 
      DEDD.results.10$mean.fdr$'HM-HM.DEDD', 
      type='l', col=col_vector[2])
lines(DEDD.results.10$mean.discoveries$'MDSeq-MDSeq.DEDD', 
      DEDD.results.10$mean.fdr$'MDSeq-MDSeq.DEDD', 
      type='l', col=col_vector[3])
lines(DEDD.results.10$mean.discoveries$'edgeR-HM.DEDD', 
      DEDD.results.10$mean.fdr$'edgeR-HM.DEDD', 
      type='l', col=col_vector[4])
legend("topleft", bty='n', col=col_vector[1:4], lwd=1,cex=1.2,  
       legend=c("Mixture model", "HM/HM", "MDSeq/MDSeq", "edgeR/HM"))
abline(h=0.05, col='lightgrey')
plot(DEDD.results.20$mean.discoveries$HMM.DEDD, 
     DEDD.results.20$mean.fdr$HMM.DEDD, 
     type='l', ylim=c(0,0.1), xlim=c(1000,1700), col=col_vector[1], 
     main=paste0("20 samples per group"))
lines(DEDD.results.20$mean.discoveries$'HM-HM.DEDD', 
      DEDD.results.20$mean.fdr$'HM-HM.DEDD', 
      type='l', col=col_vector[2])
lines(DEDD.results.20$mean.discoveries$'MDSeq-MDSeq.DEDD', 
      DEDD.results.20$mean.fdr$'MDSeq-MDSeq.DEDD', 
      type='l', col=col_vector[3])
lines(DEDD.results.20$mean.discoveries$'edgeR-HM.DEDD', 
      DEDD.results.20$mean.fdr$'edgeR-HM.DEDD', 
      type='l', col=col_vector[4])
legend("topleft", bty='n', col=col_vector[1:4], lwd=1, cex=1.2, 
       legend=c("Mixture model", "HM/HM", "MDSeq/MDSeq", "edgeR/HM"))
abline(h=0.05, col='lightgrey')
plot(DEDD.results.50$mean.discoveries$HMM.DEDD, 
     DEDD.results.50$mean.fdr$HMM.DEDD, 
     type='l', ylim=c(0,0.1), xlim=c(1400,2000), col=col_vector[1], 
     main=paste0("50 samples per group"))
lines(DEDD.results.50$mean.discoveries$'HM-HM.DEDD', 
      DEDD.results.50$mean.fdr$'HM-HM.DEDD', 
      type='l', col=col_vector[2])
lines(DEDD.results.50$mean.discoveries$'MDSeq-MDSeq.DEDD', 
      DEDD.results.50$mean.fdr$'MDSeq-MDSeq.DEDD', 
      type='l', col=col_vector[3])
lines(DEDD.results.50$mean.discoveries$'edgeR-HM.DEDD', 
      DEDD.results.50$mean.fdr$'edgeR-HM.DEDD', 
      type='l', col=col_vector[4])
legend("topleft", bty='n', col=col_vector[1:4], lwd=1, cex=1.2, 
       legend=c("Mixture model", "HM/HM", "MDSeq/MDSeq", "edgeR/HM"))
abline(h=0.05, col='lightgrey')

