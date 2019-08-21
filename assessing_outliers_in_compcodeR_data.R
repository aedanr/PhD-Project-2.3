library(here)
library(compcodeR)
library(ggplot2)
library(tidyr)
library(ROCR)
library(caret)
library(edgeR)

###################
# DEDD2, 5, 10 ####
for (i in c('DE', 'DE.lfc1', 'DE.lfc2', 'DD', 'DD.lfc1', 'DD.lfc2', 'DEDD')) {
  for (j in c('DEDD2', 'DEDD5', 'DEDD10')) {
    assign(paste0(i,'.results.',j), readRDS(here('Results/DEDD compcodeR data results July-Aug 2019', 
                                                 paste0(i,'.results.',j,'.rds'))))
  }
}
rm('i','j')

# DE ####
DE.results.DEDD2.auc <- gather(DE.results.DEDD2$auc, 'method', 'auc')
DE.results.DEDD2.fdr <- gather(DE.results.DEDD2$fdr[-grep('raw',names(DE.results.DEDD2$fdr))], 'method', 'fdr')
DE.results.DEDD2.tpr <- gather(DE.results.DEDD2$tpr[-grep('raw',names(DE.results.DEDD2$tpr))], 'method', 'tpr')
DE.results.DEDD5.auc <- gather(DE.results.DEDD5$auc, 'method', 'auc')
DE.results.DEDD5.fdr <- gather(DE.results.DEDD5$fdr[-grep('raw',names(DE.results.DEDD5$fdr))], 'method', 'fdr')
DE.results.DEDD5.tpr <- gather(DE.results.DEDD5$tpr[-grep('raw',names(DE.results.DEDD5$tpr))], 'method', 'tpr')
DE.results.DEDD10.auc <- gather(DE.results.DEDD10$auc, 'method', 'auc')
DE.results.DEDD10.fdr <- gather(DE.results.DEDD10$fdr[-grep('raw',names(DE.results.DEDD10$fdr))], 'method', 'fdr')
DE.results.DEDD10.tpr <- gather(DE.results.DEDD10$tpr[-grep('raw',names(DE.results.DEDD10$tpr))], 'method', 'tpr')

outliers.DE.DEDD2.auc <- which(DE.results.DEDD2$auc$baySeq %in% sort(DE.results.DEDD2$auc$baySeq)[1:4])
colours.DE.DEDD2 <- rep('grey',50); colours.DE.DEDD2[outliers.DE.DEDD2.auc] = c('darkgreen','blue','yellow','darkorange')
outliers.DE.DEDD5.auc <- which(DE.results.DEDD5$auc$baySeq %in% sort(DE.results.DEDD5$auc$baySeq)[1:4])
colours.DE.DEDD5 <- rep('grey',50); colours.DE.DEDD5[outliers.DE.DEDD5.auc] = c('darkgreen','blue','yellow','darkorange')
outliers.DE.DEDD10.auc <- which(DE.results.DEDD10$auc$baySeq %in% sort(DE.results.DEDD10$auc$baySeq)[1:4])
colours.DE.DEDD10 <- rep('grey',50); colours.DE.DEDD10[outliers.DE.DEDD10.auc] = c('darkgreen','blue','yellow','darkorange')

ggplot(data = DE.results.DEDD2.auc, aes(x=method, y=auc)) + geom_point(col=rep(colours.DE.DEDD2,nrow(DE.results.DEDD2.auc)/50))
ggplot(data = DE.results.DEDD2.fdr, aes(x=method, y=fdr)) + geom_point(col=rep(colours.DE.DEDD2,nrow(DE.results.DEDD2.fdr)/50))
ggplot(data = DE.results.DEDD2.tpr, aes(x=method, y=tpr)) + geom_point(col=rep(colours.DE.DEDD2,nrow(DE.results.DEDD2.tpr)/50))
# Two obvious outliers for AUC for all methods except edgeR, voom and ShrinkBayes: 29, 32.
# Often but not always among highest FDR and TPR.
ggplot(data = DE.results.DEDD5.auc, aes(x=method, y=auc)) + geom_point(col=rep(colours.DE.DEDD5,nrow(DE.results.DEDD5.auc)/50))
ggplot(data = DE.results.DEDD5.fdr, aes(x=method, y=fdr)) + geom_point(col=rep(colours.DE.DEDD5,nrow(DE.results.DEDD5.fdr)/50))
ggplot(data = DE.results.DEDD5.tpr, aes(x=method, y=tpr)) + geom_point(col=rep(colours.DE.DEDD5,nrow(DE.results.DEDD5.tpr)/50))
# One extreme outlier for AUC and two less extreme, for all methods except edgeR and voom: 26, 23, 1.
# Almost always most extreme FDRs and TPRs as well.
ggplot(data = DE.results.DEDD10.auc, aes(x=method, y=auc)) + geom_point(col=rep(colours.DE.DEDD10,nrow(DE.results.DEDD10.auc)/50))
ggplot(data = DE.results.DEDD10.fdr, aes(x=method, y=fdr)) + geom_point(col=rep(colours.DE.DEDD10,nrow(DE.results.DEDD10.fdr)/50))
ggplot(data = DE.results.DEDD10.tpr, aes(x=method, y=tpr)) + geom_point(col=rep(colours.DE.DEDD10,nrow(DE.results.DEDD10.tpr)/50))
# One extreme, three moderate outliers for AUC for all methods except edgeR and voom: 45, 9, 22, 3.
# Also most extreme FDRs, but not TPRs.

# DE.lfc1 ####
DE.lfc1.results.DEDD2.auc <- gather(DE.lfc1.results.DEDD2$auc, 'method', 'auc')
DE.lfc1.results.DEDD2.fdr <- gather(DE.lfc1.results.DEDD2$fdr[-grep('raw',names(DE.lfc1.results.DEDD2$fdr))], 'method', 'fdr')
DE.lfc1.results.DEDD2.tpr <- gather(DE.lfc1.results.DEDD2$tpr[-grep('raw',names(DE.lfc1.results.DEDD2$tpr))], 'method', 'tpr')
DE.lfc1.results.DEDD5.auc <- gather(DE.lfc1.results.DEDD5$auc, 'method', 'auc')
DE.lfc1.results.DEDD5.fdr <- gather(DE.lfc1.results.DEDD5$fdr[-grep('raw',names(DE.lfc1.results.DEDD5$fdr))], 'method', 'fdr')
DE.lfc1.results.DEDD5.tpr <- gather(DE.lfc1.results.DEDD5$tpr[-grep('raw',names(DE.lfc1.results.DEDD5$tpr))], 'method', 'tpr')
DE.lfc1.results.DEDD10.auc <- gather(DE.lfc1.results.DEDD10$auc, 'method', 'auc')
DE.lfc1.results.DEDD10.fdr <- gather(DE.lfc1.results.DEDD10$fdr[-grep('raw',names(DE.lfc1.results.DEDD10$fdr))], 'method', 'fdr')
DE.lfc1.results.DEDD10.tpr <- gather(DE.lfc1.results.DEDD10$tpr[-grep('raw',names(DE.lfc1.results.DEDD10$tpr))], 'method', 'tpr')

outliers.DE.lfc1.DEDD2.auc <- which(DE.lfc1.results.DEDD2$auc$lfc1.DESeq %in% sort(DE.lfc1.results.DEDD2$auc$lfc1.DESeq)[1:4])
colours.DE.lfc1.DEDD2 <- rep('grey',50); colours.DE.lfc1.DEDD2[outliers.DE.lfc1.DEDD2.auc] = c('darkgreen','blue','yellow','darkorange')
outliers.DE.lfc1.DEDD5.auc <- which(DE.lfc1.results.DEDD5$auc$lfc1.DESeq %in% sort(DE.lfc1.results.DEDD5$auc$lfc1.DESeq)[1:4])
colours.DE.lfc1.DEDD5 <- rep('grey',50); colours.DE.lfc1.DEDD5[outliers.DE.lfc1.DEDD5.auc] = c('darkgreen','blue','yellow','darkorange')
outliers.DE.lfc1.DEDD10.auc <- which(DE.lfc1.results.DEDD10$auc$lfc1.DESeq %in% sort(DE.lfc1.results.DEDD10$auc$lfc1.DESeq)[1:4])
colours.DE.lfc1.DEDD10 <- rep('grey',50); colours.DE.lfc1.DEDD10[outliers.DE.lfc1.DEDD10.auc] = c('darkgreen','blue','yellow','darkorange')

ggplot(data = DE.lfc1.results.DEDD2.auc, aes(x=method, y=auc)) + geom_point(col=rep(colours.DE.lfc1.DEDD2,nrow(DE.lfc1.results.DEDD2.auc)/50))
ggplot(data = DE.lfc1.results.DEDD2.fdr, aes(x=method, y=fdr)) + geom_point(col=rep(colours.DE.lfc1.DEDD2,nrow(DE.lfc1.results.DEDD2.fdr)/50))
ggplot(data = DE.lfc1.results.DEDD2.tpr, aes(x=method, y=tpr)) + geom_point(col=rep(colours.DE.lfc1.DEDD2,nrow(DE.lfc1.results.DEDD2.tpr)/50))
# Outliers for AUC for DESeq2 match MDSeq and mostly match ShrinkBayes, but no real outliers for HMs - just lower AUCs 
# generally. AUCs for the DESeq2/MDSeq outliers are higher for HMs than for those methods. Outliers: 29, 32, 40, 3; 
# two most extreme are same as outliers for DE. No real outliers for FDR, TPR.
ggplot(data = DE.lfc1.results.DEDD5.auc, aes(x=method, y=auc)) + geom_point(col=rep(colours.DE.lfc1.DEDD5,nrow(DE.lfc1.results.DEDD5.auc)/50))
ggplot(data = DE.lfc1.results.DEDD5.fdr, aes(x=method, y=fdr)) + geom_point(col=rep(colours.DE.lfc1.DEDD5,nrow(DE.lfc1.results.DEDD5.fdr)/50))
ggplot(data = DE.lfc1.results.DEDD5.tpr, aes(x=method, y=tpr)) + geom_point(col=rep(colours.DE.lfc1.DEDD5,nrow(DE.lfc1.results.DEDD5.tpr)/50))
# Again, outliers for AUC for DESeq2 match MDSeq, and the datasets with lowest AUCs for those methods have highest 
# AUCs for HMs: 26, 23, 1, 42; same as for DE. No consistent outliers for FDR, TPR.
ggplot(data = DE.lfc1.results.DEDD10.auc, aes(x=method, y=auc)) + geom_point(col=rep(colours.DE.lfc1.DEDD10,nrow(DE.lfc1.results.DEDD10.auc)/50))
ggplot(data = DE.lfc1.results.DEDD10.fdr, aes(x=method, y=fdr)) + geom_point(col=rep(colours.DE.lfc1.DEDD10,nrow(DE.lfc1.results.DEDD10.fdr)/50))
ggplot(data = DE.lfc1.results.DEDD10.tpr, aes(x=method, y=tpr)) + geom_point(col=rep(colours.DE.lfc1.DEDD10,nrow(DE.lfc1.results.DEDD10.tpr)/50))
# One outlier for AUC, this time consistent between DESEq2, MDSeq and HMs: 45, which is the extreme outlier for DE. No 
# outliers for FDR, TPR.

# DE.lfc2 ####
DE.lfc2.results.DEDD2.auc <- gather(DE.lfc2.results.DEDD2$auc, 'method', 'auc')
DE.lfc2.results.DEDD2.fdr <- gather(DE.lfc2.results.DEDD2$fdr[-grep('raw',names(DE.lfc2.results.DEDD2$fdr))], 'method', 'fdr')
DE.lfc2.results.DEDD2.tpr <- gather(DE.lfc2.results.DEDD2$tpr[-grep('raw',names(DE.lfc2.results.DEDD2$tpr))], 'method', 'tpr')
DE.lfc2.results.DEDD5.auc <- gather(DE.lfc2.results.DEDD5$auc, 'method', 'auc')
DE.lfc2.results.DEDD5.fdr <- gather(DE.lfc2.results.DEDD5$fdr[-grep('raw',names(DE.lfc2.results.DEDD5$fdr))], 'method', 'fdr')
DE.lfc2.results.DEDD5.tpr <- gather(DE.lfc2.results.DEDD5$tpr[-grep('raw',names(DE.lfc2.results.DEDD5$tpr))], 'method', 'tpr')
DE.lfc2.results.DEDD10.auc <- gather(DE.lfc2.results.DEDD10$auc, 'method', 'auc')
DE.lfc2.results.DEDD10.fdr <- gather(DE.lfc2.results.DEDD10$fdr[-grep('raw',names(DE.lfc2.results.DEDD10$fdr))], 'method', 'fdr')
DE.lfc2.results.DEDD10.tpr <- gather(DE.lfc2.results.DEDD10$tpr[-grep('raw',names(DE.lfc2.results.DEDD10$tpr))], 'method', 'tpr')

outliers.DE.lfc2.DEDD2.auc <- which(DE.lfc2.results.DEDD2$auc$lfc2.DESeq %in% sort(DE.lfc2.results.DEDD2$auc$lfc2.DESeq)[1:4])
colours.DE.lfc2.DEDD2 <- rep('grey',50); colours.DE.lfc2.DEDD2[outliers.DE.lfc2.DEDD2.auc] = c('darkgreen','blue','yellow','darkorange')
outliers.DE.lfc2.DEDD5.auc <- which(DE.lfc2.results.DEDD5$auc$lfc2.DESeq %in% sort(DE.lfc2.results.DEDD5$auc$lfc2.DESeq)[1:4])
colours.DE.lfc2.DEDD5 <- rep('grey',50); colours.DE.lfc2.DEDD5[outliers.DE.lfc2.DEDD5.auc] = c('darkgreen','blue','yellow','darkorange')
outliers.DE.lfc2.DEDD10.auc <- which(DE.lfc2.results.DEDD10$auc$lfc2.DESeq %in% sort(DE.lfc2.results.DEDD10$auc$lfc2.DESeq)[1:4])
colours.DE.lfc2.DEDD10 <- rep('grey',50); colours.DE.lfc2.DEDD10[outliers.DE.lfc2.DEDD10.auc] = c('darkgreen','blue','yellow','darkorange')

ggplot(data = DE.lfc2.results.DEDD2.auc, aes(x=method, y=auc)) + geom_point(col=rep(colours.DE.lfc2.DEDD2,nrow(DE.lfc2.results.DEDD2.auc)/50))
ggplot(data = DE.lfc2.results.DEDD2.fdr, aes(x=method, y=fdr)) + geom_point(col=rep(colours.DE.lfc2.DEDD2,nrow(DE.lfc2.results.DEDD2.fdr)/50))
ggplot(data = DE.lfc2.results.DEDD2.tpr, aes(x=method, y=tpr)) + geom_point(col=rep(colours.DE.lfc2.DEDD2,nrow(DE.lfc2.results.DEDD2.tpr)/50))
ggplot(data = DE.lfc2.results.DEDD5.auc, aes(x=method, y=auc)) + geom_point(col=rep(colours.DE.lfc2.DEDD5,nrow(DE.lfc2.results.DEDD5.auc)/50))
ggplot(data = DE.lfc2.results.DEDD5.fdr, aes(x=method, y=fdr)) + geom_point(col=rep(colours.DE.lfc2.DEDD5,nrow(DE.lfc2.results.DEDD5.fdr)/50))
ggplot(data = DE.lfc2.results.DEDD5.tpr, aes(x=method, y=tpr)) + geom_point(col=rep(colours.DE.lfc2.DEDD5,nrow(DE.lfc2.results.DEDD5.tpr)/50))
ggplot(data = DE.lfc2.results.DEDD10.auc, aes(x=method, y=auc)) + geom_point(col=rep(colours.DE.lfc2.DEDD10,nrow(DE.lfc2.results.DEDD10.auc)/50))
ggplot(data = DE.lfc2.results.DEDD10.fdr, aes(x=method, y=fdr)) + geom_point(col=rep(colours.DE.lfc2.DEDD10,nrow(DE.lfc2.results.DEDD10.fdr)/50))
ggplot(data = DE.lfc2.results.DEDD10.tpr, aes(x=method, y=tpr)) + geom_point(col=rep(colours.DE.lfc2.DEDD10,nrow(DE.lfc2.results.DEDD10.tpr)/50))
# Basically same pattern as for lfc1 - outliers match between DESeq2 and MDSeq, but not HMs.

# DD ####
DD.results.DEDD2.auc <- gather(DD.results.DEDD2$auc, 'method', 'auc')
DD.results.DEDD2.fdr <- gather(DD.results.DEDD2$fdr[-grep('raw',names(DD.results.DEDD2$fdr))], 'method', 'fdr')
DD.results.DEDD2.tpr <- gather(DD.results.DEDD2$tpr[-grep('raw',names(DD.results.DEDD2$tpr))], 'method', 'tpr')
DD.results.DEDD5.auc <- gather(DD.results.DEDD5$auc, 'method', 'auc')
DD.results.DEDD5.fdr <- gather(DD.results.DEDD5$fdr[-grep('raw',names(DD.results.DEDD5$fdr))], 'method', 'fdr')
DD.results.DEDD5.tpr <- gather(DD.results.DEDD5$tpr[-grep('raw',names(DD.results.DEDD5$tpr))], 'method', 'tpr')
DD.results.DEDD10.auc <- gather(DD.results.DEDD10$auc, 'method', 'auc')
DD.results.DEDD10.fdr <- gather(DD.results.DEDD10$fdr[-grep('raw',names(DD.results.DEDD10$fdr))], 'method', 'fdr')
DD.results.DEDD10.tpr <- gather(DD.results.DEDD10$tpr[-grep('raw',names(DD.results.DEDD10$tpr))], 'method', 'tpr')

outliers.DD.DEDD2.auc <- which(DD.results.DEDD2$auc$disp.zi.MDSeq %in% sort(DD.results.DEDD2$auc$disp.zi.MDSeq)[1:4])
colours.DD.DEDD2 <- rep('grey',50); colours.DD.DEDD2[outliers.DD.DEDD2.auc] = c('darkgreen','blue','yellow','darkorange')
outliers.DD.DEDD5.auc <- which(DD.results.DEDD5$auc$disp.zi.MDSeq %in% sort(DD.results.DEDD5$auc$disp.zi.MDSeq)[1:4])
colours.DD.DEDD5 <- rep('grey',50); colours.DD.DEDD5[outliers.DD.DEDD5.auc] = c('darkgreen','blue','yellow','darkorange')
outliers.DD.DEDD10.auc <- which(DD.results.DEDD10$auc$disp.zi.MDSeq %in% sort(DD.results.DEDD10$auc$disp.zi.MDSeq)[1:4])
colours.DD.DEDD10 <- rep('grey',50); colours.DD.DEDD10[outliers.DD.DEDD10.auc] = c('darkgreen','blue','yellow','darkorange')

ggplot(data = DD.results.DEDD2.auc, aes(x=method, y=auc)) + geom_point(col=rep(colours.DD.DEDD2,nrow(DD.results.DEDD2.auc)/50))
ggplot(data = DD.results.DEDD2.fdr, aes(x=method, y=fdr)) + geom_point(col=rep(colours.DD.DEDD2,nrow(DD.results.DEDD2.fdr)/50))
ggplot(data = DD.results.DEDD2.tpr, aes(x=method, y=tpr)) + geom_point(col=rep(colours.DD.DEDD2,nrow(DD.results.DEDD2.tpr)/50))
# No consistent outliers for AUC. Some overlap between datasets with lowest few AUCs between MDSeq, expHM, lnHM, but 
# not really consistent, and no outliers for FDR, TPR.
ggplot(data = DD.results.DEDD5.auc, aes(x=method, y=auc)) + geom_point(col=rep(colours.DD.DEDD5,nrow(DD.results.DEDD5.auc)/50))
ggplot(data = DD.results.DEDD5.fdr, aes(x=method, y=fdr)) + geom_point(col=rep(colours.DD.DEDD5,nrow(DD.results.DEDD5.fdr)/50))
ggplot(data = DD.results.DEDD5.tpr, aes(x=method, y=tpr)) + geom_point(col=rep(colours.DD.DEDD5,nrow(DD.results.DEDD5.tpr)/50))
# Consistent overlap between datasets with lowest AUCs, but no consistent outliers. No outliers for FDR, TPR.
ggplot(data = DD.results.DEDD10.auc, aes(x=method, y=auc)) + geom_point(col=rep(colours.DD.DEDD10,nrow(DD.results.DEDD10.auc)/50))
ggplot(data = DD.results.DEDD10.fdr, aes(x=method, y=fdr)) + geom_point(col=rep(colours.DD.DEDD10,nrow(DD.results.DEDD10.fdr)/50))
ggplot(data = DD.results.DEDD10.tpr, aes(x=method, y=tpr)) + geom_point(col=rep(colours.DD.DEDD10,nrow(DD.results.DEDD10.tpr)/50))
# One clear outlier for AUC for HMs, which also has lowest AUC for MDSeq: 41. Looks like generally AUCs for HMs are 
# consistently slightly lower than for MDSeq. No outliers for FDR, TPR.

# DD.lfc1 ####
DD.lfc1.results.DEDD2.auc <- gather(DD.lfc1.results.DEDD2$auc, 'method', 'auc')
DD.lfc1.results.DEDD2.fdr <- gather(DD.lfc1.results.DEDD2$fdr[-grep('raw',names(DD.lfc1.results.DEDD2$fdr))], 'method', 'fdr')
DD.lfc1.results.DEDD2.tpr <- gather(DD.lfc1.results.DEDD2$tpr[-grep('raw',names(DD.lfc1.results.DEDD2$tpr))], 'method', 'tpr')
DD.lfc1.results.DEDD5.auc <- gather(DD.lfc1.results.DEDD5$auc, 'method', 'auc')
DD.lfc1.results.DEDD5.fdr <- gather(DD.lfc1.results.DEDD5$fdr[-grep('raw',names(DD.lfc1.results.DEDD5$fdr))], 'method', 'fdr')
DD.lfc1.results.DEDD5.tpr <- gather(DD.lfc1.results.DEDD5$tpr[-grep('raw',names(DD.lfc1.results.DEDD5$tpr))], 'method', 'tpr')
DD.lfc1.results.DEDD10.auc <- gather(DD.lfc1.results.DEDD10$auc, 'method', 'auc')
DD.lfc1.results.DEDD10.fdr <- gather(DD.lfc1.results.DEDD10$fdr[-grep('raw',names(DD.lfc1.results.DEDD10$fdr))], 'method', 'fdr')
DD.lfc1.results.DEDD10.tpr <- gather(DD.lfc1.results.DEDD10$tpr[-grep('raw',names(DD.lfc1.results.DEDD10$tpr))], 'method', 'tpr')

outliers.DD.lfc1.DEDD2.auc <- which(DD.lfc1.results.DEDD2$auc$disp.zi.lfc1.MDSeq %in% sort(DD.lfc1.results.DEDD2$auc$disp.zi.lfc1.MDSeq)[1:4])
colours.DD.lfc1.DEDD2 <- rep('grey',50); colours.DD.lfc1.DEDD2[outliers.DD.lfc1.DEDD2.auc] = c('darkgreen','blue','yellow','darkorange')
outliers.DD.lfc1.DEDD5.auc <- which(DD.lfc1.results.DEDD5$auc$disp.zi.lfc1.MDSeq %in% sort(DD.lfc1.results.DEDD5$auc$disp.zi.lfc1.MDSeq)[1:4])
colours.DD.lfc1.DEDD5 <- rep('grey',50); colours.DD.lfc1.DEDD5[outliers.DD.lfc1.DEDD5.auc] = c('darkgreen','blue','yellow','darkorange')
outliers.DD.lfc1.DEDD10.auc <- which(DD.lfc1.results.DEDD10$auc$disp.zi.lfc1.MDSeq %in% sort(DD.lfc1.results.DEDD10$auc$disp.zi.lfc1.MDSeq)[1:4])
colours.DD.lfc1.DEDD10 <- rep('grey',50); colours.DD.lfc1.DEDD10[outliers.DD.lfc1.DEDD10.auc] = c('darkgreen','blue','yellow','darkorange')

ggplot(data = DD.lfc1.results.DEDD2.auc, aes(x=method, y=auc)) + geom_point(col=rep(colours.DD.lfc1.DEDD2,nrow(DD.lfc1.results.DEDD2.auc)/50))
ggplot(data = DD.lfc1.results.DEDD2.fdr, aes(x=method, y=fdr)) + geom_point(col=rep(colours.DD.lfc1.DEDD2,nrow(DD.lfc1.results.DEDD2.fdr)/50))
ggplot(data = DD.lfc1.results.DEDD2.tpr, aes(x=method, y=tpr)) + geom_point(col=rep(colours.DD.lfc1.DEDD2,nrow(DD.lfc1.results.DEDD2.tpr)/50))
ggplot(data = DD.lfc1.results.DEDD5.auc, aes(x=method, y=auc)) + geom_point(col=rep(colours.DD.lfc1.DEDD5,nrow(DD.lfc1.results.DEDD5.auc)/50))
ggplot(data = DD.lfc1.results.DEDD5.fdr, aes(x=method, y=fdr)) + geom_point(col=rep(colours.DD.lfc1.DEDD5,nrow(DD.lfc1.results.DEDD5.fdr)/50))
ggplot(data = DD.lfc1.results.DEDD5.tpr, aes(x=method, y=tpr)) + geom_point(col=rep(colours.DD.lfc1.DEDD5,nrow(DD.lfc1.results.DEDD5.tpr)/50))
ggplot(data = DD.lfc1.results.DEDD10.auc, aes(x=method, y=auc)) + geom_point(col=rep(colours.DD.lfc1.DEDD10,nrow(DD.lfc1.results.DEDD10.auc)/50))
ggplot(data = DD.lfc1.results.DEDD10.fdr, aes(x=method, y=fdr)) + geom_point(col=rep(colours.DD.lfc1.DEDD10,nrow(DD.lfc1.results.DEDD10.fdr)/50))
ggplot(data = DD.lfc1.results.DEDD10.tpr, aes(x=method, y=tpr)) + geom_point(col=rep(colours.DD.lfc1.DEDD10,nrow(DD.lfc1.results.DEDD10.tpr)/50))
# No real outliers. Lowest AUCs generally overlap between MDSeq and HMs. Highest AUCs are similar for MDSeq and HMS, 
# but lower AUCs are lower for HMs.

# DD.lfc2 ####
DD.lfc2.results.DEDD2.auc <- gather(DD.lfc2.results.DEDD2$auc, 'method', 'auc')
DD.lfc2.results.DEDD2.fdr <- gather(DD.lfc2.results.DEDD2$fdr[-grep('raw',names(DD.lfc2.results.DEDD2$fdr))], 'method', 'fdr')
DD.lfc2.results.DEDD2.tpr <- gather(DD.lfc2.results.DEDD2$tpr[-grep('raw',names(DD.lfc2.results.DEDD2$tpr))], 'method', 'tpr')
DD.lfc2.results.DEDD5.auc <- gather(DD.lfc2.results.DEDD5$auc, 'method', 'auc')
DD.lfc2.results.DEDD5.fdr <- gather(DD.lfc2.results.DEDD5$fdr[-grep('raw',names(DD.lfc2.results.DEDD5$fdr))], 'method', 'fdr')
DD.lfc2.results.DEDD5.tpr <- gather(DD.lfc2.results.DEDD5$tpr[-grep('raw',names(DD.lfc2.results.DEDD5$tpr))], 'method', 'tpr')
DD.lfc2.results.DEDD10.auc <- gather(DD.lfc2.results.DEDD10$auc, 'method', 'auc')
DD.lfc2.results.DEDD10.fdr <- gather(DD.lfc2.results.DEDD10$fdr[-grep('raw',names(DD.lfc2.results.DEDD10$fdr))], 'method', 'fdr')
DD.lfc2.results.DEDD10.tpr <- gather(DD.lfc2.results.DEDD10$tpr[-grep('raw',names(DD.lfc2.results.DEDD10$tpr))], 'method', 'tpr')

outliers.DD.lfc2.DEDD2.auc <- which(DD.lfc2.results.DEDD2$auc$disp.lfc2.lnHM %in% sort(DD.lfc2.results.DEDD2$auc$disp.lfc2.lnHM)[1:4])
colours.DD.lfc2.DEDD2 <- rep('grey',50); colours.DD.lfc2.DEDD2[outliers.DD.lfc2.DEDD2.auc] = c('darkgreen','blue','yellow','darkorange')
outliers.DD.lfc2.DEDD5.auc <- which(DD.lfc2.results.DEDD5$auc$disp.lfc2.lnHM %in% sort(DD.lfc2.results.DEDD5$auc$disp.lfc2.lnHM)[1:4])
colours.DD.lfc2.DEDD5 <- rep('grey',50); colours.DD.lfc2.DEDD5[outliers.DD.lfc2.DEDD5.auc] = c('darkgreen','blue','yellow','darkorange')
outliers.DD.lfc2.DEDD10.auc <- which(DD.lfc2.results.DEDD10$auc$disp.lfc2.lnHM %in% sort(DD.lfc2.results.DEDD10$auc$disp.lfc2.lnHM)[1:4])
colours.DD.lfc2.DEDD10 <- rep('grey',50); colours.DD.lfc2.DEDD10[outliers.DD.lfc2.DEDD10.auc] = c('darkgreen','blue','yellow','darkorange')

ggplot(data = DD.lfc2.results.DEDD2.auc, aes(x=method, y=auc)) + geom_point(col=rep(colours.DD.lfc2.DEDD2,nrow(DD.lfc2.results.DEDD2.auc)/50))
ggplot(data = DD.lfc2.results.DEDD2.fdr, aes(x=method, y=fdr)) + geom_point(col=rep(colours.DD.lfc2.DEDD2,nrow(DD.lfc2.results.DEDD2.fdr)/50))
ggplot(data = DD.lfc2.results.DEDD2.tpr, aes(x=method, y=tpr)) + geom_point(col=rep(colours.DD.lfc2.DEDD2,nrow(DD.lfc2.results.DEDD2.tpr)/50))
ggplot(data = DD.lfc2.results.DEDD5.auc, aes(x=method, y=auc)) + geom_point(col=rep(colours.DD.lfc2.DEDD5,nrow(DD.lfc2.results.DEDD5.auc)/50))
ggplot(data = DD.lfc2.results.DEDD5.fdr, aes(x=method, y=fdr)) + geom_point(col=rep(colours.DD.lfc2.DEDD5,nrow(DD.lfc2.results.DEDD5.fdr)/50))
ggplot(data = DD.lfc2.results.DEDD5.tpr, aes(x=method, y=tpr)) + geom_point(col=rep(colours.DD.lfc2.DEDD5,nrow(DD.lfc2.results.DEDD5.tpr)/50))
ggplot(data = DD.lfc2.results.DEDD10.auc, aes(x=method, y=auc)) + geom_point(col=rep(colours.DD.lfc2.DEDD10,nrow(DD.lfc2.results.DEDD10.auc)/50))
ggplot(data = DD.lfc2.results.DEDD10.fdr, aes(x=method, y=fdr)) + geom_point(col=rep(colours.DD.lfc2.DEDD10,nrow(DD.lfc2.results.DEDD10.fdr)/50))
ggplot(data = DD.lfc2.results.DEDD10.tpr, aes(x=method, y=tpr)) + geom_point(col=rep(colours.DD.lfc2.DEDD10,nrow(DD.lfc2.results.DEDD10.tpr)/50))
# No real outliers. Not much overlap between lowest AUCs for HMs and MDSeq. Looks like most AUCs are similar 
# between HMs and MDSeq, but lowest for HMs are quite a lot lower than for same datasets for MDSeq, particularly 
# 18, 6, 21 for DEDD10.

# DEDD ####
DEDD.results.DEDD2.auc <- gather(DEDD.results.DEDD2$auc, 'method', 'auc')
DEDD.results.DEDD2.fdr <- gather(DEDD.results.DEDD2$fdr, 'method', 'fdr')
DEDD.results.DEDD2.tpr <- gather(DEDD.results.DEDD2$tpr, 'method', 'tpr')
DEDD.results.DEDD5.auc <- gather(DEDD.results.DEDD5$auc, 'method', 'auc')
DEDD.results.DEDD5.fdr <- gather(DEDD.results.DEDD5$fdr, 'method', 'fdr')
DEDD.results.DEDD5.tpr <- gather(DEDD.results.DEDD5$tpr, 'method', 'tpr')
DEDD.results.DEDD10.auc <- gather(DEDD.results.DEDD10$auc, 'method', 'auc')
DEDD.results.DEDD10.fdr <- gather(DEDD.results.DEDD10$fdr, 'method', 'fdr')
DEDD.results.DEDD10.tpr <- gather(DEDD.results.DEDD10$tpr, 'method', 'tpr')

outliers.DEDD.DEDD2.auc <- which(DEDD.results.DEDD2$auc$expHM %in% sort(DEDD.results.DEDD2$auc$expHM)[1:4])
colours.DEDD.DEDD2 <- rep('grey',50); colours.DEDD.DEDD2[outliers.DEDD.DEDD2.auc] = c('darkgreen','blue','yellow','darkorange')
outliers.DEDD.DEDD5.auc <- which(DEDD.results.DEDD5$auc$expHM %in% sort(DEDD.results.DEDD5$auc$expHM)[1:4])
colours.DEDD.DEDD5 <- rep('grey',50); colours.DEDD.DEDD5[outliers.DEDD.DEDD5.auc] = c('darkgreen','blue','yellow','darkorange')
outliers.DEDD.DEDD10.auc <- which(DEDD.results.DEDD10$auc$expHM %in% sort(DEDD.results.DEDD10$auc$expHM)[1:4])
colours.DEDD.DEDD10 <- rep('grey',50); colours.DEDD.DEDD10[outliers.DEDD.DEDD10.auc] = c('darkgreen','blue','yellow','darkorange')

ggplot(data = DEDD.results.DEDD2.auc, aes(x=method, y=auc)) + geom_point(col=rep(colours.DEDD.DEDD2,nrow(DEDD.results.DEDD2.auc)/50))
ggplot(data = DEDD.results.DEDD2.fdr, aes(x=method, y=fdr)) + geom_point(col=rep(colours.DEDD.DEDD2,nrow(DEDD.results.DEDD2.fdr)/50))
ggplot(data = DEDD.results.DEDD2.tpr, aes(x=method, y=tpr)) + geom_point(col=rep(colours.DEDD.DEDD2,nrow(DEDD.results.DEDD2.tpr)/50))
# No real outliers, datasets with lowest AUCs generally match.
ggplot(data = DEDD.results.DEDD5.auc, aes(x=method, y=auc)) + geom_point(col=rep(colours.DEDD.DEDD5,nrow(DEDD.results.DEDD5.auc)/50))
ggplot(data = DEDD.results.DEDD5.fdr, aes(x=method, y=fdr)) + geom_point(col=rep(colours.DEDD.DEDD5,nrow(DEDD.results.DEDD5.fdr)/50))
ggplot(data = DEDD.results.DEDD5.tpr, aes(x=method, y=tpr)) + geom_point(col=rep(colours.DEDD.DEDD5,nrow(DEDD.results.DEDD5.tpr)/50))
# One extreme outlier for AUC: 26, which was also an extreme outlier for DE, and an outlier for DE.lfc1 for DESeq2 
# and MDSeq but not for HMs. Same dataset has highest (too high) FDRs for thr, but not for .5, and highest TPR for 
# thr and close to highest for .5.
ggplot(data = DEDD.results.DEDD10.auc, aes(x=method, y=auc)) + geom_point(col=rep(colours.DEDD.DEDD10,nrow(DEDD.results.DEDD10.auc)/50))
ggplot(data = DEDD.results.DEDD10.fdr, aes(x=method, y=fdr)) + geom_point(col=rep(colours.DEDD.DEDD10,nrow(DEDD.results.DEDD10.fdr)/50))
ggplot(data = DEDD.results.DEDD10.tpr, aes(x=method, y=tpr)) + geom_point(col=rep(colours.DEDD.DEDD10,nrow(DEDD.results.DEDD10.tpr)/50))
# One extreme outlier for AUC: 45, which was also an extreme outlier for DE and an outlier for DE.lfc1. Same dataset 
# has highest (but still conservative) FDR for thr and 0.5, but is around average for TPR.


###################
# DE2, 5, 10 ####
for (i in c('DE', 'DE.lfc1', 'DE.lfc2', 'DEDD')) {
  for (j in c('DE2', 'DE5', 'DE10')) {
    assign(paste0(i,'.results.',j), readRDS(here('Results/DE compcodeR data results July-Aug 2019', 
                                                 paste0(i,'.results.',j,'.rds'))))
  }
}
rm('i','j')

# DE ####
DE.results.DE2.auc <- gather(DE.results.DE2$auc, 'method', 'auc')
DE.results.DE2.fdr <- gather(DE.results.DE2$fdr[-grep('raw',names(DE.results.DE2$fdr))], 'method', 'fdr')
DE.results.DE2.tpr <- gather(DE.results.DE2$tpr[-grep('raw',names(DE.results.DE2$tpr))], 'method', 'tpr')
DE.results.DE5.auc <- gather(DE.results.DE5$auc, 'method', 'auc')
DE.results.DE5.fdr <- gather(DE.results.DE5$fdr[-grep('raw',names(DE.results.DE5$fdr))], 'method', 'fdr')
DE.results.DE5.tpr <- gather(DE.results.DE5$tpr[-grep('raw',names(DE.results.DE5$tpr))], 'method', 'tpr')
DE.results.DE10.auc <- gather(DE.results.DE10$auc, 'method', 'auc')
DE.results.DE10.fdr <- gather(DE.results.DE10$fdr[-grep('raw',names(DE.results.DE10$fdr))], 'method', 'fdr')
DE.results.DE10.tpr <- gather(DE.results.DE10$tpr[-grep('raw',names(DE.results.DE10$tpr))], 'method', 'tpr')

outliers.DE.DE2.auc <- which(DE.results.DE2$auc$baySeq %in% sort(DE.results.DE2$auc$baySeq)[1:4])
colours.DE.DE2 <- rep('grey',50); colours.DE.DE2[outliers.DE.DE2.auc] = c('darkgreen','blue','yellow','darkorange')
outliers.DE.DE5.auc <- which(DE.results.DE5$auc$baySeq %in% sort(DE.results.DE5$auc$baySeq)[1:4])
colours.DE.DE5 <- rep('grey',50); colours.DE.DE5[outliers.DE.DE5.auc] = c('darkgreen','blue','yellow','darkorange')
outliers.DE.DE10.auc <- which(DE.results.DE10$auc$baySeq %in% sort(DE.results.DE10$auc$baySeq)[1:4])
colours.DE.DE10 <- rep('grey',50); colours.DE.DE10[outliers.DE.DE10.auc] = c('darkgreen','blue','yellow','darkorange')

ggplot(data = DE.results.DE2.auc, aes(x=method, y=auc)) + geom_point(col=rep(colours.DE.DE2,nrow(DE.results.DE2.auc)/50))
ggplot(data = DE.results.DE2.fdr, aes(x=method, y=fdr)) + geom_point(col=rep(colours.DE.DE2,nrow(DE.results.DE2.fdr)/50))
ggplot(data = DE.results.DE2.tpr, aes(x=method, y=tpr)) + geom_point(col=rep(colours.DE.DE2,nrow(DE.results.DE2.tpr)/50))
# Three obvious outliers for AUC for all methods except edgeR and voom, and a fourth less obvious: 47, 15, 23, 36.
# Often but not always among highest FDR and TPR.
ggplot(data = DE.results.DE5.auc, aes(x=method, y=auc)) + geom_point(col=rep(colours.DE.DE5,nrow(DE.results.DE5.auc)/50))
ggplot(data = DE.results.DE5.fdr, aes(x=method, y=fdr)) + geom_point(col=rep(colours.DE.DE5,nrow(DE.results.DE5.fdr)/50))
ggplot(data = DE.results.DE5.tpr, aes(x=method, y=tpr)) + geom_point(col=rep(colours.DE.DE5,nrow(DE.results.DE5.tpr)/50))
# Two extreme outliers for AUC for all methods except edgeR and voom: 45, 32. Usually most extreme FDRs, but not 
# always TPRs.
ggplot(data = DE.results.DE10.auc, aes(x=method, y=auc)) + geom_point(col=rep(colours.DE.DE10,nrow(DE.results.DE10.auc)/50))
ggplot(data = DE.results.DE10.fdr, aes(x=method, y=fdr)) + geom_point(col=rep(colours.DE.DE10,nrow(DE.results.DE10.fdr)/50))
ggplot(data = DE.results.DE10.tpr, aes(x=method, y=tpr)) + geom_point(col=rep(colours.DE.DE10,nrow(DE.results.DE10.tpr)/50))
# One extreme and three more obvious outliers for AUC for all methods except edgeR and voom: 14, 15, 16, 7.
# Also most extreme FDRs, but not TPRs.

# DE.lfc1 ####
DE.lfc1.results.DE2.auc <- gather(DE.lfc1.results.DE2$auc, 'method', 'auc')
DE.lfc1.results.DE2.fdr <- gather(DE.lfc1.results.DE2$fdr[-grep('raw',names(DE.lfc1.results.DE2$fdr))], 'method', 'fdr')
DE.lfc1.results.DE2.tpr <- gather(DE.lfc1.results.DE2$tpr[-grep('raw',names(DE.lfc1.results.DE2$tpr))], 'method', 'tpr')
DE.lfc1.results.DE5.auc <- gather(DE.lfc1.results.DE5$auc, 'method', 'auc')
DE.lfc1.results.DE5.fdr <- gather(DE.lfc1.results.DE5$fdr[-grep('raw',names(DE.lfc1.results.DE5$fdr))], 'method', 'fdr')
DE.lfc1.results.DE5.tpr <- gather(DE.lfc1.results.DE5$tpr[-grep('raw',names(DE.lfc1.results.DE5$tpr))], 'method', 'tpr')
DE.lfc1.results.DE10.auc <- gather(DE.lfc1.results.DE10$auc, 'method', 'auc')
DE.lfc1.results.DE10.fdr <- gather(DE.lfc1.results.DE10$fdr[-grep('raw',names(DE.lfc1.results.DE10$fdr))], 'method', 'fdr')
DE.lfc1.results.DE10.tpr <- gather(DE.lfc1.results.DE10$tpr[-grep('raw',names(DE.lfc1.results.DE10$tpr))], 'method', 'tpr')

outliers.DE.lfc1.DE2.auc <- which(DE.lfc1.results.DE2$auc$lfc1.DESeq %in% sort(DE.lfc1.results.DE2$auc$lfc1.DESeq)[1:4])
colours.DE.lfc1.DE2 <- rep('grey',50); colours.DE.lfc1.DE2[outliers.DE.lfc1.DE2.auc] = c('darkgreen','blue','yellow','darkorange')
outliers.DE.lfc1.DE5.auc <- which(DE.lfc1.results.DE5$auc$lfc1.DESeq %in% sort(DE.lfc1.results.DE5$auc$lfc1.DESeq)[1:4])
colours.DE.lfc1.DE5 <- rep('grey',50); colours.DE.lfc1.DE5[outliers.DE.lfc1.DE5.auc] = c('darkgreen','blue','yellow','darkorange')
outliers.DE.lfc1.DE10.auc <- which(DE.lfc1.results.DE10$auc$lfc1.DESeq %in% sort(DE.lfc1.results.DE10$auc$lfc1.DESeq)[1:4])
colours.DE.lfc1.DE10 <- rep('grey',50); colours.DE.lfc1.DE10[outliers.DE.lfc1.DE10.auc] = c('darkgreen','blue','yellow','darkorange')

ggplot(data = DE.lfc1.results.DE2.auc, aes(x=method, y=auc)) + geom_point(col=rep(colours.DE.lfc1.DE2,nrow(DE.lfc1.results.DE2.auc)/50))
ggplot(data = DE.lfc1.results.DE2.fdr, aes(x=method, y=fdr)) + geom_point(col=rep(colours.DE.lfc1.DE2,nrow(DE.lfc1.results.DE2.fdr)/50))
ggplot(data = DE.lfc1.results.DE2.tpr, aes(x=method, y=tpr)) + geom_point(col=rep(colours.DE.lfc1.DE2,nrow(DE.lfc1.results.DE2.tpr)/50))
# Outliers for AUC for DESeq2 match MDSeq but generally not HMs; 3/4 of datasets with lowest AUCs for DESeq2 have 
# higher AUCs for HMs than for DESeq2 or MDSeq, while datasets with lowest AUCs for HMs are around average for 
# DESeq2 and MDSeq. Most extreme outlier for HMs, which is also an outlier for DESeq2 and MDSeq, is 23. Most extreme 
# outlier for DESeq2 and MDSeq, which has higher AUC for HMs, is 47. Outlier for HMs which is average for others is 
# 9. No obvious outliers for FDR, TPR.
ggplot(data = DE.lfc1.results.DE5.auc, aes(x=method, y=auc)) + geom_point(col=rep(colours.DE.lfc1.DE5,nrow(DE.lfc1.results.DE5.auc)/50))
ggplot(data = DE.lfc1.results.DE5.fdr, aes(x=method, y=fdr)) + geom_point(col=rep(colours.DE.lfc1.DE5,nrow(DE.lfc1.results.DE5.fdr)/50))
ggplot(data = DE.lfc1.results.DE5.tpr, aes(x=method, y=tpr)) + geom_point(col=rep(colours.DE.lfc1.DE5,nrow(DE.lfc1.results.DE5.tpr)/50))
# One extreme outliers for AUC for DESeq2 and MDSeq, which has lowest AUC for edgeR and among highest AUCs for HMs: 
# 45. Next lowest for DESeq2 and MDSeq is also low for HMs: 32. Lowest for HMs is average for others: 33. No 
# consistent outliers for FDR, TPR.
ggplot(data = DE.lfc1.results.DE10.auc, aes(x=method, y=auc)) + geom_point(col=rep(colours.DE.lfc1.DE10,nrow(DE.lfc1.results.DE10.auc)/50))
ggplot(data = DE.lfc1.results.DE10.fdr, aes(x=method, y=fdr)) + geom_point(col=rep(colours.DE.lfc1.DE10,nrow(DE.lfc1.results.DE10.fdr)/50))
ggplot(data = DE.lfc1.results.DE10.tpr, aes(x=method, y=tpr)) + geom_point(col=rep(colours.DE.lfc1.DE10,nrow(DE.lfc1.results.DE10.tpr)/50))
# Same pattern again for AUCs. Three obvious outliers for DESeq2 and MDSeq, which are among the highest for HMs: 16, 
# 7, 14. Lowest AUC for HMs is also among lowest for others: 49. No real outliers for FDR, TPR.

# DE.lfc2 ####
DE.lfc2.results.DE2.auc <- gather(DE.lfc2.results.DE2$auc, 'method', 'auc')
DE.lfc2.results.DE2.fdr <- gather(DE.lfc2.results.DE2$fdr[-grep('raw',names(DE.lfc2.results.DE2$fdr))], 'method', 'fdr')
DE.lfc2.results.DE2.tpr <- gather(DE.lfc2.results.DE2$tpr[-grep('raw',names(DE.lfc2.results.DE2$tpr))], 'method', 'tpr')
DE.lfc2.results.DE5.auc <- gather(DE.lfc2.results.DE5$auc, 'method', 'auc')
DE.lfc2.results.DE5.fdr <- gather(DE.lfc2.results.DE5$fdr[-grep('raw',names(DE.lfc2.results.DE5$fdr))], 'method', 'fdr')
DE.lfc2.results.DE5.tpr <- gather(DE.lfc2.results.DE5$tpr[-grep('raw',names(DE.lfc2.results.DE5$tpr))], 'method', 'tpr')
DE.lfc2.results.DE10.auc <- gather(DE.lfc2.results.DE10$auc, 'method', 'auc')
DE.lfc2.results.DE10.fdr <- gather(DE.lfc2.results.DE10$fdr[-grep('raw',names(DE.lfc2.results.DE10$fdr))], 'method', 'fdr')
DE.lfc2.results.DE10.tpr <- gather(DE.lfc2.results.DE10$tpr[-grep('raw',names(DE.lfc2.results.DE10$tpr))], 'method', 'tpr')

outliers.DE.lfc2.DE2.auc <- which(DE.lfc2.results.DE2$auc$lfc2.DESeq %in% sort(DE.lfc2.results.DE2$auc$lfc2.DESeq)[1:4])
colours.DE.lfc2.DE2 <- rep('grey',50); colours.DE.lfc2.DE2[outliers.DE.lfc2.DE2.auc] = c('darkgreen','blue','yellow','darkorange')
outliers.DE.lfc2.DE5.auc <- which(DE.lfc2.results.DE5$auc$lfc2.DESeq %in% sort(DE.lfc2.results.DE5$auc$lfc2.DESeq)[1:4])
colours.DE.lfc2.DE5 <- rep('grey',50); colours.DE.lfc2.DE5[outliers.DE.lfc2.DE5.auc] = c('darkgreen','blue','yellow','darkorange')
outliers.DE.lfc2.DE10.auc <- which(DE.lfc2.results.DE10$auc$lfc2.DESeq %in% sort(DE.lfc2.results.DE10$auc$lfc2.DESeq)[1:4])
colours.DE.lfc2.DE10 <- rep('grey',50); colours.DE.lfc2.DE10[outliers.DE.lfc2.DE10.auc] = c('darkgreen','blue','yellow','darkorange')

ggplot(data = DE.lfc2.results.DE2.auc, aes(x=method, y=auc)) + geom_point(col=rep(colours.DE.lfc2.DE2,nrow(DE.lfc2.results.DE2.auc)/50))
ggplot(data = DE.lfc2.results.DE2.fdr, aes(x=method, y=fdr)) + geom_point(col=rep(colours.DE.lfc2.DE2,nrow(DE.lfc2.results.DE2.fdr)/50))
ggplot(data = DE.lfc2.results.DE2.tpr, aes(x=method, y=tpr)) + geom_point(col=rep(colours.DE.lfc2.DE2,nrow(DE.lfc2.results.DE2.tpr)/50))
ggplot(data = DE.lfc2.results.DE5.auc, aes(x=method, y=auc)) + geom_point(col=rep(colours.DE.lfc2.DE5,nrow(DE.lfc2.results.DE5.auc)/50))
ggplot(data = DE.lfc2.results.DE5.fdr, aes(x=method, y=fdr)) + geom_point(col=rep(colours.DE.lfc2.DE5,nrow(DE.lfc2.results.DE5.fdr)/50))
ggplot(data = DE.lfc2.results.DE5.tpr, aes(x=method, y=tpr)) + geom_point(col=rep(colours.DE.lfc2.DE5,nrow(DE.lfc2.results.DE5.tpr)/50))
ggplot(data = DE.lfc2.results.DE10.auc, aes(x=method, y=auc)) + geom_point(col=rep(colours.DE.lfc2.DE10,nrow(DE.lfc2.results.DE10.auc)/50))
ggplot(data = DE.lfc2.results.DE10.fdr, aes(x=method, y=fdr)) + geom_point(col=rep(colours.DE.lfc2.DE10,nrow(DE.lfc2.results.DE10.fdr)/50))
ggplot(data = DE.lfc2.results.DE10.tpr, aes(x=method, y=tpr)) + geom_point(col=rep(colours.DE.lfc2.DE10,nrow(DE.lfc2.results.DE10.tpr)/50))
# Basically same pattern as for lfc1 - outliers match between DESeq2 and MDSeq, but not HMs.

# DEDD ####
DEDD.results.DE2.auc <- gather(DEDD.results.DE2$auc, 'method', 'auc')
DEDD.results.DE2.fdr <- gather(DEDD.results.DE2$fdr, 'method', 'fdr')
DEDD.results.DE2.tpr <- gather(DEDD.results.DE2$tpr, 'method', 'tpr')
DEDD.results.DE5.auc <- gather(DEDD.results.DE5$auc, 'method', 'auc')
DEDD.results.DE5.fdr <- gather(DEDD.results.DE5$fdr, 'method', 'fdr')
DEDD.results.DE5.tpr <- gather(DEDD.results.DE5$tpr, 'method', 'tpr')
DEDD.results.DE10.auc <- gather(DEDD.results.DE10$auc, 'method', 'auc')
DEDD.results.DE10.fdr <- gather(DEDD.results.DE10$fdr, 'method', 'fdr')
DEDD.results.DE10.tpr <- gather(DEDD.results.DE10$tpr, 'method', 'tpr')

outliers.DEDD.DE2.auc <- which(DEDD.results.DE2$auc$expHM %in% sort(DEDD.results.DE2$auc$expHM)[1:4])
colours.DEDD.DE2 <- rep('grey',50); colours.DEDD.DE2[outliers.DEDD.DE2.auc] = c('darkgreen','blue','yellow','darkorange')
outliers.DEDD.DE5.auc <- which(DEDD.results.DE5$auc$expHM %in% sort(DEDD.results.DE5$auc$expHM)[1:4])
colours.DEDD.DE5 <- rep('grey',50); colours.DEDD.DE5[outliers.DEDD.DE5.auc] = c('darkgreen','blue','yellow','darkorange')
outliers.DEDD.DE10.auc <- which(DEDD.results.DE10$auc$expHM %in% sort(DEDD.results.DE10$auc$expHM)[1:4])
colours.DEDD.DE10 <- rep('grey',50); colours.DEDD.DE10[outliers.DEDD.DE10.auc] = c('darkgreen','blue','yellow','darkorange')

ggplot(data = DEDD.results.DE2.auc, aes(x=method, y=auc)) + geom_point(col=rep(colours.DEDD.DE2,nrow(DEDD.results.DE2.auc)/50))
ggplot(data = DEDD.results.DE2.fdr, aes(x=method, y=fdr)) + geom_point(col=rep(colours.DEDD.DE2,nrow(DEDD.results.DE2.fdr)/50))
ggplot(data = DEDD.results.DE2.tpr, aes(x=method, y=tpr)) + geom_point(col=rep(colours.DEDD.DE2,nrow(DEDD.results.DE2.tpr)/50))
# No real outliers, datasets with lowest AUCs roughly match.
ggplot(data = DEDD.results.DE5.auc, aes(x=method, y=auc)) + geom_point(col=rep(colours.DEDD.DE5,nrow(DEDD.results.DE5.auc)/50))
ggplot(data = DEDD.results.DE5.fdr, aes(x=method, y=fdr)) + geom_point(col=rep(colours.DEDD.DE5,nrow(DEDD.results.DE5.fdr)/50))
ggplot(data = DEDD.results.DE5.tpr, aes(x=method, y=tpr)) + geom_point(col=rep(colours.DEDD.DE5,nrow(DEDD.results.DE5.tpr)/50))
# Three extreme and one other outlier for AUC, same for exp and ln: 45, 6, 38, 5. 45 also an outlier for DE.lfc1 for 
# DESeq2 and MDSeq, but not HMs, and for all methods except edgeR and voom for DE. No consistent outliers for FDR, 
# TPR.
ggplot(data = DEDD.results.DE10.auc, aes(x=method, y=auc)) + geom_point(col=rep(colours.DEDD.DE10,nrow(DEDD.results.DE10.auc)/50))
ggplot(data = DEDD.results.DE10.fdr, aes(x=method, y=fdr)) + geom_point(col=rep(colours.DEDD.DE10,nrow(DEDD.results.DE10.fdr)/50))
ggplot(data = DEDD.results.DE10.tpr, aes(x=method, y=tpr)) + geom_point(col=rep(colours.DEDD.DE10,nrow(DEDD.results.DE10.tpr)/50))
# One extreme and three clear outliers for AUC, same for exp and ln: 14, 7, 16, 15, all also outliers for DESeq2 and 
# MDSeq but not for HMs for DE.lfc1, and outliers for all methods except edgeR and voom for DE. 14 is also an extreme 
# outlier for FDR (still conservative for 0.5, but not for thr), and has highest TPR in most cases.



######################################
# Investigating outlier datasets ####

## DEDD10, DE ####
results.DEDD10.45 <- readRDS(here('Results/DEDD compcodeR data results July-Aug 2019', 'results.DEDD10.45.all.rds'))
names(results.DEDD10.45)
data.DEDD10.45 <- results.DEDD10.45$data
DE.DEDD10.45 <- results.DEDD10.45$DE
DD.DEDD10.45 <- results.DEDD10.45$DD
c(sum(DE.DEDD10.45==0 & DD.DEDD10.45==0), sum(DE.DEDD10.45==0 & DD.DEDD10.45==1), 
  sum(DE.DEDD10.45==1 & DD.DEDD10.45==0), sum(DE.DEDD10.45==1 & DD.DEDD10.45==1))
c(sum(DE.DEDD10.45==0 & DD.DEDD10.45==0), sum(DE.DEDD10.45==0 & DD.DEDD10.45==1), 
  sum(DE.DEDD10.45==1 & DD.DEDD10.45==0), sum(DE.DEDD10.45==1 & DD.DEDD10.45==1)) / length(DE.DEDD10.45)
# 13850 no DE or DD (85.0%)
# 815 DD only (5.0%)
# 829 DE only (5.1%)
# 804 both (4.9%)
# Total 16298 genes after filtering

# AUCs
p.edgeR.ql.DEDD10.45 <- results.DEDD10.45$p.ql.edgeR
p.edgeR.lr.DEDD10.45 <- results.DEDD10.45$p.lr.edgeR
p.edgeR.et.DEDD10.45 <- results.DEDD10.45$p.et.edgeR
p.baySeq.DEDD10.45 <- 1 - results.DEDD10.45$prob.baySeq
p.MDSeq.zi.DEDD10.45 <- results.DEDD10.45$p.mean.zi.MDSeq
p.MDSeq.nozi.DEDD10.45 <- results.DEDD10.45$p.mean.nozi.MDSeq
p.DESeq2.if.DEDD10.45 <- results.DEDD10.45$p.if.DESeq
p.DESeq2.noif.DEDD10.45 <- results.DEDD10.45$p.noif.DESeq
p.voom.DEDD10.45 <- results.DEDD10.45$p.voom
p.DSS.notrend.DEDD10.45 <- results.DEDD10.45$p.notrend.DSS
p.expHM.untr.DEDD10.45 <- results.DEDD10.45$p.mean.expHM
p.expHM.tr.DEDD10.45 <- results.DEDD10.45$p.lmean.expHM
p.lnHM.untr.DEDD10.45 <- results.DEDD10.45$p.mean.lnHM
p.lnHM.tr.DEDD10.45 <- results.DEDD10.45$p.lmean.lnHM
pred.edgeR.ql.DEDD10.45 <- prediction(1-p.edgeR.ql.DEDD10.45, DE.DEDD10.45)
pred.edgeR.lr.DEDD10.45 <- prediction(1-p.edgeR.lr.DEDD10.45, DE.DEDD10.45)
pred.edgeR.et.DEDD10.45 <- prediction(1-p.edgeR.et.DEDD10.45, DE.DEDD10.45)
pred.baySeq.DEDD10.45 <- prediction(1-p.baySeq.DEDD10.45, DE.DEDD10.45)
pred.MDSeq.zi.DEDD10.45 <- prediction(1-p.MDSeq.zi.DEDD10.45, DE.DEDD10.45)
pred.MDSeq.nozi.DEDD10.45 <- prediction(1-p.MDSeq.nozi.DEDD10.45, DE.DEDD10.45)
pred.DESeq2.if.DEDD10.45 <- prediction(1-p.DESeq2.if.DEDD10.45, DE.DEDD10.45)
pred.DESeq2.noif.DEDD10.45 <- prediction(1-p.DESeq2.noif.DEDD10.45, DE.DEDD10.45)
pred.voom.DEDD10.45 <- prediction(1-p.voom.DEDD10.45, DE.DEDD10.45)
pred.DSS.notrend.DEDD10.45 <- prediction(1-p.DSS.notrend.DEDD10.45, DE.DEDD10.45)
pred.expHM.untr.DEDD10.45 <- prediction(1-p.expHM.untr.DEDD10.45, DE.DEDD10.45)
pred.expHM.tr.DEDD10.45 <- prediction(1-p.expHM.tr.DEDD10.45, DE.DEDD10.45)
pred.lnHM.untr.DEDD10.45 <- prediction(1-p.lnHM.untr.DEDD10.45, DE.DEDD10.45)
pred.lnHM.tr.DEDD10.45 <- prediction(1-p.lnHM.tr.DEDD10.45, DE.DEDD10.45)
auc.edgeR.ql.DEDD10.45 <- performance(pred.edgeR.ql.DEDD10.45, measure='auc')@y.values[[1]]
auc.edgeR.lr.DEDD10.45 <- performance(pred.edgeR.lr.DEDD10.45, measure='auc')@y.values[[1]]
auc.edgeR.et.DEDD10.45 <- performance(pred.edgeR.et.DEDD10.45, measure='auc')@y.values[[1]]
auc.baySeq.DEDD10.45 <- performance(pred.baySeq.DEDD10.45, measure='auc')@y.values[[1]]
auc.MDSeq.zi.DEDD10.45 <- performance(pred.MDSeq.zi.DEDD10.45, measure='auc')@y.values[[1]]
auc.MDSeq.nozi.DEDD10.45 <- performance(pred.MDSeq.nozi.DEDD10.45, measure='auc')@y.values[[1]]
auc.DESeq2.if.DEDD10.45 <- performance(pred.DESeq2.if.DEDD10.45, measure='auc')@y.values[[1]]
auc.DESeq2.noif.DEDD10.45 <- performance(pred.DESeq2.noif.DEDD10.45, measure='auc')@y.values[[1]]
auc.voom.DEDD10.45 <- performance(pred.voom.DEDD10.45, measure='auc')@y.values[[1]]
auc.DSS.notrend.DEDD10.45 <- performance(pred.DSS.notrend.DEDD10.45, measure='auc')@y.values[[1]]
auc.expHM.untr.DEDD10.45 <- performance(pred.expHM.untr.DEDD10.45, measure='auc')@y.values[[1]]
auc.expHM.tr.DEDD10.45 <- performance(pred.expHM.tr.DEDD10.45, measure='auc')@y.values[[1]]
auc.lnHM.untr.DEDD10.45 <- performance(pred.lnHM.untr.DEDD10.45, measure='auc')@y.values[[1]]
auc.lnHM.tr.DEDD10.45 <- performance(pred.lnHM.tr.DEDD10.45, measure='auc')@y.values[[1]]
c(auc.edgeR.ql.DEDD10.45, auc.edgeR.lr.DEDD10.45, auc.edgeR.et.DEDD10.45, auc.baySeq.DEDD10.45, 
  auc.MDSeq.zi.DEDD10.45, auc.MDSeq.nozi.DEDD10.45, auc.DESeq2.if.DEDD10.45, auc.DESeq2.noif.DEDD10.45, 
  auc.voom.DEDD10.45, auc.DSS.notrend.DEDD10.45, auc.expHM.untr.DEDD10.45, auc.expHM.tr.DEDD10.45, 
  auc.lnHM.untr.DEDD10.45, auc.lnHM.tr.DEDD10.45)
# edgeR 0.94; voom 0.93; HMs, baySeq, DSS, DESeq2 0.86; MDSeq 0.85

# Confusion matrices
q.edgeR.ql.DEDD10.45 <- results.DEDD10.45$q.ql.edgeR
q.edgeR.lr.DEDD10.45 <- results.DEDD10.45$q.lr.edgeR
q.edgeR.et.DEDD10.45 <- results.DEDD10.45$q.et.edgeR
q.baySeq.DEDD10.45 <- results.DEDD10.45$q.baySeq
q.MDSeq.zi.DEDD10.45 <- results.DEDD10.45$q.mean.zi.MDSeq
q.MDSeq.nozi.DEDD10.45 <- results.DEDD10.45$q.mean.nozi.MDSeq
q.DESeq2.if.DEDD10.45 <- results.DEDD10.45$q.if.DESeq
q.DESeq2.noif.DEDD10.45 <- results.DEDD10.45$q.noif.DESeq
q.voom.DEDD10.45 <- results.DEDD10.45$q.voom
q.DSS.notrend.DEDD10.45 <- results.DEDD10.45$q.notrend.DSS
q.expHM.untr.DEDD10.45 <- as.vector(p.adjust(results.DEDD10.45$p.mean.expHM, method='BH'))
q.expHM.tr.DEDD10.45 <- as.vector(p.adjust(results.DEDD10.45$p.lmean.expHM, method='BH'))
q.lnHM.untr.DEDD10.45 <- as.vector(p.adjust(results.DEDD10.45$p.mean.lnHM, method='BH'))
q.lnHM.tr.DEDD10.45 <- as.vector(p.adjust(results.DEDD10.45$p.lmean.lnHM, method='BH'))
c(mean(q.edgeR.ql.DEDD10.45), mean(q.edgeR.lr.DEDD10.45), mean(q.edgeR.et.DEDD10.45), mean(q.baySeq.DEDD10.45), 
  mean(q.MDSeq.zi.DEDD10.45, na.rm=T), mean(q.MDSeq.nozi.DEDD10.45, na.rm=T), mean(q.DESeq2.if.DEDD10.45), 
  mean(q.DESeq2.noif.DEDD10.45),mean(q.voom.DEDD10.45), mean(q.DSS.notrend.DEDD10.45), mean(q.expHM.untr.DEDD10.45), 
  mean(q.expHM.tr.DEDD10.45), mean(q.lnHM.untr.DEDD10.45), mean(q.lnHM.tr.DEDD10.45))
# NA for MDSeq but all others have a mean, so there aren't genes being excluded by edgeR and voom that are kept in 
# by the others. Mean q-values are greater for edgeR and voom than the rest, but also higher for baySeq than the 
# others, and DSS in between. With NAs removed, MDSeq higher than edgeR and voom.
c(sum(is.na(q.MDSeq.zi.DEDD10.45)), sum(is.na(q.MDSeq.nozi.DEDD10.45)))
# 19 NAs for zi, 38 for nozi
for (i in c('edgeR.ql.DEDD10.45', 'edgeR.lr.DEDD10.45', 'edgeR.et.DEDD10.45', 'baySeq.DEDD10.45', 
            'MDSeq.zi.DEDD10.45', 'MDSeq.nozi.DEDD10.45', 'DESeq2.if.DEDD10.45', 'DESeq2.noif.DEDD10.45', 
            'voom.DEDD10.45', 'DSS.notrend.DEDD10.45', 'expHM.untr.DEDD10.45', 'expHM.tr.DEDD10.45', 
            'lnHM.untr.DEDD10.45', 'lnHM.tr.DEDD10.45')) {
  assign(paste0('call.',i), get(paste0('q.',i)) < 0.05)
}
c(mean(call.edgeR.ql.DEDD10.45), mean(call.edgeR.lr.DEDD10.45), mean(call.edgeR.et.DEDD10.45), 
  mean(call.baySeq.DEDD10.45),   mean(call.MDSeq.zi.DEDD10.45, na.rm=T), mean(call.MDSeq.nozi.DEDD10.45, na.rm=T), 
  mean(call.DESeq2.if.DEDD10.45), mean(call.DESeq2.noif.DEDD10.45),mean(call.voom.DEDD10.45), 
  mean(call.DSS.notrend.DEDD10.45), mean(call.expHM.untr.DEDD10.45), mean(call.expHM.tr.DEDD10.45), 
  mean(call.lnHM.untr.DEDD10.45), mean(call.lnHM.tr.DEDD10.45))
# Better performance of edgeR and voom doesn't seem to be explained by generally higher or lower calling of DE - some 
# of other methods call more genes DE (DESeq2), some less (baySeq); DSS and HMs similar to edgeR and voom.
for (i in c('edgeR.ql.DEDD10.45', 'edgeR.lr.DEDD10.45', 'edgeR.et.DEDD10.45', 'baySeq.DEDD10.45', 
            'MDSeq.zi.DEDD10.45', 'MDSeq.nozi.DEDD10.45', 'DESeq2.if.DEDD10.45', 'DESeq2.noif.DEDD10.45', 
            'voom.DEDD10.45', 'DSS.notrend.DEDD10.45', 'expHM.untr.DEDD10.45', 'expHM.tr.DEDD10.45', 
            'lnHM.untr.DEDD10.45', 'lnHM.tr.DEDD10.45')) {
  assign(paste0('conf.',i), confusionMatrix(factor(as.numeric(get(paste0('call.',i)))), factor(DE.DEDD10.45), 
                                            positive='1'))
}
c(unname(conf.edgeR.ql.DEDD10.45$overall['Accuracy']), unname(conf.edgeR.lr.DEDD10.45$overall['Accuracy']), 
  unname(conf.edgeR.et.DEDD10.45$overall['Accuracy']), unname(conf.baySeq.DEDD10.45$overall['Accuracy']), 
  unname(conf.MDSeq.zi.DEDD10.45$overall['Accuracy']), unname(conf.MDSeq.nozi.DEDD10.45$overall['Accuracy']), 
  unname(conf.DESeq2.if.DEDD10.45$overall['Accuracy']), unname(conf.DESeq2.noif.DEDD10.45$overall['Accuracy']), 
  unname(conf.voom.DEDD10.45$overall['Accuracy']), unname(conf.DSS.notrend.DEDD10.45$overall['Accuracy']), 
  unname(conf.expHM.untr.DEDD10.45$overall['Accuracy']), unname(conf.expHM.tr.DEDD10.45$overall['Accuracy']), 
  unname(conf.lnHM.untr.DEDD10.45$overall['Accuracy']), unname(conf.lnHM.tr.DEDD10.45$overall['Accuracy']))
# No big differences in overall accuracy. edgeR, voom 0.96; DSS 0.95; baySeq, HMs, MDSeq 0.94; DESeq2 0.93.
c(unname(conf.edgeR.ql.DEDD10.45$byClass['Sensitivity']), unname(conf.edgeR.lr.DEDD10.45$byClass['Sensitivity']), 
  unname(conf.edgeR.et.DEDD10.45$byClass['Sensitivity']), unname(conf.baySeq.DEDD10.45$byClass['Sensitivity']), 
  unname(conf.MDSeq.zi.DEDD10.45$byClass['Sensitivity']), unname(conf.MDSeq.nozi.DEDD10.45$byClass['Sensitivity']), 
  unname(conf.DESeq2.if.DEDD10.45$byClass['Sensitivity']), unname(conf.DESeq2.noif.DEDD10.45$byClass['Sensitivity']), 
  unname(conf.voom.DEDD10.45$byClass['Sensitivity']), unname(conf.DSS.notrend.DEDD10.45$byClass['Sensitivity']), 
  unname(conf.expHM.untr.DEDD10.45$byClass['Sensitivity']), unname(conf.expHM.tr.DEDD10.45$byClass['Sensitivity']), 
  unname(conf.lnHM.untr.DEDD10.45$byClass['Sensitivity']), unname(conf.lnHM.tr.DEDD10.45$byClass['Sensitivity']))
# TPR: edgeR 0.61-0.67; DESeq2 0.63; voom 0.60; DSS 0.59; MDSeq 0.55; HMs 0.50-0.56; baySeq 0.49.
c(unname(1-conf.edgeR.ql.DEDD10.45$byClass['Precision']), unname(1-conf.edgeR.lr.DEDD10.45$byClass['Precision']), 
  unname(1-conf.edgeR.et.DEDD10.45$byClass['Precision']), unname(1-conf.baySeq.DEDD10.45$byClass['Precision']), 
  unname(1-conf.MDSeq.zi.DEDD10.45$byClass['Precision']), unname(1-conf.MDSeq.nozi.DEDD10.45$byClass['Precision']), 
  unname(1-conf.DESeq2.if.DEDD10.45$byClass['Precision']), unname(1-conf.DESeq2.noif.DEDD10.45$byClass['Precision']), 
  unname(1-conf.voom.DEDD10.45$byClass['Precision']), unname(1-conf.DSS.notrend.DEDD10.45$byClass['Precision']), 
  unname(1-conf.expHM.untr.DEDD10.45$byClass['Precision']), unname(1-conf.expHM.tr.DEDD10.45$byClass['Precision']), 
  unname(1-conf.lnHM.untr.DEDD10.45$byClass['Precision']), unname(1-conf.lnHM.tr.DEDD10.45$byClass['Precision']))
# FDRs below 0.1 for edgeR, voom, DSS, baySeq; HMs 0.14-0.18; MDSeq 0.26; DESeq2 0.36.
c(unname(conf.edgeR.ql.DEDD10.45$byClass['Specificity']), unname(conf.edgeR.lr.DEDD10.45$byClass['Specificity']), 
  unname(conf.edgeR.et.DEDD10.45$byClass['Specificity']), unname(conf.baySeq.DEDD10.45$byClass['Specificity']), 
  unname(conf.MDSeq.zi.DEDD10.45$byClass['Specificity']), unname(conf.MDSeq.nozi.DEDD10.45$byClass['Specificity']), 
  unname(conf.DESeq2.if.DEDD10.45$byClass['Specificity']), unname(conf.DESeq2.noif.DEDD10.45$byClass['Specificity']), 
  unname(conf.voom.DEDD10.45$byClass['Specificity']), unname(conf.DSS.notrend.DEDD10.45$byClass['Specificity']), 
  unname(conf.expHM.untr.DEDD10.45$byClass['Specificity']), unname(conf.expHM.tr.DEDD10.45$byClass['Specificity']), 
  unname(conf.lnHM.untr.DEDD10.45$byClass['Specificity']), unname(conf.lnHM.tr.DEDD10.45$byClass['Specificity']))
# Specificity very high (i.e. FPR very low) for all. edgeR, voom, DSS, baySeq > 0.99; HMs 0.99; MDSeq 0.98; 
# DESeq2 0.96.
c(unname(conf.edgeR.ql.DEDD10.45$byClass['F1']), unname(conf.edgeR.lr.DEDD10.45$byClass['F1']), 
  unname(conf.edgeR.et.DEDD10.45$byClass['F1']), unname(conf.baySeq.DEDD10.45$byClass['F1']), 
  unname(conf.MDSeq.zi.DEDD10.45$byClass['F1']), unname(conf.MDSeq.nozi.DEDD10.45$byClass['F1']), 
  unname(conf.DESeq2.if.DEDD10.45$byClass['F1']), unname(conf.DESeq2.noif.DEDD10.45$byClass['F1']), 
  unname(conf.voom.DEDD10.45$byClass['F1']), unname(conf.DSS.notrend.DEDD10.45$byClass['F1']), 
  unname(conf.expHM.untr.DEDD10.45$byClass['F1']), unname(conf.expHM.tr.DEDD10.45$byClass['F1']), 
  unname(conf.lnHM.untr.DEDD10.45$byClass['F1']), unname(conf.lnHM.tr.DEDD10.45$byClass['F1']))
# F1: edgeR 0.75-0.77; voom 0.74; DSS 0.72; HMs 0.63-0.67; baySeq 0.64; MDSeq, DESeq2 0.63.
c(unname(conf.edgeR.ql.DEDD10.45$byClass['Balanced Accuracy']), unname(conf.edgeR.lr.DEDD10.45$byClass['Balanced Accuracy']), 
  unname(conf.edgeR.et.DEDD10.45$byClass['Balanced Accuracy']), unname(conf.baySeq.DEDD10.45$byClass['Balanced Accuracy']), 
  unname(conf.MDSeq.zi.DEDD10.45$byClass['Balanced Accuracy']), unname(conf.MDSeq.nozi.DEDD10.45$byClass['Balanced Accuracy']), 
  unname(conf.DESeq2.if.DEDD10.45$byClass['Balanced Accuracy']), unname(conf.DESeq2.noif.DEDD10.45$byClass['Balanced Accuracy']), 
  unname(conf.voom.DEDD10.45$byClass['Balanced Accuracy']), unname(conf.DSS.notrend.DEDD10.45$byClass['Balanced Accuracy']), 
  unname(conf.expHM.untr.DEDD10.45$byClass['Balanced Accuracy']), unname(conf.expHM.tr.DEDD10.45$byClass['Balanced Accuracy']), 
  unname(conf.lnHM.untr.DEDD10.45$byClass['Balanced Accuracy']), unname(conf.lnHM.tr.DEDD10.45$byClass['Balanced Accuracy']))
# Balanced accuracy: edgeR 0.80-0.83; voom 0.80; DESeq2, DSS 0.79; MDSeq 0.77; HMs 0.75-0.78; baySeq 0.74.

calls <- data.frame(edgeR.ql=call.edgeR.ql.DEDD10.45, edgeR.lr=call.edgeR.lr.DEDD10.45, 
                    edgeR.et=call.edgeR.et.DEDD10.45, baySeq=call.baySeq.DEDD10.45, MDSeq.zi=call.MDSeq.zi.DEDD10.45, 
                    MDSeq.nozi=call.MDSeq.nozi.DEDD10.45, DESeq2.if=call.DESeq2.if.DEDD10.45, 
                    DESeq2.noif=call.DESeq2.noif.DEDD10.45, voom=call.voom.DEDD10.45, 
                    DSS.notrend=call.DSS.notrend.DEDD10.45, expHM.untr=call.expHM.untr.DEDD10.45, 
                    expHM.tr=call.expHM.tr.DEDD10.45, lnHM.untr=call.lnHM.untr.DEDD10.45, 
                    lnHM.tr=call.lnHM.tr.DEDD10.45)
round(cor(calls, use='complete.obs'), 2)
# edgeR and voom most highly correlated. DESeq2.if and DESeq2.noif identical; other different versions of 
# the same methods > 0.9.
# Next most highly correlated are DSS with edgeR and voom, and MDSeq with HMs; then baySeq with MDSeq and 
# HMs, and MDSeq with DESeq2.


differing.results <- which((call.baySeq.DEDD10.45==call.MDSeq.zi.DEDD10.45) & (call.baySeq.DEDD10.45==call.MDSeq.nozi.DEDD10.45) &
                             (call.baySeq.DEDD10.45==call.DESeq2.if.DEDD10.45) & (call.baySeq.DEDD10.45==call.DESeq2.noif.DEDD10.45) & 
                             (call.baySeq.DEDD10.45==call.DSS.notrend.DEDD10.45) & (call.baySeq.DEDD10.45==call.expHM.untr.DEDD10.45) & 
                             (call.baySeq.DEDD10.45==call.expHM.tr.DEDD10.45) & (call.baySeq.DEDD10.45==call.lnHM.untr.DEDD10.45) & 
                             (call.baySeq.DEDD10.45==call.lnHM.tr.DEDD10.45) & 
                             (call.baySeq.DEDD10.45!=call.edgeR.ql.DEDD10.45) & (call.baySeq.DEDD10.45!=call.edgeR.lr.DEDD10.45) & 
                             (call.baySeq.DEDD10.45!=call.edgeR.et.DEDD10.45) & (call.baySeq.DEDD10.45!=call.voom.DEDD10.45))
sum(call.baySeq.DEDD10.45[differing.results]==DE.DEDD10.45[differing.results])
sum(call.baySeq.DEDD10.45[differing.results]!=DE.DEDD10.45[differing.results])
# 6 instances where others all right and edgeR and voom wrong, 25 other way around.
table(DE.DEDD10.45[differing.results])
# Overall 11 non-DE and 20 DE genes have differing results between the groups of methods.
differing.results.tab <- rbind(DE.DEDD10.45[differing.results], DD.DEDD10.45[differing.results], 
                               call.baySeq.DEDD10.45[differing.results], call.edgeR.ql.DEDD10.45[differing.results])
rownames(differing.results.tab) <- c('DE','DD','baySeq','edgeR')
colnames(differing.results.tab) <- 1:31
differing.results.tab
#  3    DE, no DD, baySeq right
# 10    DE, no DD,  edgeR right
# Total 13/829 instances = 0.016
#  0    DE,    DD, baySeq right
#  7    DE,    DD,  edgeR right
# Total 7/804 instances = 0.009
#  3 no DE, no DD, baySeq right
#  8 no DE, no DD,  edgeR right
# Total 11/13850 instances = 0.001
#  0 no DE,    DD, baySeq right
#  0 no DE,    DD,  edgeR right
# Total 0/815
# Most likely to disagree for true DE genes, whether DD or not.

# Look at genes correctly called DE by edgeR and vooom but not by others, where no DD.
# Any obvious similarities? Zeros? Outliers? Sample means different from true means?
DE.edgeR.not.others.noDD <- differing.results[4:13]
data.DEDD10.45@variable.annotations$truedispersions.S1[DE.edgeR.not.others.noDD]
data.DEDD10.45@variable.annotations$truedispersions.S2[DE.edgeR.not.others.noDD]
data.DEDD10.45@variable.annotations$truemeans.S1[DE.edgeR.not.others.noDD]
data.DEDD10.45@variable.annotations$truemeans.S2[DE.edgeR.not.others.noDD]
data.DEDD10.45@variable.annotations$truelog2foldchanges[DE.edgeR.not.others.noDD]
# All negative log FCs, i.e. expression lower in second group.
rbind(differing.results.tab, data.DEDD10.45@variable.annotations$downregulation[differing.results])
# All 17 genes  correctly called DE by edgeR and voom but not others are down in group 2.
# All 3 genes correctly called DE by others but not edgeR and voom are up in group 2.
# Should look at performances restricted to upregulated and downregulated genes, initially 
# only for genes with no DD.
# But since FDRs are high for most other methods, maybe I should be looking more at genes 
# that are incorrectly called DE by the other methods but not by edgeR. There are 8 of 
# those, all with no DE.


subset.DEup.DEDD10.45 <- which((data.DEDD10.45@variable.annotations$upregulation == 1 | 
                                 data.DEDD10.45@variable.annotations$differential.expression == 0) & 
                                 data.DEDD10.45@variable.annotations$differential.dispersion == 0)
subset.DEdown.DEDD10.45 <- which((data.DEDD10.45@variable.annotations$downregulation == 1 | 
                                    data.DEDD10.45@variable.annotations$differential.expression == 0) & 
                                   data.DEDD10.45@variable.annotations$differential.dispersion == 0)

table(DE.DEDD10.45[subset.DEup.DEDD10.45]) # 428 DE, total 14278; 3% DE
table(DE.DEDD10.45[subset.DEdown.DEDD10.45]) # 401 DE, total 14251; 3% DE
# Make sure subsetting has worked properly:
table(DD.DEDD10.45[subset.DEup.DEDD10.45])
table(DD.DEDD10.45[subset.DEdown.DEDD10.45])
# No DD as expected
table(data.DEDD10.45@variable.annotations$upregulation[subset.DEup.DEDD10.45])
table(data.DEDD10.45@variable.annotations$downregulation[subset.DEup.DEDD10.45])
# Upregulation matches total DE, no downregulation, as expected
table(data.DEDD10.45@variable.annotations$downregulation[subset.DEdown.DEDD10.45])
table(data.DEDD10.45@variable.annotations$upregulation[subset.DEdown.DEDD10.45])
# Downregulation matches total DE, no upregulation, as expected
DE.DEDD10.45.DEup <- DE.DEDD10.45[subset.DEup.DEDD10.45]
DE.DEDD10.45.DEdown <- DE.DEDD10.45[subset.DEdown.DEDD10.45]
for (i in c('edgeR.ql.DEDD10.45', 'edgeR.lr.DEDD10.45', 'edgeR.et.DEDD10.45', 'baySeq.DEDD10.45', 
            'MDSeq.zi.DEDD10.45', 'MDSeq.nozi.DEDD10.45', 'DESeq2.if.DEDD10.45', 'DESeq2.noif.DEDD10.45', 
            'voom.DEDD10.45', 'DSS.notrend.DEDD10.45', 'expHM.untr.DEDD10.45', 'expHM.tr.DEDD10.45', 
            'lnHM.untr.DEDD10.45', 'lnHM.tr.DEDD10.45')) {
  assign(paste0('p.',i,'.DEup'), get(paste0('p.',i))[subset.DEup.DEDD10.45])
  assign(paste0('p.',i,'.DEdown'), get(paste0('p.',i))[subset.DEdown.DEDD10.45])
  assign(paste0('q.',i,'.DEup'), get(paste0('q.',i))[subset.DEup.DEDD10.45])
  assign(paste0('q.',i,'.DEdown'), get(paste0('q.',i))[subset.DEdown.DEDD10.45])
}

# AUCs
for (i in c('edgeR.ql.DEDD10.45', 'edgeR.lr.DEDD10.45', 'edgeR.et.DEDD10.45', 'baySeq.DEDD10.45', 
            'MDSeq.zi.DEDD10.45', 'MDSeq.nozi.DEDD10.45', 'DESeq2.if.DEDD10.45', 'DESeq2.noif.DEDD10.45', 
            'voom.DEDD10.45', 'DSS.notrend.DEDD10.45', 'expHM.untr.DEDD10.45', 'expHM.tr.DEDD10.45', 
            'lnHM.untr.DEDD10.45', 'lnHM.tr.DEDD10.45')) {
  assign(paste0('pred.',i,'.DEup'), prediction(1-get(paste0('p.',i,'.DEup')), DE.DEDD10.45.DEup))
  assign(paste0('pred.',i,'.DEdown'), prediction(1-get(paste0('p.',i,'.DEdown')), DE.DEDD10.45.DEdown))
  assign(paste0('auc.',i,'.DEup'), performance(get(paste0('pred.',i,'.DEup')), measure='auc')@y.values[[1]])
  assign(paste0('auc.',i,'.DEdown'), performance(get(paste0('pred.',i,'.DEdown')), measure='auc')@y.values[[1]])
}
c(auc.edgeR.ql.DEDD10.45.DEup, auc.edgeR.lr.DEDD10.45.DEup, auc.edgeR.et.DEDD10.45.DEup, auc.baySeq.DEDD10.45.DEup, 
  auc.MDSeq.zi.DEDD10.45.DEup, auc.MDSeq.nozi.DEDD10.45.DEup, auc.DESeq2.if.DEDD10.45.DEup, 
  auc.DESeq2.noif.DEDD10.45.DEup, auc.voom.DEDD10.45.DEup, auc.DSS.notrend.DEDD10.45.DEup, 
  auc.expHM.untr.DEDD10.45.DEup, auc.expHM.tr.DEDD10.45.DEup, auc.lnHM.untr.DEDD10.45.DEup, 
  auc.lnHM.tr.DEDD10.45.DEup)
c(auc.edgeR.ql.DEDD10.45.DEdown, auc.edgeR.lr.DEDD10.45.DEdown, auc.edgeR.et.DEDD10.45.DEdown, 
  auc.baySeq.DEDD10.45.DEdown, auc.MDSeq.zi.DEDD10.45.DEdown, auc.MDSeq.nozi.DEDD10.45.DEdown, 
  auc.DESeq2.if.DEDD10.45.DEdown, auc.DESeq2.noif.DEDD10.45.DEdown, auc.voom.DEDD10.45.DEdown, 
  auc.DSS.notrend.DEDD10.45.DEdown, auc.expHM.untr.DEDD10.45.DEdown, auc.expHM.tr.DEDD10.45.DEdown, 
  auc.lnHM.untr.DEDD10.45.DEdown, auc.lnHM.tr.DEDD10.45.DEdown)
# All very similar for upregulated genes, all except edgeR and voom much worse but similar to each 
# other for downregulated genes.

# Test whether pattern for up/downregulated genes is general ####
results.DEDD10.49 <- readRDS(here('Results/DEDD compcodeR data results July-Aug 2019', 'results.DEDD10.49.all.rds'))
data.DEDD10.49 <- results.DEDD10.49$data
DE.DEDD10.49 <- results.DEDD10.49$DE
DD.DEDD10.49 <- results.DEDD10.49$DD
p.edgeR.ql.DEDD10.49 <- results.DEDD10.49$p.ql.edgeR
p.edgeR.lr.DEDD10.49 <- results.DEDD10.49$p.lr.edgeR
p.edgeR.et.DEDD10.49 <- results.DEDD10.49$p.et.edgeR
p.baySeq.DEDD10.49 <- 1 - results.DEDD10.49$prob.baySeq
p.MDSeq.zi.DEDD10.49 <- results.DEDD10.49$p.mean.zi.MDSeq
p.MDSeq.nozi.DEDD10.49 <- results.DEDD10.49$p.mean.nozi.MDSeq
p.DESeq2.if.DEDD10.49 <- results.DEDD10.49$p.if.DESeq
p.DESeq2.noif.DEDD10.49 <- results.DEDD10.49$p.noif.DESeq
p.voom.DEDD10.49 <- results.DEDD10.49$p.voom
p.DSS.notrend.DEDD10.49 <- results.DEDD10.49$p.notrend.DSS
p.expHM.untr.DEDD10.49 <- results.DEDD10.49$p.mean.expHM
p.expHM.tr.DEDD10.49 <- results.DEDD10.49$p.lmean.expHM
p.lnHM.untr.DEDD10.49 <- results.DEDD10.49$p.mean.lnHM
p.lnHM.tr.DEDD10.49 <- results.DEDD10.49$p.lmean.lnHM
subset.DEup.DEDD10.49 <- which((data.DEDD10.49@variable.annotations$upregulation == 1 | 
                                  data.DEDD10.49@variable.annotations$differential.expression == 0) & 
                                 data.DEDD10.49@variable.annotations$differential.dispersion == 0)
subset.DEdown.DEDD10.49 <- which((data.DEDD10.49@variable.annotations$downregulation == 1 | 
                                    data.DEDD10.49@variable.annotations$differential.expression == 0) & 
                                   data.DEDD10.49@variable.annotations$differential.dispersion == 0)
DE.DEDD10.49.DEup <- DE.DEDD10.49[subset.DEup.DEDD10.49]
DE.DEDD10.49.DEdown <- DE.DEDD10.49[subset.DEdown.DEDD10.49]
for (i in c('edgeR.ql.DEDD10.49', 'edgeR.lr.DEDD10.49', 'edgeR.et.DEDD10.49', 'baySeq.DEDD10.49', 
            'MDSeq.zi.DEDD10.49', 'MDSeq.nozi.DEDD10.49', 'DESeq2.if.DEDD10.49', 'DESeq2.noif.DEDD10.49', 
            'voom.DEDD10.49', 'DSS.notrend.DEDD10.49', 'expHM.untr.DEDD10.49', 'expHM.tr.DEDD10.49', 
            'lnHM.untr.DEDD10.49', 'lnHM.tr.DEDD10.49')) {
  assign(paste0('p.',i,'.DEup'), get(paste0('p.',i))[subset.DEup.DEDD10.49])
  assign(paste0('p.',i,'.DEdown'), get(paste0('p.',i))[subset.DEdown.DEDD10.49])
}
for (i in c('edgeR.ql.DEDD10.49', 'edgeR.lr.DEDD10.49', 'edgeR.et.DEDD10.49', 'baySeq.DEDD10.49', 
            'MDSeq.zi.DEDD10.49', 'MDSeq.nozi.DEDD10.49', 'DESeq2.if.DEDD10.49', 'DESeq2.noif.DEDD10.49', 
            'voom.DEDD10.49', 'DSS.notrend.DEDD10.49', 'expHM.untr.DEDD10.49', 'expHM.tr.DEDD10.49', 
            'lnHM.untr.DEDD10.49', 'lnHM.tr.DEDD10.49')) {
  assign(paste0('pred.',i,'.DEup'), prediction(1-get(paste0('p.',i,'.DEup')), DE.DEDD10.49.DEup))
  assign(paste0('pred.',i,'.DEdown'), prediction(1-get(paste0('p.',i,'.DEdown')), DE.DEDD10.49.DEdown))
  assign(paste0('auc.',i,'.DEup'), performance(get(paste0('pred.',i,'.DEup')), measure='auc')@y.values[[1]])
  assign(paste0('auc.',i,'.DEdown'), performance(get(paste0('pred.',i,'.DEdown')), measure='auc')@y.values[[1]])
}
c(auc.edgeR.ql.DEDD10.49.DEup, auc.edgeR.lr.DEDD10.49.DEup, auc.edgeR.et.DEDD10.49.DEup, auc.baySeq.DEDD10.49.DEup, 
  auc.MDSeq.zi.DEDD10.49.DEup, auc.MDSeq.nozi.DEDD10.49.DEup, auc.DESeq2.if.DEDD10.49.DEup, 
  auc.DESeq2.noif.DEDD10.49.DEup, auc.voom.DEDD10.49.DEup, auc.DSS.notrend.DEDD10.49.DEup, 
  auc.expHM.untr.DEDD10.49.DEup, auc.expHM.tr.DEDD10.49.DEup, auc.lnHM.untr.DEDD10.49.DEup, 
  auc.lnHM.tr.DEDD10.49.DEup)
c(auc.edgeR.ql.DEDD10.49.DEdown, auc.edgeR.lr.DEDD10.49.DEdown, auc.edgeR.et.DEDD10.49.DEdown, 
  auc.baySeq.DEDD10.49.DEdown, auc.MDSeq.zi.DEDD10.49.DEdown, auc.MDSeq.nozi.DEDD10.49.DEdown, 
  auc.DESeq2.if.DEDD10.49.DEdown, auc.DESeq2.noif.DEDD10.49.DEdown, auc.voom.DEDD10.49.DEdown, 
  auc.DSS.notrend.DEDD10.49.DEdown, auc.expHM.untr.DEDD10.49.DEdown, auc.expHM.tr.DEDD10.49.DEdown, 
  auc.lnHM.untr.DEDD10.49.DEdown, auc.lnHM.tr.DEDD10.49.DEdown)
# Doesn't seem to be a general property but should investigate more thoroughly to make sure - i.e. check whether 
# AUC averaged over the 50 datasets is roughly the same for upregulated and downregulated genes.
# If it's not a general pattern, then it must just be that the downregulated genes in DEDD10.45 happen to have 
# whatever properties it is that makes edgeR and voom work much better than all the other methods.


for (j in 1:50) {
  assign(paste0('results.DEDD10.',j), readRDS(here('Results/DEDD compcodeR data results July-Aug 2019',
                                                   paste0('results.DEDD10.',j,'.all.rds'))))
  assign(paste0('data.DEDD10.',j), get('data', get(paste0('results.DEDD10.',j))))
  assign(paste0('DE.DEDD10.',j), get('DE', get(paste0('results.DEDD10.',j))))
  assign(paste0('DD.DEDD10.',j), get('DD', get(paste0('results.DEDD10.',j))))
  assign(paste0('p.edgeR.ql.DEDD10.',j), get('p.ql.edgeR', get(paste0('results.DEDD10.',j))))
  assign(paste0('p.edgeR.lr.DEDD10.',j), get('p.lr.edgeR', get(paste0('results.DEDD10.',j))))
  assign(paste0('p.edgeR.et.DEDD10.',j), get('p.et.edgeR', get(paste0('results.DEDD10.',j))))
  assign(paste0('p.baySeq.DEDD10.',j), 1 - get('prob.baySeq', get(paste0('results.DEDD10.',j))))
  assign(paste0('p.MDSeq.zi.DEDD10.',j), get('p.mean.zi.MDSeq', get(paste0('results.DEDD10.',j))))
  assign(paste0('p.MDSeq.nozi.DEDD10.',j), get('p.mean.nozi.MDSeq', get(paste0('results.DEDD10.',j))))
  assign(paste0('p.DESeq2.if.DEDD10.',j), get('p.if.DESeq', get(paste0('results.DEDD10.',j))))
  assign(paste0('p.DESeq2.noif.DEDD10.',j), get('p.noif.DESeq', get(paste0('results.DEDD10.',j))))
  assign(paste0('p.voom.DEDD10.',j), get('p.voom', get(paste0('results.DEDD10.',j))))
  assign(paste0('p.DSS.notrend.DEDD10.',j), get('p.notrend.DSS', get(paste0('results.DEDD10.',j))))
  assign(paste0('p.expHM.untr.DEDD10.',j), get('p.mean.expHM', get(paste0('results.DEDD10.',j))))
  assign(paste0('p.expHM.tr.DEDD10.',j), get('p.lmean.expHM', get(paste0('results.DEDD10.',j))))
  assign(paste0('p.lnHM.untr.DEDD10.',j), get('p.mean.lnHM', get(paste0('results.DEDD10.',j))))
  assign(paste0('p.lnHM.tr.DEDD10.',j), get('p.lmean.lnHM', get(paste0('results.DEDD10.',j))))
  assign(paste0('subset.DEup.DEDD10.',j), which((get('upregulation',
                                                     slot(get(paste0('data.DEDD10.',j)), 'variable.annotations')) == 1 |
                                                   get('differential.expression',
                                                       slot(get(paste0('data.DEDD10.',j)), 'variable.annotations')) == 0) &
                                                  get('differential.dispersion',
                                                      slot(get(paste0('data.DEDD10.',j)), 'variable.annotations')) == 0))
  assign(paste0('subset.DEdown.DEDD10.',j), which((get('downregulation',
                                                     slot(get(paste0('data.DEDD10.',j)), 'variable.annotations')) == 1 |
                                                   get('differential.expression',
                                                       slot(get(paste0('data.DEDD10.',j)), 'variable.annotations')) == 0) &
                                                  get('differential.dispersion',
                                                      slot(get(paste0('data.DEDD10.',j)), 'variable.annotations')) == 0))
  assign(paste0('DE.DEDD10.',j,'.DEup'), get(paste0('DE.DEDD10.',j))[get(paste0('subset.DEup.DEDD10.',j))])
  assign(paste0('DE.DEDD10.',j,'.DEdown'), get(paste0('DE.DEDD10.',j))[get(paste0('subset.DEdown.DEDD10.',j))])
  for (i in c(paste0('edgeR.ql.DEDD10.',j), paste0('edgeR.lr.DEDD10.',j),paste0('edgeR.et.DEDD10.',j),
              paste0('baySeq.DEDD10.',j), paste0('MDSeq.zi.DEDD10.',j), paste0('MDSeq.nozi.DEDD10.',j),
              paste0('DESeq2.if.DEDD10.',j), paste0('DESeq2.noif.DEDD10.',j), paste0('voom.DEDD10.',j),
              paste0('DSS.notrend.DEDD10.',j), paste0('expHM.untr.DEDD10.',j), paste0('expHM.tr.DEDD10.',j),
              paste0('lnHM.untr.DEDD10.',j), paste0('lnHM.tr.DEDD10.',j))) {
    assign(paste0('p.',i,'.DEup'), get(paste0('p.',i))[get(paste0('subset.DEup.DEDD10.',j))])
    assign(paste0('p.',i,'.DEdown'), get(paste0('p.',i))[get(paste0('subset.DEdown.DEDD10.',j))])
  }
  for (i in c(paste0('edgeR.ql.DEDD10.',j), paste0('edgeR.lr.DEDD10.',j),paste0('edgeR.et.DEDD10.',j),
              paste0('baySeq.DEDD10.',j), paste0('MDSeq.zi.DEDD10.',j), paste0('MDSeq.nozi.DEDD10.',j),
              paste0('DESeq2.if.DEDD10.',j), paste0('DESeq2.noif.DEDD10.',j), paste0('voom.DEDD10.',j),
              paste0('DSS.notrend.DEDD10.',j), paste0('expHM.untr.DEDD10.',j), paste0('expHM.tr.DEDD10.',j),
              paste0('lnHM.untr.DEDD10.',j), paste0('lnHM.tr.DEDD10.',j))) {
    assign(paste0('pred.',i,'.DEup'), prediction(1-get(paste0('p.',i,'.DEup')), get(paste0('DE.DEDD10.',j,'.DEup'))))
    assign(paste0('pred.',i,'.DEdown'), prediction(1-get(paste0('p.',i,'.DEdown')), get(paste0('DE.DEDD10.',j,'.DEdown'))))
    assign(paste0('auc.',i,'.DEup'), performance(get(paste0('pred.',i,'.DEup')), measure='auc')@y.values[[1]])
    assign(paste0('auc.',i,'.DEdown'), performance(get(paste0('pred.',i,'.DEdown')), measure='auc')@y.values[[1]])
    rm(list=c(paste0('pred.',i,'.DEup'), paste0('pred.',i,'.DEdown')))
  }
  rm(list=c(paste0('results.DEDD10.',j), paste0('data.DEDD10.',j)))
}
auc.DEDD10.up <- data.frame(matrix(nrow=50, ncol=14))
names(auc.DEDD10.up) <- c('edgeR.ql', 'edgeR.lr', 'edgeR.et', 'baySeq', 'MDSeq.zi', 'MDSeq.nozi', 'DESeq2.if', 
                          'DESeq2.noif', 'voom', 'DSS.notrend', 'expHM.untr', 'expHM.tr', 'lnHM.untr', 'lnHM.tr')
auc.DEDD10.down <- auc.DEDD10.up
names(auc.DEDD10.down) <- names(auc.DEDD10.up)
for (i in 1:50) {
  for (j in names(auc.DEDD10.up)) {
    auc.DEDD10.up[[j]][i] <- get(paste0('auc.',j,'.DEDD10.',i,'.DEup'))
    auc.DEDD10.down[[j]][i] <- get(paste0('auc.',j,'.DEDD10.',i,'.DEdown'))
  }
}
for (i in c('edgeR.ql', 'edgeR.lr', 'edgeR.et', 'baySeq', 'MDSeq.zi', 'MDSeq.nozi', 'DESeq2.if', 
            'DESeq2.noif', 'voom', 'DSS.notrend', 'expHM.untr', 'expHM.tr', 'lnHM.untr', 'lnHM.tr', 
            'DE.', 'DD.', 'subset')) {
  rm(list=ls()[grep(i, ls(), fixed=T)])
}
par(mfrow=c(2,1), mar=c(3,3,1,1))
boxplot(auc.DEDD10.up, ylim=c(0.75,1)); abline(h=c(0.8,0.85,0.9), col='grey')
boxplot(auc.DEDD10.down, ylim=c(0.75,1)); abline(h=c(0.8,0.85,0.9), col='grey')
round(rbind(colMeans(auc.DEDD10.up), colMeans(auc.DEDD10.down)),3)
round(rbind(apply(auc.DEDD10.up,2,median), apply(auc.DEDD10.down,2,median)),3)
round(rbind(apply(auc.DEDD10.up,2,max), apply(auc.DEDD10.down,2,max)),3)
round(rbind(apply(auc.DEDD10.up,2,min), apply(auc.DEDD10.down,2,min)),3)
# Generally means, medians and maxima similar between all methods, but minima higher for edgeR 
# and voom, and a lot higher for downregulated genes. Look to be four dataests that give much 
# lower AUCs for the other methods with downregulated genes.
outliers.DEDD10.up <- which(auc.DEDD10.up$baySeq < 0.9)
outliers.DEDD10.down <- which(auc.DEDD10.down$baySeq < 0.9)
extreme.outliers.DEDD10.up <- which(auc.DEDD10.up$baySeq < 0.87)
extreme.outliers.DEDD10.down <- which(auc.DEDD10.down$baySeq < 0.87)

# Need to see whether the difference in performance between upregulated and downregulated genes 
# is coincidence or not. Repeating on other datasets should tell.


# DEDD5 ####
for (j in 1:50) {
  assign(paste0('results.DEDD5.',j), readRDS(here('Results/DEDD compcodeR data results July-Aug 2019',
                                                   paste0('results.DEDD5.',j,'.all.rds'))))
  assign(paste0('data.DEDD5.',j), get('data', get(paste0('results.DEDD5.',j))))
  assign(paste0('DE.DEDD5.',j), get('DE', get(paste0('results.DEDD5.',j))))
  assign(paste0('DD.DEDD5.',j), get('DD', get(paste0('results.DEDD5.',j))))
  assign(paste0('p.edgeR.ql.DEDD5.',j), get('p.ql.edgeR', get(paste0('results.DEDD5.',j))))
  assign(paste0('p.edgeR.lr.DEDD5.',j), get('p.lr.edgeR', get(paste0('results.DEDD5.',j))))
  assign(paste0('p.edgeR.et.DEDD5.',j), get('p.et.edgeR', get(paste0('results.DEDD5.',j))))
  assign(paste0('p.baySeq.DEDD5.',j), 1 - get('prob.baySeq', get(paste0('results.DEDD5.',j))))
  assign(paste0('p.MDSeq.zi.DEDD5.',j), get('p.mean.zi.MDSeq', get(paste0('results.DEDD5.',j))))
  assign(paste0('p.MDSeq.nozi.DEDD5.',j), get('p.mean.nozi.MDSeq', get(paste0('results.DEDD5.',j))))
  assign(paste0('p.DESeq2.if.DEDD5.',j), get('p.if.DESeq', get(paste0('results.DEDD5.',j))))
  assign(paste0('p.DESeq2.noif.DEDD5.',j), get('p.noif.DESeq', get(paste0('results.DEDD5.',j))))
  assign(paste0('p.voom.DEDD5.',j), get('p.voom', get(paste0('results.DEDD5.',j))))
  assign(paste0('p.DSS.notrend.DEDD5.',j), get('p.notrend.DSS', get(paste0('results.DEDD5.',j))))
  assign(paste0('p.expHM.untr.DEDD5.',j), get('p.mean.expHM', get(paste0('results.DEDD5.',j))))
  assign(paste0('p.expHM.tr.DEDD5.',j), get('p.lmean.expHM', get(paste0('results.DEDD5.',j))))
  assign(paste0('p.lnHM.untr.DEDD5.',j), get('p.mean.lnHM', get(paste0('results.DEDD5.',j))))
  assign(paste0('p.lnHM.tr.DEDD5.',j), get('p.lmean.lnHM', get(paste0('results.DEDD5.',j))))
  assign(paste0('subset.DEup.DEDD5.',j), which((get('upregulation',
                                                     slot(get(paste0('data.DEDD5.',j)), 'variable.annotations')) == 1 |
                                                   get('differential.expression',
                                                       slot(get(paste0('data.DEDD5.',j)), 'variable.annotations')) == 0) &
                                                  get('differential.dispersion',
                                                      slot(get(paste0('data.DEDD5.',j)), 'variable.annotations')) == 0))
  assign(paste0('subset.DEdown.DEDD5.',j), which((get('downregulation',
                                                       slot(get(paste0('data.DEDD5.',j)), 'variable.annotations')) == 1 |
                                                     get('differential.expression',
                                                         slot(get(paste0('data.DEDD5.',j)), 'variable.annotations')) == 0) &
                                                    get('differential.dispersion',
                                                        slot(get(paste0('data.DEDD5.',j)), 'variable.annotations')) == 0))
  assign(paste0('DE.DEDD5.',j,'.DEup'), get(paste0('DE.DEDD5.',j))[get(paste0('subset.DEup.DEDD5.',j))])
  assign(paste0('DE.DEDD5.',j,'.DEdown'), get(paste0('DE.DEDD5.',j))[get(paste0('subset.DEdown.DEDD5.',j))])
  for (i in c(paste0('edgeR.ql.DEDD5.',j), paste0('edgeR.lr.DEDD5.',j),paste0('edgeR.et.DEDD5.',j),
              paste0('baySeq.DEDD5.',j), paste0('MDSeq.zi.DEDD5.',j), paste0('MDSeq.nozi.DEDD5.',j),
              paste0('DESeq2.if.DEDD5.',j), paste0('DESeq2.noif.DEDD5.',j), paste0('voom.DEDD5.',j),
              paste0('DSS.notrend.DEDD5.',j), paste0('expHM.untr.DEDD5.',j), paste0('expHM.tr.DEDD5.',j),
              paste0('lnHM.untr.DEDD5.',j), paste0('lnHM.tr.DEDD5.',j))) {
    assign(paste0('p.',i,'.DEup'), get(paste0('p.',i))[get(paste0('subset.DEup.DEDD5.',j))])
    assign(paste0('p.',i,'.DEdown'), get(paste0('p.',i))[get(paste0('subset.DEdown.DEDD5.',j))])
  }
  for (i in c(paste0('edgeR.ql.DEDD5.',j), paste0('edgeR.lr.DEDD5.',j),paste0('edgeR.et.DEDD5.',j),
              paste0('baySeq.DEDD5.',j), paste0('MDSeq.zi.DEDD5.',j), paste0('MDSeq.nozi.DEDD5.',j),
              paste0('DESeq2.if.DEDD5.',j), paste0('DESeq2.noif.DEDD5.',j), paste0('voom.DEDD5.',j),
              paste0('DSS.notrend.DEDD5.',j), paste0('expHM.untr.DEDD5.',j), paste0('expHM.tr.DEDD5.',j),
              paste0('lnHM.untr.DEDD5.',j), paste0('lnHM.tr.DEDD5.',j))) {
    assign(paste0('pred.',i,'.DEup'), prediction(1-get(paste0('p.',i,'.DEup')), get(paste0('DE.DEDD5.',j,'.DEup'))))
    assign(paste0('pred.',i,'.DEdown'), prediction(1-get(paste0('p.',i,'.DEdown')), get(paste0('DE.DEDD5.',j,'.DEdown'))))
    assign(paste0('auc.',i,'.DEup'), performance(get(paste0('pred.',i,'.DEup')), measure='auc')@y.values[[1]])
    assign(paste0('auc.',i,'.DEdown'), performance(get(paste0('pred.',i,'.DEdown')), measure='auc')@y.values[[1]])
    rm(list=c(paste0('pred.',i,'.DEup'), paste0('pred.',i,'.DEdown')))
  }
  rm(list=c(paste0('results.DEDD5.',j), paste0('data.DEDD5.',j)))
}
auc.DEDD5.up <- data.frame(matrix(nrow=50, ncol=14))
names(auc.DEDD5.up) <- c('edgeR.ql', 'edgeR.lr', 'edgeR.et', 'baySeq', 'MDSeq.zi', 'MDSeq.nozi', 'DESeq2.if', 
                          'DESeq2.noif', 'voom', 'DSS.notrend', 'expHM.untr', 'expHM.tr', 'lnHM.untr', 'lnHM.tr')
auc.DEDD5.down <- auc.DEDD5.up
names(auc.DEDD5.down) <- names(auc.DEDD5.up)
for (i in 1:50) {
  for (j in names(auc.DEDD5.up)) {
    auc.DEDD5.up[[j]][i] <- get(paste0('auc.',j,'.DEDD5.',i,'.DEup'))
    auc.DEDD5.down[[j]][i] <- get(paste0('auc.',j,'.DEDD5.',i,'.DEdown'))
  }
}
for (i in c('edgeR.ql', 'edgeR.lr', 'edgeR.et', 'baySeq', 'MDSeq.zi', 'MDSeq.nozi', 'DESeq2.if', 
            'DESeq2.noif', 'voom', 'DSS.notrend', 'expHM.untr', 'expHM.tr', 'lnHM.untr', 'lnHM.tr', 
            'DE.', 'DD.', 'subset')) {
  rm(list=ls()[grep(i, ls(), fixed=T)])
}
par(mfrow=c(2,1), mar=c(3,3,1,1))
boxplot(auc.DEDD5.up, ylim=c(0.6,1)); abline(h=c(0.7,0.75,0.8), col='grey')
boxplot(auc.DEDD5.down, ylim=c(0.6,1)); abline(h=c(0.7,0.75,0.8), col='grey')
round(rbind(colMeans(auc.DEDD5.up), colMeans(auc.DEDD5.down)),3)
round(rbind(apply(auc.DEDD5.up,2,median), apply(auc.DEDD5.down,2,median)),3)
round(rbind(apply(auc.DEDD5.up,2,max), apply(auc.DEDD5.down,2,max)),3)
round(rbind(apply(auc.DEDD5.up,2,min), apply(auc.DEDD5.down,2,min)),3)
# If anything the opposite pattern for DEDD5, so seems that it's just coincidence that all the 
# outliers are caused by downregulated genes in DEDD10.
outliers.DEDD5.up <- which(auc.DEDD5.up$baySeq < 0.8)
outliers.DEDD5.down <- which(auc.DEDD5.down$baySeq < 0.8)
extreme.outliers.DEDD5.up <- which(auc.DEDD5.up$baySeq < 0.75)
extreme.outliers.DEDD5.down <- which(auc.DEDD5.down$baySeq < 0.75)


# DEDD2 ####
for (j in 1:50) {
  assign(paste0('results.DEDD2.',j), readRDS(here('Results/DEDD compcodeR data results July-Aug 2019',
                                                  paste0('results.DEDD2.',j,'.all.rds'))))
  assign(paste0('data.DEDD2.',j), get('data', get(paste0('results.DEDD2.',j))))
  assign(paste0('DE.DEDD2.',j), get('DE', get(paste0('results.DEDD2.',j))))
  assign(paste0('DD.DEDD2.',j), get('DD', get(paste0('results.DEDD2.',j))))
  assign(paste0('p.edgeR.ql.DEDD2.',j), get('p.ql.edgeR', get(paste0('results.DEDD2.',j))))
  assign(paste0('p.edgeR.lr.DEDD2.',j), get('p.lr.edgeR', get(paste0('results.DEDD2.',j))))
  assign(paste0('p.edgeR.et.DEDD2.',j), get('p.et.edgeR', get(paste0('results.DEDD2.',j))))
  assign(paste0('p.baySeq.DEDD2.',j), 1 - get('prob.baySeq', get(paste0('results.DEDD2.',j))))
  assign(paste0('p.MDSeq.zi.DEDD2.',j), get('p.mean.zi.MDSeq', get(paste0('results.DEDD2.',j))))
  assign(paste0('p.MDSeq.nozi.DEDD2.',j), get('p.mean.nozi.MDSeq', get(paste0('results.DEDD2.',j))))
  assign(paste0('p.DESeq2.if.DEDD2.',j), get('p.if.DESeq', get(paste0('results.DEDD2.',j))))
  assign(paste0('p.DESeq2.noif.DEDD2.',j), get('p.noif.DESeq', get(paste0('results.DEDD2.',j))))
  assign(paste0('p.voom.DEDD2.',j), get('p.voom', get(paste0('results.DEDD2.',j))))
  assign(paste0('p.DSS.notrend.DEDD2.',j), get('p.notrend.DSS', get(paste0('results.DEDD2.',j))))
  assign(paste0('p.expHM.untr.DEDD2.',j), get('p.mean.expHM', get(paste0('results.DEDD2.',j))))
  assign(paste0('p.expHM.tr.DEDD2.',j), get('p.lmean.expHM', get(paste0('results.DEDD2.',j))))
  assign(paste0('p.lnHM.untr.DEDD2.',j), get('p.mean.lnHM', get(paste0('results.DEDD2.',j))))
  assign(paste0('p.lnHM.tr.DEDD2.',j), get('p.lmean.lnHM', get(paste0('results.DEDD2.',j))))
  assign(paste0('subset.DEup.DEDD2.',j), which((get('upregulation',
                                                    slot(get(paste0('data.DEDD2.',j)), 'variable.annotations')) == 1 |
                                                  get('differential.expression',
                                                      slot(get(paste0('data.DEDD2.',j)), 'variable.annotations')) == 0) &
                                                 get('differential.dispersion',
                                                     slot(get(paste0('data.DEDD2.',j)), 'variable.annotations')) == 0))
  assign(paste0('subset.DEdown.DEDD2.',j), which((get('downregulation',
                                                      slot(get(paste0('data.DEDD2.',j)), 'variable.annotations')) == 1 |
                                                    get('differential.expression',
                                                        slot(get(paste0('data.DEDD2.',j)), 'variable.annotations')) == 0) &
                                                   get('differential.dispersion',
                                                       slot(get(paste0('data.DEDD2.',j)), 'variable.annotations')) == 0))
  assign(paste0('DE.DEDD2.',j,'.DEup'), get(paste0('DE.DEDD2.',j))[get(paste0('subset.DEup.DEDD2.',j))])
  assign(paste0('DE.DEDD2.',j,'.DEdown'), get(paste0('DE.DEDD2.',j))[get(paste0('subset.DEdown.DEDD2.',j))])
  for (i in c(paste0('edgeR.ql.DEDD2.',j), paste0('edgeR.lr.DEDD2.',j),paste0('edgeR.et.DEDD2.',j),
              paste0('baySeq.DEDD2.',j), paste0('MDSeq.zi.DEDD2.',j), paste0('MDSeq.nozi.DEDD2.',j),
              paste0('DESeq2.if.DEDD2.',j), paste0('DESeq2.noif.DEDD2.',j), paste0('voom.DEDD2.',j),
              paste0('DSS.notrend.DEDD2.',j), paste0('expHM.untr.DEDD2.',j), paste0('expHM.tr.DEDD2.',j),
              paste0('lnHM.untr.DEDD2.',j), paste0('lnHM.tr.DEDD2.',j))) {
    assign(paste0('p.',i,'.DEup'), get(paste0('p.',i))[get(paste0('subset.DEup.DEDD2.',j))])
    assign(paste0('p.',i,'.DEdown'), get(paste0('p.',i))[get(paste0('subset.DEdown.DEDD2.',j))])
  }
  for (i in c(paste0('edgeR.ql.DEDD2.',j), paste0('edgeR.lr.DEDD2.',j),paste0('edgeR.et.DEDD2.',j),
              paste0('baySeq.DEDD2.',j), paste0('MDSeq.zi.DEDD2.',j), paste0('MDSeq.nozi.DEDD2.',j),
              paste0('DESeq2.if.DEDD2.',j), paste0('DESeq2.noif.DEDD2.',j), paste0('voom.DEDD2.',j),
              paste0('DSS.notrend.DEDD2.',j), paste0('expHM.untr.DEDD2.',j), paste0('expHM.tr.DEDD2.',j),
              paste0('lnHM.untr.DEDD2.',j), paste0('lnHM.tr.DEDD2.',j))) {
    assign(paste0('pred.',i,'.DEup'), prediction(1-get(paste0('p.',i,'.DEup')), get(paste0('DE.DEDD2.',j,'.DEup'))))
    assign(paste0('pred.',i,'.DEdown'), prediction(1-get(paste0('p.',i,'.DEdown')), get(paste0('DE.DEDD2.',j,'.DEdown'))))
    assign(paste0('auc.',i,'.DEup'), performance(get(paste0('pred.',i,'.DEup')), measure='auc')@y.values[[1]])
    assign(paste0('auc.',i,'.DEdown'), performance(get(paste0('pred.',i,'.DEdown')), measure='auc')@y.values[[1]])
    rm(list=c(paste0('pred.',i,'.DEup'), paste0('pred.',i,'.DEdown')))
  }
  rm(list=c(paste0('results.DEDD2.',j), paste0('data.DEDD2.',j)))
}
auc.DEDD2.up <- data.frame(matrix(nrow=50, ncol=14))
names(auc.DEDD2.up) <- c('edgeR.ql', 'edgeR.lr', 'edgeR.et', 'baySeq', 'MDSeq.zi', 'MDSeq.nozi', 'DESeq2.if', 
                         'DESeq2.noif', 'voom', 'DSS.notrend', 'expHM.untr', 'expHM.tr', 'lnHM.untr', 'lnHM.tr')
auc.DEDD2.down <- auc.DEDD2.up
names(auc.DEDD2.down) <- names(auc.DEDD2.up)
for (i in 1:50) {
  for (j in names(auc.DEDD2.up)) {
    auc.DEDD2.up[[j]][i] <- get(paste0('auc.',j,'.DEDD2.',i,'.DEup'))
    auc.DEDD2.down[[j]][i] <- get(paste0('auc.',j,'.DEDD2.',i,'.DEdown'))
  }
}
for (i in c('edgeR.ql', 'edgeR.lr', 'edgeR.et', 'baySeq', 'MDSeq.zi', 'MDSeq.nozi', 'DESeq2.if', 
            'DESeq2.noif', 'voom', 'DSS.notrend', 'expHM.untr', 'expHM.tr', 'lnHM.untr', 'lnHM.tr', 
            'DE.', 'DD.', 'subset')) {
  rm(list=ls()[grep(i, ls(), fixed=T)])
}
par(mfrow=c(2,1), mar=c(3,3,1,1))
boxplot(auc.DEDD2.up, ylim=c(0.5,0.9)); abline(h=c(0.6,0.65,0.7), col='grey')
boxplot(auc.DEDD2.down, ylim=c(0.5,0.9)); abline(h=c(0.6,0.65,0.7), col='grey')
round(rbind(colMeans(auc.DEDD2.up), colMeans(auc.DEDD2.down)),3)
round(rbind(apply(auc.DEDD2.up,2,median), apply(auc.DEDD2.down,2,median)),3)
round(rbind(apply(auc.DEDD2.up,2,max), apply(auc.DEDD2.down,2,max)),3)
round(rbind(apply(auc.DEDD2.up,2,min), apply(auc.DEDD2.down,2,min)),3)
# No difference for DEDD2.
outliers.DEDD2.up <- which(auc.DEDD2.up$baySeq < 0.7)
outliers.DEDD2.down <- which(auc.DEDD2.down$baySeq < 0.7)
extreme.outliers.DEDD2.up <- which(auc.DEDD2.up$baySeq < 0.6)
extreme.outliers.DEDD2.down <- which(auc.DEDD2.down$baySeq < 0.6)


# DE10 ####
for (j in 1:50) {
  assign(paste0('results.DE10.',j), readRDS(here('Results/DE compcodeR data results July-Aug 2019',
                                                   paste0('results.DE10.',j,'.all.rds'))))
  assign(paste0('data.DE10.',j), get('data', get(paste0('results.DE10.',j))))
  assign(paste0('DE.DE10.',j), get('DE', get(paste0('results.DE10.',j))))
  assign(paste0('p.edgeR.ql.DE10.',j), get('p.ql.edgeR', get(paste0('results.DE10.',j))))
  assign(paste0('p.edgeR.lr.DE10.',j), get('p.lr.edgeR', get(paste0('results.DE10.',j))))
  assign(paste0('p.edgeR.et.DE10.',j), get('p.et.edgeR', get(paste0('results.DE10.',j))))
  assign(paste0('p.baySeq.DE10.',j), 1 - get('prob.baySeq', get(paste0('results.DE10.',j))))
  assign(paste0('p.MDSeq.zi.DE10.',j), get('p.mean.zi.MDSeq', get(paste0('results.DE10.',j))))
  assign(paste0('p.MDSeq.nozi.DE10.',j), get('p.mean.nozi.MDSeq', get(paste0('results.DE10.',j))))
  assign(paste0('p.DESeq2.if.DE10.',j), get('p.if.DESeq', get(paste0('results.DE10.',j))))
  assign(paste0('p.DESeq2.noif.DE10.',j), get('p.noif.DESeq', get(paste0('results.DE10.',j))))
  assign(paste0('p.voom.DE10.',j), get('p.voom', get(paste0('results.DE10.',j))))
  assign(paste0('p.DSS.notrend.DE10.',j), get('p.notrend.DSS', get(paste0('results.DE10.',j))))
  assign(paste0('p.expHM.untr.DE10.',j), get('p.mean.expHM', get(paste0('results.DE10.',j))))
  assign(paste0('p.expHM.tr.DE10.',j), get('p.lmean.expHM', get(paste0('results.DE10.',j))))
  assign(paste0('p.lnHM.untr.DE10.',j), get('p.mean.lnHM', get(paste0('results.DE10.',j))))
  assign(paste0('p.lnHM.tr.DE10.',j), get('p.lmean.lnHM', get(paste0('results.DE10.',j))))
  assign(paste0('subset.DEup.DE10.',j), which((get('upregulation',
                                                     slot(get(paste0('data.DE10.',j)), 'variable.annotations')) == 1 |
                                                   get('differential.expression',
                                                       slot(get(paste0('data.DE10.',j)), 'variable.annotations')) == 0)))
  assign(paste0('subset.DEdown.DE10.',j), which((get('downregulation',
                                                       slot(get(paste0('data.DE10.',j)), 'variable.annotations')) == 1 |
                                                     get('differential.expression',
                                                         slot(get(paste0('data.DE10.',j)), 'variable.annotations')) == 0)))
  assign(paste0('DE.DE10.',j,'.DEup'), get(paste0('DE.DE10.',j))[get(paste0('subset.DEup.DE10.',j))])
  assign(paste0('DE.DE10.',j,'.DEdown'), get(paste0('DE.DE10.',j))[get(paste0('subset.DEdown.DE10.',j))])
  for (i in c(paste0('edgeR.ql.DE10.',j), paste0('edgeR.lr.DE10.',j),paste0('edgeR.et.DE10.',j),
              paste0('baySeq.DE10.',j), paste0('MDSeq.zi.DE10.',j), paste0('MDSeq.nozi.DE10.',j),
              paste0('DESeq2.if.DE10.',j), paste0('DESeq2.noif.DE10.',j), paste0('voom.DE10.',j),
              paste0('DSS.notrend.DE10.',j), paste0('expHM.untr.DE10.',j), paste0('expHM.tr.DE10.',j),
              paste0('lnHM.untr.DE10.',j), paste0('lnHM.tr.DE10.',j))) {
    assign(paste0('p.',i,'.DEup'), get(paste0('p.',i))[get(paste0('subset.DEup.DE10.',j))])
    assign(paste0('p.',i,'.DEdown'), get(paste0('p.',i))[get(paste0('subset.DEdown.DE10.',j))])
  }
  for (i in c(paste0('edgeR.ql.DE10.',j), paste0('edgeR.lr.DE10.',j),paste0('edgeR.et.DE10.',j),
              paste0('baySeq.DE10.',j), paste0('MDSeq.zi.DE10.',j), paste0('MDSeq.nozi.DE10.',j),
              paste0('DESeq2.if.DE10.',j), paste0('DESeq2.noif.DE10.',j), paste0('voom.DE10.',j),
              paste0('DSS.notrend.DE10.',j), paste0('expHM.untr.DE10.',j), paste0('expHM.tr.DE10.',j),
              paste0('lnHM.untr.DE10.',j), paste0('lnHM.tr.DE10.',j))) {
    assign(paste0('pred.',i,'.DEup'), prediction(1-get(paste0('p.',i,'.DEup')), get(paste0('DE.DE10.',j,'.DEup'))))
    assign(paste0('pred.',i,'.DEdown'), prediction(1-get(paste0('p.',i,'.DEdown')), get(paste0('DE.DE10.',j,'.DEdown'))))
    assign(paste0('auc.',i,'.DEup'), performance(get(paste0('pred.',i,'.DEup')), measure='auc')@y.values[[1]])
    assign(paste0('auc.',i,'.DEdown'), performance(get(paste0('pred.',i,'.DEdown')), measure='auc')@y.values[[1]])
    rm(list=c(paste0('pred.',i,'.DEup'), paste0('pred.',i,'.DEdown')))
  }
  rm(list=c(paste0('results.DE10.',j), paste0('data.DE10.',j)))
}
auc.DE10.up <- data.frame(matrix(nrow=50, ncol=14))
names(auc.DE10.up) <- c('edgeR.ql', 'edgeR.lr', 'edgeR.et', 'baySeq', 'MDSeq.zi', 'MDSeq.nozi', 'DESeq2.if', 
                          'DESeq2.noif', 'voom', 'DSS.notrend', 'expHM.untr', 'expHM.tr', 'lnHM.untr', 'lnHM.tr')
auc.DE10.down <- auc.DE10.up
names(auc.DE10.down) <- names(auc.DE10.up)
for (i in 1:50) {
  for (j in names(auc.DE10.up)) {
    auc.DE10.up[[j]][i] <- get(paste0('auc.',j,'.DE10.',i,'.DEup'))
    auc.DE10.down[[j]][i] <- get(paste0('auc.',j,'.DE10.',i,'.DEdown'))
  }
}
for (i in c('edgeR.ql', 'edgeR.lr', 'edgeR.et', 'baySeq', 'MDSeq.zi', 'MDSeq.nozi', 'DESeq2.if', 
            'DESeq2.noif', 'voom', 'DSS.notrend', 'expHM.untr', 'expHM.tr', 'lnHM.untr', 'lnHM.tr', 
            'DE.', 'DD.', 'subset')) {
  rm(list=ls()[grep(i, ls(), fixed=T)])
}
par(mfrow=c(2,1), mar=c(3,3,1,1))
boxplot(auc.DE10.up, ylim=c(0.7,1)); abline(h=c(0.8,0.85,0.9), col='grey')
boxplot(auc.DE10.down, ylim=c(0.7,1)); abline(h=c(0.8,0.85,0.9), col='grey')
round(rbind(colMeans(auc.DE10.up), colMeans(auc.DE10.down)),3)
round(rbind(apply(auc.DE10.up,2,median), apply(auc.DE10.down,2,median)),3)
round(rbind(apply(auc.DE10.up,2,max), apply(auc.DE10.down,2,max)),3)
round(rbind(apply(auc.DE10.up,2,min), apply(auc.DE10.down,2,min)),3)
# More outliers for upregulated for DE10 but not much difference.
outliers.DE10.up <- which(auc.DE10.up$baySeq < 0.9)
outliers.DE10.down <- which(auc.DE10.down$baySeq < 0.9)
extreme.outliers.DE10.up <- which(auc.DE10.up$baySeq < 0.85)
extreme.outliers.DE10.down <- which(auc.DE10.down$baySeq < 0.85)


# DE5 ####
for (j in 1:50) {
  assign(paste0('results.DE5.',j), readRDS(here('Results/DE compcodeR data results July-Aug 2019',
                                                 paste0('results.DE5.',j,'.all.rds'))))
  assign(paste0('data.DE5.',j), get('data', get(paste0('results.DE5.',j))))
  assign(paste0('DE.DE5.',j), get('DE', get(paste0('results.DE5.',j))))
  assign(paste0('p.edgeR.ql.DE5.',j), get('p.ql.edgeR', get(paste0('results.DE5.',j))))
  assign(paste0('p.edgeR.lr.DE5.',j), get('p.lr.edgeR', get(paste0('results.DE5.',j))))
  assign(paste0('p.edgeR.et.DE5.',j), get('p.et.edgeR', get(paste0('results.DE5.',j))))
  assign(paste0('p.baySeq.DE5.',j), 1 - get('prob.baySeq', get(paste0('results.DE5.',j))))
  assign(paste0('p.MDSeq.zi.DE5.',j), get('p.mean.zi.MDSeq', get(paste0('results.DE5.',j))))
  assign(paste0('p.MDSeq.nozi.DE5.',j), get('p.mean.nozi.MDSeq', get(paste0('results.DE5.',j))))
  assign(paste0('p.DESeq2.if.DE5.',j), get('p.if.DESeq', get(paste0('results.DE5.',j))))
  assign(paste0('p.DESeq2.noif.DE5.',j), get('p.noif.DESeq', get(paste0('results.DE5.',j))))
  assign(paste0('p.voom.DE5.',j), get('p.voom', get(paste0('results.DE5.',j))))
  assign(paste0('p.DSS.notrend.DE5.',j), get('p.notrend.DSS', get(paste0('results.DE5.',j))))
  assign(paste0('p.expHM.untr.DE5.',j), get('p.mean.expHM', get(paste0('results.DE5.',j))))
  assign(paste0('p.expHM.tr.DE5.',j), get('p.lmean.expHM', get(paste0('results.DE5.',j))))
  assign(paste0('p.lnHM.untr.DE5.',j), get('p.mean.lnHM', get(paste0('results.DE5.',j))))
  assign(paste0('p.lnHM.tr.DE5.',j), get('p.lmean.lnHM', get(paste0('results.DE5.',j))))
  assign(paste0('subset.DEup.DE5.',j), which((get('upregulation',
                                                   slot(get(paste0('data.DE5.',j)), 'variable.annotations')) == 1 |
                                                 get('differential.expression',
                                                     slot(get(paste0('data.DE5.',j)), 'variable.annotations')) == 0)))
  assign(paste0('subset.DEdown.DE5.',j), which((get('downregulation',
                                                     slot(get(paste0('data.DE5.',j)), 'variable.annotations')) == 1 |
                                                   get('differential.expression',
                                                       slot(get(paste0('data.DE5.',j)), 'variable.annotations')) == 0)))
  assign(paste0('DE.DE5.',j,'.DEup'), get(paste0('DE.DE5.',j))[get(paste0('subset.DEup.DE5.',j))])
  assign(paste0('DE.DE5.',j,'.DEdown'), get(paste0('DE.DE5.',j))[get(paste0('subset.DEdown.DE5.',j))])
  for (i in c(paste0('edgeR.ql.DE5.',j), paste0('edgeR.lr.DE5.',j),paste0('edgeR.et.DE5.',j),
              paste0('baySeq.DE5.',j), paste0('MDSeq.zi.DE5.',j), paste0('MDSeq.nozi.DE5.',j),
              paste0('DESeq2.if.DE5.',j), paste0('DESeq2.noif.DE5.',j), paste0('voom.DE5.',j),
              paste0('DSS.notrend.DE5.',j), paste0('expHM.untr.DE5.',j), paste0('expHM.tr.DE5.',j),
              paste0('lnHM.untr.DE5.',j), paste0('lnHM.tr.DE5.',j))) {
    assign(paste0('p.',i,'.DEup'), get(paste0('p.',i))[get(paste0('subset.DEup.DE5.',j))])
    assign(paste0('p.',i,'.DEdown'), get(paste0('p.',i))[get(paste0('subset.DEdown.DE5.',j))])
  }
  for (i in c(paste0('edgeR.ql.DE5.',j), paste0('edgeR.lr.DE5.',j),paste0('edgeR.et.DE5.',j),
              paste0('baySeq.DE5.',j), paste0('MDSeq.zi.DE5.',j), paste0('MDSeq.nozi.DE5.',j),
              paste0('DESeq2.if.DE5.',j), paste0('DESeq2.noif.DE5.',j), paste0('voom.DE5.',j),
              paste0('DSS.notrend.DE5.',j), paste0('expHM.untr.DE5.',j), paste0('expHM.tr.DE5.',j),
              paste0('lnHM.untr.DE5.',j), paste0('lnHM.tr.DE5.',j))) {
    assign(paste0('pred.',i,'.DEup'), prediction(1-get(paste0('p.',i,'.DEup')), get(paste0('DE.DE5.',j,'.DEup'))))
    assign(paste0('pred.',i,'.DEdown'), prediction(1-get(paste0('p.',i,'.DEdown')), get(paste0('DE.DE5.',j,'.DEdown'))))
    assign(paste0('auc.',i,'.DEup'), performance(get(paste0('pred.',i,'.DEup')), measure='auc')@y.values[[1]])
    assign(paste0('auc.',i,'.DEdown'), performance(get(paste0('pred.',i,'.DEdown')), measure='auc')@y.values[[1]])
    rm(list=c(paste0('pred.',i,'.DEup'), paste0('pred.',i,'.DEdown')))
  }
  rm(list=c(paste0('results.DE5.',j), paste0('data.DE5.',j)))
}
auc.DE5.up <- data.frame(matrix(nrow=50, ncol=14))
names(auc.DE5.up) <- c('edgeR.ql', 'edgeR.lr', 'edgeR.et', 'baySeq', 'MDSeq.zi', 'MDSeq.nozi', 'DESeq2.if', 
                        'DESeq2.noif', 'voom', 'DSS.notrend', 'expHM.untr', 'expHM.tr', 'lnHM.untr', 'lnHM.tr')
auc.DE5.down <- auc.DE5.up
names(auc.DE5.down) <- names(auc.DE5.up)
for (i in 1:50) {
  for (j in names(auc.DE5.up)) {
    auc.DE5.up[[j]][i] <- get(paste0('auc.',j,'.DE5.',i,'.DEup'))
    auc.DE5.down[[j]][i] <- get(paste0('auc.',j,'.DE5.',i,'.DEdown'))
  }
}
for (i in c('edgeR.ql', 'edgeR.lr', 'edgeR.et', 'baySeq', 'MDSeq.zi', 'MDSeq.nozi', 'DESeq2.if', 
            'DESeq2.noif', 'voom', 'DSS.notrend', 'expHM.untr', 'expHM.tr', 'lnHM.untr', 'lnHM.tr', 
            'DE.', 'DD.', 'subset')) {
  rm(list=ls()[grep(i, ls(), fixed=T)])
}
par(mfrow=c(2,1), mar=c(3,3,1,1))
boxplot(auc.DE5.up, ylim=c(0.65,1)); abline(h=c(0.7,0.75,0.8), col='grey')
boxplot(auc.DE5.down, ylim=c(0.65,1)); abline(h=c(0.7,0.75,0.8), col='grey')
round(rbind(colMeans(auc.DE5.up), colMeans(auc.DE5.down)),3)
round(rbind(apply(auc.DE5.up,2,median), apply(auc.DE5.down,2,median)),3)
round(rbind(apply(auc.DE5.up,2,max), apply(auc.DE5.down,2,max)),3)
round(rbind(apply(auc.DE5.up,2,min), apply(auc.DE5.down,2,min)),3)
# No difference for DE5.
outliers.DE5.up <- which(auc.DE5.up$baySeq < 0.8)
outliers.DE5.down <- which(auc.DE5.down$baySeq < 0.8)
extreme.outliers.DE5.up <- which(auc.DE5.up$baySeq < 0.75)
extreme.outliers.DE5.down <- which(auc.DE5.down$baySeq < 0.75)


# DE2 ####
for (j in 1:50) {
  assign(paste0('results.DE2.',j), readRDS(here('Results/DE compcodeR data results July-Aug 2019',
                                                paste0('results.DE2.',j,'.all.rds'))))
  assign(paste0('data.DE2.',j), get('data', get(paste0('results.DE2.',j))))
  assign(paste0('DE.DE2.',j), get('DE', get(paste0('results.DE2.',j))))
  assign(paste0('p.edgeR.ql.DE2.',j), get('p.ql.edgeR', get(paste0('results.DE2.',j))))
  assign(paste0('p.edgeR.lr.DE2.',j), get('p.lr.edgeR', get(paste0('results.DE2.',j))))
  assign(paste0('p.edgeR.et.DE2.',j), get('p.et.edgeR', get(paste0('results.DE2.',j))))
  assign(paste0('p.baySeq.DE2.',j), 1 - get('prob.baySeq', get(paste0('results.DE2.',j))))
  assign(paste0('p.MDSeq.zi.DE2.',j), get('p.mean.zi.MDSeq', get(paste0('results.DE2.',j))))
  assign(paste0('p.MDSeq.nozi.DE2.',j), get('p.mean.nozi.MDSeq', get(paste0('results.DE2.',j))))
  assign(paste0('p.DESeq2.if.DE2.',j), get('p.if.DESeq', get(paste0('results.DE2.',j))))
  assign(paste0('p.DESeq2.noif.DE2.',j), get('p.noif.DESeq', get(paste0('results.DE2.',j))))
  assign(paste0('p.voom.DE2.',j), get('p.voom', get(paste0('results.DE2.',j))))
  assign(paste0('p.DSS.notrend.DE2.',j), get('p.notrend.DSS', get(paste0('results.DE2.',j))))
  assign(paste0('p.expHM.untr.DE2.',j), get('p.mean.expHM', get(paste0('results.DE2.',j))))
  assign(paste0('p.expHM.tr.DE2.',j), get('p.lmean.expHM', get(paste0('results.DE2.',j))))
  assign(paste0('p.lnHM.untr.DE2.',j), get('p.mean.lnHM', get(paste0('results.DE2.',j))))
  assign(paste0('p.lnHM.tr.DE2.',j), get('p.lmean.lnHM', get(paste0('results.DE2.',j))))
  assign(paste0('subset.DEup.DE2.',j), which((get('upregulation',
                                                  slot(get(paste0('data.DE2.',j)), 'variable.annotations')) == 1 |
                                                get('differential.expression',
                                                    slot(get(paste0('data.DE2.',j)), 'variable.annotations')) == 0)))
  assign(paste0('subset.DEdown.DE2.',j), which((get('downregulation',
                                                    slot(get(paste0('data.DE2.',j)), 'variable.annotations')) == 1 |
                                                  get('differential.expression',
                                                      slot(get(paste0('data.DE2.',j)), 'variable.annotations')) == 0)))
  assign(paste0('DE.DE2.',j,'.DEup'), get(paste0('DE.DE2.',j))[get(paste0('subset.DEup.DE2.',j))])
  assign(paste0('DE.DE2.',j,'.DEdown'), get(paste0('DE.DE2.',j))[get(paste0('subset.DEdown.DE2.',j))])
  for (i in c(paste0('edgeR.ql.DE2.',j), paste0('edgeR.lr.DE2.',j),paste0('edgeR.et.DE2.',j),
              paste0('baySeq.DE2.',j), paste0('MDSeq.zi.DE2.',j), paste0('MDSeq.nozi.DE2.',j),
              paste0('DESeq2.if.DE2.',j), paste0('DESeq2.noif.DE2.',j), paste0('voom.DE2.',j),
              paste0('DSS.notrend.DE2.',j), paste0('expHM.untr.DE2.',j), paste0('expHM.tr.DE2.',j),
              paste0('lnHM.untr.DE2.',j), paste0('lnHM.tr.DE2.',j))) {
    assign(paste0('p.',i,'.DEup'), get(paste0('p.',i))[get(paste0('subset.DEup.DE2.',j))])
    assign(paste0('p.',i,'.DEdown'), get(paste0('p.',i))[get(paste0('subset.DEdown.DE2.',j))])
  }
  for (i in c(paste0('edgeR.ql.DE2.',j), paste0('edgeR.lr.DE2.',j),paste0('edgeR.et.DE2.',j),
              paste0('baySeq.DE2.',j), paste0('MDSeq.zi.DE2.',j), paste0('MDSeq.nozi.DE2.',j),
              paste0('DESeq2.if.DE2.',j), paste0('DESeq2.noif.DE2.',j), paste0('voom.DE2.',j),
              paste0('DSS.notrend.DE2.',j), paste0('expHM.untr.DE2.',j), paste0('expHM.tr.DE2.',j),
              paste0('lnHM.untr.DE2.',j), paste0('lnHM.tr.DE2.',j))) {
    assign(paste0('pred.',i,'.DEup'), prediction(1-get(paste0('p.',i,'.DEup')), get(paste0('DE.DE2.',j,'.DEup'))))
    assign(paste0('pred.',i,'.DEdown'), prediction(1-get(paste0('p.',i,'.DEdown')), get(paste0('DE.DE2.',j,'.DEdown'))))
    assign(paste0('auc.',i,'.DEup'), performance(get(paste0('pred.',i,'.DEup')), measure='auc')@y.values[[1]])
    assign(paste0('auc.',i,'.DEdown'), performance(get(paste0('pred.',i,'.DEdown')), measure='auc')@y.values[[1]])
    rm(list=c(paste0('pred.',i,'.DEup'), paste0('pred.',i,'.DEdown')))
  }
  rm(list=c(paste0('results.DE2.',j), paste0('data.DE2.',j)))
}
auc.DE2.up <- data.frame(matrix(nrow=50, ncol=14))
names(auc.DE2.up) <- c('edgeR.ql', 'edgeR.lr', 'edgeR.et', 'baySeq', 'MDSeq.zi', 'MDSeq.nozi', 'DESeq2.if', 
                       'DESeq2.noif', 'voom', 'DSS.notrend', 'expHM.untr', 'expHM.tr', 'lnHM.untr', 'lnHM.tr')
auc.DE2.down <- auc.DE2.up
names(auc.DE2.down) <- names(auc.DE2.up)
for (i in 1:50) {
  for (j in names(auc.DE2.up)) {
    auc.DE2.up[[j]][i] <- get(paste0('auc.',j,'.DE2.',i,'.DEup'))
    auc.DE2.down[[j]][i] <- get(paste0('auc.',j,'.DE2.',i,'.DEdown'))
  }
}
for (i in c('edgeR.ql', 'edgeR.lr', 'edgeR.et', 'baySeq', 'MDSeq.zi', 'MDSeq.nozi', 'DESeq2.if', 
            'DESeq2.noif', 'voom', 'DSS.notrend', 'expHM.untr', 'expHM.tr', 'lnHM.untr', 'lnHM.tr', 
            'DE.', 'DD.', 'subset')) {
  rm(list=ls()[grep(i, ls(), fixed=T)])
}
par(mfrow=c(2,1), mar=c(3,3,1,1))
boxplot(auc.DE2.up, ylim=c(0.44,0.9)); abline(h=c(0.55,0.6,0.65), col='grey')
boxplot(auc.DE2.down, ylim=c(0.44,0.9)); abline(h=c(0.55,0.6,0.65), col='grey')
round(rbind(colMeans(auc.DE2.up), colMeans(auc.DE2.down)),3)
round(rbind(apply(auc.DE2.up,2,median), apply(auc.DE2.down,2,median)),3)
round(rbind(apply(auc.DE2.up,2,max), apply(auc.DE2.down,2,max)),3)
round(rbind(apply(auc.DE2.up,2,min), apply(auc.DE2.down,2,min)),3)
# No difference for DE2.
outliers.DE2.up <- which(auc.DE2.up$baySeq < 0.65)
outliers.DE2.down <- which(auc.DE2.down$baySeq < 0.65)
extreme.outliers.DE2.up <- which(auc.DE2.up$baySeq < 0.59)
extreme.outliers.DE2.down <- which(auc.DE2.down$baySeq < 0.59)



outliers.DEDD10.down[which(outliers.DEDD10.down %in% outliers.DEDD10.up)]
outliers.DEDD5.down[which(outliers.DEDD5.down %in% outliers.DEDD5.up)]
outliers.DEDD2.down[which(outliers.DEDD2.down %in% outliers.DEDD2.up)]
outliers.DE10.down[which(outliers.DE10.down %in% outliers.DE10.up)]
outliers.DE5.down[which(outliers.DE5.down %in% outliers.DE5.up)]
outliers.DE2.down[which(outliers.DE2.down %in% outliers.DE2.up)]
# No overlap between outliers for upregulated and downregulated genes in any dataset.
par(mfrow=c(4,1), mar=c(3,3,1,1))
boxplot(auc.DEDD10.up, ylim=c(0.75,1)); abline(h=c(0.8,0.85,0.9), col='grey')
boxplot(auc.DEDD10.up[outliers.DEDD10.down,], ylim=c(0.75,1)); abline(h=c(0.8,0.85,0.9), col='grey')
boxplot(auc.DEDD10.down, ylim=c(0.75,1)); abline(h=c(0.8,0.85,0.9), col='grey')
boxplot(auc.DEDD10.down[outliers.DEDD10.up,], ylim=c(0.75,1)); abline(h=c(0.8,0.85,0.9), col='grey')
par(mfrow=c(4,1), mar=c(3,3,1,1))
boxplot(auc.DEDD5.up, ylim=c(0.6,1)); abline(h=c(0.7,0.75,0.8), col='grey')
boxplot(auc.DEDD5.up[outliers.DEDD5.down,], ylim=c(0.6,1)); abline(h=c(0.7,0.75,0.8), col='grey')
boxplot(auc.DEDD5.down, ylim=c(0.6,1)); abline(h=c(0.7,0.75,0.8), col='grey')
boxplot(auc.DEDD5.down[outliers.DEDD5.up,], ylim=c(0.6,1)); abline(h=c(0.7,0.75,0.8), col='grey')
par(mfrow=c(4,1), mar=c(3,3,1,1))
boxplot(auc.DEDD2.up, ylim=c(0.5,0.9)); abline(h=c(0.6,0.65,0.7), col='grey')
boxplot(auc.DEDD2.up[outliers.DEDD2.down,], ylim=c(0.5,0.9)); abline(h=c(0.6,0.65,0.7), col='grey')
boxplot(auc.DEDD2.down, ylim=c(0.5,0.9)); abline(h=c(0.6,0.65,0.7), col='grey')
boxplot(auc.DEDD2.down[outliers.DEDD2.up,], ylim=c(0.5,0.9)); abline(h=c(0.6,0.65,0.7), col='grey')
par(mfrow=c(4,1), mar=c(3,3,1,1))
boxplot(auc.DE10.up, ylim=c(0.7,1)); abline(h=c(0.8,0.85,0.9), col='grey')
boxplot(auc.DE10.up[outliers.DE10.down,], ylim=c(0.7,1)); abline(h=c(0.8,0.85,0.9), col='grey')
boxplot(auc.DE10.down, ylim=c(0.7,1)); abline(h=c(0.8,0.85,0.9), col='grey')
boxplot(auc.DE10.down[outliers.DE10.up,], ylim=c(0.7,1)); abline(h=c(0.8,0.85,0.9), col='grey')
par(mfrow=c(4,1), mar=c(3,3,1,1))
boxplot(auc.DE5.up, ylim=c(0.65,1)); abline(h=c(0.7,0.75,0.8), col='grey')
boxplot(auc.DE5.up[outliers.DE5.down,], ylim=c(0.65,1)); abline(h=c(0.7,0.75,0.8), col='grey')
boxplot(auc.DE5.down, ylim=c(0.65,1)); abline(h=c(0.7,0.75,0.8), col='grey')
boxplot(auc.DE5.down[outliers.DE5.up,], ylim=c(0.65,1)); abline(h=c(0.7,0.75,0.8), col='grey')
par(mfrow=c(4,1), mar=c(3,3,1,1))
boxplot(auc.DE2.up, ylim=c(0.44,0.9)); abline(h=c(0.55,0.6,0.65), col='grey')
boxplot(auc.DE2.up[outliers.DE2.down,], ylim=c(0.44,0.9)); abline(h=c(0.55,0.6,0.65), col='grey')
boxplot(auc.DE2.down, ylim=c(0.44,0.9)); abline(h=c(0.55,0.6,0.65), col='grey')
boxplot(auc.DE2.down[outliers.DE2.up,], ylim=c(0.44,0.9)); abline(h=c(0.55,0.6,0.65), col='grey')
# Outliers for upregulated genes will make good controls for outliers for downregulated genes, 
# and vice versa; i.e. can use each as an example of data that all methods perform well on to 
# compare to data with DE in the opposite direction for which methods other than edgeR and voom 
# perform badly.

# First look at DE only datasets to allow analysis on basis of DE only without having to add the 
# extra step of removing DD genes.


# Load data for outliers for DE only datasets ####
for (i in outliers.DE10.up) {
  assign(paste0('results.DE10.',i,'_outlier.up'), readRDS(here('Results/DE compcodeR data results July-Aug 2019',
                                                               paste0('results.DE10.',i,'.all.rds'))))
  assign(paste0('DE.DE10.',i,'_outlier.up'), get('DE', get(paste0('results.DE10.',i,'_outlier.up'))))
  assign(paste0('data.DE10.',i,'_outlier.up'), get('data', get(paste0('results.DE10.',i,'_outlier.up'))))
  assign(paste0('DEup.DE10.',i,'_outlier.up'), get('upregulation', slot(get(paste0('data.DE10.',i,'_outlier.up')),
                                                                        'variable.annotations')))
  assign(paste0('DEdown.DE10.',i,'_outlier.up'), get('downregulation', slot(get(paste0('data.DE10.',i,'_outlier.up')),
                                                                            'variable.annotations')))
  assign(paste0('lfc.DE10.',i,'_outlier.up'), get('truelog2foldchanges', slot(get(paste0('data.DE10.',i,'_outlier.up')),
                                                                              'variable.annotations')))
  assign(paste0('truemeans1.DE10.',i,'_outlier.up'), get('truemeans.S1', slot(get(paste0('data.DE10.',i,'_outlier.up')),
                                                                              'variable.annotations')))
  assign(paste0('truemeans2.DE10.',i,'_outlier.up'), get('truemeans.S2', slot(get(paste0('data.DE10.',i,'_outlier.up')),
                                                                              'variable.annotations')))
  assign(paste0('truedisps.DE10.',i,'_outlier.up'), get('truedispersions.S1', slot(get(paste0('data.DE10.',i,'_outlier.up')),
                                                                                   'variable.annotations')))
  assign(paste0('M.DE10.',i,'_outlier.up'), get('M.value', slot(get(paste0('data.DE10.',i,'_outlier.up')),
                                                                                   'variable.annotations')))
  assign(paste0('A.DE10.',i,'_outlier.up'), get('A.value', slot(get(paste0('data.DE10.',i,'_outlier.up')),
                                                                'variable.annotations')))
  assign(paste0('counts.DE10.',i,'_outlier.up'), slot(get(paste0('data.DE10.',i,'_outlier.up')),'count.matrix'))
}
for (i in outliers.DE10.down) {
  assign(paste0('results.DE10.',i,'_outlier.down'), readRDS(here('Results/DE compcodeR data results July-Aug 2019',
                                                               paste0('results.DE10.',i,'.all.rds'))))
  assign(paste0('DE.DE10.',i,'_outlier.down'), get('DE', get(paste0('results.DE10.',i,'_outlier.down'))))
  assign(paste0('data.DE10.',i,'_outlier.down'), get('data', get(paste0('results.DE10.',i,'_outlier.down'))))
  assign(paste0('DEup.DE10.',i,'_outlier.down'), get('upregulation', slot(get(paste0('data.DE10.',i,'_outlier.down')),
                                                                        'variable.annotations')))
  assign(paste0('DEdown.DE10.',i,'_outlier.down'), get('downregulation', slot(get(paste0('data.DE10.',i,'_outlier.down')),
                                                                            'variable.annotations')))
  assign(paste0('lfc.DE10.',i,'_outlier.down'), get('truelog2foldchanges', slot(get(paste0('data.DE10.',i,'_outlier.down')),
                                                                              'variable.annotations')))
  assign(paste0('truemeans1.DE10.',i,'_outlier.down'), get('truemeans.S1', slot(get(paste0('data.DE10.',i,'_outlier.down')),
                                                                              'variable.annotations')))
  assign(paste0('truemeans2.DE10.',i,'_outlier.down'), get('truemeans.S2', slot(get(paste0('data.DE10.',i,'_outlier.down')),
                                                                              'variable.annotations')))
  assign(paste0('truedisps.DE10.',i,'_outlier.down'), get('truedispersions.S1', slot(get(paste0('data.DE10.',i,'_outlier.down')),
                                                                                   'variable.annotations')))
  assign(paste0('M.DE10.',i,'_outlier.down'), get('M.value', slot(get(paste0('data.DE10.',i,'_outlier.down')),
                                                                'variable.annotations')))
  assign(paste0('A.DE10.',i,'_outlier.down'), get('A.value', slot(get(paste0('data.DE10.',i,'_outlier.down')),
                                                                'variable.annotations')))
  assign(paste0('counts.DE10.',i,'_outlier.down'), slot(get(paste0('data.DE10.',i,'_outlier.down')),'count.matrix'))
}
for (i in outliers.DE5.up) {
  assign(paste0('results.DE5.',i,'_outlier.up'), readRDS(here('Results/DE compcodeR data results July-Aug 2019',
                                                               paste0('results.DE5.',i,'.all.rds'))))
  assign(paste0('DE.DE5.',i,'_outlier.up'), get('DE', get(paste0('results.DE5.',i,'_outlier.up'))))
  assign(paste0('data.DE5.',i,'_outlier.up'), get('data', get(paste0('results.DE5.',i,'_outlier.up'))))
  assign(paste0('DEup.DE5.',i,'_outlier.up'), get('upregulation', slot(get(paste0('data.DE5.',i,'_outlier.up')),
                                                                        'variable.annotations')))
  assign(paste0('DEdown.DE5.',i,'_outlier.up'), get('downregulation', slot(get(paste0('data.DE5.',i,'_outlier.up')),
                                                                            'variable.annotations')))
  assign(paste0('lfc.DE5.',i,'_outlier.up'), get('truelog2foldchanges', slot(get(paste0('data.DE5.',i,'_outlier.up')),
                                                                              'variable.annotations')))
  assign(paste0('truemeans1.DE5.',i,'_outlier.up'), get('truemeans.S1', slot(get(paste0('data.DE5.',i,'_outlier.up')),
                                                                              'variable.annotations')))
  assign(paste0('truemeans2.DE5.',i,'_outlier.up'), get('truemeans.S2', slot(get(paste0('data.DE5.',i,'_outlier.up')),
                                                                              'variable.annotations')))
  assign(paste0('truedisps.DE5.',i,'_outlier.up'), get('truedispersions.S1', slot(get(paste0('data.DE5.',i,'_outlier.up')),
                                                                                   'variable.annotations')))
  assign(paste0('M.DE5.',i,'_outlier.up'), get('M.value', slot(get(paste0('data.DE5.',i,'_outlier.up')),
                                                                'variable.annotations')))
  assign(paste0('A.DE5.',i,'_outlier.up'), get('A.value', slot(get(paste0('data.DE5.',i,'_outlier.up')),
                                                                'variable.annotations')))
  assign(paste0('counts.DE5.',i,'_outlier.up'), slot(get(paste0('data.DE5.',i,'_outlier.up')),'count.matrix'))
}
for (i in outliers.DE5.down) {
  assign(paste0('results.DE5.',i,'_outlier.down'), readRDS(here('Results/DE compcodeR data results July-Aug 2019',
                                                                 paste0('results.DE5.',i,'.all.rds'))))
  assign(paste0('DE.DE5.',i,'_outlier.down'), get('DE', get(paste0('results.DE5.',i,'_outlier.down'))))
  assign(paste0('data.DE5.',i,'_outlier.down'), get('data', get(paste0('results.DE5.',i,'_outlier.down'))))
  assign(paste0('DEup.DE5.',i,'_outlier.down'), get('upregulation', slot(get(paste0('data.DE5.',i,'_outlier.down')),
                                                                          'variable.annotations')))
  assign(paste0('DEdown.DE5.',i,'_outlier.down'), get('downregulation', slot(get(paste0('data.DE5.',i,'_outlier.down')),
                                                                              'variable.annotations')))
  assign(paste0('lfc.DE5.',i,'_outlier.down'), get('truelog2foldchanges', slot(get(paste0('data.DE5.',i,'_outlier.down')),
                                                                                'variable.annotations')))
  assign(paste0('truemeans1.DE5.',i,'_outlier.down'), get('truemeans.S1', slot(get(paste0('data.DE5.',i,'_outlier.down')),
                                                                                'variable.annotations')))
  assign(paste0('truemeans2.DE5.',i,'_outlier.down'), get('truemeans.S2', slot(get(paste0('data.DE5.',i,'_outlier.down')),
                                                                                'variable.annotations')))
  assign(paste0('truedisps.DE5.',i,'_outlier.down'), get('truedispersions.S1', slot(get(paste0('data.DE5.',i,'_outlier.down')),
                                                                                     'variable.annotations')))
  assign(paste0('M.DE5.',i,'_outlier.down'), get('M.value', slot(get(paste0('data.DE5.',i,'_outlier.down')),
                                                                  'variable.annotations')))
  assign(paste0('A.DE5.',i,'_outlier.down'), get('A.value', slot(get(paste0('data.DE5.',i,'_outlier.down')),
                                                                  'variable.annotations')))
  assign(paste0('counts.DE5.',i,'_outlier.down'), slot(get(paste0('data.DE5.',i,'_outlier.down')),'count.matrix'))
}
for (i in outliers.DE2.up) {
  assign(paste0('results.DE2.',i,'_outlier.up'), readRDS(here('Results/DE compcodeR data results July-Aug 2019',
                                                               paste0('results.DE2.',i,'.all.rds'))))
  assign(paste0('DE.DE2.',i,'_outlier.up'), get('DE', get(paste0('results.DE2.',i,'_outlier.up'))))
  assign(paste0('data.DE2.',i,'_outlier.up'), get('data', get(paste0('results.DE2.',i,'_outlier.up'))))
  assign(paste0('DEup.DE2.',i,'_outlier.up'), get('upregulation', slot(get(paste0('data.DE2.',i,'_outlier.up')),
                                                                        'variable.annotations')))
  assign(paste0('DEdown.DE2.',i,'_outlier.up'), get('downregulation', slot(get(paste0('data.DE2.',i,'_outlier.up')),
                                                                            'variable.annotations')))
  assign(paste0('lfc.DE2.',i,'_outlier.up'), get('truelog2foldchanges', slot(get(paste0('data.DE2.',i,'_outlier.up')),
                                                                              'variable.annotations')))
  assign(paste0('truemeans1.DE2.',i,'_outlier.up'), get('truemeans.S1', slot(get(paste0('data.DE2.',i,'_outlier.up')),
                                                                              'variable.annotations')))
  assign(paste0('truemeans2.DE2.',i,'_outlier.up'), get('truemeans.S2', slot(get(paste0('data.DE2.',i,'_outlier.up')),
                                                                              'variable.annotations')))
  assign(paste0('truedisps.DE2.',i,'_outlier.up'), get('truedispersions.S1', slot(get(paste0('data.DE2.',i,'_outlier.up')),
                                                                                   'variable.annotations')))
  assign(paste0('M.DE2.',i,'_outlier.up'), get('M.value', slot(get(paste0('data.DE2.',i,'_outlier.up')),
                                                                'variable.annotations')))
  assign(paste0('A.DE2.',i,'_outlier.up'), get('A.value', slot(get(paste0('data.DE2.',i,'_outlier.up')),
                                                                'variable.annotations')))
  assign(paste0('counts.DE2.',i,'_outlier.up'), slot(get(paste0('data.DE2.',i,'_outlier.up')),'count.matrix'))
}
for (i in outliers.DE2.down) {
  assign(paste0('results.DE2.',i,'_outlier.down'), readRDS(here('Results/DE compcodeR data results July-Aug 2019',
                                                                 paste0('results.DE2.',i,'.all.rds'))))
  assign(paste0('DE.DE2.',i,'_outlier.down'), get('DE', get(paste0('results.DE2.',i,'_outlier.down'))))
  assign(paste0('data.DE2.',i,'_outlier.down'), get('data', get(paste0('results.DE2.',i,'_outlier.down'))))
  assign(paste0('DEup.DE2.',i,'_outlier.down'), get('upregulation', slot(get(paste0('data.DE2.',i,'_outlier.down')),
                                                                          'variable.annotations')))
  assign(paste0('DEdown.DE2.',i,'_outlier.down'), get('downregulation', slot(get(paste0('data.DE2.',i,'_outlier.down')),
                                                                              'variable.annotations')))
  assign(paste0('lfc.DE2.',i,'_outlier.down'), get('truelog2foldchanges', slot(get(paste0('data.DE2.',i,'_outlier.down')),
                                                                                'variable.annotations')))
  assign(paste0('truemeans1.DE2.',i,'_outlier.down'), get('truemeans.S1', slot(get(paste0('data.DE2.',i,'_outlier.down')),
                                                                                'variable.annotations')))
  assign(paste0('truemeans2.DE2.',i,'_outlier.down'), get('truemeans.S2', slot(get(paste0('data.DE2.',i,'_outlier.down')),
                                                                                'variable.annotations')))
  assign(paste0('truedisps.DE2.',i,'_outlier.down'), get('truedispersions.S1', slot(get(paste0('data.DE2.',i,'_outlier.down')),
                                                                                     'variable.annotations')))
  assign(paste0('M.DE2.',i,'_outlier.down'), get('M.value', slot(get(paste0('data.DE2.',i,'_outlier.down')),
                                                                  'variable.annotations')))
  assign(paste0('A.DE2.',i,'_outlier.down'), get('A.value', slot(get(paste0('data.DE2.',i,'_outlier.down')),
                                                                  'variable.annotations')))
  assign(paste0('counts.DE2.',i,'_outlier.down'), slot(get(paste0('data.DE2.',i,'_outlier.down')),'count.matrix'))
}

# Now have, for each DE only dataset and for outliers for upregulation and downregulation:
# - differentially expressed genes (DE.)
# - upregulated genes (DEup.)
# - downregulated genes (DEdown.)
# - log fold changes (lfc.)
# - true means in both groups (truemeans1., truemeans2.)
# - true dispersions (truedisps.)
# - M values (M.)
# - A values (A.)
# - count matrix (counts.)

# Number and proportion of up/down-regulated genes ####
for (i in outliers.DE10.up) {
  print(c(sum(get(paste0('DEup.DE10.',i,'_outlier.up'))), sum(get(paste0('DEdown.DE10.',i,'_outlier.up')))))
}
for (i in outliers.DE10.down) {
  print(c(sum(get(paste0('DEup.DE10.',i,'_outlier.down'))), sum(get(paste0('DEdown.DE10.',i,'_outlier.down')))))
}
for (i in outliers.DE5.up) {
  print(c(sum(get(paste0('DEup.DE5.',i,'_outlier.up'))), sum(get(paste0('DEdown.DE5.',i,'_outlier.up')))))
}
for (i in outliers.DE5.down) {
  print(c(sum(get(paste0('DEup.DE5.',i,'_outlier.down'))), sum(get(paste0('DEdown.DE5.',i,'_outlier.down')))))
}
for (i in outliers.DE2.up) {
  print(c(sum(get(paste0('DEup.DE2.',i,'_outlier.up'))), sum(get(paste0('DEdown.DE2.',i,'_outlier.up')))))
}
for (i in outliers.DE2.down) {
  print(c(sum(get(paste0('DEup.DE2.',i,'_outlier.down'))), sum(get(paste0('DEdown.DE2.',i,'_outlier.down')))))
}
# 415-457 upregulated, 359-411 downregulated
# Always more upregulated than downregulated

for (i in outliers.DE10.up) {
  print(c(sum(get(paste0('DEup.DE10.',i,'_outlier.up'))) / sum(get(paste0('DE.DE10.',i,'_outlier.up')))))
}
for (i in outliers.DE10.down) {
  print(c(sum(get(paste0('DEup.DE10.',i,'_outlier.down'))) / sum(get(paste0('DE.DE10.',i,'_outlier.down')))))
}
for (i in outliers.DE5.up) {
  print(c(sum(get(paste0('DEup.DE5.',i,'_outlier.up'))) / sum(get(paste0('DE.DE5.',i,'_outlier.up')))))
}
for (i in outliers.DE5.down) {
  print(c(sum(get(paste0('DEup.DE5.',i,'_outlier.down'))) / sum(get(paste0('DE.DE5.',i,'_outlier.down')))))
}
for (i in outliers.DE2.up) {
  print(c(sum(get(paste0('DEup.DE2.',i,'_outlier.up'))) / sum(get(paste0('DE.DE2.',i,'_outlier.up')))))
}
for (i in outliers.DE2.down) {
  print(c(sum(get(paste0('DEup.DE2.',i,'_outlier.down'))) / sum(get(paste0('DE.DE2.',i,'_outlier.down')))))
}
# 51.2-53.4% upregulated for datasets that are outliers for upregulation
# 51.5-54.4% upregulated for datasets that are outliers for downregulation

# Can't conclude any differences in numbers or proportions of up v downregulated genes, but looks like 
# downregulated genes generally more likely to be filtered out; presumably because, given random 
# distribution of means in group 1, those that are downregulated in group 2 are more likely to have low 
# counts and therefore more likely to be filtered out. The way the data is simulated also means that in 
# general, the distribution of means for DE genes will be different in group 2 from group 1 - greater 
# variance.


# Distribution of log fold changes ####
# Plot densities of LFC for each dataset. Non-outliers blue, outliers red.
length(outliers.DE10.up); length(outliers.DE10.down) # 12, 6
length(outliers.DE5.up); length(outliers.DE5.down) # 4, 7
length(outliers.DE2.up); length(outliers.DE2.down) # 11, 4

par(mfrow=c(8,6), mar=c(2,2,1,1))
for (i in outliers.DE10.up) {
  plot(density(0 - get(paste0('lfc.DE10.',i,'_outlier.up'))[which(
    get(paste0('DEdown.DE10.',i,'_outlier.up')) == 1)]), col='blue')
  lines(density(get(paste0('lfc.DE10.',i,'_outlier.up'))[which(
    get(paste0('DEup.DE10.',i,'_outlier.up')) == 1)]), col='red')
}
for (i in outliers.DE10.down) {
  plot(density(0 - get(paste0('lfc.DE10.',i,'_outlier.down'))[which(
    get(paste0('DEup.DE10.',i,'_outlier.down')) == 1)]), col='blue')
  lines(density(get(paste0('lfc.DE10.',i,'_outlier.down'))[which(
    get(paste0('DEdown.DE10.',i,'_outlier.down')) == 1)]), col='red')
}
for (i in outliers.DE5.up) {
  plot(density(0 - get(paste0('lfc.DE5.',i,'_outlier.up'))[which(
    get(paste0('DEdown.DE5.',i,'_outlier.up')) == 1)]), col='blue')
  lines(density(get(paste0('lfc.DE5.',i,'_outlier.up'))[which(
    get(paste0('DEup.DE5.',i,'_outlier.up')) == 1)]), col='red')
}
for (i in outliers.DE5.down) {
  plot(density(0 - get(paste0('lfc.DE5.',i,'_outlier.down'))[which(
    get(paste0('DEup.DE5.',i,'_outlier.down')) == 1)]), col='blue')
  lines(density(get(paste0('lfc.DE5.',i,'_outlier.down'))[which(
    get(paste0('DEdown.DE5.',i,'_outlier.down')) == 1)]), col='red')
}
for (i in outliers.DE2.up) {
  plot(density(0 - get(paste0('lfc.DE2.',i,'_outlier.up'))[which(
    get(paste0('DEdown.DE2.',i,'_outlier.up')) == 1)]), col='blue')
  lines(density(get(paste0('lfc.DE2.',i,'_outlier.up'))[which(
    get(paste0('DEup.DE2.',i,'_outlier.up')) == 1)]), col='red')
}
for (i in outliers.DE2.down) {
  plot(density(0 - get(paste0('lfc.DE2.',i,'_outlier.down'))[which(
    get(paste0('DEup.DE2.',i,'_outlier.down')) == 1)]), col='blue')
  lines(density(get(paste0('lfc.DE2.',i,'_outlier.down'))[which(
    get(paste0('DEdown.DE2.',i,'_outlier.down')) == 1)]), col='red')
}
# No noticeable differences at all between outlier and non-outlier datasets.


# M-A plots (actual LFC v log mean) - look for differences in patterns of DE ####
# Non-DE genes grey, genes DE in direction for which dataset is not outlier blue, genes 
# DE in direction for which dataset is an outlier red
par(mfrow=c(4,4), mar=c(1,1,0.5,0.5), mgp=c(3,0,0))
for (i in outliers.DE10.up) {
  plot(get(paste0('A.DE10.',i,'_outlier.up')), get(paste0('M.DE10.',i,'_outlier.up')), 
       pch=20, col='grey', ylim=c(-4.5,4.5), xlim=c(1,17))
  points(get(paste0('A.DE10.',i,'_outlier.up'))[which(get(paste0('DEdown.DE10.',i,'_outlier.up')) == 1)], 
         get(paste0('M.DE10.',i,'_outlier.up'))[which(get(paste0('DEdown.DE10.',i,'_outlier.up')) == 1)], 
         pch=20, col='blue')
  points(get(paste0('A.DE10.',i,'_outlier.up'))[which(get(paste0('DEup.DE10.',i,'_outlier.up')) == 1)], 
         get(paste0('M.DE10.',i,'_outlier.up'))[which(get(paste0('DEup.DE10.',i,'_outlier.up')) == 1)], 
         pch=20, col='red')
  abline(h=seq(-8,8,2), col='grey')
}
par(mfrow=c(4,4), mar=c(1,1,0.5,0.5), mgp=c(3,0,0))
for (i in outliers.DE10.down) {
  plot(get(paste0('A.DE10.',i,'_outlier.down')), get(paste0('M.DE10.',i,'_outlier.down')), 
       pch=20, col='grey', ylim=c(-4.5,4.5), xlim=c(1,17))
  points(get(paste0('A.DE10.',i,'_outlier.down'))[which(get(paste0('DEdown.DE10.',i,'_outlier.down')) == 1)], 
         get(paste0('M.DE10.',i,'_outlier.down'))[which(get(paste0('DEdown.DE10.',i,'_outlier.down')) == 1)], 
         pch=20, col='red')
  points(get(paste0('A.DE10.',i,'_outlier.down'))[which(get(paste0('DEup.DE10.',i,'_outlier.down')) == 1)], 
         get(paste0('M.DE10.',i,'_outlier.down'))[which(get(paste0('DEup.DE10.',i,'_outlier.down')) == 1)], 
         pch=20, col='blue')
  abline(h=seq(-8,8,2), col='grey')
}
par(mfrow=c(4,4), mar=c(1,1,0.5,0.5), mgp=c(3,0,0))
for (i in outliers.DE5.up) {
  plot(get(paste0('A.DE5.',i,'_outlier.up')), get(paste0('M.DE5.',i,'_outlier.up')), 
       pch=20, col='grey', ylim=c(-6,6), xlim=c(1,17))
  points(get(paste0('A.DE5.',i,'_outlier.up'))[which(get(paste0('DEdown.DE5.',i,'_outlier.up')) == 1)], 
         get(paste0('M.DE5.',i,'_outlier.up'))[which(get(paste0('DEdown.DE5.',i,'_outlier.up')) == 1)], 
         pch=20, col='blue')
  points(get(paste0('A.DE5.',i,'_outlier.up'))[which(get(paste0('DEup.DE5.',i,'_outlier.up')) == 1)], 
         get(paste0('M.DE5.',i,'_outlier.up'))[which(get(paste0('DEup.DE5.',i,'_outlier.up')) == 1)], 
         pch=20, col='red')
  abline(h=seq(-8,8,2), col='grey')
}
for (i in outliers.DE5.down) {
  plot(get(paste0('A.DE5.',i,'_outlier.down')), get(paste0('M.DE5.',i,'_outlier.down')), 
       pch=20, col='grey', ylim=c(-6,6), xlim=c(1,17))
  points(get(paste0('A.DE5.',i,'_outlier.down'))[which(get(paste0('DEdown.DE5.',i,'_outlier.down')) == 1)], 
         get(paste0('M.DE5.',i,'_outlier.down'))[which(get(paste0('DEdown.DE5.',i,'_outlier.down')) == 1)], 
         pch=20, col='red')
  points(get(paste0('A.DE5.',i,'_outlier.down'))[which(get(paste0('DEup.DE5.',i,'_outlier.down')) == 1)], 
         get(paste0('M.DE5.',i,'_outlier.down'))[which(get(paste0('DEup.DE5.',i,'_outlier.down')) == 1)], 
         pch=20, col='blue')
  abline(h=seq(-8,8,2), col='grey')
}
par(mfrow=c(4,4), mar=c(1,1,0.5,0.5), mgp=c(3,0,0))
for (i in outliers.DE2.up) {
  plot(get(paste0('A.DE2.',i,'_outlier.up')), get(paste0('M.DE2.',i,'_outlier.up')), 
       pch=20, col='grey', ylim=c(-10,10), xlim=c(1,17))
  points(get(paste0('A.DE2.',i,'_outlier.up'))[which(get(paste0('DEdown.DE2.',i,'_outlier.up')) == 1)], 
         get(paste0('M.DE2.',i,'_outlier.up'))[which(get(paste0('DEdown.DE2.',i,'_outlier.up')) == 1)], 
         pch=20, col='blue')
  points(get(paste0('A.DE2.',i,'_outlier.up'))[which(get(paste0('DEup.DE2.',i,'_outlier.up')) == 1)], 
         get(paste0('M.DE2.',i,'_outlier.up'))[which(get(paste0('DEup.DE2.',i,'_outlier.up')) == 1)], 
         pch=20, col='red')
  abline(h=seq(-8,8,2), col='grey')
}
for (i in outliers.DE2.down) {
  plot(get(paste0('A.DE2.',i,'_outlier.down')), get(paste0('M.DE2.',i,'_outlier.down')), 
       pch=20, col='grey', ylim=c(-10,10), xlim=c(1,17))
  points(get(paste0('A.DE2.',i,'_outlier.down'))[which(get(paste0('DEdown.DE2.',i,'_outlier.down')) == 1)], 
         get(paste0('M.DE2.',i,'_outlier.down'))[which(get(paste0('DEdown.DE2.',i,'_outlier.down')) == 1)], 
         pch=20, col='red')
  points(get(paste0('A.DE2.',i,'_outlier.down'))[which(get(paste0('DEup.DE2.',i,'_outlier.down')) == 1)], 
         get(paste0('M.DE2.',i,'_outlier.down'))[which(get(paste0('DEup.DE2.',i,'_outlier.down')) == 1)], 
         pch=20, col='blue')
  abline(h=seq(-8,8,2), col='grey')
}
# No obvious differences. Would expect maybe to see more non-DE genes with high observed 
# LFC or more DE genes with low observed LFC in outlier datasets, but that doesn't seem 
# to be the case.


# Look at differences between true and observed LFCs ####
# May expect outlier datasets to have less correlation between true and observed LFCs

# Plot true v (observed - true) LFC
par(mfrow=c(4,4), mar=c(1,1,0.5,0.5), mgp=c(3,0,0))
for (i in outliers.DE10.up) {
  plot(get(paste0('lfc.DE10.',i,'_outlier.up')), get(paste0('M.DE10.',i,'_outlier.up')) - 
         get(paste0('lfc.DE10.',i,'_outlier.up')), pch=20, col='grey', main='')
  points(get(paste0('lfc.DE10.',i,'_outlier.up'))[which(get(paste0('DEdown.DE10.',i,'_outlier.up')) == 1)], 
         get(paste0('M.DE10.',i,'_outlier.up'))[which(get(paste0('DEdown.DE10.',i,'_outlier.up')) == 1)] - 
           get(paste0('lfc.DE10.',i,'_outlier.up'))[which(get(paste0('DEdown.DE10.',i,'_outlier.up')) == 1)], 
         pch=20, col='blue')
  points(get(paste0('lfc.DE10.',i,'_outlier.up'))[which(get(paste0('DEup.DE10.',i,'_outlier.up')) == 1)], 
         get(paste0('M.DE10.',i,'_outlier.up'))[which(get(paste0('DEup.DE10.',i,'_outlier.up')) == 1)] - 
           get(paste0('lfc.DE10.',i,'_outlier.up'))[which(get(paste0('DEup.DE10.',i,'_outlier.up')) == 1)], 
         pch=20, col='red')
  lines(c(-4,4),c(0,0))
}
par(mfrow=c(4,4), mar=c(1,1,0.5,0.5), mgp=c(3,0,0))
for (i in outliers.DE10.down) {
  plot(get(paste0('lfc.DE10.',i,'_outlier.down')), get(paste0('M.DE10.',i,'_outlier.down')) - 
         get(paste0('lfc.DE10.',i,'_outlier.down')), pch=20, col='grey', main='')
  points(get(paste0('lfc.DE10.',i,'_outlier.down'))[which(get(paste0('DEdown.DE10.',i,'_outlier.down')) == 1)], 
         get(paste0('M.DE10.',i,'_outlier.down'))[which(get(paste0('DEdown.DE10.',i,'_outlier.down')) == 1)] - 
           get(paste0('lfc.DE10.',i,'_outlier.down'))[which(get(paste0('DEdown.DE10.',i,'_outlier.down')) == 1)], 
         pch=20, col='red')
  points(get(paste0('lfc.DE10.',i,'_outlier.down'))[which(get(paste0('DEup.DE10.',i,'_outlier.down')) == 1)], 
         get(paste0('M.DE10.',i,'_outlier.down'))[which(get(paste0('DEup.DE10.',i,'_outlier.down')) == 1)] - 
           get(paste0('lfc.DE10.',i,'_outlier.down'))[which(get(paste0('DEup.DE10.',i,'_outlier.down')) == 1)], 
         pch=20, col='blue')
  lines(c(-4,4),c(0,0))
}
par(mfrow=c(4,4), mar=c(1,1,0.5,0.5), mgp=c(3,0,0))
for (i in outliers.DE5.down) {
  plot(get(paste0('lfc.DE5.',i,'_outlier.down')), get(paste0('M.DE5.',i,'_outlier.down')) - 
         get(paste0('lfc.DE5.',i,'_outlier.down')), pch=20, col='grey', main='')
  points(get(paste0('lfc.DE5.',i,'_outlier.down'))[which(get(paste0('DEdown.DE5.',i,'_outlier.down')) == 1)], 
         get(paste0('M.DE5.',i,'_outlier.down'))[which(get(paste0('DEdown.DE5.',i,'_outlier.down')) == 1)] - 
           get(paste0('lfc.DE5.',i,'_outlier.down'))[which(get(paste0('DEdown.DE5.',i,'_outlier.down')) == 1)], 
         pch=20, col='red')
  points(get(paste0('lfc.DE5.',i,'_outlier.down'))[which(get(paste0('DEup.DE5.',i,'_outlier.down')) == 1)], 
         get(paste0('M.DE5.',i,'_outlier.down'))[which(get(paste0('DEup.DE5.',i,'_outlier.down')) == 1)] - 
           get(paste0('lfc.DE5.',i,'_outlier.down'))[which(get(paste0('DEup.DE5.',i,'_outlier.down')) == 1)], 
         pch=20, col='blue')
  lines(c(-4,4),c(0,0))
}
for (i in outliers.DE5.down) {
  plot(get(paste0('lfc.DE5.',i,'_outlier.down')), get(paste0('M.DE5.',i,'_outlier.down')) - 
         get(paste0('lfc.DE5.',i,'_outlier.down')), pch=20, col='grey', main='')
  points(get(paste0('lfc.DE5.',i,'_outlier.down'))[which(get(paste0('DEdown.DE5.',i,'_outlier.down')) == 1)], 
         get(paste0('M.DE5.',i,'_outlier.down'))[which(get(paste0('DEdown.DE5.',i,'_outlier.down')) == 1)] - 
           get(paste0('lfc.DE5.',i,'_outlier.down'))[which(get(paste0('DEdown.DE5.',i,'_outlier.down')) == 1)], 
         pch=20, col='blue')
  points(get(paste0('lfc.DE5.',i,'_outlier.down'))[which(get(paste0('DEup.DE5.',i,'_outlier.down')) == 1)], 
         get(paste0('M.DE5.',i,'_outlier.down'))[which(get(paste0('DEup.DE5.',i,'_outlier.down')) == 1)] - 
           get(paste0('lfc.DE5.',i,'_outlier.down'))[which(get(paste0('DEup.DE5.',i,'_outlier.down')) == 1)], 
         pch=20, col='red')
  lines(c(-4,4),c(0,0))
}
par(mfrow=c(4,4), mar=c(1,1,0.5,0.5), mgp=c(3,0,0))
for (i in outliers.DE2.down) {
  plot(get(paste0('lfc.DE2.',i,'_outlier.down')), get(paste0('M.DE2.',i,'_outlier.down')) - 
         get(paste0('lfc.DE2.',i,'_outlier.down')), pch=20, col='grey', main='')
  points(get(paste0('lfc.DE2.',i,'_outlier.down'))[which(get(paste0('DEdown.DE2.',i,'_outlier.down')) == 1)], 
         get(paste0('M.DE2.',i,'_outlier.down'))[which(get(paste0('DEdown.DE2.',i,'_outlier.down')) == 1)] - 
           get(paste0('lfc.DE2.',i,'_outlier.down'))[which(get(paste0('DEdown.DE2.',i,'_outlier.down')) == 1)], 
         pch=20, col='red')
  points(get(paste0('lfc.DE2.',i,'_outlier.down'))[which(get(paste0('DEup.DE2.',i,'_outlier.down')) == 1)], 
         get(paste0('M.DE2.',i,'_outlier.down'))[which(get(paste0('DEup.DE2.',i,'_outlier.down')) == 1)] - 
           get(paste0('lfc.DE2.',i,'_outlier.down'))[which(get(paste0('DEup.DE2.',i,'_outlier.down')) == 1)], 
         pch=20, col='blue')
  lines(c(-4,4),c(0,0))
}
for (i in outliers.DE2.down) {
  plot(get(paste0('lfc.DE2.',i,'_outlier.down')), get(paste0('M.DE2.',i,'_outlier.down')) - 
         get(paste0('lfc.DE2.',i,'_outlier.down')), pch=20, col='grey', main='')
  points(get(paste0('lfc.DE2.',i,'_outlier.down'))[which(get(paste0('DEdown.DE2.',i,'_outlier.down')) == 1)], 
         get(paste0('M.DE2.',i,'_outlier.down'))[which(get(paste0('DEdown.DE2.',i,'_outlier.down')) == 1)] - 
           get(paste0('lfc.DE2.',i,'_outlier.down'))[which(get(paste0('DEdown.DE2.',i,'_outlier.down')) == 1)], 
         pch=20, col='blue')
  points(get(paste0('lfc.DE2.',i,'_outlier.down'))[which(get(paste0('DEup.DE2.',i,'_outlier.down')) == 1)], 
         get(paste0('M.DE2.',i,'_outlier.down'))[which(get(paste0('DEup.DE2.',i,'_outlier.down')) == 1)] - 
           get(paste0('lfc.DE2.',i,'_outlier.down'))[which(get(paste0('DEup.DE2.',i,'_outlier.down')) == 1)], 
         pch=20, col='red')
  lines(c(-4,4),c(0,0))
}
# Doesn't seem to be any more deviation from 0 for red genes than for blue genes. If anything the other 
# way around, but even then not really.

# Plot densities of differences between true and observed LFCs
par(mfrow=c(4,4), mar=c(1,1,0.5,0.5), mgp=c(3,0,0))
for (i in outliers.DE10.up) {
  plot(density(-get(paste0('lfc.DE10.',i,'_outlier.up'))[which(get(paste0('DEdown.DE10.',i,'_outlier.up')) == 1)] + 
                 get(paste0('M.DE10.',i,'_outlier.up'))[which(get(paste0('DEdown.DE10.',i,'_outlier.up')) == 1)]), 
       col='blue', main='')
  lines(density(get(paste0('lfc.DE10.',i,'_outlier.up'))[which(get(paste0('DEup.DE10.',i,'_outlier.up')) == 1)] - 
                  get(paste0('M.DE10.',i,'_outlier.up'))[which(get(paste0('DEup.DE10.',i,'_outlier.up')) == 1)]), 
        col='red')
}
par(mfrow=c(4,4), mar=c(1,1,0.5,0.5), mgp=c(3,0,0))
for (i in outliers.DE10.down) {
  plot(density(-get(paste0('lfc.DE10.',i,'_outlier.down'))[which(get(paste0('DEdown.DE10.',i,'_outlier.down')) == 1)] + 
                 get(paste0('M.DE10.',i,'_outlier.down'))[which(get(paste0('DEdown.DE10.',i,'_outlier.down')) == 1)]), 
       col='red', main='')
  lines(density(get(paste0('lfc.DE10.',i,'_outlier.down'))[which(get(paste0('DEup.DE10.',i,'_outlier.down')) == 1)] - 
                  get(paste0('M.DE10.',i,'_outlier.down'))[which(get(paste0('DEup.DE10.',i,'_outlier.down')) == 1)]), 
        col='blue')
}
par(mfrow=c(4,4), mar=c(1,1,0.5,0.5), mgp=c(3,0,0))
for (i in outliers.DE5.up) {
  plot(density(-get(paste0('lfc.DE5.',i,'_outlier.up'))[which(get(paste0('DEdown.DE5.',i,'_outlier.up')) == 1)] + 
                 get(paste0('M.DE5.',i,'_outlier.up'))[which(get(paste0('DEdown.DE5.',i,'_outlier.up')) == 1)]), 
       col='blue', main='')
  lines(density(get(paste0('lfc.DE5.',i,'_outlier.up'))[which(get(paste0('DEup.DE5.',i,'_outlier.up')) == 1)] - 
                  get(paste0('M.DE5.',i,'_outlier.up'))[which(get(paste0('DEup.DE5.',i,'_outlier.up')) == 1)]), 
        col='red')
}
for (i in outliers.DE5.down) {
  plot(density(-get(paste0('lfc.DE5.',i,'_outlier.down'))[which(get(paste0('DEdown.DE5.',i,'_outlier.down')) == 1)] + 
                 get(paste0('M.DE5.',i,'_outlier.down'))[which(get(paste0('DEdown.DE5.',i,'_outlier.down')) == 1)]), 
       col='red', main='')
  lines(density(get(paste0('lfc.DE5.',i,'_outlier.down'))[which(get(paste0('DEup.DE5.',i,'_outlier.down')) == 1)] - 
                  get(paste0('M.DE5.',i,'_outlier.down'))[which(get(paste0('DEup.DE5.',i,'_outlier.down')) == 1)]), 
        col='blue')
}
par(mfrow=c(4,4), mar=c(1,1,0.5,0.5), mgp=c(3,0,0))
for (i in outliers.DE2.up) {
  plot(density(-get(paste0('lfc.DE2.',i,'_outlier.up'))[which(get(paste0('DEdown.DE2.',i,'_outlier.up')) == 1)] + 
                 get(paste0('M.DE2.',i,'_outlier.up'))[which(get(paste0('DEdown.DE2.',i,'_outlier.up')) == 1)]), 
       col='blue', main='')
  lines(density(get(paste0('lfc.DE2.',i,'_outlier.up'))[which(get(paste0('DEup.DE2.',i,'_outlier.up')) == 1)] - 
                  get(paste0('M.DE2.',i,'_outlier.up'))[which(get(paste0('DEup.DE2.',i,'_outlier.up')) == 1)]), 
        col='red')
}
for (i in outliers.DE2.down) {
  plot(density(-get(paste0('lfc.DE2.',i,'_outlier.down'))[which(get(paste0('DEdown.DE2.',i,'_outlier.down')) == 1)] + 
                 get(paste0('M.DE2.',i,'_outlier.down'))[which(get(paste0('DEdown.DE2.',i,'_outlier.down')) == 1)]), 
       col='red', main='')
  lines(density(get(paste0('lfc.DE2.',i,'_outlier.down'))[which(get(paste0('DEup.DE2.',i,'_outlier.down')) == 1)] - 
                  get(paste0('M.DE2.',i,'_outlier.down'))[which(get(paste0('DEup.DE2.',i,'_outlier.down')) == 1)]), 
        col='blue')
}
# No obvious differences.


# Look at distributions of non-DE genes ####
# If they are skewed in one direction that could cause false positives which 
# could explain poor performance.
par(mfrow=c(4,4), mar=c(1,1,0.5,0.5), mgp=c(3,0,0))
for (i in outliers.DE10.up) {
  plot(density(get(paste0('lfc.DE10.',i,'_outlier.up'))[which(get(paste0('DE.DE10.',i,'_outlier.up')) == 0)]), 
       main='', lwd=2); abline(v=seq(-0.4,0.4,0.1), h=seq(0,3,0.5), col='grey')
}
par(mfrow=c(4,4), mar=c(1,1,0.5,0.5), mgp=c(3,0,0))
for (i in outliers.DE10.down) {
  plot(density(get(paste0('lfc.DE10.',i,'_outlier.down'))[which(get(paste0('DE.DE10.',i,'_outlier.down')) == 0)]), 
       main='', lwd=2); abline(v=seq(-0.4,0.4,0.1), h=seq(0,3,0.5), col='grey')
}
par(mfrow=c(4,4), mar=c(1,1,0.5,0.5), mgp=c(3,0,0))
for (i in outliers.DE5.up) {
  plot(density(get(paste0('lfc.DE5.',i,'_outlier.up'))[which(get(paste0('DE.DE5.',i,'_outlier.up')) == 0)]), 
       main='', lwd=2); abline(v=seq(-0.4,0.4,0.1), h=seq(0,3,0.5), col='grey')
}
for (i in outliers.DE5.down) {
  plot(density(get(paste0('lfc.DE5.',i,'_outlier.down'))[which(get(paste0('DE.DE5.',i,'_outlier.down')) == 0)]), 
       main='', lwd=2); abline(v=seq(-0.4,0.4,0.1), h=seq(0,3,0.5), col='grey')
}
par(mfrow=c(4,4), mar=c(1,1,0.5,0.5), mgp=c(3,0,0))
for (i in outliers.DE2.up) {
  plot(density(get(paste0('lfc.DE2.',i,'_outlier.up'))[which(get(paste0('DE.DE2.',i,'_outlier.up')) == 0)]), 
       main='', lwd=2); abline(v=seq(-0.4,0.4,0.1), h=seq(0,3,0.5), col='grey')
}
for (i in outliers.DE2.down) {
  plot(density(get(paste0('lfc.DE2.',i,'_outlier.down'))[which(get(paste0('DE.DE2.',i,'_outlier.down')) == 0)]), 
       main='', lwd=2); abline(v=seq(-0.4,0.4,0.1), h=seq(0,3,0.5), col='grey')
}
# All look identical and perfectly symmetrical.




# Go back to looking at genes called differently by different methods ####
# Particularly genes called correctly by edgeR and voom, and wrongly by others

# DE10 ####
for (j in 1:50) {
  assign(paste0('results.DE10.',j), readRDS(here('Results/DE compcodeR data results July-Aug 2019',
                                                 paste0('results.DE10.',j,'.all.rds'))))
  assign(paste0('data.DE10.',j), get('data', get(paste0('results.DE10.',j))))
  assign(paste0('DE.DE10.',j), get('DE', get(paste0('results.DE10.',j))))
  assign(paste0('p.edgeR.ql.DE10.',j), get('p.ql.edgeR', get(paste0('results.DE10.',j))))
  assign(paste0('p.edgeR.lr.DE10.',j), get('p.lr.edgeR', get(paste0('results.DE10.',j))))
  assign(paste0('p.edgeR.et.DE10.',j), get('p.et.edgeR', get(paste0('results.DE10.',j))))
  assign(paste0('p.baySeq.DE10.',j), 1 - get('prob.baySeq', get(paste0('results.DE10.',j))))
  assign(paste0('p.MDSeq.zi.DE10.',j), get('p.mean.zi.MDSeq', get(paste0('results.DE10.',j))))
  assign(paste0('p.MDSeq.nozi.DE10.',j), get('p.mean.nozi.MDSeq', get(paste0('results.DE10.',j))))
  assign(paste0('p.DESeq2.if.DE10.',j), get('p.if.DESeq', get(paste0('results.DE10.',j))))
  assign(paste0('p.DESeq2.noif.DE10.',j), get('p.noif.DESeq', get(paste0('results.DE10.',j))))
  assign(paste0('p.voom.DE10.',j), get('p.voom', get(paste0('results.DE10.',j))))
  assign(paste0('p.DSS.notrend.DE10.',j), get('p.notrend.DSS', get(paste0('results.DE10.',j))))
  assign(paste0('p.expHM.untr.DE10.',j), get('p.mean.expHM', get(paste0('results.DE10.',j))))
  assign(paste0('p.expHM.tr.DE10.',j), get('p.lmean.expHM', get(paste0('results.DE10.',j))))
  assign(paste0('p.lnHM.untr.DE10.',j), get('p.mean.lnHM', get(paste0('results.DE10.',j))))
  assign(paste0('p.lnHM.tr.DE10.',j), get('p.lmean.lnHM', get(paste0('results.DE10.',j))))
  assign(paste0('subset.DEup.DE10.',j), 
         which((get('upregulation',slot(get(paste0('data.DE10.',j)), 
                                        'variable.annotations')) == 1 |
                  get('differential.expression', slot(get(paste0('data.DE10.',j)), 
                                                      'variable.annotations')) == 0)))
  assign(paste0('subset.DEdown.DE10.',j), 
         which((get('downregulation',slot(get(paste0('data.DE10.',j)), 
                                          'variable.annotations')) == 1 |
                  get('differential.expression',slot(get(paste0('data.DE10.',j)), 
                                                     'variable.annotations')) == 0)))
  assign(paste0('DE.DE10.',j,'.DEup'), get(paste0('DE.DE10.',j))[get(paste0('subset.DEup.DE10.',j))])
  assign(paste0('DE.DE10.',j,'.DEdown'), get(paste0('DE.DE10.',j))[get(paste0('subset.DEdown.DE10.',j))])
  for (i in c(paste0('edgeR.ql.DE10.',j), paste0('edgeR.lr.DE10.',j),paste0('edgeR.et.DE10.',j),
              paste0('baySeq.DE10.',j), paste0('MDSeq.zi.DE10.',j), paste0('MDSeq.nozi.DE10.',j),
              paste0('DESeq2.if.DE10.',j), paste0('DESeq2.noif.DE10.',j), paste0('voom.DE10.',j),
              paste0('DSS.notrend.DE10.',j), paste0('expHM.untr.DE10.',j), paste0('expHM.tr.DE10.',j),
              paste0('lnHM.untr.DE10.',j), paste0('lnHM.tr.DE10.',j))) {
    assign(paste0('p.',i,'.DEup'), get(paste0('p.',i))[get(paste0('subset.DEup.DE10.',j))])
    assign(paste0('p.',i,'.DEdown'), get(paste0('p.',i))[get(paste0('subset.DEdown.DE10.',j))])
  }
  for (i in c(paste0('edgeR.ql.DE10.',j), paste0('edgeR.lr.DE10.',j),paste0('edgeR.et.DE10.',j),
              paste0('baySeq.DE10.',j), paste0('MDSeq.zi.DE10.',j), paste0('MDSeq.nozi.DE10.',j),
              paste0('DESeq2.if.DE10.',j), paste0('DESeq2.noif.DE10.',j), paste0('voom.DE10.',j),
              paste0('DSS.notrend.DE10.',j), paste0('expHM.untr.DE10.',j), paste0('expHM.tr.DE10.',j),
              paste0('lnHM.untr.DE10.',j), paste0('lnHM.tr.DE10.',j))) {
    assign(paste0('pred.',i,'.DEup'), prediction(1-get(paste0('p.',i,'.DEup')), 
                                                 get(paste0('DE.DE10.',j,'.DEup'))))
    assign(paste0('pred.',i,'.DEdown'), prediction(1-get(paste0('p.',i,'.DEdown')), 
                                                   get(paste0('DE.DE10.',j,'.DEdown'))))
    assign(paste0('pred.',i,'.overall'), prediction(1-get(paste0('p.',i)), 
                                                    get(paste0('DE.DE10.',j))))
    assign(paste0('auc.',i,'.DEup'), performance(get(paste0('pred.',i,'.DEup')), 
                                                 measure='auc')@y.values[[1]])
    assign(paste0('auc.',i,'.DEdown'), performance(get(paste0('pred.',i,'.DEdown')), 
                                                   measure='auc')@y.values[[1]])
    assign(paste0('auc.',i,'.overall'), performance(get(paste0('pred.',i,'.overall')), 
                                                    measure='auc')@y.values[[1]])
    # rm(list=c(paste0('pred.',i,'.DEup'), paste0('pred.',i,'.DEdown')))
  }
  rm(list=c(paste0('results.DE10.',j), paste0('data.DE10.',j)))
}
auc.DE10.up <- data.frame(matrix(nrow=50, ncol=14))
names(auc.DE10.up) <- c('edgeR.ql', 'edgeR.lr', 'edgeR.et', 'baySeq', 'MDSeq.zi', 'MDSeq.nozi', 'DESeq2.if', 
                        'DESeq2.noif', 'voom', 'DSS.notrend', 'expHM.untr', 'expHM.tr', 'lnHM.untr', 'lnHM.tr')
auc.DE10.down <- auc.DE10.up
names(auc.DE10.down) <- names(auc.DE10.up)
auc.DE10.overall <- auc.DE10.up
names(auc.DE10.overall) <- names(auc.DE10.up)
for (i in 1:50) {
  for (j in names(auc.DE10.up)) {
    auc.DE10.up[[j]][i] <- get(paste0('auc.',j,'.DE10.',i,'.DEup'))
    auc.DE10.down[[j]][i] <- get(paste0('auc.',j,'.DE10.',i,'.DEdown'))
    auc.DE10.overall[[j]][i] <- get(paste0('auc.',j,'.DE10.',i,'.overall'))
  }
}
for (i in c('edgeR.ql', 'edgeR.lr', 'edgeR.et', 'baySeq', 'MDSeq.zi', 'MDSeq.nozi', 'DESeq2.if', 
            'DESeq2.noif', 'voom', 'DSS.notrend', 'expHM.untr', 'expHM.tr', 'lnHM.untr', 'lnHM.tr', 
            'DE.', 'DD.', 'subset')) {
  rm(list=ls()[grep(i, ls(), fixed=T)])
}
par(mfrow=c(3,1), mar=c(3,3,1,1))
boxplot(auc.DE10.up, ylim=c(0.7,1)); abline(h=c(0.8,0.85,0.9), col='grey')
boxplot(auc.DE10.down, ylim=c(0.7,1)); abline(h=c(0.8,0.85,0.9), col='grey')
boxplot(auc.DE10.overall, ylim=c(0.7,1)); abline(h=c(0.8,0.85,0.9), col='grey')
round(rbind(colMeans(auc.DE10.up), colMeans(auc.DE10.down), colMeans(auc.DE10.overall)),3)
round(rbind(apply(auc.DE10.up,2,median), apply(auc.DE10.down,2,median), apply(auc.DE10.overall,2,median)),3)
round(rbind(apply(auc.DE10.up,2,max), apply(auc.DE10.down,2,max), apply(auc.DE10.overall,2,max)),3)
round(rbind(apply(auc.DE10.up,2,min), apply(auc.DE10.down,2,min), apply(auc.DE10.overall,2,min)),3)
# More outliers for upregulated for DE10 but not much difference.
outliers.DE10.up <- which(auc.DE10.up$baySeq < 0.9)
outliers.DE10.down <- which(auc.DE10.down$baySeq < 0.9)
outliers.DE10 <- which(auc.DE10.overall$baySeq < 0.88)
auc.DE10.overall[outliers.DE10,]
# edgeR 0.95-0.96; voom 0.93-0.95; baySeq 0.83-0.88; MDSeq 0.84-0.88; DESeq2, DSS 0.85-0.89; HMs 0.84-0.90
# 14 by far the biggest gap so should pay particular attention to that
non.outliers.DE10 <- which(auc.DE10.overall$baySeq>0.943)
# selects top 4 for baySeq
auc.DE10.overall[non.outliers.DE10,]
# DSS, HMs 0.95-0.97; edgeR, MDSeq, DESeq2 0.95-0.96; voom, baySeq 0.94-0.95
# All very close so good controls

# Load data
for (i in outliers.DE10) {
  assign(paste0('results.DE10.',i,'_outlier'), readRDS(here('Results/DE compcodeR data results July-Aug 2019',
                                                               paste0('results.DE10.',i,'.all.rds'))))
  assign(paste0('DE.DE10.',i,'_outlier'), get('DE', get(paste0('results.DE10.',i,'_outlier'))))
  assign(paste0('data.DE10.',i,'_outlier'), get('data', get(paste0('results.DE10.',i,'_outlier'))))
  assign(paste0('lfc.DE10.',i,'_outlier'), get('truelog2foldchanges', slot(get(paste0('data.DE10.',i,'_outlier')),
                                                                              'variable.annotations')))
  assign(paste0('truemeans1.DE10.',i,'_outlier'), get('truemeans.S1', slot(get(paste0('data.DE10.',i,'_outlier')),
                                                                              'variable.annotations')))
  assign(paste0('truemeans2.DE10.',i,'_outlier'), get('truemeans.S2', slot(get(paste0('data.DE10.',i,'_outlier')),
                                                                              'variable.annotations')))
  assign(paste0('truedisps.DE10.',i,'_outlier'), get('truedispersions.S1', slot(get(paste0('data.DE10.',i,'_outlier')),
                                                                                   'variable.annotations')))
  assign(paste0('M.DE10.',i,'_outlier'), get('M.value', slot(get(paste0('data.DE10.',i,'_outlier')),
                                                                'variable.annotations')))
  assign(paste0('A.DE10.',i,'_outlier'), get('A.value', slot(get(paste0('data.DE10.',i,'_outlier')),
                                                                'variable.annotations')))
  assign(paste0('counts.DE10.',i,'_outlier'), slot(get(paste0('data.DE10.',i,'_outlier')),'count.matrix'))
}
for (i in non.outliers.DE10) {
  assign(paste0('results.DE10.',i,'_non.outlier'), readRDS(here('Results/DE compcodeR data results July-Aug 2019',
                                                                 paste0('results.DE10.',i,'.all.rds'))))
  assign(paste0('DE.DE10.',i,'_non.outlier'), get('DE', get(paste0('results.DE10.',i,'_non.outlier'))))
  assign(paste0('data.DE10.',i,'_non.outlier'), get('data', get(paste0('results.DE10.',i,'_non.outlier'))))
  assign(paste0('lfc.DE10.',i,'_non.outlier'), get('truelog2foldchanges', slot(get(paste0('data.DE10.',i,'_non.outlier')),
                                                                                'variable.annotations')))
  assign(paste0('truemeans1.DE10.',i,'_non.outlier'), get('truemeans.S1', slot(get(paste0('data.DE10.',i,'_non.outlier')),
                                                                                'variable.annotations')))
  assign(paste0('truemeans2.DE10.',i,'_non.outlier'), get('truemeans.S2', slot(get(paste0('data.DE10.',i,'_non.outlier')),
                                                                                'variable.annotations')))
  assign(paste0('truedisps.DE10.',i,'_non.outlier'), get('truedispersions.S1', slot(get(paste0('data.DE10.',i,'_non.outlier')),
                                                                                     'variable.annotations')))
  assign(paste0('M.DE10.',i,'_non.outlier'), get('M.value', slot(get(paste0('data.DE10.',i,'_non.outlier')),
                                                                  'variable.annotations')))
  assign(paste0('A.DE10.',i,'_non.outlier'), get('A.value', slot(get(paste0('data.DE10.',i,'_non.outlier')),
                                                                  'variable.annotations')))
  assign(paste0('counts.DE10.',i,'_non.outlier'), slot(get(paste0('data.DE10.',i,'_non.outlier')),'count.matrix'))
}
# Now have, for DE10 and for outliers non-outliers:
# - differentially expressed genes (DE.)
# - log fold changes (lfc.)
# - true means in both groups (truemeans1., truemeans2.)
# - true dispersions (truedisps.)
# - M values (M.)
# - A values (A.)
# - count matrix (counts.)

# Distribution of log fold changes ####
# Plot densities of LFC for each dataset. Non-outliers blue, outliers red.
par(mfrow=c(4,2), mar=c(2,2,1,1))
for (i in outliers.DE10) {
  plot(density(get(paste0('lfc.DE10.',i,'_outlier'))[which(get(paste0('DE.DE10.',i,'_outlier')) == 1)]), 
       main='', xlim=c(-4,4), ylim=c(0,0.4), lwd=2)
  abline(h=seq(0.05,0.35,0.05), v=seq(-3.5,3.5,0.5), col='grey')
}
for (i in non.outliers.DE10) {
  plot(density(get(paste0('lfc.DE10.',i,'_non.outlier'))[which(get(paste0('DE.DE10.',i,'_non.outlier')) == 1)]), 
       main='', xlim=c(-4,4), ylim=c(0,0.4), lwd=2)
  abline(h=seq(0.05,0.35,0.05), v=seq(-3.5,3.5,0.5), col='grey')
}
# No noticeable differences at all between outlier and non-outlier datasets.

# M-A plots (actual LFC v log mean) - look for differences in patterns of DE ####
par(mfrow=c(4,2), mar=c(1,1,0.5,0.5), mgp=c(3,0,0))
for (i in outliers.DE10) {
  plot(get(paste0('A.DE10.',i,'_outlier')), get(paste0('M.DE10.',i,'_outlier')), 
       pch=20, col='grey', ylim=c(-4.5,4.5), xlim=c(1,17))
  points(get(paste0('A.DE10.',i,'_outlier'))[which(get(paste0('DE.DE10.',i,'_outlier')) == 1)], 
         get(paste0('M.DE10.',i,'_outlier'))[which(get(paste0('DE.DE10.',i,'_outlier')) == 1)], 
         pch=20, col='blue')
  abline(h=seq(-8,8,2), col='grey')
}
for (i in non.outliers.DE10) {
  plot(get(paste0('A.DE10.',i,'_non.outlier')), get(paste0('M.DE10.',i,'_non.outlier')), 
       pch=20, col='grey', ylim=c(-4.5,4.5), xlim=c(1,17))
  points(get(paste0('A.DE10.',i,'_non.outlier'))[which(get(paste0('DE.DE10.',i,'_non.outlier')) == 1)], 
         get(paste0('M.DE10.',i,'_non.outlier'))[which(get(paste0('DE.DE10.',i,'_non.outlier')) == 1)], 
         pch=20, col='red')
  abline(h=seq(-8,8,2), col='grey')
}
# No obvious differences. Would expect maybe to see more non-DE genes with high observed LFC or more DE 
# genes with low observed LFC in outlier datasets, but that doesn't seem to be the case.
for (i in outliers.DE10) {
  print(c(range(get(paste0('M.DE10.',i,'_outlier'))[which(get(paste0('DE.DE10.',i,'_outlier')) == 1)]), 
        range(get(paste0('M.DE10.',i,'_outlier'))[which(get(paste0('DE.DE10.',i,'_outlier')) == 0)])))
}
for (i in non.outliers.DE10) {
  print(c(range(get(paste0('M.DE10.',i,'_non.outlier'))[which(get(paste0('DE.DE10.',i,'_non.outlier')) == 1)]), 
          range(get(paste0('M.DE10.',i,'_non.outlier'))[which(get(paste0('DE.DE10.',i,'_non.outlier')) == 0)])))
}
# Actually looks like a wider range of observed LFC for DE genes in outlier datasets compared to non-outlier
for (i in outliers.DE10) {
  print(c(mean(abs(get(paste0('M.DE10.',i,'_outlier'))[which(get(paste0('DE.DE10.',i,'_outlier')) == 1)]) > 2), 
          mean(abs(get(paste0('M.DE10.',i,'_outlier'))[which(get(paste0('DE.DE10.',i,'_outlier')) == 0)]) > 2)))
}
for (i in non.outliers.DE10) {
  print(c(mean(abs(get(paste0('M.DE10.',i,'_non.outlier'))[which(get(paste0('DE.DE10.',i,'_non.outlier')) == 1)]) > 2), 
          mean(abs(get(paste0('M.DE10.',i,'_non.outlier'))[which(get(paste0('DE.DE10.',i,'_non.outlier')) == 0)]) > 2)))
}
# For the biggest outlier (14; index 2 in outliers.DE10), clearly lower proportion of DE genes with extreme 
# observed LFCs, but not for the other outliers, and also a lower proportion of non-DE genes with extreme 
# observed LFCs. If look at LFC > 1 rather than 2, no obvious differences between outliers and non-outliers.

# Look at differences between true and observed LFCs ####
# May expect outlier datasets to have less correlation between true and observed LFCs

# Plot (observed - true) v true LFC
par(mfrow=c(4,2), mar=c(1,1,0.5,0.5), mgp=c(3,0,0))
for (i in outliers.DE10) {
  plot(get(paste0('lfc.DE10.',i,'_outlier')), get(paste0('M.DE10.',i,'_outlier')) - 
         get(paste0('lfc.DE10.',i,'_outlier')), pch=20, main='', ylim=c(-6,6), xlim=c(-3.5,3.5))
  abline(h=seq(-5,5,1), col='grey')
}
for (i in non.outliers.DE10) {
  plot(get(paste0('lfc.DE10.',i,'_non.outlier')), get(paste0('M.DE10.',i,'_non.outlier')) - 
         get(paste0('lfc.DE10.',i,'_non.outlier')), pch=20, main='', ylim=c(-6,6), xlim=c(-3.5,3.5))
  abline(h=seq(-5,5,1), col='grey')
}
# Doesn't seem to be any more deviation from 0 for outliers than non-outliers.

# Plot densities of differences between true and observed LFCs
par(mfrow=c(4,2), mar=c(1,1,0.5,0.5), mgp=c(3,0,0))
for (i in outliers.DE10) {
  plot(density(get(paste0('lfc.DE10.',i,'_outlier'))[which(get(paste0('DE.DE10.',i,'_outlier')) == 1)] - 
                 get(paste0('M.DE10.',i,'_outlier'))[which(get(paste0('DE.DE10.',i,'_outlier')) == 1)]), 
       main='', xlim=c(-4,4)); abline(v=seq(-3,3,1), col='grey')
}
for (i in non.outliers.DE10) {
  plot(density(get(paste0('lfc.DE10.',i,'_non.outlier'))[which(get(paste0('DE.DE10.',i,'_non.outlier')) == 1)] - 
                 get(paste0('M.DE10.',i,'_non.outlier'))[which(get(paste0('DE.DE10.',i,'_non.outlier')) == 1)]), 
       main='', xlim=c(-4,4)); abline(v=seq(-3,3,1), col='grey')
}
# No obvious differences.

# Look at distributions of non-DE genes ####
# If they are skewed in one direction that could cause false positives which 
# could explain poor performance.
par(mfrow=c(4,2), mar=c(1,1,0.5,0.5), mgp=c(3,0,0))
for (i in outliers.DE10) {
  plot(density(get(paste0('lfc.DE10.',i,'_outlier'))[which(get(paste0('DE.DE10.',i,'_outlier')) == 0)]), 
       main='', lwd=2); abline(v=seq(-0.4,0.4,0.1), h=seq(0,3,0.5), col='grey')
}
for (i in outliers.DE10) {
  plot(density(get(paste0('lfc.DE10.',i,'_outlier'))[which(get(paste0('DE.DE10.',i,'_outlier')) == 0)]), 
       main='', lwd=2); abline(v=seq(-0.4,0.4,0.1), h=seq(0,3,0.5), col='grey')
}
# All look identical and perfectly symmetrical.






# Look at genes called differently by good and bad methods ####
# Look at DE10.14 as that is the most extreme outlier
results.DE10.14 <- readRDS(here('Results/DE compcodeR data results July-Aug 2019', 'results.DE10.14.all.rds'))
data.DE10.14 <- results.DE10.14$data
norm.DE10.14 <- t(t(data.DE10.14@count.matrix) / calcNormFactors(data.DE10.14@count.matrix))
DE.DE10.14 <- results.DE10.14$DE

# AUCs
p.edgeR.ql.DE10.14 <- results.DE10.14$p.ql.edgeR
p.edgeR.lr.DE10.14 <- results.DE10.14$p.lr.edgeR
p.edgeR.et.DE10.14 <- results.DE10.14$p.et.edgeR
p.baySeq.DE10.14 <- 1 - results.DE10.14$prob.baySeq
p.MDSeq.zi.DE10.14 <- results.DE10.14$p.mean.zi.MDSeq
p.MDSeq.nozi.DE10.14 <- results.DE10.14$p.mean.nozi.MDSeq
p.DESeq2.if.DE10.14 <- results.DE10.14$p.if.DESeq
p.DESeq2.noif.DE10.14 <- results.DE10.14$p.noif.DESeq
p.voom.DE10.14 <- results.DE10.14$p.voom
p.DSS.notrend.DE10.14 <- results.DE10.14$p.notrend.DSS
p.expHM.untr.DE10.14 <- results.DE10.14$p.mean.expHM
p.expHM.tr.DE10.14 <- results.DE10.14$p.lmean.expHM
p.lnHM.untr.DE10.14 <- results.DE10.14$p.mean.lnHM
p.lnHM.tr.DE10.14 <- results.DE10.14$p.lmean.lnHM
pred.edgeR.ql.DE10.14 <- prediction(1-p.edgeR.ql.DE10.14, DE.DE10.14)
pred.edgeR.lr.DE10.14 <- prediction(1-p.edgeR.lr.DE10.14, DE.DE10.14)
pred.edgeR.et.DE10.14 <- prediction(1-p.edgeR.et.DE10.14, DE.DE10.14)
pred.baySeq.DE10.14 <- prediction(1-p.baySeq.DE10.14, DE.DE10.14)
pred.MDSeq.zi.DE10.14 <- prediction(1-p.MDSeq.zi.DE10.14, DE.DE10.14)
pred.MDSeq.nozi.DE10.14 <- prediction(1-p.MDSeq.nozi.DE10.14, DE.DE10.14)
pred.DESeq2.if.DE10.14 <- prediction(1-p.DESeq2.if.DE10.14, DE.DE10.14)
pred.DESeq2.noif.DE10.14 <- prediction(1-p.DESeq2.noif.DE10.14, DE.DE10.14)
pred.voom.DE10.14 <- prediction(1-p.voom.DE10.14, DE.DE10.14)
pred.DSS.notrend.DE10.14 <- prediction(1-p.DSS.notrend.DE10.14, DE.DE10.14)
pred.expHM.untr.DE10.14 <- prediction(1-p.expHM.untr.DE10.14, DE.DE10.14)
pred.expHM.tr.DE10.14 <- prediction(1-p.expHM.tr.DE10.14, DE.DE10.14)
pred.lnHM.untr.DE10.14 <- prediction(1-p.lnHM.untr.DE10.14, DE.DE10.14)
pred.lnHM.tr.DE10.14 <- prediction(1-p.lnHM.tr.DE10.14, DE.DE10.14)
auc.edgeR.ql.DE10.14 <- performance(pred.edgeR.ql.DE10.14, measure='auc')@y.values[[1]]
auc.edgeR.lr.DE10.14 <- performance(pred.edgeR.lr.DE10.14, measure='auc')@y.values[[1]]
auc.edgeR.et.DE10.14 <- performance(pred.edgeR.et.DE10.14, measure='auc')@y.values[[1]]
auc.baySeq.DE10.14 <- performance(pred.baySeq.DE10.14, measure='auc')@y.values[[1]]
auc.MDSeq.zi.DE10.14 <- performance(pred.MDSeq.zi.DE10.14, measure='auc')@y.values[[1]]
auc.MDSeq.nozi.DE10.14 <- performance(pred.MDSeq.nozi.DE10.14, measure='auc')@y.values[[1]]
auc.DESeq2.if.DE10.14 <- performance(pred.DESeq2.if.DE10.14, measure='auc')@y.values[[1]]
auc.DESeq2.noif.DE10.14 <- performance(pred.DESeq2.noif.DE10.14, measure='auc')@y.values[[1]]
auc.voom.DE10.14 <- performance(pred.voom.DE10.14, measure='auc')@y.values[[1]]
auc.DSS.notrend.DE10.14 <- performance(pred.DSS.notrend.DE10.14, measure='auc')@y.values[[1]]
auc.expHM.untr.DE10.14 <- performance(pred.expHM.untr.DE10.14, measure='auc')@y.values[[1]]
auc.expHM.tr.DE10.14 <- performance(pred.expHM.tr.DE10.14, measure='auc')@y.values[[1]]
auc.lnHM.untr.DE10.14 <- performance(pred.lnHM.untr.DE10.14, measure='auc')@y.values[[1]]
auc.lnHM.tr.DE10.14 <- performance(pred.lnHM.tr.DE10.14, measure='auc')@y.values[[1]]
rm(list=ls()[which(grepl('pred.',ls()))])
c(auc.edgeR.ql.DE10.14, auc.edgeR.lr.DE10.14, auc.edgeR.et.DE10.14, auc.baySeq.DE10.14, 
  auc.MDSeq.zi.DE10.14, auc.MDSeq.nozi.DE10.14, auc.DESeq2.if.DE10.14, auc.DESeq2.noif.DE10.14, 
  auc.voom.DE10.14, auc.DSS.notrend.DE10.14, auc.expHM.untr.DE10.14, auc.expHM.tr.DE10.14, 
  auc.lnHM.untr.DE10.14, auc.lnHM.tr.DE10.14)
# edgeR 0.95-0.96; voom 0.94; HMs 0.84-0.86; DESeq2, DSS 0.85; MDSeq 0.84; baySeq 0.83

# Confusion matrices
q.edgeR.ql.DE10.14 <- results.DE10.14$q.ql.edgeR
q.edgeR.lr.DE10.14 <- results.DE10.14$q.lr.edgeR
q.edgeR.et.DE10.14 <- results.DE10.14$q.et.edgeR
q.baySeq.DE10.14 <- results.DE10.14$q.baySeq
q.MDSeq.zi.DE10.14 <- results.DE10.14$q.mean.zi.MDSeq
q.MDSeq.nozi.DE10.14 <- results.DE10.14$q.mean.nozi.MDSeq
q.DESeq2.if.DE10.14 <- results.DE10.14$q.if.DESeq
q.DESeq2.noif.DE10.14 <- results.DE10.14$q.noif.DESeq
q.voom.DE10.14 <- results.DE10.14$q.voom
q.DSS.notrend.DE10.14 <- results.DE10.14$q.notrend.DSS
q.expHM.untr.DE10.14 <- as.vector(p.adjust(results.DE10.14$p.mean.expHM, method='BH'))
q.expHM.tr.DE10.14 <- as.vector(p.adjust(results.DE10.14$p.lmean.expHM, method='BH'))
q.lnHM.untr.DE10.14 <- as.vector(p.adjust(results.DE10.14$p.mean.lnHM, method='BH'))
q.lnHM.tr.DE10.14 <- as.vector(p.adjust(results.DE10.14$p.lmean.lnHM, method='BH'))
c(mean(q.edgeR.ql.DE10.14), mean(q.edgeR.lr.DE10.14), mean(q.edgeR.et.DE10.14), mean(q.baySeq.DE10.14), 
  mean(q.MDSeq.zi.DE10.14, na.rm=T), mean(q.MDSeq.nozi.DE10.14, na.rm=T), mean(q.DESeq2.if.DE10.14), 
  mean(q.DESeq2.noif.DE10.14),mean(q.voom.DE10.14), mean(q.DSS.notrend.DE10.14), mean(q.expHM.untr.DE10.14), 
  mean(q.expHM.tr.DE10.14), mean(q.lnHM.untr.DE10.14), mean(q.lnHM.tr.DE10.14))
# As for DEDD10.45, only MDSeq has removed genes.
c(sum(is.na(q.MDSeq.zi.DE10.14)), sum(is.na(q.MDSeq.nozi.DE10.14)))
# 10 genes removed for for zi, 20 for nozi; interesting that nozi removed exactly double again
# edgeR mean q-value around 0.89; voom 0.86; MDSeq 0.84; DSS 0.80; baySeq 0.78; HMs 0.47; DESeq2 0.43
# Doesn't seem to correlate with performance; probably says something about distribution of values 
# for non-DE genes, but looks like says nothing about why some methods perform badly for this

for (i in c('edgeR.ql.DE10.14', 'edgeR.lr.DE10.14', 'edgeR.et.DE10.14', 'baySeq.DE10.14', 
            'MDSeq.zi.DE10.14', 'MDSeq.nozi.DE10.14', 'DESeq2.if.DE10.14', 'DESeq2.noif.DE10.14', 
            'voom.DE10.14', 'DSS.notrend.DE10.14', 'expHM.untr.DE10.14', 'expHM.tr.DE10.14', 
            'lnHM.untr.DE10.14', 'lnHM.tr.DE10.14')) {
  assign(paste0('call.',i), get(paste0('q.',i)) < 0.05)
}
c(mean(call.edgeR.ql.DE10.14), mean(call.edgeR.lr.DE10.14), mean(call.edgeR.et.DE10.14), 
  mean(call.baySeq.DE10.14),   mean(call.MDSeq.zi.DE10.14, na.rm=T), mean(call.MDSeq.nozi.DE10.14, na.rm=T), 
  mean(call.DESeq2.if.DE10.14), mean(call.DESeq2.noif.DE10.14),mean(call.voom.DE10.14), 
  mean(call.DSS.notrend.DE10.14), mean(call.expHM.untr.DE10.14), mean(call.expHM.tr.DE10.14), 
  mean(call.lnHM.untr.DE10.14), mean(call.lnHM.tr.DE10.14))
# DESeq2 calls 0.08 genes DE; MDSeq 0.06; HMs, edgeR 0.03-0.04; voom, baySeq 0.03
# Similar patter to DEDD2.45. Proportion of genes called DE doesn't seem to correlate with performance.

for (i in c('edgeR.ql.DE10.14', 'edgeR.lr.DE10.14', 'edgeR.et.DE10.14', 'baySeq.DE10.14', 
            'MDSeq.zi.DE10.14', 'MDSeq.nozi.DE10.14', 'DESeq2.if.DE10.14', 'DESeq2.noif.DE10.14', 
            'voom.DE10.14', 'DSS.notrend.DE10.14', 'expHM.untr.DE10.14', 'expHM.tr.DE10.14', 
            'lnHM.untr.DE10.14', 'lnHM.tr.DE10.14')) {
  assign(paste0('conf.',i), confusionMatrix(factor(as.numeric(get(paste0('call.',i)))), factor(DE.DE10.14), 
                                            positive='1'))
}
c(unname(conf.edgeR.ql.DE10.14$overall['Accuracy']), unname(conf.edgeR.lr.DE10.14$overall['Accuracy']), 
  unname(conf.edgeR.et.DE10.14$overall['Accuracy']), unname(conf.baySeq.DE10.14$overall['Accuracy']), 
  unname(conf.MDSeq.zi.DE10.14$overall['Accuracy']), unname(conf.MDSeq.nozi.DE10.14$overall['Accuracy']), 
  unname(conf.DESeq2.if.DE10.14$overall['Accuracy']), unname(conf.DESeq2.noif.DE10.14$overall['Accuracy']), 
  unname(conf.voom.DE10.14$overall['Accuracy']), unname(conf.DSS.notrend.DE10.14$overall['Accuracy']), 
  unname(conf.expHM.untr.DE10.14$overall['Accuracy']), unname(conf.expHM.tr.DE10.14$overall['Accuracy']), 
  unname(conf.lnHM.untr.DE10.14$overall['Accuracy']), unname(conf.lnHM.tr.DE10.14$overall['Accuracy']))
# No real differences in overall accuracy. edgeR, voom, DSS 0.98; baySeq 0.97; HMs 0.96; MDSeq 0.95; DESeq2 0.93

c(unname(conf.edgeR.ql.DE10.14$byClass['Sensitivity']), unname(conf.edgeR.lr.DE10.14$byClass['Sensitivity']), 
  unname(conf.edgeR.et.DE10.14$byClass['Sensitivity']), unname(conf.baySeq.DE10.14$byClass['Sensitivity']), 
  unname(conf.MDSeq.zi.DE10.14$byClass['Sensitivity']), unname(conf.MDSeq.nozi.DE10.14$byClass['Sensitivity']), 
  unname(conf.DESeq2.if.DE10.14$byClass['Sensitivity']), unname(conf.DESeq2.noif.DE10.14$byClass['Sensitivity']), 
  unname(conf.voom.DE10.14$byClass['Sensitivity']), unname(conf.DSS.notrend.DE10.14$byClass['Sensitivity']), 
  unname(conf.expHM.untr.DE10.14$byClass['Sensitivity']), unname(conf.expHM.tr.DE10.14$byClass['Sensitivity']), 
  unname(conf.lnHM.untr.DE10.14$byClass['Sensitivity']), unname(conf.lnHM.tr.DE10.14$byClass['Sensitivity']))
# TPR: edgeR 0.61-0.66; DESeq2 0.64; voom, DSS 0.60; MDSeq 0.57; HMs 0.48-0.56; baySeq 0.46.
# Same pattern as for DEDD10.45. Doesn't really correlate with AUC.

c(unname(1-conf.edgeR.ql.DE10.14$byClass['Precision']), unname(1-conf.edgeR.lr.DE10.14$byClass['Precision']), 
  unname(1-conf.edgeR.et.DE10.14$byClass['Precision']), unname(1-conf.baySeq.DE10.14$byClass['Precision']), 
  unname(1-conf.MDSeq.zi.DE10.14$byClass['Precision']), unname(1-conf.MDSeq.nozi.DE10.14$byClass['Precision']), 
  unname(1-conf.DESeq2.if.DE10.14$byClass['Precision']), unname(1-conf.DESeq2.noif.DE10.14$byClass['Precision']), 
  unname(1-conf.voom.DE10.14$byClass['Precision']), unname(1-conf.DSS.notrend.DE10.14$byClass['Precision']), 
  unname(1-conf.expHM.untr.DE10.14$byClass['Precision']), unname(1-conf.expHM.tr.DE10.14$byClass['Precision']), 
  unname(1-conf.lnHM.untr.DE10.14$byClass['Precision']), unname(1-conf.lnHM.tr.DE10.14$byClass['Precision']))
# Good FDR control only for voom, edgeR and DSS (0.03-0.09); baySeq, HMs 0.2-0.3; MDSeq, DESeq2 > 0.5.
# Apart from voom, edgeR and DSS, much worse than for DEDD10.45.

c(unname(conf.edgeR.ql.DE10.14$byClass['Specificity']), unname(conf.edgeR.lr.DE10.14$byClass['Specificity']), 
  unname(conf.edgeR.et.DE10.14$byClass['Specificity']), unname(conf.baySeq.DE10.14$byClass['Specificity']), 
  unname(conf.MDSeq.zi.DE10.14$byClass['Specificity']), unname(conf.MDSeq.nozi.DE10.14$byClass['Specificity']), 
  unname(conf.DESeq2.if.DE10.14$byClass['Specificity']), unname(conf.DESeq2.noif.DE10.14$byClass['Specificity']), 
  unname(conf.voom.DE10.14$byClass['Specificity']), unname(conf.DSS.notrend.DE10.14$byClass['Specificity']), 
  unname(conf.expHM.untr.DE10.14$byClass['Specificity']), unname(conf.expHM.tr.DE10.14$byClass['Specificity']), 
  unname(conf.lnHM.untr.DE10.14$byClass['Specificity']), unname(conf.lnHM.tr.DE10.14$byClass['Specificity']))
# Specificity very high (i.e. FPR very low) for all. edgeR, voom, DSS, baySeq > 0.99; HMs 0.99; MDSeq 0.97; 
# DESeq2 0.94.

c(unname(conf.edgeR.ql.DE10.14$byClass['F1']), unname(conf.edgeR.lr.DE10.14$byClass['F1']), 
  unname(conf.edgeR.et.DE10.14$byClass['F1']), unname(conf.baySeq.DE10.14$byClass['F1']), 
  unname(conf.MDSeq.zi.DE10.14$byClass['F1']), unname(conf.MDSeq.nozi.DE10.14$byClass['F1']), 
  unname(conf.DESeq2.if.DE10.14$byClass['F1']), unname(conf.DESeq2.noif.DE10.14$byClass['F1']), 
  unname(conf.voom.DE10.14$byClass['F1']), unname(conf.DSS.notrend.DE10.14$byClass['F1']), 
  unname(conf.expHM.untr.DE10.14$byClass['F1']), unname(conf.expHM.tr.DE10.14$byClass['F1']), 
  unname(conf.lnHM.untr.DE10.14$byClass['F1']), unname(conf.lnHM.tr.DE10.14$byClass['F1']))
# F1: edgeR 0.75-0.77; voom 0.73; DSS 0.72; HMs 0.57-0.62; baySeq 0.58; MDSeq 0.52; DESeq2 0.47.

c(unname(conf.edgeR.ql.DE10.14$byClass['Balanced Accuracy']), unname(conf.edgeR.lr.DE10.14$byClass['Balanced Accuracy']), 
  unname(conf.edgeR.et.DE10.14$byClass['Balanced Accuracy']), unname(conf.baySeq.DE10.14$byClass['Balanced Accuracy']), 
  unname(conf.MDSeq.zi.DE10.14$byClass['Balanced Accuracy']), unname(conf.MDSeq.nozi.DE10.14$byClass['Balanced Accuracy']), 
  unname(conf.DESeq2.if.DE10.14$byClass['Balanced Accuracy']), unname(conf.DESeq2.noif.DE10.14$byClass['Balanced Accuracy']), 
  unname(conf.voom.DE10.14$byClass['Balanced Accuracy']), unname(conf.DSS.notrend.DE10.14$byClass['Balanced Accuracy']), 
  unname(conf.expHM.untr.DE10.14$byClass['Balanced Accuracy']), unname(conf.expHM.tr.DE10.14$byClass['Balanced Accuracy']), 
  unname(conf.lnHM.untr.DE10.14$byClass['Balanced Accuracy']), unname(conf.lnHM.tr.DE10.14$byClass['Balanced Accuracy']))
# Balanced accuracy: edgeR 0.81-0.83; voom, DSS 0.80; DESeq2 0.79; MDSeq 0.77; HMs 0.74-0.777; baySeq 0.73.

calls <- data.frame(edgeR.ql=call.edgeR.ql.DE10.14, edgeR.lr=call.edgeR.lr.DE10.14, 
                    edgeR.et=call.edgeR.et.DE10.14, baySeq=call.baySeq.DE10.14, MDSeq.zi=call.MDSeq.zi.DE10.14, 
                    MDSeq.nozi=call.MDSeq.nozi.DE10.14, DESeq2.if=call.DESeq2.if.DE10.14, 
                    DESeq2.noif=call.DESeq2.noif.DE10.14, voom=call.voom.DE10.14, 
                    DSS.notrend=call.DSS.notrend.DE10.14, expHM.untr=call.expHM.untr.DE10.14, 
                    expHM.tr=call.expHM.tr.DE10.14, lnHM.untr=call.lnHM.untr.DE10.14, 
                    lnHM.tr=call.lnHM.tr.DE10.14)
round(cor(calls, use='complete.obs'), 2)
# edgeR and voom most highly correlated. different versions of MDSeq and DESeq2 identical.

differing.results <- which((call.baySeq.DE10.14==call.MDSeq.zi.DE10.14) & (call.baySeq.DE10.14==call.MDSeq.nozi.DE10.14) &
                             (call.baySeq.DE10.14==call.DESeq2.if.DE10.14) & (call.baySeq.DE10.14==call.DESeq2.noif.DE10.14) & 
                             (call.baySeq.DE10.14==call.DSS.notrend.DE10.14) & (call.baySeq.DE10.14==call.expHM.untr.DE10.14) & 
                             (call.baySeq.DE10.14==call.expHM.tr.DE10.14) & (call.baySeq.DE10.14==call.lnHM.untr.DE10.14) & 
                             (call.baySeq.DE10.14==call.lnHM.tr.DE10.14) & 
                             (call.baySeq.DE10.14!=call.edgeR.ql.DE10.14) & (call.baySeq.DE10.14!=call.edgeR.lr.DE10.14) & 
                             (call.baySeq.DE10.14!=call.edgeR.et.DE10.14) & (call.baySeq.DE10.14!=call.voom.DE10.14))
sum(call.baySeq.DE10.14[differing.results]==DE.DE10.14[differing.results])
sum(call.baySeq.DE10.14[differing.results]!=DE.DE10.14[differing.results])
# 2 instances where others all right and edgeR and voom wrong, 15 other way around.
table(DE.DE10.14[differing.results])
# Overall 5 non-DE and 12 DE genes have differing results between the groups of methods.
differing.results.tab <- rbind(DE.DE10.14[differing.results], call.baySeq.DE10.14[differing.results], 
                               call.edgeR.ql.DE10.14[differing.results])
rownames(differing.results.tab) <- c('DE','baySeq','edgeR')
colnames(differing.results.tab) <- 1:17
differing.results.tab
#  1    DE, baySeq right
# 11    DE,  edgeR right
# Total 12/811 instances = 0.015
#  1 no DE, baySeq right
#  4 no DE,  edgeR right
# Total 5/15547 = 0.0003
# Most likely to disagree for true DE genes.

# Look at genes correctly called DE by edgeR and vooom but not by others.
# Any obvious similarities? Zeros? Outliers? Sample means different from true means?
DE.edgeR.not.others <- differing.results[which(differing.results %in% which(DE.DE10.14==1 & call.edgeR.lr.DE10.14==1))]
data.DE10.14@variable.annotations$truedispersions.S1[DE.edgeR.not.others]
rbind(data.DE10.14@variable.annotations$truemeans.S1[DE.edgeR.not.others], 
      rowMeans(data.DE10.14@count.matrix[DE.edgeR.not.others,1:10]), 
      rowMeans(norm.DE10.14[DE.edgeR.not.others,1:10]))
rbind(data.DE10.14@variable.annotations$truemeans.S2[DE.edgeR.not.others], 
      rowMeans(data.DE10.14@count.matrix[DE.edgeR.not.others,11:20]), 
      rowMeans(norm.DE10.14[DE.edgeR.not.others,11:20]))
rbind(data.DE10.14@variable.annotations$truelog2foldchanges[DE.edgeR.not.others], 
      data.DE10.14@variable.annotations$M.value[DE.edgeR.not.others])
mean(abs(data.DE10.14@variable.annotations$truelog2foldchanges[which(DE.DE10.14==1)]))
# True LFCs for these genes all smaller than average, but not extremely small; roughly half are in 
# bottom 10% of LFCs.
plot(data.DE10.14@variable.annotations$truemeans.S1, rowMeans(data.DE10.14@count.matrix[,1:10]), 
     pch=20, col='grey', xlim=c(0,3000), ylim=c(0,3000))
points(data.DE10.14@variable.annotations$truemeans.S2, rowMeans(data.DE10.14@count.matrix[,11:20]), 
       pch=20, col='grey')
points(data.DE10.14@variable.annotations$truemeans.S1[which(DE.DE10.14==1)], 
       rowMeans(data.DE10.14@count.matrix[which(DE.DE10.14==1),1:10]), pch=20, col='red')
points(data.DE10.14@variable.annotations$truemeans.S2[which(DE.DE10.14==1)], 
       rowMeans(data.DE10.14@count.matrix[which(DE.DE10.14==1),11:20]), pch=20, col='red')
points(data.DE10.14@variable.annotations$truemeans.S1[DE.edgeR.not.others], 
     rowMeans(data.DE10.14@count.matrix[DE.edgeR.not.others,1:10]), pch=20, col='blue')
points(data.DE10.14@variable.annotations$truemeans.S2[DE.edgeR.not.others], 
     rowMeans(data.DE10.14@count.matrix[DE.edgeR.not.others,11:20]), pch=20, col='blue')
lines(c(0,3000), c(0,3000))
# Genes called correctly by edgeR and voom and wrongly by others don't have any noticeable 
# differences from other DE genes in terms of difference between true and samples means.

# Seems very strange that sample means consistently overestimate true means. Seems to be 
# caused by simulation of sequencing depth factors, which are multuplicative, not additive. 
# Following a comment in compcodeR code at 
# https://github.com/csoneson/compcodeR/blob/master/R/generateSyntheticData.R, can obtain 
# means used in NB simulation for each gene and sample by multiplying each 'true mean' 
# (@variable.annotations$truemeans.S1 or S2) by seq depth (default 1e7; 
# @info.parameters$seqdepth) and the depth factor for that sample 
# (@sample.annotations$depth.factor), and dividing by the sum of the 'true means'. Doing 
# this in a test dataset gave 'actual' means that matched the range of the sample means.
# Thought that using normalised counts rather than raw counts might remove this effect, but 
# it doesn't. Plots look nearly identical:
plot(data.DE10.14@variable.annotations$truemeans.S1, rowMeans(norm.DE10.14[,1:10]), 
     pch=20, col='grey', xlim=c(0,3000), ylim=c(0,3000))
points(data.DE10.14@variable.annotations$truemeans.S2, rowMeans(norm.DE10.14[,11:20]), 
       pch=20, col='grey')
points(data.DE10.14@variable.annotations$truemeans.S1[which(DE.DE10.14==1)], 
       rowMeans(norm.DE10.14[which(DE.DE10.14==1),1:10]), pch=20, col='red')
points(data.DE10.14@variable.annotations$truemeans.S2[which(DE.DE10.14==1)], 
       rowMeans(norm.DE10.14[which(DE.DE10.14==1),11:20]), pch=20, col='red')
points(data.DE10.14@variable.annotations$truemeans.S1[DE.edgeR.not.others], 
       rowMeans(norm.DE10.14[DE.edgeR.not.others,1:10]), pch=20, col='blue')
points(data.DE10.14@variable.annotations$truemeans.S2[DE.edgeR.not.others], 
       rowMeans(norm.DE10.14[DE.edgeR.not.others,11:20]), pch=20, col='blue')
lines(c(0,3000), c(0,3000))
# Looking at 'actual' rather than 'true' means might be useful in terms of seeing if there 
# is anything unusual about outlier datasets, but suspect not, because methods work on 
# the data that is simulated, not on the underlying model.
DE10.14.actualmeans.S1 <- matrix(nrow=nrow(data.DE10.14@count.matrix), ncol=ncol(data.DE10.14@count.matrix))
DE10.14.actualmeans.S2 <- matrix(nrow=nrow(data.DE10.14@count.matrix), ncol=ncol(data.DE10.14@count.matrix))
for (i in 1:ncol(data.DE10.14@count.matrix)) {
  for (j in 1:nrow(data.DE10.14@count.matrix)) {
    DE10.14.actualmeans.S1[j,i] <- data.DE10.14@variable.annotations$truemeans.S1[j] * 
      data.DE10.14@info.parameters$seqdepth * 
      data.DE10.14@sample.annotations$depth.factor[i] / 
      sum(data.DE10.14@variable.annotations$truemeans.S1)
    DE10.14.actualmeans.S2[j,i] <- data.DE10.14@variable.annotations$truemeans.S2[j] * 
      data.DE10.14@info.parameters$seqdepth * 
      data.DE10.14@sample.annotations$depth.factor[i] / 
      sum(data.DE10.14@variable.annotations$truemeans.S2)
  }
}
plot(rowMeans(DE10.14.actualmeans.S1), rowMeans(data.DE10.14@count.matrix[,1:10]), 
     pch=20, col='grey', xlim=c(0,3000), ylim=c(0,3000))
points(rowMeans(DE10.14.actualmeans.S2), rowMeans(data.DE10.14@count.matrix[,11:20]), 
       pch=20, col='grey')
points(rowMeans(DE10.14.actualmeans.S1)[which(DE.DE10.14==1)], 
       rowMeans(data.DE10.14@count.matrix[which(DE.DE10.14==1),1:10]), pch=20, col='red')
points(rowMeans(DE10.14.actualmeans.S2)[which(DE.DE10.14==1)], 
       rowMeans(data.DE10.14@count.matrix[which(DE.DE10.14==1),11:20]), pch=20, col='red')
points(rowMeans(DE10.14.actualmeans.S1)[DE.edgeR.not.others], 
       rowMeans(data.DE10.14@count.matrix[DE.edgeR.not.others,1:10]), pch=20, col='blue')
points(rowMeans(DE10.14.actualmeans.S2)[DE.edgeR.not.others], 
       rowMeans(data.DE10.14@count.matrix[DE.edgeR.not.others,11:20]), pch=20, col='blue')
lines(c(0,3000), c(0,3000))
plot(rowMeans(DE10.14.actualmeans.S1), rowMeans(norm.DE10.14[,1:10]), 
     pch=20, col='grey', xlim=c(0,3000), ylim=c(0,3000))
points(rowMeans(DE10.14.actualmeans.S2), rowMeans(norm.DE10.14[,11:20]), 
       pch=20, col='grey')
points(rowMeans(DE10.14.actualmeans.S1)[which(DE.DE10.14==1)], 
       rowMeans(norm.DE10.14[which(DE.DE10.14==1),1:10]), pch=20, col='red')
points(rowMeans(DE10.14.actualmeans.S2)[which(DE.DE10.14==1)], 
       rowMeans(norm.DE10.14[which(DE.DE10.14==1),11:20]), pch=20, col='red')
points(rowMeans(DE10.14.actualmeans.S1)[DE.edgeR.not.others], 
       rowMeans(norm.DE10.14[DE.edgeR.not.others,1:10]), pch=20, col='blue')
points(rowMeans(DE10.14.actualmeans.S2)[DE.edgeR.not.others], 
       rowMeans(norm.DE10.14[DE.edgeR.not.others,11:20]), pch=20, col='blue')
lines(c(0,3000), c(0,3000))
# Genes called correctly by edgeR and voom and wrongly by others still don't stand out from 
# other DE genes.

# Look at distributions of counts and sample variances
data.DE10.14@count.matrix[DE.edgeR.not.others,]
par(mfrow=c(4,2), mar=c(2,2,1,1), mgp=c(2,0.5,0))
boxplot(t(data.DE10.14@count.matrix[DE.edgeR.not.others,1:10]))
boxplot(t(data.DE10.14@count.matrix[DE.edgeR.not.others,11:20]))
# Is a lot of overlap between groups, but unusually so?
DE.all.correct <- which((call.baySeq.DE10.14==call.MDSeq.zi.DE10.14) & (call.baySeq.DE10.14==call.MDSeq.nozi.DE10.14) &
                          (call.baySeq.DE10.14==call.DESeq2.if.DE10.14) & (call.baySeq.DE10.14==call.DESeq2.noif.DE10.14) & 
                          (call.baySeq.DE10.14==call.DSS.notrend.DE10.14) & (call.baySeq.DE10.14==call.expHM.untr.DE10.14) & 
                          (call.baySeq.DE10.14==call.expHM.tr.DE10.14) & (call.baySeq.DE10.14==call.lnHM.untr.DE10.14) & 
                          (call.baySeq.DE10.14==call.lnHM.tr.DE10.14) & (call.baySeq.DE10.14==call.edgeR.ql.DE10.14) & 
                          (call.baySeq.DE10.14==call.edgeR.lr.DE10.14) & (call.baySeq.DE10.14==call.edgeR.et.DE10.14) & 
                          (call.baySeq.DE10.14==call.voom.DE10.14) & DE.DE10.14==1 & call.baySeq.DE10.14==1)
DE.all.wrong <- which((call.baySeq.DE10.14==call.MDSeq.zi.DE10.14) & (call.baySeq.DE10.14==call.MDSeq.nozi.DE10.14) &
                          (call.baySeq.DE10.14==call.DESeq2.if.DE10.14) & (call.baySeq.DE10.14==call.DESeq2.noif.DE10.14) & 
                          (call.baySeq.DE10.14==call.DSS.notrend.DE10.14) & (call.baySeq.DE10.14==call.expHM.untr.DE10.14) & 
                          (call.baySeq.DE10.14==call.expHM.tr.DE10.14) & (call.baySeq.DE10.14==call.lnHM.untr.DE10.14) & 
                          (call.baySeq.DE10.14==call.lnHM.tr.DE10.14) & (call.baySeq.DE10.14==call.edgeR.ql.DE10.14) & 
                          (call.baySeq.DE10.14==call.edgeR.lr.DE10.14) & (call.baySeq.DE10.14==call.edgeR.et.DE10.14) & 
                          (call.baySeq.DE10.14==call.voom.DE10.14) & DE.DE10.14==1 & call.baySeq.DE10.14==0)
boxplot(t(data.DE10.14@count.matrix[DE.all.correct[1:11],1:10]))
boxplot(t(data.DE10.14@count.matrix[DE.all.correct[1:11],11:20]))
boxplot(t(data.DE10.14@count.matrix[DE.all.wrong[1:11],1:10]))
boxplot(t(data.DE10.14@count.matrix[DE.all.wrong[1:11],11:20]))
# No obvious differences between genes that edgeR and voom call correctly and others wrongly, 
# genes that all call correctly, and genes that all call wrongly
for (i in 1:11) {
  print(round(c(t.test(data.DE10.14@count.matrix[DE.edgeR.not.others[i],1:10], 
                 data.DE10.14@count.matrix[DE.edgeR.not.others[i],11:20])$p.value, 
          t.test(data.DE10.14@count.matrix[DE.all.correct[i],1:10], 
                 data.DE10.14@count.matrix[DE.all.correct[i],11:20])$p.value, 
          t.test(data.DE10.14@count.matrix[DE.all.wrong[i],1:10], 
                 data.DE10.14@count.matrix[DE.all.wrong[i],11:20])$p.value), 3))
}
# Genes called correctly by edgeR and voom and wrongly by others all fail to be called by 
# t-test; genes called correctly by all are all called by t-test; genes called wrongly by 
# all all fail to be called by t-test.
# This says that for these genes at least, edgeR and voom do a better job of correctly 
# calling DE, but says nothing about how this differs between datsets.
# Could see if this explains the outliers by seeing whether there are more truly DE genes
# that are difficult to call in outlier datasets than in others. This is really what I've 
# been trying to do already, but so far haven't really looked at the distributions of counts
# for each gene. p-values from t-tests might give an easy proxy for difficulty in calling DE.


# Look at distribution of p-values for DE and non-DE genes in outlier an non-outlier datasets ####
for (i in outliers.DE10) {
  assign(paste0('ttest.DE10.',i,'_outlier'), numeric(length(get(paste0('DE.DE10.',i,'_outlier')))))
  for (j in 1:length(get(paste0('DE.DE10.',i,'_outlier')))) {
    assign(paste0('ttest.DE10.',i,'_outlier'), 
           `[<-`(get(paste0('ttest.DE10.',i,'_outlier')), j, 
                 value=t.test(get(paste0('counts.DE10.',i,'_outlier'))[j,1:10], 
                              get(paste0('counts.DE10.',i,'_outlier'))[j,11:20])$p.value))
  }
}
for (i in non.outliers.DE10) {
  assign(paste0('ttest.DE10.',i,'_non.outlier'), numeric(length(get(paste0('DE.DE10.',i,'_non.outlier')))))
  for (j in 1:length(get(paste0('DE.DE10.',i,'_non.outlier')))) {
    assign(paste0('ttest.DE10.',i,'_non.outlier'), 
           `[<-`(get(paste0('ttest.DE10.',i,'_non.outlier')), j, 
                 value=t.test(get(paste0('counts.DE10.',i,'_non.outlier'))[j,1:10], 
                              get(paste0('counts.DE10.',i,'_non.outlier'))[j,11:20])$p.value))
  }
}
par(mfrow=c(4,2), mar=c(2,2,1,1), mgp=c(2,0.5,0))
plot(density(ttest.DE10.7_outlier[which(DE.DE10.7_outlier == 1)]), main='')
plot(density(ttest.DE10.14_outlier[which(DE.DE10.14_outlier == 1)]), main='')
plot(density(ttest.DE10.15_outlier[which(DE.DE10.15_outlier == 1)]), main='')
plot(density(ttest.DE10.16_outlier[which(DE.DE10.16_outlier == 1)]), main='')
plot(density(ttest.DE10.1_non.outlier[which(DE.DE10.1_non.outlier == 1)]), main='')
plot(density(ttest.DE10.6_non.outlier[which(DE.DE10.6_non.outlier == 1)]), main='')
plot(density(ttest.DE10.20_non.outlier[which(DE.DE10.20_non.outlier == 1)]), main='')
plot(density(ttest.DE10.27_non.outlier[which(DE.DE10.27_non.outlier == 1)]), main='')
par(mfrow=c(4,2), mar=c(2,2,1,1), mgp=c(2,0.5,0))
plot(density(log(ttest.DE10.7_outlier)[which(DE.DE10.7_outlier == 1)]), main='')
plot(density(log(ttest.DE10.14_outlier)[which(DE.DE10.14_outlier == 1)]), main='')
plot(density(log(ttest.DE10.15_outlier)[which(DE.DE10.15_outlier == 1)]), main='')
plot(density(log(ttest.DE10.16_outlier)[which(DE.DE10.16_outlier == 1)]), main='')
plot(density(log(ttest.DE10.1_non.outlier)[which(DE.DE10.1_non.outlier == 1)]), main='')
plot(density(log(ttest.DE10.6_non.outlier)[which(DE.DE10.6_non.outlier == 1)]), main='')
plot(density(log(ttest.DE10.20_non.outlier)[which(DE.DE10.20_non.outlier == 1)]), main='')
plot(density(log(ttest.DE10.27_non.outlier)[which(DE.DE10.27_non.outlier == 1)]), main='')
c(mean(ttest.DE10.7_outlier[which(DE.DE10.7_outlier == 1)] > 0.05), 
  mean(ttest.DE10.14_outlier[which(DE.DE10.14_outlier == 1)] > 0.05), 
  mean(ttest.DE10.15_outlier[which(DE.DE10.15_outlier == 1)] > 0.05), 
  mean(ttest.DE10.16_outlier[which(DE.DE10.16_outlier == 1)] > 0.05), 
  mean(ttest.DE10.1_non.outlier[which(DE.DE10.1_non.outlier == 1)] > 0.05), 
  mean(ttest.DE10.6_non.outlier[which(DE.DE10.6_non.outlier == 1)] > 0.05), 
  mean(ttest.DE10.20_non.outlier[which(DE.DE10.20_non.outlier == 1)] > 0.05), 
  mean(ttest.DE10.27_non.outlier[which(DE.DE10.27_non.outlier == 1)] > 0.05))
# [1] 0.2732843 0.2872996 0.2628993 0.2506234 0.1908302 0.2148148 0.2121588 0.1961023
# No obvious differences in densities but clear difference in t-test p-values for
# outlier datasets. Seems to suggest that edgeR and voom are just better for difficult
# data.
# There is a question about the simulation though - seems odd that the mean t-test p-values
# should vary so much between datasets when there are around 800 DE genes. Should loook at 
# the distribution over all 50 datasets (and for other settings) and see whether there is 
# something odd about the values in the outlier datasets.
par(mfrow=c(4,2), mar=c(2,2,1,1), mgp=c(2,0.5,0))
plot(density(ttest.DE10.7_outlier[which(DE.DE10.7_outlier == 0)]), main='')
plot(density(ttest.DE10.14_outlier[which(DE.DE10.14_outlier == 0)]), main='')
plot(density(ttest.DE10.15_outlier[which(DE.DE10.15_outlier == 0)]), main='')
plot(density(ttest.DE10.16_outlier[which(DE.DE10.16_outlier == 0)]), main='')
plot(density(ttest.DE10.1_non.outlier[which(DE.DE10.1_non.outlier == 0)]), main='')
plot(density(ttest.DE10.6_non.outlier[which(DE.DE10.6_non.outlier == 0)]), main='')
plot(density(ttest.DE10.20_non.outlier[which(DE.DE10.20_non.outlier == 0)]), main='')
plot(density(ttest.DE10.27_non.outlier[which(DE.DE10.27_non.outlier == 0)]), main='')
c(mean(ttest.DE10.7_outlier[which(DE.DE10.7_outlier == 0)] < 0.05), 
  mean(ttest.DE10.14_outlier[which(DE.DE10.14_outlier == 0)] < 0.05), 
  mean(ttest.DE10.15_outlier[which(DE.DE10.15_outlier == 0)] < 0.05), 
  mean(ttest.DE10.16_outlier[which(DE.DE10.16_outlier == 0)] < 0.05), 
  mean(ttest.DE10.1_non.outlier[which(DE.DE10.1_non.outlier == 0)] < 0.05), 
  mean(ttest.DE10.6_non.outlier[which(DE.DE10.6_non.outlier == 0)] < 0.05), 
  mean(ttest.DE10.20_non.outlier[which(DE.DE10.20_non.outlier == 0)] < 0.05), 
  mean(ttest.DE10.27_non.outlier[which(DE.DE10.27_non.outlier == 0)] < 0.05))
# Very clear differences in distributions of t-test p-values for non-DE genes between 
# outlier and non-outlier datasets. Should also look at how this varies over all 
# datasets.

for (j in 1:50) {
  assign(paste0('results.DE10.',j), readRDS(here('Results/DE compcodeR data results July-Aug 2019',
                                                 paste0('results.DE10.',j,'.all.rds'))))
  assign(paste0('data.DE10.',j), get('data', get(paste0('results.DE10.',j))))
  assign(paste0('counts.DE10.',j), slot(get(paste0('data.DE10.',j)),'count.matrix'))
  assign(paste0('DE.DE10.',j), get('DE', get(paste0('results.DE10.',j))))
  assign(paste0('p.edgeR.ql.DE10.',j), get('p.ql.edgeR', get(paste0('results.DE10.',j))))
  assign(paste0('p.edgeR.lr.DE10.',j), get('p.lr.edgeR', get(paste0('results.DE10.',j))))
  assign(paste0('p.edgeR.et.DE10.',j), get('p.et.edgeR', get(paste0('results.DE10.',j))))
  assign(paste0('p.baySeq.DE10.',j), 1 - get('prob.baySeq', get(paste0('results.DE10.',j))))
  assign(paste0('p.MDSeq.zi.DE10.',j), get('p.mean.zi.MDSeq', get(paste0('results.DE10.',j))))
  assign(paste0('p.MDSeq.nozi.DE10.',j), get('p.mean.nozi.MDSeq', get(paste0('results.DE10.',j))))
  assign(paste0('p.DESeq2.if.DE10.',j), get('p.if.DESeq', get(paste0('results.DE10.',j))))
  assign(paste0('p.DESeq2.noif.DE10.',j), get('p.noif.DESeq', get(paste0('results.DE10.',j))))
  assign(paste0('p.voom.DE10.',j), get('p.voom', get(paste0('results.DE10.',j))))
  assign(paste0('p.DSS.notrend.DE10.',j), get('p.notrend.DSS', get(paste0('results.DE10.',j))))
  assign(paste0('p.expHM.untr.DE10.',j), get('p.mean.expHM', get(paste0('results.DE10.',j))))
  assign(paste0('p.expHM.tr.DE10.',j), get('p.lmean.expHM', get(paste0('results.DE10.',j))))
  assign(paste0('p.lnHM.untr.DE10.',j), get('p.mean.lnHM', get(paste0('results.DE10.',j))))
  assign(paste0('p.lnHM.tr.DE10.',j), get('p.lmean.lnHM', get(paste0('results.DE10.',j))))
  for (i in c(paste0('edgeR.ql.DE10.',j), paste0('edgeR.lr.DE10.',j),paste0('edgeR.et.DE10.',j),
              paste0('baySeq.DE10.',j), paste0('MDSeq.zi.DE10.',j), paste0('MDSeq.nozi.DE10.',j),
              paste0('DESeq2.if.DE10.',j), paste0('DESeq2.noif.DE10.',j), paste0('voom.DE10.',j),
              paste0('DSS.notrend.DE10.',j), paste0('expHM.untr.DE10.',j), paste0('expHM.tr.DE10.',j),
              paste0('lnHM.untr.DE10.',j), paste0('lnHM.tr.DE10.',j))) {
    assign(paste0('pred.',i), prediction(1-get(paste0('p.',i)), 
                                                    get(paste0('DE.DE10.',j))))
    assign(paste0('auc.',i), performance(get(paste0('pred.',i)), 
                                                    measure='auc')@y.values[[1]])
    rm(list=paste0('pred.',i))
  }
  rm(list=paste0('results.DE10.',j))
  assign(paste0('ttest.DE10.',j), numeric(length(get(paste0('DE.DE10.',j)))))
  for (k in 1:length(get(paste0('DE.DE10.',j)))) {
    assign(paste0('ttest.DE10.',j), 
           `[<-`(get(paste0('ttest.DE10.',j)), k, 
                 value=t.test(get(paste0('counts.DE10.',j))[k,1:10], 
                              get(paste0('counts.DE10.',j))[k,11:20])$p.value))
  }
}
prop.nonsig.ttest.DE <- numeric(50)
prop.sig.ttest.noDE <- numeric(50)
for (i in 1:50) {
  prop.nonsig.ttest.DE[i] <- mean(get(paste0('ttest.DE10.',i))[
    which(get(paste0('DE.DE10.',i)) == 1)] > 0.05)
  prop.sig.ttest.noDE[i] <- mean(get(paste0('ttest.DE10.',i))[
    which(get(paste0('DE.DE10.',i)) == 0)] < 0.05)
}
col=rep('black',50); col[outliers.DE10] <- 'red'; col[non.outliers.DE10] <- 'blue'
plot(prop.nonsig.ttest.DE, col=col, pch=16)
plot(prop.sig.ttest.noDE, col=col, pch=16)
# Outlier datasets are among those with highest proportions of non-significant 
# t-tests for DE genes, but not really clearly outside the normal range, but 
# for proportions of significant t-tests for non-DE genes, outliers are the 
# highest, and two of them are clearly outside what seems to be the normal range.
auc.DE10 <- data.frame(matrix(nrow=50, ncol=14))
names(auc.DE10) <- c('edgeR.ql', 'edgeR.lr', 'edgeR.et', 'baySeq', 'MDSeq.zi', 'MDSeq.nozi', 'DESeq2.if', 
                        'DESeq2.noif', 'voom', 'DSS.notrend', 'expHM.untr', 'expHM.tr', 'lnHM.untr', 'lnHM.tr')
for (i in 1:50) {
  for (j in names(auc.DE10.up)) {
    auc.DE10[[j]][i] <- get(paste0('auc.',j,'.DE10.',i))
  }
}
par(mfrow=c(7,4), mar=c(1.5,1.5,0.5,0.5), mgp=c(2,0,0))
for (i in names(auc.DE10)) {
  plot(prop.nonsig.ttest.DE, auc.DE10[[i]], pch=20, main=i, ylim=c(0.83,1))
  plot(prop.sig.ttest.noDE, auc.DE10[[i]], pch=20, main=i, ylim=c(0.83,1))
}
# Clear trend of AUCs decreasing with proportion of nonsignificant t-tests among 
# DE genes and proportion of significant t-tests among non-DE genes for all methods 
# except edgeR and voom.
# Wonder if TMM normalisation may be favouring edgeR and voom. Doubt it but should 
# rule it out.
library(limma)
library(edgeR)
library(DESeq2)
library(DSS)
group <- factor(c(rep(1,10), rep(2,10)))
design <- model.matrix(~group)
for (i in 1:50) {
  # Load data
  filename <- paste0('DE10.', i)
  counts <- readRDS(here('Simulated data',paste0(filename,'.rds')))
  DE <- counts@variable.annotations$differential.expression
  
  # Normalise and create DGEList object
  # nf <- calcNormFactors(counts@count.matrix)
  # nf <- rep(1,20)
  dat.DESeq <- DESeqDataSetFromMatrix(countData=counts@count.matrix, colData=data.frame(group), 
                                      design=~group)
  dat.DESeq <- estimateSizeFactors(dat.DESeq)
  nf <- dat.DESeq$sizeFactor
  norm <- t(t(counts@count.matrix) / nf)
  dge <- DGEList(counts=counts@count.matrix, norm.factors=nf, group=group)
  
  
  ## edgeR
  dat.edgeR <- estimateDisp(dge, design)
  qlfit.edgeR <- glmQLFit(dat.edgeR, design)
  qltest.edgeR <- glmQLFTest(qlfit.edgeR)
  lrfit.edgeR <- glmFit(dat.edgeR, design)
  lrtest.edgeR <- glmLRT(lrfit.edgeR)
  et.edgeR <- exactTest(dat.edgeR)
  p.ql.edgeR <- qltest.edgeR$table$PValue
  p.lr.edgeR <- lrtest.edgeR$table$PValue
  p.et.edgeR <- et.edgeR$table$PValue
  pred.ql.edgeR <- prediction(1-p.ql.edgeR, DE)
  pred.lr.edgeR <- prediction(1-p.lr.edgeR, DE)
  pred.et.edgeR <- prediction(1-p.et.edgeR, DE)
  auc.ql.edgeR <- performance(pred.ql.edgeR, measure='auc')@y.values[[1]]
  auc.lr.edgeR <- performance(pred.lr.edgeR, measure='auc')@y.values[[1]]
  auc.et.edgeR <- performance(pred.et.edgeR, measure='auc')@y.values[[1]]
  rm(list=c('dat.edgeR', 'qlfit.edgeR', 'qltest.edgeR', 'lrfit.edgeR', 'lrtest.edgeR', 'et.edgeR', 
            'pred.ql.edgeR', 'pred.lr.edgeR', 'pred.et.edgeR'))
  
  ## DESeq2
  # dat.DESeq <- DESeqDataSetFromMatrix(countData=counts@count.matrix, colData=data.frame(group), 
  #                                     design=~group)
  # dat.DESeq <- estimateSizeFactors(dat.DESeq)
  # dat.DESeq$sizeFactor <- nf
  dat.DESeq <- DESeq(dat.DESeq, minReplicatesForReplace=Inf)
  res.noif.DESeq <- results(dat.DESeq, independentFiltering=F, cooksCutoff=F)
  res.if.DESeq <- results(dat.DESeq, cooksCutoff=F, alpha=0.05)
  p.noif.DESeq <- res.noif.DESeq$pvalue
  p.if.DESeq <- res.noif.DESeq$pvalue
  pred.noif.DESeq <- prediction(1-p.noif.DESeq, DE)
  pred.if.DESeq <- prediction(1-p.if.DESeq, DE)
  auc.noif.DESeq <- performance(pred.noif.DESeq, measure='auc')@y.values[[1]]
  auc.if.DESeq <- performance(pred.if.DESeq, measure='auc')@y.values[[1]]
  rm(list=c('dat.DESeq', 'res.noif.DESeq', 'res.if.DESeq', 'pred.noif.DESeq', 'pred.if.DESeq'))
  
  ## limma-voom
  dat.voom <- voom(dge)
  fit.voom <- lmFit(dat.voom, design)
  res.voom <- eBayes(fit.voom)
  p.voom <- topTable(res.voom, number=Inf, sort.by='none')$P.Value
  pred.voom <- prediction(1-p.voom, DE)
  auc.voom <- performance(pred.voom, measure='auc')@y.values[[1]]
  rm(list=c('dat.voom', 'fit.voom', 'res.voom', 'pred.voom'))
  
  ## DSS
  dat.DSS <- newSeqCountSet(counts=matrix(counts@count.matrix, ncol=length(group)), 
                            designs=as.numeric(group), normalizationFactor=nf)
  notrend.DSS <- estDispersion(dat.DSS)
  res.notrend.DSS <- waldTest(notrend.DSS, 1, 2)[order(waldTest(notrend.DSS, 1, 2)$geneIndex),]
  p.notrend.DSS <- res.notrend.DSS$pval
  pred.notrend.DSS <- prediction(1-p.notrend.DSS, DE)
  auc.notrend.DSS <- performance(pred.notrend.DSS, measure='auc')@y.values[[1]]
  rm(list=c('dat.DSS', 'notrend.DSS', 'res.notrend.DSS', 'pred.notrend.DSS'))
  
  results <- list(data = counts, 
                  DE = DE, 
                  p.ql.edgeR = p.ql.edgeR, 
                  p.lr.edgeR = p.lr.edgeR, 
                  p.et.edgeR = p.et.edgeR, 
                  p.noif.DESeq = p.noif.DESeq, 
                  p.if.DESeq = p.if.DESeq, 
                  p.voom = p.voom, 
                  p.notrend.DSS = p.notrend.DSS, 
                  auc.ql.edgeR = auc.ql.edgeR, 
                  auc.lr.edgeR = auc.lr.edgeR, 
                  auc.et.edgeR = auc.et.edgeR, 
                  auc.noif.DESeq = auc.noif.DESeq, 
                  auc.if.DESeq = auc.if.DESeq, 
                  auc.voom = auc.voom, 
                  auc.notrend.DSS = auc.notrend.DSS)
  # assign(paste0('res.no.norm.',i), results)
  assign(paste0('res.DESeq.norm.',i), results)
  rm(list=c('counts', 'DE', 'nf', 'norm', 'dge', 'results'))
}
auc.DE10.no.norm <- data.frame(matrix(nrow=50, ncol=7))
names(auc.DE10.no.norm) <- c('ql.edgeR', 'lr.edgeR', 'et.edgeR', 'if.DESeq', 'noif.DESeq', 
                             'voom', 'notrend.DSS')
for (i in 1:50) {
  for (j in names(auc.DE10.no.norm)) {
    auc.DE10.no.norm[[j]][i] <- unlist(get(paste0('res.no.norm.',i))[paste0('auc.',j)])
  }
}
boxplot(auc.DE10.no.norm)
boxplot(auc.DE10)
# No normalisation gives almost identical results to TMM
auc.DE10.DESeq.norm <- data.frame(matrix(nrow=50, ncol=7))
names(auc.DE10.DESeq.norm) <- c('ql.edgeR', 'lr.edgeR', 'et.edgeR', 'if.DESeq', 'noif.DESeq', 
                             'voom', 'notrend.DSS')
for (i in 1:50) {
  for (j in names(auc.DE10.DESeq.norm)) {
    auc.DE10.DESeq.norm[[j]][i] <- unlist(get(paste0('res.DESeq.norm.',i))[paste0('auc.',j)])
  }
}
boxplot(auc.DE10.DESeq.norm)
boxplot(auc.DE10)
# Normalisation!!!!!!!!!!!!!!!!!
# Outlier datasets look to generally now be outliers for edgeR and voom, but not for 
# DESeq2 and DSS.



#### ####



