library(here)
library(compcodeR)
library(ggplot2)
library(tidyr)

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
library(ROCR)
library(caret)
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






