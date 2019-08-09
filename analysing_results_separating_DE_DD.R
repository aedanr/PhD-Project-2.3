library(here)

for (i in c('DE_noDD', 'DE_noDD.lfc1', 'DE_noDD.lfc2', 
            'DE_withDD', 'DE_withDD.lfc1', 'DE_withDD.lfc2', 
            'DD_noDE', 'DD_noDE.lfc1', 'DD_noDE.lfc2', 
            'DD_withDE', 'DD_withDE.lfc1', 'DD_withDE.lfc2', 
            'DEDD_DEandDD', 'DEDD_DEnoDD', 'DEDD_DDnoDE')) {
  for (j in c('DEDD2', 'DEDD5', 'DEDD10')) {
    assign(paste0(i,'.results.',j), readRDS(here('Results/DEDD compcodeR data results July-Aug 2019', 
                                     paste0(i,'.results.',j,'.rds'))))
  }
}
rm('i','j')


# DE ####
# AUC
par(mfcol=c(3,2), mar=c(3,3,1,1))
boxplot(DE_noDD.results.DEDD2$auc, names=1:17, ylim=range(c(DE_noDD.results.DEDD2$auc, DE_withDD.results.DEDD2$auc)))
abline(h=c(0.6,0.7,0.8), col='grey')
boxplot(DE_noDD.results.DEDD5$auc, names=1:14, ylim=range(c(DE_noDD.results.DEDD5$auc, DE_withDD.results.DEDD5$auc)))
abline(h=c(0.8,0.85,0.9), col='grey')
boxplot(DE_noDD.results.DEDD10$auc, names=1:14, ylim=range(c(DE_noDD.results.DEDD10$auc, DE_withDD.results.DEDD10$auc)))
abline(h=c(0.9,0.92,0.94,0.96), col='grey')
boxplot(DE_withDD.results.DEDD2$auc, names=1:17, ylim=range(c(DE_noDD.results.DEDD2$auc, DE_withDD.results.DEDD2$auc)))
abline(h=c(0.6,0.7,0.8), col='grey')
boxplot(DE_withDD.results.DEDD5$auc, names=1:14, ylim=range(c(DE_noDD.results.DEDD5$auc, DE_withDD.results.DEDD5$auc)))
abline(h=c(0.8,0.85,0.9), col='grey')
boxplot(DE_withDD.results.DEDD10$auc, names=1:14, ylim=range(c(DE_noDD.results.DEDD10$auc, DE_withDD.results.DEDD10$auc)))
abline(h=c(0.9,0.92,0.94,0.96), col='grey')
cbind(round(colMeans(DE_noDD.results.DEDD2$auc),3),
      round(colMeans(DE_withDD.results.DEDD2$auc),3))
cbind(round(colMeans(DE_noDD.results.DEDD5$auc),3),
      round(colMeans(DE_withDD.results.DEDD5$auc),3))
cbind(round(colMeans(DE_noDD.results.DEDD10$auc),3),
      round(colMeans(DE_withDD.results.DEDD10$auc),3))

# FDR curves
par(mfrow=c(6,3), mar=c(2.5,3,1,1), mgp=c(1.8,0.7,0))
for (i in 1:17) {
  plot(DE_noDD.results.DEDD2$mean.discoveries[which(DE_noDD.results.DEDD2$mean.discoveries[,i]<1000),i], 
       DE_noDD.results.DEDD2$mean.fdr[which(DE_noDD.results.DEDD2$mean.discoveries[,i]<1000),i], 
       type='l', xlab='', ylab=names(DE_noDD.results.DEDD2$mean.fdr)[i], ylim=c(0,0.8))
  lines(DE_withDD.results.DEDD2$mean.discoveries[which(DE_withDD.results.DEDD2$mean.discoveries[,i]<1000),i], 
       DE_withDD.results.DEDD2$mean.fdr[which(DE_withDD.results.DEDD2$mean.discoveries[,i]<1000),i], 
       xlab='', col='grey', ylab=names(DE_withDD.results.DEDD2$mean.fdr)[i], ylim=c(0,0.8))
}
par(mfrow=c(6,3), mar=c(2.5,3,1,1), mgp=c(1.8,0.7,0))
for (i in 1:14) {
  plot(DE_noDD.results.DEDD5$mean.discoveries[which(DE_noDD.results.DEDD5$mean.discoveries[,i]<1000),i], 
       DE_noDD.results.DEDD5$mean.fdr[which(DE_noDD.results.DEDD5$mean.discoveries[,i]<1000),i], 
       type='l',xlab='', ylab=names(DE_noDD.results.DEDD5$mean.fdr)[i], ylim=c(0,0.4))
  lines(DE_withDD.results.DEDD5$mean.discoveries[which(DE_withDD.results.DEDD5$mean.discoveries[,i]<1000),i], 
       DE_withDD.results.DEDD5$mean.fdr[which(DE_withDD.results.DEDD5$mean.discoveries[,i]<1000),i], 
       xlab='', col='grey', ylab=names(DE_withDD.results.DEDD5$mean.fdr)[i], ylim=c(0,0.4))
}
par(mfrow=c(6,3), mar=c(2.5,3,1,1), mgp=c(1.8,0.7,0))
for (i in 1:14) {
  plot(DE_noDD.results.DEDD10$mean.discoveries[which(DE_noDD.results.DEDD10$mean.discoveries[,i]<1000),i], 
       DE_noDD.results.DEDD10$mean.fdr[which(DE_noDD.results.DEDD10$mean.discoveries[,i]<1000),i], 
       type='l', xlab='', ylab=names(DE_noDD.results.DEDD10$mean.fdr)[i], ylim=c(0,0.2))
  lines(DE_withDD.results.DEDD10$mean.discoveries[which(DE_withDD.results.DEDD10$mean.discoveries[,i]<1000),i], 
       DE_withDD.results.DEDD10$mean.fdr[which(DE_withDD.results.DEDD10$mean.discoveries[,i]<1000),i], 
       xlab='', col='grey', ylab=names(DE_withDD.results.DEDD10$mean.fdr)[i], ylim=c(0,0.2))
}

# FDR
par(mfcol=c(3,2), mar=c(2.5,2.5,0.5,0.5))
boxplot(DE_noDD.results.DEDD2$fdr[-grep('raw',names(DE_noDD.results.DEDD2$fdr))], 
        names=1:32, ylim=c(0,range(c(DE_noDD.results.DEDD2$fdr,DE_withDD.results.DEDD2$fdr), 
                                   na.rm=T)[2])); abline(h=c(0.05,0.1), col='grey')
boxplot(DE_noDD.results.DEDD5$fdr[-grep('raw',names(DE_noDD.results.DEDD5$fdr))], 
        names=1:25, ylim=c(0,range(c(DE_noDD.results.DEDD5$fdr,DE_withDD.results.DEDD5$fdr), 
                                   na.rm=T)[2])); abline(h=c(0.05,0.1), col='grey')
boxplot(DE_noDD.results.DEDD10$fdr[-grep('raw',names(DE_noDD.results.DEDD10$fdr))], 
        names=1:25, ylim=c(0,range(c(DE_noDD.results.DEDD10$fdr,DE_withDD.results.DEDD10$fdr), 
                                   na.rm=T)[2])); abline(h=c(0.05,0.1), col='grey')
boxplot(DE_withDD.results.DEDD2$fdr[-grep('raw',names(DE_withDD.results.DEDD2$fdr))], 
        names=1:32, ylim=c(0,range(c(DE_noDD.results.DEDD2$fdr,DE_withDD.results.DEDD2$fdr), 
                                   na.rm=T)[2])); abline(h=c(0.05,0.1), col='grey')
boxplot(DE_withDD.results.DEDD5$fdr[-grep('raw',names(DE_withDD.results.DEDD5$fdr))], 
        names=1:25, ylim=c(0,range(c(DE_noDD.results.DEDD5$fdr,DE_withDD.results.DEDD5$fdr), 
                                   na.rm=T)[2])); abline(h=c(0.05,0.1), col='grey')
boxplot(DE_withDD.results.DEDD10$fdr[-grep('raw',names(DE_withDD.results.DEDD10$fdr))], 
        names=1:25, ylim=c(0,range(c(DE_noDD.results.DEDD10$fdr,DE_withDD.results.DEDD10$fdr), 
                                   na.rm=T)[2])); abline(h=c(0.05,0.1), col='grey')
cbind(round(colMeans(DE_noDD.results.DEDD2$fdr[-grep('raw',names(DE_noDD.results.DEDD2$fdr))]),3),
      round(colMeans(DE_withDD.results.DEDD2$fdr[-grep('raw',names(DE_withDD.results.DEDD2$fdr))]),3))
cbind(round(colMeans(DE_noDD.results.DEDD5$fdr[-grep('raw',names(DE_noDD.results.DEDD5$fdr))]),3),
      round(colMeans(DE_withDD.results.DEDD5$fdr[-grep('raw',names(DE_withDD.results.DEDD5$fdr))]),3))
cbind(round(colMeans(DE_noDD.results.DEDD10$fdr[-grep('raw',names(DE_noDD.results.DEDD10$fdr))]),3),
      round(colMeans(DE_withDD.results.DEDD10$fdr[-grep('raw',names(DE_withDD.results.DEDD10$fdr))]),3))

# TPR
par(mfcol=c(3,2), mar=c(2.5,2.5,0.5,0.5))
boxplot(DE_noDD.results.DEDD2$tpr[-grep('raw',names(DE_noDD.results.DEDD2$tpr))], 
        names=1:32, ylim=c(0,range(c(DE_noDD.results.DEDD2$tpr,DE_withDD.results.DEDD2$tpr), 
                                   na.rm=T)[2])); abline(h=c(0.1,0.2,0.3), col='grey')
boxplot(DE_noDD.results.DEDD5$tpr[-grep('raw',names(DE_noDD.results.DEDD5$tpr))], 
        names=1:25, ylim=c(0,range(c(DE_noDD.results.DEDD5$tpr,DE_withDD.results.DEDD5$tpr), 
                                   na.rm=T)[2])); abline(h=c(0.1,0.3,0.5), col='grey')
boxplot(DE_noDD.results.DEDD10$tpr[-grep('raw',names(DE_noDD.results.DEDD10$tpr))], 
        names=1:25, ylim=c(0,range(c(DE_noDD.results.DEDD10$tpr,DE_withDD.results.DEDD10$tpr), 
                                   na.rm=T)[2])); abline(h=c(0.4,0.55,0.7), col='grey')
boxplot(DE_withDD.results.DEDD2$tpr[-grep('raw',names(DE_withDD.results.DEDD2$tpr))], 
        names=1:32, ylim=c(0,range(c(DE_noDD.results.DEDD2$tpr,DE_withDD.results.DEDD2$tpr), 
                                   na.rm=T)[2])); abline(h=c(0.1,0.2,0.3), col='grey')
boxplot(DE_withDD.results.DEDD5$tpr[-grep('raw',names(DE_withDD.results.DEDD5$tpr))], 
        names=1:25, ylim=c(0,range(c(DE_noDD.results.DEDD5$tpr,DE_withDD.results.DEDD5$tpr), 
                                   na.rm=T)[2])); abline(h=c(0.1,0.3,0.5), col='grey')
boxplot(DE_withDD.results.DEDD10$tpr[-grep('raw',names(DE_withDD.results.DEDD10$tpr))], 
        names=1:25, ylim=c(0,range(c(DE_noDD.results.DEDD10$tpr,DE_withDD.results.DEDD10$tpr), 
                                   na.rm=T)[2])); abline(h=c(0.4,0.55,0.7), col='grey')
cbind(round(colMeans(DE_noDD.results.DEDD2$tpr[-grep('raw',names(DE_noDD.results.DEDD2$tpr))]),3),
      round(colMeans(DE_withDD.results.DEDD2$tpr[-grep('raw',names(DE_withDD.results.DEDD2$tpr))]),3))
cbind(round(colMeans(DE_noDD.results.DEDD5$tpr[-grep('raw',names(DE_noDD.results.DEDD5$tpr))]),3),
      round(colMeans(DE_withDD.results.DEDD5$tpr[-grep('raw',names(DE_withDD.results.DEDD5$tpr))]),3))
cbind(round(colMeans(DE_noDD.results.DEDD10$tpr[-grep('raw',names(DE_noDD.results.DEDD10$tpr))]),3),
      round(colMeans(DE_withDD.results.DEDD10$tpr[-grep('raw',names(DE_withDD.results.DEDD10$tpr))]),3))


# DE lfc1 ####
# AUC
par(mfcol=c(3,2), mar=c(2.5,2.5,0.5,0.5))
boxplot(DE_noDD.lfc1.results.DEDD2$auc, names=1:9, ylim=c(0.5,1))
abline(h=c(0.7,0.8,0.9), col='grey')
boxplot(DE_noDD.lfc1.results.DEDD5$auc, names=1:8, ylim=c(0.5,1))
abline(h=c(0.7,0.8,0.9), col='grey')
boxplot(DE_noDD.lfc1.results.DEDD10$auc, names=1:8, ylim=c(0.5,1))
abline(h=c(0.7,0.8,0.9), col='grey')
boxplot(DE_withDD.lfc1.results.DEDD2$auc, names=1:9, ylim=c(0.5,1))
abline(h=c(0.7,0.8,0.9), col='grey')
boxplot(DE_withDD.lfc1.results.DEDD5$auc, names=1:8, ylim=c(0.5,1))
abline(h=c(0.7,0.8,0.9), col='grey')
boxplot(DE_withDD.lfc1.results.DEDD10$auc, names=1:8, ylim=c(0.5,1))
abline(h=c(0.7,0.8,0.9), col='grey')
cbind(round(colMeans(DE_noDD.lfc1.results.DEDD2$auc),3),
      round(colMeans(DE_withDD.lfc1.results.DEDD2$auc),3))
cbind(round(colMeans(DE_noDD.lfc1.results.DEDD5$auc),3),
      round(colMeans(DE_withDD.lfc1.results.DEDD5$auc),3))
cbind(round(colMeans(DE_noDD.lfc1.results.DEDD10$auc),3),
      round(colMeans(DE_withDD.lfc1.results.DEDD10$auc),3))

# FDR
par(mfcol=c(3,2), mar=c(2.5,2.5,0.5,0.5))
boxplot(DE_noDD.lfc1.results.DEDD2$fdr[-grep('raw',names(DE_noDD.lfc1.results.DEDD2$fdr))], 
        names=1:14, ylim=c(0,1)); abline(h=c(0.05,0.1), col='grey')
boxplot(DE_noDD.lfc1.results.DEDD5$fdr[-grep('raw',names(DE_noDD.lfc1.results.DEDD5$fdr))], 
        names=1:12, ylim=c(0,0.4)); abline(h=c(0.05,0.1), col='grey')
boxplot(DE_noDD.lfc1.results.DEDD10$fdr[-grep('raw',names(DE_noDD.lfc1.results.DEDD10$fdr))], 
        names=1:12, ylim=c(0,0.1)); abline(h=c(0.05,0.1), col='grey')
boxplot(DE_withDD.lfc1.results.DEDD2$fdr[-grep('raw',names(DE_withDD.lfc1.results.DEDD2$fdr))], 
        names=1:14, ylim=c(0,1)); abline(h=c(0.05,0.1), col='grey')
boxplot(DE_withDD.lfc1.results.DEDD5$fdr[-grep('raw',names(DE_withDD.lfc1.results.DEDD5$fdr))], 
        names=1:12, ylim=c(0,0.4)); abline(h=c(0.05,0.1), col='grey')
boxplot(DE_withDD.lfc1.results.DEDD10$fdr[-grep('raw',names(DE_withDD.lfc1.results.DEDD10$fdr))], 
        names=1:12, ylim=c(0,0.1)); abline(h=c(0.05,0.1), col='grey')
cbind(round(colMeans(DE_noDD.lfc1.results.DEDD2$fdr[-grep('raw',names(DE_noDD.lfc1.results.DEDD2$fdr))]),3),
      round(colMeans(DE_withDD.lfc1.results.DEDD2$fdr[-grep('raw',names(DE_withDD.lfc1.results.DEDD2$fdr))]),3))
cbind(round(colMeans(DE_noDD.lfc1.results.DEDD5$fdr[-grep('raw',names(DE_noDD.lfc1.results.DEDD5$fdr))]),3),
      round(colMeans(DE_withDD.lfc1.results.DEDD5$fdr[-grep('raw',names(DE_withDD.lfc1.results.DEDD5$fdr))]),3))
cbind(round(colMeans(DE_noDD.lfc1.results.DEDD10$fdr[-grep('raw',names(DE_noDD.lfc1.results.DEDD10$fdr))]),3),
      round(colMeans(DE_withDD.lfc1.results.DEDD10$fdr[-grep('raw',names(DE_withDD.lfc1.results.DEDD10$fdr))]),3))

# TPR
par(mfcol=c(3,2), mar=c(2.5,2.5,0.5,0.5))
boxplot(DE_noDD.lfc1.results.DEDD2$tpr[-grep('raw',names(DE_noDD.lfc1.results.DEDD2$tpr))], 
        names=1:14, ylim=c(0,0.15)); abline(h=c(0.05,0.1), col='grey')
boxplot(DE_noDD.lfc1.results.DEDD5$tpr[-grep('raw',names(DE_noDD.lfc1.results.DEDD5$tpr))], 
        names=1:12, ylim=c(0,0.15)); abline(h=c(0.05,0.1), col='grey')
boxplot(DE_noDD.lfc1.results.DEDD10$tpr[-grep('raw',names(DE_noDD.lfc1.results.DEDD10$tpr))], 
        names=1:12, ylim=c(0,0.25)); abline(h=c(0.05,0.1), col='grey')
boxplot(DE_withDD.lfc1.results.DEDD2$tpr[-grep('raw',names(DE_withDD.lfc1.results.DEDD2$tpr))], 
        names=1:14, ylim=c(0,0.15)); abline(h=c(0.05,0.1), col='grey')
boxplot(DE_withDD.lfc1.results.DEDD5$tpr[-grep('raw',names(DE_withDD.lfc1.results.DEDD5$tpr))], 
        names=1:12, ylim=c(0,0.15)); abline(h=c(0.05,0.1), col='grey')
boxplot(DE_withDD.lfc1.results.DEDD10$tpr[-grep('raw',names(DE_withDD.lfc1.results.DEDD10$tpr))], 
        names=1:12, ylim=c(0,0.25)); abline(h=c(0.05,0.1), col='grey')
cbind(round(colMeans(DE_noDD.lfc1.results.DEDD2$tpr[-grep('raw',names(DE_noDD.lfc1.results.DEDD2$tpr))]),3),
      round(colMeans(DE_withDD.lfc1.results.DEDD2$tpr[-grep('raw',names(DE_withDD.lfc1.results.DEDD2$tpr))]),3))
cbind(round(colMeans(DE_noDD.lfc1.results.DEDD5$tpr[-grep('raw',names(DE_noDD.lfc1.results.DEDD5$tpr))]),3),
      round(colMeans(DE_withDD.lfc1.results.DEDD5$tpr[-grep('raw',names(DE_withDD.lfc1.results.DEDD5$tpr))]),3))
cbind(round(colMeans(DE_noDD.lfc1.results.DEDD10$tpr[-grep('raw',names(DE_noDD.lfc1.results.DEDD10$tpr))]),3),
      round(colMeans(DE_withDD.lfc1.results.DEDD10$tpr[-grep('raw',names(DE_withDD.lfc1.results.DEDD10$tpr))]),3))


# DE lfc2 ####
# AUC
par(mfcol=c(3,2), mar=c(2.5,2.5,0.5,0.5))
boxplot(DE_noDD.lfc2.results.DEDD2$auc, names=1:9, ylim=c(0.5,1))
abline(h=c(0.6,0.7,0.8), col='grey')
boxplot(DE_noDD.lfc2.results.DEDD5$auc, names=1:8, ylim=c(0.5,1))
abline(h=c(0.6,0.7,0.8), col='grey')
boxplot(DE_noDD.lfc2.results.DEDD10$auc, names=1:8, ylim=c(0.5,1))
abline(h=c(0.6,0.7,0.8), col='grey')
boxplot(DE_withDD.lfc2.results.DEDD2$auc, names=1:9, ylim=c(0.5,1))
abline(h=c(0.6,0.7,0.8), col='grey')
boxplot(DE_withDD.lfc2.results.DEDD5$auc, names=1:8, ylim=c(0.5,1))
abline(h=c(0.6,0.7,0.8), col='grey')
boxplot(DE_withDD.lfc2.results.DEDD10$auc, names=1:8, ylim=c(0.5,1))
abline(h=c(0.6,0.7,0.8), col='grey')
cbind(round(colMeans(DE_noDD.lfc2.results.DEDD2$auc),3),
      round(colMeans(DE_withDD.lfc2.results.DEDD2$auc),3))
cbind(round(colMeans(DE_noDD.lfc2.results.DEDD5$auc),3),
      round(colMeans(DE_withDD.lfc2.results.DEDD5$auc),3))
cbind(round(colMeans(DE_noDD.lfc2.results.DEDD10$auc),3),
      round(colMeans(DE_withDD.lfc2.results.DEDD10$auc),3))

# FDR
par(mfcol=c(3,2), mar=c(2.5,2.5,0.5,0.5))
boxplot(DE_noDD.lfc2.results.DEDD2$fdr[-grep('raw',names(DE_noDD.lfc2.results.DEDD2$fdr))], 
        names=1:14, ylim=c(0,1)); abline(h=c(0.05,0.1), col='grey')
boxplot(DE_noDD.lfc2.results.DEDD5$fdr[-grep('raw',names(DE_noDD.lfc2.results.DEDD5$fdr))], 
        names=1:12, ylim=c(0,1)); abline(h=c(0.05,0.1), col='grey')
boxplot(DE_noDD.lfc2.results.DEDD10$fdr[-grep('raw',names(DE_noDD.lfc2.results.DEDD10$fdr))], 
        names=1:12, ylim=c(0,1)); abline(h=c(0.05,0.1), col='grey')
boxplot(DE_withDD.lfc2.results.DEDD2$fdr[-grep('raw',names(DE_withDD.lfc2.results.DEDD2$fdr))], 
        names=1:14, ylim=c(0,1)); abline(h=c(0.05,0.1), col='grey')
boxplot(DE_withDD.lfc2.results.DEDD5$fdr[-grep('raw',names(DE_withDD.lfc2.results.DEDD5$fdr))], 
        names=1:12, ylim=c(0,1)); abline(h=c(0.05,0.1), col='grey')
boxplot(DE_withDD.lfc2.results.DEDD10$fdr[-grep('raw',names(DE_withDD.lfc2.results.DEDD10$fdr))], 
        names=1:12, ylim=c(0,1)); abline(h=c(0.05,0.1), col='grey')
cbind(round(colMeans(DE_noDD.lfc2.results.DEDD2$fdr[-grep('raw',names(DE_noDD.lfc2.results.DEDD2$fdr))]),3),
      round(colMeans(DE_withDD.lfc2.results.DEDD2$fdr[-grep('raw',names(DE_withDD.lfc2.results.DEDD2$fdr))]),3))
cbind(round(colMeans(DE_noDD.lfc2.results.DEDD5$fdr[-grep('raw',names(DE_noDD.lfc2.results.DEDD5$fdr))]),3),
      round(colMeans(DE_withDD.lfc2.results.DEDD5$fdr[-grep('raw',names(DE_withDD.lfc2.results.DEDD5$fdr))]),3))
cbind(round(colMeans(DE_noDD.lfc2.results.DEDD10$fdr[-grep('raw',names(DE_noDD.lfc2.results.DEDD10$fdr))]),3),
      round(colMeans(DE_withDD.lfc2.results.DEDD10$fdr[-grep('raw',names(DE_withDD.lfc2.results.DEDD10$fdr))]),3))

# TPR
par(mfcol=c(3,2), mar=c(2.5,2.5,0.5,0.5))
boxplot(DE_noDD.lfc2.results.DEDD2$tpr[-grep('raw',names(DE_noDD.lfc2.results.DEDD2$tpr))], 
        names=1:14, ylim=c(0,0.05)); abline(h=c(0.05,0.1), col='grey')
boxplot(DE_noDD.lfc2.results.DEDD5$tpr[-grep('raw',names(DE_noDD.lfc2.results.DEDD5$tpr))], 
        names=1:12, ylim=c(0,0.05)); abline(h=c(0.05,0.1), col='grey')
boxplot(DE_noDD.lfc2.results.DEDD10$tpr[-grep('raw',names(DE_noDD.lfc2.results.DEDD10$tpr))], 
        names=1:12, ylim=c(0,0.05)); abline(h=c(0.05,0.1), col='grey')
boxplot(DE_withDD.lfc2.results.DEDD2$tpr[-grep('raw',names(DE_withDD.lfc2.results.DEDD2$tpr))], 
        names=1:14, ylim=c(0,0.05)); abline(h=c(0.05,0.1), col='grey')
boxplot(DE_withDD.lfc2.results.DEDD5$tpr[-grep('raw',names(DE_withDD.lfc2.results.DEDD5$tpr))], 
        names=1:12, ylim=c(0,0.05)); abline(h=c(0.05,0.1), col='grey')
boxplot(DE_withDD.lfc2.results.DEDD10$tpr[-grep('raw',names(DE_withDD.lfc2.results.DEDD10$tpr))], 
        names=1:12, ylim=c(0,0.05)); abline(h=c(0.05,0.1), col='grey')
cbind(round(colMeans(DE_noDD.lfc2.results.DEDD2$tpr[-grep('raw',names(DE_noDD.lfc2.results.DEDD2$tpr))]),3),
      round(colMeans(DE_withDD.lfc2.results.DEDD2$tpr[-grep('raw',names(DE_withDD.lfc2.results.DEDD2$tpr))]),3))
cbind(round(colMeans(DE_noDD.lfc2.results.DEDD5$tpr[-grep('raw',names(DE_noDD.lfc2.results.DEDD5$tpr))]),3),
      round(colMeans(DE_withDD.lfc2.results.DEDD5$tpr[-grep('raw',names(DE_withDD.lfc2.results.DEDD5$tpr))]),3))
cbind(round(colMeans(DE_noDD.lfc2.results.DEDD10$tpr[-grep('raw',names(DE_noDD.lfc2.results.DEDD10$tpr))]),3),
      round(colMeans(DE_withDD.lfc2.results.DEDD10$tpr[-grep('raw',names(DE_withDD.lfc2.results.DEDD10$tpr))]),3))


# DD ####
# AUC
par(mfcol=c(3,2), mar=c(2.5,2.5,0.5,0.5))
boxplot(DD_noDE.results.DEDD2$auc, names=1:6, ylim=c(0.2,0.8))
abline(h=c(0.45,0.5,0.55), col='grey')
boxplot(DD_noDE.results.DEDD5$auc, names=1:6, ylim=c(0.2,0.8))
abline(h=c(0.5,0.55,0.6,0.65), col='grey')
boxplot(DD_noDE.results.DEDD10$auc, names=1:6, ylim=c(0.2,0.8))
abline(h=c(0.55,0.6,0.65,0.7), col='grey')
boxplot(DD_withDE.results.DEDD2$auc, names=1:6, ylim=c(0.2,0.8))
abline(h=c(0.45,0.5,0.55), col='grey')
boxplot(DD_withDE.results.DEDD5$auc, names=1:6, ylim=c(0.2,0.8))
abline(h=c(0.5,0.55,0.6,0.65), col='grey')
boxplot(DD_withDE.results.DEDD10$auc, names=1:6, ylim=c(0.2,0.8))
abline(h=c(0.55,0.6,0.65,0.7), col='grey')
cbind(round(colMeans(DD_noDE.results.DEDD2$auc),3),
      round(colMeans(DD_withDE.results.DEDD2$auc),3))
cbind(round(colMeans(DD_noDE.results.DEDD5$auc),3),
      round(colMeans(DD_withDE.results.DEDD5$auc),3))
cbind(round(colMeans(DD_noDE.results.DEDD10$auc),3),
      round(colMeans(DD_withDE.results.DEDD10$auc),3))

# FDR
par(mfcol=c(3,2), mar=c(2.5,2.5,0.5,0.5))
boxplot(DD_noDE.results.DEDD2$fdr[-grep('raw',names(DD_noDE.results.DEDD2$fdr))], 
        names=1:14, ylim=c(0,1)); abline(h=c(0.05,0.1), col='grey')
boxplot(DD_noDE.results.DEDD5$fdr[-grep('raw',names(DD_noDE.results.DEDD5$fdr))], 
        names=1:14, ylim=c(0,1)); abline(h=c(0.05,0.1), col='grey')
boxplot(DD_noDE.results.DEDD10$fdr[-grep('raw',names(DD_noDE.results.DEDD10$fdr))], 
        names=1:14, ylim=c(0,1)); abline(h=c(0.05,0.1), col='grey')
boxplot(DD_withDE.results.DEDD2$fdr[-grep('raw',names(DD_withDE.results.DEDD2$fdr))], 
        names=1:14, ylim=c(0,1)); abline(h=c(0.05,0.1), col='grey')
boxplot(DD_withDE.results.DEDD5$fdr[-grep('raw',names(DD_withDE.results.DEDD5$fdr))], 
        names=1:14, ylim=c(0,1)); abline(h=c(0.05,0.1), col='grey')
boxplot(DD_withDE.results.DEDD10$fdr[-grep('raw',names(DD_withDE.results.DEDD10$fdr))], 
        names=1:14, ylim=c(0,1)); abline(h=c(0.05,0.1), col='grey')
cbind(round(colMeans(DD_noDE.results.DEDD2$fdr[-grep('raw',names(DD_noDE.results.DEDD2$fdr))]),3),
      round(colMeans(DD_withDE.results.DEDD2$fdr[-grep('raw',names(DD_withDE.results.DEDD2$fdr))]),3))
cbind(round(colMeans(DD_noDE.results.DEDD5$fdr[-grep('raw',names(DD_noDE.results.DEDD5$fdr))]),3),
      round(colMeans(DD_withDE.results.DEDD5$fdr[-grep('raw',names(DD_withDE.results.DEDD5$fdr))]),3))
cbind(round(colMeans(DD_noDE.results.DEDD10$fdr[-grep('raw',names(DD_noDE.results.DEDD10$fdr))]),3),
      round(colMeans(DD_withDE.results.DEDD10$fdr[-grep('raw',names(DD_withDE.results.DEDD10$fdr))]),3))

# TPR
par(mfcol=c(3,2), mar=c(2.5,2.5,0.5,0.5))
boxplot(DD_noDE.results.DEDD2$tpr[-grep('raw',names(DD_noDE.results.DEDD2$tpr))], 
        names=1:14, ylim=c(0,0.1))
boxplot(DD_noDE.results.DEDD5$tpr[-grep('raw',names(DD_noDE.results.DEDD5$tpr))], 
        names=1:14, ylim=c(0,0.1))
boxplot(DD_noDE.results.DEDD10$tpr[-grep('raw',names(DD_noDE.results.DEDD10$tpr))], 
        names=1:14, ylim=c(0,0.1))
boxplot(DD_withDE.results.DEDD2$tpr[-grep('raw',names(DD_withDE.results.DEDD2$tpr))], 
        names=1:14, ylim=c(0,0.1))
boxplot(DD_withDE.results.DEDD5$tpr[-grep('raw',names(DD_withDE.results.DEDD5$tpr))], 
        names=1:14, ylim=c(0,0.1))
boxplot(DD_withDE.results.DEDD10$tpr[-grep('raw',names(DD_withDE.results.DEDD10$tpr))], 
        names=1:14, ylim=c(0,0.1))
cbind(round(colMeans(DD_noDE.results.DEDD2$tpr[-grep('raw',names(DD_noDE.results.DEDD2$tpr))]),3),
      round(colMeans(DD_withDE.results.DEDD2$tpr[-grep('raw',names(DD_withDE.results.DEDD2$tpr))]),3))
cbind(round(colMeans(DD_noDE.results.DEDD5$tpr[-grep('raw',names(DD_noDE.results.DEDD5$tpr))]),3),
      round(colMeans(DD_withDE.results.DEDD5$tpr[-grep('raw',names(DD_withDE.results.DEDD5$tpr))]),3))
cbind(round(colMeans(DD_noDE.results.DEDD10$tpr[-grep('raw',names(DD_noDE.results.DEDD10$tpr))]),3),
      round(colMeans(DD_withDE.results.DEDD10$tpr[-grep('raw',names(DD_withDE.results.DEDD10$tpr))]),3))


# DE lfc1 ####
# AUC
par(mfcol=c(3,2), mar=c(2.5,2.5,0.5,0.5))
boxplot(DD_noDE.lfc1.results.DEDD2$auc, ylim=c(0.2,0.8)); abline(h=c(0.45,0.5,0.55), col='grey')
boxplot(DD_noDE.lfc1.results.DEDD5$auc, ylim=c(0.2,0.8)); abline(h=c(0.5,0.55,0.6), col='grey')
boxplot(DD_noDE.lfc1.results.DEDD10$auc, ylim=c(0.2,0.8)); abline(h=c(0.55,0.6,0.65,0.7), col='grey')
boxplot(DD_withDE.lfc1.results.DEDD2$auc, ylim=c(0.2,0.8)); abline(h=c(0.45,0.5,0.55), col='grey')
boxplot(DD_withDE.lfc1.results.DEDD5$auc, ylim=c(0.2,0.8)); abline(h=c(0.5,0.55,0.6), col='grey')
boxplot(DD_withDE.lfc1.results.DEDD10$auc, ylim=c(0.2,0.8)); abline(h=c(0.55,0.6,0.65,0.7), col='grey')
cbind(round(colMeans(DD_noDE.lfc1.results.DEDD2$auc),3),
      round(colMeans(DD_withDE.lfc1.results.DEDD2$auc),3))
cbind(round(colMeans(DD_noDE.lfc1.results.DEDD5$auc),3),
      round(colMeans(DD_withDE.lfc1.results.DEDD5$auc),3))
cbind(round(colMeans(DD_noDE.lfc1.results.DEDD10$auc),3),
      round(colMeans(DD_withDE.lfc1.results.DEDD10$auc),3))

# FDR
par(mfcol=c(3,2), mar=c(2.5,2.5,0.5,0.5))
boxplot(DD_noDE.lfc1.results.DEDD2$fdr[-grep('raw',names(DD_noDE.lfc1.results.DEDD2$fdr))], 
        names=1:8, ylim=c(0,1)); abline(h=c(0.05,0.1), col='grey')
boxplot(DD_noDE.lfc1.results.DEDD5$fdr[-grep('raw',names(DD_noDE.lfc1.results.DEDD5$fdr))], 
        names=1:8, ylim=c(0,1)); abline(h=c(0.05,0.1), col='grey')
boxplot(DD_noDE.lfc1.results.DEDD10$fdr[-grep('raw',names(DD_noDE.lfc1.results.DEDD10$fdr))], 
        names=1:8, ylim=c(0,1)); abline(h=c(0.05,0.1), col='grey')
boxplot(DD_withDE.lfc1.results.DEDD2$fdr[-grep('raw',names(DD_withDE.lfc1.results.DEDD2$fdr))], 
        names=1:8, ylim=c(0,1)); abline(h=c(0.05,0.1), col='grey')
boxplot(DD_withDE.lfc1.results.DEDD5$fdr[-grep('raw',names(DD_withDE.lfc1.results.DEDD5$fdr))], 
        names=1:8, ylim=c(0,1)); abline(h=c(0.05,0.1), col='grey')
boxplot(DD_withDE.lfc1.results.DEDD10$fdr[-grep('raw',names(DD_withDE.lfc1.results.DEDD10$fdr))], 
        names=1:8, ylim=c(0,1)); abline(h=c(0.05,0.1), col='grey')
cbind(round(colMeans(DD_noDE.lfc1.results.DEDD2$fdr[-grep('raw',names(DD_noDE.lfc1.results.DEDD2$fdr))], na.rm=T),3),
      round(colMeans(DD_withDE.lfc1.results.DEDD2$fdr[-grep('raw',names(DD_withDE.lfc1.results.DEDD2$fdr))], na.rm=T),3))
cbind(round(colMeans(DD_noDE.lfc1.results.DEDD5$fdr[-grep('raw',names(DD_noDE.lfc1.results.DEDD5$fdr))], na.rm=T),3),
      round(colMeans(DD_withDE.lfc1.results.DEDD5$fdr[-grep('raw',names(DD_withDE.lfc1.results.DEDD5$fdr))], na.rm=T),3))
cbind(round(colMeans(DD_noDE.lfc1.results.DEDD10$fdr[-grep('raw',names(DD_noDE.lfc1.results.DEDD10$fdr))], na.rm=T),3),
      round(colMeans(DD_withDE.lfc1.results.DEDD10$fdr[-grep('raw',names(DD_withDE.lfc1.results.DEDD10$fdr))], na.rm=T),3))

# TPR
par(mfcol=c(3,2), mar=c(2.5,2.5,0.5,0.5))
boxplot(DD_noDE.lfc1.results.DEDD2$tpr[-grep('raw',names(DD_noDE.lfc1.results.DEDD2$tpr))], 
        names=1:8, ylim=c(0,0.1))
boxplot(DD_noDE.lfc1.results.DEDD5$tpr[-grep('raw',names(DD_noDE.lfc1.results.DEDD5$tpr))], 
        names=1:8, ylim=c(0,0.1))
boxplot(DD_noDE.lfc1.results.DEDD10$tpr[-grep('raw',names(DD_noDE.lfc1.results.DEDD10$tpr))], 
        names=1:8, ylim=c(0,0.1))
boxplot(DD_withDE.lfc1.results.DEDD2$tpr[-grep('raw',names(DD_withDE.lfc1.results.DEDD2$tpr))], 
        names=1:8, ylim=c(0,0.1))
boxplot(DD_withDE.lfc1.results.DEDD5$tpr[-grep('raw',names(DD_withDE.lfc1.results.DEDD5$tpr))], 
        names=1:8, ylim=c(0,0.1))
boxplot(DD_withDE.lfc1.results.DEDD10$tpr[-grep('raw',names(DD_withDE.lfc1.results.DEDD10$tpr))], 
        names=1:8, ylim=c(0,0.1))
cbind(round(colMeans(DD_noDE.lfc1.results.DEDD2$tpr[-grep('raw',names(DD_noDE.lfc1.results.DEDD2$tpr))]),3),
      round(colMeans(DD_withDE.lfc1.results.DEDD2$tpr[-grep('raw',names(DD_withDE.lfc1.results.DEDD2$tpr))]),3))
cbind(round(colMeans(DD_noDE.lfc1.results.DEDD5$tpr[-grep('raw',names(DD_noDE.lfc1.results.DEDD5$tpr))]),3),
      round(colMeans(DD_withDE.lfc1.results.DEDD5$tpr[-grep('raw',names(DD_withDE.lfc1.results.DEDD5$tpr))]),3))
cbind(round(colMeans(DD_noDE.lfc1.results.DEDD10$tpr[-grep('raw',names(DD_noDE.lfc1.results.DEDD10$tpr))]),3),
      round(colMeans(DD_withDE.lfc1.results.DEDD10$tpr[-grep('raw',names(DD_withDE.lfc1.results.DEDD10$tpr))]),3))


# DE lfc2 ####
# AUC
par(mfcol=c(3,2), mar=c(2.5,2.5,0.5,0.5))
boxplot(DD_noDE.lfc2.results.DEDD2$auc, ylim=c(0.2,0.8))
abline(h=c(0.45,0.5,0.55), col='grey')
boxplot(DD_noDE.lfc2.results.DEDD5$auc, ylim=c(0.2,0.8))
abline(h=c(0.5,0.55,0.6), col='grey')
boxplot(DD_noDE.lfc2.results.DEDD10$auc, ylim=c(0.2,0.8))
abline(h=c(0.55,0.5,0.65), col='grey')
boxplot(DD_withDE.lfc2.results.DEDD2$auc, ylim=c(0.2,0.8))
abline(h=c(0.45,0.5,0.55), col='grey')
boxplot(DD_withDE.lfc2.results.DEDD5$auc, ylim=c(0.2,0.8))
abline(h=c(0.5,0.55,0.6), col='grey')
boxplot(DD_withDE.lfc2.results.DEDD10$auc, ylim=c(0.2,0.8))
abline(h=c(0.55,0.5,0.65), col='grey')
cbind(round(colMeans(DD_noDE.lfc2.results.DEDD2$auc),3),
      round(colMeans(DD_withDE.lfc2.results.DEDD2$auc),3))
cbind(round(colMeans(DD_noDE.lfc2.results.DEDD5$auc),3),
      round(colMeans(DD_withDE.lfc2.results.DEDD5$auc),3))
cbind(round(colMeans(DD_noDE.lfc2.results.DEDD10$auc),3),
      round(colMeans(DD_withDE.lfc2.results.DEDD10$auc),3))

# FDR
par(mfcol=c(3,2), mar=c(2.5,2.5,0.5,0.5))
boxplot(DD_noDE.lfc2.results.DEDD2$fdr[-grep('raw',names(DD_noDE.lfc2.results.DEDD2$fdr))], 
        names=1:8, ylim=c(0,1)); abline(h=c(0.05,0.1), col='grey')
boxplot(DD_noDE.lfc2.results.DEDD5$fdr[-grep('raw',names(DD_noDE.lfc2.results.DEDD5$fdr))], 
        names=1:8, ylim=c(0,1)); abline(h=c(0.05,0.1), col='grey')
boxplot(DD_noDE.lfc2.results.DEDD10$fdr[-grep('raw',names(DD_noDE.lfc2.results.DEDD10$fdr))], 
        names=1:8, ylim=c(0,1)); abline(h=c(0.05,0.1), col='grey')
boxplot(DD_withDE.lfc2.results.DEDD2$fdr[-grep('raw',names(DD_withDE.lfc2.results.DEDD2$fdr))], 
        names=1:8, ylim=c(0,1)); abline(h=c(0.05,0.1), col='grey')
boxplot(DD_withDE.lfc2.results.DEDD5$fdr[-grep('raw',names(DD_withDE.lfc2.results.DEDD5$fdr))], 
        names=1:8, ylim=c(0,1)); abline(h=c(0.05,0.1), col='grey')
boxplot(DD_withDE.lfc2.results.DEDD10$fdr[-grep('raw',names(DD_withDE.lfc2.results.DEDD10$fdr))], 
        names=1:8, ylim=c(0,1)); abline(h=c(0.05,0.1), col='grey')
cbind(round(colMeans(DD_noDE.lfc2.results.DEDD2$fdr[-grep('raw',names(DD_noDE.lfc2.results.DEDD2$fdr))], na.rm=T),3),
      round(colMeans(DD_withDE.lfc2.results.DEDD2$fdr[-grep('raw',names(DD_withDE.lfc2.results.DEDD2$fdr))], na.rm=T),3))
cbind(round(colMeans(DD_noDE.lfc2.results.DEDD5$fdr[-grep('raw',names(DD_noDE.lfc2.results.DEDD5$fdr))], na.rm=T),3),
      round(colMeans(DD_withDE.lfc2.results.DEDD5$fdr[-grep('raw',names(DD_withDE.lfc2.results.DEDD5$fdr))], na.rm=T),3))
cbind(round(colMeans(DD_noDE.lfc2.results.DEDD10$fdr[-grep('raw',names(DD_noDE.lfc2.results.DEDD10$fdr))], na.rm=T),3),
      round(colMeans(DD_withDE.lfc2.results.DEDD10$fdr[-grep('raw',names(DD_withDE.lfc2.results.DEDD10$fdr))], na.rm=T),3))

# TPR
par(mfcol=c(3,2), mar=c(2.5,2.5,0.5,0.5))
boxplot(DD_noDE.lfc2.results.DEDD2$tpr[-grep('raw',names(DD_noDE.lfc2.results.DEDD2$tpr))], 
        names=1:8, ylim=c(0,0.05))
boxplot(DD_noDE.lfc2.results.DEDD5$tpr[-grep('raw',names(DD_noDE.lfc2.results.DEDD5$tpr))], 
        names=1:8, ylim=c(0,0.05))
boxplot(DD_noDE.lfc2.results.DEDD10$tpr[-grep('raw',names(DD_noDE.lfc2.results.DEDD10$tpr))], 
        names=1:8, ylim=c(0,0.05))
boxplot(DD_withDE.lfc2.results.DEDD2$tpr[-grep('raw',names(DD_withDE.lfc2.results.DEDD2$tpr))], 
        names=1:8, ylim=c(0,0.05))
boxplot(DD_withDE.lfc2.results.DEDD5$tpr[-grep('raw',names(DD_withDE.lfc2.results.DEDD5$tpr))], 
        names=1:8, ylim=c(0,0.05))
boxplot(DD_withDE.lfc2.results.DEDD10$tpr[-grep('raw',names(DD_withDE.lfc2.results.DEDD10$tpr))], 
        names=1:8, ylim=c(0,0.05))
cbind(round(colMeans(DD_noDE.lfc2.results.DEDD2$tpr[-grep('raw',names(DD_noDE.lfc2.results.DEDD2$tpr))]),3),
      round(colMeans(DD_withDE.lfc2.results.DEDD2$tpr[-grep('raw',names(DD_withDE.lfc2.results.DEDD2$tpr))]),3))
cbind(round(colMeans(DD_noDE.lfc2.results.DEDD5$tpr[-grep('raw',names(DD_noDE.lfc2.results.DEDD5$tpr))]),3),
      round(colMeans(DD_withDE.lfc2.results.DEDD5$tpr[-grep('raw',names(DD_withDE.lfc2.results.DEDD5$tpr))]),3))
cbind(round(colMeans(DD_noDE.lfc2.results.DEDD10$tpr[-grep('raw',names(DD_noDE.lfc2.results.DEDD10$tpr))]),3),
      round(colMeans(DD_withDE.lfc2.results.DEDD10$tpr[-grep('raw',names(DD_withDE.lfc2.results.DEDD10$tpr))]),3))


# DEDD ####
# AUC
par(mfcol=c(3,3), mar=c(2.5,2.5,0.5,0.5))
boxplot(DEDD_DEandDD.results.DEDD2$auc, ylim=c(0.4,1))
abline(h=c(0.5,0.6,0.7), col='grey')
boxplot(DEDD_DEandDD.results.DEDD5$auc, ylim=c(0.4,1))
abline(h=c(0.5,0.6,0.7,0.8), col='grey')
boxplot(DEDD_DEandDD.results.DEDD10$auc, ylim=c(0.4,1))
abline(h=c(0.6,0.7,0.8,0.9), col='grey')
boxplot(DEDD_DEnoDD.results.DEDD2$auc, ylim=c(0.4,1))
abline(h=c(0.5,0.6,0.7), col='grey')
boxplot(DEDD_DEnoDD.results.DEDD5$auc, ylim=c(0.4,1))
abline(h=c(0.5,0.6,0.7,0.8), col='grey')
boxplot(DEDD_DEnoDD.results.DEDD10$auc, ylim=c(0.4,1))
abline(h=c(0.6,0.7,0.8,0.9), col='grey')
boxplot(DEDD_DDnoDE.results.DEDD2$auc, ylim=c(0.4,1))
abline(h=c(0.5,0.6,0.7), col='grey')
boxplot(DEDD_DDnoDE.results.DEDD5$auc, ylim=c(0.4,1))
abline(h=c(0.5,0.6,0.7,0.8), col='grey')
boxplot(DEDD_DDnoDE.results.DEDD10$auc, ylim=c(0.4,1))
abline(h=c(0.6,0.7,0.8,0.9), col='grey')
cbind(round(colMeans(DEDD_DEandDD.results.DEDD2$auc),3), 
      round(colMeans(DEDD_DEnoDD.results.DEDD2$auc),3), 
      round(colMeans(DEDD_DDnoDE.results.DEDD2$auc),3))
cbind(round(colMeans(DEDD_DEandDD.results.DEDD5$auc),3), 
      round(colMeans(DEDD_DEnoDD.results.DEDD5$auc),3), 
      round(colMeans(DEDD_DDnoDE.results.DEDD5$auc),3))
cbind(round(colMeans(DEDD_DEandDD.results.DEDD10$auc),3), 
      round(colMeans(DEDD_DEnoDD.results.DEDD10$auc),3), 
      round(colMeans(DEDD_DDnoDE.results.DEDD10$auc),3))

# FDR
par(mfcol=c(3,3), mar=c(2.5,2.5,0.5,0.5))
boxplot(DEDD_DEandDD.results.DEDD2$fdr, ylim=c(0,1), names=1:6); abline(h=c(0.05,0.1), col='grey')
boxplot(DEDD_DEandDD.results.DEDD5$fdr, ylim=c(0,1), names=1:6); abline(h=c(0.05,0.1), col='grey')
boxplot(DEDD_DEandDD.results.DEDD10$fdr, ylim=c(0,1), names=1:6); abline(h=c(0.05,0.1), col='grey')
boxplot(DEDD_DEnoDD.results.DEDD2$fdr, ylim=c(0,1), names=1:6); abline(h=c(0.05,0.1), col='grey')
boxplot(DEDD_DEnoDD.results.DEDD5$fdr, ylim=c(0,1), names=1:6); abline(h=c(0.05,0.1), col='grey')
boxplot(DEDD_DEnoDD.results.DEDD10$fdr, ylim=c(0,1), names=1:6); abline(h=c(0.05,0.1), col='grey')
boxplot(DEDD_DDnoDE.results.DEDD2$fdr, ylim=c(0,1), names=1:6); abline(h=c(0.05,0.1), col='grey')
boxplot(DEDD_DDnoDE.results.DEDD5$fdr, ylim=c(0,1), names=1:6); abline(h=c(0.05,0.1), col='grey')
boxplot(DEDD_DDnoDE.results.DEDD10$fdr, ylim=c(0,1), names=1:6); abline(h=c(0.05,0.1), col='grey')
cbind(round(colMeans(DEDD_DEandDD.results.DEDD2$fdr, na.rm=T),3), 
      round(colMeans(DEDD_DEnoDD.results.DEDD2$fdr, na.rm=T),3), 
      round(colMeans(DEDD_DDnoDE.results.DEDD2$fdr, na.rm=T),3))
cbind(round(colMeans(DEDD_DEandDD.results.DEDD5$fdr, na.rm=T),3), 
      round(colMeans(DEDD_DEnoDD.results.DEDD5$fdr, na.rm=T),3), 
      round(colMeans(DEDD_DDnoDE.results.DEDD5$fdr, na.rm=T),3))
cbind(round(colMeans(DEDD_DEandDD.results.DEDD10$fdr, na.rm=T),3), 
      round(colMeans(DEDD_DEnoDD.results.DEDD10$fdr, na.rm=T),3), 
      round(colMeans(DEDD_DDnoDE.results.DEDD10$fdr, na.rm=T),3))

# TPR
par(mfcol=c(3,3), mar=c(2.5,2.5,0.5,0.5))
boxplot(DEDD_DEandDD.results.DEDD2$tpr, ylim=c(0,0.1), names=1:6)
boxplot(DEDD_DEandDD.results.DEDD5$tpr, ylim=c(0,0.15), names=1:6)
boxplot(DEDD_DEandDD.results.DEDD10$tpr, ylim=c(0,0.5), names=1:6)
boxplot(DEDD_DEnoDD.results.DEDD2$tpr, ylim=c(0,0.1), names=1:6)
boxplot(DEDD_DEnoDD.results.DEDD5$tpr, ylim=c(0,0.15), names=1:6)
boxplot(DEDD_DEnoDD.results.DEDD10$tpr, ylim=c(0,0.5), names=1:6)
boxplot(DEDD_DDnoDE.results.DEDD2$tpr, ylim=c(0,0.1), names=1:6)
boxplot(DEDD_DDnoDE.results.DEDD5$tpr, ylim=c(0,0.15), names=1:6)
boxplot(DEDD_DDnoDE.results.DEDD10$tpr, ylim=c(0,0.5), names=1:6)
cbind(round(colMeans(DEDD_DEandDD.results.DEDD2$tpr, na.rm=T),3), 
      round(colMeans(DEDD_DEnoDD.results.DEDD2$tpr, na.rm=T),3), 
      round(colMeans(DEDD_DDnoDE.results.DEDD2$tpr, na.rm=T),3))
cbind(round(colMeans(DEDD_DEandDD.results.DEDD5$tpr, na.rm=T),3), 
      round(colMeans(DEDD_DEnoDD.results.DEDD5$tpr, na.rm=T),3), 
      round(colMeans(DEDD_DDnoDE.results.DEDD5$tpr, na.rm=T),3))
cbind(round(colMeans(DEDD_DEandDD.results.DEDD10$tpr, na.rm=T),3), 
      round(colMeans(DEDD_DEnoDD.results.DEDD10$tpr, na.rm=T),3), 
      round(colMeans(DEDD_DDnoDE.results.DEDD10$tpr, na.rm=T),3))

