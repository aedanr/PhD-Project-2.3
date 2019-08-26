library(here)

###################
# DEDD2, 5, 10 ####
for (i in c('DE', 'DE.lfc1', 'DE.lfc2', 'DD', 'DD.lfc1', 'DD.lfc2', 'DEDD')) {
  for (j in c('DEDD2', 'DEDD5', 'DEDD10', 'DEDD20')) {
    assign(paste0(i,'.results.',j), readRDS(here('Results/DEDD compcodeR data results July-Aug 2019', 
                                                 paste0(i,'.results.',j,'.rds'))))
  }
}
rm('i','j')


# DE ####
par(mfrow=c(4,1), mar=c(3,3,1,1))
boxplot(DE.results.DEDD2$auc, names=1:17)
boxplot(DE.results.DEDD5$auc, names=1:14)
boxplot(DE.results.DEDD10$auc, names=1:14)
boxplot(DE.results.DEDD20$auc, names=1:14)
round(colMeans(DE.results.DEDD2$auc),3)
round(colMeans(DE.results.DEDD5$auc),3)
round(colMeans(DE.results.DEDD10$auc),3)
round(colMeans(DE.results.DEDD20$auc),3)

par(mfrow=c(6,3), mar=c(2.5,3,1,1), mgp=c(1.8,0.7,0))
for (i in 1:17) {
  plot(DE.results.DEDD2$mean.discoveries[which(DE.results.DEDD2$mean.discoveries[,i]<1000),i], 
       DE.results.DEDD2$mean.fdr[which(DE.results.DEDD2$mean.discoveries[,i]<1000),i], type='l', 
       xlab='', ylab=names(DE.results.DEDD2$mean.fdr)[i], ylim=c(0,0.8))
}
par(mfrow=c(6,3), mar=c(2.5,3,1,1), mgp=c(1.8,0.7,0))
for (i in 1:14) {
  plot(DE.results.DEDD5$mean.discoveries[which(DE.results.DEDD5$mean.discoveries[,i]<1000),i], 
       DE.results.DEDD5$mean.fdr[which(DE.results.DEDD5$mean.discoveries[,i]<1000),i], type='l', 
       xlab='', ylab=names(DE.results.DEDD5$mean.fdr)[i], ylim=c(0,0.2))
}
par(mfrow=c(6,3), mar=c(2.5,3,1,1), mgp=c(1.8,0.7,0))
for (i in 1:14) {
  plot(DE.results.DEDD10$mean.discoveries[which(DE.results.DEDD10$mean.discoveries[,i]<1000),i], 
       DE.results.DEDD10$mean.fdr[which(DE.results.DEDD10$mean.discoveries[,i]<1000),i], type='l', 
       xlab='', ylab=names(DE.results.DEDD10$mean.fdr)[i], ylim=c(0,0.05))
}
par(mfrow=c(6,3), mar=c(2.5,3,1,1), mgp=c(1.8,0.7,0))
for (i in 1:14) {
  plot(DE.results.DEDD20$mean.discoveries[which(DE.results.DEDD20$mean.discoveries[,i]<2000),i], 
       DE.results.DEDD20$mean.fdr[which(DE.results.DEDD20$mean.discoveries[,i]<2000),i], type='l', 
       xlab='', ylab=names(DE.results.DEDD20$mean.fdr)[i], ylim=c(0,0.05))
}

par(mfrow=c(4,1), mar=c(3,3,1,1))
boxplot(DE.results.DEDD2$fdr[-grep('raw',names(DE.results.DEDD2$fdr))], 
        names=1:32); abline(h=0.05)
boxplot(DE.results.DEDD5$fdr[-grep('raw',names(DE.results.DEDD5$fdr))], 
        names=1:25); abline(h=0.05)
boxplot(DE.results.DEDD10$fdr[-grep('raw',names(DE.results.DEDD10$fdr))], 
        names=1:25); abline(h=0.05)
boxplot(DE.results.DEDD20$fdr[-grep('raw',names(DE.results.DEDD20$fdr))], 
        names=1:25); abline(h=0.05)
round(colMeans(DE.results.DEDD2$fdr[-grep('raw',names(DE.results.DEDD2$fdr))], na.rm=T),3)
round(colMeans(DE.results.DEDD5$fdr[-grep('raw',names(DE.results.DEDD5$fdr))], na.rm=T),3)
round(colMeans(DE.results.DEDD10$fdr[-grep('raw',names(DE.results.DEDD10$fdr))], na.rm=T),3)
round(colMeans(DE.results.DEDD20$fdr[-grep('raw',names(DE.results.DEDD20$fdr))], na.rm=T),3)

par(mfrow=c(4,1), mar=c(3,3,1,1))
boxplot(DE.results.DEDD2$tpr[-grep('raw',names(DE.results.DEDD2$tpr))], 
        names=1:32)
boxplot(DE.results.DEDD5$tpr[-grep('raw',names(DE.results.DEDD5$tpr))], 
        names=1:25)
boxplot(DE.results.DEDD10$tpr[-grep('raw',names(DE.results.DEDD10$tpr))], 
        names=1:25)
boxplot(DE.results.DEDD20$tpr[-grep('raw',names(DE.results.DEDD20$tpr))], 
        names=1:25)
round(colMeans(DE.results.DEDD2$tpr[-grep('raw',names(DE.results.DEDD2$tpr))], na.rm=T),3)
round(colMeans(DE.results.DEDD5$tpr[-grep('raw',names(DE.results.DEDD5$tpr))], na.rm=T),3)
round(colMeans(DE.results.DEDD10$tpr[-grep('raw',names(DE.results.DEDD10$tpr))], na.rm=T),3)
round(colMeans(DE.results.DEDD20$tpr[-grep('raw',names(DE.results.DEDD20$tpr))], na.rm=T),3)


# DE lfc1 ####
par(mfrow=c(4,1), mar=c(3,3,1,1))
boxplot(DE.lfc1.results.DEDD2$auc, names=1:9)
boxplot(DE.lfc1.results.DEDD5$auc, names=1:8)
boxplot(DE.lfc1.results.DEDD10$auc, names=1:8)
boxplot(DE.lfc1.results.DEDD20$auc, names=1:8)
round(colMeans(DE.lfc1.results.DEDD2$auc),3)
round(colMeans(DE.lfc1.results.DEDD5$auc),3)
round(colMeans(DE.lfc1.results.DEDD10$auc),3)
round(colMeans(DE.lfc1.results.DEDD20$auc),3)

par(mfrow=c(4,1), mar=c(3,3,1,1))
boxplot(DE.lfc1.results.DEDD2$fdr[-grep('raw',names(DE.lfc1.results.DEDD2$fdr))], 
        names=1:14); abline(h=0.05)
boxplot(DE.lfc1.results.DEDD5$fdr[-grep('raw',names(DE.lfc1.results.DEDD5$fdr))], 
        names=1:12); abline(h=0.05)
boxplot(DE.lfc1.results.DEDD10$fdr[-grep('raw',names(DE.lfc1.results.DEDD10$fdr))], 
        names=1:12); abline(h=0.05)
boxplot(DE.lfc1.results.DEDD20$fdr[-grep('raw',names(DE.lfc1.results.DEDD20$fdr))], 
        names=1:12); abline(h=0.05)
round(colMeans(DE.lfc1.results.DEDD2$fdr[-grep('raw',names(DE.lfc1.results.DEDD2$fdr))], na.rm=T),3)
round(colMeans(DE.lfc1.results.DEDD5$fdr[-grep('raw',names(DE.lfc1.results.DEDD5$fdr))], na.rm=T),3)
round(colMeans(DE.lfc1.results.DEDD10$fdr[-grep('raw',names(DE.lfc1.results.DEDD10$fdr))], na.rm=T),3)
round(colMeans(DE.lfc1.results.DEDD20$fdr[-grep('raw',names(DE.lfc1.results.DEDD20$fdr))], na.rm=T),3)

par(mfrow=c(4,1), mar=c(3,3,1,1))
boxplot(DE.lfc1.results.DEDD2$tpr[-grep('raw',names(DE.lfc1.results.DEDD2$tpr))], 
        names=1:14)
boxplot(DE.lfc1.results.DEDD5$tpr[-grep('raw',names(DE.lfc1.results.DEDD5$tpr))], 
        names=1:12)
boxplot(DE.lfc1.results.DEDD10$tpr[-grep('raw',names(DE.lfc1.results.DEDD10$tpr))], 
        names=1:12)
boxplot(DE.lfc1.results.DEDD20$tpr[-grep('raw',names(DE.lfc1.results.DEDD20$tpr))], 
        names=1:12)
round(colMeans(DE.lfc1.results.DEDD2$tpr[-grep('raw',names(DE.lfc1.results.DEDD2$tpr))], na.rm=T),3)
round(colMeans(DE.lfc1.results.DEDD5$tpr[-grep('raw',names(DE.lfc1.results.DEDD5$tpr))], na.rm=T),3)
round(colMeans(DE.lfc1.results.DEDD10$tpr[-grep('raw',names(DE.lfc1.results.DEDD10$tpr))], na.rm=T),3)
round(colMeans(DE.lfc1.results.DEDD20$tpr[-grep('raw',names(DE.lfc1.results.DEDD20$tpr))], na.rm=T),3)


# DE lfc2 ####
par(mfrow=c(4,1), mar=c(3,3,1,1))
boxplot(DE.lfc2.results.DEDD2$auc, names=1:9)
boxplot(DE.lfc2.results.DEDD5$auc, names=1:8)
boxplot(DE.lfc2.results.DEDD10$auc, names=1:8)
boxplot(DE.lfc2.results.DEDD20$auc, names=1:8)
round(colMeans(DE.lfc2.results.DEDD2$auc),3)
round(colMeans(DE.lfc2.results.DEDD5$auc),3)
round(colMeans(DE.lfc2.results.DEDD10$auc),3)
round(colMeans(DE.lfc2.results.DEDD20$auc),3)

par(mfrow=c(4,1), mar=c(3,3,1,1))
boxplot(DE.lfc2.results.DEDD2$fdr[-grep('raw',names(DE.lfc2.results.DEDD2$fdr))], 
        names=1:14); abline(h=0.05)
boxplot(DE.lfc2.results.DEDD5$fdr[-grep('raw',names(DE.lfc2.results.DEDD5$fdr))], 
        names=1:12); abline(h=0.05)
boxplot(DE.lfc2.results.DEDD10$fdr[-grep('raw',names(DE.lfc2.results.DEDD10$fdr))], 
        names=1:12); abline(h=0.05)
boxplot(DE.lfc2.results.DEDD20$fdr[-grep('raw',names(DE.lfc2.results.DEDD20$fdr))], 
        names=1:12); abline(h=0.05)
round(colMeans(DE.lfc2.results.DEDD2$fdr[-grep('raw',names(DE.lfc2.results.DEDD2$fdr))], na.rm=T),3)
round(colMeans(DE.lfc2.results.DEDD5$fdr[-grep('raw',names(DE.lfc2.results.DEDD5$fdr))], na.rm=T),3)
round(colMeans(DE.lfc2.results.DEDD10$fdr[-grep('raw',names(DE.lfc2.results.DEDD10$fdr))], na.rm=T),3)
round(colMeans(DE.lfc2.results.DEDD20$fdr[-grep('raw',names(DE.lfc2.results.DEDD20$fdr))], na.rm=T),3)

par(mfrow=c(4,1), mar=c(3,3,1,1))
boxplot(DE.lfc2.results.DEDD2$tpr[-grep('raw',names(DE.lfc2.results.DEDD2$tpr))], 
        names=1:14)
boxplot(DE.lfc2.results.DEDD5$tpr[-grep('raw',names(DE.lfc2.results.DEDD5$tpr))], 
        names=1:12)
boxplot(DE.lfc2.results.DEDD10$tpr[-grep('raw',names(DE.lfc2.results.DEDD10$tpr))], 
        names=1:12)
boxplot(DE.lfc2.results.DEDD20$tpr[-grep('raw',names(DE.lfc2.results.DEDD20$tpr))], 
        names=1:12)
round(colMeans(DE.lfc2.results.DEDD2$tpr[-grep('raw',names(DE.lfc2.results.DEDD2$tpr))], na.rm=T),3)
round(colMeans(DE.lfc2.results.DEDD5$tpr[-grep('raw',names(DE.lfc2.results.DEDD5$tpr))], na.rm=T),3)
round(colMeans(DE.lfc2.results.DEDD10$tpr[-grep('raw',names(DE.lfc2.results.DEDD10$tpr))], na.rm=T),3)
round(colMeans(DE.lfc2.results.DEDD20$tpr[-grep('raw',names(DE.lfc2.results.DEDD20$tpr))], na.rm=T),3)


# DD ####
par(mfrow=c(4,1), mar=c(3,3,1,1))
boxplot(DD.results.DEDD2$auc)
boxplot(DD.results.DEDD5$auc)
boxplot(DD.results.DEDD10$auc)
boxplot(DD.results.DEDD20$auc)
round(colMeans(DD.results.DEDD2$auc),3)
round(colMeans(DD.results.DEDD5$auc),3)
round(colMeans(DD.results.DEDD10$auc),3)
round(colMeans(DD.results.DEDD20$auc),3)

par(mfrow=c(4,1), mar=c(3,3,1,1))
boxplot(DD.results.DEDD2$fdr[-grep('raw',names(DD.results.DEDD2$fdr))], 
        names=1:14); abline(h=0.05)
boxplot(DD.results.DEDD5$fdr[-grep('raw',names(DD.results.DEDD5$fdr))], 
        names=1:14); abline(h=0.05)
boxplot(DD.results.DEDD10$fdr[-grep('raw',names(DD.results.DEDD10$fdr))], 
        names=1:14); abline(h=0.05)
boxplot(DD.results.DEDD20$fdr[-grep('raw',names(DD.results.DEDD20$fdr))], 
        names=1:14); abline(h=0.05)
round(colMeans(DD.results.DEDD2$fdr[-grep('raw',names(DD.results.DEDD2$fdr))], na.rm=T),3)
round(colMeans(DD.results.DEDD5$fdr[-grep('raw',names(DD.results.DEDD5$fdr))], na.rm=T),3)
round(colMeans(DD.results.DEDD10$fdr[-grep('raw',names(DD.results.DEDD10$fdr))], na.rm=T),3)
round(colMeans(DD.results.DEDD20$fdr[-grep('raw',names(DD.results.DEDD20$fdr))], na.rm=T),3)

par(mfrow=c(4,1), mar=c(3,3,1,1))
boxplot(DD.results.DEDD2$tpr[-grep('raw',names(DD.results.DEDD2$tpr))], 
        names=1:14)
boxplot(DD.results.DEDD5$tpr[-grep('raw',names(DD.results.DEDD5$tpr))], 
        names=1:14)
boxplot(DD.results.DEDD10$tpr[-grep('raw',names(DD.results.DEDD10$tpr))], 
        names=1:14)
boxplot(DD.results.DEDD20$tpr[-grep('raw',names(DD.results.DEDD20$tpr))], 
        names=1:14)
round(colMeans(DD.results.DEDD2$tpr[-grep('raw',names(DD.results.DEDD2$tpr))], na.rm=T),3)
round(colMeans(DD.results.DEDD5$tpr[-grep('raw',names(DD.results.DEDD5$tpr))], na.rm=T),3)
round(colMeans(DD.results.DEDD10$tpr[-grep('raw',names(DD.results.DEDD10$tpr))], na.rm=T),3)
round(colMeans(DD.results.DEDD20$tpr[-grep('raw',names(DD.results.DEDD20$tpr))], na.rm=T),3)


# DD lfc1 ####
par(mfrow=c(4,1), mar=c(3,3,1,1))
boxplot(DD.lfc1.results.DEDD2$auc)
boxplot(DD.lfc1.results.DEDD5$auc)
boxplot(DD.lfc1.results.DEDD10$auc)
boxplot(DD.lfc1.results.DEDD20$auc)
round(colMeans(DD.lfc1.results.DEDD2$auc),3)
round(colMeans(DD.lfc1.results.DEDD5$auc),3)
round(colMeans(DD.lfc1.results.DEDD10$auc),3)
round(colMeans(DD.lfc1.results.DEDD20$auc),3)

par(mfrow=c(4,1), mar=c(3,3,1,1))
boxplot(DD.lfc1.results.DEDD2$fdr[-grep('raw',names(DD.lfc1.results.DEDD2$fdr))], 
        names=1:8); abline(h=0.05)
boxplot(DD.lfc1.results.DEDD5$fdr[-grep('raw',names(DD.lfc1.results.DEDD5$fdr))], 
        names=1:8); abline(h=0.05)
boxplot(DD.lfc1.results.DEDD10$fdr[-grep('raw',names(DD.lfc1.results.DEDD10$fdr))], 
        names=1:8); abline(h=0.05)
boxplot(DD.lfc1.results.DEDD20$fdr[-grep('raw',names(DD.lfc1.results.DEDD20$fdr))], 
        names=1:8); abline(h=0.05)
round(colMeans(DD.lfc1.results.DEDD2$fdr[-grep('raw',names(DD.lfc1.results.DEDD2$fdr))], na.rm=T),3)
round(colMeans(DD.lfc1.results.DEDD5$fdr[-grep('raw',names(DD.lfc1.results.DEDD5$fdr))], na.rm=T),3)
round(colMeans(DD.lfc1.results.DEDD10$fdr[-grep('raw',names(DD.lfc1.results.DEDD10$fdr))], na.rm=T),3)
round(colMeans(DD.lfc1.results.DEDD20$fdr[-grep('raw',names(DD.lfc1.results.DEDD20$fdr))], na.rm=T),3)

par(mfrow=c(4,1), mar=c(3,3,1,1))
boxplot(DD.lfc1.results.DEDD2$tpr[-grep('raw',names(DD.lfc1.results.DEDD2$tpr))], 
        names=1:8)
boxplot(DD.lfc1.results.DEDD5$tpr[-grep('raw',names(DD.lfc1.results.DEDD5$tpr))], 
        names=1:8)
boxplot(DD.lfc1.results.DEDD10$tpr[-grep('raw',names(DD.lfc1.results.DEDD10$tpr))], 
        names=1:8)
boxplot(DD.lfc1.results.DEDD20$tpr[-grep('raw',names(DD.lfc1.results.DEDD20$tpr))], 
        names=1:8)
round(colMeans(DD.lfc1.results.DEDD2$tpr[-grep('raw',names(DD.lfc1.results.DEDD2$tpr))], na.rm=T),3)
round(colMeans(DD.lfc1.results.DEDD5$tpr[-grep('raw',names(DD.lfc1.results.DEDD5$tpr))], na.rm=T),3)
round(colMeans(DD.lfc1.results.DEDD10$tpr[-grep('raw',names(DD.lfc1.results.DEDD10$tpr))], na.rm=T),3)
round(colMeans(DD.lfc1.results.DEDD20$tpr[-grep('raw',names(DD.lfc1.results.DEDD20$tpr))], na.rm=T),3)


# DD lfc2 ####
par(mfrow=c(4,1), mar=c(3,3,1,1))
boxplot(DD.lfc2.results.DEDD2$auc)
boxplot(DD.lfc2.results.DEDD5$auc)
boxplot(DD.lfc2.results.DEDD10$auc)
boxplot(DD.lfc2.results.DEDD20$auc)
round(colMeans(DD.lfc2.results.DEDD2$auc),3)
round(colMeans(DD.lfc2.results.DEDD5$auc),3)
round(colMeans(DD.lfc2.results.DEDD10$auc),3)
round(colMeans(DD.lfc2.results.DEDD20$auc),3)

par(mfrow=c(4,1), mar=c(3,3,1,1))
boxplot(DD.lfc2.results.DEDD2$fdr[-grep('raw',names(DD.lfc2.results.DEDD2$fdr))], 
        names=1:8); abline(h=0.05)
boxplot(DD.lfc2.results.DEDD5$fdr[-grep('raw',names(DD.lfc2.results.DEDD5$fdr))], 
        names=1:8); abline(h=0.05)
boxplot(DD.lfc2.results.DEDD10$fdr[-grep('raw',names(DD.lfc2.results.DEDD10$fdr))], 
        names=1:8); abline(h=0.05)
boxplot(DD.lfc2.results.DEDD20$fdr[-grep('raw',names(DD.lfc2.results.DEDD20$fdr))], 
        names=1:8); abline(h=0.05)
round(colMeans(DD.lfc2.results.DEDD2$fdr[-grep('raw',names(DD.lfc2.results.DEDD2$fdr))], na.rm=T),3)
round(colMeans(DD.lfc2.results.DEDD5$fdr[-grep('raw',names(DD.lfc2.results.DEDD5$fdr))], na.rm=T),3)
round(colMeans(DD.lfc2.results.DEDD10$fdr[-grep('raw',names(DD.lfc2.results.DEDD10$fdr))], na.rm=T),3)
round(colMeans(DD.lfc2.results.DEDD20$fdr[-grep('raw',names(DD.lfc2.results.DEDD20$fdr))], na.rm=T),3)

par(mfrow=c(4,1), mar=c(3,3,1,1))
boxplot(DD.lfc2.results.DEDD2$tpr[-grep('raw',names(DD.lfc2.results.DEDD2$tpr))], 
        names=1:8)
boxplot(DD.lfc2.results.DEDD5$tpr[-grep('raw',names(DD.lfc2.results.DEDD5$tpr))], 
        names=1:8)
boxplot(DD.lfc2.results.DEDD10$tpr[-grep('raw',names(DD.lfc2.results.DEDD10$tpr))], 
        names=1:8)
boxplot(DD.lfc2.results.DEDD20$tpr[-grep('raw',names(DD.lfc2.results.DEDD20$tpr))], 
        names=1:8)
round(colMeans(DD.lfc2.results.DEDD2$tpr[-grep('raw',names(DD.lfc2.results.DEDD2$tpr))], na.rm=T),3)
round(colMeans(DD.lfc2.results.DEDD5$tpr[-grep('raw',names(DD.lfc2.results.DEDD5$tpr))], na.rm=T),3)
round(colMeans(DD.lfc2.results.DEDD10$tpr[-grep('raw',names(DD.lfc2.results.DEDD10$tpr))], na.rm=T),3)
round(colMeans(DD.lfc2.results.DEDD20$tpr[-grep('raw',names(DD.lfc2.results.DEDD20$tpr))], na.rm=T),3)


# DEDD ####
par(mfrow=c(4,1), mar=c(3,3,1,1))
boxplot(DEDD.results.DEDD2$auc)
boxplot(DEDD.results.DEDD5$auc)
boxplot(DEDD.results.DEDD10$auc)
boxplot(DEDD.results.DEDD20$auc)
round(colMeans(DEDD.results.DEDD2$auc),3)
round(colMeans(DEDD.results.DEDD5$auc),3)
round(colMeans(DEDD.results.DEDD10$auc),3)
round(colMeans(DEDD.results.DEDD20$auc),3)

par(mfrow=c(4,1), mar=c(3,3,1,1))
boxplot(DEDD.results.DEDD2$fdr); abline(h=0.05)
boxplot(DEDD.results.DEDD5$fdr); abline(h=0.05)
boxplot(DEDD.results.DEDD10$fdr); abline(h=0.05)
boxplot(DEDD.results.DEDD20$fdr); abline(h=0.05)
round(colMeans(DEDD.results.DEDD2$fdr, na.rm=T),3)
round(colMeans(DEDD.results.DEDD5$fdr, na.rm=T),3)
round(colMeans(DEDD.results.DEDD10$fdr, na.rm=T),3)
round(colMeans(DEDD.results.DEDD20$fdr, na.rm=T),3)

par(mfrow=c(4,1), mar=c(3,3,1,1))
boxplot(DEDD.results.DEDD2$tpr)
boxplot(DEDD.results.DEDD5$tpr)
boxplot(DEDD.results.DEDD10$tpr)
boxplot(DEDD.results.DEDD20$tpr)
round(colMeans(DEDD.results.DEDD2$tpr, na.rm=T),3)
round(colMeans(DEDD.results.DEDD5$tpr, na.rm=T),3)
round(colMeans(DEDD.results.DEDD10$tpr, na.rm=T),3)
round(colMeans(DEDD.results.DEDD20$tpr, na.rm=T),3)


#################
# DE2, 5, 10 ####
for (i in c('DE', 'DE.lfc1', 'DE.lfc2', 'DEDD')) {
  for (j in c('DE2', 'DE5', 'DE10', 'DE20')) {
    assign(paste0(i,'.results.',j), readRDS(here('Results/DE compcodeR data results July-Aug 2019', 
                                                 paste0(i,'.results.',j,'.rds'))))
  }
}
rm('i','j')


# DE ####
par(mfrow=c(4,1), mar=c(3,3,1,1))
boxplot(DE.results.DE2$auc, names=1:15)
boxplot(DE.results.DE5$auc, names=1:15)
boxplot(DE.results.DE10$auc, names=1:14)
boxplot(DE.results.DE20$auc, names=1:14)
round(colMeans(DE.results.DE2$auc),3)
round(colMeans(DE.results.DE5$auc),3)
round(colMeans(DE.results.DE10$auc),3)
round(colMeans(DE.results.DE20$auc),3)

par(mfrow=c(5,3), mar=c(2.5,3,1,1), mgp=c(1.8,0.7,0))
for (i in 1:15) {
  plot(DE.results.DE2$mean.discoveries[which(DE.results.DE2$mean.discoveries[,i]<1000),i], 
       DE.results.DE2$mean.fdr[which(DE.results.DE2$mean.discoveries[,i]<1000),i], type='l', 
       xlab='', ylab=names(DE.results.DE2$mean.fdr)[i], ylim=c(0,0.8))
}
par(mfrow=c(5,3), mar=c(2.5,3,1,1), mgp=c(1.8,0.7,0))
for (i in 1:14) {
  plot(DE.results.DE5$mean.discoveries[which(DE.results.DE5$mean.discoveries[,i]<1000),i], 
       DE.results.DE5$mean.fdr[which(DE.results.DE5$mean.discoveries[,i]<1000),i], type='l', 
       xlab='', ylab=names(DE.results.DE5$mean.fdr)[i], ylim=c(0,0.2))
}
par(mfrow=c(5,3), mar=c(2.5,3,1,1), mgp=c(1.8,0.7,0))
for (i in 1:14) {
  plot(DE.results.DE10$mean.discoveries[which(DE.results.DE10$mean.discoveries[,i]<1000),i], 
       DE.results.DE10$mean.fdr[which(DE.results.DE10$mean.discoveries[,i]<1000),i], type='l', 
       xlab='', ylab=names(DE.results.DE10$mean.fdr)[i], ylim=c(0,0.05))
}
par(mfrow=c(5,3), mar=c(2.5,3,1,1), mgp=c(1.8,0.7,0))
for (i in 1:14) {
  plot(DE.results.DE20$mean.discoveries[which(DE.results.DE20$mean.discoveries[,i]<1000),i], 
       DE.results.DE20$mean.fdr[which(DE.results.DE20$mean.discoveries[,i]<1000),i], type='l', 
       xlab='', ylab=names(DE.results.DE10$mean.fdr)[i], ylim=c(0,0.05))
}

par(mfrow=c(4,1), mar=c(3,3,1,1))
boxplot(DE.results.DE2$fdr[-grep('raw',names(DE.results.DE2$fdr))], 
        names=1:28); abline(h=0.05)
boxplot(DE.results.DE5$fdr[-grep('raw',names(DE.results.DE5$fdr))], 
        names=1:28); abline(h=0.05)
boxplot(DE.results.DE10$fdr[-grep('raw',names(DE.results.DE10$fdr))], 
        names=1:25); abline(h=0.05)
boxplot(DE.results.DE20$fdr[-grep('raw',names(DE.results.DE20$fdr))], 
        names=1:25); abline(h=0.05)
round(colMeans(DE.results.DE2$fdr[-grep('raw',names(DE.results.DE2$fdr))], na.rm=T),3)
round(colMeans(DE.results.DE5$fdr[-grep('raw',names(DE.results.DE5$fdr))], na.rm=T),3)
round(colMeans(DE.results.DE10$fdr[-grep('raw',names(DE.results.DE10$fdr))], na.rm=T),3)
round(colMeans(DE.results.DE20$fdr[-grep('raw',names(DE.results.DE20$fdr))], na.rm=T),3)

par(mfrow=c(4,1), mar=c(3,3,1,1))
boxplot(DE.results.DE2$tpr[-grep('raw',names(DE.results.DE2$tpr))], 
        names=1:28)
boxplot(DE.results.DE5$tpr[-grep('raw',names(DE.results.DE5$tpr))], 
        names=1:28)
boxplot(DE.results.DE10$tpr[-grep('raw',names(DE.results.DE10$tpr))], 
        names=1:25)
boxplot(DE.results.DE20$tpr[-grep('raw',names(DE.results.DE20$tpr))], 
        names=1:25)
round(colMeans(DE.results.DE2$tpr[-grep('raw',names(DE.results.DE2$tpr))], na.rm=T),3)
round(colMeans(DE.results.DE5$tpr[-grep('raw',names(DE.results.DE5$tpr))], na.rm=T),3)
round(colMeans(DE.results.DE10$tpr[-grep('raw',names(DE.results.DE10$tpr))], na.rm=T),3)
round(colMeans(DE.results.DE20$tpr[-grep('raw',names(DE.results.DE20$tpr))], na.rm=T),3)


# DE lfc1 ####
par(mfrow=c(4,1), mar=c(3,3,1,1))
boxplot(DE.lfc1.results.DE2$auc, names=1:8)
boxplot(DE.lfc1.results.DE5$auc, names=1:8)
boxplot(DE.lfc1.results.DE10$auc, names=1:8)
boxplot(DE.lfc1.results.DE20$auc, names=1:8)
round(colMeans(DE.lfc1.results.DE2$auc),3)
round(colMeans(DE.lfc1.results.DE5$auc),3)
round(colMeans(DE.lfc1.results.DE10$auc),3)
round(colMeans(DE.lfc1.results.DE20$auc),3)

par(mfrow=c(4,1), mar=c(3,3,1,1))
boxplot(DE.lfc1.results.DE2$fdr[-grep('raw',names(DE.lfc1.results.DE2$fdr))], 
        names=1:12); abline(h=0.05)
boxplot(DE.lfc1.results.DE5$fdr[-grep('raw',names(DE.lfc1.results.DE5$fdr))], 
        names=1:12); abline(h=0.05)
boxplot(DE.lfc1.results.DE10$fdr[-grep('raw',names(DE.lfc1.results.DE10$fdr))], 
        names=1:12); abline(h=0.05)
boxplot(DE.lfc1.results.DE20$fdr[-grep('raw',names(DE.lfc1.results.DE20$fdr))], 
        names=1:12); abline(h=0.05)
round(colMeans(DE.lfc1.results.DE2$fdr[-grep('raw',names(DE.lfc1.results.DE2$fdr))], na.rm=T),3)
round(colMeans(DE.lfc1.results.DE5$fdr[-grep('raw',names(DE.lfc1.results.DE5$fdr))], na.rm=T),3)
round(colMeans(DE.lfc1.results.DE10$fdr[-grep('raw',names(DE.lfc1.results.DE10$fdr))], na.rm=T),3)
round(colMeans(DE.lfc1.results.DE20$fdr[-grep('raw',names(DE.lfc1.results.DE20$fdr))], na.rm=T),3)

par(mfrow=c(4,1), mar=c(3,3,1,1))
boxplot(DE.lfc1.results.DE2$tpr[-grep('raw',names(DE.lfc1.results.DE2$tpr))], 
        names=1:12)
boxplot(DE.lfc1.results.DE5$tpr[-grep('raw',names(DE.lfc1.results.DE5$tpr))], 
        names=1:12)
boxplot(DE.lfc1.results.DE10$tpr[-grep('raw',names(DE.lfc1.results.DE10$tpr))], 
        names=1:12)
boxplot(DE.lfc1.results.DE20$tpr[-grep('raw',names(DE.lfc1.results.DE20$tpr))], 
        names=1:12)
round(colMeans(DE.lfc1.results.DE2$tpr[-grep('raw',names(DE.lfc1.results.DE2$tpr))], na.rm=T),3)
round(colMeans(DE.lfc1.results.DE5$tpr[-grep('raw',names(DE.lfc1.results.DE5$tpr))], na.rm=T),3)
round(colMeans(DE.lfc1.results.DE10$tpr[-grep('raw',names(DE.lfc1.results.DE10$tpr))], na.rm=T),3)
round(colMeans(DE.lfc1.results.DE20$tpr[-grep('raw',names(DE.lfc1.results.DE20$tpr))], na.rm=T),3)


# DE lfc2 ####
par(mfrow=c(4,1), mar=c(3,3,1,1))
boxplot(DE.lfc2.results.DE2$auc, names=1:8)
boxplot(DE.lfc2.results.DE5$auc, names=1:8)
boxplot(DE.lfc2.results.DE10$auc, names=1:8)
boxplot(DE.lfc2.results.DE20$auc, names=1:8)
round(colMeans(DE.lfc2.results.DE2$auc),3)
round(colMeans(DE.lfc2.results.DE5$auc),3)
round(colMeans(DE.lfc2.results.DE10$auc),3)
round(colMeans(DE.lfc2.results.DE20$auc),3)

par(mfrow=c(4,1), mar=c(3,3,1,1))
boxplot(DE.lfc2.results.DE2$fdr[-grep('raw',names(DE.lfc2.results.DE2$fdr))], 
        names=1:12); abline(h=0.05)
boxplot(DE.lfc2.results.DE5$fdr[-grep('raw',names(DE.lfc2.results.DE5$fdr))], 
        names=1:12); abline(h=0.05)
boxplot(DE.lfc2.results.DE10$fdr[-grep('raw',names(DE.lfc2.results.DE10$fdr))], 
        names=1:12); abline(h=0.05)
boxplot(DE.lfc2.results.DE20$fdr[-grep('raw',names(DE.lfc2.results.DE20$fdr))], 
        names=1:12); abline(h=0.05)
round(colMeans(DE.lfc2.results.DE2$fdr[-grep('raw',names(DE.lfc2.results.DE2$fdr))], na.rm=T),3)
round(colMeans(DE.lfc2.results.DE5$fdr[-grep('raw',names(DE.lfc2.results.DE5$fdr))], na.rm=T),3)
round(colMeans(DE.lfc2.results.DE10$fdr[-grep('raw',names(DE.lfc2.results.DE10$fdr))], na.rm=T),3)
round(colMeans(DE.lfc2.results.DE20$fdr[-grep('raw',names(DE.lfc2.results.DE20$fdr))], na.rm=T),3)

par(mfrow=c(4,1), mar=c(3,3,1,1))
boxplot(DE.lfc2.results.DE2$tpr[-grep('raw',names(DE.lfc2.results.DE2$tpr))], 
        names=1:12)
boxplot(DE.lfc2.results.DE5$tpr[-grep('raw',names(DE.lfc2.results.DE5$tpr))], 
        names=1:12)
boxplot(DE.lfc2.results.DE10$tpr[-grep('raw',names(DE.lfc2.results.DE10$tpr))], 
        names=1:12)
boxplot(DE.lfc2.results.DE20$tpr[-grep('raw',names(DE.lfc2.results.DE20$tpr))], 
        names=1:12)
round(colMeans(DE.lfc2.results.DE2$tpr[-grep('raw',names(DE.lfc2.results.DE2$tpr))], na.rm=T),3)
round(colMeans(DE.lfc2.results.DE5$tpr[-grep('raw',names(DE.lfc2.results.DE5$tpr))], na.rm=T),3)
round(colMeans(DE.lfc2.results.DE10$tpr[-grep('raw',names(DE.lfc2.results.DE10$tpr))], na.rm=T),3)
round(colMeans(DE.lfc2.results.DE20$tpr[-grep('raw',names(DE.lfc2.results.DE20$tpr))], na.rm=T),3)


# DEDD ####
par(mfrow=c(4,1), mar=c(3,3,1,1))
boxplot(DEDD.results.DE2$auc)
boxplot(DEDD.results.DE5$auc)
boxplot(DEDD.results.DE10$auc)
boxplot(DEDD.results.DE20$auc)
round(colMeans(DEDD.results.DE2$auc),3)
round(colMeans(DEDD.results.DE5$auc),3)
round(colMeans(DEDD.results.DE10$auc),3)
round(colMeans(DEDD.results.DE20$auc),3)

par(mfrow=c(4,1), mar=c(3,3,1,1))
boxplot(DEDD.results.DE2$fdr); abline(h=0.05)
boxplot(DEDD.results.DE5$fdr); abline(h=0.05)
boxplot(DEDD.results.DE10$fdr); abline(h=0.05)
boxplot(DEDD.results.DE20$fdr); abline(h=0.05)
round(colMeans(DEDD.results.DE2$fdr, na.rm=T),3)
round(colMeans(DEDD.results.DE5$fdr, na.rm=T),3)
round(colMeans(DEDD.results.DE10$fdr, na.rm=T),3)
round(colMeans(DEDD.results.DE20$fdr, na.rm=T),3)

par(mfrow=c(4,1), mar=c(3,3,1,1))
boxplot(DEDD.results.DE2$tpr)
boxplot(DEDD.results.DE5$tpr)
boxplot(DEDD.results.DE10$tpr)
boxplot(DEDD.results.DE20$tpr)
round(colMeans(DEDD.results.DE2$tpr),3)
round(colMeans(DEDD.results.DE5$tpr),3)
round(colMeans(DEDD.results.DE10$tpr),3)
round(colMeans(DEDD.results.DE20$tpr),3)


#################
# DD2, 5, 10 ####
for (i in c('DD', 'DD.lfc1', 'DD.lfc2', 'DEDD')) {
  for (j in c('DD2', 'DD5', 'DD10', 'DD20')) {
    assign(paste0(i,'.results.',j), readRDS(here('Results/DD compcodeR data results July-Aug 2019', 
                                                 paste0(i,'.results.',j,'.rds'))))
  }
}
rm('i','j')


# DD ####
par(mfrow=c(4,1), mar=c(3,3,1,1))
boxplot(DD.results.DD2$auc)
boxplot(DD.results.DD5$auc)
boxplot(DD.results.DD10$auc)
boxplot(DD.results.DD20$auc)
round(colMeans(DD.results.DD2$auc),3)
round(colMeans(DD.results.DD5$auc),3)
round(colMeans(DD.results.DD10$auc),3)
round(colMeans(DD.results.DD20$auc),3)

par(mfcol=c(6,4), mar=c(2.5,3,1,1), mgp=c(1.8,0.7,0))
for (i in 1:6) {
  plot(DD.results.DD2$mean.discoveries[which(DD.results.DD2$mean.discoveries[,i]<1000),i], 
       DD.results.DD2$mean.fdr[which(DD.results.DD2$mean.discoveries[,i]<1000),i], type='l', 
       xlab='', ylab=names(DD.results.DD2$mean.fdr)[i], ylim=c(0,1))
}
for (i in 1:6) {
  plot(DD.results.DD5$mean.discoveries[which(DD.results.DD5$mean.discoveries[,i]<1000),i], 
       DD.results.DD5$mean.fdr[which(DD.results.DD5$mean.discoveries[,i]<1000),i], type='l', 
       xlab='', ylab=names(DD.results.DD5$mean.fdr)[i], ylim=c(0,1))
}
for (i in 1:6) {
  plot(DD.results.DD10$mean.discoveries[which(DD.results.DD10$mean.discoveries[,i]<1000),i], 
       DD.results.DD10$mean.fdr[which(DD.results.DD10$mean.discoveries[,i]<1000),i], type='l', 
       xlab='', ylab=names(DD.results.DD10$mean.fdr)[i], ylim=c(0,1))
}
for (i in 1:6) {
  plot(DD.results.DD20$mean.discoveries[which(DD.results.DD20$mean.discoveries[,i]<1000),i], 
       DD.results.DD20$mean.fdr[which(DD.results.DD20$mean.discoveries[,i]<1000),i], type='l', 
       xlab='', ylab=names(DD.results.DD10$mean.fdr)[i], ylim=c(0,1))
}

par(mfrow=c(4,1), mar=c(3,3,1,1))
boxplot(DD.results.DD2$fdr[-grep('raw',names(DD.results.DD2$fdr))], 
        names=1:14); abline(h=0.05)
boxplot(DD.results.DD5$fdr[-grep('raw',names(DD.results.DD5$fdr))], 
        names=1:14); abline(h=0.05)
boxplot(DD.results.DD10$fdr[-grep('raw',names(DD.results.DD10$fdr))], 
        names=1:14); abline(h=0.05)
boxplot(DD.results.DD20$fdr[-grep('raw',names(DD.results.DD20$fdr))], 
        names=1:14); abline(h=0.05)
round(colMeans(DD.results.DD2$fdr[-grep('raw',names(DD.results.DD2$fdr))], na.rm=T),3)
round(colMeans(DD.results.DD5$fdr[-grep('raw',names(DD.results.DD5$fdr))], na.rm=T),3)
round(colMeans(DD.results.DD10$fdr[-grep('raw',names(DD.results.DD10$fdr))], na.rm=T),3)
round(colMeans(DD.results.DD20$fdr[-grep('raw',names(DD.results.DD20$fdr))], na.rm=T),3)

par(mfrow=c(4,1), mar=c(3,3,1,1))
boxplot(DD.results.DD2$tpr[-grep('raw',names(DD.results.DD2$tpr))], 
        names=1:14)
boxplot(DD.results.DD5$tpr[-grep('raw',names(DD.results.DD5$tpr))], 
        names=1:14)
boxplot(DD.results.DD10$tpr[-grep('raw',names(DD.results.DD10$tpr))], 
        names=1:14)
boxplot(DD.results.DD20$tpr[-grep('raw',names(DD.results.DD20$tpr))], 
        names=1:14)
round(colMeans(DD.results.DD2$tpr[-grep('raw',names(DD.results.DD2$tpr))], na.rm=T),3)
round(colMeans(DD.results.DD5$tpr[-grep('raw',names(DD.results.DD5$tpr))], na.rm=T),3)
round(colMeans(DD.results.DD10$tpr[-grep('raw',names(DD.results.DD10$tpr))], na.rm=T),3)
round(colMeans(DD.results.DD20$tpr[-grep('raw',names(DD.results.DD20$tpr))], na.rm=T),3)


# DD lfc1 ####
par(mfrow=c(4,1), mar=c(3,3,1,1))
boxplot(DD.lfc1.results.DD2$auc)
boxplot(DD.lfc1.results.DD5$auc)
boxplot(DD.lfc1.results.DD10$auc)
boxplot(DD.lfc1.results.DD20$auc)
round(colMeans(DD.lfc1.results.DD2$auc),3)
round(colMeans(DD.lfc1.results.DD5$auc),3)
round(colMeans(DD.lfc1.results.DD10$auc),3)
round(colMeans(DD.lfc1.results.DD20$auc),3)

par(mfrow=c(4,1), mar=c(3,3,1,1))
boxplot(DD.lfc1.results.DD2$fdr[-grep('raw',names(DD.lfc1.results.DD2$fdr))], 
        names=1:8); abline(h=0.05)
boxplot(DD.lfc1.results.DD5$fdr[-grep('raw',names(DD.lfc1.results.DD5$fdr))], 
        names=1:8); abline(h=0.05)
boxplot(DD.lfc1.results.DD10$fdr[-grep('raw',names(DD.lfc1.results.DD10$fdr))], 
        names=1:8); abline(h=0.05)
boxplot(DD.lfc1.results.DD20$fdr[-grep('raw',names(DD.lfc1.results.DD20$fdr))], 
        names=1:8); abline(h=0.05)
round(colMeans(DD.lfc1.results.DD2$fdr[-grep('raw',names(DD.lfc1.results.DD2$fdr))], na.rm=T),3)
round(colMeans(DD.lfc1.results.DD5$fdr[-grep('raw',names(DD.lfc1.results.DD5$fdr))], na.rm=T),3)
round(colMeans(DD.lfc1.results.DD10$fdr[-grep('raw',names(DD.lfc1.results.DD10$fdr))], na.rm=T),3)
round(colMeans(DD.lfc1.results.DD20$fdr[-grep('raw',names(DD.lfc1.results.DD20$fdr))], na.rm=T),3)

par(mfrow=c(4,1), mar=c(3,3,1,1))
boxplot(DD.lfc1.results.DD2$tpr[-grep('raw',names(DD.lfc1.results.DD2$tpr))], 
        names=1:8)
boxplot(DD.lfc1.results.DD5$tpr[-grep('raw',names(DD.lfc1.results.DD5$tpr))], 
        names=1:8)
boxplot(DD.lfc1.results.DD10$tpr[-grep('raw',names(DD.lfc1.results.DD10$tpr))], 
        names=1:8)
boxplot(DD.lfc1.results.DD20$tpr[-grep('raw',names(DD.lfc1.results.DD20$tpr))], 
        names=1:8)
round(colMeans(DD.lfc1.results.DD2$tpr[-grep('raw',names(DD.lfc1.results.DD2$tpr))], na.rm=T),3)
round(colMeans(DD.lfc1.results.DD5$tpr[-grep('raw',names(DD.lfc1.results.DD5$tpr))], na.rm=T),3)
round(colMeans(DD.lfc1.results.DD10$tpr[-grep('raw',names(DD.lfc1.results.DD10$tpr))], na.rm=T),3)
round(colMeans(DD.lfc1.results.DD20$tpr[-grep('raw',names(DD.lfc1.results.DD20$tpr))], na.rm=T),3)


# DD lfc2 ####
par(mfrow=c(4,1), mar=c(3,3,1,1))
boxplot(DD.lfc2.results.DD2$auc)
boxplot(DD.lfc2.results.DD5$auc)
boxplot(DD.lfc2.results.DD10$auc)
boxplot(DD.lfc2.results.DD10$auc)
round(colMeans(DD.lfc2.results.DD2$auc),3)
round(colMeans(DD.lfc2.results.DD5$auc),3)
round(colMeans(DD.lfc2.results.DD10$auc),3)
round(colMeans(DD.lfc2.results.DD10$auc),3)

par(mfrow=c(4,1), mar=c(3,3,1,1))
boxplot(DD.lfc2.results.DD2$fdr[-grep('raw',names(DD.lfc2.results.DD2$fdr))], 
        names=1:8); abline(h=0.05)
boxplot(DD.lfc2.results.DD5$fdr[-grep('raw',names(DD.lfc2.results.DD5$fdr))], 
        names=1:8); abline(h=0.05)
boxplot(DD.lfc2.results.DD10$fdr[-grep('raw',names(DD.lfc2.results.DD10$fdr))], 
        names=1:8); abline(h=0.05)
boxplot(DD.lfc2.results.DD20$fdr[-grep('raw',names(DD.lfc2.results.DD20$fdr))], 
        names=1:8); abline(h=0.05)
round(colMeans(DD.lfc2.results.DD2$fdr[-grep('raw',names(DD.lfc2.results.DD2$fdr))], na.rm=T),3)
round(colMeans(DD.lfc2.results.DD5$fdr[-grep('raw',names(DD.lfc2.results.DD5$fdr))], na.rm=T),3)
round(colMeans(DD.lfc2.results.DD10$fdr[-grep('raw',names(DD.lfc2.results.DD10$fdr))], na.rm=T),3)
round(colMeans(DD.lfc2.results.DD20$fdr[-grep('raw',names(DD.lfc2.results.DD20$fdr))], na.rm=T),3)

par(mfrow=c(4,1), mar=c(3,3,1,1))
boxplot(DD.lfc2.results.DD2$tpr[-grep('raw',names(DD.lfc2.results.DD2$tpr))], 
        names=1:8)
boxplot(DD.lfc2.results.DD5$tpr[-grep('raw',names(DD.lfc2.results.DD5$tpr))], 
        names=1:8)
boxplot(DD.lfc2.results.DD10$tpr[-grep('raw',names(DD.lfc2.results.DD10$tpr))], 
        names=1:8)
boxplot(DD.lfc2.results.DD20$tpr[-grep('raw',names(DD.lfc2.results.DD20$tpr))], 
        names=1:8)
round(colMeans(DD.lfc2.results.DD2$tpr[-grep('raw',names(DD.lfc2.results.DD2$tpr))], na.rm=T),3)
round(colMeans(DD.lfc2.results.DD5$tpr[-grep('raw',names(DD.lfc2.results.DD5$tpr))], na.rm=T),3)
round(colMeans(DD.lfc2.results.DD10$tpr[-grep('raw',names(DD.lfc2.results.DD10$tpr))], na.rm=T),3)
round(colMeans(DD.lfc2.results.DD20$tpr[-grep('raw',names(DD.lfc2.results.DD20$tpr))], na.rm=T),3)


# DEDD ####
par(mfrow=c(4,1), mar=c(3,3,1,1))
boxplot(DEDD.results.DD2$auc)
boxplot(DEDD.results.DD5$auc)
boxplot(DEDD.results.DD10$auc)
boxplot(DEDD.results.DD20$auc)
round(colMeans(DEDD.results.DD2$auc),3)
round(colMeans(DEDD.results.DD5$auc),3)
round(colMeans(DEDD.results.DD10$auc),3)
round(colMeans(DEDD.results.DD20$auc),3)

par(mfrow=c(4,1), mar=c(3,3,1,1))
boxplot(DEDD.results.DD2$fdr); abline(h=0.05)
boxplot(DEDD.results.DD5$fdr); abline(h=0.05)
boxplot(DEDD.results.DD10$fdr); abline(h=0.05)
boxplot(DEDD.results.DD20$fdr); abline(h=0.05)
round(colMeans(DEDD.results.DD2$fdr, na.rm=T),3)
round(colMeans(DEDD.results.DD5$fdr, na.rm=T),3)
round(colMeans(DEDD.results.DD10$fdr, na.rm=T),3)
round(colMeans(DEDD.results.DD20$fdr, na.rm=T),3)

par(mfrow=c(4,1), mar=c(3,3,1,1))
boxplot(DEDD.results.DD2$tpr)
boxplot(DEDD.results.DD5$tpr)
boxplot(DEDD.results.DD10$tpr)
boxplot(DEDD.results.DD20$tpr)
round(colMeans(DEDD.results.DD2$tpr),3)
round(colMeans(DEDD.results.DD5$tpr),3)
round(colMeans(DEDD.results.DD10$tpr),3)
round(colMeans(DEDD.results.DD20$tpr),3)
