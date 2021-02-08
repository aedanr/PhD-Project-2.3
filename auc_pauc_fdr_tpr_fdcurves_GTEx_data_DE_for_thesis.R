library(here)
library(RColorBrewer)

## Load data, create colour vector ####
folder <- "Results/GTEx combined results Sept 2020"
for (i in c("2", "5", "10", "20", "50")) {
  assign(paste0("DE.results.DEDD", i), 
         readRDS(here(folder, paste0("DE.results.DEDD", i, ".rds"))))
}
rm(i,folder)

qual_col_pals = brewer.pal.info[brewer.pal.info$category == "qual" & 
                                  brewer.pal.info$colorblind == T,]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector <- col_vector[c(1:3,9,11,26)]
rm(qual_col_pals)


## Keep only one version of each method: edgeR QL, DESeq2 IF, MDSeq ZI, lnHM log 
## (long chain for muscle).
## Separate objects for each metric and for blood and muscle
for (metric in c("pauc", "auc", "fpr", "fdr", "tpr", "mean.fdr", "mean.discoveries")) {
  for (size in c("2", "5", "10", "20", "50")) {
    assign(
      paste0(metric, "_blood_", size), 
      get(paste0("DE.results.DEDD", size))[[metric]][
        , c("blood_eR.ql", "blood_DES.if", "blood_voom", "blood_baySeq", "blood_MD.zi", "blood_lnHM.log")
      ]
    )
  }
}
for (metric in c("pauc", "auc", "fpr", "fdr", "tpr", "mean.fdr", "mean.discoveries")) {
  for (size in c("2", "5", "10", "20", "50")) {
    assign(
      paste0(metric, "_muscle_", size), 
      get(paste0("DE.results.DEDD", size))[[metric]][
        , c("muscle_eR.ql", "muscle_DES.if", "muscle_voom", "muscle_baySeq", "muscle_MD.zi", "muscle_lnHM.log.lc")
      ]
    )
  }
}

positions = rep(1:30)
offsets <- c(rep(0, 6), rep(1, 6), rep(2, 6), rep(3, 6), rep(4, 6))

### AUC ####
par(mfrow=c(2,1), mar=c(2,2.5,2,0.5), mgp=c(3,0.7,0))

plot(NA, xlim=c(1,34), ylim=c(0.65,1), xaxt="n", yaxt="n", 
     main="AUC - blood data", cex.main=2)
boxplot(cbind(auc_blood_2, auc_blood_5, auc_blood_10, auc_blood_20, auc_blood_50), 
        col=col_vector[1:6], pch=20, cex.axis=2, xaxt='n', ylim=c(0,0.8), 
        at=positions + offsets, lwd=0.2, cex=0.5, 
        add=T)
axis(side=1, at=c(3.5,10.5,17.5,24.5,31.5), labels=c(2,5,10,20,50), tick=F, cex.axis=2)
legend("bottomright", fill=col_vector[1:6], bty='n', cex=1.5, ncol=1, 
       legend=c("edgeR", "DESeq2", "voom", "baySeq", "MDSeq", "HM"))

plot(NA, xlim=c(1,34), ylim=c(0.65,1), xaxt="n", yaxt="n", 
     main="AUC - muscle data", cex.main=2)
boxplot(cbind(auc_muscle_2, auc_muscle_5, auc_muscle_10, auc_muscle_20, auc_muscle_50), 
        col=col_vector[1:6], pch=20, cex.axis=2, xaxt='n', ylim=c(0,0.8), 
        at=positions + offsets, lwd=0.2, cex=0.5, 
        add=T)
axis(side=1, at=c(3.5,10.5,17.5,24.5,31.5), labels=c(2,5,10,20,50), tick=F, cex.axis=2)

rbind(colMeans(auc_blood_2), 
      colMeans(auc_blood_5), 
      colMeans(auc_blood_10), 
      colMeans(auc_blood_20), 
      colMeans(auc_blood_50))
rbind(colMeans(auc_muscle_2), 
      colMeans(auc_muscle_5), 
      colMeans(auc_muscle_10), 
      colMeans(auc_muscle_20), 
      colMeans(auc_muscle_50))


### pAUC ####
par(mfrow=c(2,1), mar=c(2,2.5,2,0.5), mgp=c(3,0.7,0))

plot(NA, xlim=c(1,34), ylim=c(0,0.05), xaxt="n", yaxt="n", 
     main="Partial AUC - blood data", cex.main=2)
boxplot(cbind(pauc_blood_2, pauc_blood_5, pauc_blood_10, pauc_blood_20, pauc_blood_50), 
        col=col_vector[1:6], pch=20, cex.axis=2, xaxt='n', ylim=c(0,0.8), 
        at=positions + offsets, lwd=0.2, cex=0.5, 
        add=T)
axis(side=1, at=c(3.5,10.5,17.5,24.5,31.5), labels=c(2,5,10,20,50), tick=F, cex.axis=2)
legend("bottomright", fill=col_vector[1:6], bty='n', cex=1.5, ncol=1, 
       legend=c("edgeR", "DESeq2", "voom", "baySeq", "MDSeq", "HM"))

plot(NA, xlim=c(1,34), ylim=c(0,0.05), xaxt="n", yaxt="n", 
     main="Partial AUC - muscle data", cex.main=2)
boxplot(cbind(pauc_muscle_2, pauc_muscle_5, pauc_muscle_10, pauc_muscle_20, pauc_muscle_50), 
        col=col_vector[1:6], pch=20, cex.axis=2, xaxt='n', ylim=c(0,0.8), 
        at=positions + offsets, lwd=0.2, cex=0.5, 
        add=T)
axis(side=1, at=c(3.5,10.5,17.5,24.5,31.5), labels=c(2,5,10,20,50), tick=F, cex.axis=2)

rbind(colMeans(pauc_blood_2), 
      colMeans(pauc_blood_5), 
      colMeans(pauc_blood_10), 
      colMeans(pauc_blood_20), 
      colMeans(pauc_blood_50))
rbind(colMeans(pauc_muscle_2), 
      colMeans(pauc_muscle_5), 
      colMeans(pauc_muscle_10), 
      colMeans(pauc_muscle_20), 
      colMeans(pauc_muscle_50))


### FDR ####
par(mfrow=c(2,1), mar=c(2,2.5,2,0.5), mgp=c(3,0.7,0))

plot(NA, xlim=c(1,34), ylim=c(0,1), xaxt="n", yaxt="n", 
     main="FDR - blood data", cex.main=2)
lines(c(0.5,6.5), c(0.05,0.05), col="lightgrey")
lines(c(7.5,13.5), c(0.05,0.05), col="lightgrey")
lines(c(14.5,20.5), c(0.05,0.05), col="lightgrey")
lines(c(21.5,27.5), c(0.05,0.05), col="lightgrey")
lines(c(28.5,34.5), c(0.05,0.05), col="lightgrey")
boxplot(cbind(fdr_blood_2, fdr_blood_5, fdr_blood_10, fdr_blood_20, fdr_blood_50), 
        col=col_vector[1:6], pch=20, cex.axis=2, xaxt='n', ylim=c(0,0.8), 
        at=positions + offsets, lwd=0.2, cex=0.5, 
        add=T)
axis(side=1, at=c(3.5,10.5,17.5,24.5,31.5), labels=c(2,5,10,20,50), tick=F, cex.axis=2)

plot(NA, xlim=c(1,34), ylim=c(0,1), xaxt="n", yaxt="n", 
     main="FDR - muscle data", cex.main=2)
lines(c(0.5,6.5), c(0.05,0.05), col="lightgrey")
lines(c(7.5,13.5), c(0.05,0.05), col="lightgrey")
lines(c(14.5,20.5), c(0.05,0.05), col="lightgrey")
lines(c(21.5,27.5), c(0.05,0.05), col="lightgrey")
lines(c(28.5,34.5), c(0.05,0.05), col="lightgrey")
boxplot(cbind(fdr_muscle_2, fdr_muscle_5, fdr_muscle_10, fdr_muscle_20, fdr_muscle_50), 
        col=col_vector[1:6], pch=20, cex.axis=2, xaxt='n', ylim=c(0,0.8), 
        at=positions + offsets, lwd=0.2, cex=0.5, 
        add=T)
axis(side=1, at=c(3.5,10.5,17.5,24.5,31.5), labels=c(2,5,10,20,50), tick=F, cex.axis=2)
legend("topright", fill=col_vector[1:6], bty='n', cex=1.5, ncol=1, 
       legend=c("edgeR", "DESeq2", "voom", "baySeq", "MDSeq", "HM"))

rbind(colMeans(fdr_blood_2), 
      colMeans(fdr_blood_5), 
      colMeans(fdr_blood_10), 
      colMeans(fdr_blood_20), 
      colMeans(fdr_blood_50))
rbind(colMeans(fdr_muscle_2), 
      colMeans(fdr_muscle_5), 
      colMeans(fdr_muscle_10), 
      colMeans(fdr_muscle_20), 
      colMeans(fdr_muscle_50))


### TPR ####
par(mfrow=c(2,1), mar=c(2,2.5,2,0.5), mgp=c(3,0.7,0))

plot(NA, xlim=c(1,34), ylim=c(0,1), xaxt="n", yaxt="n", 
     main="TPR - blood data", cex.main=2)
boxplot(cbind(tpr_blood_2, tpr_blood_5, tpr_blood_10, tpr_blood_20, tpr_blood_50), 
        col=col_vector[1:6], pch=20, cex.axis=2, xaxt='n', ylim=c(0,0.8), 
        at=positions + offsets, lwd=0.2, cex=0.5, 
        add=T)
axis(side=1, at=c(3.5,10.5,17.5,24.5,31.5), labels=c(2,5,10,20,50), tick=F, cex.axis=2)
legend("bottomright", fill=col_vector[1:6], bty='n', cex=1.5, ncol=1, 
       legend=c("edgeR", "DESeq2", "voom", "baySeq", "MDSeq", "HM"))

plot(NA, xlim=c(1,34), ylim=c(0,1), xaxt="n", yaxt="n", 
     main="TPR - muscle data", cex.main=2)
boxplot(cbind(tpr_muscle_2, tpr_muscle_5, tpr_muscle_10, tpr_muscle_20, tpr_muscle_50), 
        col=col_vector[1:6], pch=20, cex.axis=2, xaxt='n', ylim=c(0,0.8), 
        at=positions + offsets, lwd=0.2, cex=0.5, 
        add=T)
axis(side=1, at=c(3.5,10.5,17.5,24.5,31.5), labels=c(2,5,10,20,50), tick=F, cex.axis=2)

rbind(colMeans(tpr_blood_2), 
      colMeans(tpr_blood_5), 
      colMeans(tpr_blood_10), 
      colMeans(tpr_blood_20), 
      colMeans(tpr_blood_50))
rbind(colMeans(tpr_muscle_2), 
      colMeans(tpr_muscle_5), 
      colMeans(tpr_muscle_10), 
      colMeans(tpr_muscle_20), 
      colMeans(tpr_muscle_50))


### FDR curves ####
par(mfrow=c(2,5), mar=c(1.5,3.5,2.5,1), mgp=c(3,1,0))

plot(mean.discoveries_blood_2$blood_eR.ql, mean.fdr_blood_2$blood_eR.ql, 
     type='l', ylim=c(0,0.6), xlim=c(0,2000), lwd=2, col=col_vector[1], cex.axis=3, 
     xlab=NA, ylab=NA, xaxt='n', main=paste0("2"), cex.main=3)
axis(side=1, at=c(0,500,1000,1500,2000), tick=T, cex.axis=3, mgp=c(3,2,0))
text("Blood", x=100, y=0.1, cex=3, adj=0)
lines(mean.discoveries_blood_2$blood_DES.if, mean.fdr_blood_2$blood_DES.if, 
      lwd=2, col=col_vector[2], cex.axis=3)
lines(mean.discoveries_blood_2$blood_voom, mean.fdr_blood_2$blood_voom, 
      lwd=2, col=col_vector[3], cex.axis=3)
lines(mean.discoveries_blood_2$blood_baySeq, mean.fdr_blood_2$blood_baySeq, 
      lwd=2, col=col_vector[4], cex.axis=3)
lines(mean.discoveries_blood_2$blood_MD.zi, mean.fdr_blood_2$blood_MD.zi, 
      lwd=2, col=col_vector[5], cex.axis=3)
lines(mean.discoveries_blood_2$blood_lnHM.log, mean.fdr_blood_2$blood_lnHM.log, 
      lwd=2, col=col_vector[6], cex.axis=3)
abline(h=0.05, col='lightgrey')

plot(mean.discoveries_blood_5$blood_eR.ql, mean.fdr_blood_5$blood_eR.ql, 
     type='l', ylim=c(0,0.6), xlim=c(0,2000), lwd=2, col=col_vector[1], cex.axis=3, 
     xlab=NA, ylab=NA, xaxt='n', yaxt='n', main=paste0("5"), cex.main=3)
axis(side=1, at=c(0,500,1000,1500,2000), tick=T, cex.axis=3, mgp=c(3,2,0))
lines(mean.discoveries_blood_5$blood_DES.if, mean.fdr_blood_5$blood_DES.if, 
      lwd=2, col=col_vector[2], cex.axis=3)
lines(mean.discoveries_blood_5$blood_voom, mean.fdr_blood_5$blood_voom, 
      lwd=2, col=col_vector[3], cex.axis=3)
lines(mean.discoveries_blood_5$blood_baySeq, mean.fdr_blood_5$blood_baySeq, 
      lwd=2, col=col_vector[4], cex.axis=3)
lines(mean.discoveries_blood_5$blood_MD.zi, mean.fdr_blood_5$blood_MD.zi, 
      lwd=2, col=col_vector[5], cex.axis=3)
lines(mean.discoveries_blood_5$blood_lnHM.log, mean.fdr_blood_5$blood_lnHM.log, 
      lwd=2, col=col_vector[6], cex.axis=3)
abline(h=0.05, col='lightgrey')

plot(mean.discoveries_blood_10$blood_eR.ql, mean.fdr_blood_10$blood_eR.ql, 
     type='l', ylim=c(0,0.6), xlim=c(0,2000), lwd=2, col=col_vector[1], cex.axis=3, 
     xlab=NA, ylab=NA, xaxt='n', yaxt='n', main=paste0("10"), cex.main=3)
axis(side=1, at=c(0,500,1000,1500,2000), tick=T, cex.axis=3, mgp=c(3,2,0))
lines(mean.discoveries_blood_10$blood_DES.if, mean.fdr_blood_10$blood_DES.if, 
      lwd=2, col=col_vector[2], cex.axis=3)
lines(mean.discoveries_blood_10$blood_voom, mean.fdr_blood_10$blood_voom, 
      lwd=2, col=col_vector[3], cex.axis=3)
lines(mean.discoveries_blood_10$blood_baySeq, mean.fdr_blood_10$blood_baySeq, 
      lwd=2, col=col_vector[4], cex.axis=3)
lines(mean.discoveries_blood_10$blood_MD.zi, mean.fdr_blood_10$blood_MD.zi, 
      lwd=2, col=col_vector[5], cex.axis=3)
lines(mean.discoveries_blood_10$blood_lnHM.log, mean.fdr_blood_10$blood_lnHM.log, 
      lwd=2, col=col_vector[6], cex.axis=3)
abline(h=0.05, col='lightgrey')

plot(mean.discoveries_blood_20$blood_eR.ql, mean.fdr_blood_20$blood_eR.ql, 
     type='l', ylim=c(0,0.6), xlim=c(0,2000), lwd=2, col=col_vector[1], cex.axis=3, 
     xlab=NA, ylab=NA, xaxt='n', yaxt='n', main=paste0("20"), cex.main=3)
axis(side=1, at=c(0,500,1000,1500,2000), tick=T, cex.axis=3, mgp=c(3,2,0))
lines(mean.discoveries_blood_20$blood_DES.if, mean.fdr_blood_20$blood_DES.if, 
      lwd=2, col=col_vector[2], cex.axis=3)
lines(mean.discoveries_blood_20$blood_voom, mean.fdr_blood_20$blood_voom, 
      lwd=2, col=col_vector[3], cex.axis=3)
lines(mean.discoveries_blood_20$blood_baySeq, mean.fdr_blood_20$blood_baySeq, 
      lwd=2, col=col_vector[4], cex.axis=3)
lines(mean.discoveries_blood_20$blood_MD.zi, mean.fdr_blood_20$blood_MD.zi, 
      lwd=2, col=col_vector[5], cex.axis=3)
lines(mean.discoveries_blood_20$blood_lnHM.log, mean.fdr_blood_20$blood_lnHM.log, 
      lwd=2, col=col_vector[6], cex.axis=3)
abline(h=0.05, col='lightgrey')

plot(mean.discoveries_blood_50$blood_eR.ql, mean.fdr_blood_50$blood_eR.ql, 
     type='l', ylim=c(0,0.6), xlim=c(0,2000), lwd=2, col=col_vector[1], cex.axis=3, 
     xlab=NA, ylab=NA, xaxt='n', yaxt='n', main=paste0("50"), cex.main=3)
axis(side=1, at=c(0,500,1000,1500,2000), tick=T, cex.axis=3, mgp=c(3,2,0))
lines(mean.discoveries_blood_50$blood_DES.if, mean.fdr_blood_50$blood_DES.if, 
      lwd=2, col=col_vector[2], cex.axis=3)
lines(mean.discoveries_blood_50$blood_voom, mean.fdr_blood_50$blood_voom, 
      lwd=2, col=col_vector[3], cex.axis=3)
lines(mean.discoveries_blood_50$blood_baySeq, mean.fdr_blood_50$blood_baySeq, 
      lwd=2, col=col_vector[4], cex.axis=3)
lines(mean.discoveries_blood_50$blood_MD.zi, mean.fdr_blood_50$blood_MD.zi, 
      lwd=2, col=col_vector[5], cex.axis=3)
lines(mean.discoveries_blood_50$blood_lnHM.log, mean.fdr_blood_50$blood_lnHM.log, 
      lwd=2, col=col_vector[6], cex.axis=3)
abline(h=0.05, col='lightgrey')

plot(mean.discoveries_muscle_2$muscle_eR.ql, mean.fdr_muscle_2$muscle_eR.ql, 
     type='l', ylim=c(0,0.6), xlim=c(0,2000), lwd=2, col=col_vector[1], cex.axis=3, 
     xlab=NA, ylab=NA, xaxt='n')
text("Muscle", x=100, y=0.1, cex=3, adj=0)
lines(mean.discoveries_muscle_2$muscle_DES.if, mean.fdr_muscle_2$muscle_DES.if, 
      lwd=2, col=col_vector[2], cex.axis=3)
lines(mean.discoveries_muscle_2$muscle_voom, mean.fdr_muscle_2$muscle_voom, 
      lwd=2, col=col_vector[3], cex.axis=3)
lines(mean.discoveries_muscle_2$muscle_baySeq, mean.fdr_muscle_2$muscle_baySeq, 
      lwd=2, col=col_vector[4], cex.axis=3)
lines(mean.discoveries_muscle_2$muscle_MD.zi, mean.fdr_muscle_2$muscle_MD.zi, 
      lwd=2, col=col_vector[5], cex.axis=3)
lines(mean.discoveries_muscle_2$muscle_lnHM.log.lc, mean.fdr_muscle_2$muscle_lnHM.log.lc, 
      lwd=2, col=col_vector[6], cex.axis=3)
abline(h=0.05, col='lightgrey')

plot(mean.discoveries_muscle_5$muscle_eR.ql, mean.fdr_muscle_5$muscle_eR.ql, 
     type='l', ylim=c(0,0.6), xlim=c(0,2000), lwd=2, col=col_vector[1], cex.axis=3, 
     xlab=NA, ylab=NA, xaxt='n', yaxt='n')
lines(mean.discoveries_muscle_5$muscle_DES.if, mean.fdr_muscle_5$muscle_DES.if, 
      lwd=2, col=col_vector[2], cex.axis=3)
lines(mean.discoveries_muscle_5$muscle_voom, mean.fdr_muscle_5$muscle_voom, 
      lwd=2, col=col_vector[3], cex.axis=3)
lines(mean.discoveries_muscle_5$muscle_baySeq, mean.fdr_muscle_5$muscle_baySeq, 
      lwd=2, col=col_vector[4], cex.axis=3)
lines(mean.discoveries_muscle_5$muscle_MD.zi, mean.fdr_muscle_5$muscle_MD.zi, 
      lwd=2, col=col_vector[5], cex.axis=3)
lines(mean.discoveries_muscle_5$muscle_lnHM.log.lc, mean.fdr_muscle_5$muscle_lnHM.log.lc, 
      lwd=2, col=col_vector[6], cex.axis=3)
abline(h=0.05, col='lightgrey')

plot(mean.discoveries_muscle_10$muscle_eR.ql, mean.fdr_muscle_10$muscle_eR.ql, 
     type='l', ylim=c(0,0.6), xlim=c(0,2000), lwd=2, col=col_vector[1], cex.axis=3, 
     xlab=NA, ylab=NA, xaxt='n', yaxt='n')
lines(mean.discoveries_muscle_10$muscle_DES.if, mean.fdr_muscle_10$muscle_DES.if, 
      lwd=2, col=col_vector[2], cex.axis=3)
lines(mean.discoveries_muscle_10$muscle_voom, mean.fdr_muscle_10$muscle_voom, 
      lwd=2, col=col_vector[3], cex.axis=3)
lines(mean.discoveries_muscle_10$muscle_baySeq, mean.fdr_muscle_10$muscle_baySeq, 
      lwd=2, col=col_vector[4], cex.axis=3)
lines(mean.discoveries_muscle_10$muscle_MD.zi, mean.fdr_muscle_10$muscle_MD.zi, 
      lwd=2, col=col_vector[5], cex.axis=3)
lines(mean.discoveries_muscle_10$muscle_lnHM.log.lc, mean.fdr_muscle_10$muscle_lnHM.log.lc, 
      lwd=2, col=col_vector[6], cex.axis=3)
abline(h=0.05, col='lightgrey')

plot(mean.discoveries_muscle_20$muscle_eR.ql, mean.fdr_muscle_20$muscle_eR.ql, 
     type='l', ylim=c(0,0.6), xlim=c(0,2000), lwd=2, col=col_vector[1], cex.axis=3, 
     xlab=NA, ylab=NA, xaxt='n', yaxt='n')
lines(mean.discoveries_muscle_20$muscle_DES.if, mean.fdr_muscle_20$muscle_DES.if, 
      lwd=2, col=col_vector[2], cex.axis=3)
lines(mean.discoveries_muscle_20$muscle_voom, mean.fdr_muscle_20$muscle_voom, 
      lwd=2, col=col_vector[3], cex.axis=3)
lines(mean.discoveries_muscle_20$muscle_baySeq, mean.fdr_muscle_20$muscle_baySeq, 
      lwd=2, col=col_vector[4], cex.axis=3)
lines(mean.discoveries_muscle_20$muscle_MD.zi, mean.fdr_muscle_20$muscle_MD.zi, 
      lwd=2, col=col_vector[5], cex.axis=3)
lines(mean.discoveries_muscle_20$muscle_lnHM.log.lc, mean.fdr_muscle_20$muscle_lnHM.log.lc, 
      lwd=2, col=col_vector[6], cex.axis=3)
abline(h=0.05, col='lightgrey')

plot(mean.discoveries_muscle_50$muscle_eR.ql, mean.fdr_muscle_50$muscle_eR.ql, 
     type='l', ylim=c(0,0.6), xlim=c(0,2000), lwd=2, col=col_vector[1], cex.axis=3, 
     xlab=NA, ylab=NA, xaxt='n', yaxt='n')
lines(mean.discoveries_muscle_50$muscle_DES.if, mean.fdr_muscle_50$muscle_DES.if, 
      lwd=2, col=col_vector[2], cex.axis=3)
lines(mean.discoveries_muscle_50$muscle_voom, mean.fdr_muscle_50$muscle_voom, 
      lwd=2, col=col_vector[3], cex.axis=3)
lines(mean.discoveries_muscle_50$muscle_baySeq, mean.fdr_muscle_50$muscle_baySeq, 
      lwd=2, col=col_vector[4], cex.axis=3)
lines(mean.discoveries_muscle_50$muscle_MD.zi, mean.fdr_muscle_50$muscle_MD.zi, 
      lwd=2, col=col_vector[5], cex.axis=3)
lines(mean.discoveries_muscle_50$muscle_lnHM.log.lc, mean.fdr_muscle_50$muscle_lnHM.log.lc, 
      lwd=2, col=col_vector[6], cex.axis=3)
abline(h=0.05, col='lightgrey')

legend("topleft", col=col_vector[1:6], bty='n', cex=2.5, ncol=1, lty=1, lwd=2, 
       legend=c("edgeR", "DESeq2", "voom", "baySeq", "MDSeq", "HM"))

