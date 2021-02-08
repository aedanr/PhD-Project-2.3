library(here)
library(RColorBrewer)

## Load data, create colour vector ####
folder <- "Results/GTEx combined results Sept 2020"
for (i in c("2", "5", "10", "20", "50")) {
  assign(paste0("DEDD.results.DEDD", i), 
         readRDS(here(folder, paste0("DEDD.results.DEDD", i, ".rds"))))
}
rm(i,folder)

qual_col_pals = brewer.pal.info[brewer.pal.info$category == "qual" & 
                                  brewer.pal.info$colorblind == T,]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector <- col_vector[c(1:3,9,11,26)]
rm(qual_col_pals)


## Keep only diffVar, voom/lnHM.log, lnHMM (posterior threshold and BFDR; long chain for muscle).
## Separate objects for each metric and for blood and muscle
for (metric in c("pauc", "auc", "mean.fdr", "mean.discoveries")) {
  for (size in c("2", "5", "10", "20", "50")) {
    assign(
      paste0(metric, "_blood_", size), 
      get(paste0("DEDD.results.DEDD", size))[[metric]][
        , c("blood_dV", "blood_voom_HM", "blood_lnHMM")
      ]
    )
  }
}
for (metric in c("fpr", "fdr", "tpr")) {
  for (size in c("2", "5", "10", "20", "50")) {
    assign(
      paste0(metric, "_blood_", size), 
      get(paste0("DEDD.results.DEDD", size))[[metric]][
        , c("blood_dV", "blood_voom_HM", "blood_lnHMM.thr", "blood_lnHMM.bfdr")
      ]
    )
  }
}
for (metric in c("pauc", "auc", "mean.fdr", "mean.discoveries")) {
  for (size in c("2", "5", "10", "20", "50")) {
    assign(
      paste0(metric, "_muscle_", size), 
      get(paste0("DEDD.results.DEDD", size))[[metric]][
        , c("muscle_dV", "muscle_voom_HM.lc", "muscle_lnHMM.lc")
      ]
    )
  }
}
for (metric in c("fpr", "fdr", "tpr")) {
  for (size in c("2", "5", "10", "20", "50")) {
    assign(
      paste0(metric, "_muscle_", size), 
      get(paste0("DEDD.results.DEDD", size))[[metric]][
        , c("muscle_dV", "muscle_voom_HM.lc", "muscle_lnHMM.thr.lc", "muscle_lnHMM.bfdr.lc")
      ]
    )
  }
}


### AUC ####
positions <- rep(1:15)
offsets <- c(rep(0, 3), rep(0.5, 3), rep(1, 3), rep(1.5, 3), rep(2, 3))

par(mfrow=c(2,1), mar=c(2,2.5,2,0.5), mgp=c(3,0.7,0))

plot(NA, xlim=c(1,17), ylim=c(0.42,0.98), xaxt="n", yaxt="n", 
     main="AUC - blood data", cex.main=2)
boxplot(cbind(auc_blood_2, auc_blood_5, auc_blood_10, auc_blood_20, auc_blood_50), 
        col=col_vector[1:3], pch=20, cex.axis=2, xaxt='n', ylim=c(0,0.8), 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(2,5.5,9,12.5,16), labels=c(2,5,10,20,50), tick=F, cex.axis=2)
legend("bottomright", fill=col_vector[1:6], bty='n', cex=1.8, ncol=1, 
       legend=c("diffVar", "Hybrid", "HMM"))

plot(NA, xlim=c(1,17), ylim=c(0.42,0.98), xaxt="n", yaxt="n", 
     main="AUC - muscle data", cex.main=2)
boxplot(cbind(auc_muscle_2, auc_muscle_5, auc_muscle_10, auc_muscle_20, auc_muscle_50), 
        col=col_vector[1:3], pch=20, cex.axis=2, xaxt='n', ylim=c(0,0.8), 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(2,5.5,9,12.5,16), labels=c(2,5,10,20,50), tick=F, cex.axis=2)

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
positions <- rep(1:15)
offsets <- c(rep(0, 3), rep(0.5, 3), rep(1, 3), rep(1.5, 3), rep(2, 3))

par(mfrow=c(2,1), mar=c(2,2.5,2,0.5), mgp=c(3,0.7,0))

plot(NA, xlim=c(1,17), ylim=c(0,0.044), xaxt="n", yaxt="n", 
     main="Partial AUC - blood data", cex.main=2)
boxplot(cbind(pauc_blood_2, pauc_blood_5, pauc_blood_10, pauc_blood_20, pauc_blood_50), 
        col=col_vector[1:3], pch=20, cex.axis=2, xaxt='n', ylim=c(0,0.8), 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(2,5.5,9,12.5,16), labels=c(2,5,10,20,50), tick=F, cex.axis=2)
legend("bottomright", fill=col_vector[1:6], bty='n', cex=1.8, ncol=1, 
       legend=c("diffVar", "Hybrid", "HMM"))

plot(NA, xlim=c(1,17), ylim=c(0,0.044), xaxt="n", yaxt="n", 
     main="Partial AUC - muscle data", cex.main=2)
boxplot(cbind(pauc_muscle_2, pauc_muscle_5, pauc_muscle_10, pauc_muscle_20, pauc_muscle_50), 
        col=col_vector[1:3], pch=20, cex.axis=2, xaxt='n', ylim=c(0,0.8), 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(2,5.5,9,12.5,16), labels=c(2,5,10,20,50), tick=F, cex.axis=2)


### FDR ####
positions = rep(1:20)
offsets <- c(rep(0, 4), rep(1, 4), rep(2, 4), rep(3, 4), rep(4, 4))

par(mfrow=c(2,1), mar=c(2,2.5,2,0.5), mgp=c(3,0.7,0))

plot(NA, xlim=c(1,24), ylim=c(0,1), xaxt="n", yaxt="n", 
     main="FDR - blood data", cex.main=2)
lines(c(0.5,4.5), c(0.05,0.05), col="lightgrey")
lines(c(5.5,9.5), c(0.05,0.05), col="lightgrey")
lines(c(10.5,14.5), c(0.05,0.05), col="lightgrey")
lines(c(15.5,19.5), c(0.05,0.05), col="lightgrey")
lines(c(20.5,24.5), c(0.05,0.05), col="lightgrey")
boxplot(cbind(fdr_blood_2, fdr_blood_5, fdr_blood_10, fdr_blood_20, fdr_blood_50), 
        col=col_vector[1:4], pch=20, cex.axis=2, xaxt='n', ylim=c(0,0.8), 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(2.5,7.5,12.5,17.5,22.5), labels=c(2,5,10,20,50), tick=F, cex.axis=2)

plot(NA, xlim=c(1,24), ylim=c(0,1), xaxt="n", yaxt="n", 
     main="FDR - muscle data", cex.main=2)
lines(c(0.5,4.5), c(0.05,0.05), col="lightgrey")
lines(c(5.5,9.5), c(0.05,0.05), col="lightgrey")
lines(c(10.5,14.5), c(0.05,0.05), col="lightgrey")
lines(c(15.5,19.5), c(0.05,0.05), col="lightgrey")
lines(c(20.5,24.5), c(0.05,0.05), col="lightgrey")
boxplot(cbind(fdr_muscle_2, fdr_muscle_5, fdr_muscle_10, fdr_muscle_20, fdr_muscle_50), 
        col=col_vector[1:4], pch=20, cex.axis=2, xaxt='n', ylim=c(0,0.8), 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(2.5,7.5,12.5,17.5,22.5), labels=c(2,5,10,20,50), tick=F, cex.axis=2)
legend("topright", fill=col_vector[1:4], bty='n', cex=2, ncol=1, 
       legend=c("diffVar", "Hybrid", "HMM, posterior threshold", "HMM, BFDR"))

rbind(colMeans(fdr_blood_2, na.rm=T), 
      colMeans(fdr_blood_5, na.rm=T), 
      colMeans(fdr_blood_10, na.rm=T), 
      colMeans(fdr_blood_20, na.rm=T), 
      colMeans(fdr_blood_50, na.rm=T))
rbind(colMeans(fdr_muscle_2, na.rm=T), 
      colMeans(fdr_muscle_5, na.rm=T), 
      colMeans(fdr_muscle_10, na.rm=T), 
      colMeans(fdr_muscle_20, na.rm=T), 
      colMeans(fdr_muscle_50, na.rm=T))


### TPR ####
positions = rep(1:20)
offsets <- c(rep(0, 4), rep(1, 4), rep(2, 4), rep(3, 4), rep(4, 4))

par(mfrow=c(2,1), mar=c(2,2.5,2,0.5), mgp=c(3,0.7,0))

plot(NA, xlim=c(1,24), ylim=c(0,1), xaxt="n", yaxt="n", 
     main="TPR - blood data", cex.main=2)
boxplot(cbind(tpr_blood_2, tpr_blood_5, tpr_blood_10, tpr_blood_20, tpr_blood_50), 
        col=col_vector[1:4], pch=20, cex.axis=2, xaxt='n', ylim=c(0,0.8), 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(2.5,7.5,12.5,17.5,22.5), labels=c(2,5,10,20,50), tick=F, cex.axis=2)
legend("topleft", fill=col_vector[1:4], bty='n', cex=2, ncol=1, 
       legend=c("diffVar", "Hybrid", "HMM, posterior threshold", "HMM, BFDR"))

plot(NA, xlim=c(1,24), ylim=c(0,1), xaxt="n", yaxt="n", 
     main="TPR - muscle data", cex.main=2)
boxplot(cbind(tpr_muscle_2, tpr_muscle_5, tpr_muscle_10, tpr_muscle_20, tpr_muscle_50), 
        col=col_vector[1:4], pch=20, cex.axis=2, xaxt='n', ylim=c(0,0.8), 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(2.5,7.5,12.5,17.5,22.5), labels=c(2,5,10,20,50), tick=F, cex.axis=2)

rbind(colMeans(tpr_blood_2, na.rm=T), 
      colMeans(tpr_blood_5, na.rm=T), 
      colMeans(tpr_blood_10, na.rm=T), 
      colMeans(tpr_blood_20, na.rm=T), 
      colMeans(tpr_blood_50, na.rm=T))
rbind(colMeans(tpr_muscle_2, na.rm=T), 
      colMeans(tpr_muscle_5, na.rm=T), 
      colMeans(tpr_muscle_10, na.rm=T), 
      colMeans(tpr_muscle_20, na.rm=T), 
      colMeans(tpr_muscle_50, na.rm=T))


### FDR curves blood and muscle ####
par(mfrow=c(2,5), mar=c(1.5,3.5,2.5,1), mgp=c(3,1,0))

plot(mean.discoveries_blood_2$blood_dV, mean.fdr_blood_2$blood_dV, 
     type='l', ylim=c(0,0.95), xlim=c(0,2000), lwd=2, col=col_vector[1], cex.axis=3, 
     xlab=NA, ylab=NA, xaxt='n', main=paste0("2"), cex.main=3)
axis(side=1, at=c(0,500,1000,1500,2000), tick=T, cex.axis=3, mgp=c(3,2,0))
text("Blood", x=0, y=0.9, cex=3, adj=0)
lines(mean.discoveries_blood_2$blood_voom_HM, mean.fdr_blood_2$blood_voom_HM, 
      lwd=2, col=col_vector[2], cex.axis=3)
lines(mean.discoveries_blood_2$blood_lnHMM, mean.fdr_blood_2$blood_lnHMM, 
      lwd=2, col=col_vector[3], cex.axis=3)
abline(h=0.05, col='lightgrey')

plot(mean.discoveries_blood_5$blood_dV, mean.fdr_blood_5$blood_dV, 
     type='l', ylim=c(0,0.95), xlim=c(0,2000), lwd=2, col=col_vector[1], cex.axis=3, 
     xlab=NA, ylab=NA, xaxt='n', yaxt='n', main=paste0("5"), cex.main=3)
axis(side=1, at=c(0,500,1000,1500,2000), tick=T, cex.axis=3, mgp=c(3,2,0))
lines(mean.discoveries_blood_5$blood_voom_HM, mean.fdr_blood_5$blood_voom_HM, 
      lwd=2, col=col_vector[2], cex.axis=3)
lines(mean.discoveries_blood_5$blood_lnHMM, mean.fdr_blood_5$blood_lnHMM, 
      lwd=2, col=col_vector[3], cex.axis=3)
abline(h=0.05, col='lightgrey')

plot(mean.discoveries_blood_10$blood_dV, mean.fdr_blood_10$blood_dV, 
     type='l', ylim=c(0,0.95), xlim=c(0,2000), lwd=2, col=col_vector[1], cex.axis=3, 
     xlab=NA, ylab=NA, xaxt='n', yaxt='n', main=paste0("10"), cex.main=3)
axis(side=1, at=c(0,500,1000,1500,2000), tick=T, cex.axis=3, mgp=c(3,2,0))
lines(mean.discoveries_blood_10$blood_voom_HM, mean.fdr_blood_10$blood_voom_HM, 
      lwd=2, col=col_vector[2], cex.axis=3)
lines(mean.discoveries_blood_10$blood_lnHMM, mean.fdr_blood_10$blood_lnHMM, 
      lwd=2, col=col_vector[3], cex.axis=3)
abline(h=0.05, col='lightgrey')

plot(mean.discoveries_blood_20$blood_dV, mean.fdr_blood_20$blood_dV, 
     type='l', ylim=c(0,0.95), xlim=c(0,2000), lwd=2, col=col_vector[1], cex.axis=3, 
     xlab=NA, ylab=NA, xaxt='n', yaxt='n', main=paste0("20"), cex.main=3)
axis(side=1, at=c(0,500,1000,1500,2000), tick=T, cex.axis=3, mgp=c(3,2,0))
lines(mean.discoveries_blood_20$blood_voom_HM, mean.fdr_blood_20$blood_voom_HM, 
      lwd=2, col=col_vector[2], cex.axis=3)
lines(mean.discoveries_blood_20$blood_lnHMM, mean.fdr_blood_20$blood_lnHMM, 
      lwd=2, col=col_vector[3], cex.axis=3)
abline(h=0.05, col='lightgrey')

plot(mean.discoveries_blood_50$blood_dV, mean.fdr_blood_50$blood_dV, 
     type='l', ylim=c(0,0.95), xlim=c(0,2000), lwd=2, col=col_vector[1], cex.axis=3, 
     xlab=NA, ylab=NA, xaxt='n', yaxt='n', main=paste0("50"), cex.main=3)
axis(side=1, at=c(0,500,1000,1500,2000), tick=T, cex.axis=3, mgp=c(3,2,0))
lines(mean.discoveries_blood_50$blood_voom_HM, mean.fdr_blood_50$blood_voom_HM, 
      lwd=2, col=col_vector[2], cex.axis=3)
lines(mean.discoveries_blood_50$blood_lnHMM, mean.fdr_blood_50$blood_lnHMM, 
      lwd=2, col=col_vector[3], cex.axis=3)
abline(h=0.05, col='lightgrey')

plot(mean.discoveries_muscle_2$muscle_dV, mean.fdr_muscle_2$muscle_dV, 
     type='l', ylim=c(0,0.95), xlim=c(0,2000), lwd=2, col=col_vector[1], cex.axis=3, 
     xlab=NA, ylab=NA, xaxt='n')
text("Muscle", x=0, y=0.9, cex=3, adj=0)
lines(mean.discoveries_muscle_2$muscle_voom_HM.lc, mean.fdr_muscle_2$muscle_voom_HM.lc, 
      lwd=2, col=col_vector[2], cex.axis=3)
lines(mean.discoveries_muscle_2$muscle_lnHMM.lc, mean.fdr_muscle_2$muscle_lnHMM.lc, 
      lwd=2, col=col_vector[3], cex.axis=3)
abline(h=0.05, col='lightgrey')

plot(mean.discoveries_muscle_5$muscle_dV, mean.fdr_muscle_5$muscle_dV, 
     type='l', ylim=c(0,0.95), xlim=c(0,2000), lwd=2, col=col_vector[1], cex.axis=3, 
     xlab=NA, ylab=NA, xaxt='n', yaxt='n')
lines(mean.discoveries_muscle_5$muscle_voom_HM.lc, mean.fdr_muscle_5$muscle_voom_HM.lc, 
      lwd=2, col=col_vector[2], cex.axis=3)
lines(mean.discoveries_muscle_5$muscle_lnHMM.lc, mean.fdr_muscle_5$muscle_lnHMM.lc, 
      lwd=2, col=col_vector[3], cex.axis=3)
abline(h=0.05, col='lightgrey')

plot(mean.discoveries_muscle_10$muscle_dV, mean.fdr_muscle_10$muscle_dV, 
     type='l', ylim=c(0,0.95), xlim=c(0,2000), lwd=2, col=col_vector[1], cex.axis=3, 
     xlab=NA, ylab=NA, xaxt='n', yaxt='n')
lines(mean.discoveries_muscle_10$muscle_voom_HM.lc, mean.fdr_muscle_10$muscle_voom_HM.lc, 
      lwd=2, col=col_vector[2], cex.axis=3)
lines(mean.discoveries_muscle_10$muscle_lnHMM.lc, mean.fdr_muscle_10$muscle_lnHMM.lc, 
      lwd=2, col=col_vector[3], cex.axis=3)
abline(h=0.05, col='lightgrey')

plot(mean.discoveries_muscle_20$muscle_dV, mean.fdr_muscle_20$muscle_dV, 
     type='l', ylim=c(0,0.95), xlim=c(0,2000), lwd=2, col=col_vector[1], cex.axis=3, 
     xlab=NA, ylab=NA, xaxt='n', yaxt='n')
lines(mean.discoveries_muscle_20$muscle_voom_HM.lc, mean.fdr_muscle_20$muscle_voom_HM.lc, 
      lwd=2, col=col_vector[2], cex.axis=3)
lines(mean.discoveries_muscle_20$muscle_lnHMM.lc, mean.fdr_muscle_20$muscle_lnHMM.lc, 
      lwd=2, col=col_vector[3], cex.axis=3)
abline(h=0.05, col='lightgrey')

plot(mean.discoveries_muscle_50$muscle_dV, mean.fdr_muscle_50$muscle_dV, 
     type='l', ylim=c(0,0.95), xlim=c(0,2000), lwd=2, col=col_vector[1], cex.axis=3, 
     xlab=NA, ylab=NA, xaxt='n', yaxt='n')
lines(mean.discoveries_muscle_50$muscle_voom_HM.lc, mean.fdr_muscle_50$muscle_voom_HM.lc, 
      lwd=2, col=col_vector[2], cex.axis=3)
lines(mean.discoveries_muscle_50$muscle_lnHMM.lc, mean.fdr_muscle_50$muscle_lnHMM.lc, 
      lwd=2, col=col_vector[3], cex.axis=3)
abline(h=0.05, col='lightgrey')

legend("topleft", col=col_vector[1:6], bty='n', cex=2.5, ncol=1, lty=1, lwd=2, 
       legend=c("diffVar", "Hybrid", "HMM"))

