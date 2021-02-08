library(here)
library(RColorBrewer)

## Load data, create colour vector ####
folder <- "Results/GTEx DD lfc1 results Jan 2021"
for (i in c("2", "5", "10", "20", "50")) {
  for (j in c("blood", "muscle")) {
    assign(paste0("DD.results.", j, ".DEDD", i), 
           readRDS(here(folder, paste0("DD.lfc1.results.", j, "_", i, "_DEDD.rds"))))
  }
}
rm(i,j,folder)

qual_col_pals = brewer.pal.info[brewer.pal.info$category == "qual" & 
                                  brewer.pal.info$colorblind == T,]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector <- col_vector[c(1:3,9,11,26)]
rm(qual_col_pals)


positions = rep(1:10)
offsets <- c(rep(0, 2), rep(1, 2), rep(2, 2), rep(3, 2), rep(4, 2))

### AUC ####
par(mfrow=c(2,1), mar=c(2,2.5,2,0.5), mgp=c(3,0.7,0))

plot(NA, xlim=c(0.75,14.25), ylim=c(0.5,0.86), xaxt="n", yaxt="n", 
     main="AUC - blood data", cex.main=2)
boxplot(cbind(DD.results.blood.DEDD2$auc, 
              DD.results.blood.DEDD5$auc, 
              DD.results.blood.DEDD10$auc, 
              DD.results.blood.DEDD20$auc, 
              DD.results.blood.DEDD50$auc), 
        col=col_vector[1:2], pch=20, cex.axis=2, xaxt='n', 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(1.5,4.5,7.5,10.5,13.5), labels=c(2,5,10,20,50), tick=F, cex.axis=2)
legend("bottomright", fill=col_vector[1:2], bty='n', cex=2, ncol=1, 
       legend=c("MDSeq", "HM"))

plot(NA, xlim=c(0.75,14.25), ylim=c(0.5,0.86), xaxt="n", yaxt="n", 
     main="AUC - muscle data", cex.main=2)
boxplot(cbind(DD.results.muscle.DEDD2$auc, 
              DD.results.muscle.DEDD5$auc, 
              DD.results.muscle.DEDD10$auc, 
              DD.results.muscle.DEDD20$auc, 
              DD.results.muscle.DEDD50$auc), 
        col=col_vector[1:2], pch=20, cex.axis=2, xaxt='n', 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(1.5,4.5,7.5,10.5,13.5), labels=c(2,5,10,20,50), tick=F, cex.axis=2)

rbind(colMeans(DD.results.blood.DEDD2$auc), 
      colMeans(DD.results.blood.DEDD5$auc), 
      colMeans(DD.results.blood.DEDD10$auc), 
      colMeans(DD.results.blood.DEDD20$auc), 
      colMeans(DD.results.blood.DEDD50$auc))
rbind(colMeans(DD.results.muscle.DEDD2$auc), 
      colMeans(DD.results.muscle.DEDD5$auc), 
      colMeans(DD.results.muscle.DEDD10$auc), 
      colMeans(DD.results.muscle.DEDD20$auc), 
      colMeans(DD.results.muscle.DEDD50$auc))


### pAUC ####
par(mfrow=c(2,1), mar=c(2,2.5,2,0.5), mgp=c(3,0.7,0))

plot(NA, xlim=c(0.75,14.25), ylim=c(0,0.032), xaxt="n", yaxt="n", 
     main="Partial AUC - blood data", cex.main=2)
boxplot(cbind(DD.results.blood.DEDD2$pauc, 
              DD.results.blood.DEDD5$pauc, 
              DD.results.blood.DEDD10$pauc, 
              DD.results.blood.DEDD20$pauc, 
              DD.results.blood.DEDD50$pauc), 
        col=col_vector[1:2], pch=20, cex.axis=2, xaxt='n', 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(1.5,4.5,7.5,10.5,13.5), labels=c(2,5,10,20,50), tick=F, cex.axis=2)

plot(NA, xlim=c(0.75,14.25), ylim=c(0,0.032), xaxt="n", yaxt="n", 
     main="Partial AUC - muscle data", cex.main=2)
boxplot(cbind(DD.results.muscle.DEDD2$pauc, 
              DD.results.muscle.DEDD5$pauc, 
              DD.results.muscle.DEDD10$pauc, 
              DD.results.muscle.DEDD20$pauc, 
              DD.results.muscle.DEDD50$pauc), 
        col=col_vector[1:2], pch=20, cex.axis=2, xaxt='n', 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(1.5,4.5,7.5,10.5,13.5), labels=c(2,5,10,20,50), tick=F, cex.axis=2)
legend("topleft", fill=col_vector[1:2], bty='n', cex=2, ncol=1, 
       legend=c("MDSeq", "HM"))

rbind(colMeans(DD.results.blood.DEDD2$pauc), 
      colMeans(DD.results.blood.DEDD5$pauc), 
      colMeans(DD.results.blood.DEDD10$pauc), 
      colMeans(DD.results.blood.DEDD20$pauc), 
      colMeans(DD.results.blood.DEDD50$pauc))
rbind(colMeans(DD.results.muscle.DEDD2$pauc), 
      colMeans(DD.results.muscle.DEDD5$pauc), 
      colMeans(DD.results.muscle.DEDD10$pauc), 
      colMeans(DD.results.muscle.DEDD20$pauc), 
      colMeans(DD.results.muscle.DEDD50$pauc))


### FDR ####
par(mfrow=c(2,1), mar=c(2,2.5,2,0.5), mgp=c(3,0.7,0))

plot(NA, xlim=c(0.75,14.25), ylim=c(0,0.95), xaxt="n", yaxt="n", 
     main="FDR - blood data", cex.main=2)
lines(c(0.5,2.5), c(0.05,0.05), col="lightgrey")
lines(c(3.5,5.5), c(0.05,0.05), col="lightgrey")
lines(c(6.5,8.5), c(0.05,0.05), col="lightgrey")
lines(c(9.5,11.5), c(0.05,0.05), col="lightgrey")
lines(c(12.5,14.5), c(0.05,0.05), col="lightgrey")
boxplot(cbind(DD.results.blood.DEDD2$fdr, 
              DD.results.blood.DEDD5$fdr, 
              DD.results.blood.DEDD10$fdr, 
              DD.results.blood.DEDD20$fdr, 
              DD.results.blood.DEDD50$fdr), 
        col=col_vector[1:2], pch=20, cex.axis=2, xaxt='n', 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(1.5,4.5,7.5,10.5,13.5), labels=c(2,5,10,20,50), tick=F, cex.axis=2)

plot(NA, xlim=c(0.75,14.25), ylim=c(0,0.95), xaxt="n", yaxt="n", 
     main="FDR - muscle data", cex.main=2)
lines(c(0.5,2.5), c(0.05,0.05), col="lightgrey")
lines(c(3.5,5.5), c(0.05,0.05), col="lightgrey")
lines(c(6.5,8.5), c(0.05,0.05), col="lightgrey")
lines(c(9.5,11.5), c(0.05,0.05), col="lightgrey")
lines(c(12.5,14.5), c(0.05,0.05), col="lightgrey")
boxplot(cbind(DD.results.muscle.DEDD2$fdr, 
              DD.results.muscle.DEDD5$fdr, 
              DD.results.muscle.DEDD10$fdr, 
              DD.results.muscle.DEDD20$fdr, 
              DD.results.muscle.DEDD50$fdr), 
        col=col_vector[1:2], pch=20, cex.axis=2, xaxt='n', 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(1.5,4.5,7.5,10.5,13.5), labels=c(2,5,10,20,50), tick=F, cex.axis=2)
legend("topright", fill=col_vector[1:2], bty='n', cex=2, ncol=1, 
       legend=c("MDSeq", "HM"))

rbind(colMeans(DD.results.blood.DEDD2$fdr), 
      colMeans(DD.results.blood.DEDD5$fdr), 
      colMeans(DD.results.blood.DEDD10$fdr), 
      colMeans(DD.results.blood.DEDD20$fdr), 
      colMeans(DD.results.blood.DEDD50$fdr))
rbind(colMeans(DD.results.muscle.DEDD2$fdr), 
      colMeans(DD.results.muscle.DEDD5$fdr), 
      colMeans(DD.results.muscle.DEDD10$fdr), 
      colMeans(DD.results.muscle.DEDD20$fdr), 
      colMeans(DD.results.muscle.DEDD50$fdr))


### TPR ####
par(mfrow=c(2,1), mar=c(2,2.5,2,0.5), mgp=c(3,0.7,0))

plot(NA, xlim=c(0.75,14.25), ylim=c(0,0.41), xaxt="n", yaxt="n", 
     main="TPR - blood data", cex.main=2)
boxplot(cbind(DD.results.blood.DEDD2$tpr, 
              DD.results.blood.DEDD5$tpr, 
              DD.results.blood.DEDD10$tpr, 
              DD.results.blood.DEDD20$tpr, 
              DD.results.blood.DEDD50$tpr), 
        col=col_vector[1:2], pch=20, cex.axis=2, xaxt='n', 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(1.5,4.5,7.5,10.5,13.5), labels=c(2,5,10,20,50), tick=F, cex.axis=2)

plot(NA, xlim=c(0.75,14.25), ylim=c(0,0.41), xaxt="n", yaxt="n", 
     main="TPR - muscle data", cex.main=2)
boxplot(cbind(DD.results.muscle.DEDD2$tpr, 
              DD.results.muscle.DEDD5$tpr, 
              DD.results.muscle.DEDD10$tpr, 
              DD.results.muscle.DEDD20$tpr, 
              DD.results.muscle.DEDD50$tpr), 
        col=col_vector[1:2], pch=20, cex.axis=2, xaxt='n', 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(1.5,4.5,7.5,10.5,13.5), labels=c(2,5,10,20,50), tick=F, cex.axis=2)
legend("topleft", fill=col_vector[1:2], bty='n', cex=2, ncol=1, 
       legend=c("MDSeq", "HM"))

rbind(colMeans(DD.results.blood.DEDD2$tpr), 
      colMeans(DD.results.blood.DEDD5$tpr), 
      colMeans(DD.results.blood.DEDD10$tpr), 
      colMeans(DD.results.blood.DEDD20$tpr), 
      colMeans(DD.results.blood.DEDD50$tpr))
rbind(colMeans(DD.results.muscle.DEDD2$tpr), 
      colMeans(DD.results.muscle.DEDD5$tpr), 
      colMeans(DD.results.muscle.DEDD10$tpr), 
      colMeans(DD.results.muscle.DEDD20$tpr), 
      colMeans(DD.results.muscle.DEDD50$tpr))


### FDR curves blood and muscle ####
par(mfrow=c(2,5), mar=c(1.5,3.5,2.5,1), mgp=c(3,1,0))

plot(DD.results.blood.DEDD2$mean.discoveries$disp.lfc1.MDSeq.zi, 
     DD.results.blood.DEDD2$mean.fdr$disp.lfc1.MDSeq.zi, 
     type='l', ylim=c(0,0.95), xlim=c(0,2000), lwd=2, col=col_vector[1], cex.axis=3, 
     xlab=NA, ylab=NA, xaxt='n', main=paste0("2"), cex.main=3)
axis(side=1, at=c(0,500,1000,1500,2000), tick=T, cex.axis=3, mgp=c(3,2,0))
text("Blood", x=100, y=0.1, cex=3, adj=0)
lines(DD.results.blood.DEDD2$mean.discoveries$disp.lfc1.lnHM, 
      DD.results.blood.DEDD2$mean.fdr$disp.lfc1.lnHM, 
      lwd=2, col=col_vector[2], cex.axis=3)
abline(h=0.05, col='lightgrey')

plot(DD.results.blood.DEDD5$mean.discoveries$disp.lfc1.MDSeq.zi, 
     DD.results.blood.DEDD5$mean.fdr$disp.lfc1.MDSeq.zi, 
     type='l', ylim=c(0,0.95), xlim=c(0,2000), lwd=2, col=col_vector[1], cex.axis=3, 
     xlab=NA, ylab=NA, xaxt='n', yaxt='n', main=paste0("5"), cex.main=3)
axis(side=1, at=c(0,500,1000,1500,2000), tick=T, cex.axis=3, mgp=c(3,2,0))
lines(DD.results.blood.DEDD5$mean.discoveries$disp.lfc1.lnHM, 
      DD.results.blood.DEDD5$mean.fdr$disp.lfc1.lnHM, 
      lwd=2, col=col_vector[2], cex.axis=3)
abline(h=0.05, col='lightgrey')

plot(DD.results.blood.DEDD10$mean.discoveries$disp.lfc1.MDSeq.zi, 
     DD.results.blood.DEDD10$mean.fdr$disp.lfc1.MDSeq.zi, 
     type='l', ylim=c(0,0.95), xlim=c(0,2000), lwd=2, col=col_vector[1], cex.axis=3, 
     xlab=NA, ylab=NA, xaxt='n', yaxt='n', main=paste0("10"), cex.main=3)
axis(side=1, at=c(0,500,1000,1500,2000), tick=T, cex.axis=3, mgp=c(3,2,0))
lines(DD.results.blood.DEDD10$mean.discoveries$disp.lfc1.lnHM, 
      DD.results.blood.DEDD10$mean.fdr$disp.lfc1.lnHM, 
      lwd=2, col=col_vector[2], cex.axis=3)
abline(h=0.05, col='lightgrey')

plot(DD.results.blood.DEDD20$mean.discoveries$disp.lfc1.MDSeq.zi, 
     DD.results.blood.DEDD20$mean.fdr$disp.lfc1.MDSeq.zi, 
     type='l', ylim=c(0,0.95), xlim=c(0,2000), lwd=2, col=col_vector[1], cex.axis=3, 
     xlab=NA, ylab=NA, xaxt='n', yaxt='n', main=paste0("20"), cex.main=3)
axis(side=1, at=c(0,500,1000,1500,2000), tick=T, cex.axis=3, mgp=c(3,2,0))
lines(DD.results.blood.DEDD20$mean.discoveries$disp.lfc1.lnHM, 
      DD.results.blood.DEDD20$mean.fdr$disp.lfc1.lnHM, 
      lwd=2, col=col_vector[2], cex.axis=3)
abline(h=0.05, col='lightgrey')

plot(DD.results.blood.DEDD50$mean.discoveries$disp.lfc1.MDSeq.zi, 
     DD.results.blood.DEDD50$mean.fdr$disp.lfc1.MDSeq.zi, 
     type='l', ylim=c(0,0.95), xlim=c(0,2000), lwd=2, col=col_vector[1], cex.axis=3, 
     xlab=NA, ylab=NA, xaxt='n', yaxt='n', main=paste0("50"), cex.main=3)
axis(side=1, at=c(0,500,1000,1500,2000), tick=T, cex.axis=3, mgp=c(3,2,0))
lines(DD.results.blood.DEDD50$mean.discoveries$disp.lfc1.lnHM, 
      DD.results.blood.DEDD50$mean.fdr$disp.lfc1.lnHM, 
      lwd=2, col=col_vector[2], cex.axis=3)
abline(h=0.05, col='lightgrey')

plot(DD.results.muscle.DEDD2$mean.discoveries$disp.lfc1.MDSeq.zi, 
     DD.results.muscle.DEDD2$mean.fdr$disp.lfc1.MDSeq.zi, 
     type='l', ylim=c(0,0.95), xlim=c(0,2000), lwd=2, col=col_vector[1], cex.axis=3, 
     xlab=NA, ylab=NA, xaxt='n')
text("Muscle", x=100, y=0.1, cex=3, adj=0)
lines(DD.results.muscle.DEDD2$mean.discoveries$disp.lfc1.lnHM, 
      DD.results.muscle.DEDD2$mean.fdr$disp.lfc1.lnHM, 
      lwd=2, col=col_vector[2], cex.axis=3)
abline(h=0.05, col='lightgrey')

plot(DD.results.muscle.DEDD5$mean.discoveries$disp.lfc1.MDSeq.zi, 
     DD.results.muscle.DEDD5$mean.fdr$disp.lfc1.MDSeq.zi, 
     type='l', ylim=c(0,0.95), xlim=c(0,2000), lwd=2, col=col_vector[1], cex.axis=3, 
     xlab=NA, ylab=NA, xaxt='n', yaxt='n')
lines(DD.results.muscle.DEDD5$mean.discoveries$disp.lfc1.lnHM, 
      DD.results.muscle.DEDD5$mean.fdr$disp.lfc1.lnHM, 
      lwd=2, col=col_vector[2], cex.axis=3)
abline(h=0.05, col='lightgrey')

plot(DD.results.muscle.DEDD10$mean.discoveries$disp.lfc1.MDSeq.zi, 
     DD.results.muscle.DEDD10$mean.fdr$disp.lfc1.MDSeq.zi, 
     type='l', ylim=c(0,0.95), xlim=c(0,2000), lwd=2, col=col_vector[1], cex.axis=3, 
     xlab=NA, ylab=NA, xaxt='n', yaxt='n')
lines(DD.results.muscle.DEDD10$mean.discoveries$disp.lfc1.lnHM, 
      DD.results.muscle.DEDD10$mean.fdr$disp.lfc1.lnHM, 
      lwd=2, col=col_vector[2], cex.axis=3)
abline(h=0.05, col='lightgrey')

plot(DD.results.muscle.DEDD20$mean.discoveries$disp.lfc1.MDSeq.zi, 
     DD.results.muscle.DEDD20$mean.fdr$disp.lfc1.MDSeq.zi, 
     type='l', ylim=c(0,0.95), xlim=c(0,2000), lwd=2, col=col_vector[1], cex.axis=3, 
     xlab=NA, ylab=NA, xaxt='n', yaxt='n')
lines(DD.results.muscle.DEDD20$mean.discoveries$disp.lfc1.lnHM, 
      DD.results.muscle.DEDD20$mean.fdr$disp.lfc1.lnHM, 
      lwd=2, col=col_vector[2], cex.axis=3)
abline(h=0.05, col='lightgrey')

plot(DD.results.muscle.DEDD50$mean.discoveries$disp.lfc1.MDSeq.zi, 
     DD.results.muscle.DEDD50$mean.fdr$disp.lfc1.MDSeq.zi, 
     type='l', ylim=c(0,0.95), xlim=c(0,2000), lwd=2, col=col_vector[1], cex.axis=3, 
     xlab=NA, ylab=NA, xaxt='n', yaxt='n')
lines(DD.results.muscle.DEDD50$mean.discoveries$disp.lfc1.lnHM, 
      DD.results.muscle.DEDD50$mean.fdr$disp.lfc1.lnHM, 
      lwd=2, col=col_vector[2], cex.axis=3)
abline(h=0.05, col='lightgrey')

legend("topleft", col=col_vector[1:2], bty='n', cex=2.5, ncol=1, lty=1, lwd=2, 
       legend=c("MDSeq", "HM"))

