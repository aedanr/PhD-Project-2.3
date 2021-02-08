library(here)
library(RColorBrewer)

## Load data, create colour vector ####
folder <- "Results/compcodeR DE, DD, DEDD results Feb 2020"
for (i in c("DEDD2", "DEDD5", "DEDD10", "DEDD20", "DEDD50")) {
  assign(paste0("DD.results.", i), 
         readRDS(here(folder, paste0("DD.lfc1.results.", i, ".rds"))))
}
for (i in c("DD2", "DD5", "DD10", "DD20", "DD50")) {
  assign(paste0("DD.results.", i), 
         readRDS(here(folder, paste0("DD.lfc1.results.", i, ".rds"))))
}
rm(i,folder)

qual_col_pals = brewer.pal.info[brewer.pal.info$category == "qual" & 
                                  brewer.pal.info$colorblind == T,]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector <- col_vector[c(1:2,9,11,26)]
rm(qual_col_pals)

## Keep only TMM and one version of each method ####
# MDSeq with ZI, lnHM.
{
  DD.results.DD2$pauc <- DD.results.DD2$pauc[, c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  DD.results.DD2$auc <- DD.results.DD2$auc[, c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  DD.results.DD2$fpr <- DD.results.DD2$fpr[, c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  DD.results.DD2$fdr <- DD.results.DD2$fdr[, c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  DD.results.DD2$tpr <- DD.results.DD2$tpr[, c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  DD.results.DD2$disp.lfc1.fdr <- DD.results.DD2$disp.lfc1.fdr[
    , c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")
    ]
  DD.results.DD2$mean.discoveries <- DD.results.DD2$mean.discoveries[
    , c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")
    ]
  
  DD.results.DD5$pauc <- DD.results.DD5$pauc[, c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  DD.results.DD5$auc <- DD.results.DD5$auc[, c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  DD.results.DD5$fpr <- DD.results.DD5$fpr[, c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  DD.results.DD5$fdr <- DD.results.DD5$fdr[, c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  DD.results.DD5$tpr <- DD.results.DD5$tpr[, c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  DD.results.DD5$disp.lfc1.fdr <- DD.results.DD5$disp.lfc1.fdr[
    , c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  DD.results.DD5$mean.discoveries <- DD.results.DD5$mean.discoveries[
    , c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  
  DD.results.DD10$pauc <- DD.results.DD10$pauc[, c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  DD.results.DD10$auc <- DD.results.DD10$auc[, c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  DD.results.DD10$fpr <- DD.results.DD10$fpr[, c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  DD.results.DD10$fdr <- DD.results.DD10$fdr[, c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  DD.results.DD10$tpr <- DD.results.DD10$tpr[, c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  DD.results.DD10$disp.lfc1.fdr <- DD.results.DD10$disp.lfc1.fdr[
    , c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  DD.results.DD10$mean.discoveries <- DD.results.DD10$mean.discoveries[
    , c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  
  DD.results.DD20$pauc <- DD.results.DD20$pauc[, c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  DD.results.DD20$auc <- DD.results.DD20$auc[, c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  DD.results.DD20$fpr <- DD.results.DD20$fpr[, c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  DD.results.DD20$fdr <- DD.results.DD20$fdr[, c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  DD.results.DD20$tpr <- DD.results.DD20$tpr[, c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  DD.results.DD20$disp.lfc1.fdr <- DD.results.DD20$disp.lfc1.fdr[
    , c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  DD.results.DD20$mean.discoveries <- DD.results.DD20$mean.discoveries[
    , c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  
  DD.results.DD50$pauc <- DD.results.DD50$pauc[, c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  DD.results.DD50$auc <- DD.results.DD50$auc[, c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  DD.results.DD50$fpr <- DD.results.DD50$fpr[, c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  DD.results.DD50$fdr <- DD.results.DD50$fdr[, c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  DD.results.DD50$tpr <- DD.results.DD50$tpr[, c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  DD.results.DD50$disp.lfc1.fdr <- DD.results.DD50$disp.lfc1.fdr[
    , c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  DD.results.DD50$mean.discoveries <- DD.results.DD50$mean.discoveries[
    , c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  
  DD.results.DEDD2$pauc <- DD.results.DEDD2$pauc[
    , c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  DD.results.DEDD2$auc <- DD.results.DEDD2$auc[
    , c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  DD.results.DEDD2$fpr <- DD.results.DEDD2$fpr[
    , c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  DD.results.DEDD2$fdr <- DD.results.DEDD2$fdr[
    , c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  DD.results.DEDD2$tpr <- DD.results.DEDD2$tpr[
    , c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  DD.results.DEDD2$disp.lfc1.fdr <- DD.results.DEDD2$disp.lfc1.fdr[
    , c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  DD.results.DEDD2$mean.discoveries <- DD.results.DEDD2$mean.discoveries[
    , c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  
  DD.results.DEDD5$pauc <- DD.results.DEDD5$pauc[
    , c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  DD.results.DEDD5$auc <- DD.results.DEDD5$auc[
    , c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  DD.results.DEDD5$fpr <- DD.results.DEDD5$fpr[
    , c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  DD.results.DEDD5$fdr <- DD.results.DEDD5$fdr[
    , c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  DD.results.DEDD5$tpr <- DD.results.DEDD5$tpr[
    , c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  DD.results.DEDD5$disp.lfc1.fdr <- DD.results.DEDD5$disp.lfc1.fdr[
    , c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  DD.results.DEDD5$mean.discoveries <- DD.results.DEDD5$mean.discoveries[
    , c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  
  DD.results.DEDD10$pauc <- DD.results.DEDD10$pauc[
    , c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  DD.results.DEDD10$auc <- DD.results.DEDD10$auc[
    , c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  DD.results.DEDD10$fpr <- DD.results.DEDD10$fpr[
    , c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  DD.results.DEDD10$fdr <- DD.results.DEDD10$fdr[
    , c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  DD.results.DEDD10$tpr <- DD.results.DEDD10$tpr[
    , c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  DD.results.DEDD10$disp.lfc1.fdr <- DD.results.DEDD10$disp.lfc1.fdr[
    , c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  DD.results.DEDD10$mean.discoveries <- DD.results.DEDD10$mean.discoveries[
    , c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  
  DD.results.DEDD20$pauc <- DD.results.DEDD20$pauc[
    , c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  DD.results.DEDD20$auc <- DD.results.DEDD20$auc[
    , c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  DD.results.DEDD20$fpr <- DD.results.DEDD20$fpr[
    , c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  DD.results.DEDD20$fdr <- DD.results.DEDD20$fdr[
    , c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  DD.results.DEDD20$tpr <- DD.results.DEDD20$tpr[
    , c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  DD.results.DEDD20$disp.lfc1.fdr <- DD.results.DEDD20$disp.lfc1.fdr[
    , c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  DD.results.DEDD20$mean.discoveries <- DD.results.DEDD20$mean.discoveries[
    , c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  
  DD.results.DEDD50$pauc <- DD.results.DEDD50$pauc[
    , c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  DD.results.DEDD50$auc <- DD.results.DEDD50$auc[
    , c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  DD.results.DEDD50$fpr <- DD.results.DEDD50$fpr[
    , c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  DD.results.DEDD50$fdr <- DD.results.DEDD50$fdr[
    , c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  DD.results.DEDD50$tpr <- DD.results.DEDD50$tpr[
    , c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  DD.results.DEDD50$disp.lfc1.fdr <- DD.results.DEDD50$disp.lfc1.fdr[
    , c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  DD.results.DEDD50$mean.discoveries <- DD.results.DEDD50$mean.discoveries[
    , c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
}

positions = rep(1:10)
offsets <- c(rep(0, 2), rep(1, 2), rep(2, 2), rep(3, 2), rep(4, 2))

### AUC ####
par(mfrow=c(2,1), mar=c(2.5,2.5,2.5,0.5), mgp=c(3,0.7,0))
plot(NA, xlim=c(0.75,14.25), ylim=c(0.5,0.86), xaxt="n", yaxt="n", 
     main="AUC - differences in dispersion only", cex.main=2)
boxplot(cbind(DD.results.DD2$auc, 
              DD.results.DD5$auc, 
              DD.results.DD10$auc, 
              DD.results.DD20$auc, 
              DD.results.DD50$auc), 
        col=col_vector[1:2], pch=20, cex.axis=2, xaxt='n', ylim=c(0,0.8), 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(1.5,4.5,7.5,10.5,13.5), labels=c(2,5,10,20,50), tick=F, cex.axis=2)
legend("topleft", fill=col_vector[1:2], bty='n', cex=2, ncol=1, 
       legend=c("MDSeq", "HM"))
plot(NA, xlim=c(0.75,14.25), ylim=c(0.5,0.86), xaxt="n", yaxt="n", 
     main="AUC - differences in mean and dispersion", cex.main=2)
boxplot(cbind(DD.results.DEDD2$auc, 
              DD.results.DEDD5$auc, 
              DD.results.DEDD10$auc, 
              DD.results.DEDD20$auc, 
              DD.results.DEDD50$auc), 
        col=col_vector[1:2], pch=20, cex.axis=2, xaxt='n', ylim=c(0,0.8), 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(1.5,4.5,7.5,10.5,13.5), labels=c(2,5,10,20,50), tick=F, cex.axis=2)

rbind(colMeans(DD.results.DD2$auc), 
      colMeans(DD.results.DD5$auc), 
      colMeans(DD.results.DD10$auc), 
      colMeans(DD.results.DD20$auc), 
      colMeans(DD.results.DD50$auc))
rbind(colMeans(DD.results.DEDD2$auc), 
      colMeans(DD.results.DEDD5$auc), 
      colMeans(DD.results.DEDD10$auc), 
      colMeans(DD.results.DEDD20$auc), 
      colMeans(DD.results.DEDD50$auc))

### pAUC ####
par(mfrow=c(2,1), mar=c(2.5,2.5,2.5,0.5), mgp=c(3,0.7,0))
plot(NA, xlim=c(0.75,14.25), ylim=c(0,0.034), xaxt="n", yaxt="n", 
     main="Partial AUC - differences in dispersion only", cex.main=2)
boxplot(cbind(DD.results.DD2$pauc, 
              DD.results.DD5$pauc, 
              DD.results.DD10$pauc, 
              DD.results.DD20$pauc, 
              DD.results.DD50$pauc), 
        col=col_vector[1:2], pch=20, cex.axis=2, xaxt='n', ylim=c(0,0.8), 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(1.5,4.5,7.5,10.5,13.5), labels=c(2,5,10,20,50), tick=F, cex.axis=2)
legend("topleft", fill=col_vector[1:2], bty='n', cex=2, ncol=1, 
       legend=c("MDSeq", "HM"))
plot(NA, xlim=c(0.75,14.25), ylim=c(0,0.034), xaxt="n", yaxt="n", 
     main="Partial AUC - differences in mean and dispersion", cex.main=2)
boxplot(cbind(DD.results.DEDD2$pauc, 
              DD.results.DEDD5$pauc, 
              DD.results.DEDD10$pauc, 
              DD.results.DEDD20$pauc, 
              DD.results.DEDD50$pauc), 
        col=col_vector[1:2], pch=20, cex.axis=2, xaxt='n', ylim=c(0,0.8), 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(1.5,4.5,7.5,10.5,13.5), labels=c(2,5,10,20,50), tick=F, cex.axis=2)

rbind(colMeans(DD.results.DD2$pauc), 
      colMeans(DD.results.DD5$pauc), 
      colMeans(DD.results.DD10$pauc), 
      colMeans(DD.results.DD20$pauc), 
      colMeans(DD.results.DD50$pauc))
rbind(colMeans(DD.results.DEDD2$pauc), 
      colMeans(DD.results.DEDD5$pauc), 
      colMeans(DD.results.DEDD10$pauc), 
      colMeans(DD.results.DEDD20$pauc), 
      colMeans(DD.results.DEDD50$pauc))

### FDR ####
par(mfrow=c(2,1), mar=c(2.5,2.5,2.5,0.5), mgp=c(3,0.7,0))
plot(NA, xlim=c(0.75,14.25), ylim=c(0,1), xaxt="n", yaxt="n", 
     main="FDR - differences in dispersion only", cex.main=2)
lines(c(0.5,2.5), c(0.05,0.05), col="lightgrey")
lines(c(3.5,5.5), c(0.05,0.05), col="lightgrey")
lines(c(6.5,8.5), c(0.05,0.05), col="lightgrey")
lines(c(9.5,11.5), c(0.05,0.05), col="lightgrey")
lines(c(12.5,14.5), c(0.05,0.05), col="lightgrey")
boxplot(cbind(DD.results.DD2$fdr, 
              DD.results.DD5$fdr, 
              DD.results.DD10$fdr, 
              DD.results.DD20$fdr, 
              DD.results.DD50$fdr), 
        col=col_vector[1:2], pch=20, cex.axis=2, xaxt='n', ylim=c(0,0.8), 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(1.5,4.5,7.5,10.5,13.5), labels=c(2,5,10,20,50), tick=F, cex.axis=2)
legend("topright", fill=col_vector[1:2], bty='n', cex=2, ncol=1, 
       legend=c("MDSeq", "HM"))
plot(NA, xlim=c(0.75,14.25), ylim=c(0,1), xaxt="n", yaxt="n", 
     main="FDR - differences in mean and dispersion", cex.main=2)
lines(c(0.5,2.5), c(0.05,0.05), col="lightgrey")
lines(c(3.5,5.5), c(0.05,0.05), col="lightgrey")
lines(c(6.5,8.5), c(0.05,0.05), col="lightgrey")
lines(c(9.5,11.5), c(0.05,0.05), col="lightgrey")
lines(c(12.5,14.5), c(0.05,0.05), col="lightgrey")
boxplot(cbind(DD.results.DEDD2$fdr, 
              DD.results.DEDD5$fdr, 
              DD.results.DEDD10$fdr, 
              DD.results.DEDD20$fdr, 
              DD.results.DEDD50$fdr), 
        col=col_vector[1:2], pch=20, cex.axis=2, xaxt='n', ylim=c(0,0.8), 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(1.5,4.5,7.5,10.5,13.5), labels=c(2,5,10,20,50), tick=F, cex.axis=2)

rbind(colMeans(DD.results.DD2$fdr, na.rm=T), 
      colMeans(DD.results.DD5$fdr, na.rm=T), 
      colMeans(DD.results.DD10$fdr, na.rm=T), 
      colMeans(DD.results.DD20$fdr, na.rm=T), 
      colMeans(DD.results.DD50$fdr, na.rm=T))
rbind(colMeans(DD.results.DEDD2$fdr, na.rm=T), 
      colMeans(DD.results.DEDD5$fdr, na.rm=T), 
      colMeans(DD.results.DEDD10$fdr, na.rm=T), 
      colMeans(DD.results.DEDD20$fdr, na.rm=T), 
      colMeans(DD.results.DEDD50$fdr, na.rm=T))

### TPR ####
par(mfrow=c(2,1), mar=c(2.5,2.5,2.5,0.5), mgp=c(3,0.7,0))
plot(NA, xlim=c(0.75,14.25), ylim=c(0,0.115), xaxt="n", yaxt="n", 
     main="TPR - differences in dispersion only", cex.main=2)
boxplot(cbind(DD.results.DD2$tpr, 
              DD.results.DD5$tpr, 
              DD.results.DD10$tpr, 
              DD.results.DD20$tpr, 
              DD.results.DD50$tpr), 
        col=col_vector[1:2], pch=20, cex.axis=2, xaxt='n', ylim=c(0,0.8), 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        main="Differences in dispersion only", cex.main=2, 
        add=T)
axis(side=1, at=c(1.5,4.5,7.5,10.5,13.5), labels=c(2,5,10,20,50), tick=F, cex.axis=2)
plot(NA, xlim=c(0.75,14.25), ylim=c(0,0.115), xaxt="n", yaxt="n", 
     main="TPR - differences in mean and dispersion", cex.main=2)
boxplot(cbind(DD.results.DEDD2$tpr, 
              DD.results.DEDD5$tpr, 
              DD.results.DEDD10$tpr, 
              DD.results.DEDD20$tpr, 
              DD.results.DEDD50$tpr), 
        col=col_vector[1:2], pch=20, cex.axis=2, xaxt='n', ylim=c(0,0.8), 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(1.5,4.5,7.5,10.5,13.5), labels=c(2,5,10,20,50), tick=F, cex.axis=2)
legend("topleft", fill=col_vector[1:2], bty='n', cex=2, ncol=1, 
       legend=c("MDSeq", "HM"))

rbind(colMeans(DD.results.DD2$tpr), 
      colMeans(DD.results.DD5$tpr), 
      colMeans(DD.results.DD10$tpr), 
      colMeans(DD.results.DD20$tpr), 
      colMeans(DD.results.DD50$tpr))
rbind(colMeans(DD.results.DEDD2$tpr), 
      colMeans(DD.results.DEDD5$tpr), 
      colMeans(DD.results.DEDD10$tpr), 
      colMeans(DD.results.DEDD20$tpr), 
      colMeans(DD.results.DEDD50$tpr))


### FDR curves ####
par(mfrow=c(2,5), mar=c(1.5,3.5,2.5,1), mgp=c(3,1,0))

plot(DD.results.DD2$mean.discoveries$disp.lfc1.MDSeq.zi.tmm, 
     DD.results.DD2$mean.fdr$disp.lfc1.MDSeq.zi.tmm, 
     type='l', ylim=c(0,0.95), xlim=c(0,2000), lwd=2, col=col_vector[1], cex.axis=3, 
     xlab=NA, ylab=NA, xaxt='n', main=paste0("2"), cex.main=3)
axis(side=1, at=c(0,500,1000,1500,2000), tick=T, cex.axis=3, mgp=c(3,2,0))
text("Differences in\ndispersion only", x=0, y=0.15, cex=2.5, adj=0)
lines(DD.results.DD2$mean.discoveries$disp.lfc1.lnHM.tmm, 
      DD.results.DD2$mean.fdr$disp.lfc1.lnHM.tmm, 
      lwd=2, col=col_vector[2], cex.axis=3)
abline(h=0.05, col='lightgrey')

plot(DD.results.DD5$mean.discoveries$disp.lfc1.MDSeq.zi.tmm, 
     DD.results.DD5$mean.fdr$disp.lfc1.MDSeq.zi.tmm, 
     type='l', ylim=c(0,0.95), xlim=c(0,2000), lwd=2, col=col_vector[1], cex.axis=3, 
     xlab=NA, ylab=NA, xaxt='n', yaxt='n', main=paste0("5"), cex.main=3)
axis(side=1, at=c(0,500,1000,1500,2000), tick=T, cex.axis=3, mgp=c(3,2,0))
lines(DD.results.DD5$mean.discoveries$disp.lfc1.lnHM.tmm, 
      DD.results.DD5$mean.fdr$disp.lfc1.lnHM.tmm, 
      lwd=2, col=col_vector[2], cex.axis=3)
abline(h=0.05, col='lightgrey')

plot(DD.results.DD10$mean.discoveries$disp.lfc1.MDSeq.zi.tmm, 
     DD.results.DD10$mean.fdr$disp.lfc1.MDSeq.zi.tmm, 
     type='l', ylim=c(0,0.95), xlim=c(0,2000), lwd=2, col=col_vector[1], cex.axis=3, 
     xlab=NA, ylab=NA, xaxt='n', yaxt='n', main=paste0("10"), cex.main=3)
axis(side=1, at=c(0,500,1000,1500,2000), tick=T, cex.axis=3, mgp=c(3,2,0))
lines(DD.results.DD10$mean.discoveries$disp.lfc1.lnHM.tmm, 
      DD.results.DD10$mean.fdr$disp.lfc1.lnHM.tmm, 
      lwd=2, col=col_vector[2], cex.axis=3)
abline(h=0.05, col='lightgrey')

plot(DD.results.DD20$mean.discoveries$disp.lfc1.MDSeq.zi.tmm, 
     DD.results.DD20$mean.fdr$disp.lfc1.MDSeq.zi.tmm, 
     type='l', ylim=c(0,0.95), xlim=c(0,2000), lwd=2, col=col_vector[1], cex.axis=3, 
     xlab=NA, ylab=NA, xaxt='n', yaxt='n', main=paste0("20"), cex.main=3)
axis(side=1, at=c(0,500,1000,1500,2000), tick=T, cex.axis=3, mgp=c(3,2,0))
lines(DD.results.DD20$mean.discoveries$disp.lfc1.lnHM.tmm, 
      DD.results.DD20$mean.fdr$disp.lfc1.lnHM.tmm, 
      lwd=2, col=col_vector[2], cex.axis=3)
abline(h=0.05, col='lightgrey')

plot(DD.results.DD50$mean.discoveries$disp.lfc1.MDSeq.zi.tmm, 
     DD.results.DD50$mean.fdr$disp.lfc1.MDSeq.zi.tmm, 
     type='l', ylim=c(0,0.95), xlim=c(0,2000), lwd=2, col=col_vector[1], cex.axis=3, 
     xlab=NA, ylab=NA, xaxt='n', yaxt='n', main=paste0("50"), cex.main=3)
axis(side=1, at=c(0,500,1000,1500,2000), tick=T, cex.axis=3, mgp=c(3,2,0))
lines(DD.results.DD50$mean.discoveries$disp.lfc1.lnHM.tmm, 
      DD.results.DD50$mean.fdr$disp.lfc1.lnHM.tmm, 
      lwd=2, col=col_vector[2], cex.axis=3)
abline(h=0.05, col='lightgrey')

plot(DD.results.DEDD2$mean.discoveries$disp.lfc1.MDSeq.zi.tmm, 
     DD.results.DEDD2$mean.fdr$disp.lfc1.MDSeq.zi.tmm, 
     type='l', ylim=c(0,0.95), xlim=c(0,2000), lwd=2, col=col_vector[1], cex.axis=3, 
     xlab=NA, ylab=NA, xaxt='n')
text("Differences in\nmean and dispersion", x=0, y=0.15, cex=2.5, adj=0)
lines(DD.results.DEDD2$mean.discoveries$disp.lfc1.lnHM.tmm, 
      DD.results.DEDD2$mean.fdr$disp.lfc1.lnHM.tmm, 
      lwd=2, col=col_vector[2], cex.axis=3)
abline(h=0.05, col='lightgrey')

plot(DD.results.DEDD5$mean.discoveries$disp.lfc1.MDSeq.zi.tmm, 
     DD.results.DEDD5$mean.fdr$disp.lfc1.MDSeq.zi.tmm, 
     type='l', ylim=c(0,0.95), xlim=c(0,2000), lwd=2, col=col_vector[1], cex.axis=3, 
     xlab=NA, ylab=NA, xaxt='n', yaxt='n')
lines(DD.results.DEDD5$mean.discoveries$disp.lfc1.lnHM.tmm, 
      DD.results.DEDD5$mean.fdr$disp.lfc1.lnHM.tmm, 
      lwd=2, col=col_vector[2], cex.axis=3)
abline(h=0.05, col='lightgrey')

plot(DD.results.DEDD10$mean.discoveries$disp.lfc1.MDSeq.zi.tmm, 
     DD.results.DEDD10$mean.fdr$disp.lfc1.MDSeq.zi.tmm, 
     type='l', ylim=c(0,0.95), xlim=c(0,2000), lwd=2, col=col_vector[1], cex.axis=3, 
     xlab=NA, ylab=NA, xaxt='n', yaxt='n')
lines(DD.results.DEDD10$mean.discoveries$disp.lfc1.lnHM.tmm, 
      DD.results.DEDD10$mean.fdr$disp.lfc1.lnHM.tmm, 
      lwd=2, col=col_vector[2], cex.axis=3)
abline(h=0.05, col='lightgrey')

plot(DD.results.DEDD20$mean.discoveries$disp.lfc1.MDSeq.zi.tmm, 
     DD.results.DEDD20$mean.fdr$disp.lfc1.MDSeq.zi.tmm, 
     type='l', ylim=c(0,0.95), xlim=c(0,2000), lwd=2, col=col_vector[1], cex.axis=3, 
     xlab=NA, ylab=NA, xaxt='n', yaxt='n')
lines(DD.results.DEDD20$mean.discoveries$disp.lfc1.lnHM.tmm, 
      DD.results.DEDD20$mean.fdr$disp.lfc1.lnHM.tmm, 
      lwd=2, col=col_vector[2], cex.axis=3)
abline(h=0.05, col='lightgrey')

plot(DD.results.DEDD50$mean.discoveries$disp.lfc1.MDSeq.zi.tmm, 
     DD.results.DEDD50$mean.fdr$disp.lfc1.MDSeq.zi.tmm, 
     type='l', ylim=c(0,0.95), xlim=c(0,2000), lwd=2, col=col_vector[1], cex.axis=3, 
     xlab=NA, ylab=NA, xaxt='n', yaxt='n')
lines(DD.results.DEDD50$mean.discoveries$disp.lfc1.lnHM.tmm, 
      DD.results.DEDD50$mean.fdr$disp.lfc1.lnHM.tmm, 
      lwd=2, col=col_vector[2], cex.axis=3)
abline(h=0.05, col='lightgrey')

legend("topleft", col=col_vector[1:2], bty='n', cex=2.5, ncol=1, lty=1, lwd=2, 
       legend=c("MDSeq", "HM"))
