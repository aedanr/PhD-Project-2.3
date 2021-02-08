library(here)
library(RColorBrewer)

## Load data, create colour vector ####
folder <- "Results/compcodeR DE, DD, DEDD results Feb 2020"
for (i in c("DE2", "DE5", "DE10", "DE20", "DE50")) {
  assign(paste0("DEDD.results.", i), 
         readRDS(here(folder, paste0("DEDD.results.", i, ".rds"))))
}
for (i in c("DD2", "DD5", "DD10", "DD20", "DD50")) {
  assign(paste0("DEDD.results.", i), 
         readRDS(here(folder, paste0("DEDD.results.", i, ".rds"))))
}
folder <- "Results/compcodeR combined results Dec 2020 including diffVar"
for (i in c("DEDD2", "DEDD5", "DEDD10", "DEDD20", "DEDD50")) {
  assign(paste0("DEDD.results.", i), 
         readRDS(here(folder, paste0("DEDD.results.", i, ".rds"))))
}
rm(i,folder)

qual_col_pals = brewer.pal.info[brewer.pal.info$category == "qual" & 
                                  brewer.pal.info$colorblind == T,]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector <- col_vector[c(1:3,9,11,26)]
rm(qual_col_pals)


## Keep only TMM, for diffVar, edgeR/lnHM.log, lnHMM (posterior threshold and BFDR) ####
# edgeR/HM hybrid since edgeR was best DE for compcodeR data; switch to voom for GTEx 
# since voom was best there.
{
  DEDD.results.DD2$pauc <- DEDD.results.DD2$pauc$lnHMM.tmm
  DEDD.results.DD2$auc <- DEDD.results.DD2$auc$lnHMM.tmm
  DEDD.results.DD2$fpr <- DEDD.results.DD2$fpr[, c("thr.lnHMM.tmm", "bfdr.lnHMM.tmm")]
  DEDD.results.DD2$fdr <- DEDD.results.DD2$fdr[, c("thr.lnHMM.tmm", "bfdr.lnHMM.tmm")]
  DEDD.results.DD2$tpr <- DEDD.results.DD2$tpr[, c("thr.lnHMM.tmm", "bfdr.lnHMM.tmm")]
  DEDD.results.DD2$mean.fdr <- DEDD.results.DD2$mean.fdr$lnHMM.tmm
  DEDD.results.DD2$mean.discoveries <- DEDD.results.DD2$mean.discoveries$lnHMM.tmm
  DEDD.results.DD5$pauc <- DEDD.results.DD5$pauc$lnHMM.tmm
  DEDD.results.DD5$auc <- DEDD.results.DD5$auc$lnHMM.tmm
  DEDD.results.DD5$fpr <- DEDD.results.DD5$fpr[, c("thr.lnHMM.tmm", "bfdr.lnHMM.tmm")]
  DEDD.results.DD5$fdr <- DEDD.results.DD5$fdr[, c("thr.lnHMM.tmm", "bfdr.lnHMM.tmm")]
  DEDD.results.DD5$tpr <- DEDD.results.DD5$tpr[, c("thr.lnHMM.tmm", "bfdr.lnHMM.tmm")]
  DEDD.results.DD5$mean.fdr <- DEDD.results.DD5$mean.fdr$lnHMM.tmm
  DEDD.results.DD5$mean.discoveries <- DEDD.results.DD5$mean.discoveries$lnHMM.tmm
  DEDD.results.DD10$pauc <- DEDD.results.DD10$pauc$lnHMM.tmm
  DEDD.results.DD10$auc <- DEDD.results.DD10$auc$lnHMM.tmm
  DEDD.results.DD10$fpr <- DEDD.results.DD10$fpr[, c("thr.lnHMM.tmm", "bfdr.lnHMM.tmm")]
  DEDD.results.DD10$fdr <- DEDD.results.DD10$fdr[, c("thr.lnHMM.tmm", "bfdr.lnHMM.tmm")]
  DEDD.results.DD10$tpr <- DEDD.results.DD10$tpr[, c("thr.lnHMM.tmm", "bfdr.lnHMM.tmm")]
  DEDD.results.DD10$mean.fdr <- DEDD.results.DD10$mean.fdr$lnHMM.tmm
  DEDD.results.DD10$mean.discoveries <- DEDD.results.DD10$mean.discoveries$lnHMM.tmm
  DEDD.results.DD20$pauc <- DEDD.results.DD20$pauc$lnHMM.tmm
  DEDD.results.DD20$auc <- DEDD.results.DD20$auc$lnHMM.tmm
  DEDD.results.DD20$fpr <- DEDD.results.DD20$fpr[, c("thr.lnHMM.tmm", "bfdr.lnHMM.tmm")]
  DEDD.results.DD20$fdr <- DEDD.results.DD20$fdr[, c("thr.lnHMM.tmm", "bfdr.lnHMM.tmm")]
  DEDD.results.DD20$tpr <- DEDD.results.DD20$tpr[, c("thr.lnHMM.tmm", "bfdr.lnHMM.tmm")]
  DEDD.results.DD20$mean.fdr <- DEDD.results.DD20$mean.fdr$lnHMM.tmm
  DEDD.results.DD20$mean.discoveries <- DEDD.results.DD20$mean.discoveries$lnHMM.tmm
  DEDD.results.DD50$pauc <- DEDD.results.DD50$pauc$lnHMM.tmm
  DEDD.results.DD50$auc <- DEDD.results.DD50$auc$lnHMM.tmm
  DEDD.results.DD50$fpr <- DEDD.results.DD50$fpr[, c("thr.lnHMM.tmm", "bfdr.lnHMM.tmm")]
  DEDD.results.DD50$fdr <- DEDD.results.DD50$fdr[, c("thr.lnHMM.tmm", "bfdr.lnHMM.tmm")]
  DEDD.results.DD50$tpr <- DEDD.results.DD50$tpr[, c("thr.lnHMM.tmm", "bfdr.lnHMM.tmm")]
  DEDD.results.DD50$mean.fdr <- DEDD.results.DD50$mean.fdr$lnHMM.tmm
  DEDD.results.DD50$mean.discoveries <- DEDD.results.DD50$mean.discoveries$lnHMM.tmm
  
  DEDD.results.DE2$pauc <- DEDD.results.DE2$pauc$lnHMM.tmm
  DEDD.results.DE2$auc <- DEDD.results.DE2$auc$lnHMM.tmm
  DEDD.results.DE2$fpr <- DEDD.results.DE2$fpr[, c("thr.lnHMM.tmm", "bfdr.lnHMM.tmm")]
  DEDD.results.DE2$fdr <- DEDD.results.DE2$fdr[, c("thr.lnHMM.tmm", "bfdr.lnHMM.tmm")]
  DEDD.results.DE2$tpr <- DEDD.results.DE2$tpr[, c("thr.lnHMM.tmm", "bfdr.lnHMM.tmm")]
  DEDD.results.DE2$mean.fdr <- DEDD.results.DE2$mean.fdr$lnHMM.tmm
  DEDD.results.DE2$mean.discoveries <- DEDD.results.DE2$mean.discoveries$lnHMM.tmm
  DEDD.results.DE5$pauc <- DEDD.results.DE5$pauc$lnHMM.tmm
  DEDD.results.DE5$auc <- DEDD.results.DE5$auc$lnHMM.tmm
  DEDD.results.DE5$fpr <- DEDD.results.DE5$fpr[, c("thr.lnHMM.tmm", "bfdr.lnHMM.tmm")]
  DEDD.results.DE5$fdr <- DEDD.results.DE5$fdr[, c("thr.lnHMM.tmm", "bfdr.lnHMM.tmm")]
  DEDD.results.DE5$tpr <- DEDD.results.DE5$tpr[, c("thr.lnHMM.tmm", "bfdr.lnHMM.tmm")]
  DEDD.results.DE5$mean.fdr <- DEDD.results.DE5$mean.fdr$lnHMM.tmm
  DEDD.results.DE5$mean.discoveries <- DEDD.results.DE5$mean.discoveries$lnHMM.tmm
  DEDD.results.DE10$pauc <- DEDD.results.DE10$pauc$lnHMM.tmm
  DEDD.results.DE10$auc <- DEDD.results.DE10$auc$lnHMM.tmm
  DEDD.results.DE10$fpr <- DEDD.results.DE10$fpr[, c("thr.lnHMM.tmm", "bfdr.lnHMM.tmm")]
  DEDD.results.DE10$fdr <- DEDD.results.DE10$fdr[, c("thr.lnHMM.tmm", "bfdr.lnHMM.tmm")]
  DEDD.results.DE10$tpr <- DEDD.results.DE10$tpr[, c("thr.lnHMM.tmm", "bfdr.lnHMM.tmm")]
  DEDD.results.DE10$mean.fdr <- DEDD.results.DE10$mean.fdr$lnHMM.tmm
  DEDD.results.DE10$mean.discoveries <- DEDD.results.DE10$mean.discoveries$lnHMM.tmm
  DEDD.results.DE20$pauc <- DEDD.results.DE20$pauc$lnHMM.tmm
  DEDD.results.DE20$auc <- DEDD.results.DE20$auc$lnHMM.tmm
  DEDD.results.DE20$fpr <- DEDD.results.DE20$fpr[, c("thr.lnHMM.tmm", "bfdr.lnHMM.tmm")]
  DEDD.results.DE20$fdr <- DEDD.results.DE20$fdr[, c("thr.lnHMM.tmm", "bfdr.lnHMM.tmm")]
  DEDD.results.DE20$tpr <- DEDD.results.DE20$tpr[, c("thr.lnHMM.tmm", "bfdr.lnHMM.tmm")]
  DEDD.results.DE20$mean.fdr <- DEDD.results.DE20$mean.fdr$lnHMM.tmm
  DEDD.results.DE20$mean.discoveries <- DEDD.results.DE20$mean.discoveries$lnHMM.tmm
  DEDD.results.DE50$pauc <- DEDD.results.DE50$pauc$lnHMM.tmm
  DEDD.results.DE50$auc <- DEDD.results.DE50$auc$lnHMM.tmm
  DEDD.results.DE50$fpr <- DEDD.results.DE50$fpr[, c("thr.lnHMM.tmm", "bfdr.lnHMM.tmm")]
  DEDD.results.DE50$fdr <- DEDD.results.DE50$fdr[, c("thr.lnHMM.tmm", "bfdr.lnHMM.tmm")]
  DEDD.results.DE50$tpr <- DEDD.results.DE50$tpr[, c("thr.lnHMM.tmm", "bfdr.lnHMM.tmm")]
  DEDD.results.DE50$mean.fdr <- DEDD.results.DE50$mean.fdr$lnHMM.tmm
  DEDD.results.DE50$mean.discoveries <- DEDD.results.DE50$mean.discoveries$lnHMM.tmm
  
  DEDD.results.DEDD2$pauc <- DEDD.results.DEDD2$pauc[, c("dV.tmm", "edgeR_HM.tmm", "lnHMM.tmm")]
  DEDD.results.DEDD2$auc <- DEDD.results.DEDD2$auc[, c("dV.tmm", "edgeR_HM.tmm", "lnHMM.tmm")]
  DEDD.results.DEDD2$fpr <- DEDD.results.DEDD2$fpr[
    , c("dV.tmm", "edgeR_HM.tmm", "thr.lnHMM.tmm", "bfdr.lnHMM.tmm")]
  DEDD.results.DEDD2$fdr <- DEDD.results.DEDD2$fdr[
    , c("dV.tmm", "edgeR_HM.tmm", "thr.lnHMM.tmm", "bfdr.lnHMM.tmm")]
  DEDD.results.DEDD2$tpr <- DEDD.results.DEDD2$tpr[
    , c("dV.tmm", "edgeR_HM.tmm", "thr.lnHMM.tmm", "bfdr.lnHMM.tmm")]
  DEDD.results.DEDD2$mean.fdr <- DEDD.results.DEDD2$mean.fdr[
    , c("dV.tmm", "edgeR_HM.tmm", "lnHMM.tmm")]
  DEDD.results.DEDD2$mean.discoveries <- DEDD.results.DEDD2$mean.discoveries[
    , c("dV.tmm", "edgeR_HM.tmm", "lnHMM.tmm")]
  
  DEDD.results.DEDD5$pauc <- DEDD.results.DEDD5$pauc[, c("dV.tmm", "edgeR_HM.tmm", "lnHMM.tmm")]
  DEDD.results.DEDD5$auc <- DEDD.results.DEDD5$auc[, c("dV.tmm", "edgeR_HM.tmm", "lnHMM.tmm")]
  DEDD.results.DEDD5$fpr <- DEDD.results.DEDD5$fpr[
    , c("dV.tmm", "edgeR_HM.tmm", "thr.lnHMM.tmm", "bfdr.lnHMM.tmm")]
  DEDD.results.DEDD5$fdr <- DEDD.results.DEDD5$fdr[
    , c("dV.tmm", "edgeR_HM.tmm", "thr.lnHMM.tmm", "bfdr.lnHMM.tmm")]
  DEDD.results.DEDD5$tpr <- DEDD.results.DEDD5$tpr[
    , c("dV.tmm", "edgeR_HM.tmm", "thr.lnHMM.tmm", "bfdr.lnHMM.tmm")]
  DEDD.results.DEDD5$mean.fdr <- DEDD.results.DEDD5$mean.fdr[
    , c("dV.tmm", "edgeR_HM.tmm", "lnHMM.tmm")]
  DEDD.results.DEDD5$mean.discoveries <- DEDD.results.DEDD5$mean.discoveries[
    , c("dV.tmm", "edgeR_HM.tmm", "lnHMM.tmm")]
  
  DEDD.results.DEDD10$pauc <- DEDD.results.DEDD10$pauc[, c("dV.tmm", "edgeR_HM.tmm", "lnHMM.tmm")]
  DEDD.results.DEDD10$auc <- DEDD.results.DEDD10$auc[, c("dV.tmm", "edgeR_HM.tmm", "lnHMM.tmm")]
  DEDD.results.DEDD10$fpr <- DEDD.results.DEDD10$fpr[
    , c("dV.tmm", "edgeR_HM.tmm", "thr.lnHMM.tmm", "bfdr.lnHMM.tmm")]
  DEDD.results.DEDD10$fdr <- DEDD.results.DEDD10$fdr[
    , c("dV.tmm", "edgeR_HM.tmm", "thr.lnHMM.tmm", "bfdr.lnHMM.tmm")]
  DEDD.results.DEDD10$tpr <- DEDD.results.DEDD10$tpr[
    , c("dV.tmm", "edgeR_HM.tmm", "thr.lnHMM.tmm", "bfdr.lnHMM.tmm")]
  DEDD.results.DEDD10$mean.fdr <- DEDD.results.DEDD10$mean.fdr[
    , c("dV.tmm", "edgeR_HM.tmm", "lnHMM.tmm")]
  DEDD.results.DEDD10$mean.discoveries <- DEDD.results.DEDD10$mean.discoveries[
    , c("dV.tmm", "edgeR_HM.tmm", "lnHMM.tmm")]
  
  DEDD.results.DEDD20$pauc <- DEDD.results.DEDD20$pauc[, c("dV.tmm", "edgeR_HM.tmm", "lnHMM.tmm")]
  DEDD.results.DEDD20$auc <- DEDD.results.DEDD20$auc[, c("dV.tmm", "edgeR_HM.tmm", "lnHMM.tmm")]
  DEDD.results.DEDD20$fpr <- DEDD.results.DEDD20$fpr[
    , c("dV.tmm", "edgeR_HM.tmm", "thr.lnHMM.tmm", "bfdr.lnHMM.tmm")]
  DEDD.results.DEDD20$fdr <- DEDD.results.DEDD20$fdr[
    , c("dV.tmm", "edgeR_HM.tmm", "thr.lnHMM.tmm", "bfdr.lnHMM.tmm")]
  DEDD.results.DEDD20$tpr <- DEDD.results.DEDD20$tpr[
    , c("dV.tmm", "edgeR_HM.tmm", "thr.lnHMM.tmm", "bfdr.lnHMM.tmm")]
  DEDD.results.DEDD20$mean.fdr <- DEDD.results.DEDD20$mean.fdr[
    , c("dV.tmm", "edgeR_HM.tmm", "lnHMM.tmm")]
  DEDD.results.DEDD20$mean.discoveries <- DEDD.results.DEDD20$mean.discoveries[
    , c("dV.tmm", "edgeR_HM.tmm", "lnHMM.tmm")]
  
  DEDD.results.DEDD50$pauc <- DEDD.results.DEDD50$pauc[, c("dV.tmm", "edgeR_HM.tmm", "lnHMM.tmm")]
  DEDD.results.DEDD50$auc <- DEDD.results.DEDD50$auc[, c("dV.tmm", "edgeR_HM.tmm", "lnHMM.tmm")]
  DEDD.results.DEDD50$fpr <- DEDD.results.DEDD50$fpr[
    , c("dV.tmm", "edgeR_HM.tmm", "thr.lnHMM.tmm", "bfdr.lnHMM.tmm")]
  DEDD.results.DEDD50$fdr <- DEDD.results.DEDD50$fdr[
    , c("dV.tmm", "edgeR_HM.tmm", "thr.lnHMM.tmm", "bfdr.lnHMM.tmm")]
  DEDD.results.DEDD50$tpr <- DEDD.results.DEDD50$tpr[
    , c("dV.tmm", "edgeR_HM.tmm", "thr.lnHMM.tmm", "bfdr.lnHMM.tmm")]
  DEDD.results.DEDD50$mean.fdr <- DEDD.results.DEDD50$mean.fdr[
    , c("dV.tmm", "edgeR_HM.tmm", "lnHMM.tmm")]
  DEDD.results.DEDD50$mean.discoveries <- DEDD.results.DEDD50$mean.discoveries[
    , c("dV.tmm", "edgeR_HM.tmm", "lnHMM.tmm")]
}

### AUC ####
positions = rep(1:15)
offsets <- c(rep(0, 3), rep(1, 3), rep(2, 3), rep(3, 3), rep(4, 3))

par(mfrow=c(3,1), mar=c(3.5,3.5,3.5,0.5), mgp=c(3,1,0))

plot(NA, xlim=c(0.7,5.3), ylim=c(0.48,0.98), xaxt="n", yaxt="n", ylab=NA, xlab=NA, 
     main="AUC - differences in mean only", cex.main=3)
boxplot(cbind(DEDD.results.DE5$auc, 
              DEDD.results.DE5$auc, 
              DEDD.results.DE10$auc, 
              DEDD.results.DE20$auc, 
              DEDD.results.DE50$auc), 
        col=col_vector[3], border=col_vector[3], pch=20, cex.axis=3, xaxt='n', 
        lwd=0.2, cex=0.5, boxwex=0.6, 
        add=T)
axis(side=1, at=1:5, labels=c(2,5,10,20,50), tick=F, cex.axis=3, mgp=c(3,1.5,0))

plot(NA, xlim=c(0.7,5.3), ylim=c(0.48,0.98), xaxt="n", yaxt="n", ylab=NA, xlab=NA, 
     main="AUC - differences in dispersion only", cex.main=3)
boxplot(cbind(DEDD.results.DD2$auc, 
              DEDD.results.DD5$auc, 
              DEDD.results.DD10$auc, 
              DEDD.results.DD20$auc, 
              DEDD.results.DD50$auc), 
        col=col_vector[3], border=col_vector[3], pch=20, cex.axis=3, xaxt='n',
        lwd=0.2, cex=0.5, boxwex=0.6, 
        add=T)
axis(side=1, at=1:5, labels=c(2,5,10,20,50), tick=F, cex.axis=3, mgp=c(3,1.5,0))

plot(NA, xlim=c(1,19), ylim=c(0.48,0.98), xaxt="n", yaxt="n", ylab=NA, xlab=NA, 
     main="AUC - differences in mean and dispersion", cex.main=3)
boxplot(cbind(DEDD.results.DEDD2$auc, 
              DEDD.results.DEDD5$auc, 
              DEDD.results.DEDD10$auc, 
              DEDD.results.DEDD20$auc, 
              DEDD.results.DEDD50$auc), 
        col=col_vector[1:3], border=col_vector[1:3], pch=20, cex.axis=3, xaxt='n', 
        at=positions + offsets, lwd=0.2, cex=0.5, 
        add=T)
axis(side=1, at=c(2,6,10,14,18), labels=c(2,5,10,20,50), tick=F, cex.axis=3, mgp=c(3,1.5,0))
legend("bottomright", fill=col_vector[1:3], border=col_vector[1:3], bty='n', cex=2.5, ncol=1, 
       legend=c("diffVar", "Hybrid", "HMM"))

rbind(mean(DEDD.results.DE2$auc), 
      mean(DEDD.results.DE5$auc), 
      mean(DEDD.results.DE10$auc), 
      mean(DEDD.results.DE20$auc), 
      mean(DEDD.results.DE50$auc))
rbind(mean(DEDD.results.DD2$auc), 
      mean(DEDD.results.DD5$auc), 
      mean(DEDD.results.DD10$auc), 
      mean(DEDD.results.DD20$auc), 
      mean(DEDD.results.DD50$auc))
rbind(colMeans(DEDD.results.DEDD2$auc), 
      colMeans(DEDD.results.DEDD5$auc), 
      colMeans(DEDD.results.DEDD10$auc), 
      colMeans(DEDD.results.DEDD20$auc), 
      colMeans(DEDD.results.DEDD50$auc))


### pAUC ####
positions = rep(1:15)
offsets <- c(rep(0, 3), rep(1, 3), rep(2, 3), rep(3, 3), rep(4, 3))

par(mfrow=c(3,1), mar=c(3.5,3.5,3.5,0.5), mgp=c(3,1,0))

plot(NA, xlim=c(0.7,5.3), ylim=c(0.001,0.048), xaxt="n", yaxt="n", ylab=NA, xlab=NA, 
     main="Partial AUC - differences in mean only", cex.main=3)
boxplot(cbind(DEDD.results.DE5$pauc, 
              DEDD.results.DE5$pauc, 
              DEDD.results.DE10$pauc, 
              DEDD.results.DE20$pauc, 
              DEDD.results.DE50$pauc), 
        col=col_vector[3], border=col_vector[3], pch=20, cex.axis=3, xaxt='n', 
        lwd=0.2, cex=0.5, boxwex=0.6, 
        add=T)
axis(side=1, at=1:5, labels=c(2,5,10,20,50), tick=F, cex.axis=3, mgp=c(3,1.5,0))

plot(NA, xlim=c(0.7,5.3), ylim=c(0.001,0.048), xaxt="n", yaxt="n", ylab=NA, xlab=NA, 
     main="Partial AUC - differences in dispersion only", cex.main=3)
boxplot(cbind(DEDD.results.DD2$pauc, 
              DEDD.results.DD5$pauc, 
              DEDD.results.DD10$pauc, 
              DEDD.results.DD20$pauc, 
              DEDD.results.DD50$pauc), 
        col=col_vector[3], border=col_vector[3], pch=20, cex.axis=3, xaxt='n', 
        lwd=0.2, cex=0.5, boxwex=0.6, 
        add=T)
axis(side=1, at=1:5, labels=c(2,5,10,20,50), tick=F, cex.axis=3, mgp=c(3,1.5,0))

plot(NA, xlim=c(1,19), ylim=c(0.001,0.048), xaxt="n", yaxt="n", ylab=NA, xlab=NA, 
     main="Partial AUC - differences in mean and dispersion", cex.main=3)
boxplot(cbind(DEDD.results.DEDD2$pauc, 
              DEDD.results.DEDD5$pauc, 
              DEDD.results.DEDD10$pauc, 
              DEDD.results.DEDD20$pauc, 
              DEDD.results.DEDD50$pauc), 
        col=col_vector[1:3], border=col_vector[1:3], pch=20, cex.axis=3, xaxt='n', 
        at=positions + offsets, lwd=0.2, cex=0.5, 
        add=T)
axis(side=1, at=c(2,6,10,14,18), labels=c(2,5,10,20,50), tick=F, cex.axis=3, mgp=c(3,1.5,0))
legend("bottomright", fill=col_vector[1:3], border=col_vector[1:3], bty='n', cex=2.5, ncol=1, 
       legend=c("diffVar", "Hybrid", "HMM"))

rbind(mean(DEDD.results.DE2$pauc), 
      mean(DEDD.results.DE5$pauc), 
      mean(DEDD.results.DE10$pauc), 
      mean(DEDD.results.DE20$pauc), 
      mean(DEDD.results.DE50$pauc))
rbind(mean(DEDD.results.DD2$pauc), 
      mean(DEDD.results.DD5$pauc), 
      mean(DEDD.results.DD10$pauc), 
      mean(DEDD.results.DD20$pauc), 
      mean(DEDD.results.DD50$pauc))
rbind(colMeans(DEDD.results.DEDD2$pauc), 
      colMeans(DEDD.results.DEDD5$pauc), 
      colMeans(DEDD.results.DEDD10$pauc), 
      colMeans(DEDD.results.DEDD20$pauc), 
      colMeans(DEDD.results.DEDD50$pauc))


### FDR ####
par(mfrow=c(3,1), mar=c(3.5,3.5,3.5,0.5), mgp=c(3,1,0))

positions = rep(1:10)
offsets <- c(rep(0, 2), rep(0.5, 2), rep(1, 2), rep(1.5, 2), rep(2, 2))
plot(NA, xlim=c(0.75,12.25), ylim=c(0,1), xaxt="n", yaxt="n", ylab=NA, xlab=NA, 
     main="FDR - differences in mean only", cex.main=3)
lines(c(0.5,2.5), c(0.05,0.05), col="lightgrey")
lines(c(3,5), c(0.05,0.05), col="lightgrey")
lines(c(5.5,7.5), c(0.05,0.05), col="lightgrey")
lines(c(8,10), c(0.05,0.05), col="lightgrey")
lines(c(10.5,12.5), c(0.05,0.05), col="lightgrey")
boxplot(cbind(DEDD.results.DE2$fdr, 
              DEDD.results.DE5$fdr, 
              DEDD.results.DE10$fdr, 
              DEDD.results.DE20$fdr, 
              DEDD.results.DE50$fdr), 
        col=col_vector[3:4], border=col_vector[3:4], pch=20, cex.axis=3, xaxt='n', 
        at=positions + offsets, lwd=0.2, cex=0.5, 
        add=T)
axis(side=1, at=c(1.5,4,6.5,9,11.5), labels=c(2,5,10,20,50), tick=F, cex.axis=3, mgp=c(3,1.5,0))
legend("topright", fill=col_vector[3:4], border=col_vector[3:4], bty='n', cex=2.5, ncol=1, 
       legend=c("Posterior threshold", "BFDR"))

plot(NA, xlim=c(0.75,12.25), ylim=c(0,1), xaxt="n", yaxt="n", ylab=NA, xlab=NA, 
     main="FDR - differences in dispersion only", cex.main=3)
lines(c(0.5,2.5), c(0.05,0.05), col="lightgrey")
lines(c(3,5), c(0.05,0.05), col="lightgrey")
lines(c(5.5,7.5), c(0.05,0.05), col="lightgrey")
lines(c(8,10), c(0.05,0.05), col="lightgrey")
lines(c(10.5,12.5), c(0.05,0.05), col="lightgrey")
boxplot(cbind(DEDD.results.DD2$fdr, 
              DEDD.results.DD5$fdr, 
              DEDD.results.DD10$fdr, 
              DEDD.results.DD20$fdr, 
              DEDD.results.DD50$fdr), 
        col=col_vector[3:4], border=col_vector[3:4], pch=20, cex.axis=3, xaxt='n', 
        at=positions + offsets, lwd=0.2, cex=0.5, 
        add=T)
axis(side=1, at=c(1.5,4,6.5,9,11.5), labels=c(2,5,10,20,50), tick=F, cex.axis=3, mgp=c(3,1.5,0))
legend("topright", fill=col_vector[3:4], border=col_vector[3:4], bty='n', cex=2.5, ncol=1, 
       legend=c("Posterior threshold", "BFDR"))

positions = rep(1:20)
offsets <- c(rep(0, 4), rep(1, 4), rep(2, 4), rep(3, 4), rep(4, 4))
plot(NA, xlim=c(1,24), ylim=c(0,1), xaxt="n", yaxt="n", ylab=NA, xlab=NA, 
     main="FDR - differences in mean and dispersion", cex.main=3)
lines(c(0.5,4.5), c(0.05,0.05), col="lightgrey")
lines(c(5.5,9.5), c(0.05,0.05), col="lightgrey")
lines(c(10.5,14.5), c(0.05,0.05), col="lightgrey")
lines(c(15.5,19.5), c(0.05,0.05), col="lightgrey")
lines(c(20.5,24.5), c(0.05,0.05), col="lightgrey")
boxplot(cbind(DEDD.results.DEDD2$fdr, 
              DEDD.results.DEDD5$fdr, 
              DEDD.results.DEDD10$fdr, 
              DEDD.results.DEDD20$fdr, 
              DEDD.results.DEDD50$fdr), 
        col=col_vector[1:4], border=col_vector[1:4], pch=20, cex.axis=3, xaxt='n', 
        at=positions + offsets, lwd=0.2, cex=0.5, 
        add=T)
axis(side=1, at=c(2.5,7.5,12.5,17.5,22.5), labels=c(2,5,10,20,50), tick=F, cex.axis=3, mgp=c(3,1.5,0))
legend("topright", fill=col_vector[1:4], border=col_vector[1:4], bty='n', cex=2.5, ncol=1, 
       legend=c("diffVar", "Hybrid", "HMM, posterior threshold", "HMM, BFDR"))

rbind(colMeans(DEDD.results.DE2$fdr, na.rm=T), 
      colMeans(DEDD.results.DE5$fdr, na.rm=T), 
      colMeans(DEDD.results.DE10$fdr, na.rm=T), 
      colMeans(DEDD.results.DE20$fdr, na.rm=T), 
      colMeans(DEDD.results.DE50$fdr, na.rm=T))
rbind(colMeans(DEDD.results.DD2$fdr, na.rm=T), 
      colMeans(DEDD.results.DD5$fdr, na.rm=T), 
      colMeans(DEDD.results.DD10$fdr, na.rm=T), 
      colMeans(DEDD.results.DD20$fdr, na.rm=T), 
      colMeans(DEDD.results.DD50$fdr, na.rm=T))
rbind(colMeans(DEDD.results.DEDD2$fdr, na.rm=T), 
      colMeans(DEDD.results.DEDD5$fdr, na.rm=T), 
      colMeans(DEDD.results.DEDD10$fdr, na.rm=T), 
      colMeans(DEDD.results.DEDD20$fdr, na.rm=T), 
      colMeans(DEDD.results.DEDD50$fdr, na.rm=T))


### TPR ####
par(mfrow=c(3,1), mar=c(3.5,3.5,3.5,0.5), mgp=c(3,1,0))

positions = rep(1:10)
offsets <- c(rep(0, 2), rep(0.5, 2), rep(1, 2), rep(1.5, 2), rep(2, 2))
plot(NA, xlim=c(0.75,12.25), ylim=c(0,0.92), xaxt="n", yaxt="n", ylab=NA, xlab=NA, 
     main="TPR - differences in mean only", cex.main=3)
boxplot(cbind(DEDD.results.DE2$tpr, 
              DEDD.results.DE5$tpr, 
              DEDD.results.DE10$tpr, 
              DEDD.results.DE20$tpr, 
              DEDD.results.DE50$tpr), 
        col=col_vector[3:4], border=col_vector[3:4], pch=20, cex.axis=3, xaxt='n', 
        at=positions + offsets, lwd=0.2, cex=0.5, 
        add=T)
axis(side=1, at=c(1.5,4,6.5,9,11.5), labels=c(2,5,10,20,50), tick=F, cex.axis=3, mgp=c(3,1.5,0))
legend("topleft", fill=col_vector[3:4],border=col_vector[3:4],  bty='n', cex=2.5, ncol=1, 
       legend=c("Posterior threshold", "BFDR"))

plot(NA, xlim=c(0.75,12.25), ylim=c(0,0.92), xaxt="n", yaxt="n", ylab=NA, xlab=NA, 
     main="TPR - differences in dispersion only", cex.main=3)
boxplot(cbind(DEDD.results.DD2$tpr, 
              DEDD.results.DD5$tpr, 
              DEDD.results.DD10$tpr, 
              DEDD.results.DD20$tpr, 
              DEDD.results.DD50$tpr), 
        col=col_vector[3:4], border=col_vector[3:4], pch=20, cex.axis=3, xaxt='n', 
        at=positions + offsets, lwd=0.2, cex=0.5, 
        add=T)
axis(side=1, at=c(1.5,4,6.5,9,11.5), labels=c(2,5,10,20,50), tick=F, cex.axis=3, mgp=c(3,1.5,0))
legend("topleft", fill=col_vector[3:4], border=col_vector[3:4], bty='n', cex=2.5, ncol=1, 
       legend=c("Posterior threshold", "BFDR"))

positions = rep(1:20)
offsets <- c(rep(0, 4), rep(1, 4), rep(2, 4), rep(3, 4), rep(4, 4))
plot(NA, xlim=c(1,24), ylim=c(0,0.92), xaxt="n", yaxt="n", ylab=NA, xlab=NA, 
     main="TPR - differences in mean and dispersion", cex.main=3)
boxplot(cbind(DEDD.results.DEDD2$tpr, 
              DEDD.results.DEDD5$tpr, 
              DEDD.results.DEDD10$tpr, 
              DEDD.results.DEDD20$tpr, 
              DEDD.results.DEDD50$tpr), 
        col=col_vector[1:4], border=col_vector[1:4], pch=20, cex.axis=3, xaxt='n', 
        at=positions + offsets, lwd=0.2, cex=0.5, 
        add=T)
axis(side=1, at=c(2.5,7.5,12.5,17.5,22.5), labels=c(2,5,10,20,50), tick=F, cex.axis=3, mgp=c(3,1.5,0))
legend("topleft", fill=col_vector[1:4], border=col_vector[1:4], bty='n', cex=2.5, ncol=1, 
       legend=c("diffVar", "Hybrid", "HMM, posterior threshold", "HMM, BFDR"))

rbind(colMeans(DEDD.results.DE2$tpr, na.rm=T), 
      colMeans(DEDD.results.DE5$tpr, na.rm=T), 
      colMeans(DEDD.results.DE10$tpr, na.rm=T), 
      colMeans(DEDD.results.DE20$tpr, na.rm=T), 
      colMeans(DEDD.results.DE50$tpr, na.rm=T))
rbind(colMeans(DEDD.results.DD2$tpr, na.rm=T), 
      colMeans(DEDD.results.DD5$tpr, na.rm=T), 
      colMeans(DEDD.results.DD10$tpr, na.rm=T), 
      colMeans(DEDD.results.DD20$tpr, na.rm=T), 
      colMeans(DEDD.results.DD50$tpr, na.rm=T))
rbind(colMeans(DEDD.results.DEDD2$tpr, na.rm=T), 
      colMeans(DEDD.results.DEDD5$tpr, na.rm=T), 
      colMeans(DEDD.results.DEDD10$tpr, na.rm=T), 
      colMeans(DEDD.results.DEDD20$tpr, na.rm=T), 
      colMeans(DEDD.results.DEDD50$tpr, na.rm=T))


### FDR curves ####
par(mfrow=c(3,5), mar=c(1.5,3.5,2.5,1), mgp=c(3,1,0))

plot(DEDD.results.DE2$mean.discoveries, DEDD.results.DE2$mean.fdr, 
     type='l', ylim=c(0,0.95), xlim=c(0,2000), lwd=2, col=col_vector[3], cex.axis=3, 
     xlab=NA, ylab=NA, xaxt='n', main=paste0("2"), cex.main=3)
axis(side=1, at=c(0,500,1000,1500,2000), tick=T, cex.axis=3, mgp=c(3,2,0))
text("Differences in\nmean only", x=0, y=0.2, cex=2.5, adj=0)
abline(h=0.05, col='lightgrey')
plot(DEDD.results.DE5$mean.discoveries, DEDD.results.DE5$mean.fdr, 
     type='l', ylim=c(0,0.95), xlim=c(0,2000), lwd=2, col=col_vector[3], cex.axis=3, 
     xlab=NA, ylab=NA, xaxt='n', yaxt='n', main=paste0("5"), cex.main=3)
axis(side=1, at=c(0,500,1000,1500,2000), tick=T, cex.axis=3, mgp=c(3,2,0))
abline(h=0.05, col='lightgrey')
plot(DEDD.results.DE10$mean.discoveries, DEDD.results.DE10$mean.fdr, 
     type='l', ylim=c(0,0.95), xlim=c(0,2000), lwd=2, col=col_vector[3], cex.axis=3, 
     xlab=NA, ylab=NA, xaxt='n', yaxt='n', main=paste0("10"), cex.main=3)
axis(side=1, at=c(0,500,1000,1500,2000), tick=T, cex.axis=3, mgp=c(3,2,0))
abline(h=0.05, col='lightgrey')
plot(DEDD.results.DE20$mean.discoveries, DEDD.results.DE20$mean.fdr, 
     type='l', ylim=c(0,0.95), xlim=c(0,2000), lwd=2, col=col_vector[3], cex.axis=3, 
     xlab=NA, ylab=NA, xaxt='n', yaxt='n', main=paste0("20"), cex.main=3)
axis(side=1, at=c(0,500,1000,1500,2000), tick=T, cex.axis=3, mgp=c(3,2,0))
abline(h=0.05, col='lightgrey')
plot(DEDD.results.DE50$mean.discoveries, DEDD.results.DE50$mean.fdr, 
     type='l', ylim=c(0,0.95), xlim=c(0,2000), lwd=2, col=col_vector[3], cex.axis=3, 
     xlab=NA, ylab=NA, xaxt='n', yaxt='n', main=paste0("50"), cex.main=3)
axis(side=1, at=c(0,500,1000,1500,2000), tick=T, cex.axis=3, mgp=c(3,2,0))
abline(h=0.05, col='lightgrey')

plot(DEDD.results.DD2$mean.discoveries, DEDD.results.DD2$mean.fdr, 
     type='l', ylim=c(0,0.95), xlim=c(0,2000), lwd=2, col=col_vector[3], cex.axis=3, 
     xlab=NA, ylab=NA, xaxt='n')
text("Differences in\ndispersion only", x=0, y=0.2, cex=2.5, adj=0)
abline(h=0.05, col='lightgrey')
plot(DEDD.results.DD5$mean.discoveries, DEDD.results.DD5$mean.fdr, 
     type='l', ylim=c(0,0.95), xlim=c(0,2000), lwd=2, col=col_vector[3], cex.axis=3, 
     xlab=NA, ylab=NA, xaxt='n', yaxt='n')
abline(h=0.05, col='lightgrey')
plot(DEDD.results.DD10$mean.discoveries, DEDD.results.DD10$mean.fdr, 
     type='l', ylim=c(0,0.95), xlim=c(0,2000), lwd=2, col=col_vector[3], cex.axis=3, 
     xlab=NA, ylab=NA, xaxt='n', yaxt='n')
abline(h=0.05, col='lightgrey')
plot(DEDD.results.DD20$mean.discoveries, DEDD.results.DD20$mean.fdr, 
     type='l', ylim=c(0,0.95), xlim=c(0,2000), lwd=2, col=col_vector[3], cex.axis=3, 
     xlab=NA, ylab=NA, xaxt='n', yaxt='n')
abline(h=0.05, col='lightgrey')
plot(DEDD.results.DD50$mean.discoveries, DEDD.results.DD50$mean.fdr, 
     type='l', ylim=c(0,0.95), xlim=c(0,2000), lwd=2, col=col_vector[3], cex.axis=3, 
     xlab=NA, ylab=NA, xaxt='n', yaxt='n')
abline(h=0.05, col='lightgrey')

plot(DEDD.results.DEDD2$mean.discoveries$dV.tmm, DEDD.results.DEDD2$mean.fdr$dV.tmm, 
     type='l', ylim=c(0,0.95), xlim=c(0,2000), lwd=2, col=col_vector[1], cex.axis=3, 
     xlab=NA, ylab=NA, xaxt='n')
text("Differences in\nmean and dispersion", x=0, y=0.8, cex=2.5, adj=0)
lines(DEDD.results.DEDD2$mean.discoveries$edgeR_HM.tmm, DEDD.results.DEDD2$mean.fdr$edgeR_HM.tmm, 
      lwd=2, col=col_vector[2], cex.axis=3)
lines(DEDD.results.DEDD2$mean.discoveries$lnHMM.tmm, DEDD.results.DEDD2$mean.fdr$lnHMM.tmm, 
      lwd=2, col=col_vector[3], cex.axis=3)
abline(h=0.05, col='lightgrey')
plot(DEDD.results.DEDD5$mean.discoveries$dV.tmm, DEDD.results.DEDD5$mean.fdr$dV.tmm, 
     type='l', ylim=c(0,0.95), xlim=c(0,2000), lwd=2, col=col_vector[1], cex.axis=3, 
     xlab=NA, ylab=NA, xaxt='n', yaxt='n')
lines(DEDD.results.DEDD5$mean.discoveries$edgeR_HM.tmm, DEDD.results.DEDD5$mean.fdr$edgeR_HM.tmm, 
      lwd=2, col=col_vector[2], cex.axis=3)
lines(DEDD.results.DEDD5$mean.discoveries$lnHMM.tmm, DEDD.results.DEDD5$mean.fdr$lnHMM.tmm, 
      lwd=2, col=col_vector[3], cex.axis=3)
abline(h=0.05, col='lightgrey')
plot(DEDD.results.DEDD10$mean.discoveries$dV.tmm, DEDD.results.DEDD10$mean.fdr$dV.tmm, 
     type='l', ylim=c(0,0.95), xlim=c(0,2000), lwd=2, col=col_vector[1], cex.axis=3, 
     xlab=NA, ylab=NA, xaxt='n', yaxt='n')
lines(DEDD.results.DEDD10$mean.discoveries$edgeR_HM.tmm, DEDD.results.DEDD10$mean.fdr$edgeR_HM.tmm, 
      lwd=2, col=col_vector[2], cex.axis=3)
lines(DEDD.results.DEDD10$mean.discoveries$lnHMM.tmm, DEDD.results.DEDD10$mean.fdr$lnHMM.tmm, 
      lwd=2, col=col_vector[3], cex.axis=3)
abline(h=0.05, col='lightgrey')
plot(DEDD.results.DEDD20$mean.discoveries$dV.tmm, DEDD.results.DEDD20$mean.fdr$dV.tmm, 
     type='l', ylim=c(0,0.95), xlim=c(0,2000), lwd=2, col=col_vector[1], cex.axis=3, 
     xlab=NA, ylab=NA, xaxt='n', yaxt='n')
lines(DEDD.results.DEDD20$mean.discoveries$edgeR_HM.tmm, DEDD.results.DEDD20$mean.fdr$edgeR_HM.tmm, 
      lwd=2, col=col_vector[2], cex.axis=3)
lines(DEDD.results.DEDD20$mean.discoveries$lnHMM.tmm, DEDD.results.DEDD20$mean.fdr$lnHMM.tmm, 
      lwd=2, col=col_vector[3], cex.axis=3)
abline(h=0.05, col='lightgrey')
plot(DEDD.results.DEDD50$mean.discoveries$dV.tmm, DEDD.results.DEDD50$mean.fdr$dV.tmm, 
     type='l', ylim=c(0,0.95), xlim=c(0,2000), lwd=2, col=col_vector[1], cex.axis=3, 
     xlab=NA, ylab=NA, xaxt='n', yaxt='n')
lines(DEDD.results.DEDD50$mean.discoveries$edgeR_HM.tmm, DEDD.results.DEDD50$mean.fdr$edgeR_HM.tmm, 
      lwd=2, col=col_vector[2], cex.axis=3)
lines(DEDD.results.DEDD50$mean.discoveries$lnHMM.tmm, DEDD.results.DEDD50$mean.fdr$lnHMM.tmm, 
      lwd=2, col=col_vector[3], cex.axis=3)
abline(h=0.05, col='lightgrey')

legend("topleft", col=col_vector[1:3], bty='n', cex=2.5, ncol=1, lty=1, lwd=2, 
       legend=c("diffVar", "Hybrid", "HMM"))
