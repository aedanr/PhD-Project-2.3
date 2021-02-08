library(here)
library(RColorBrewer)

## Load data, create colour vector ####
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
  DEDD.results.DEDD2$fdr <- DEDD.results.DEDD2$fdr[
    , c("dV.tmm", "edgeR_HM.tmm", "thr.lnHMM.tmm", "bfdr.lnHMM.tmm")]
  DEDD.results.DEDD2$tpr <- DEDD.results.DEDD2$tpr[
    , c("dV.tmm", "edgeR_HM.tmm", "thr.lnHMM.tmm", "bfdr.lnHMM.tmm")]
  DEDD.results.DEDD5$fdr <- DEDD.results.DEDD5$fdr[
    , c("dV.tmm", "edgeR_HM.tmm", "thr.lnHMM.tmm", "bfdr.lnHMM.tmm")]
  DEDD.results.DEDD5$tpr <- DEDD.results.DEDD5$tpr[
    , c("dV.tmm", "edgeR_HM.tmm", "thr.lnHMM.tmm", "bfdr.lnHMM.tmm")]
  DEDD.results.DEDD10$fdr <- DEDD.results.DEDD10$fdr[
    , c("dV.tmm", "edgeR_HM.tmm", "thr.lnHMM.tmm", "bfdr.lnHMM.tmm")]
  DEDD.results.DEDD10$tpr <- DEDD.results.DEDD10$tpr[
    , c("dV.tmm", "edgeR_HM.tmm", "thr.lnHMM.tmm", "bfdr.lnHMM.tmm")]
  DEDD.results.DEDD20$fdr <- DEDD.results.DEDD20$fdr[
    , c("dV.tmm", "edgeR_HM.tmm", "thr.lnHMM.tmm", "bfdr.lnHMM.tmm")]
  DEDD.results.DEDD20$tpr <- DEDD.results.DEDD20$tpr[
    , c("dV.tmm", "edgeR_HM.tmm", "thr.lnHMM.tmm", "bfdr.lnHMM.tmm")]
  DEDD.results.DEDD50$fdr <- DEDD.results.DEDD50$fdr[
    , c("dV.tmm", "edgeR_HM.tmm", "thr.lnHMM.tmm", "bfdr.lnHMM.tmm")]
  DEDD.results.DEDD50$tpr <- DEDD.results.DEDD50$tpr[
    , c("dV.tmm", "edgeR_HM.tmm", "thr.lnHMM.tmm", "bfdr.lnHMM.tmm")]
}


## Plot FDR and TPR ####
positions = rep(1:20)
offsets <- c(rep(0, 4), rep(2, 4), rep(4, 4), rep(6, 4), rep(8, 4))

par(mfrow=c(2,1), mar=c(2,2,2,0.5), mgp=c(3,0.7,0))
plot(NA, xlim=c(1,28), ylim=c(0,1), xaxt="n", yaxt="n", 
     main="False discovery rate", cex.main=1.5)
lines(c(0.5,4.5), c(0.05,0.05), col="lightgrey")
lines(c(6.5,10.5), c(0.05,0.05), col="lightgrey")
lines(c(12.5,16.5), c(0.05,0.05), col="lightgrey")
lines(c(18.5,22.5), c(0.05,0.05), col="lightgrey")
lines(c(24.5,28.5), c(0.05,0.05), col="lightgrey")
boxplot(cbind(DEDD.results.DEDD2$fdr, 
              DEDD.results.DEDD5$fdr, 
              DEDD.results.DEDD10$fdr, 
              DEDD.results.DEDD20$fdr, 
              DEDD.results.DEDD50$fdr), 
        col=col_vector[1:4], pch=20, cex.axis=1.5, xaxt='n',
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(2.5,8.5,14.5,20.5,26.5), labels=c(2,5,10,20,50), tick=F, cex.axis=1.5)
legend("topright", fill=col_vector[1:4], bty='n', cex=1.5, ncol=1, 
       legend=c("diffVar", "Hybrid", "HMM, posterior threshold", "HMM, BFDR"))
plot(NA, xlim=c(1,28), ylim=c(0,1), xaxt="n", yaxt="n", 
     main="Sensitivity", cex.main=1.5)
boxplot(cbind(DEDD.results.DEDD2$tpr, 
              DEDD.results.DEDD5$tpr, 
              DEDD.results.DEDD10$tpr, 
              DEDD.results.DEDD20$tpr, 
              DEDD.results.DEDD50$tpr), 
        col=col_vector[1:4], pch=20, cex.axis=1.5, xaxt='n', 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(2.5,8.5,14.5,20.5,26.5), labels=c(2,5,10,20,50), tick=F, cex.axis=1.5)
