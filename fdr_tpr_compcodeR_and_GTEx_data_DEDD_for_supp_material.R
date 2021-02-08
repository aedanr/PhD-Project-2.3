library(here)
library(RColorBrewer)

# Create colour vector
qual_col_pals = brewer.pal.info[brewer.pal.info$category == "qual" & 
                                  brewer.pal.info$colorblind == T,]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector <- col_vector[c(1:3,9,11,26)]
rm(qual_col_pals)

## Load compcodeR data, keep only diffVar, edgeR/lnHM log, lnHMM (posterior threshold and BFDR) ####
folder <- "Results/compcodeR combined results Dec 2020 including diffVar"
for (i in c("DEDD2", "DEDD5", "DEDD10", "DEDD20", "DEDD50")) {
  assign(paste0("DEDD.results.compcodeR.", i), 
         readRDS(here(folder, paste0("DEDD.results.", i, ".rds"))))
}
rm(i,folder)

{
  DEDD.results.compcodeR.DEDD2$fdr <- DEDD.results.compcodeR.DEDD2$fdr[
    , c("dV.tmm", "edgeR_HM.tmm", "thr.lnHMM.tmm", "bfdr.lnHMM.tmm")]
  DEDD.results.compcodeR.DEDD2$tpr <- DEDD.results.compcodeR.DEDD2$tpr[
    , c("dV.tmm", "edgeR_HM.tmm", "thr.lnHMM.tmm", "bfdr.lnHMM.tmm")]
  DEDD.results.compcodeR.DEDD5$fdr <- DEDD.results.compcodeR.DEDD5$fdr[
    , c("dV.tmm", "edgeR_HM.tmm", "thr.lnHMM.tmm", "bfdr.lnHMM.tmm")]
  DEDD.results.compcodeR.DEDD5$tpr <- DEDD.results.compcodeR.DEDD5$tpr[
    , c("dV.tmm", "edgeR_HM.tmm", "thr.lnHMM.tmm", "bfdr.lnHMM.tmm")]
  DEDD.results.compcodeR.DEDD10$fdr <- DEDD.results.compcodeR.DEDD10$fdr[
    , c("dV.tmm", "edgeR_HM.tmm", "thr.lnHMM.tmm", "bfdr.lnHMM.tmm")]
  DEDD.results.compcodeR.DEDD10$tpr <- DEDD.results.compcodeR.DEDD10$tpr[
    , c("dV.tmm", "edgeR_HM.tmm", "thr.lnHMM.tmm", "bfdr.lnHMM.tmm")]
  DEDD.results.compcodeR.DEDD20$fdr <- DEDD.results.compcodeR.DEDD20$fdr[
    , c("dV.tmm", "edgeR_HM.tmm", "thr.lnHMM.tmm", "bfdr.lnHMM.tmm")]
  DEDD.results.compcodeR.DEDD20$tpr <- DEDD.results.compcodeR.DEDD20$tpr[
    , c("dV.tmm", "edgeR_HM.tmm", "thr.lnHMM.tmm", "bfdr.lnHMM.tmm")]
  DEDD.results.compcodeR.DEDD50$fdr <- DEDD.results.compcodeR.DEDD50$fdr[
    , c("dV.tmm", "edgeR_HM.tmm", "thr.lnHMM.tmm", "bfdr.lnHMM.tmm")]
  DEDD.results.compcodeR.DEDD50$tpr <- DEDD.results.compcodeR.DEDD50$tpr[
    , c("dV.tmm", "edgeR_HM.tmm", "thr.lnHMM.tmm", "bfdr.lnHMM.tmm")]
}

## Load GTEx data, keep only diffVar, voom/lnHM log, lnHMM (posterior threshold and BFDR) ####
folder <- "Results/GTEx combined results Sept 2020"
for (i in c("2", "5", "10", "20", "50")) {
  assign(paste0("DEDD.results.GTEx.DEDD", i), 
         readRDS(here(folder, paste0("DEDD.results.DEDD", i, ".rds"))))
}
rm(i,folder)

for (metric in c("fdr", "tpr")) {
  for (size in c("2", "5", "10", "20", "50")) {
    assign(
      paste0(metric, "_blood_", size), 
      get(paste0("DEDD.results.GTEx.DEDD", size))[[metric]][
        , c("blood_dV", "blood_voom_HM", "blood_lnHMM.thr", "blood_lnHMM.bfdr")
      ]
    )
    assign(
      paste0(metric, "_muscle_", size), 
      get(paste0("DEDD.results.GTEx.DEDD", size))[[metric]][
        , c("muscle_dV", "muscle_voom_HM", "muscle_lnHMM.thr", "muscle_lnHMM.bfdr")
      ]
    )
  }
}
rm(metric, size)


## Plot FDR and TPR for each dataset ####
positions = rep(1:20)
offsets <- c(rep(0, 4), rep(2, 4), rep(4, 4), rep(6, 4), rep(8, 4))

par(mfrow=c(3,2), mar=c(2,2.5,2.5,0.5), mgp=c(3,0.7,0))

# First row - compcodeR
plot(NA, xlim=c(1,28), ylim=c(0,1), xaxt="n", yaxt="n", 
     main="False discovery rate", cex.main=2.5)
lines(c(0.5,4.5), c(0.05,0.05), col="lightgrey")
lines(c(6.5,10.5), c(0.05,0.05), col="lightgrey")
lines(c(12.5,16.5), c(0.05,0.05), col="lightgrey")
lines(c(18.5,22.5), c(0.05,0.05), col="lightgrey")
lines(c(24.5,28.5), c(0.05,0.05), col="lightgrey")
boxplot(cbind(DEDD.results.compcodeR.DEDD2$fdr, 
              DEDD.results.compcodeR.DEDD5$fdr, 
              DEDD.results.compcodeR.DEDD10$fdr, 
              DEDD.results.compcodeR.DEDD20$fdr, 
              DEDD.results.compcodeR.DEDD50$fdr), 
        col=col_vector[1:4], pch=20, cex.axis=2, xaxt='n', 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(2.5,8.5,14.5,20.5,26.5), labels=c(2,5,10,20,50), tick=F, cex.axis=2.5, mgp=c(3,1,0))
mtext('Simulated data', cex=1.5, side=3, line=-2)
legend("topright", fill=col_vector[1:6], bty='n', cex=2, ncol=1, 
       legend=c("diffVar", "Hybrid", "HMM, posterior threshold", "HMM, BFDR"))
plot(NA, xlim=c(1,28), ylim=c(0,1), xaxt="n", yaxt="n", 
     main="Sensitivity", cex.main=2.5)
boxplot(cbind(DEDD.results.compcodeR.DEDD2$tpr, 
              DEDD.results.compcodeR.DEDD5$tpr, 
              DEDD.results.compcodeR.DEDD10$tpr, 
              DEDD.results.compcodeR.DEDD20$tpr, 
              DEDD.results.compcodeR.DEDD50$tpr), 
        col=col_vector[1:4], pch=20, cex.axis=2, xaxt='n', 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(2.5,8.5,14.5,20.5,26.5), labels=c(2,5,10,20,50), tick=F, cex.axis=2.5, mgp=c(3,1,0))
mtext('Simulated data', cex=1.5, side=3, line=-2)

# Second row - blood
plot(NA, xlim=c(1,28), ylim=c(0,1), xaxt="n", yaxt="n")
lines(c(0.5,4.5), c(0.05,0.05), col="lightgrey")
lines(c(6.5,10.5), c(0.05,0.05), col="lightgrey")
lines(c(12.5,16.5), c(0.05,0.05), col="lightgrey")
lines(c(18.5,22.5), c(0.05,0.05), col="lightgrey")
lines(c(24.5,28.5), c(0.05,0.05), col="lightgrey")
boxplot(cbind(fdr_blood_2, fdr_blood_5, fdr_blood_10, fdr_blood_20, fdr_blood_50), 
        col=col_vector[1:4], pch=20, cex.axis=2, xaxt='n', 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(2.5,8.5,14.5,20.5,26.5), labels=c(2,5,10,20,50), tick=F, cex.axis=2.5, mgp=c(3,1,0))
mtext('GTEx blood data', cex=1.5, side=3, line=-2)
plot(NA, xlim=c(1,28), ylim=c(0,1), xaxt="n", yaxt="n")
boxplot(cbind(tpr_blood_2, tpr_blood_5, tpr_blood_10, tpr_blood_20, tpr_blood_50), 
        col=col_vector[1:4], pch=20, cex.axis=2, xaxt='n', 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(2.5,8.5,14.5,20.5,26.5), labels=c(2,5,10,20,50), tick=F, cex.axis=2.5, mgp=c(3,1,0))
mtext('GTEx blood data', cex=1.5, side=3, line=-2)

# Third row - muscle
plot(NA, xlim=c(1,28), ylim=c(0,1), xaxt="n", yaxt="n")
lines(c(0.5,4.5), c(0.05,0.05), col="lightgrey")
lines(c(6.5,10.5), c(0.05,0.05), col="lightgrey")
lines(c(12.5,16.5), c(0.05,0.05), col="lightgrey")
lines(c(18.5,22.5), c(0.05,0.05), col="lightgrey")
lines(c(24.5,28.5), c(0.05,0.05), col="lightgrey")
boxplot(cbind(fdr_muscle_2, fdr_muscle_5, fdr_muscle_10, fdr_muscle_20, fdr_muscle_50), 
        col=col_vector[1:4], pch=20, cex.axis=2, xaxt='n', 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(2.5,8.5,14.5,20.5,26.5), labels=c(2,5,10,20,50), tick=F, cex.axis=2.5, mgp=c(3,1,0))
mtext('GTEx muscle data', cex=1.5, side=3, line=-2)
plot(NA, xlim=c(1,28), ylim=c(0,1), xaxt="n", yaxt="n")
boxplot(cbind(tpr_muscle_2, tpr_muscle_5, tpr_muscle_10, tpr_muscle_20, tpr_muscle_50), 
        col=col_vector[1:4], pch=20, cex.axis=2, xaxt='n', 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(2.5,8.5,14.5,20.5,26.5), labels=c(2,5,10,20,50), tick=F, cex.axis=2.5, mgp=c(3,1,0))
mtext('GTEx muscle data', cex=1.5, side=3, line=-2)


## Only 5 and 50 samples per group to match FDR curves (haven't decided whether to use) ####
positions = rep(1:8)
offsets <- c(rep(0, 4), rep(1, 4))
par(mfrow=c(3,2), mar=c(2,2.5,2.5,0.5), mgp=c(3,0.7,0))

# First row - compcodeR
plot(NA, xlim=c(0,10), ylim=c(0,1), xaxt="n", yaxt="n", 
     main="False discovery rate", cex.main=2.5)
lines(c(0.5,4.5), c(0.05,0.05), col="lightgrey")
lines(c(5.5,9.5), c(0.05,0.05), col="lightgrey")
boxplot(cbind(DEDD.results.compcodeR.DEDD5$fdr, 
              DEDD.results.compcodeR.DEDD50$fdr), 
        col=col_vector[1:4], pch=20, cex.axis=2, xaxt='n', 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(2.5,8.5), labels=c(5,50), tick=F, cex.axis=2.5, mgp=c(3,1,0))
mtext('Simulated data', cex=1.5, side=3, line=-2)
plot(NA, xlim=c(0,10), ylim=c(0,1), xaxt="n", yaxt="n", 
     main="Sensitivity", cex.main=2.5)
boxplot(cbind(DEDD.results.compcodeR.DEDD5$tpr, 
              DEDD.results.compcodeR.DEDD50$tpr), 
        col=col_vector[1:4], pch=20, cex.axis=2, xaxt='n', 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(2.5,8.5), labels=c(5,50), tick=F, cex.axis=2.5, mgp=c(3,1,0))
mtext('Simulated data', cex=1.5, side=3, line=-2)
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.5, ncol=1, 
       legend=c("diffVar", "Hybrid", "HMM, posterior threshold", "HMM, BFDR"))

# Second row - blood
plot(NA, xlim=c(0,10), ylim=c(0,1), xaxt="n", yaxt="n")
lines(c(0.5,4.5), c(0.05,0.05), col="lightgrey")
lines(c(5.5,9.5), c(0.05,0.05), col="lightgrey")
boxplot(cbind(fdr_blood_5, fdr_blood_50), 
        col=col_vector[1:4], pch=20, cex.axis=2, xaxt='n', 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(2.5,8.5), labels=c(5,50), tick=F, cex.axis=2.5, mgp=c(3,1,0))
mtext('GTEx blood data', cex=1.5, side=3, line=-2)
plot(NA, xlim=c(0,10), ylim=c(0,1), xaxt="n", yaxt="n")
boxplot(cbind(tpr_blood_5, tpr_blood_50), 
        col=col_vector[1:4], pch=20, cex.axis=2, xaxt='n', 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(2.5,8.5), labels=c(5,50), tick=F, cex.axis=2.5, mgp=c(3,1,0))
mtext('GTEx blood data', cex=1.5, side=3, line=-2)

# Third row - muscle
plot(NA, xlim=c(0,10), ylim=c(0,1), xaxt="n", yaxt="n")
lines(c(0.5,4.5), c(0.05,0.05), col="lightgrey")
lines(c(5.5,9.5), c(0.05,0.05), col="lightgrey")
boxplot(cbind(fdr_muscle_5,  fdr_muscle_50), 
        col=col_vector[1:4], pch=20, cex.axis=2, xaxt='n', 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(2.5,8.5), labels=c(5,50), tick=F, cex.axis=2.5, mgp=c(3,1,0))
mtext('GTEx muscle data', cex=1.5, side=3, line=-2)
plot(NA, xlim=c(0,10), ylim=c(0,1), xaxt="n", yaxt="n")
boxplot(cbind(tpr_muscle_5, tpr_muscle_50), 
        col=col_vector[1:4], pch=20, cex.axis=2, xaxt='n', 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(2.5,8.5), labels=c(5,50), tick=F, cex.axis=2.5, mgp=c(3,1,0))
mtext('GTEx muscle data', cex=1.5, side=3, line=-2)
