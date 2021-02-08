library(here)
library(RColorBrewer)

# Make colour vector
qual_col_pals = brewer.pal.info[brewer.pal.info$category == "qual" & 
                                  brewer.pal.info$colorblind == T,]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector <- col_vector[c(1:3,9,11,26)]
rm(qual_col_pals)

## Load compcodeR data and keep only MDSeq ZI, ExpVarQuant and lnHM log ####
folder <- "Results/compcodeR combined results Dec 2020 including diffVar"
for (i in c("DEDD5", "DEDD50")) {
  assign(paste0("DD.results.compcodeR.", i), 
         readRDS(here(folder, paste0("DD.results.", i, ".rds"))))
}
rm(i,folder)

{
  DD.results.compcodeR.DEDD5$fdr <- DD.results.compcodeR.DEDD5$fdr[, c("disp.MDSeq.zi.tmm", "EVQ.tmm", "disp.lnHM.log.tmm")]
  DD.results.compcodeR.DEDD5$tpr <- DD.results.compcodeR.DEDD5$tpr[, c("disp.MDSeq.zi.tmm", "EVQ.tmm", "disp.lnHM.log.tmm")]
  DD.results.compcodeR.DEDD50$fdr <- DD.results.compcodeR.DEDD50$fdr[, c("disp.MDSeq.zi.tmm", "EVQ.tmm", "disp.lnHM.log.tmm")]
  DD.results.compcodeR.DEDD50$tpr <- DD.results.compcodeR.DEDD50$tpr[, c("disp.MDSeq.zi.tmm", "EVQ.tmm", "disp.lnHM.log.tmm")]
}

## Load GTEx data, rearrange and keep only MDSeq ZI, ExpVarQuant and lnHM log ####
folder <- "Results/GTEx combined results Dec 2020"
for (i in c("5", "50")) {
  assign(paste0("DD.results.GTEx.DEDD", i), 
         readRDS(here(folder, paste0("DD.results.DEDD", i, ".rds"))))
}
rm(i,folder)

for (metric in c("fdr", "tpr")) {
  for (size in c("5", "50")) {
    assign(
      paste0(metric, "_blood_", size), 
      get(paste0("DD.results.GTEx.DEDD", size))[[metric]][
        , c("blood_MD.zi", "blood_EVQ", "blood_lnHM.log")
      ]
    )
    assign(
      paste0(metric, "_muscle_", size), 
      get(paste0("DD.results.GTEx.DEDD", size))[[metric]][
        , c("muscle_MD.zi", "muscle_EVQ", "muscle_lnHM.log")
      ]
    )
  }
}
rm(size, metric)

## Plot FDR and TPR for each dataset ####
positions = rep(1:6)
offsets <- c(rep(0, 3), rep(0.5, 3))

par(mfrow=c(3,2), mar=c(2,2.5,2.5,0.5), mgp=c(3,0.7,0))

# First row - compcodeR
plot(NA, xlim=c(0.5,7), ylim=c(0,1.1), xaxt="n", yaxt="n", 
     main="False discovery rate", cex.main=2.5)
lines(c(0.5,3.5), c(0.05,0.05), col="lightgrey")
lines(c(4,7), c(0.05,0.05), col="lightgrey")
boxplot(cbind(DD.results.compcodeR.DEDD5$fdr, 
              DD.results.compcodeR.DEDD50$fdr), 
        col=col_vector[1:3], pch=20, cex.axis=2, xaxt='n', 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(2,6), labels=c(5,50), tick=F, cex.axis=2.5, mgp=c(3,1,0))
mtext('Simulated data', cex=1.5, side=3, line=-2)
legend("topright", fill=col_vector[1:3], bty='n', cex=2, ncol=1, 
       legend=c("MDSeq", "GAMLSS", "HM"))
plot(NA, xlim=c(0.5,7), ylim=c(0,0.65), xaxt="n", yaxt="n", 
     main="Sensitivity", cex.main=2.5)
boxplot(cbind(DD.results.compcodeR.DEDD5$tpr, 
              DD.results.compcodeR.DEDD50$tpr), 
        col=col_vector[1:3], pch=20, cex.axis=2, xaxt='n', 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(2,6), labels=c(5,50), tick=F, cex.axis=2.5, mgp=c(3,1,0))
mtext('Simulated data', cex=1.5, side=3, line=-2)

# Second row - blood
plot(NA, xlim=c(0.5,7), ylim=c(0,1.1), xaxt="n", yaxt="n")
lines(c(0.5,3.5), c(0.05,0.05), col="lightgrey")
lines(c(4,7), c(0.05,0.05), col="lightgrey")
boxplot(cbind(fdr_blood_5, fdr_blood_50), 
        col=col_vector[1:3], pch=20, cex.axis=2, xaxt='n', 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(2,6), labels=c(5,50), tick=F, cex.axis=2.5, mgp=c(3,1,0))
mtext('GTEx blood data', cex=1.5, side=3, line=-2)
plot(NA, xlim=c(0.5,7), ylim=c(0,0.65), xaxt="n", yaxt="n")
boxplot(cbind(tpr_blood_5, tpr_blood_50), 
        col=col_vector[1:3], pch=20, cex.axis=2, xaxt='n', 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(2,6), labels=c(5,50), tick=F, cex.axis=2.5, mgp=c(3,1,0))
mtext('GTEx blood data', cex=1.5, side=3, line=-2)

# Third row - muscle
plot(NA, xlim=c(0.5,7), ylim=c(0,1.1), xaxt="n", yaxt="n")
lines(c(0.5,3.5), c(0.05,0.05), col="lightgrey")
lines(c(4,7), c(0.05,0.05), col="lightgrey")
boxplot(cbind(fdr_muscle_5, fdr_muscle_50), 
        col=col_vector[1:3], pch=20, cex.axis=2, xaxt='n', 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(2,6), labels=c(5,50), tick=F, cex.axis=2.5, mgp=c(3,1,0))
mtext('GTEx muscle data', cex=1.5, side=3, line=-2)
plot(NA, xlim=c(0.5,7), ylim=c(0,0.65), xaxt="n", yaxt="n")
boxplot(cbind(tpr_muscle_5, tpr_muscle_50), 
        col=col_vector[1:3], pch=20, cex.axis=2, xaxt='n', 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(2,6), labels=c(5,50), tick=F, cex.axis=2.5, mgp=c(3,1,0))
mtext('GTEx muscle data', cex=1.5, side=3, line=-2)
