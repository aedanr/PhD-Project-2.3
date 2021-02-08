library(here)
library(RColorBrewer)

## Colour vector ####
qual_col_pals = brewer.pal.info[brewer.pal.info$category == "qual" & 
                                  brewer.pal.info$colorblind == T,]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# pie(rep(1,length(col_vector)), col=col_vector)
col_vector <- col_vector[c(1:3,9,11,26)]
rm(qual_col_pals)

## Load GTEx data ####
folder <- "Results/GTEx combined results Sept 2020"
for (i in c("2", "5", "10", "20", "50")) {
  assign(paste0("DEDD.results.GTEx.DEDD", i), 
         readRDS(here(folder, paste0("DEDD.results.DEDD", i, ".rds"))))
}
rm(i,folder)

## Load compcodeR data ####
folder <- "Results/compcodeR combined results Dec 2020 including diffVar"
for (i in c("2", "5", "10", "20", "50")) {
  assign(paste0("DEDD.results.compcodeR.DEDD", i), 
         readRDS(here(folder, paste0("DEDD.results.DEDD", i, ".rds"))))
}
rm(i,folder)


## FDR curves - diffVar, hybrid, lnHMM ####
## Hybrid: edgeR/lnHM log for compcodeR, voom/lnHM log for GTEx
par(mfrow=c(3,2), mar=c(2.5,3.5,3,1), mgp=c(3,1.5,0))

# First row - compcodeR
plot(DEDD.results.compcodeR.DEDD5$mean.discoveries$dV.tmm, 
     DEDD.results.compcodeR.DEDD5$mean.fdr$dV.tmm, 
     type='l', ylim=c(0,1.05), xlim=c(0,2000), lwd=2, col=col_vector[1], 
     main=paste0("5 samples per group"), cex.main=2.5, yaxt='n', xaxt='n', ylab=NA, xlab=NA)
axis(side=1, mgp=c(3,1.5,0), cex.axis=2.5)
axis(side=2, mgp=c(3,1,0), cex.axis=2.5)
mtext('Simulated data', cex=1.5, side=3, line=-2.5)
lines(DEDD.results.compcodeR.DEDD5$mean.discoveries$edgeR_HM.tmm, 
      DEDD.results.compcodeR.DEDD5$mean.fdr$edgeR_HM.tmm, 
      lwd=2, col=col_vector[2])
lines(DEDD.results.compcodeR.DEDD5$mean.discoveries$lnHMM.tmm, 
      DEDD.results.compcodeR.DEDD5$mean.fdr$lnHMM.tmm, 
      lwd=2, col=col_vector[3])
plot(DEDD.results.compcodeR.DEDD50$mean.discoveries$dV.tmm, 
     DEDD.results.compcodeR.DEDD50$mean.fdr$dV.tmm, 
     type='l', ylim=c(0,0.25), xlim=c(0,2000), lwd=2, col=col_vector[1], 
     main=paste0("50 samples per group"), cex.main=2.5, yaxt='n', xaxt='n', ylab=NA, xlab=NA)
axis(side=1, mgp=c(3,1.5,0), cex.axis=2.5)
axis(side=2, mgp=c(3,1,0), cex.axis=2.5)
mtext('Simulated data', cex=1.5, side=3, line=-2.5)
lines(DEDD.results.compcodeR.DEDD50$mean.discoveries$edgeR_HM.tmm, 
      DEDD.results.compcodeR.DEDD50$mean.fdr$edgeR_HM.tmm, 
      lwd=2, col=col_vector[2])
lines(DEDD.results.compcodeR.DEDD50$mean.discoveries$lnHMM.tmm, 
      DEDD.results.compcodeR.DEDD50$mean.fdr$lnHMM.tmm, 
      lwd=2, col=col_vector[3])
legend("left", col=col_vector[1:3], bty='n', cex=2.5, ncol=1, lty=1, lwd=2, 
       legend=c("diffVar", "Hybrid", "HMM"))

# Second row - blood
plot(DEDD.results.GTEx.DEDD5$mean.discoveries$blood_dV, 
     DEDD.results.GTEx.DEDD5$mean.fdr$blood_dV, 
     type='l', ylim=c(0,1.05), xlim=c(0,2000), lwd=2, col=col_vector[1], 
     yaxt='n', xaxt='n', ylab=NA, xlab=NA)
axis(side=1, mgp=c(3,1.5,0), cex.axis=2.5)
axis(side=2, mgp=c(3,1,0), cex.axis=2.5)
mtext('GTEx blood data', cex=1.5, side=3, line=-2.5)
lines(DEDD.results.GTEx.DEDD5$mean.discoveries$blood_voom_HM, 
      DEDD.results.GTEx.DEDD5$mean.fdr$blood_voom_HM, 
      lwd=2, col=col_vector[2], cex.axis=2.5)
lines(DEDD.results.GTEx.DEDD5$mean.discoveries$blood_lnHMM, 
      DEDD.results.GTEx.DEDD5$mean.fdr$blood_lnHMM, 
      lwd=2, col=col_vector[3], cex.axis=2.5)
plot(DEDD.results.GTEx.DEDD50$mean.discoveries$blood_dV, 
     DEDD.results.GTEx.DEDD50$mean.fdr$blood_dV, 
     type='l', ylim=c(0,0.25), xlim=c(0,2000), lwd=2, col=col_vector[1], 
     yaxt='n', xaxt='n', ylab=NA, xlab=NA)
axis(side=1, mgp=c(3,1.5,0), cex.axis=2.5)
axis(side=2, mgp=c(3,1,0), cex.axis=2.5)
mtext('GTEx blood data', cex=1.5, side=3, line=-2.5)
lines(DEDD.results.GTEx.DEDD50$mean.discoveries$blood_voom_HM, 
      DEDD.results.GTEx.DEDD50$mean.fdr$blood_voom_HM, 
      lwd=2, col=col_vector[2], cex.axis=2.5)
lines(DEDD.results.GTEx.DEDD50$mean.discoveries$blood_lnHMM, 
      DEDD.results.GTEx.DEDD50$mean.fdr$blood_lnHMM, 
      lwd=2, col=col_vector[3], cex.axis=2.5)

# Third row - muscle
plot(DEDD.results.GTEx.DEDD5$mean.discoveries$muscle_dV, 
     DEDD.results.GTEx.DEDD5$mean.fdr$muscle_dV, 
     type='l', ylim=c(0,1.05), xlim=c(0,2000), lwd=2, col=col_vector[1], 
     yaxt='n', xaxt='n', ylab=NA, xlab=NA)
axis(side=1, mgp=c(3,1.5,0), cex.axis=2.5)
axis(side=2, mgp=c(3,1,0), cex.axis=2.5)
mtext('GTEx muscle data', cex=1.5, side=3, line=-2.5)
lines(DEDD.results.GTEx.DEDD5$mean.discoveries$muscle_voom_HM, 
      DEDD.results.GTEx.DEDD5$mean.fdr$muscle_voom_HM, 
      lwd=2, col=col_vector[2], cex.axis=2.5)
lines(DEDD.results.GTEx.DEDD5$mean.discoveries$muscle_lnHMM, 
      DEDD.results.GTEx.DEDD5$mean.fdr$muscle_lnHMM, 
      lwd=2, col=col_vector[3], cex.axis=2.5)
plot(DEDD.results.GTEx.DEDD50$mean.discoveries$muscle_dV, 
     DEDD.results.GTEx.DEDD50$mean.fdr$muscle_dV, 
     type='l', ylim=c(0,0.25), xlim=c(0,2000), lwd=2, col=col_vector[1], 
     yaxt='n', xaxt='n', ylab=NA, xlab=NA)
axis(side=1, mgp=c(3,1.5,0), cex.axis=2.5)
axis(side=2, mgp=c(3,1,0), cex.axis=2.5)
mtext('GTEx muscle data', cex=1.5, side=3, line=-2.5)
lines(DEDD.results.GTEx.DEDD50$mean.discoveries$muscle_voom_HM, 
      DEDD.results.GTEx.DEDD50$mean.fdr$muscle_voom_HM, 
      lwd=2, col=col_vector[2], cex.axis=2.5)
lines(DEDD.results.GTEx.DEDD50$mean.discoveries$muscle_lnHMM, 
      DEDD.results.GTEx.DEDD50$mean.fdr$muscle_lnHMM, 
      lwd=2, col=col_vector[3], cex.axis=2.5)
