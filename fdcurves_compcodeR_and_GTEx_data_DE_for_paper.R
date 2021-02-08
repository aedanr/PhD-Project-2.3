library(here)
library(RColorBrewer)

## Colour vector ####
qual_col_pals = brewer.pal.info[brewer.pal.info$category == "qual" & 
                                  brewer.pal.info$colorblind == T,]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector <- col_vector[c(1:3,9,11,26)]
rm(qual_col_pals)

## Load compcodeR data ####
folder <- "Results/compcodeR DE, DD, DEDD results Feb 2020"
for (i in c("2", "5", "10", "20", "50")) {
  assign(paste0("DE.results.compcodeR.DEDD", i), 
         readRDS(here(folder, paste0("DE.results.DEDD", i, ".rds"))))
}
rm(i,folder)

## Load GTEx data ####
folder <- "Results/GTEx combined results Sept 2020"
for (i in c("2", "5", "10", "20", "50")) {
  assign(paste0("DE.results.GTEx.DEDD", i), 
         readRDS(here(folder, paste0("DE.results.DEDD", i, ".rds"))))
}
rm(i,folder)


## FDR curves - edgeR QL, DESeq2 IF, voom, baySeq, MDSeq ZI, lnHM log ####
par(mfrow=c(3,2), mar=c(2.5,3.5,3,1), mgp=c(3,1.5,0))

# First row - compcodeR
plot(DE.results.compcodeR.DEDD5$mean.discoveries$edgeR.ql.tmm, 
     DE.results.compcodeR.DEDD5$mean.fdr$edgeR.ql.tmm, 
     type='l', ylim=c(0,0.5), xlim=c(0,2000), lwd=2, col=col_vector[1], 
     main=paste0("5 samples per group"), cex.main=2.5, yaxt='n', xaxt='n', ylab=NA, xlab=NA)
axis(side=1, mgp=c(3,1.5,0), cex.axis=2.5)
axis(side=2, mgp=c(3,1,0), cex.axis=2.5)
mtext('Simulated data', cex=1.5, side=3, line=-2.5)
lines(DE.results.compcodeR.DEDD5$mean.discoveries$DESeq2.if.tmm, 
     DE.results.compcodeR.DEDD5$mean.fdr$DESeq2.if.tmm, 
     lwd=2, col=col_vector[2])
lines(DE.results.compcodeR.DEDD5$mean.discoveries$voom.tmm, 
      DE.results.compcodeR.DEDD5$mean.fdr$voom.tmm, 
      lwd=2, col=col_vector[3])
lines(DE.results.compcodeR.DEDD5$mean.discoveries$baySeq.tmm, 
      DE.results.compcodeR.DEDD5$mean.fdr$baySeq.tmm, 
      lwd=2, col=col_vector[4])
lines(DE.results.compcodeR.DEDD5$mean.discoveries$mean.MDSeq.zi.tmm, 
      DE.results.compcodeR.DEDD5$mean.fdr$mean.MDSeq.zi.tmm, 
      lwd=2, col=col_vector[5])
lines(DE.results.compcodeR.DEDD5$mean.discoveries$mean.lnHM.log.tmm, 
      DE.results.compcodeR.DEDD5$mean.fdr$mean.lnHM.log.tmm, 
      lwd=2, col=col_vector[6])
plot(DE.results.compcodeR.DEDD50$mean.discoveries$edgeR.ql.tmm, 
     DE.results.compcodeR.DEDD50$mean.fdr$edgeR.ql.tmm, 
     type='l', ylim=c(0,0.25), xlim=c(0,2000), lwd=2, col=col_vector[1], 
     main=paste0("50 samples per group"), cex.main=2.5, yaxt='n', xaxt='n', ylab=NA, xlab=NA)
axis(side=1, mgp=c(3,1.5,0), cex.axis=2.5)
axis(side=2, mgp=c(3,1,0), cex.axis=2.5)
mtext('Simulated data', cex=1.5, side=3, line=-2.5)
lines(DE.results.compcodeR.DEDD50$mean.discoveries$DESeq2.if.tmm, 
      DE.results.compcodeR.DEDD50$mean.fdr$DESeq2.if.tmm, 
      lwd=2, col=col_vector[2])
lines(DE.results.compcodeR.DEDD50$mean.discoveries$voom.tmm, 
      DE.results.compcodeR.DEDD50$mean.fdr$voom.tmm, 
      lwd=2, col=col_vector[3])
lines(DE.results.compcodeR.DEDD50$mean.discoveries$baySeq.tmm, 
      DE.results.compcodeR.DEDD50$mean.fdr$baySeq.tmm, 
      lwd=2, col=col_vector[4])
lines(DE.results.compcodeR.DEDD50$mean.discoveries$mean.MDSeq.zi.tmm, 
      DE.results.compcodeR.DEDD50$mean.fdr$mean.MDSeq.zi.tmm, 
      lwd=2, col=col_vector[5])
lines(DE.results.compcodeR.DEDD50$mean.discoveries$mean.lnHM.log.tmm, 
      DE.results.compcodeR.DEDD50$mean.fdr$mean.lnHM.log.tmm, 
      lwd=2, col=col_vector[6])
legend("topleft", col=col_vector[1:6], bty='n', cex=2.5, ncol=1, lty=1, lwd=2, 
       legend=c("edgeR", "DESeq2", "voom", "baySeq", "MDSeq", "HM"))

# Second row - blood
plot(DE.results.GTEx.DEDD5$mean.discoveries$blood_eR.ql, 
     DE.results.GTEx.DEDD5$mean.fdr$blood_eR.ql, 
     type='l', ylim=c(0,0.5), xlim=c(0,2000), lwd=2, col=col_vector[1], 
     yaxt='n', xaxt='n', ylab=NA, xlab=NA)
axis(side=1, mgp=c(3,1.5,0), cex.axis=2.5)
axis(side=2, mgp=c(3,1,0), cex.axis=2.5)
mtext('GTEx blood data', cex=1.5, side=3, line=-2.5)
lines(DE.results.GTEx.DEDD5$mean.discoveries$blood_DES.if, 
      DE.results.GTEx.DEDD5$mean.fdr$blood_DES.if, 
      lwd=2, col=col_vector[2], cex.axis=2.5)
lines(DE.results.GTEx.DEDD5$mean.discoveries$blood_voom, 
      DE.results.GTEx.DEDD5$mean.fdr$blood_voom, 
      lwd=2, col=col_vector[3], cex.axis=2.5)
lines(DE.results.GTEx.DEDD5$mean.discoveries$blood_baySeq, 
      DE.results.GTEx.DEDD5$mean.fdr$blood_baySeq, 
      lwd=2, col=col_vector[4], cex.axis=2.5)
lines(DE.results.GTEx.DEDD5$mean.discoveries$blood_lnHM.log, 
      DE.results.GTEx.DEDD5$mean.fdr$blood_lnHM.log, 
      lwd=2, col=col_vector[5], cex.axis=2.5)
lines(DE.results.GTEx.DEDD5$mean.discoveries$blood_MD.zi, 
      DE.results.GTEx.DEDD5$mean.fdr$blood_MD.zi, 
      lwd=2, col=col_vector[6], cex.axis=2.5)
plot(DE.results.GTEx.DEDD50$mean.discoveries$blood_eR.ql, 
     DE.results.GTEx.DEDD50$mean.fdr$blood_eR.ql, 
     type='l', ylim=c(0,0.25), xlim=c(0,2000), lwd=2, col=col_vector[1], 
     yaxt='n', xaxt='n', ylab=NA, xlab=NA)
axis(side=1, mgp=c(3,1.5,0), cex.axis=2.5)
axis(side=2, mgp=c(3,1,0), cex.axis=2.5)
mtext('GTEx blood data', cex=1.5, side=3, line=-2.5)
lines(DE.results.GTEx.DEDD50$mean.discoveries$blood_DES.if, 
      DE.results.GTEx.DEDD50$mean.fdr$blood_DES.if, 
      lwd=2, col=col_vector[2], cex.axis=2.5)
lines(DE.results.GTEx.DEDD50$mean.discoveries$blood_voom, 
      DE.results.GTEx.DEDD50$mean.fdr$blood_voom, 
      lwd=2, col=col_vector[3], cex.axis=2.5)
lines(DE.results.GTEx.DEDD50$mean.discoveries$blood_baySeq, 
      DE.results.GTEx.DEDD50$mean.fdr$blood_baySeq, 
      lwd=2, col=col_vector[4], cex.axis=2.5)
lines(DE.results.GTEx.DEDD50$mean.discoveries$blood_MD.zi, 
      DE.results.GTEx.DEDD50$mean.fdr$blood_MD.zi, 
      lwd=2, col=col_vector[5], cex.axis=2.5)
lines(DE.results.GTEx.DEDD50$mean.discoveries$blood_lnHM.log, 
      DE.results.GTEx.DEDD50$mean.fdr$blood_lnHM.log, 
      lwd=2, col=col_vector[6], cex.axis=2.5)

# Third row - muscle
plot(DE.results.GTEx.DEDD5$mean.discoveries$muscle_eR.ql, 
     DE.results.GTEx.DEDD5$mean.fdr$muscle_eR.ql, 
     type='l', ylim=c(0,0.5), xlim=c(0,2000), lwd=2, col=col_vector[1], 
     yaxt='n', xaxt='n', ylab=NA, xlab=NA)
axis(side=1, mgp=c(3,1.5,0), cex.axis=2.5)
axis(side=2, mgp=c(3,1,0), cex.axis=2.5)
mtext('GTEx muscle data', cex=1.5, side=3, line=-2.5)
lines(DE.results.GTEx.DEDD5$mean.discoveries$muscle_DES.if, 
      DE.results.GTEx.DEDD5$mean.fdr$muscle_DES.if, 
      lwd=2, col=col_vector[2], cex.axis=2.5)
lines(DE.results.GTEx.DEDD5$mean.discoveries$muscle_voom, 
      DE.results.GTEx.DEDD5$mean.fdr$muscle_voom, 
      lwd=2, col=col_vector[3], cex.axis=2.5)
lines(DE.results.GTEx.DEDD5$mean.discoveries$muscle_baySeq, 
      DE.results.GTEx.DEDD5$mean.fdr$muscle_baySeq, 
      lwd=2, col=col_vector[4], cex.axis=2.5)
lines(DE.results.GTEx.DEDD5$mean.discoveries$muscle_lnHM.log, 
      DE.results.GTEx.DEDD5$mean.fdr$muscle_lnHM.log, 
      lwd=2, col=col_vector[5], cex.axis=2.5)
lines(DE.results.GTEx.DEDD5$mean.discoveries$muscle_MD.zi, 
      DE.results.GTEx.DEDD5$mean.fdr$muscle_MD.zi, 
      lwd=2, col=col_vector[6], cex.axis=2.5)
plot(DE.results.GTEx.DEDD50$mean.discoveries$muscle_eR.ql, 
     DE.results.GTEx.DEDD50$mean.fdr$muscle_eR.ql, 
     type='l', ylim=c(0,0.25), xlim=c(0,2000), lwd=2, col=col_vector[1], 
     yaxt='n', xaxt='n', ylab=NA, xlab=NA)
axis(side=1, mgp=c(3,1.5,0), cex.axis=2.5)
axis(side=2, mgp=c(3,1,0), cex.axis=2.5)
mtext('GTEx muscle data', cex=1.5, side=3, line=-2.5)
lines(DE.results.GTEx.DEDD50$mean.discoveries$muscle_DES.if, 
      DE.results.GTEx.DEDD50$mean.fdr$muscle_DES.if, 
      lwd=2, col=col_vector[2], cex.axis=2.5)
lines(DE.results.GTEx.DEDD50$mean.discoveries$muscle_voom, 
      DE.results.GTEx.DEDD50$mean.fdr$muscle_voom, 
      lwd=2, col=col_vector[3], cex.axis=2.5)
lines(DE.results.GTEx.DEDD50$mean.discoveries$muscle_baySeq, 
      DE.results.GTEx.DEDD50$mean.fdr$muscle_baySeq, 
      lwd=2, col=col_vector[4], cex.axis=2.5)
lines(DE.results.GTEx.DEDD50$mean.discoveries$muscle_MD.zi, 
      DE.results.GTEx.DEDD50$mean.fdr$muscle_MD.zi, 
      lwd=2, col=col_vector[5], cex.axis=2.5)
lines(DE.results.GTEx.DEDD50$mean.discoveries$muscle_lnHM.log, 
      DE.results.GTEx.DEDD50$mean.fdr$muscle_lnHM.log, 
      lwd=2, col=col_vector[6], cex.axis=2.5)
