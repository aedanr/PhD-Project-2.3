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
folder <- "Results/GTEx combined results Dec 2020"
for (i in c("2", "5", "10", "20", "50")) {
  assign(paste0("DD.results.GTEx.DEDD", i), 
         readRDS(here(folder, paste0("DD.results.DEDD", i, ".rds"))))
}
rm(i,folder)

## Load compcodeR data ####
folder <- "Results/compcodeR combined results Dec 2020 including diffVar"
for (i in c("2", "5", "10", "20", "50")) {
  assign(paste0("DD.results.compcodeR.DEDD", i), 
         readRDS(here(folder, paste0("DD.results.DEDD", i, ".rds"))))
}
rm(i,folder)


## FDR curves - MDSeq ZI, ExpVarQuant, lnHM log ####
par(mfrow=c(3,2), mar=c(2.5,3.5,3,1), mgp=c(3,1.5,0))

# First row - compcodeR
plot(DD.results.compcodeR.DEDD5$mean.discoveries$disp.MDSeq.zi.tmm, 
     DD.results.compcodeR.DEDD5$mean.fdr$disp.MDSeq.zi.tmm, 
     type='l', ylim=c(0,0.95), xlim=c(0,2000), lwd=2, col=col_vector[1], 
     main=paste0("5 samples per group"), cex.main=2.5, yaxt='n', xaxt='n', ylab=NA, xlab=NA)
axis(side=1, mgp=c(3,1.5,0), cex.axis=2.5)
axis(side=2, mgp=c(3,1,0), cex.axis=2.5)
mtext('Simulated data', cex=1.5, side=3, line=-2.5)
lines(DD.results.compcodeR.DEDD5$mean.discoveries$EVQ.tmm, 
     DD.results.compcodeR.DEDD5$mean.fdr$EVQ.tmm, 
     lwd=2, col=col_vector[2])
lines(DD.results.compcodeR.DEDD5$mean.discoveries$disp.lnHM.log.tmm, 
      DD.results.compcodeR.DEDD5$mean.fdr$disp.lnHM.log.tmm, 
      lwd=2, col=col_vector[3])
legend("bottomright", col=col_vector[1:3], bty='n', cex=2.5, ncol=1, lty=1, lwd=2, 
       legend=c("MDSeq", "GAMLSS", "HM"))
plot(DD.results.compcodeR.DEDD50$mean.discoveries$disp.MDSeq.zi.tmm, 
     DD.results.compcodeR.DEDD50$mean.fdr$disp.MDSeq.zi.tmm, 
     type='l', ylim=c(0,0.6), xlim=c(0,2000), lwd=2, col=col_vector[1], 
     main=paste0("50 samples per group"), cex.main=2.5, yaxt='n', xaxt='n', ylab=NA, xlab=NA)
axis(side=1, mgp=c(3,1.5,0), cex.axis=2.5)
axis(side=2, mgp=c(3,1,0), cex.axis=2.5)
mtext('Simulated data', cex=1.5, side=3, line=-2.5)
lines(DD.results.compcodeR.DEDD50$mean.discoveries$EVQ.tmm, 
      DD.results.compcodeR.DEDD50$mean.fdr$EVQ.tmm, 
      lwd=2, col=col_vector[2])
lines(DD.results.compcodeR.DEDD50$mean.discoveries$disp.lnHM.log.tmm, 
      DD.results.compcodeR.DEDD50$mean.fdr$disp.lnHM.log.tmm, 
      lwd=2, col=col_vector[3])

# Second row - blood
plot(DD.results.GTEx.DEDD5$mean.discoveries$blood_MD.zi, 
     DD.results.GTEx.DEDD5$mean.fdr$blood_MD.zi, 
     type='l', ylim=c(0,0.95), xlim=c(0,2000), lwd=2, col=col_vector[1], 
     yaxt='n', xaxt='n', ylab=NA, xlab=NA)
axis(side=1, mgp=c(3,1.5,0), cex.axis=2.5)
axis(side=2, mgp=c(3,1,0), cex.axis=2.5)
mtext('GTEx blood data', cex=1.5, side=3, line=-2.5)
lines(DD.results.GTEx.DEDD5$mean.discoveries$blood_EVQ, 
      DD.results.GTEx.DEDD5$mean.fdr$blood_EVQ, 
      lwd=2, col=col_vector[2], cex.axis=2.5)
lines(DD.results.GTEx.DEDD5$mean.discoveries$blood_lnHM.log, 
      DD.results.GTEx.DEDD5$mean.fdr$blood_lnHM.log, 
      lwd=2, col=col_vector[3], cex.axis=2.5)
plot(DD.results.GTEx.DEDD50$mean.discoveries$blood_MD.zi, 
     DD.results.GTEx.DEDD50$mean.fdr$blood_MD.zi, 
     type='l', ylim=c(0,0.6), xlim=c(0,2000), lwd=2, col=col_vector[1], 
     yaxt='n', xaxt='n', ylab=NA, xlab=NA)
axis(side=1, mgp=c(3,1.5,0), cex.axis=2.5)
axis(side=2, mgp=c(3,1,0), cex.axis=2.5)
mtext('GTEx blood data', cex=1.5, side=3, line=-2.5)
lines(DD.results.GTEx.DEDD50$mean.discoveries$blood_EVQ, 
      DD.results.GTEx.DEDD50$mean.fdr$blood_EVQ, 
      lwd=2, col=col_vector[2], cex.axis=2.5)
lines(DD.results.GTEx.DEDD50$mean.discoveries$blood_lnHM.log, 
      DD.results.GTEx.DEDD50$mean.fdr$blood_lnHM.log, 
      lwd=2, col=col_vector[3], cex.axis=2.5)

# Third row - muscle
plot(DD.results.GTEx.DEDD5$mean.discoveries$muscle_MD.zi, 
     DD.results.GTEx.DEDD5$mean.fdr$muscle_MD.zi, 
     type='l', ylim=c(0,0.95), xlim=c(0,2000), lwd=2, col=col_vector[1], 
     yaxt='n', xaxt='n', ylab=NA, xlab=NA)
axis(side=1, mgp=c(3,1.5,0), cex.axis=2.5)
axis(side=2, mgp=c(3,1,0), cex.axis=2.5)
mtext('GTEx muscle data', cex=1.5, side=3, line=-2.5)
lines(DD.results.GTEx.DEDD5$mean.discoveries$muscle_EVQ, 
      DD.results.GTEx.DEDD5$mean.fdr$muscle_EVQ, 
      lwd=2, col=col_vector[2])
lines(DD.results.GTEx.DEDD5$mean.discoveries$muscle_lnHM.log, 
      DD.results.GTEx.DEDD5$mean.fdr$muscle_lnHM.log, 
      lwd=2, col=col_vector[3])
plot(DD.results.GTEx.DEDD50$mean.discoveries$muscle_MD.zi, 
     DD.results.GTEx.DEDD50$mean.fdr$muscle_MD.zi, 
     type='l', ylim=c(0,0.6), xlim=c(0,2000), lwd=2, col=col_vector[1], 
     yaxt='n', xaxt='n', ylab=NA, xlab=NA)
axis(side=1, mgp=c(3,1.5,0), cex.axis=2.5)
axis(side=2, mgp=c(3,1,0), cex.axis=2.5)
mtext('GTEx muscle data', cex=1.5, side=3, line=-2.5)
lines(DD.results.GTEx.DEDD50$mean.discoveries$muscle_EVQ, 
      DD.results.GTEx.DEDD50$mean.fdr$muscle_EVQ, 
      lwd=2, col=col_vector[2])
lines(DD.results.GTEx.DEDD50$mean.discoveries$muscle_lnHM.log, 
      DD.results.GTEx.DEDD50$mean.fdr$muscle_lnHM.log, 
      lwd=2, col=col_vector[3])
