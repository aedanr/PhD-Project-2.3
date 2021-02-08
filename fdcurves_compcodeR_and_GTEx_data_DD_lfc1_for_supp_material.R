library(here)
library(RColorBrewer)

# Make colour vector
qual_col_pals = brewer.pal.info[brewer.pal.info$category == "qual" & 
                                  brewer.pal.info$colorblind == T,]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector <- col_vector[c(1:3,9,11,26)]
rm(qual_col_pals)

## Load compcodeR data and keep only MDSeq ZI and lnHM log ####
folder <- "Results/compcodeR DE, DD, DEDD results Feb 2020"
for (i in c("DEDD2", "DEDD5", "DEDD10", "DEDD20", "DEDD50")) {
  assign(paste0("DD.results.compcodeR.", i), 
         readRDS(here(folder, paste0("DD.lfc1.results.", i, ".rds"))))
}
rm(i,folder)

## Load GTEx data ####
folder <- "Results/GTEx DD lfc1 results Jan 2021"
for (i in c("2", "5", "10", "20", "50")) {
  for (j in c("blood", "muscle")) {
    assign(paste0("DD.results.", j, ".DEDD", i), 
           readRDS(here(folder, paste0("DD.lfc1.results.", j, "_", i, "_DEDD.rds"))))
  }
}
rm(i,j,folder)


par(mfrow=c(3,5), mar=c(2.5,3.5,3,1), mgp=c(3,1.5,0))

# First row - compcodeR
plot(DD.results.compcodeR.DEDD2$mean.discoveries$disp.lfc1.MDSeq.zi.tmm, 
     DD.results.compcodeR.DEDD2$mean.fdr$disp.lfc1.MDSeq.zi.tmm, 
     type='l', ylim=c(0,0.95), xlim=c(0,2000), lwd=2, col=col_vector[1], 
     main=paste0("2 samples per group"), cex.main=2.5, yaxt='n', xaxt='n', ylab=NA, xlab=NA)
axis(side=1, mgp=c(3,1.5,0), cex.axis=2.5)
axis(side=2, mgp=c(3,1,0), cex.axis=2.5)
lines(DD.results.compcodeR.DEDD2$mean.discoveries$disp.lfc1.lnHM.tmm, 
      DD.results.compcodeR.DEDD2$mean.fdr$disp.lfc1.lnHM.tmm, 
      lwd=2, col=col_vector[3])
legend("bottomright", col=col_vector[c(1,3)], bty='n', cex=2, ncol=1, lty=1, lwd=2, 
       legend=c("MDSeq", "HM"))

plot(DD.results.compcodeR.DEDD5$mean.discoveries$disp.lfc1.MDSeq.zi.tmm, 
     DD.results.compcodeR.DEDD5$mean.fdr$disp.lfc1.MDSeq.zi.tmm, 
     type='l', ylim=c(0,0.95), xlim=c(0,2000), lwd=2, col=col_vector[1], 
     main=paste0("5 samples per group"), cex.main=2.5, yaxt='n', xaxt='n', ylab=NA, xlab=NA)
axis(side=1, mgp=c(3,1.5,0), cex.axis=2.5)
axis(side=2, mgp=c(3,1,0), cex.axis=2.5)
lines(DD.results.compcodeR.DEDD5$mean.discoveries$disp.lfc1.lnHM.tmm, 
      DD.results.compcodeR.DEDD5$mean.fdr$disp.lfc1.lnHM.tmm, 
      lwd=2, col=col_vector[3])

plot(DD.results.compcodeR.DEDD10$mean.discoveries$disp.lfc1.MDSeq.zi.tmm, 
     DD.results.compcodeR.DEDD10$mean.fdr$disp.lfc1.MDSeq.zi.tmm, 
     type='l', ylim=c(0,0.95), xlim=c(0,2000), lwd=2, col=col_vector[1], 
     main=paste0("10 samples per group"), cex.main=2.5, yaxt='n', xaxt='n', ylab=NA, xlab=NA)
axis(side=1, mgp=c(3,1.5,0), cex.axis=2.5)
axis(side=2, mgp=c(3,1,0), cex.axis=2.5)
mtext('Simulated data', cex=1.5, side=3, line=-2.5)
lines(DD.results.compcodeR.DEDD10$mean.discoveries$disp.lfc1.lnHM.tmm, 
      DD.results.compcodeR.DEDD10$mean.fdr$disp.lfc1.lnHM.tmm, 
      lwd=2, col=col_vector[3])

plot(DD.results.compcodeR.DEDD20$mean.discoveries$disp.lfc1.MDSeq.zi.tmm, 
     DD.results.compcodeR.DEDD20$mean.fdr$disp.lfc1.MDSeq.zi.tmm, 
     type='l', ylim=c(0,0.95), xlim=c(0,2000), lwd=2, col=col_vector[1], 
     main=paste0("20 samples per group"), cex.main=2.5, yaxt='n', xaxt='n', ylab=NA, xlab=NA)
axis(side=1, mgp=c(3,1.5,0), cex.axis=2.5)
axis(side=2, mgp=c(3,1,0), cex.axis=2.5)
lines(DD.results.compcodeR.DEDD20$mean.discoveries$disp.lfc1.lnHM.tmm, 
      DD.results.compcodeR.DEDD20$mean.fdr$disp.lfc1.lnHM.tmm, 
      lwd=2, col=col_vector[3])

plot(DD.results.compcodeR.DEDD50$mean.discoveries$disp.lfc1.MDSeq.zi.tmm, 
     DD.results.compcodeR.DEDD50$mean.fdr$disp.lfc1.MDSeq.zi.tmm, 
     type='l', ylim=c(0,0.95), xlim=c(0,2000), lwd=2, col=col_vector[1], 
     main=paste0("50 samples per group"), cex.main=2.5, yaxt='n', xaxt='n', ylab=NA, xlab=NA)
axis(side=1, mgp=c(3,1.5,0), cex.axis=2.5)
axis(side=2, mgp=c(3,1,0), cex.axis=2.5)
lines(DD.results.compcodeR.DEDD50$mean.discoveries$disp.lfc1.lnHM.tmm, 
      DD.results.compcodeR.DEDD50$mean.fdr$disp.lfc1.lnHM.tmm, 
      lwd=2, col=col_vector[3])

# Second row - blood
plot(DD.results.blood.DEDD2$mean.discoveries$disp.lfc1.MDSeq.zi, 
     DD.results.blood.DEDD2$mean.fdr$disp.lfc1.MDSeq.zi, 
     type='l', ylim=c(0,0.95), xlim=c(0,2000), lwd=2, col=col_vector[1], 
     yaxt='n', xaxt='n', ylab=NA, xlab=NA)
axis(side=1, mgp=c(3,1.5,0), cex.axis=2.5)
axis(side=2, mgp=c(3,1,0), cex.axis=2.5)
lines(DD.results.blood.DEDD2$mean.discoveries$disp.lfc1.lnHM, 
      DD.results.blood.DEDD2$mean.fdr$disp.lfc1.lnHM, 
      lwd=2, col=col_vector[3])

plot(DD.results.blood.DEDD5$mean.discoveries$disp.lfc1.MDSeq.zi, 
     DD.results.blood.DEDD5$mean.fdr$disp.lfc1.MDSeq.zi, 
     type='l', ylim=c(0,0.95), xlim=c(0,2000), lwd=2, col=col_vector[1], 
     yaxt='n', xaxt='n', ylab=NA, xlab=NA)
axis(side=1, mgp=c(3,1.5,0), cex.axis=2.5)
axis(side=2, mgp=c(3,1,0), cex.axis=2.5)
lines(DD.results.blood.DEDD5$mean.discoveries$disp.lfc1.lnHM, 
      DD.results.blood.DEDD5$mean.fdr$disp.lfc1.lnHM, 
      lwd=2, col=col_vector[3])

plot(DD.results.blood.DEDD10$mean.discoveries$disp.lfc1.MDSeq.zi, 
     DD.results.blood.DEDD10$mean.fdr$disp.lfc1.MDSeq.zi, 
     type='l', ylim=c(0,0.95), xlim=c(0,2000), lwd=2, col=col_vector[1], 
     yaxt='n', xaxt='n', ylab=NA, xlab=NA)
axis(side=1, mgp=c(3,1.5,0), cex.axis=2.5)
axis(side=2, mgp=c(3,1,0), cex.axis=2.5)
mtext('GTEx blood data', cex=1.5, side=3, line=-2.5)
lines(DD.results.blood.DEDD10$mean.discoveries$disp.lfc1.lnHM, 
      DD.results.blood.DEDD10$mean.fdr$disp.lfc1.lnHM, 
      lwd=2, col=col_vector[3])

plot(DD.results.blood.DEDD20$mean.discoveries$disp.lfc1.MDSeq.zi, 
     DD.results.blood.DEDD20$mean.fdr$disp.lfc1.MDSeq.zi, 
     type='l', ylim=c(0,0.95), xlim=c(0,2000), lwd=2, col=col_vector[1], 
     yaxt='n', xaxt='n', ylab=NA, xlab=NA)
axis(side=1, mgp=c(3,1.5,0), cex.axis=2.5)
axis(side=2, mgp=c(3,1,0), cex.axis=2.5)
lines(DD.results.blood.DEDD20$mean.discoveries$disp.lfc1.lnHM, 
      DD.results.blood.DEDD20$mean.fdr$disp.lfc1.lnHM, 
      lwd=2, col=col_vector[3])

plot(DD.results.blood.DEDD50$mean.discoveries$disp.lfc1.MDSeq.zi, 
     DD.results.blood.DEDD50$mean.fdr$disp.lfc1.MDSeq.zi, 
     type='l', ylim=c(0,0.95), xlim=c(0,2000), lwd=2, col=col_vector[1], 
     yaxt='n', xaxt='n', ylab=NA, xlab=NA)
axis(side=1, mgp=c(3,1.5,0), cex.axis=2.5)
axis(side=2, mgp=c(3,1,0), cex.axis=2.5)
lines(DD.results.blood.DEDD50$mean.discoveries$disp.lfc1.lnHM, 
      DD.results.blood.DEDD50$mean.fdr$disp.lfc1.lnHM, 
      lwd=2, col=col_vector[3])

# Third row - muscle
plot(DD.results.muscle.DEDD2$mean.discoveries$disp.lfc1.MDSeq.zi, 
     DD.results.muscle.DEDD2$mean.fdr$disp.lfc1.MDSeq.zi, 
     type='l', ylim=c(0,0.95), xlim=c(0,2000), lwd=2, col=col_vector[1], 
     yaxt='n', xaxt='n', ylab=NA, xlab=NA)
axis(side=1, mgp=c(3,1.5,0), cex.axis=2.5)
axis(side=2, mgp=c(3,1,0), cex.axis=2.5)
lines(DD.results.muscle.DEDD2$mean.discoveries$disp.lfc1.lnHM, 
      DD.results.muscle.DEDD2$mean.fdr$disp.lfc1.lnHM, 
      lwd=2, col=col_vector[3])

plot(DD.results.muscle.DEDD5$mean.discoveries$disp.lfc1.MDSeq.zi, 
     DD.results.muscle.DEDD5$mean.fdr$disp.lfc1.MDSeq.zi, 
     type='l', ylim=c(0,0.95), xlim=c(0,2000), lwd=2, col=col_vector[1], 
     yaxt='n', xaxt='n', ylab=NA, xlab=NA)
axis(side=1, mgp=c(3,1.5,0), cex.axis=2.5)
axis(side=2, mgp=c(3,1,0), cex.axis=2.5)
lines(DD.results.muscle.DEDD5$mean.discoveries$disp.lfc1.lnHM, 
      DD.results.muscle.DEDD5$mean.fdr$disp.lfc1.lnHM, 
      lwd=2, col=col_vector[3])

plot(DD.results.muscle.DEDD10$mean.discoveries$disp.lfc1.MDSeq.zi, 
     DD.results.muscle.DEDD10$mean.fdr$disp.lfc1.MDSeq.zi, 
     type='l', ylim=c(0,0.95), xlim=c(0,2000), lwd=2, col=col_vector[1], 
     yaxt='n', xaxt='n', ylab=NA, xlab=NA)
axis(side=1, mgp=c(3,1.5,0), cex.axis=2.5)
axis(side=2, mgp=c(3,1,0), cex.axis=2.5)
mtext('GTEx muscle data', cex=1.5, side=3, line=-2.5)
lines(DD.results.muscle.DEDD10$mean.discoveries$disp.lfc1.lnHM, 
      DD.results.muscle.DEDD10$mean.fdr$disp.lfc1.lnHM, 
      lwd=2, col=col_vector[3])

plot(DD.results.muscle.DEDD20$mean.discoveries$disp.lfc1.MDSeq.zi, 
     DD.results.muscle.DEDD20$mean.fdr$disp.lfc1.MDSeq.zi, 
     type='l', ylim=c(0,0.95), xlim=c(0,2000), lwd=2, col=col_vector[1], 
     yaxt='n', xaxt='n', ylab=NA, xlab=NA)
axis(side=1, mgp=c(3,1.5,0), cex.axis=2.5)
axis(side=2, mgp=c(3,1,0), cex.axis=2.5)
lines(DD.results.muscle.DEDD20$mean.discoveries$disp.lfc1.lnHM, 
      DD.results.muscle.DEDD20$mean.fdr$disp.lfc1.lnHM, 
      lwd=2, col=col_vector[3])

plot(DD.results.muscle.DEDD50$mean.discoveries$disp.lfc1.MDSeq.zi, 
     DD.results.muscle.DEDD50$mean.fdr$disp.lfc1.MDSeq.zi, 
     type='l', ylim=c(0,0.95), xlim=c(0,2000), lwd=2, col=col_vector[1], 
     yaxt='n', xaxt='n', ylab=NA, xlab=NA)
axis(side=1, mgp=c(3,1.5,0), cex.axis=2.5)
axis(side=2, mgp=c(3,1,0), cex.axis=2.5)
lines(DD.results.muscle.DEDD50$mean.discoveries$disp.lfc1.lnHM, 
      DD.results.muscle.DEDD50$mean.fdr$disp.lfc1.lnHM, 
      lwd=2, col=col_vector[3])
