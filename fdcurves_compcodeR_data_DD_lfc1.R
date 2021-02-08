library(here)
library(RColorBrewer)

## Colour vector ####
qual_col_pals = brewer.pal.info[brewer.pal.info$category == "qual" & 
                                  brewer.pal.info$colorblind == T,]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# pie(rep(1,length(col_vector)), col=col_vector)
col_vector <- col_vector[c(1:3,9,11,26)]
rm(qual_col_pals)

## Load compcodeR data ####
folder <- "Results/compcodeR DE, DD, DEDD results Feb 2020"
for (i in c("2", "5", "10", "20", "50")) {
  assign(paste0("DD.results.compcodeR.DEDD", i), 
         readRDS(here(folder, paste0("DD.results.DEDD", i, ".rds"))))
}
folder <- "Results/compcodeR DE, DD, DEDD results Feb 2020"
for (i in c("2", "5", "10", "20", "50")) {
  assign(paste0("DD.lfc1.results.compcodeR.DEDD", i), 
         readRDS(here(folder, paste0("DD.lfc1.results.DEDD", i, ".rds"))))
}
rm(i,folder)


## FDR curves - MDSeq ZI, lnHM log ####
par(mfcol=c(5,2), mar=c(2.5,3.5,3,1), mgp=c(3,1.5,0))

# Min lfc 0
plot(DD.results.compcodeR.DEDD2$mean.discoveries$disp.MDSeq.zi.tmm, 
     DD.results.compcodeR.DEDD2$mean.fdr$disp.MDSeq.zi.tmm, 
     type='l', ylim=c(0,1), xlim=c(0,2000), lwd=2, col=col_vector[1], 
     main=paste0("5 samples per group"), cex.main=2.5, yaxt='n', xaxt='n', ylab=NA, xlab=NA)
axis(side=1, mgp=c(3,1.5,0), cex.axis=2.5)
axis(side=2, mgp=c(3,1,0), cex.axis=2.5)
lines(DD.results.compcodeR.DEDD2$mean.discoveries$disp.lnHM.log.tmm, 
      DD.results.compcodeR.DEDD2$mean.fdr$disp.lnHM.log.tmm, 
      lwd=2, col=col_vector[2])
legend("bottomright", col=col_vector[1:3], bty='n', cex=2, ncol=1, lty=1, lwd=2, 
       legend=c("MDSeq", "HM"))

plot(DD.results.compcodeR.DEDD5$mean.discoveries$disp.MDSeq.zi.tmm, 
     DD.results.compcodeR.DEDD5$mean.fdr$disp.MDSeq.zi.tmm, 
     type='l', ylim=c(0,1), xlim=c(0,2000), lwd=2, col=col_vector[1], 
     main=paste0("5 samples per group"), cex.main=2.5, yaxt='n', xaxt='n', ylab=NA, xlab=NA)
axis(side=1, mgp=c(3,1.5,0), cex.axis=2.5)
axis(side=2, mgp=c(3,1,0), cex.axis=2.5)
lines(DD.results.compcodeR.DEDD5$mean.discoveries$disp.lnHM.log.tmm, 
      DD.results.compcodeR.DEDD5$mean.fdr$disp.lnHM.log.tmm, 
      lwd=2, col=col_vector[2])

plot(DD.results.compcodeR.DEDD10$mean.discoveries$disp.MDSeq.zi.tmm, 
     DD.results.compcodeR.DEDD10$mean.fdr$disp.MDSeq.zi.tmm, 
     type='l', ylim=c(0,1), xlim=c(0,2000), lwd=2, col=col_vector[1], 
     main=paste0("10 samples per group"), cex.main=2.5, yaxt='n', xaxt='n', ylab=NA, xlab=NA)
axis(side=1, mgp=c(3,1.5,0), cex.axis=2.5)
axis(side=2, mgp=c(3,1,0), cex.axis=2.5)
lines(DD.results.compcodeR.DEDD10$mean.discoveries$disp.lnHM.log.tmm, 
      DD.results.compcodeR.DEDD10$mean.fdr$disp.lnHM.log.tmm, 
      lwd=2, col=col_vector[2])

plot(DD.results.compcodeR.DEDD20$mean.discoveries$disp.MDSeq.zi.tmm, 
     DD.results.compcodeR.DEDD20$mean.fdr$disp.MDSeq.zi.tmm, 
     type='l', ylim=c(0,1), xlim=c(0,2000), lwd=2, col=col_vector[1], 
     main=paste0("20 samples per group"), cex.main=2.5, yaxt='n', xaxt='n', ylab=NA, xlab=NA)
axis(side=1, mgp=c(3,1.5,0), cex.axis=2.5)
axis(side=2, mgp=c(3,1,0), cex.axis=2.5)
lines(DD.results.compcodeR.DEDD20$mean.discoveries$disp.lnHM.log.tmm, 
      DD.results.compcodeR.DEDD20$mean.fdr$disp.lnHM.log.tmm, 
      lwd=2, col=col_vector[2])

plot(DD.results.compcodeR.DEDD50$mean.discoveries$disp.MDSeq.zi.tmm, 
     DD.results.compcodeR.DEDD50$mean.fdr$disp.MDSeq.zi.tmm, 
     type='l', ylim=c(0,1), xlim=c(0,2000), lwd=2, col=col_vector[1], 
     main=paste0("50 samples per group"), cex.main=2.5, yaxt='n', xaxt='n', ylab=NA, xlab=NA)
axis(side=1, mgp=c(3,1.5,0), cex.axis=2.5)
axis(side=2, mgp=c(3,1,0), cex.axis=2.5)
lines(DD.results.compcodeR.DEDD50$mean.discoveries$disp.lnHM.log.tmm, 
      DD.results.compcodeR.DEDD50$mean.fdr$disp.lnHM.log.tmm, 
      lwd=2, col=col_vector[2])

# Min lfc 1
plot(DD.lfc1.results.compcodeR.DEDD2$mean.discoveries$disp.lfc1.MDSeq.zi.tmm, 
     DD.lfc1.results.compcodeR.DEDD2$mean.fdr$disp.lfc1.MDSeq.zi.tmm, 
     type='l', ylim=c(0,1), xlim=c(0,2000), lwd=2, col=col_vector[1], 
     main=paste0("5 samples per group"), cex.main=2.5, yaxt='n', xaxt='n', ylab=NA, xlab=NA)
axis(side=1, mgp=c(3,1.5,0), cex.axis=2.5)
axis(side=2, mgp=c(3,1,0), cex.axis=2.5)
lines(DD.lfc1.results.compcodeR.DEDD2$mean.discoveries$disp.lfc1.lnHM.tmm, 
      DD.lfc1.results.compcodeR.DEDD2$mean.fdr$disp.lfc1.lnHM.tmm, 
      lwd=2, col=col_vector[2])
legend("bottomright", col=col_vector[1:3], bty='n', cex=2, ncol=1, lty=1, lwd=2, 
       legend=c("MDSeq", "HM"))

plot(DD.lfc1.results.compcodeR.DEDD5$mean.discoveries$disp.lfc1.MDSeq.zi.tmm, 
     DD.lfc1.results.compcodeR.DEDD5$mean.fdr$disp.lfc1.MDSeq.zi.tmm, 
     type='l', ylim=c(0,1), xlim=c(0,2000), lwd=2, col=col_vector[1], 
     main=paste0("5 samples per group"), cex.main=2.5, yaxt='n', xaxt='n', ylab=NA, xlab=NA)
axis(side=1, mgp=c(3,1.5,0), cex.axis=2.5)
axis(side=2, mgp=c(3,1,0), cex.axis=2.5)
lines(DD.lfc1.results.compcodeR.DEDD5$mean.discoveries$disp.lfc1.lnHM.tmm, 
      DD.lfc1.results.compcodeR.DEDD5$mean.fdr$disp.lfc1.lnHM.tmm, 
      lwd=2, col=col_vector[2])

plot(DD.lfc1.results.compcodeR.DEDD10$mean.discoveries$disp.lfc1.MDSeq.zi.tmm, 
     DD.lfc1.results.compcodeR.DEDD10$mean.fdr$disp.lfc1.MDSeq.zi.tmm, 
     type='l', ylim=c(0,1), xlim=c(0,2000), lwd=2, col=col_vector[1], 
     main=paste0("10 samples per group"), cex.main=2.5, yaxt='n', xaxt='n', ylab=NA, xlab=NA)
axis(side=1, mgp=c(3,1.5,0), cex.axis=2.5)
axis(side=2, mgp=c(3,1,0), cex.axis=2.5)
lines(DD.lfc1.results.compcodeR.DEDD10$mean.discoveries$disp.lfc1.lnHM.tmm, 
      DD.lfc1.results.compcodeR.DEDD10$mean.fdr$disp.lfc1.lnHM.tmm, 
      lwd=2, col=col_vector[2])

plot(DD.lfc1.results.compcodeR.DEDD20$mean.discoveries$disp.lfc1.MDSeq.zi.tmm, 
     DD.lfc1.results.compcodeR.DEDD20$mean.fdr$disp.lfc1.MDSeq.zi.tmm, 
     type='l', ylim=c(0,1), xlim=c(0,2000), lwd=2, col=col_vector[1], 
     main=paste0("20 samples per group"), cex.main=2.5, yaxt='n', xaxt='n', ylab=NA, xlab=NA)
axis(side=1, mgp=c(3,1.5,0), cex.axis=2.5)
axis(side=2, mgp=c(3,1,0), cex.axis=2.5)
lines(DD.lfc1.results.compcodeR.DEDD20$mean.discoveries$disp.lfc1.lnHM.tmm, 
      DD.lfc1.results.compcodeR.DEDD20$mean.fdr$disp.lfc1.lnHM.tmm, 
      lwd=2, col=col_vector[2])

plot(DD.lfc1.results.compcodeR.DEDD50$mean.discoveries$disp.lfc1.MDSeq.zi.tmm, 
     DD.lfc1.results.compcodeR.DEDD50$mean.fdr$disp.lfc1.MDSeq.zi.tmm, 
     type='l', ylim=c(0,1), xlim=c(0,2000), lwd=2, col=col_vector[1], 
     main=paste0("50 samples per group"), cex.main=2.5, yaxt='n', xaxt='n', ylab=NA, xlab=NA)
axis(side=1, mgp=c(3,1.5,0), cex.axis=2.5)
axis(side=2, mgp=c(3,1,0), cex.axis=2.5)
lines(DD.lfc1.results.compcodeR.DEDD50$mean.discoveries$disp.lfc1.lnHM.tmm, 
      DD.lfc1.results.compcodeR.DEDD50$mean.fdr$disp.lfc1.lnHM.tmm, 
      lwd=2, col=col_vector[2])

