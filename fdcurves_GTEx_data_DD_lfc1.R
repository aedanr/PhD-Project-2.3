library(here)
library(RColorBrewer)

## Load data, create colour vector ####
folder <- "Results/GTEx DD lfc1 results Jan 2021"
for (i in c("2", "5", "10", "20", "50")) {
  for (j in c("blood", "muscle")) {
    assign(paste0("DD.lfc1.results.", j, ".DEDD", i), 
           readRDS(here(folder, paste0("DD.lfc1.results.", j, "_", i, "_DEDD.rds"))))
  }
}
rm(i,j,folder)

qual_col_pals = brewer.pal.info[brewer.pal.info$category == "qual" & 
                                  brewer.pal.info$colorblind == T,]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector <- col_vector[c(1:3,9,11,26)]
rm(qual_col_pals)


par(mfcol=c(5,2), mar=c(2.5,3.5,3,1), mgp=c(3,1.5,0))

# Blood
plot(DD.lfc1.results.blood.DEDD2$mean.discoveries$disp.lfc1.MDSeq.zi, 
     DD.lfc1.results.blood.DEDD2$mean.fdr$disp.lfc1.MDSeq.zi, 
     type='l', ylim=c(0,1), xlim=c(0,2000), lwd=2, col=col_vector[1], 
     main=paste0("2 samples per group"), cex.main=2.5, yaxt='n', xaxt='n', ylab=NA, xlab=NA)
axis(side=1, mgp=c(3,1.5,0), cex.axis=2.5)
axis(side=2, mgp=c(3,1,0), cex.axis=2.5)
lines(DD.lfc1.results.blood.DEDD2$mean.discoveries$disp.lfc1.lnHM, 
      DD.lfc1.results.blood.DEDD2$mean.fdr$disp.lfc1.lnHM, 
      lwd=2, col=col_vector[2])
legend("bottomright", col=col_vector[1:3], bty='n', cex=2, ncol=1, lty=1, lwd=2, 
       legend=c("MDSeq", "HM"))

plot(DD.lfc1.results.blood.DEDD5$mean.discoveries$disp.lfc1.MDSeq.zi, 
     DD.lfc1.results.blood.DEDD5$mean.fdr$disp.lfc1.MDSeq.zi, 
     type='l', ylim=c(0,1), xlim=c(0,2000), lwd=2, col=col_vector[1], 
     main=paste0("5 samples per group"), cex.main=2.5, yaxt='n', xaxt='n', ylab=NA, xlab=NA)
axis(side=1, mgp=c(3,1.5,0), cex.axis=2.5)
axis(side=2, mgp=c(3,1,0), cex.axis=2.5)
lines(DD.lfc1.results.blood.DEDD5$mean.discoveries$disp.lfc1.lnHM, 
      DD.lfc1.results.blood.DEDD5$mean.fdr$disp.lfc1.lnHM, 
      lwd=2, col=col_vector[2])

plot(DD.lfc1.results.blood.DEDD10$mean.discoveries$disp.lfc1.MDSeq.zi, 
     DD.lfc1.results.blood.DEDD10$mean.fdr$disp.lfc1.MDSeq.zi, 
     type='l', ylim=c(0,1), xlim=c(0,2000), lwd=2, col=col_vector[1], 
     main=paste0("10 samples per group"), cex.main=2.5, yaxt='n', xaxt='n', ylab=NA, xlab=NA)
axis(side=1, mgp=c(3,1.5,0), cex.axis=2.5)
axis(side=2, mgp=c(3,1,0), cex.axis=2.5)
lines(DD.lfc1.results.blood.DEDD10$mean.discoveries$disp.lfc1.lnHM, 
      DD.lfc1.results.blood.DEDD10$mean.fdr$disp.lfc1.lnHM, 
      lwd=2, col=col_vector[2])

plot(DD.lfc1.results.blood.DEDD20$mean.discoveries$disp.lfc1.MDSeq.zi, 
     DD.lfc1.results.blood.DEDD20$mean.fdr$disp.lfc1.MDSeq.zi, 
     type='l', ylim=c(0,1), xlim=c(0,2000), lwd=2, col=col_vector[1], 
     main=paste0("20 samples per group"), cex.main=2.5, yaxt='n', xaxt='n', ylab=NA, xlab=NA)
axis(side=1, mgp=c(3,1.5,0), cex.axis=2.5)
axis(side=2, mgp=c(3,1,0), cex.axis=2.5)
lines(DD.lfc1.results.blood.DEDD20$mean.discoveries$disp.lfc1.lnHM, 
      DD.lfc1.results.blood.DEDD20$mean.fdr$disp.lfc1.lnHM, 
      lwd=2, col=col_vector[2])

plot(DD.lfc1.results.blood.DEDD50$mean.discoveries$disp.lfc1.MDSeq.zi, 
     DD.lfc1.results.blood.DEDD50$mean.fdr$disp.lfc1.MDSeq.zi, 
     type='l', ylim=c(0,1), xlim=c(0,2000), lwd=2, col=col_vector[1], 
     main=paste0("50 samples per group"), cex.main=2.5, yaxt='n', xaxt='n', ylab=NA, xlab=NA)
axis(side=1, mgp=c(3,1.5,0), cex.axis=2.5)
axis(side=2, mgp=c(3,1,0), cex.axis=2.5)
lines(DD.lfc1.results.blood.DEDD50$mean.discoveries$disp.lfc1.lnHM, 
      DD.lfc1.results.blood.DEDD50$mean.fdr$disp.lfc1.lnHM, 
      lwd=2, col=col_vector[2])

# Muscle
plot(DD.lfc1.results.muscle.DEDD2$mean.discoveries$disp.lfc1.MDSeq.zi, 
     DD.lfc1.results.muscle.DEDD2$mean.fdr$disp.lfc1.MDSeq.zi, 
     type='l', ylim=c(0,1), xlim=c(0,2000), lwd=2, col=col_vector[1], 
     main=paste0("5 samples per group"), cex.main=2.5, yaxt='n', xaxt='n', ylab=NA, xlab=NA)
axis(side=1, mgp=c(3,1.5,0), cex.axis=2.5)
axis(side=2, mgp=c(3,1,0), cex.axis=2.5)
lines(DD.lfc1.results.muscle.DEDD2$mean.discoveries$disp.lfc1.lnHM, 
      DD.lfc1.results.muscle.DEDD2$mean.fdr$disp.lfc1.lnHM, 
      lwd=2, col=col_vector[2])
legend("bottomright", col=col_vector[1:3], bty='n', cex=2, ncol=1, lty=1, lwd=2, 
       legend=c("MDSeq", "HM"))

plot(DD.lfc1.results.muscle.DEDD5$mean.discoveries$disp.lfc1.MDSeq.zi, 
     DD.lfc1.results.muscle.DEDD5$mean.fdr$disp.lfc1.MDSeq.zi, 
     type='l', ylim=c(0,1), xlim=c(0,2000), lwd=2, col=col_vector[1], 
     main=paste0("5 samples per group"), cex.main=2.5, yaxt='n', xaxt='n', ylab=NA, xlab=NA)
axis(side=1, mgp=c(3,1.5,0), cex.axis=2.5)
axis(side=2, mgp=c(3,1,0), cex.axis=2.5)
lines(DD.lfc1.results.muscle.DEDD5$mean.discoveries$disp.lfc1.lnHM, 
      DD.lfc1.results.muscle.DEDD5$mean.fdr$disp.lfc1.lnHM, 
      lwd=2, col=col_vector[2])

plot(DD.lfc1.results.muscle.DEDD10$mean.discoveries$disp.lfc1.MDSeq.zi, 
     DD.lfc1.results.muscle.DEDD10$mean.fdr$disp.lfc1.MDSeq.zi, 
     type='l', ylim=c(0,1), xlim=c(0,2000), lwd=2, col=col_vector[1], 
     main=paste0("10 samples per group"), cex.main=2.5, yaxt='n', xaxt='n', ylab=NA, xlab=NA)
axis(side=1, mgp=c(3,1.5,0), cex.axis=2.5)
axis(side=2, mgp=c(3,1,0), cex.axis=2.5)
lines(DD.lfc1.results.muscle.DEDD10$mean.discoveries$disp.lfc1.lnHM, 
      DD.lfc1.results.muscle.DEDD10$mean.fdr$disp.lfc1.lnHM, 
      lwd=2, col=col_vector[2])

plot(DD.lfc1.results.muscle.DEDD20$mean.discoveries$disp.lfc1.MDSeq.zi, 
     DD.lfc1.results.muscle.DEDD20$mean.fdr$disp.lfc1.MDSeq.zi, 
     type='l', ylim=c(0,1), xlim=c(0,2000), lwd=2, col=col_vector[1], 
     main=paste0("20 samples per group"), cex.main=2.5, yaxt='n', xaxt='n', ylab=NA, xlab=NA)
axis(side=1, mgp=c(3,1.5,0), cex.axis=2.5)
axis(side=2, mgp=c(3,1,0), cex.axis=2.5)
lines(DD.lfc1.results.muscle.DEDD20$mean.discoveries$disp.lfc1.lnHM, 
      DD.lfc1.results.muscle.DEDD20$mean.fdr$disp.lfc1.lnHM, 
      lwd=2, col=col_vector[2])

plot(DD.lfc1.results.muscle.DEDD50$mean.discoveries$disp.lfc1.MDSeq.zi, 
     DD.lfc1.results.muscle.DEDD50$mean.fdr$disp.lfc1.MDSeq.zi, 
     type='l', ylim=c(0,1), xlim=c(0,2000), lwd=2, col=col_vector[1], 
     main=paste0("50 samples per group"), cex.main=2.5, yaxt='n', xaxt='n', ylab=NA, xlab=NA)
axis(side=1, mgp=c(3,1.5,0), cex.axis=2.5)
axis(side=2, mgp=c(3,1,0), cex.axis=2.5)
lines(DD.lfc1.results.muscle.DEDD50$mean.discoveries$disp.lfc1.lnHM, 
      DD.lfc1.results.muscle.DEDD50$mean.fdr$disp.lfc1.lnHM, 
      lwd=2, col=col_vector[2])

