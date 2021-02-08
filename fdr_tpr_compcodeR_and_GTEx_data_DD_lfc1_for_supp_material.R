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

{
  DD.results.compcodeR.DEDD2$fdr <- DD.results.compcodeR.DEDD2$fdr[, c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  DD.results.compcodeR.DEDD2$tpr <- DD.results.compcodeR.DEDD2$tpr[, c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  DD.results.compcodeR.DEDD5$fdr <- DD.results.compcodeR.DEDD5$fdr[, c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  DD.results.compcodeR.DEDD5$tpr <- DD.results.compcodeR.DEDD5$tpr[, c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  DD.results.compcodeR.DEDD10$fdr <- DD.results.compcodeR.DEDD10$fdr[, c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  DD.results.compcodeR.DEDD10$tpr <- DD.results.compcodeR.DEDD10$tpr[, c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  DD.results.compcodeR.DEDD20$fdr <- DD.results.compcodeR.DEDD20$fdr[, c("disp.lfc1.MDSeq.zi.tmm",  "disp.lfc1.lnHM.tmm")]
  DD.results.compcodeR.DEDD20$tpr <- DD.results.compcodeR.DEDD20$tpr[, c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  DD.results.compcodeR.DEDD50$fdr <- DD.results.compcodeR.DEDD50$fdr[, c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
  DD.results.compcodeR.DEDD50$tpr <- DD.results.compcodeR.DEDD50$tpr[, c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.lnHM.tmm")]
}

## Load GTEx data ####
folder <- "Results/GTEx DD lfc1 results Jan 2021"
for (i in c("2", "5", "10", "20", "50")) {
  for (j in c("blood", "muscle")) {
    assign(paste0("DD.results.", j, ".DEDD", i), 
           readRDS(here(folder, paste0("DD.lfc1.results.", j, "_", i, "_DEDD.rds"))))
  }
}
rm(i,j,folder)

## Plot FDR and TPR for each dataset ####
positions = rep(1:10)
offsets <- c(rep(0, 2), rep(1, 2), rep(2, 2), rep(3, 2), rep(4, 2))

par(mfrow=c(3,2), mar=c(2,2.5,2.5,0.5), mgp=c(3,0.7,0))

# First row - compcodeR
plot(NA, xlim=c(0.75,14.25), ylim=c(0,1), xaxt="n", yaxt="n", 
     main="False discovery rate", cex.main=2.5)
lines(c(0.5,2.5), c(0.05,0.05), col="lightgrey")
lines(c(3.5,5.5), c(0.05,0.05), col="lightgrey")
lines(c(6.5,8.5), c(0.05,0.05), col="lightgrey")
lines(c(9.5,11.5), c(0.05,0.05), col="lightgrey")
lines(c(12.5,14.5), c(0.05,0.05), col="lightgrey")
boxplot(cbind(DD.results.compcodeR.DEDD2$fdr, 
              DD.results.compcodeR.DEDD5$fdr, 
              DD.results.compcodeR.DEDD10$fdr, 
              DD.results.compcodeR.DEDD20$fdr, 
              DD.results.compcodeR.DEDD50$fdr), 
        col=col_vector[c(1,3)], pch=20, cex.axis=2, xaxt='n', 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(1.5,4.5,7.5,10.5,13.5), labels=c(2,5,10,20,50), tick=F, cex.axis=2.5, mgp=c(3,1,0))
mtext('Simulated data', cex=1.5, side=3, line=-2)
legend("topright", fill=col_vector[c(1,3)], bty='n', cex=2, ncol=1, 
       legend=c("MDSeq", "HM"))
plot(NA, xlim=c(0.75,14.25), ylim=c(0,0.42), xaxt="n", yaxt="n", 
     main="Sensitivity", cex.main=2.5)
boxplot(cbind(DD.results.compcodeR.DEDD2$tpr, 
              DD.results.compcodeR.DEDD5$tpr, 
              DD.results.compcodeR.DEDD10$tpr, 
              DD.results.compcodeR.DEDD20$tpr, 
              DD.results.compcodeR.DEDD50$tpr), 
        col=col_vector[c(1,3)], pch=20, cex.axis=2, xaxt='n', 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(1.5,4.5,7.5,10.5,13.5), labels=c(2,5,10,20,50), tick=F, cex.axis=2.5, mgp=c(3,1,0))
mtext('Simulated data', cex=1.5, side=3, line=-2)

# Second row - blood
plot(NA, xlim=c(0.75,14.25), ylim=c(0,1), xaxt="n", yaxt="n", 
     main="False discovery rate", cex.main=2.5)
lines(c(0.5,2.5), c(0.05,0.05), col="lightgrey")
lines(c(3.5,5.5), c(0.05,0.05), col="lightgrey")
lines(c(6.5,8.5), c(0.05,0.05), col="lightgrey")
lines(c(9.5,11.5), c(0.05,0.05), col="lightgrey")
lines(c(12.5,14.5), c(0.05,0.05), col="lightgrey")
boxplot(cbind(DD.results.blood.DEDD2$fdr, 
              DD.results.blood.DEDD5$fdr, 
              DD.results.blood.DEDD10$fdr, 
              DD.results.blood.DEDD20$fdr, 
              DD.results.blood.DEDD50$fdr), 
        col=col_vector[c(1,3)], pch=20, cex.axis=2, xaxt='n', 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(1.5,4.5,7.5,10.5,13.5), labels=c(2,5,10,20,50), tick=F, cex.axis=2.5, mgp=c(3,1,0))
mtext('Simulated data', cex=1.5, side=3, line=-2)
plot(NA, xlim=c(0.75,14.25), ylim=c(0,0.42), xaxt="n", yaxt="n", 
     main="Sensitivity", cex.main=2.5)
boxplot(cbind(DD.results.blood.DEDD2$tpr, 
              DD.results.blood.DEDD5$tpr, 
              DD.results.blood.DEDD10$tpr, 
              DD.results.blood.DEDD20$tpr, 
              DD.results.blood.DEDD50$tpr), 
        col=col_vector[c(1,3)], pch=20, cex.axis=2, xaxt='n', 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(1.5,4.5,7.5,10.5,13.5), labels=c(2,5,10,20,50), tick=F, cex.axis=2.5, mgp=c(3,1,0))
mtext('GTEx blood data', cex=1.5, side=3, line=-2)

# Third row - muscle
plot(NA, xlim=c(0.75,14.25), ylim=c(0,1), xaxt="n", yaxt="n", 
     main="False discovery rate", cex.main=2.5)
lines(c(0.5,2.5), c(0.05,0.05), col="lightgrey")
lines(c(3.5,5.5), c(0.05,0.05), col="lightgrey")
lines(c(6.5,8.5), c(0.05,0.05), col="lightgrey")
lines(c(9.5,11.5), c(0.05,0.05), col="lightgrey")
lines(c(12.5,14.5), c(0.05,0.05), col="lightgrey")
boxplot(cbind(DD.results.muscle.DEDD2$fdr, 
              DD.results.muscle.DEDD5$fdr, 
              DD.results.muscle.DEDD10$fdr, 
              DD.results.muscle.DEDD20$fdr, 
              DD.results.muscle.DEDD50$fdr), 
        col=col_vector[c(1,3)], pch=20, cex.axis=2, xaxt='n', 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(1.5,4.5,7.5,10.5,13.5), labels=c(2,5,10,20,50), tick=F, cex.axis=2.5, mgp=c(3,1,0))
mtext('Simulated data', cex=1.5, side=3, line=-2)
plot(NA, xlim=c(0.75,14.25), ylim=c(0,0.42), xaxt="n", yaxt="n", 
     main="Sensitivity", cex.main=2.5)
boxplot(cbind(DD.results.muscle.DEDD2$tpr, 
              DD.results.muscle.DEDD5$tpr, 
              DD.results.muscle.DEDD10$tpr, 
              DD.results.muscle.DEDD20$tpr, 
              DD.results.muscle.DEDD50$tpr), 
        col=col_vector[c(1,3)], pch=20, cex.axis=2, xaxt='n', 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(1.5,4.5,7.5,10.5,13.5), labels=c(2,5,10,20,50), tick=F, cex.axis=2.5, mgp=c(3,1,0))
mtext('GTEx muscle data', cex=1.5, side=3, line=-2)
