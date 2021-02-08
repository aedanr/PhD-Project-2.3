library(here)
library(RColorBrewer)

## Load data, create colour vector ####
folder <- "Results/compcodeR combined results Dec 2020 including diffVar"
for (i in c("DEDD2", "DEDD5", "DEDD10", "DEDD20", "DEDD50")) {
  assign(paste0("DD.results.", i), 
         readRDS(here(folder, paste0("DD.results.", i, ".rds"))))
}
rm(i,folder)

qual_col_pals = brewer.pal.info[brewer.pal.info$category == "qual" & 
                                  brewer.pal.info$colorblind == T,]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector <- col_vector[c(1:3,9,11,26)]
rm(qual_col_pals)


## Keep only TMM and one version of each method ####
# MDSeq with ZI, lnHM with log transformation.
{
  DD.results.DEDD2$fdr <- DD.results.DEDD2$fdr[, c("disp.MDSeq.zi.tmm", "EVQ.tmm", "disp.lnHM.log.tmm")]
  DD.results.DEDD2$tpr <- DD.results.DEDD2$tpr[, c("disp.MDSeq.zi.tmm", "EVQ.tmm", "disp.lnHM.log.tmm")]
  DD.results.DEDD5$fdr <- DD.results.DEDD5$fdr[, c("disp.MDSeq.zi.tmm", "EVQ.tmm", "disp.lnHM.log.tmm")]
  DD.results.DEDD5$tpr <- DD.results.DEDD5$tpr[, c("disp.MDSeq.zi.tmm", "EVQ.tmm", "disp.lnHM.log.tmm")]
  DD.results.DEDD10$fdr <- DD.results.DEDD10$fdr[, c("disp.MDSeq.zi.tmm", "EVQ.tmm", "disp.lnHM.log.tmm")]
  DD.results.DEDD10$tpr <- DD.results.DEDD10$tpr[, c("disp.MDSeq.zi.tmm", "EVQ.tmm", "disp.lnHM.log.tmm")]
  DD.results.DEDD20$fdr <- DD.results.DEDD20$fdr[, c("disp.MDSeq.zi.tmm", "EVQ.tmm", "disp.lnHM.log.tmm")]
  DD.results.DEDD20$tpr <- DD.results.DEDD20$tpr[, c("disp.MDSeq.zi.tmm", "EVQ.tmm", "disp.lnHM.log.tmm")]
  DD.results.DEDD50$fdr <- DD.results.DEDD50$fdr[, c("disp.MDSeq.zi.tmm", "EVQ.tmm", "disp.lnHM.log.tmm")]
  DD.results.DEDD50$tpr <- DD.results.DEDD50$tpr[, c("disp.MDSeq.zi.tmm", "EVQ.tmm", "disp.lnHM.log.tmm")]
}

positions = rep(1:15)
offsets <- c(rep(0, 3), rep(1, 3), rep(2, 3), rep(3, 3), rep(4, 3))

par(mfrow=c(2,1), mar=c(2,2,2,0.5), mgp=c(3,0.7,0))
plot(NA, xlim=c(0.75,19.25), ylim=c(0,1), xaxt="n", yaxt="n", 
     main="False discovery rate", cex.main=1.5)
lines(c(0.25,3.75), c(0.05,0.05), col="lightgrey")
lines(c(4.25,7.75), c(0.05,0.05), col="lightgrey")
lines(c(8.25,11.75), c(0.05,0.05), col="lightgrey")
lines(c(12.25,15.75), c(0.05,0.05), col="lightgrey")
lines(c(16.25,19.75), c(0.05,0.05), col="lightgrey")
boxplot(cbind(DD.results.DEDD2$fdr, 
              DD.results.DEDD5$fdr, 
              DD.results.DEDD10$fdr, 
              DD.results.DEDD20$fdr, 
              DD.results.DEDD50$fdr), 
        col=col_vector[1:3], pch=20, cex.axis=1.5, xaxt='n', 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(2,6,10,14,18), labels=c(2,5,10,20,50), tick=F, cex.axis=1.5)
plot(NA, xlim=c(0.75,19.25), ylim=c(0,0.43), xaxt="n", yaxt="n", 
     main="Sensitivity", cex.main=1.5)
boxplot(cbind(DD.results.DEDD2$tpr, 
              DD.results.DEDD5$tpr, 
              DD.results.DEDD10$tpr, 
              DD.results.DEDD20$tpr, 
              DD.results.DEDD50$tpr), 
        col=col_vector[1:3], pch=20, cex.axis=1.5, xaxt='n', 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(2,6,10,14,18), labels=c(2,5,10,20,50), tick=F, cex.axis=1.5)
legend("topleft", fill=col_vector[1:3], bty='n', cex=1.5, ncol=1, 
       legend=c("MDSeq", "GAMLSS", "HM"))
