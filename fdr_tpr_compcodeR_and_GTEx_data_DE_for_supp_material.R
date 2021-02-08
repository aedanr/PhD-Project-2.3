library(here)
library(RColorBrewer)

# Make colour vector
qual_col_pals = brewer.pal.info[brewer.pal.info$category == "qual" & 
                                  brewer.pal.info$colorblind == T,]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# pie(rep(1,length(col_vector)), col=col_vector)
col_vector <- col_vector[c(1:3,9,11,26)]
rm(qual_col_pals)

## Load compcodeR data and keep only edgeR QL, DESeq2 IF, voom, baySeq, MDSeq ZI, lnHM log ####
folder <- "Results/compcodeR DE, DD, DEDD results Feb 2020"
for (i in c("DEDD2", "DEDD5", "DEDD10", "DEDD20", "DEDD50")) {
  assign(paste0("DE.results.compcodeR.", i), 
         readRDS(here(folder, paste0("DE.results.", i, ".rds"))))
}
rm(i,folder)

{
  DE.results.compcodeR.DEDD2$fdr <- DE.results.compcodeR.DEDD2$fdr[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", "mean.lnHM.log.tmm")
  ]
  DE.results.compcodeR.DEDD2$tpr <- DE.results.compcodeR.DEDD2$tpr[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", "mean.lnHM.log.tmm")
  ]
  DE.results.compcodeR.DEDD5$fdr <- DE.results.compcodeR.DEDD5$fdr[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", "mean.lnHM.log.tmm")
  ]
  DE.results.compcodeR.DEDD5$tpr <- DE.results.compcodeR.DEDD5$tpr[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", "mean.lnHM.log.tmm")
  ]
  DE.results.compcodeR.DEDD10$fdr <- DE.results.compcodeR.DEDD10$fdr[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", "mean.lnHM.log.tmm")
  ]
  DE.results.compcodeR.DEDD10$tpr <- DE.results.compcodeR.DEDD10$tpr[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", "mean.lnHM.log.tmm")
  ]
  DE.results.compcodeR.DEDD20$fdr <- DE.results.compcodeR.DEDD20$fdr[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", "mean.lnHM.log.tmm")
  ]
  DE.results.compcodeR.DEDD20$tpr <- DE.results.compcodeR.DEDD20$tpr[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", "mean.lnHM.log.tmm")
  ]
  DE.results.compcodeR.DEDD50$fdr <- DE.results.compcodeR.DEDD50$fdr[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", "mean.lnHM.log.tmm")
  ]
  DE.results.compcodeR.DEDD50$tpr <- DE.results.compcodeR.DEDD50$tpr[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", "mean.lnHM.log.tmm")
  ]
}

## Load GTEx data, rearrange and keep only edgeR QL, DESeq2 IF, voom, baySeq, MDSeq ZI, lnHM log ####
folder <- "Results/GTEx combined results Sept 2020"
for (i in c("2", "5", "10", "20", "50")) {
  assign(paste0("DE.results.DEDD", i), 
         readRDS(here(folder, paste0("DE.results.DEDD", i, ".rds"))))
}
rm(i,folder)

for (metric in c("fdr", "tpr")) {
  for (size in c("2", "5", "10", "20", "50")) {
    assign(
      paste0(metric, "_blood_", size), 
      get(paste0("DE.results.DEDD", size))[[metric]][
        , c("blood_eR.ql", "blood_DES.if", "blood_voom", "blood_baySeq", "blood_MD.zi", "blood_lnHM.log")
      ]
    )
    assign(
      paste0(metric, "_muscle_", size), 
      get(paste0("DE.results.DEDD", size))[[metric]][
        , c("muscle_eR.ql", "muscle_DES.if", "muscle_voom", "muscle_baySeq", "muscle_MD.zi", "muscle_lnHM.log")
      ]
    )
  }
}
rm(metric, size)


## Plot FDR and TPR for each dataset ####
positions = rep(1:30)
offsets <- c(rep(0, 6), rep(2, 6), rep(4, 6), rep(6, 6), rep(8, 6))

par(mfrow=c(3,2), mar=c(2,2.5,2.5,0.5), mgp=c(3,0.7,0))

# First row - compcodeR
plot(NA, xlim=c(1,38), ylim=c(0,1.1), xaxt="n", yaxt="n", 
     main="False discovery rate", cex.main=2.5)
lines(c(0,7), c(0.05,0.05), col="lightgrey")
lines(c(8,15), c(0.05,0.05), col="lightgrey")
lines(c(16,23), c(0.05,0.05), col="lightgrey")
lines(c(24,31), c(0.05,0.05), col="lightgrey")
lines(c(32,39), c(0.05,0.05), col="lightgrey")
boxplot(cbind(DE.results.compcodeR.DEDD2$fdr, 
              DE.results.compcodeR.DEDD5$fdr, 
              DE.results.compcodeR.DEDD10$fdr, 
              DE.results.compcodeR.DEDD20$fdr, 
              DE.results.compcodeR.DEDD50$fdr), 
        col=col_vector[1:6], pch=20, cex.axis=2, xaxt='n', border=col_vector[1:6], 
        at=positions + offsets, lwd=0.5, cex=0.5,
        add=T)
axis(side=1, at=c(3.5,11.5,19.5,27.5,35.5), labels=c(2,5,10,20,50), tick=F, cex.axis=2.5, mgp=c(3,1,0))
mtext('Simulated data', cex=1.5, side=3, line=-2)
legend("topright", fill=col_vector[1:6], bty='n', cex=2, ncol=1, 
       legend=c("edgeR", "DESeq2", "voom", "baySeq", "MDSeq", "HM"))
plot(NA, xlim=c(1,38), ylim=c(0,1), xaxt="n", yaxt="n", 
     main="Sensitivity", cex.main=2.5)
boxplot(cbind(DE.results.compcodeR.DEDD2$tpr, 
              DE.results.compcodeR.DEDD5$tpr, 
              DE.results.compcodeR.DEDD10$tpr, 
              DE.results.compcodeR.DEDD20$tpr, 
              DE.results.compcodeR.DEDD50$tpr), 
        col=col_vector[1:6], pch=20, cex.axis=2, xaxt='n', border=col_vector[1:6], 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(3.5,11.5,19.5,27.5,35.5), labels=c(2,5,10,20,50), tick=F, cex.axis=2.5, mgp=c(3,1,0))
mtext('Simulated data', cex=1.5, side=3, line=-2)

# Second row - blood
plot(NA, xlim=c(1,38), ylim=c(0,1.1), xaxt="n", yaxt="n")
lines(c(0,7), c(0.05,0.05), col="lightgrey")
lines(c(8,15), c(0.05,0.05), col="lightgrey")
lines(c(16,23), c(0.05,0.05), col="lightgrey")
lines(c(24,31), c(0.05,0.05), col="lightgrey")
lines(c(32,39), c(0.05,0.05), col="lightgrey")
boxplot(cbind(fdr_blood_2, fdr_blood_5, fdr_blood_10, fdr_blood_20, fdr_blood_50), 
        col=col_vector[1:6], pch=20, cex.axis=2, xaxt='n', border=col_vector[1:6], 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(3.5,11.5,19.5,27.5,35.5), labels=c(2,5,10,20,50), tick=F, cex.axis=2.5, mgp=c(3,1,0))
mtext('GTEx blood data', cex=1.5, side=3, line=-2)
plot(NA, xlim=c(1,38), ylim=c(0,1), xaxt="n", yaxt="n")
boxplot(cbind(tpr_blood_2, tpr_blood_5, tpr_blood_10, tpr_blood_20, tpr_blood_50), 
        col=col_vector[1:6], pch=20, cex.axis=2, xaxt='n', border=col_vector[1:6], 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(3.5,11.5,19.5,27.5,35.5), labels=c(2,5,10,20,50), tick=F, cex.axis=2.5, mgp=c(3,1,0))
mtext('GTEx blood data', cex=1.5, side=3, line=-2)

# Third row - muscle
plot(NA, xlim=c(1,38), ylim=c(0,1.1), xaxt="n", yaxt="n")
lines(c(0,7), c(0.05,0.05), col="lightgrey")
lines(c(8,15), c(0.05,0.05), col="lightgrey")
lines(c(16,23), c(0.05,0.05), col="lightgrey")
lines(c(24,31), c(0.05,0.05), col="lightgrey")
lines(c(32,39), c(0.05,0.05), col="lightgrey")
boxplot(cbind(fdr_muscle_2, fdr_muscle_5, fdr_muscle_10, fdr_muscle_20, fdr_muscle_50), 
        col=col_vector[1:6], pch=20, cex.axis=2, xaxt='n', border=col_vector[1:6], 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(3.5,11.5,19.5,27.5,35.5), labels=c(2,5,10,20,50), tick=F, cex.axis=2.5, mgp=c(3,1,0))
mtext('GTEx muscle data', cex=1.5, side=3, line=-2)
plot(NA, xlim=c(1,38), ylim=c(0,1), xaxt="n", yaxt="n")
boxplot(cbind(tpr_muscle_2, tpr_muscle_5, tpr_muscle_10, tpr_muscle_20, tpr_muscle_50), 
        col=col_vector[1:6], pch=20, cex.axis=2, xaxt='n', border=col_vector[1:6], 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(3.5,11.5,19.5,27.5,35.5), labels=c(2,5,10,20,50), tick=F, cex.axis=2.5, mgp=c(3,1,0))
mtext('GTEx muscle data', cex=1.5, side=3, line=-2)
