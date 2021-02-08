library(here)
library(RColorBrewer)

## Load data, create colour vector ####
folder <- "Results/compcodeR DE, DD, DEDD results Feb 2020"
for (i in c("DEDD2", "DEDD5", "DEDD10", "DEDD20", "DEDD50")) {
  assign(paste0("DE.results.", i), 
         readRDS(here(folder, paste0("DE.results.", i, ".rds"))))
}
rm(i,folder)

qual_col_pals = brewer.pal.info[brewer.pal.info$category == "qual" & 
                                  brewer.pal.info$colorblind == T,]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector <- col_vector[c(1:3,9,11,26)]
rm(qual_col_pals)


## Keep only TMM and one version of each method ####
# Previously ran 2020-01-17_auc_pauc_fdr_tpr_false_discovery_plots_compcodeR_data_all_methods 
# to choose best version of each method (or default/recommended if there is one): 
# edgeR QL, DESeq2 IF, MDSeq ZI, lnHM log; only one for voom, baySeq.
{
  DE.results.DEDD2$fdr <- DE.results.DEDD2$fdr[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", "mean.lnHM.log.tmm")
  ]
  DE.results.DEDD2$tpr <- DE.results.DEDD2$tpr[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", "mean.lnHM.log.tmm")
  ]
  DE.results.DEDD5$fdr <- DE.results.DEDD5$fdr[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", "mean.lnHM.log.tmm")
  ]
  DE.results.DEDD5$tpr <- DE.results.DEDD5$tpr[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", "mean.lnHM.log.tmm")
  ]
  DE.results.DEDD10$fdr <- DE.results.DEDD10$fdr[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", "mean.lnHM.log.tmm")
  ]
  DE.results.DEDD10$tpr <- DE.results.DEDD10$tpr[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", "mean.lnHM.log.tmm")
  ]
  DE.results.DEDD20$fdr <- DE.results.DEDD20$fdr[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", "mean.lnHM.log.tmm")
  ]
  DE.results.DEDD20$tpr <- DE.results.DEDD20$tpr[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", "mean.lnHM.log.tmm")
  ]
  DE.results.DEDD50$fdr <- DE.results.DEDD50$fdr[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", "mean.lnHM.log.tmm")
  ]
  DE.results.DEDD50$tpr <- DE.results.DEDD50$tpr[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", "mean.lnHM.log.tmm")
  ]
}

positions = rep(1:30)
offsets <- c(rep(0, 6), rep(2, 6), rep(4, 6), rep(6, 6), rep(8, 6))

## Plot FDR and TPR ####
par(mfrow=c(2,1), mar=c(2,2,2,0.5), mgp=c(3,0.7,0))
plot(NA, xlim=c(1,38), ylim=c(0,0.7), xaxt="n", yaxt="n", 
     main="False discovery rate", cex.main=1.5)
lines(c(0,7), c(0.05,0.05), col="lightgrey")
lines(c(8,15), c(0.05,0.05), col="lightgrey")
lines(c(16,23), c(0.05,0.05), col="lightgrey")
lines(c(24,31), c(0.05,0.05), col="lightgrey")
lines(c(32,39), c(0.05,0.05), col="lightgrey")
boxplot(cbind(DE.results.DEDD2$fdr, 
              DE.results.DEDD5$fdr, 
              DE.results.DEDD10$fdr, 
              DE.results.DEDD20$fdr, 
              DE.results.DEDD50$fdr), 
        col=col_vector[1:6], pch=20, cex.axis=1.5, xaxt='n', ylim=c(0,0.8), 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(3.5,11.5,19.5,27.5,35.5), labels=c(2,5,10,20,50), tick=F, cex.axis=1.5)
legend("topright", fill=col_vector[1:6], bty='n', cex=1.5, ncol=1, 
       legend=c("edgeR", "DESeq2", "voom", "baySeq", "MDSeq", "HM"))
plot(NA, xlim=c(1,38), ylim=c(0,1), xaxt="n", yaxt="n", 
     main="Sensitivity", cex.main=1.5)
boxplot(cbind(DE.results.DEDD2$tpr, 
              DE.results.DEDD5$tpr, 
              DE.results.DEDD10$tpr, 
              DE.results.DEDD20$tpr, 
              DE.results.DEDD50$tpr), 
        col=col_vector[1:6], pch=20, cex.axis=1.5, xaxt='n', ylim=c(0,0.8), 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(3.5,11.5,19.5,27.5,35.5), labels=c(2,5,10,20,50), tick=F, cex.axis=1.5)
