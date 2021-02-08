library(here)
library(RColorBrewer)

## Load data, create colour vector ####
folder <- "Results/compcodeR DE, DD, DEDD results Feb 2020"
for (i in c("DEDD2", "DEDD5", "DEDD10", "DEDD20", "DEDD50")) {
  assign(paste0("DE.results.", i), 
         readRDS(here(folder, paste0("DE.results.", i, ".rds"))))
}
for (i in c("DE2", "DE5", "DE10", "DE20", "DE50")) {
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
# Ran 2020-01-17_auc_pauc_fdr_tpr_false_discovery_plots_compcodeR_data_all_methods to choose 
# best version of each method.
# edgeR: LR > ET > QL for AUC; same for pAUC except QL best for 2 samples; only QL has good 
# FDR control for 5 and 10 samples; QL has worst TPR; QL best FDR curves for 2 samples, 
# otherwise indistinguishable. QL seems to be generally recommended in the user guide, but 
# only recommended ahead of ET ("classic") because of additional functionality that isn't 
# relevant here. QL preferred over LR because reflects uncertainty in dispersion estimates. 
# Since more interested in results for small samples for DE, use QL.
# DESeq2: with IF since that's the default and recommended.
# DSS: lfdr gives much better FDRs for up to 10 samples, but very conservative for 20 and 
# 50; couldn't get it to run for some datasets with 50 samples, so don't include unless can 
# get it to work (may be an updated version now, unlikely but should try; hasn't been an 
# update on the issue I opened on github, and local FDR function which was causing the issue 
# hasn't been updated). Could also include but note results missing for some datasets. See 
# what situation is with GTEx data.
# MDSeq: with ZI since it's the default.
# HM: lnHM, log-transformed since lnHM generally gives better results than expHM, and not 
# much difference with and without log transformation, and it's needed if I'm going to 
# include tests at a given LFC. Might still test expHM v lnHM on same size, larger posterior 
# samples though.
{
  DE.results.DE2$pauc <- DE.results.DE2$pauc[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  DE.results.DE2$auc <- DE.results.DE2$auc[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  DE.results.DE2$fpr <- DE.results.DE2$fpr[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  DE.results.DE2$fdr <- DE.results.DE2$fdr[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  DE.results.DE2$tpr <- DE.results.DE2$tpr[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  DE.results.DE2$mean.fdr <- DE.results.DE2$mean.fdr[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  DE.results.DE2$mean.discoveries <- DE.results.DE2$mean.discoveries[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  
  DE.results.DE5$pauc <- DE.results.DE5$pauc[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  DE.results.DE5$auc <- DE.results.DE5$auc[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  DE.results.DE5$fpr <- DE.results.DE5$fpr[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  DE.results.DE5$fdr <- DE.results.DE5$fdr[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  DE.results.DE5$tpr <- DE.results.DE5$tpr[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  DE.results.DE5$mean.fdr <- DE.results.DE5$mean.fdr[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  DE.results.DE5$mean.discoveries <- DE.results.DE5$mean.discoveries[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  
  DE.results.DE10$pauc <- DE.results.DE10$pauc[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  DE.results.DE10$auc <- DE.results.DE10$auc[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  DE.results.DE10$fpr <- DE.results.DE10$fpr[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  DE.results.DE10$fdr <- DE.results.DE10$fdr[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  DE.results.DE10$tpr <- DE.results.DE10$tpr[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  DE.results.DE10$mean.fdr <- DE.results.DE10$mean.fdr[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  DE.results.DE10$mean.discoveries <- DE.results.DE10$mean.discoveries[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  
  DE.results.DE20$pauc <- DE.results.DE20$pauc[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  DE.results.DE20$auc <- DE.results.DE20$auc[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  DE.results.DE20$fpr <- DE.results.DE20$fpr[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  DE.results.DE20$fdr <- DE.results.DE20$fdr[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  DE.results.DE20$tpr <- DE.results.DE20$tpr[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  DE.results.DE20$mean.fdr <- DE.results.DE20$mean.fdr[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  DE.results.DE20$mean.discoveries <- DE.results.DE20$mean.discoveries[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  
  DE.results.DE50$pauc <- DE.results.DE50$pauc[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  DE.results.DE50$auc <- DE.results.DE50$auc[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  DE.results.DE50$fpr <- DE.results.DE50$fpr[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  DE.results.DE50$fdr <- DE.results.DE50$fdr[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  DE.results.DE50$tpr <- DE.results.DE50$tpr[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  DE.results.DE50$mean.fdr <- DE.results.DE50$mean.fdr[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  DE.results.DE50$mean.discoveries <- DE.results.DE50$mean.discoveries[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  
  DE.results.DEDD2$pauc <- DE.results.DEDD2$pauc[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  DE.results.DEDD2$auc <- DE.results.DEDD2$auc[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  DE.results.DEDD2$fpr <- DE.results.DEDD2$fpr[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  DE.results.DEDD2$fdr <- DE.results.DEDD2$fdr[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  DE.results.DEDD2$tpr <- DE.results.DEDD2$tpr[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  DE.results.DEDD2$mean.fdr <- DE.results.DEDD2$mean.fdr[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  DE.results.DEDD2$mean.discoveries <- DE.results.DEDD2$mean.discoveries[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  
  DE.results.DEDD5$pauc <- DE.results.DEDD5$pauc[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  DE.results.DEDD5$auc <- DE.results.DEDD5$auc[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  DE.results.DEDD5$fpr <- DE.results.DEDD5$fpr[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  DE.results.DEDD5$fdr <- DE.results.DEDD5$fdr[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  DE.results.DEDD5$tpr <- DE.results.DEDD5$tpr[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  DE.results.DEDD5$mean.fdr <- DE.results.DEDD5$mean.fdr[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  DE.results.DEDD5$mean.discoveries <- DE.results.DEDD5$mean.discoveries[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  
  DE.results.DEDD10$pauc <- DE.results.DEDD10$pauc[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  DE.results.DEDD10$auc <- DE.results.DEDD10$auc[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  DE.results.DEDD10$fpr <- DE.results.DEDD10$fpr[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  DE.results.DEDD10$fdr <- DE.results.DEDD10$fdr[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  DE.results.DEDD10$tpr <- DE.results.DEDD10$tpr[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  DE.results.DEDD10$mean.fdr <- DE.results.DEDD10$mean.fdr[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  DE.results.DEDD10$mean.discoveries <- DE.results.DEDD10$mean.discoveries[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  
  DE.results.DEDD20$pauc <- DE.results.DEDD20$pauc[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  DE.results.DEDD20$auc <- DE.results.DEDD20$auc[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  DE.results.DEDD20$fpr <- DE.results.DEDD20$fpr[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  DE.results.DEDD20$fdr <- DE.results.DEDD20$fdr[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  DE.results.DEDD20$tpr <- DE.results.DEDD20$tpr[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  DE.results.DEDD20$mean.fdr <- DE.results.DEDD20$mean.fdr[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  DE.results.DEDD20$mean.discoveries <- DE.results.DEDD20$mean.discoveries[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  
  DE.results.DEDD50$pauc <- DE.results.DEDD50$pauc[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  DE.results.DEDD50$auc <- DE.results.DEDD50$auc[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  DE.results.DEDD50$fpr <- DE.results.DEDD50$fpr[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  DE.results.DEDD50$fdr <- DE.results.DEDD50$fdr[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  DE.results.DEDD50$tpr <- DE.results.DEDD50$tpr[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  DE.results.DEDD50$mean.fdr <- DE.results.DEDD50$mean.fdr[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
  DE.results.DEDD50$mean.discoveries <- DE.results.DEDD50$mean.discoveries[
    , c("edgeR.ql.tmm", "DESeq2.if.tmm", "voom.tmm", "baySeq.tmm", "mean.MDSeq.zi.tmm", 
        "mean.lnHM.log.tmm")
  ]
}

positions = rep(1:30)
offsets <- c(rep(0, 6), rep(1, 6), rep(2, 6), rep(3, 6), rep(4, 6))

### AUC ####
par(mfrow=c(2,1), mar=c(2.5,2.5,2.5,0.5), mgp=c(3,0.7,0))
plot(NA, xlim=c(1,34), ylim=c(0.75,1), xaxt="n", yaxt="n", 
     main="AUC - differences in mean only", cex.main=2)
boxplot(cbind(DE.results.DE2$auc, 
              DE.results.DE5$auc, 
              DE.results.DE10$auc, 
              DE.results.DE20$auc, 
              DE.results.DE50$auc), 
        col=col_vector[1:6], pch=20, cex.axis=2, xaxt='n', ylim=c(0,0.8), 
        at=positions + offsets, lwd=0.2, cex=0.5, 
        add=T)
axis(side=1, at=c(3.5,10.5,17.5,24.5,31.5), labels=c(2,5,10,20,50), tick=F, cex.axis=2)
legend("bottomright", fill=col_vector[1:6], bty='n', cex=1.5, ncol=1, 
       legend=c("edgeR", "DESeq2", "voom", "baySeq", "MDSeq", "HM"))
plot(NA, xlim=c(1,34), ylim=c(0.75,1), xaxt="n", yaxt="n", 
     main="AUC - differences in mean and dispersion", cex.main=2)
boxplot(cbind(DE.results.DEDD2$auc, 
              DE.results.DEDD5$auc, 
              DE.results.DEDD10$auc, 
              DE.results.DEDD20$auc, 
              DE.results.DEDD50$auc), 
        col=col_vector[1:6], pch=20, cex.axis=2, xaxt='n', ylim=c(0,0.8), 
        at=positions + offsets, lwd=0.2, cex=0.5, 
        add=T)
axis(side=1, at=c(3.5,10.5,17.5,24.5,31.5), labels=c(2,5,10,20,50), tick=F, cex.axis=2)

rbind(colMeans(DE.results.DE2$auc), 
      colMeans(DE.results.DE5$auc), 
      colMeans(DE.results.DE10$auc), 
      colMeans(DE.results.DE20$auc), 
      colMeans(DE.results.DE50$auc))
rbind(colMeans(DE.results.DEDD2$auc), 
      colMeans(DE.results.DEDD5$auc), 
      colMeans(DE.results.DEDD10$auc), 
      colMeans(DE.results.DEDD20$auc), 
      colMeans(DE.results.DEDD50$auc))

### pAUC ####
par(mfrow=c(2,1), mar=c(2.5,2.5,2.5,0.5), mgp=c(3,0.7,0))
plot(NA, xlim=c(1,34), ylim=c(0.009,0.05), xaxt="n", yaxt="n", 
     main="Partial AUC - differences in mean only", cex.main=2)
boxplot(cbind(DE.results.DE2$pauc, 
              DE.results.DE5$pauc, 
              DE.results.DE10$pauc, 
              DE.results.DE20$pauc, 
              DE.results.DE50$pauc), 
        col=col_vector[1:6], pch=20, cex.axis=2, xaxt='n', ylim=c(0,0.8), 
        at=positions + offsets, lwd=0.2, cex=0.5, 
        add=T)
axis(side=1, at=c(3.5,10.5,17.5,24.5,31.5), labels=c(2,5,10,20,50), tick=F, cex.axis=2)
legend("bottomright", fill=col_vector[1:6], bty='n', cex=1.5, ncol=1, 
       legend=c("edgeR", "DESeq2", "voom", "baySeq", "MDSeq", "HM"))
plot(NA, xlim=c(1,34), ylim=c(0.009,0.05), xaxt="n", yaxt="n", 
     main="Partial AUC - differences in mean and dispersion", cex.main=2)
boxplot(cbind(DE.results.DEDD2$pauc, 
              DE.results.DEDD5$pauc, 
              DE.results.DEDD10$pauc, 
              DE.results.DEDD20$pauc, 
              DE.results.DEDD50$pauc), 
        col=col_vector[1:6], pch=20, cex.axis=2, xaxt='n', ylim=c(0,0.8), 
        at=positions + offsets, lwd=0.2, cex=0.5, 
        add=T)
axis(side=1, at=c(3.5,10.5,17.5,24.5,31.5), labels=c(2,5,10,20,50), tick=F, cex.axis=2)

rbind(colMeans(DE.results.DE2$pauc), 
      colMeans(DE.results.DE5$pauc), 
      colMeans(DE.results.DE10$pauc), 
      colMeans(DE.results.DE20$pauc), 
      colMeans(DE.results.DE50$pauc))
rbind(colMeans(DE.results.DEDD2$pauc), 
      colMeans(DE.results.DEDD5$pauc), 
      colMeans(DE.results.DEDD10$pauc), 
      colMeans(DE.results.DEDD20$pauc), 
      colMeans(DE.results.DEDD50$pauc))

### FDR ####
par(mfrow=c(2,1), mar=c(2.5,2.5,2.5,0.5), mgp=c(3,0.7,0))
plot(NA, xlim=c(1,34), ylim=c(0,0.8), xaxt="n", yaxt="n", 
     main="FDR - differences in mean only", cex.main=2)
lines(c(0.5,6.5), c(0.05,0.05), col="lightgrey")
lines(c(7.5,13.5), c(0.05,0.05), col="lightgrey")
lines(c(14.5,20.5), c(0.05,0.05), col="lightgrey")
lines(c(21.5,27.5), c(0.05,0.05), col="lightgrey")
lines(c(28.5,34.5), c(0.05,0.05), col="lightgrey")
boxplot(cbind(DE.results.DE2$fdr, 
              DE.results.DE5$fdr, 
              DE.results.DE10$fdr, 
              DE.results.DE20$fdr, 
              DE.results.DE50$fdr), 
        col=col_vector[1:6], pch=20, cex.axis=2, xaxt='n', ylim=c(0,0.8), 
        at=positions + offsets, lwd=0.2, cex=0.5, 
        add=T)
axis(side=1, at=c(3.5,10.5,17.5,24.5,31.5), labels=c(2,5,10,20,50), tick=F, cex.axis=2)
plot(NA, xlim=c(1,34), ylim=c(0,0.8), xaxt="n", yaxt="n", 
     main="FDR - differences in mean and dispersion", cex.main=2)
lines(c(0.5,6.5), c(0.05,0.05), col="lightgrey")
lines(c(7.5,13.5), c(0.05,0.05), col="lightgrey")
lines(c(14.5,20.5), c(0.05,0.05), col="lightgrey")
lines(c(21.5,27.5), c(0.05,0.05), col="lightgrey")
lines(c(28.5,34.5), c(0.05,0.05), col="lightgrey")
boxplot(cbind(DE.results.DEDD2$fdr, 
              DE.results.DEDD5$fdr, 
              DE.results.DEDD10$fdr, 
              DE.results.DEDD20$fdr, 
              DE.results.DEDD50$fdr), 
        col=col_vector[1:6], pch=20, cex.axis=2, xaxt='n', ylim=c(0,0.8), 
        at=positions + offsets, lwd=0.2, cex=0.5, 
        add=T)
axis(side=1, at=c(3.5,10.5,17.5,24.5,31.5), labels=c(2,5,10,20,50), tick=F, cex.axis=2)
legend("topright", fill=col_vector[1:6], bty='n', cex=1.5, ncol=1, 
       legend=c("edgeR", "DESeq2", "voom", "baySeq", "MDSeq", "HM"))

rbind(colMeans(DE.results.DE2$fdr), 
      colMeans(DE.results.DE5$fdr), 
      colMeans(DE.results.DE10$fdr), 
      colMeans(DE.results.DE20$fdr), 
      colMeans(DE.results.DE50$fdr))
rbind(colMeans(DE.results.DEDD2$fdr), 
      colMeans(DE.results.DEDD5$fdr), 
      colMeans(DE.results.DEDD10$fdr), 
      colMeans(DE.results.DEDD20$fdr), 
      colMeans(DE.results.DEDD50$fdr))

### TPR ####
par(mfrow=c(2,1), mar=c(2.5,2.5,2.5,0.5), mgp=c(3,0.7,0))
plot(NA, xlim=c(1,34), ylim=c(0,1), xaxt="n", yaxt="n", 
     main="TPR - differences in mean only", cex.main=2)
boxplot(cbind(DE.results.DE2$tpr, 
              DE.results.DE5$tpr, 
              DE.results.DE10$tpr, 
              DE.results.DE20$tpr, 
              DE.results.DE50$tpr), 
        col=col_vector[1:6], pch=20, cex.axis=2, xaxt='n', ylim=c(0,0.8), 
        at=positions + offsets, lwd=0.2, cex=0.5, 
        add=T)
axis(side=1, at=c(3.5,10.5,17.5,24.5,31.5), labels=c(2,5,10,20,50), tick=F, cex.axis=2)
legend("bottomright", fill=col_vector[1:6], bty='n', cex=1.5, ncol=1, 
       legend=c("edgeR", "DESeq2", "voom", "baySeq", "MDSeq", "HM"))
plot(NA, xlim=c(1,34), ylim=c(0,1), xaxt="n", yaxt="n", 
     main="TPR - differences in mean and dispersion", cex.main=2)
boxplot(cbind(DE.results.DEDD2$tpr, 
              DE.results.DEDD5$tpr, 
              DE.results.DEDD10$tpr, 
              DE.results.DEDD20$tpr, 
              DE.results.DEDD50$tpr), 
        col=col_vector[1:6], pch=20, cex.axis=2, xaxt='n', ylim=c(0,0.8), 
        at=positions + offsets, lwd=0.2, cex=0.5, 
        add=T)
axis(side=1, at=c(3.5,10.5,17.5,24.5,31.5), labels=c(2,5,10,20,50), tick=F, cex.axis=2)

rbind(colMeans(DE.results.DE2$tpr), 
      colMeans(DE.results.DE5$tpr), 
      colMeans(DE.results.DE10$tpr), 
      colMeans(DE.results.DE20$tpr), 
      colMeans(DE.results.DE50$tpr))
rbind(colMeans(DE.results.DEDD2$tpr), 
      colMeans(DE.results.DEDD5$tpr), 
      colMeans(DE.results.DEDD10$tpr), 
      colMeans(DE.results.DEDD20$tpr), 
      colMeans(DE.results.DEDD50$tpr))


### FDR curves ####
par(mfrow=c(2,5), mar=c(1.5,3.5,2.5,1), mgp=c(3,1,0))

plot(DE.results.DE2$mean.discoveries$edgeR.ql.tmm, DE.results.DE2$mean.fdr$edgeR.ql.tmm, 
     type='l', ylim=c(0,0.6), xlim=c(0,2000), lwd=2, col=col_vector[1], cex.axis=3, 
     xlab=NA, ylab=NA, xaxt='n', main=paste0("2"), cex.main=3)
axis(side=1, at=c(0,500,1000,1500,2000), tick=T, cex.axis=3, mgp=c(3,2,0))
text("Differences in\nmean only", x=0, y=0.1, cex=2.5, adj=0)
lines(DE.results.DE2$mean.discoveries$DESeq2.if.tmm, DE.results.DE2$mean.fdr$DESeq2.if.tmm, 
      lwd=2, col=col_vector[2], cex.axis=3)
lines(DE.results.DE2$mean.discoveries$voom.tmm, DE.results.DE2$mean.fdr$voom.tmm, 
      lwd=2, col=col_vector[3], cex.axis=3)
lines(DE.results.DE2$mean.discoveries$baySeq.tmm, DE.results.DE2$mean.fdr$baySeq.tmm, 
      lwd=2, col=col_vector[4], cex.axis=3)
lines(DE.results.DE2$mean.discoveries$mean.MDSeq.zi.tmm, DE.results.DE2$mean.fdr$mean.MDSeq.zi.tmm, 
      lwd=2, col=col_vector[5], cex.axis=3)
lines(DE.results.DE2$mean.discoveries$mean.lnHM.log.tmm, DE.results.DE2$mean.fdr$mean.lnHM.log.tmm, 
      lwd=2, col=col_vector[6], cex.axis=3)
abline(h=0.05, col='lightgrey')

plot(DE.results.DE5$mean.discoveries$edgeR.ql.tmm, DE.results.DE5$mean.fdr$edgeR.ql.tmm, 
     type='l', ylim=c(0,0.6), xlim=c(0,2000), lwd=2, col=col_vector[1], cex.axis=3, 
     xlab=NA, ylab=NA, xaxt='n', yaxt='n', main=paste0("5"), cex.main=3)
axis(side=1, at=c(0,500,1000,1500,2000), tick=T, cex.axis=3, mgp=c(3,2,0))
lines(DE.results.DE5$mean.discoveries$DESeq2.if.tmm, DE.results.DE5$mean.fdr$DESeq2.if.tmm, 
      lwd=2, col=col_vector[2], cex.axis=3)
lines(DE.results.DE5$mean.discoveries$voom.tmm, DE.results.DE5$mean.fdr$voom.tmm, 
      lwd=2, col=col_vector[3], cex.axis=3)
lines(DE.results.DE5$mean.discoveries$baySeq.tmm, DE.results.DE5$mean.fdr$baySeq.tmm, 
      lwd=2, col=col_vector[4], cex.axis=3)
lines(DE.results.DE5$mean.discoveries$mean.MDSeq.zi.tmm, DE.results.DE5$mean.fdr$mean.MDSeq.zi.tmm, 
      lwd=2, col=col_vector[5], cex.axis=3)
lines(DE.results.DE5$mean.discoveries$mean.lnHM.log.tmm, DE.results.DE5$mean.fdr$mean.lnHM.log.tmm, 
      lwd=2, col=col_vector[6], cex.axis=3)
abline(h=0.05, col='lightgrey')

plot(DE.results.DE10$mean.discoveries$edgeR.ql.tmm, DE.results.DE10$mean.fdr$edgeR.ql.tmm, 
     type='l', ylim=c(0,0.6), xlim=c(0,2000), lwd=2, col=col_vector[1], cex.axis=3, 
     xlab=NA, ylab=NA, xaxt='n', yaxt='n', main=paste0("10"), cex.main=3)
axis(side=1, at=c(0,500,1000,1500,2000), tick=T, cex.axis=3, mgp=c(3,2,0))
lines(DE.results.DE10$mean.discoveries$DESeq2.if.tmm, DE.results.DE10$mean.fdr$DESeq2.if.tmm, 
      lwd=2, col=col_vector[2], cex.axis=3)
lines(DE.results.DE10$mean.discoveries$voom.tmm, DE.results.DE10$mean.fdr$voom.tmm, 
      lwd=2, col=col_vector[3], cex.axis=3)
lines(DE.results.DE10$mean.discoveries$baySeq.tmm, DE.results.DE10$mean.fdr$baySeq.tmm, 
      lwd=2, col=col_vector[4], cex.axis=3)
lines(DE.results.DE10$mean.discoveries$mean.MDSeq.zi.tmm, DE.results.DE10$mean.fdr$mean.MDSeq.zi.tmm, 
      lwd=2, col=col_vector[5], cex.axis=3)
lines(DE.results.DE10$mean.discoveries$mean.lnHM.log.tmm, DE.results.DE10$mean.fdr$mean.lnHM.log.tmm, 
      lwd=2, col=col_vector[6], cex.axis=3)
abline(h=0.05, col='lightgrey')

plot(DE.results.DE20$mean.discoveries$edgeR.ql.tmm, DE.results.DE20$mean.fdr$edgeR.ql.tmm, 
     type='l', ylim=c(0,0.6), xlim=c(0,2000), lwd=2, col=col_vector[1], cex.axis=3, 
     xlab=NA, ylab=NA, xaxt='n', yaxt='n', main=paste0("20"), cex.main=3)
axis(side=1, at=c(0,500,1000,1500,2000), tick=T, cex.axis=3, mgp=c(3,2,0))
lines(DE.results.DE20$mean.discoveries$DESeq2.if.tmm, DE.results.DE20$mean.fdr$DESeq2.if.tmm, 
      lwd=2, col=col_vector[2], cex.axis=3)
lines(DE.results.DE20$mean.discoveries$voom.tmm, DE.results.DE20$mean.fdr$voom.tmm, 
      lwd=2, col=col_vector[3], cex.axis=3)
lines(DE.results.DE20$mean.discoveries$baySeq.tmm, DE.results.DE20$mean.fdr$baySeq.tmm, 
      lwd=2, col=col_vector[4], cex.axis=3)
lines(DE.results.DE20$mean.discoveries$mean.MDSeq.zi.tmm, DE.results.DE20$mean.fdr$mean.MDSeq.zi.tmm, 
      lwd=2, col=col_vector[5], cex.axis=3)
lines(DE.results.DE20$mean.discoveries$mean.lnHM.log.tmm, DE.results.DE20$mean.fdr$mean.lnHM.log.tmm, 
      lwd=2, col=col_vector[6], cex.axis=3)
abline(h=0.05, col='lightgrey')

plot(DE.results.DE50$mean.discoveries$edgeR.ql.tmm, DE.results.DE50$mean.fdr$edgeR.ql.tmm, 
     type='l', ylim=c(0,0.6), xlim=c(0,2000), lwd=2, col=col_vector[1], cex.axis=3, 
     xlab=NA, ylab=NA, xaxt='n', yaxt='n', main=paste0("50"), cex.main=3)
axis(side=1, at=c(0,500,1000,1500,2000), tick=T, cex.axis=3, mgp=c(3,2,0))
lines(DE.results.DE50$mean.discoveries$DESeq2.if.tmm, DE.results.DE50$mean.fdr$DESeq2.if.tmm, 
      lwd=2, col=col_vector[2], cex.axis=3)
lines(DE.results.DE50$mean.discoveries$voom.tmm, DE.results.DE50$mean.fdr$voom.tmm, 
      lwd=2, col=col_vector[3], cex.axis=3)
lines(DE.results.DE50$mean.discoveries$baySeq.tmm, DE.results.DE50$mean.fdr$baySeq.tmm, 
      lwd=2, col=col_vector[4], cex.axis=3)
lines(DE.results.DE50$mean.discoveries$mean.MDSeq.zi.tmm, DE.results.DE50$mean.fdr$mean.MDSeq.zi.tmm, 
      lwd=2, col=col_vector[5], cex.axis=3)
lines(DE.results.DE50$mean.discoveries$mean.lnHM.log.tmm, DE.results.DE50$mean.fdr$mean.lnHM.log.tmm, 
      lwd=2, col=col_vector[6], cex.axis=3)
abline(h=0.05, col='lightgrey')

plot(DE.results.DEDD2$mean.discoveries$edgeR.ql.tmm, DE.results.DEDD2$mean.fdr$edgeR.ql.tmm, 
     type='l', ylim=c(0,0.6), xlim=c(0,2000), lwd=2, col=col_vector[1], cex.axis=3, 
     xlab=NA, ylab=NA, xaxt='n')
text("Differences in\nmean and dispersion", x=0, y=0.1, cex=2.5, adj=0)
lines(DE.results.DEDD2$mean.discoveries$DESeq2.if.tmm, DE.results.DEDD2$mean.fdr$DESeq2.if.tmm, 
      lwd=2, col=col_vector[2], cex.axis=3)
lines(DE.results.DEDD2$mean.discoveries$voom.tmm, DE.results.DEDD2$mean.fdr$voom.tmm, 
      lwd=2, col=col_vector[3], cex.axis=3)
lines(DE.results.DEDD2$mean.discoveries$baySeq.tmm, DE.results.DEDD2$mean.fdr$baySeq.tmm, 
      lwd=2, col=col_vector[4], cex.axis=3)
lines(DE.results.DEDD2$mean.discoveries$mean.MDSeq.zi.tmm, DE.results.DEDD2$mean.fdr$mean.MDSeq.zi.tmm, 
      lwd=2, col=col_vector[5], cex.axis=3)
lines(DE.results.DEDD2$mean.discoveries$mean.lnHM.log.tmm, DE.results.DEDD2$mean.fdr$mean.lnHM.log.tmm, 
      lwd=2, col=col_vector[6], cex.axis=3)
abline(h=0.05, col='lightgrey')

plot(DE.results.DEDD5$mean.discoveries$edgeR.ql.tmm, DE.results.DEDD5$mean.fdr$edgeR.ql.tmm, 
     type='l', ylim=c(0,0.6), xlim=c(0,2000), lwd=2, col=col_vector[1], cex.axis=3, 
     xlab=NA, ylab=NA, xaxt='n', yaxt='n')
lines(DE.results.DEDD5$mean.discoveries$DESeq2.if.tmm, DE.results.DEDD5$mean.fdr$DESeq2.if.tmm, 
      lwd=2, col=col_vector[2], cex.axis=3)
lines(DE.results.DEDD5$mean.discoveries$voom.tmm, DE.results.DEDD5$mean.fdr$voom.tmm, 
      lwd=2, col=col_vector[3], cex.axis=3)
lines(DE.results.DEDD5$mean.discoveries$baySeq.tmm, DE.results.DEDD5$mean.fdr$baySeq.tmm, 
      lwd=2, col=col_vector[4], cex.axis=3)
lines(DE.results.DEDD5$mean.discoveries$mean.MDSeq.zi.tmm, DE.results.DEDD5$mean.fdr$mean.MDSeq.zi.tmm, 
      lwd=2, col=col_vector[5], cex.axis=3)
lines(DE.results.DEDD5$mean.discoveries$mean.lnHM.log.tmm, DE.results.DEDD5$mean.fdr$mean.lnHM.log.tmm, 
      lwd=2, col=col_vector[6], cex.axis=3)
abline(h=0.05, col='lightgrey')

plot(DE.results.DEDD10$mean.discoveries$edgeR.ql.tmm, DE.results.DEDD10$mean.fdr$edgeR.ql.tmm, 
     type='l', ylim=c(0,0.6), xlim=c(0,2000), lwd=2, col=col_vector[1], cex.axis=3, 
     xlab=NA, ylab=NA, xaxt='n', yaxt='n')
lines(DE.results.DEDD10$mean.discoveries$DESeq2.if.tmm, DE.results.DEDD10$mean.fdr$DESeq2.if.tmm, 
      lwd=2, col=col_vector[2], cex.axis=3)
lines(DE.results.DEDD10$mean.discoveries$voom.tmm, DE.results.DEDD10$mean.fdr$voom.tmm, 
      lwd=2, col=col_vector[3], cex.axis=3)
lines(DE.results.DEDD10$mean.discoveries$baySeq.tmm, DE.results.DEDD10$mean.fdr$baySeq.tmm, 
      lwd=2, col=col_vector[4], cex.axis=3)
lines(DE.results.DEDD10$mean.discoveries$mean.MDSeq.zi.tmm, DE.results.DEDD10$mean.fdr$mean.MDSeq.zi.tmm, 
      lwd=2, col=col_vector[5], cex.axis=3)
lines(DE.results.DEDD10$mean.discoveries$mean.lnHM.log.tmm, DE.results.DEDD10$mean.fdr$mean.lnHM.log.tmm, 
      lwd=2, col=col_vector[6], cex.axis=3)
abline(h=0.05, col='lightgrey')

plot(DE.results.DEDD20$mean.discoveries$edgeR.ql.tmm, DE.results.DEDD20$mean.fdr$edgeR.ql.tmm, 
     type='l', ylim=c(0,0.6), xlim=c(0,2000), lwd=2, col=col_vector[1], cex.axis=3, 
     xlab=NA, ylab=NA, xaxt='n', yaxt='n')
lines(DE.results.DEDD20$mean.discoveries$DESeq2.if.tmm, DE.results.DEDD20$mean.fdr$DESeq2.if.tmm, 
      lwd=2, col=col_vector[2], cex.axis=3)
lines(DE.results.DEDD20$mean.discoveries$voom.tmm, DE.results.DEDD20$mean.fdr$voom.tmm, 
      lwd=2, col=col_vector[3], cex.axis=3)
lines(DE.results.DEDD20$mean.discoveries$baySeq.tmm, DE.results.DEDD20$mean.fdr$baySeq.tmm, 
      lwd=2, col=col_vector[4], cex.axis=3)
lines(DE.results.DEDD20$mean.discoveries$mean.MDSeq.zi.tmm, DE.results.DEDD20$mean.fdr$mean.MDSeq.zi.tmm, 
      lwd=2, col=col_vector[5], cex.axis=3)
lines(DE.results.DEDD20$mean.discoveries$mean.lnHM.log.tmm, DE.results.DEDD20$mean.fdr$mean.lnHM.log.tmm, 
      lwd=2, col=col_vector[6], cex.axis=3)
abline(h=0.05, col='lightgrey')

plot(DE.results.DEDD50$mean.discoveries$edgeR.ql.tmm, DE.results.DEDD50$mean.fdr$edgeR.ql.tmm, 
     type='l', ylim=c(0,0.6), xlim=c(0,2000), lwd=2, col=col_vector[1], cex.axis=3, 
     xlab=NA, ylab=NA, xaxt='n', yaxt='n')
lines(DE.results.DEDD50$mean.discoveries$DESeq2.if.tmm, DE.results.DEDD50$mean.fdr$DESeq2.if.tmm, 
      lwd=2, col=col_vector[2], cex.axis=3)
lines(DE.results.DEDD50$mean.discoveries$voom.tmm, DE.results.DEDD50$mean.fdr$voom.tmm, 
      lwd=2, col=col_vector[3], cex.axis=3)
lines(DE.results.DEDD50$mean.discoveries$baySeq.tmm, DE.results.DEDD50$mean.fdr$baySeq.tmm, 
      lwd=2, col=col_vector[4], cex.axis=3)
lines(DE.results.DEDD50$mean.discoveries$mean.MDSeq.zi.tmm, DE.results.DEDD50$mean.fdr$mean.MDSeq.zi.tmm, 
      lwd=2, col=col_vector[5], cex.axis=3)
lines(DE.results.DEDD50$mean.discoveries$mean.lnHM.log.tmm, DE.results.DEDD50$mean.fdr$mean.lnHM.log.tmm, 
      lwd=2, col=col_vector[6], cex.axis=3)
abline(h=0.05, col='lightgrey')

legend("topleft", col=col_vector[1:6], bty='n', cex=2.5, ncol=1, lty=1, lwd=2, 
       legend=c("edgeR", "DESeq2", "voom", "baySeq", "MDSeq", "HM"))
