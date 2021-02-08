library(here)
library(RColorBrewer)

# Results initially saved by dataset, but want to assess by experiment type, 
# i.e. detecting differential expression, dispersion or distribution.

## Load data, create colour vector ####
folder <- "Results/GTEx blood simulated DE, DD, DEDD results Feb 2020"
for (i in c("DE", "DD", "DEDD")) {
  for (j in c("2_DEDD", "5_DEDD", "10_DEDD", "20_DEDD", "50_DEDD")) {
    assign(paste0(i, ".results.blood_", j), 
           readRDS(here(folder, paste0(i, ".results.blood_", j, ".rds"))))
  }
}
for (i in c("DE", "DEDD")) {
  for (j in c("2_DE", "5_DE", "10_DE", "20_DE", "50_DE")) {
    assign(paste0(i, ".results.blood_", j), 
           readRDS(here(folder, paste0(i, ".results.blood_", j, ".rds"))))
  }
}
for (i in c("DD", "DEDD")) {
  for (j in c("2_DD", "5_DD", "10_DD", "20_DD", "50_DD")) {
    assign(paste0(i, ".results.blood_", j), 
           readRDS(here(folder, paste0(i, ".results.blood_", j, ".rds"))))
  }
}
rm(i,j)

n <- 23
qual_col_pals = brewer.pal.info[brewer.pal.info$category == "qual" & 
                                  brewer.pal.info$colorblind == T,]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, 
                           rownames(qual_col_pals)))
# pie(rep(1,n), col=col_vector)
col_vector <- col_vector[-c(7,10,11,12,20)]


#################################
#### Differential dispersion ####
#################################

#############
#### AUC ####
par(mfrow=c(2,3), mar=c(1,2.5,5.5,1), mgp=c(3,0.7,0))
boxplot(DD.results.blood_2_DD$auc, names=NA, col=col_vector[1:7], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0.35,1.1), 
        main=paste0("AUC diff disp, 2 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in dispersion only"), 
        cex.main=1.7)
lines(c(7.5,7.5), c(0.35,0.95), col="lightgrey")
legend("topleft", fill=col_vector[1:7], bty='n', cex=1.5, ncol=3, 
       legend=c("diffVar", "MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", 
                "lnHM", "lnHM log"))
boxplot(DD.results.blood_5_DD$auc, names=NA, col=col_vector[1:7], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0.35,1.1), 
        main=paste0("AUC diff disp, 5 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in dispersion only"), 
        cex.main=1.7)
lines(c(7.5,7.5), c(0.35,0.95), col="lightgrey")
legend("topleft", fill=col_vector[1:7], bty='n', cex=1.5, ncol=3, 
       legend=c("diffVar", "MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", 
                "lnHM", "lnHM log"))
boxplot(DD.results.blood_10_DD$auc, names=NA, col=col_vector[1:7], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0.35,1.1), 
        main=paste0("AUC diff disp, 10 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in dispersion only"), 
        cex.main=1.7)
lines(c(7.5,7.5), c(0.35,0.95), col="lightgrey")
legend("topleft", fill=col_vector[1:7], bty='n', cex=1.5, ncol=3, 
       legend=c("diffVar", "MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", 
                "lnHM", "lnHM log"))
boxplot(DD.results.blood_20_DD$auc, names=NA, col=col_vector[1:7], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0.35,1.1), 
        main=paste0("AUC diff disp, 20 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in dispersion only"), 
        cex.main=1.7)
lines(c(7.5,7.5), c(0.35,0.95), col="lightgrey")
legend("topleft", fill=col_vector[1:7], bty='n', cex=1.5, ncol=3, 
       legend=c("diffVar", "MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", 
                "lnHM", "lnHM log"))
boxplot(DD.results.blood_50_DD$auc, names=NA, col=col_vector[1:7], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0.35,1.1), 
        main=paste0("AUC diff disp, 50 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in dispersion only"), 
        cex.main=1.7)
lines(c(7.5,7.5), c(0.35,0.95), col="lightgrey")
legend("topleft", fill=col_vector[1:7], bty='n', cex=1.5, ncol=3, 
       legend=c("diffVar", "MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", 
                "lnHM", "lnHM log"))
rbind(colMeans(DD.results.blood_2_DD$auc), 
      colMeans(DD.results.blood_5_DD$auc), 
      colMeans(DD.results.blood_10_DD$auc), 
      colMeans(DD.results.blood_20_DD$auc), 
      colMeans(DD.results.blood_50_DD$auc))

par(mfrow=c(2,3), mar=c(1,2.5,5.5,1), mgp=c(3,0.7,0))
boxplot(DD.results.blood_2_DEDD$auc, names=NA, col=col_vector[1:7], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0.35,1.1), 
        main=paste0("AUC diff disp, 2 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean and dispersion"), 
        cex.main=1.7)
lines(c(7.5,7.5), c(0.35,0.95), col="lightgrey")
legend("topleft", fill=col_vector[1:7], bty='n', cex=1.5, ncol=3, 
       legend=c("diffVar", "MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", 
                "lnHM", "lnHM log"))
boxplot(DD.results.blood_5_DEDD$auc, names=NA, col=col_vector[1:7], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0.35,1.1), 
        main=paste0("AUC diff disp, 5 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean and dispersion"), 
        cex.main=1.7)
lines(c(7.5,7.5), c(0.35,0.95), col="lightgrey")
legend("topleft", fill=col_vector[1:7], bty='n', cex=1.5, ncol=3, 
       legend=c("diffVar", "MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", 
                "lnHM", "lnHM log"))
boxplot(DD.results.blood_10_DEDD$auc, names=NA, col=col_vector[1:7], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0.35,1.1), 
        main=paste0("AUC diff disp, 10 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean and dispersion"), 
        cex.main=1.7)
lines(c(7.5,7.5), c(0.35,0.95), col="lightgrey")
legend("topleft", fill=col_vector[1:7], bty='n', cex=1.5, ncol=3, 
       legend=c("diffVar", "MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", 
                "lnHM", "lnHM log"))
boxplot(DD.results.blood_20_DEDD$auc, names=NA, col=col_vector[1:7], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0.35,1.1), 
        main=paste0("AUC diff disp, 20 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean and dispersion"), 
        cex.main=1.7)
lines(c(7.5,7.5), c(0.35,0.95), col="lightgrey")
legend("topleft", fill=col_vector[1:7], bty='n', cex=1.5, ncol=3, 
       legend=c("diffVar", "MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", 
                "lnHM", "lnHM log"))
boxplot(DD.results.blood_50_DEDD$auc, names=NA, col=col_vector[1:7], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0.35,1.1), 
        main=paste0("AUC diff disp, 50 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean and dispersion"), 
        cex.main=1.7)
lines(c(7.5,7.5), c(0.35,0.95), col="lightgrey")
legend("topleft", fill=col_vector[1:7], bty='n', cex=1.5, ncol=3, 
       legend=c("diffVar", "MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", 
                "lnHM", "lnHM log"))
rbind(colMeans(DD.results.blood_2_DEDD$auc), 
      colMeans(DD.results.blood_5_DEDD$auc), 
      colMeans(DD.results.blood_10_DEDD$auc), 
      colMeans(DD.results.blood_20_DEDD$auc), 
      colMeans(DD.results.blood_50_DEDD$auc))


##############
#### pAUC ####
par(mfrow=c(2,3), mar=c(1,2.5,5.5,1), mgp=c(3,0.7,0))
boxplot(DD.results.blood_2_DD$pauc, names=NA, col=col_vector[1:7], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,0.04), 
        main=paste0("Partial AUC diff disp, 2 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in dispersion only"), 
        cex.main=1.7)
lines(c(7.5,7.5), c(0,0.033), col="lightgrey")
legend("topleft", fill=col_vector[1:7], bty='n', cex=1.5, ncol=3, 
       legend=c("diffVar", "MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", 
                "lnHM", "lnHM log"))
boxplot(DD.results.blood_5_DD$pauc, names=NA, col=col_vector[1:7], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,0.04), 
        main=paste0("Partial AUC diff disp, 5 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in dispersion only"), 
        cex.main=1.7)
lines(c(7.5,7.5), c(0,0.033), col="lightgrey")
legend("topleft", fill=col_vector[1:7], bty='n', cex=1.5, ncol=3, 
       legend=c("diffVar", "MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", 
                "lnHM", "lnHM log"))
boxplot(DD.results.blood_10_DD$pauc, names=NA, col=col_vector[1:7], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,0.04), 
        main=paste0("Partial AUC diff disp, 10 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in dispersion only"), 
        cex.main=1.7)
lines(c(7.5,7.5), c(0,0.033), col="lightgrey")
legend("topleft", fill=col_vector[1:7], bty='n', cex=1.5, ncol=3, 
       legend=c("diffVar", "MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", 
                "lnHM", "lnHM log"))
boxplot(DD.results.blood_20_DD$pauc, names=NA, col=col_vector[1:7], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,0.04), 
        main=paste0("Partial AUC diff disp, 20 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in dispersion only"), 
        cex.main=1.7)
lines(c(7.5,7.5), c(0,0.033), col="lightgrey")
legend("topleft", fill=col_vector[1:7], bty='n', cex=1.5, ncol=3, 
       legend=c("diffVar", "MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", 
                "lnHM", "lnHM log"))
boxplot(DD.results.blood_50_DD$pauc, names=NA, col=col_vector[1:7], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,0.04), 
        main=paste0("Partial AUC diff disp, 50 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in dispersion only"), 
        cex.main=1.7)
lines(c(7.5,7.5), c(0,0.033), col="lightgrey")
legend("topleft", fill=col_vector[1:7], bty='n', cex=1.5, ncol=3, 
       legend=c("diffVar", "MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", 
                "lnHM", "lnHM log"))
rbind(colMeans(DD.results.blood_2_DD$pauc), 
      colMeans(DD.results.blood_5_DD$pauc), 
      colMeans(DD.results.blood_10_DD$pauc), 
      colMeans(DD.results.blood_20_DD$pauc), 
      colMeans(DD.results.blood_50_DD$pauc))

par(mfrow=c(2,3), mar=c(1,2.5,5.5,1), mgp=c(3,0.7,0))
boxplot(DD.results.blood_2_DEDD$pauc, names=NA, col=col_vector[1:7], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,0.04), 
        main=paste0("Partial AUC diff disp, 2 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean and dispersion"), 
        cex.main=1.7)
lines(c(7.5,7.5), c(0,0.033), col="lightgrey")
legend("topleft", fill=col_vector[1:7], bty='n', cex=1.5, ncol=3, 
       legend=c("diffVar", "MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", 
                "lnHM", "lnHM log"))
boxplot(DD.results.blood_5_DEDD$pauc, names=NA, col=col_vector[1:7], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,0.04), 
        main=paste0("Partial AUC diff disp, 5 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean and dispersion"), 
        cex.main=1.7)
lines(c(7.5,7.5), c(0,0.033), col="lightgrey")
legend("topleft", fill=col_vector[1:7], bty='n', cex=1.5, ncol=3, 
       legend=c("diffVar", "MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", 
                "lnHM", "lnHM log"))
boxplot(DD.results.blood_10_DEDD$pauc, names=NA, col=col_vector[1:7], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,0.04), 
        main=paste0("Partial AUC diff disp, 10 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean and dispersion"), 
        cex.main=1.7)
lines(c(7.5,7.5), c(0,0.033), col="lightgrey")
legend("topleft", fill=col_vector[1:7], bty='n', cex=1.5, ncol=3, 
       legend=c("diffVar", "MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", 
                "lnHM", "lnHM log"))
boxplot(DD.results.blood_20_DEDD$pauc, names=NA, col=col_vector[1:7], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,0.04), 
        main=paste0("Partial AUC diff disp, 20 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean and dispersion"), 
        cex.main=1.7)
lines(c(7.5,7.5), c(0,0.033), col="lightgrey")
legend("topleft", fill=col_vector[1:7], bty='n', cex=1.5, ncol=3, 
       legend=c("diffVar", "MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", 
                "lnHM", "lnHM log"))
boxplot(DD.results.blood_50_DEDD$pauc, names=NA, col=col_vector[1:7], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,0.04), 
        main=paste0("Partial AUC diff disp, 50 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean and dispersion"), 
        cex.main=1.7)
lines(c(7.5,7.5), c(0,0.033), col="lightgrey")
legend("topleft", fill=col_vector[1:7], bty='n', cex=1.5, ncol=3, 
       legend=c("diffVar", "MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", 
                "lnHM", "lnHM log"))
rbind(colMeans(DD.results.blood_2_DEDD$pauc), 
      colMeans(DD.results.blood_5_DEDD$pauc), 
      colMeans(DD.results.blood_10_DEDD$pauc), 
      colMeans(DD.results.blood_20_DEDD$pauc), 
      colMeans(DD.results.blood_50_DEDD$pauc))


#############
#### FDR ####
par(mfrow=c(2,3), mar=c(1,2.5,5.5,1), mgp=c(3,0.7,0))
boxplot(DD.results.blood_2_DD$fdr, names=NA, col=col_vector[1:7], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.2), 
        main=paste0("FDR diff disp, 2 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in dispersion only"), 
        cex.main=1.7)
lines(c(7.5,7.5), c(0,1), col="lightgrey"); abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[1:7], bty='n', cex=1.5, ncol=3, 
       legend=c("diffVar", "MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", 
                "lnHM", "lnHM log"))
boxplot(DD.results.blood_5_DD$fdr, names=NA, col=col_vector[1:7], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.2), 
        main=paste0("FDR diff disp, 5 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in dispersion only"), 
        cex.main=1.7)
lines(c(7.5,7.5), c(0,1), col="lightgrey"); abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[1:7], bty='n', cex=1.5, ncol=3, 
       legend=c("diffVar", "MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", 
                "lnHM", "lnHM log"))
boxplot(DD.results.blood_10_DD$fdr, names=NA, col=col_vector[1:7], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.2), 
        main=paste0("FDR diff disp, 10 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in dispersion only"), 
        cex.main=1.7)
lines(c(7.5,7.5), c(0,1), col="lightgrey"); abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[1:7], bty='n', cex=1.5, ncol=3, 
       legend=c("diffVar", "MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", 
                "lnHM", "lnHM log"))
boxplot(DD.results.blood_20_DD$fdr, names=NA, col=col_vector[1:7], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.2), 
        main=paste0("FDR diff disp, 20 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in dispersion only"), 
        cex.main=1.7)
lines(c(7.5,7.5), c(0,1), col="lightgrey"); abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[1:7], bty='n', cex=1.5, ncol=3, 
       legend=c("diffVar", "MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", 
                "lnHM", "lnHM log"))
boxplot(DD.results.blood_50_DD$fdr, names=NA, col=col_vector[1:7], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.2), 
        main=paste0("FDR diff disp, 50 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in dispersion only"), 
        cex.main=1.7)
lines(c(7.5,7.5), c(0,1), col="lightgrey"); abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[1:7], bty='n', cex=1.5, ncol=3, 
       legend=c("diffVar", "MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", 
                "lnHM", "lnHM log"))
rbind(colMeans(DD.results.blood_2_DD$fdr), 
      colMeans(DD.results.blood_5_DD$fdr), 
      colMeans(DD.results.blood_10_DD$fdr), 
      colMeans(DD.results.blood_20_DD$fdr), 
      colMeans(DD.results.blood_50_DD$fdr))

par(mfrow=c(2,3), mar=c(1,2.5,5.5,1), mgp=c(3,0.7,0))
boxplot(DD.results.blood_2_DEDD$fdr, names=NA, col=col_vector[1:7], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.2), 
        main=paste0("FDR diff disp, 2 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean and dispersion"), 
        cex.main=1.7)
lines(c(7.5,7.5), c(0,1), col="lightgrey"); abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[1:7], bty='n', cex=1.5, ncol=3, 
       legend=c("diffVar", "MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", 
                "lnHM", "lnHM log"))
boxplot(DD.results.blood_5_DEDD$fdr, names=NA, col=col_vector[1:7], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.2), 
        main=paste0("FDR diff disp, 5 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean and dispersion"), 
        cex.main=1.7)
lines(c(7.5,7.5), c(0,1), col="lightgrey"); abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[1:7], bty='n', cex=1.5, ncol=3, 
       legend=c("diffVar", "MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", 
                "lnHM", "lnHM log"))
boxplot(DD.results.blood_10_DEDD$fdr, names=NA, col=col_vector[1:7], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.2), 
        main=paste0("FDR diff disp, 10 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean and dispersion"), 
        cex.main=1.7)
lines(c(7.5,7.5), c(0,1), col="lightgrey"); abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[1:7], bty='n', cex=1.5, ncol=3, 
       legend=c("diffVar", "MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", 
                "lnHM", "lnHM log"))
boxplot(DD.results.blood_20_DEDD$fdr, names=NA, col=col_vector[1:7], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.2), 
        main=paste0("FDR diff disp, 20 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean and dispersion"), 
        cex.main=1.7)
lines(c(7.5,7.5), c(0,1), col="lightgrey"); abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[1:7], bty='n', cex=1.5, ncol=3, 
       legend=c("diffVar", "MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", 
                "lnHM", "lnHM log"))
boxplot(DD.results.blood_50_DEDD$fdr, names=NA, col=col_vector[1:7], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.2), 
        main=paste0("FDR diff disp, 50 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean and dispersion"), 
        cex.main=1.7)
lines(c(7.5,7.5), c(0,1), col="lightgrey"); abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[1:7], bty='n', cex=1.5, ncol=3, 
       legend=c("diffVar", "MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", 
                "lnHM", "lnHM log"))
rbind(colMeans(DD.results.blood_2_DEDD$fdr), 
      colMeans(DD.results.blood_5_DEDD$fdr), 
      colMeans(DD.results.blood_10_DEDD$fdr), 
      colMeans(DD.results.blood_20_DEDD$fdr), 
      colMeans(DD.results.blood_50_DEDD$fdr))


#############
#### TPR ####
par(mfrow=c(2,3), mar=c(1,2.5,5.5,1), mgp=c(3,0.7,0))
boxplot(DD.results.blood_2_DD$tpr, names=NA, col=col_vector[1:7], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,0.8), 
        main=paste0("TPR diff disp, 2 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in dispersion only"), 
        cex.main=1.7)
lines(c(7.5,7.5), c(0,0.7), col="lightgrey")
legend("topleft", fill=col_vector[1:7], bty='n', cex=1.5, ncol=3, 
       legend=c("diffVar", "MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", 
                "lnHM", "lnHM log"))
boxplot(DD.results.blood_5_DD$tpr, names=NA, col=col_vector[1:7], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,0.8), 
        main=paste0("TPR diff disp, 5 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in dispersion only"), 
        cex.main=1.7)
lines(c(7.5,7.5), c(0,0.7), col="lightgrey")
legend("topleft", fill=col_vector[1:7], bty='n', cex=1.5, ncol=3, 
       legend=c("diffVar", "MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", 
                "lnHM", "lnHM log"))
boxplot(DD.results.blood_10_DD$tpr, names=NA, col=col_vector[1:7], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,0.8), 
        main=paste0("TPR diff disp, 10 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in dispersion only"), 
        cex.main=1.7)
lines(c(7.5,7.5), c(0,0.7), col="lightgrey")
legend("topleft", fill=col_vector[1:7], bty='n', cex=1.5, ncol=3, 
       legend=c("diffVar", "MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", 
                "lnHM", "lnHM log"))
boxplot(DD.results.blood_20_DD$tpr, names=NA, col=col_vector[1:7], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,0.8), 
        main=paste0("TPR diff disp, 20 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in dispersion only"), 
        cex.main=1.7)
lines(c(7.5,7.5), c(0,0.7), col="lightgrey")
legend("topleft", fill=col_vector[1:7], bty='n', cex=1.5, ncol=3, 
       legend=c("diffVar", "MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", 
                "lnHM", "lnHM log"))
boxplot(DD.results.blood_50_DD$tpr, names=NA, col=col_vector[1:7], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,0.8), 
        main=paste0("TPR diff disp, 50 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in dispersion only"), 
        cex.main=1.7)
lines(c(7.5,7.5), c(0,0.7), col="lightgrey")
legend("topleft", fill=col_vector[1:7], bty='n', cex=1.5, ncol=3, 
       legend=c("diffVar", "MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", 
                "lnHM", "lnHM log"))
rbind(colMeans(DD.results.blood_2_DD$tpr), 
      colMeans(DD.results.blood_5_DD$tpr), 
      colMeans(DD.results.blood_10_DD$tpr), 
      colMeans(DD.results.blood_20_DD$tpr), 
      colMeans(DD.results.blood_50_DD$tpr))

par(mfrow=c(2,3), mar=c(1,2.5,5.5,1), mgp=c(3,0.7,0))
boxplot(DD.results.blood_2_DEDD$tpr, names=NA, col=col_vector[1:7], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,0.8), 
        main=paste0("TPR diff disp, 2 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean and dispersion"), 
        cex.main=1.7)
lines(c(7.5,7.5), c(0,0.7), col="lightgrey")
legend("topleft", fill=col_vector[1:7], bty='n', cex=1.5, ncol=3, 
       legend=c("diffVar", "MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", 
                "lnHM", "lnHM log"))
boxplot(DD.results.blood_5_DEDD$tpr, names=NA, col=col_vector[1:7], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,0.8), 
        main=paste0("TPR diff disp, 5 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean and dispersion"), 
        cex.main=1.7)
lines(c(7.5,7.5), c(0,0.7), col="lightgrey")
legend("topleft", fill=col_vector[1:7], bty='n', cex=1.5, ncol=3, 
       legend=c("diffVar", "MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", 
                "lnHM", "lnHM log"))
boxplot(DD.results.blood_10_DEDD$tpr, names=NA, col=col_vector[1:7], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,0.8), 
        main=paste0("TPR diff disp, 10 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean and dispersion"), 
        cex.main=1.7)
lines(c(7.5,7.5), c(0,0.7), col="lightgrey")
legend("topleft", fill=col_vector[1:7], bty='n', cex=1.5, ncol=3, 
       legend=c("diffVar", "MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", 
                "lnHM", "lnHM log"))
boxplot(DD.results.blood_20_DEDD$tpr, names=NA, col=col_vector[1:7], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,0.8), 
        main=paste0("TPR diff disp, 20 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean and dispersion"), 
        cex.main=1.7)
lines(c(7.5,7.5), c(0,0.7), col="lightgrey")
legend("topleft", fill=col_vector[1:7], bty='n', cex=1.5, ncol=3, 
       legend=c("diffVar", "MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", 
                "lnHM", "lnHM log"))
boxplot(DD.results.blood_50_DEDD$tpr, names=NA, col=col_vector[1:7], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,0.8), 
        main=paste0("TPR diff disp, 50 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean and dispersion"), 
        cex.main=1.7)
lines(c(7.5,7.5), c(0,0.7), col="lightgrey")
legend("topleft", fill=col_vector[1:7], bty='n', cex=1.5, ncol=3, 
       legend=c("diffVar", "MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", 
                "lnHM", "lnHM log"))
rbind(colMeans(DD.results.blood_2_DEDD$tpr), 
      colMeans(DD.results.blood_5_DEDD$tpr), 
      colMeans(DD.results.blood_10_DEDD$tpr), 
      colMeans(DD.results.blood_20_DEDD$tpr), 
      colMeans(DD.results.blood_50_DEDD$tpr))


###############################
#### False discovery plots ####
par(mfrow=c(4,5), mar=c(2,1.5,4,0.5), mgp=c(3,0.5,0))

plot(DD.results.blood_2_DD$mean.discoveries$diffVar.tmm, 
     DD.results.blood_2_DD$mean.fdr$diffVar.tmm, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("2 samples per group\n", "TMM normalisation\n", 
                 "Differences in dispersion only"))
lines(DD.results.blood_2_DD$mean.discoveries$disp.MDSeq.zi.tmm, 
      DD.results.blood_2_DD$mean.fdr$disp.MDSeq.zi.tmm, 
      type='l', col=col_vector[2])
lines(DD.results.blood_2_DD$mean.discoveries$disp.MDSeq.nozi.tmm, 
      DD.results.blood_2_DD$mean.fdr$disp.MDSeq.nozi.tmm, 
      type='l', col=col_vector[3])
lines(DD.results.blood_2_DD$mean.discoveries$disp.expHM.untr.tmm, 
      DD.results.blood_2_DD$mean.fdr$disp.expHM.untr.tmm, 
      type='l', col=col_vector[4])
lines(DD.results.blood_2_DD$mean.discoveries$disp.expHM.log.tmm, 
      DD.results.blood_2_DD$mean.fdr$disp.expHM.log.tmm, 
      type='l', col=col_vector[5])
lines(DD.results.blood_2_DD$mean.discoveries$disp.lnHM.untr.tmm, 
      DD.results.blood_2_DD$mean.fdr$disp.lnHM.untr.tmm, 
      type='l', col=col_vector[6])
lines(DD.results.blood_2_DD$mean.discoveries$disp.lnHM.log.tmm, 
      DD.results.blood_2_DD$mean.fdr$disp.lnHM.log.tmm, 
      type='l', col=col_vector[7])
lines(c(0,1000), c(0.05,0.05), col="lightgrey")
legend("bottomright", bty='n', col=col_vector[1:7], lty=1, ncol=1, lwd=2, 
       legend=c("diffVar", "MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", 
                "lnHM", "lnHM log"))
plot(DD.results.blood_5_DD$mean.discoveries$diffVar.tmm, 
     DD.results.blood_5_DD$mean.fdr$diffVar.tmm, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("5 samples per group\n", "TMM normalisation\n", 
                 "Differences in dispersion only"))
lines(DD.results.blood_5_DD$mean.discoveries$disp.MDSeq.zi.tmm, 
      DD.results.blood_5_DD$mean.fdr$disp.MDSeq.zi.tmm, 
      type='l', col=col_vector[2])
lines(DD.results.blood_5_DD$mean.discoveries$disp.MDSeq.nozi.tmm, 
      DD.results.blood_5_DD$mean.fdr$disp.MDSeq.nozi.tmm, 
      type='l', col=col_vector[3])
lines(DD.results.blood_5_DD$mean.discoveries$disp.expHM.untr.tmm, 
      DD.results.blood_5_DD$mean.fdr$disp.expHM.untr.tmm, 
      type='l', col=col_vector[4])
lines(DD.results.blood_5_DD$mean.discoveries$disp.expHM.log.tmm, 
      DD.results.blood_5_DD$mean.fdr$disp.expHM.log.tmm, 
      type='l', col=col_vector[5])
lines(DD.results.blood_5_DD$mean.discoveries$disp.lnHM.untr.tmm, 
      DD.results.blood_5_DD$mean.fdr$disp.lnHM.untr.tmm, 
      type='l', col=col_vector[6])
lines(DD.results.blood_5_DD$mean.discoveries$disp.lnHM.log.tmm, 
      DD.results.blood_5_DD$mean.fdr$disp.lnHM.log.tmm, 
      type='l', col=col_vector[7])
abline(h=0.05, col="lightgrey")
plot(DD.results.blood_10_DD$mean.discoveries$diffVar.tmm, 
     DD.results.blood_10_DD$mean.fdr$diffVar.tmm, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("10 samples per group\n", "TMM normalisation\n", 
                 "Differences in dispersion only"))
lines(DD.results.blood_10_DD$mean.discoveries$disp.MDSeq.zi.tmm, 
      DD.results.blood_10_DD$mean.fdr$disp.MDSeq.zi.tmm, 
      type='l', col=col_vector[2])
lines(DD.results.blood_10_DD$mean.discoveries$disp.MDSeq.nozi.tmm, 
      DD.results.blood_10_DD$mean.fdr$disp.MDSeq.nozi.tmm, 
      type='l', col=col_vector[3])
lines(DD.results.blood_10_DD$mean.discoveries$disp.expHM.untr.tmm, 
      DD.results.blood_10_DD$mean.fdr$disp.expHM.untr.tmm, 
      type='l', col=col_vector[4])
lines(DD.results.blood_10_DD$mean.discoveries$disp.expHM.log.tmm, 
      DD.results.blood_10_DD$mean.fdr$disp.expHM.log.tmm, 
      type='l', col=col_vector[5])
lines(DD.results.blood_10_DD$mean.discoveries$disp.lnHM.untr.tmm, 
      DD.results.blood_10_DD$mean.fdr$disp.lnHM.untr.tmm, 
      type='l', col=col_vector[6])
lines(DD.results.blood_10_DD$mean.discoveries$disp.lnHM.log.tmm, 
      DD.results.blood_10_DD$mean.fdr$disp.lnHM.log.tmm, 
      type='l', col=col_vector[7])
abline(h=0.05, col="lightgrey")
plot(DD.results.blood_20_DD$mean.discoveries$diffVar.tmm, 
     DD.results.blood_20_DD$mean.fdr$diffVar.tmm, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("20 samples per group\n", "TMM normalisation\n", 
                 "Differences in dispersion only"))
lines(DD.results.blood_20_DD$mean.discoveries$disp.MDSeq.zi.tmm, 
      DD.results.blood_20_DD$mean.fdr$disp.MDSeq.zi.tmm, 
      type='l', col=col_vector[2])
lines(DD.results.blood_20_DD$mean.discoveries$disp.MDSeq.nozi.tmm, 
      DD.results.blood_20_DD$mean.fdr$disp.MDSeq.nozi.tmm, 
      type='l', col=col_vector[3])
lines(DD.results.blood_20_DD$mean.discoveries$disp.expHM.untr.tmm, 
      DD.results.blood_20_DD$mean.fdr$disp.expHM.untr.tmm, 
      type='l', col=col_vector[4])
lines(DD.results.blood_20_DD$mean.discoveries$disp.expHM.log.tmm, 
      DD.results.blood_20_DD$mean.fdr$disp.expHM.log.tmm, 
      type='l', col=col_vector[5])
lines(DD.results.blood_20_DD$mean.discoveries$disp.lnHM.untr.tmm, 
      DD.results.blood_20_DD$mean.fdr$disp.lnHM.untr.tmm, 
      type='l', col=col_vector[6])
lines(DD.results.blood_20_DD$mean.discoveries$disp.lnHM.log.tmm, 
      DD.results.blood_20_DD$mean.fdr$disp.lnHM.log.tmm, 
      type='l', col=col_vector[7])
abline(h=0.05, col="lightgrey")
plot(DD.results.blood_50_DD$mean.discoveries$diffVar.tmm, 
     DD.results.blood_50_DD$mean.fdr$diffVar.tmm, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("50 samples per group\n", "TMM normalisation\n", 
                 "Differences in dispersion only"))
lines(DD.results.blood_50_DD$mean.discoveries$disp.MDSeq.zi.tmm, 
      DD.results.blood_50_DD$mean.fdr$disp.MDSeq.zi.tmm, 
      type='l', col=col_vector[2])
lines(DD.results.blood_50_DD$mean.discoveries$disp.MDSeq.nozi.tmm, 
      DD.results.blood_50_DD$mean.fdr$disp.MDSeq.nozi.tmm, 
      type='l', col=col_vector[3])
lines(DD.results.blood_50_DD$mean.discoveries$disp.expHM.untr.tmm, 
      DD.results.blood_50_DD$mean.fdr$disp.expHM.untr.tmm, 
      type='l', col=col_vector[4])
lines(DD.results.blood_50_DD$mean.discoveries$disp.expHM.log.tmm, 
      DD.results.blood_50_DD$mean.fdr$disp.expHM.log.tmm, 
      type='l', col=col_vector[5])
lines(DD.results.blood_50_DD$mean.discoveries$disp.lnHM.untr.tmm, 
      DD.results.blood_50_DD$mean.fdr$disp.lnHM.untr.tmm, 
      type='l', col=col_vector[6])
lines(DD.results.blood_50_DD$mean.discoveries$disp.lnHM.log.tmm, 
      DD.results.blood_50_DD$mean.fdr$disp.lnHM.log.tmm, 
      type='l', col=col_vector[7])
abline(h=0.05, col="lightgrey")

plot(DD.results.blood_2_DD$mean.discoveries$diffVar.rle, 
     DD.results.blood_2_DD$mean.fdr$diffVar.rle, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("2 samples per group\n", "RLE normalisation\n", 
                 "Differences in dispersion only"))
lines(DD.results.blood_2_DD$mean.discoveries$disp.MDSeq.zi.rle, 
      DD.results.blood_2_DD$mean.fdr$disp.MDSeq.zi.rle, 
      type='l', col=col_vector[2])
lines(DD.results.blood_2_DD$mean.discoveries$disp.MDSeq.nozi.rle, 
      DD.results.blood_2_DD$mean.fdr$disp.MDSeq.nozi.rle, 
      type='l', col=col_vector[3])
lines(DD.results.blood_2_DD$mean.discoveries$disp.expHM.untr.rle, 
      DD.results.blood_2_DD$mean.fdr$disp.expHM.untr.rle, 
      type='l', col=col_vector[4])
lines(DD.results.blood_2_DD$mean.discoveries$disp.expHM.log.rle, 
      DD.results.blood_2_DD$mean.fdr$disp.expHM.log.rle, 
      type='l', col=col_vector[5])
lines(DD.results.blood_2_DD$mean.discoveries$disp.lnHM.untr.rle, 
      DD.results.blood_2_DD$mean.fdr$disp.lnHM.untr.rle, 
      type='l', col=col_vector[6])
lines(DD.results.blood_2_DD$mean.discoveries$disp.lnHM.log.rle, 
      DD.results.blood_2_DD$mean.fdr$disp.lnHM.log.rle, 
      type='l', col=col_vector[7])
lines(c(0,1000), c(0.05,0.05), col="lightgrey")
legend("bottomright", bty='n', col=col_vector[1:7], lty=1, ncol=1, lwd=2, 
       legend=c("diffVar", "MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", 
                "lnHM", "lnHM log"))
plot(DD.results.blood_5_DD$mean.discoveries$diffVar.rle, 
     DD.results.blood_5_DD$mean.fdr$diffVar.rle, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("5 samples per group\n", "RLE normalisation\n", 
                 "Differences in dispersion only"))
lines(DD.results.blood_5_DD$mean.discoveries$disp.MDSeq.zi.rle, 
      DD.results.blood_5_DD$mean.fdr$disp.MDSeq.zi.rle, 
      type='l', col=col_vector[2])
lines(DD.results.blood_5_DD$mean.discoveries$disp.MDSeq.nozi.rle, 
      DD.results.blood_5_DD$mean.fdr$disp.MDSeq.nozi.rle, 
      type='l', col=col_vector[3])
lines(DD.results.blood_5_DD$mean.discoveries$disp.expHM.untr.rle, 
      DD.results.blood_5_DD$mean.fdr$disp.expHM.untr.rle, 
      type='l', col=col_vector[4])
lines(DD.results.blood_5_DD$mean.discoveries$disp.expHM.log.rle, 
      DD.results.blood_5_DD$mean.fdr$disp.expHM.log.rle, 
      type='l', col=col_vector[5])
lines(DD.results.blood_5_DD$mean.discoveries$disp.lnHM.untr.rle, 
      DD.results.blood_5_DD$mean.fdr$disp.lnHM.untr.rle, 
      type='l', col=col_vector[6])
lines(DD.results.blood_5_DD$mean.discoveries$disp.lnHM.log.rle, 
      DD.results.blood_5_DD$mean.fdr$disp.lnHM.log.rle, 
      type='l', col=col_vector[7])
abline(h=0.05, col='lightgrey')
plot(DD.results.blood_10_DD$mean.discoveries$diffVar.rle, 
     DD.results.blood_10_DD$mean.fdr$diffVar.rle, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("10 samples per group\n", "RLE normalisation\n", 
                 "Differences in dispersion only"))
lines(DD.results.blood_10_DD$mean.discoveries$disp.MDSeq.zi.rle, 
      DD.results.blood_10_DD$mean.fdr$disp.MDSeq.zi.rle, 
      type='l', col=col_vector[2])
lines(DD.results.blood_10_DD$mean.discoveries$disp.MDSeq.nozi.rle, 
      DD.results.blood_10_DD$mean.fdr$disp.MDSeq.nozi.rle, 
      type='l', col=col_vector[3])
lines(DD.results.blood_10_DD$mean.discoveries$disp.expHM.untr.rle, 
      DD.results.blood_10_DD$mean.fdr$disp.expHM.untr.rle, 
      type='l', col=col_vector[4])
lines(DD.results.blood_10_DD$mean.discoveries$disp.expHM.log.rle, 
      DD.results.blood_10_DD$mean.fdr$disp.expHM.log.rle, 
      type='l', col=col_vector[5])
lines(DD.results.blood_10_DD$mean.discoveries$disp.lnHM.untr.rle, 
      DD.results.blood_10_DD$mean.fdr$disp.lnHM.untr.rle, 
      type='l', col=col_vector[6])
lines(DD.results.blood_10_DD$mean.discoveries$disp.lnHM.log.rle, 
      DD.results.blood_10_DD$mean.fdr$disp.lnHM.log.rle, 
      type='l', col=col_vector[7])
abline(h=0.05, col='lightgrey')
plot(DD.results.blood_20_DD$mean.discoveries$diffVar.rle, 
     DD.results.blood_20_DD$mean.fdr$diffVar.rle, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("20 samples per group\n", "RLE normalisation\n", 
                 "Differences in dispersion only"))
lines(DD.results.blood_20_DD$mean.discoveries$disp.MDSeq.zi.rle, 
      DD.results.blood_20_DD$mean.fdr$disp.MDSeq.zi.rle, 
      type='l', col=col_vector[2])
lines(DD.results.blood_20_DD$mean.discoveries$disp.MDSeq.nozi.rle, 
      DD.results.blood_20_DD$mean.fdr$disp.MDSeq.nozi.rle, 
      type='l', col=col_vector[3])
lines(DD.results.blood_20_DD$mean.discoveries$disp.expHM.untr.rle, 
      DD.results.blood_20_DD$mean.fdr$disp.expHM.untr.rle, 
      type='l', col=col_vector[4])
lines(DD.results.blood_20_DD$mean.discoveries$disp.expHM.log.rle, 
      DD.results.blood_20_DD$mean.fdr$disp.expHM.log.rle, 
      type='l', col=col_vector[5])
lines(DD.results.blood_20_DD$mean.discoveries$disp.lnHM.untr.rle, 
      DD.results.blood_20_DD$mean.fdr$disp.lnHM.untr.rle, 
      type='l', col=col_vector[6])
lines(DD.results.blood_20_DD$mean.discoveries$disp.lnHM.log.rle, 
      DD.results.blood_20_DD$mean.fdr$disp.lnHM.log.rle, 
      type='l', col=col_vector[7])
abline(h=0.05, col='lightgrey')
plot(DD.results.blood_50_DD$mean.discoveries$diffVar.rle, 
     DD.results.blood_50_DD$mean.fdr$diffVar.rle, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("50 samples per group\n", "RLE normalisation\n", 
                 "Differences in dispersion only"))
lines(DD.results.blood_50_DD$mean.discoveries$disp.MDSeq.nozi.rle, 
      DD.results.blood_50_DD$mean.fdr$disp.MDSeq.nozi.rle, 
      type='l', col=col_vector[2])
lines(DD.results.blood_50_DD$mean.discoveries$disp.MDSeq.zi.rle, 
      DD.results.blood_50_DD$mean.fdr$disp.MDSeq.zi.rle, 
      type='l', col=col_vector[3])
lines(DD.results.blood_50_DD$mean.discoveries$disp.expHM.untr.rle, 
      DD.results.blood_50_DD$mean.fdr$disp.expHM.untr.rle, 
      type='l', col=col_vector[4])
lines(DD.results.blood_50_DD$mean.discoveries$disp.expHM.log.rle, 
      DD.results.blood_50_DD$mean.fdr$disp.expHM.log.rle, 
      type='l', col=col_vector[5])
lines(DD.results.blood_50_DD$mean.discoveries$disp.lnHM.untr.rle, 
      DD.results.blood_50_DD$mean.fdr$disp.lnHM.untr.rle, 
      type='l', col=col_vector[6])
lines(DD.results.blood_50_DD$mean.discoveries$disp.lnHM.log.rle, 
      DD.results.blood_50_DD$mean.fdr$disp.lnHM.log.rle, 
      type='l', col=col_vector[7])
abline(h=0.05, col='lightgrey')

plot(DD.results.blood_2_DEDD$mean.discoveries$diffVar.tmm, 
     DD.results.blood_2_DEDD$mean.fdr$diffVar.tmm, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("2 samples per group\n", "TMM normalisation\n", 
                 "Differences in mean and dispersion"))
lines(DD.results.blood_2_DEDD$mean.discoveries$disp.MDSeq.zi.tmm, 
      DD.results.blood_2_DEDD$mean.fdr$disp.MDSeq.zi.tmm, 
      type='l', col=col_vector[2])
lines(DD.results.blood_2_DEDD$mean.discoveries$disp.MDSeq.nozi.tmm, 
      DD.results.blood_2_DEDD$mean.fdr$disp.MDSeq.nozi.tmm, 
      type='l', col=col_vector[3])
lines(DD.results.blood_2_DEDD$mean.discoveries$disp.expHM.untr.tmm, 
      DD.results.blood_2_DEDD$mean.fdr$disp.expHM.untr.tmm, 
      type='l', col=col_vector[4])
lines(DD.results.blood_2_DEDD$mean.discoveries$disp.expHM.log.tmm, 
      DD.results.blood_2_DEDD$mean.fdr$disp.expHM.log.tmm, 
      type='l', col=col_vector[5])
lines(DD.results.blood_2_DEDD$mean.discoveries$disp.lnHM.untr.tmm, 
      DD.results.blood_2_DEDD$mean.fdr$disp.lnHM.untr.tmm, 
      type='l', col=col_vector[6])
lines(DD.results.blood_2_DEDD$mean.discoveries$disp.lnHM.log.tmm, 
      DD.results.blood_2_DEDD$mean.fdr$disp.lnHM.log.tmm, 
      type='l', col=col_vector[7])
lines(c(0,1000), c(0.05,0.05), col="lightgrey")
legend("bottomright", bty='n', col=col_vector[1:7], lty=1, ncol=1, lwd=2, 
       legend=c("diffVar", "MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", 
                "lnHM", "lnHM log"))
plot(DD.results.blood_5_DEDD$mean.discoveries$diffVar.tmm, 
     DD.results.blood_5_DEDD$mean.fdr$diffVar.tmm, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("5 samples per group\n", "TMM normalisation\n", 
                 "Differences in mean and dispersion"))
lines(DD.results.blood_5_DEDD$mean.discoveries$disp.MDSeq.zi.tmm, 
      DD.results.blood_5_DEDD$mean.fdr$disp.MDSeq.zi.tmm, 
      type='l', col=col_vector[2])
lines(DD.results.blood_5_DEDD$mean.discoveries$disp.MDSeq.nozi.tmm, 
      DD.results.blood_5_DEDD$mean.fdr$disp.MDSeq.nozi.tmm, 
      type='l', col=col_vector[3])
lines(DD.results.blood_5_DEDD$mean.discoveries$disp.expHM.untr.tmm, 
      DD.results.blood_5_DEDD$mean.fdr$disp.expHM.untr.tmm, 
      type='l', col=col_vector[4])
lines(DD.results.blood_5_DEDD$mean.discoveries$disp.expHM.log.tmm, 
      DD.results.blood_5_DEDD$mean.fdr$disp.expHM.log.tmm, 
      type='l', col=col_vector[5])
lines(DD.results.blood_5_DEDD$mean.discoveries$disp.lnHM.untr.tmm, 
      DD.results.blood_5_DEDD$mean.fdr$disp.lnHM.untr.tmm, 
      type='l', col=col_vector[6])
lines(DD.results.blood_5_DEDD$mean.discoveries$disp.lnHM.log.tmm, 
      DD.results.blood_5_DEDD$mean.fdr$disp.lnHM.log.tmm, 
      type='l', col=col_vector[7])
abline(h=0.05, col="lightgrey")
plot(DD.results.blood_10_DEDD$mean.discoveries$diffVar.tmm, 
     DD.results.blood_10_DEDD$mean.fdr$diffVar.tmm, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("10 samples per group\n", "TMM normalisation\n", 
                 "Differences in mean and dispersion"))
lines(DD.results.blood_10_DEDD$mean.discoveries$disp.MDSeq.zi.tmm, 
      DD.results.blood_10_DEDD$mean.fdr$disp.MDSeq.zi.tmm, 
      type='l', col=col_vector[2])
lines(DD.results.blood_10_DEDD$mean.discoveries$disp.MDSeq.nozi.tmm, 
      DD.results.blood_10_DEDD$mean.fdr$disp.MDSeq.nozi.tmm, 
      type='l', col=col_vector[3])
lines(DD.results.blood_10_DEDD$mean.discoveries$disp.expHM.untr.tmm, 
      DD.results.blood_10_DEDD$mean.fdr$disp.expHM.untr.tmm, 
      type='l', col=col_vector[4])
lines(DD.results.blood_10_DEDD$mean.discoveries$disp.expHM.log.tmm, 
      DD.results.blood_10_DEDD$mean.fdr$disp.expHM.log.tmm, 
      type='l', col=col_vector[5])
lines(DD.results.blood_10_DEDD$mean.discoveries$disp.lnHM.untr.tmm, 
      DD.results.blood_10_DEDD$mean.fdr$disp.lnHM.untr.tmm, 
      type='l', col=col_vector[6])
lines(DD.results.blood_10_DEDD$mean.discoveries$disp.lnHM.log.tmm, 
      DD.results.blood_10_DEDD$mean.fdr$disp.lnHM.log.tmm, 
      type='l', col=col_vector[7])
abline(h=0.05, col="lightgrey")
plot(DD.results.blood_20_DEDD$mean.discoveries$diffVar.tmm, 
     DD.results.blood_20_DEDD$mean.fdr$diffVar.tmm, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("20 samples per group\n", "TMM normalisation\n", 
                 "Differences in mean and dispersion"))
lines(DD.results.blood_20_DEDD$mean.discoveries$disp.MDSeq.zi.tmm, 
      DD.results.blood_20_DEDD$mean.fdr$disp.MDSeq.zi.tmm, 
      type='l', col=col_vector[2])
lines(DD.results.blood_20_DEDD$mean.discoveries$disp.MDSeq.nozi.tmm, 
      DD.results.blood_20_DEDD$mean.fdr$disp.MDSeq.nozi.tmm, 
      type='l', col=col_vector[3])
lines(DD.results.blood_20_DEDD$mean.discoveries$disp.expHM.untr.tmm, 
      DD.results.blood_20_DEDD$mean.fdr$disp.expHM.untr.tmm, 
      type='l', col=col_vector[4])
lines(DD.results.blood_20_DEDD$mean.discoveries$disp.expHM.log.tmm, 
      DD.results.blood_20_DEDD$mean.fdr$disp.expHM.log.tmm, 
      type='l', col=col_vector[5])
lines(DD.results.blood_20_DEDD$mean.discoveries$disp.lnHM.untr.tmm, 
      DD.results.blood_20_DEDD$mean.fdr$disp.lnHM.untr.tmm, 
      type='l', col=col_vector[6])
lines(DD.results.blood_20_DEDD$mean.discoveries$disp.lnHM.log.tmm, 
      DD.results.blood_20_DEDD$mean.fdr$disp.lnHM.log.tmm, 
      type='l', col=col_vector[7])
abline(h=0.05, col="lightgrey")
plot(DD.results.blood_50_DEDD$mean.discoveries$diffVar.tmm, 
     DD.results.blood_50_DEDD$mean.fdr$diffVar.tmm, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("50 samples per group\n", "TMM normalisation\n", 
                 "Differences in mean and dispersion"))
lines(DD.results.blood_50_DEDD$mean.discoveries$disp.MDSeq.zi.tmm, 
      DD.results.blood_50_DEDD$mean.fdr$disp.MDSeq.zi.tmm, 
      type='l', col=col_vector[2])
lines(DD.results.blood_50_DEDD$mean.discoveries$disp.MDSeq.nozi.tmm, 
      DD.results.blood_50_DEDD$mean.fdr$disp.MDSeq.nozi.tmm, 
      type='l', col=col_vector[3])
lines(DD.results.blood_50_DEDD$mean.discoveries$disp.expHM.untr.tmm, 
      DD.results.blood_50_DEDD$mean.fdr$disp.expHM.untr.tmm, 
      type='l', col=col_vector[4])
lines(DD.results.blood_50_DEDD$mean.discoveries$disp.expHM.log.tmm, 
      DD.results.blood_50_DEDD$mean.fdr$disp.expHM.log.tmm, 
      type='l', col=col_vector[5])
lines(DD.results.blood_50_DEDD$mean.discoveries$disp.lnHM.untr.tmm, 
      DD.results.blood_50_DEDD$mean.fdr$disp.lnHM.untr.tmm, 
      type='l', col=col_vector[6])
lines(DD.results.blood_50_DEDD$mean.discoveries$disp.lnHM.log.tmm, 
      DD.results.blood_50_DEDD$mean.fdr$disp.lnHM.log.tmm, 
      type='l', col=col_vector[7])
abline(h=0.05, col="lightgrey")

plot(DD.results.blood_2_DEDD$mean.discoveries$diffVar.rle, 
     DD.results.blood_2_DEDD$mean.fdr$diffVar.rle, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("2 samples per group\n", "RLE normalisation\n", 
                 "Differences in mean and dispersion"))
lines(DD.results.blood_2_DEDD$mean.discoveries$disp.MDSeq.zi.rle, 
      DD.results.blood_2_DEDD$mean.fdr$disp.MDSeq.zi.rle, 
      type='l', col=col_vector[2])
lines(DD.results.blood_2_DEDD$mean.discoveries$disp.MDSeq.nozi.rle, 
      DD.results.blood_2_DEDD$mean.fdr$disp.MDSeq.nozi.rle, 
      type='l', col=col_vector[3])
lines(DD.results.blood_2_DEDD$mean.discoveries$disp.expHM.untr.rle, 
      DD.results.blood_2_DEDD$mean.fdr$disp.expHM.untr.rle, 
      type='l', col=col_vector[4])
lines(DD.results.blood_2_DEDD$mean.discoveries$disp.expHM.log.rle, 
      DD.results.blood_2_DEDD$mean.fdr$disp.expHM.log.rle, 
      type='l', col=col_vector[5])
lines(DD.results.blood_2_DEDD$mean.discoveries$disp.lnHM.untr.rle, 
      DD.results.blood_2_DEDD$mean.fdr$disp.lnHM.untr.rle, 
      type='l', col=col_vector[6])
lines(DD.results.blood_2_DEDD$mean.discoveries$disp.lnHM.log.rle, 
      DD.results.blood_2_DEDD$mean.fdr$disp.lnHM.log.rle, 
      type='l', col=col_vector[7])
lines(c(0,1000), c(0.05,0.05), col="lightgrey")
legend("bottomright", bty='n', col=col_vector[1:7], lty=1, ncol=1, lwd=2, 
       legend=c("diffVar", "MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", 
                "lnHM", "lnHM log"))
plot(DD.results.blood_5_DEDD$mean.discoveries$diffVar.rle, 
     DD.results.blood_5_DEDD$mean.fdr$diffVar.rle, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("5 samples per group\n", "RLE normalisation\n", 
                 "Differences in mean and dispersion"))
lines(DD.results.blood_5_DEDD$mean.discoveries$disp.MDSeq.zi.rle, 
      DD.results.blood_5_DEDD$mean.fdr$disp.MDSeq.zi.rle, 
      type='l', col=col_vector[2])
lines(DD.results.blood_5_DEDD$mean.discoveries$disp.MDSeq.nozi.rle, 
      DD.results.blood_5_DEDD$mean.fdr$disp.MDSeq.nozi.rle, 
      type='l', col=col_vector[3])
lines(DD.results.blood_5_DEDD$mean.discoveries$disp.expHM.untr.rle, 
      DD.results.blood_5_DEDD$mean.fdr$disp.expHM.untr.rle, 
      type='l', col=col_vector[4])
lines(DD.results.blood_5_DEDD$mean.discoveries$disp.expHM.log.rle, 
      DD.results.blood_5_DEDD$mean.fdr$disp.expHM.log.rle, 
      type='l', col=col_vector[5])
lines(DD.results.blood_5_DEDD$mean.discoveries$disp.lnHM.untr.rle, 
      DD.results.blood_5_DEDD$mean.fdr$disp.lnHM.untr.rle, 
      type='l', col=col_vector[6])
lines(DD.results.blood_5_DEDD$mean.discoveries$disp.lnHM.log.rle, 
      DD.results.blood_5_DEDD$mean.fdr$disp.lnHM.log.rle, 
      type='l', col=col_vector[7])
abline(h=0.05, col='lightgrey')
plot(DD.results.blood_10_DEDD$mean.discoveries$diffVar.rle, 
     DD.results.blood_10_DEDD$mean.fdr$diffVar.rle, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("10 samples per group\n", "RLE normalisation\n", 
                 "Differences in mean and dispersion"))
lines(DD.results.blood_10_DEDD$mean.discoveries$disp.MDSeq.zi.rle, 
      DD.results.blood_10_DEDD$mean.fdr$disp.MDSeq.zi.rle, 
      type='l', col=col_vector[2])
lines(DD.results.blood_10_DEDD$mean.discoveries$disp.MDSeq.nozi.rle, 
      DD.results.blood_10_DEDD$mean.fdr$disp.MDSeq.nozi.rle, 
      type='l', col=col_vector[3])
lines(DD.results.blood_10_DEDD$mean.discoveries$disp.expHM.untr.rle, 
      DD.results.blood_10_DEDD$mean.fdr$disp.expHM.untr.rle, 
      type='l', col=col_vector[4])
lines(DD.results.blood_10_DEDD$mean.discoveries$disp.expHM.log.rle, 
      DD.results.blood_10_DEDD$mean.fdr$disp.expHM.log.rle, 
      type='l', col=col_vector[5])
lines(DD.results.blood_10_DEDD$mean.discoveries$disp.lnHM.untr.rle, 
      DD.results.blood_10_DEDD$mean.fdr$disp.lnHM.untr.rle, 
      type='l', col=col_vector[6])
lines(DD.results.blood_10_DEDD$mean.discoveries$disp.lnHM.log.rle, 
      DD.results.blood_10_DEDD$mean.fdr$disp.lnHM.log.rle, 
      type='l', col=col_vector[7])
abline(h=0.05, col='lightgrey')
plot(DD.results.blood_20_DEDD$mean.discoveries$diffVar.rle, 
     DD.results.blood_20_DEDD$mean.fdr$diffVar.rle, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("20 samples per group\n", "RLE normalisation\n", 
                 "Differences in mean and dispersion"))
lines(DD.results.blood_20_DEDD$mean.discoveries$disp.MDSeq.zi.rle, 
      DD.results.blood_20_DEDD$mean.fdr$disp.MDSeq.zi.rle, 
      type='l', col=col_vector[2])
lines(DD.results.blood_20_DEDD$mean.discoveries$disp.MDSeq.nozi.rle, 
      DD.results.blood_20_DEDD$mean.fdr$disp.MDSeq.nozi.rle, 
      type='l', col=col_vector[3])
lines(DD.results.blood_20_DEDD$mean.discoveries$disp.expHM.untr.rle, 
      DD.results.blood_20_DEDD$mean.fdr$disp.expHM.untr.rle, 
      type='l', col=col_vector[4])
lines(DD.results.blood_20_DEDD$mean.discoveries$disp.expHM.log.rle, 
      DD.results.blood_20_DEDD$mean.fdr$disp.expHM.log.rle, 
      type='l', col=col_vector[5])
lines(DD.results.blood_20_DEDD$mean.discoveries$disp.lnHM.untr.rle, 
      DD.results.blood_20_DEDD$mean.fdr$disp.lnHM.untr.rle, 
      type='l', col=col_vector[6])
lines(DD.results.blood_20_DEDD$mean.discoveries$disp.lnHM.log.rle, 
      DD.results.blood_20_DEDD$mean.fdr$disp.lnHM.log.rle, 
      type='l', col=col_vector[7])
abline(h=0.05, col='lightgrey')
plot(DD.results.blood_50_DEDD$mean.discoveries$diffVar.rle, 
     DD.results.blood_50_DEDD$mean.fdr$diffVar.rle, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("50 samples per group\n", "RLE normalisation\n", 
                 "Differences in mean and dispersion"))
lines(DD.results.blood_50_DEDD$mean.discoveries$disp.MDSeq.nozi.rle, 
      DD.results.blood_50_DEDD$mean.fdr$disp.MDSeq.nozi.rle, 
      type='l', col=col_vector[2])
lines(DD.results.blood_50_DEDD$mean.discoveries$disp.MDSeq.zi.rle, 
      DD.results.blood_50_DEDD$mean.fdr$disp.MDSeq.zi.rle, 
      type='l', col=col_vector[3])
lines(DD.results.blood_50_DEDD$mean.discoveries$disp.expHM.untr.rle, 
      DD.results.blood_50_DEDD$mean.fdr$disp.expHM.untr.rle, 
      type='l', col=col_vector[4])
lines(DD.results.blood_50_DEDD$mean.discoveries$disp.expHM.log.rle, 
      DD.results.blood_50_DEDD$mean.fdr$disp.expHM.log.rle, 
      type='l', col=col_vector[5])
lines(DD.results.blood_50_DEDD$mean.discoveries$disp.lnHM.untr.rle, 
      DD.results.blood_50_DEDD$mean.fdr$disp.lnHM.untr.rle, 
      type='l', col=col_vector[6])
lines(DD.results.blood_50_DEDD$mean.discoveries$disp.lnHM.log.rle, 
      DD.results.blood_50_DEDD$mean.fdr$disp.lnHM.log.rle, 
      type='l', col=col_vector[7])
abline(h=0.05, col='lightgrey')



#################################
#### Differential expression ####
#################################

##############
#### AUC #####
par(mfrow=c(2,3), mar=c(0.5,2,5,1), mgp=c(3,0.7,0))
boxplot(DE.results.blood_2_DE$auc, names=NA, col=col_vector[1:14], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0.65,1.1), 
        main=paste0("AUC diff exp, 2 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean only"), 
        cex.main=1.7)
lines(c(14.5,14.5), c(0.6,1), col="lightgrey")
legend("topleft", fill=col_vector[1:14], bty='n', cex=1.2, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                "voom", "DSS", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DE.results.blood_5_DE$auc, names=NA, col=col_vector[1:14], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0.65,1.1), 
        main=paste0("AUC diff exp, 5 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean only"), 
        cex.main=1.7)
lines(c(14.5,14.5), c(0.6,1), col="lightgrey")
legend("topleft", fill=col_vector[1:14], bty='n', cex=1.2, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                "voom", "DSS", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DE.results.blood_10_DE$auc, names=NA, col=col_vector[1:14], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0.65,1.1), 
        main=paste0("AUC diff exp, 10 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean only"), 
        cex.main=1.7)
lines(c(14.5,14.5), c(0.6,1), col="lightgrey")
legend("topleft", fill=col_vector[1:14], bty='n', cex=1.2, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                "voom", "DSS", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DE.results.blood_20_DE$auc, names=NA, col=col_vector[1:14], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0.65,1.1), 
        main=paste0("AUC diff exp, 20 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean only"), 
        cex.main=1.7)
lines(c(14.5,14.5), c(0.6,1), col="lightgrey")
legend("topleft", fill=col_vector[1:14], bty='n', cex=1.2, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                "voom", "DSS", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DE.results.blood_50_DE$auc, names=NA, col=col_vector[1:14], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0.65,1.1), 
        main=paste0("AUC diff exp, 50 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean only"), 
        cex.main=1.7)
lines(c(14.5,14.5), c(0.6,1), col="lightgrey")
legend("topleft", fill=col_vector[1:14], bty='n', cex=1.2, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                "voom", "DSS", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                "expHM", "expHM log", "lnHM", "lnHM log"))
rbind(colMeans(DE.results.blood_2_DE$auc), 
      colMeans(DE.results.blood_5_DE$auc), 
      colMeans(DE.results.blood_10_DE$auc), 
      colMeans(DE.results.blood_20_DE$auc), 
      colMeans(DE.results.blood_50_DE$auc))

par(mfrow=c(2,3), mar=c(0.5,2,5,1), mgp=c(3,0.7,0))
boxplot(DE.results.blood_2_DEDD$auc, names=NA, col=col_vector[1:14], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0.65,1.1), 
        main=paste0("AUC diff exp, 2 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean and dispersion"), 
        cex.main=1.7)
lines(c(14.5,14.5), c(0.6,1), col="lightgrey")
legend("topleft", fill=col_vector[1:14], bty='n', cex=1.4, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                "voom", "DSS", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DE.results.blood_5_DEDD$auc, names=NA, col=col_vector[1:14], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0.65,1.1), 
        main=paste0("AUC diff exp, 5 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean and dispersion"), 
        cex.main=1.7)
lines(c(14.5,14.5), c(0.6,1), col="lightgrey")
legend("topleft", fill=col_vector[1:14], bty='n', cex=1.4, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                "voom", "DSS", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DE.results.blood_10_DEDD$auc, names=NA, col=col_vector[1:14], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0.65,1.1), 
        main=paste0("AUC diff exp, 10 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean and dispersion"), 
        cex.main=1.7)
lines(c(14.5,14.5), c(0.6,1), col="lightgrey")
legend("topleft", fill=col_vector[1:14], bty='n', cex=1.4, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                "voom", "DSS", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DE.results.blood_20_DEDD$auc, names=NA, col=col_vector[1:14], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0.65,1.1), 
        main=paste0("AUC diff exp, 20 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean and dispersion"), 
        cex.main=1.7)
lines(c(14.5,14.5), c(0.6,1), col="lightgrey")
legend("topleft", fill=col_vector[1:14], bty='n', cex=1.4, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                "voom", "DSS", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DE.results.blood_50_DEDD$auc, names=NA, col=col_vector[1:14], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0.65,1.1), 
        main=paste0("AUC diff exp, 50 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean and dispersion"), 
        cex.main=1.7)
lines(c(14.5,14.5), c(0.6,1), col="lightgrey")
legend("topleft", fill=col_vector[1:14], bty='n', cex=1.4, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                "voom", "DSS", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                "expHM", "expHM log", "lnHM", "lnHM log"))
rbind(colMeans(DE.results.blood_2_DEDD$auc), 
      colMeans(DE.results.blood_5_DEDD$auc), 
      colMeans(DE.results.blood_10_DEDD$auc), 
      colMeans(DE.results.blood_20_DEDD$auc), 
      colMeans(DE.results.blood_50_DEDD$auc))


##############
#### pAUC ####
par(mfrow=c(2,3), mar=c(0.5,2,5,1), mgp=c(3,0.7,0))
boxplot(DE.results.blood_2_DE$pauc, names=NA, col=col_vector[1:14], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,0.065), 
        main=paste0("Partial AUC diff exp, 2 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean only"), 
        cex.main=1.7)
lines(c(14.5,14.5), c(0,0.05), col="lightgrey")
legend("topleft", fill=col_vector[1:14], bty='n', cex=1.4, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                "voom", "DSS", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DE.results.blood_5_DE$pauc, names=NA, col=col_vector[1:14], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,0.065), 
        main=paste0("Partial AUC diff exp, 5 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean only"), 
        cex.main=1.7)
lines(c(14.5,14.5), c(0,0.05), col="lightgrey")
legend("topleft", fill=col_vector[1:14], bty='n', cex=1.4, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                "voom", "DSS", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DE.results.blood_10_DE$pauc, names=NA, col=col_vector[1:14], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,0.065), 
        main=paste0("Partial AUC diff exp, 10 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean only"), 
        cex.main=1.7)
lines(c(14.5,14.5), c(0,0.05), col="lightgrey")
legend("topleft", fill=col_vector[1:14], bty='n', cex=1.4, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                "voom", "DSS", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DE.results.blood_20_DE$pauc, names=NA, col=col_vector[1:14], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,0.065), 
        main=paste0("Partial AUC diff exp, 20 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean only"), 
        cex.main=1.7)
lines(c(14.5,14.5), c(0,0.05), col="lightgrey")
legend("topleft", fill=col_vector[1:14], bty='n', cex=1.4, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                "voom", "DSS", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DE.results.blood_50_DE$pauc, names=NA, col=col_vector[1:14], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,0.065), 
        main=paste0("Partial AUC diff exp, 50 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean only"), 
        cex.main=1.7)
lines(c(14.5,14.5), c(0,0.05), col="lightgrey")
legend("topleft", fill=col_vector[1:14], bty='n', cex=1.4, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                "voom", "DSS", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                "expHM", "expHM log", "lnHM", "lnHM log"))
rbind(colMeans(DE.results.blood_2_DE$pauc), 
      colMeans(DE.results.blood_5_DE$pauc), 
      colMeans(DE.results.blood_10_DE$pauc), 
      colMeans(DE.results.blood_20_DE$pauc), 
      colMeans(DE.results.blood_50_DE$pauc))

par(mfrow=c(2,3), mar=c(0.5,2,5,1), mgp=c(3,0.7,0))
boxplot(DE.results.blood_2_DEDD$pauc, names=NA, col=col_vector[1:14], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,0.065), 
        main=paste0("Partial AUC diff exp, 2 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean and dispersion"), 
        cex.main=1.7)
lines(c(14.5,14.5), c(0,0.05), col="lightgrey")
legend("topleft", fill=col_vector[1:14], bty='n', cex=1.4, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                "voom", "DSS", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DE.results.blood_5_DEDD$pauc, names=NA, col=col_vector[1:14], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,0.065), 
        main=paste0("Partial AUC diff exp, 5 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean and dispersion"), 
        cex.main=1.7)
lines(c(14.5,14.5), c(0,0.05), col="lightgrey")
legend("topleft", fill=col_vector[1:14], bty='n', cex=1.4, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                "voom", "DSS", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DE.results.blood_10_DEDD$pauc, names=NA, col=col_vector[1:14], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,0.065), 
        main=paste0("Partial AUC diff exp, 10 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean and dispersion"), 
        cex.main=1.7)
lines(c(14.5,14.5), c(0,0.05), col="lightgrey")
legend("topleft", fill=col_vector[1:14], bty='n', cex=1.4, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                "voom", "DSS", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DE.results.blood_20_DEDD$pauc, names=NA, col=col_vector[1:14], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,0.065), 
        main=paste0("Partial AUC diff exp, 20 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean and dispersion"), 
        cex.main=1.7)
lines(c(14.5,14.5), c(0,0.05), col="lightgrey")
legend("topleft", fill=col_vector[1:14], bty='n', cex=1.4, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                "voom", "DSS", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DE.results.blood_50_DEDD$pauc, names=NA, col=col_vector[1:14], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,0.065), 
        main=paste0("Partial AUC diff exp, 50 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean and dispersion"), 
        cex.main=1.7)
lines(c(14.5,14.5), c(0,0.05), col="lightgrey")
legend("topleft", fill=col_vector[1:14], bty='n', cex=1.4, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                "voom", "DSS", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                "expHM", "expHM log", "lnHM", "lnHM log"))
rbind(colMeans(DE.results.blood_2_DEDD$pauc), 
      colMeans(DE.results.blood_5_DEDD$pauc), 
      colMeans(DE.results.blood_10_DEDD$pauc), 
      colMeans(DE.results.blood_20_DEDD$pauc), 
      colMeans(DE.results.blood_50_DEDD$pauc))


#############
#### FDR ####
par(mfrow=c(2,3), mar=c(0.5,2,5,1), mgp=c(3,0.7,0))
boxplot(DE.results.blood_2_DE$fdr, names=NA, col=col_vector[1:15], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.3), 
        main=paste0("FDR diff exp, 2 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean only"), 
        cex.main=1.7)
lines(c(15.5,15.5), c(0,1), col="lightgrey"); abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[1:15], bty='n', cex=1.4, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                "voom", "DSS", "DSS lfdr", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DE.results.blood_5_DE$fdr, names=NA, col=col_vector[1:15], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.3), 
        main=paste0("FDR diff exp, 5 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean only"), 
        cex.main=1.7)
lines(c(15.5,15.5), c(0,1), col="lightgrey"); abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[1:15], bty='n', cex=1.4, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                "voom", "DSS", "DSS lfdr", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DE.results.blood_10_DE$fdr, names=NA, col=col_vector[1:15], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.3), 
        main=paste0("FDR diff exp, 10 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean only"), 
        cex.main=1.7)
lines(c(15.5,15.5), c(0,1), col="lightgrey"); abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[1:15], bty='n', cex=1.4, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                "voom", "DSS", "DSS lfdr", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DE.results.blood_20_DE$fdr, names=NA, col=col_vector[1:15], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.3), 
        main=paste0("FDR diff exp, 20 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean only"), 
        cex.main=1.7)
lines(c(15.5,15.5), c(0,1), col="lightgrey"); abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[1:15], bty='n', cex=1.4, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                "voom", "DSS", "DSS lfdr", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DE.results.blood_50_DE$fdr, names=NA, col=col_vector[1:15], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.3), 
        main=paste0("FDR diff exp, 50 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean only"), 
        cex.main=1.7)
lines(c(15.5,15.5), c(0,1), col="lightgrey"); abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[1:15], bty='n', cex=1.4, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                "voom", "DSS", "DSS lfdr", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                "expHM", "expHM log", "lnHM", "lnHM log"))
rbind(colMeans(DE.results.blood_2_DE$fdr), 
      colMeans(DE.results.blood_5_DE$fdr), 
      colMeans(DE.results.blood_10_DE$fdr), 
      colMeans(DE.results.blood_20_DE$fdr), 
      colMeans(DE.results.blood_50_DE$fdr))

par(mfrow=c(2,3), mar=c(0.5,2,5,1), mgp=c(3,0.7,0))
boxplot(DE.results.blood_2_DEDD$fdr, names=NA, col=col_vector[1:15], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.3), 
        main=paste0("FDR diff exp, 2 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean and dispersion"), 
        cex.main=1.7)
lines(c(15.5,15.5), c(0,1), col="lightgrey"); abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[1:15], bty='n', cex=1.4, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                "voom", "DSS", "DSS lfdr", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DE.results.blood_5_DEDD$fdr, names=NA, col=col_vector[1:15], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.3), 
        main=paste0("FDR diff exp, 5 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean and dispersion"), 
        cex.main=1.7)
lines(c(15.5,15.5), c(0,1), col="lightgrey"); abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[1:15], bty='n', cex=1.4, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                "voom", "DSS", "DSS lfdr", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DE.results.blood_10_DEDD$fdr, names=NA, col=col_vector[1:15], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.3), 
        main=paste0("FDR diff exp, 10 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean and dispersion"), 
        cex.main=1.7)
lines(c(15.5,15.5), c(0,1), col="lightgrey"); abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[1:15], bty='n', cex=1.4, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                "voom", "DSS", "DSS lfdr", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DE.results.blood_20_DEDD$fdr, names=NA, col=col_vector[1:15], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.3), 
        main=paste0("FDR diff exp, 20 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean and dispersion"), 
        cex.main=1.7)
lines(c(15.5,15.5), c(0,1), col="lightgrey"); abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[1:15], bty='n', cex=1.4, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                "voom", "DSS", "DSS lfdr", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DE.results.blood_50_DEDD$fdr, names=NA, col=col_vector[1:15], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.3), 
        main=paste0("FDR diff exp, 50 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean and dispersion"), 
        cex.main=1.7)
lines(c(15.5,15.5), c(0,1), col="lightgrey"); abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[1:15], bty='n', cex=1.4, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                "voom", "DSS", "DSS lfdr", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                "expHM", "expHM log", "lnHM", "lnHM log"))
rbind(colMeans(DE.results.blood_2_DEDD$fdr), 
      colMeans(DE.results.blood_5_DEDD$fdr), 
      colMeans(DE.results.blood_10_DEDD$fdr), 
      colMeans(DE.results.blood_20_DEDD$fdr), 
      colMeans(DE.results.blood_50_DEDD$fdr))


#############
#### TPR ####
par(mfrow=c(2,3), mar=c(0.5,2,5,1), mgp=c(3,0.7,0))
boxplot(DE.results.blood_2_DE$tpr, names=NA, col=col_vector[1:15], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.25), 
        main=paste0("TPR diff exp, 2 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean only"), 
        cex.main=1.7)
lines(c(15.5,15.5), c(0,0.95), col="lightgrey")
legend("topleft", fill=col_vector[1:15], bty='n', cex=1.4, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                "voom", "DSS", "DSS lfdr", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DE.results.blood_5_DE$tpr, names=NA, col=col_vector[1:15], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.25), 
        main=paste0("TPR diff exp, 5 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean only"), 
        cex.main=1.7)
lines(c(15.5,15.5), c(0,0.95), col="lightgrey")
legend("topleft", fill=col_vector[1:15], bty='n', cex=1.4, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                "voom", "DSS", "DSS lfdr", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DE.results.blood_10_DE$tpr, names=NA, col=col_vector[1:15], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.25), 
        main=paste0("TPR diff exp, 10 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean only"), 
        cex.main=1.7)
lines(c(15.5,15.5), c(0,0.95), col="lightgrey")
legend("topleft", fill=col_vector[1:15], bty='n', cex=1.4, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                "voom", "DSS", "DSS lfdr", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DE.results.blood_20_DE$tpr, names=NA, col=col_vector[1:15], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.25), 
        main=paste0("TPR diff exp, 20 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean only"), 
        cex.main=1.7)
lines(c(15.5,15.5), c(0,0.95), col="lightgrey")
legend("topleft", fill=col_vector[1:15], bty='n', cex=1.4, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                "voom", "DSS", "DSS lfdr", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DE.results.blood_50_DE$tpr, names=NA, col=col_vector[1:15], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.25), 
        main=paste0("TPR diff exp, 50 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean only"), 
        cex.main=1.7)
lines(c(15.5,15.5), c(0,0.95), col="lightgrey")
legend("topleft", fill=col_vector[1:15], bty='n', cex=1.4, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                "voom", "DSS", "DSS lfdr", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                "expHM", "expHM log", "lnHM", "lnHM log"))
rbind(colMeans(DE.results.blood_2_DE$tpr), colMeans(DE.results.blood_5_DE$tpr), 
      colMeans(DE.results.blood_10_DE$tpr), colMeans(DE.results.blood_20_DE$tpr), 
      colMeans(DE.results.blood_50_DE$tpr))

par(mfrow=c(2,3), mar=c(0.5,2,5,1), mgp=c(3,0.7,0))
boxplot(DE.results.blood_2_DEDD$tpr, names=NA, col=col_vector[1:15], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.25), 
        main=paste0("TPR diff exp, 2 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean and dispersion"), 
        cex.main=1.7)
lines(c(15.5,15.5), c(0,0.95), col="lightgrey")
legend("topleft", fill=col_vector[1:15], bty='n', cex=1.4, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                "voom", "DSS", "DSS lfdr", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DE.results.blood_5_DEDD$tpr, names=NA, col=col_vector[1:15], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.25), 
        main=paste0("TPR diff exp, 5 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean and dispersion"), 
        cex.main=1.7)
lines(c(15.5,15.5), c(0,0.95), col="lightgrey")
legend("topleft", fill=col_vector[1:15], bty='n', cex=1.4, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                "voom", "DSS", "DSS lfdr", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DE.results.blood_10_DEDD$tpr, names=NA, col=col_vector[1:15], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.25), 
        main=paste0("TPR diff exp, 10 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean and dispersion"), 
        cex.main=1.7)
lines(c(15.5,15.5), c(0,0.95), col="lightgrey")
legend("topleft", fill=col_vector[1:15], bty='n', cex=1.4, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                "voom", "DSS", "DSS lfdr", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DE.results.blood_20_DEDD$tpr, names=NA, col=col_vector[1:15], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.25), 
        main=paste0("TPR diff exp, 20 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean and dispersion"), 
        cex.main=1.7)
lines(c(15.5,15.5), c(0,0.95), col="lightgrey")
legend("topleft", fill=col_vector[1:15], bty='n', cex=1.4, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                "voom", "DSS", "DSS lfdr", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DE.results.blood_50_DEDD$tpr, names=NA, col=col_vector[1:15], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.25), 
        main=paste0("TPR diff exp, 50 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean and dispersion"), 
        cex.main=1.7)
lines(c(15.5,15.5), c(0,0.95), col="lightgrey")
legend("topleft", fill=col_vector[1:15], bty='n', cex=1.4, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                "voom", "DSS", "DSS lfdr", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                "expHM", "expHM log", "lnHM", "lnHM log"))
rbind(colMeans(DE.results.blood_2_DEDD$tpr), 
      colMeans(DE.results.blood_5_DEDD$tpr), 
      colMeans(DE.results.blood_10_DEDD$tpr), 
      colMeans(DE.results.blood_20_DEDD$tpr), 
      colMeans(DE.results.blood_50_DEDD$tpr))


###############################
#### False discovery plots ####
par(mfrow=c(4,5), mar=c(2,1.5,2.5,0.5), mgp=c(3,0.5,0))

plot(DE.results.blood_2_DE$mean.discoveries$edgeR.ql.tmm, 
     DE.results.blood_2_DE$mean.fdr$edgeR.ql.tmm, 
     type='l', ylim=c(0,0.7), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("2 samples per group, TMM\n", 
                 "Differences in mean only"))
lines(DE.results.blood_2_DE$mean.discoveries$edgeR.lr.tmm, 
      DE.results.blood_2_DE$mean.fdr$edgeR.lr.tmm, 
      type='l', col=col_vector[2])
lines(DE.results.blood_2_DE$mean.discoveries$edgeR.et.tmm, 
      DE.results.blood_2_DE$mean.fdr$edgeR.et.tmm, 
      type='l', col=col_vector[3])
lines(DE.results.blood_2_DE$mean.discoveries$DESeq2.noif.tmm, 
      DE.results.blood_2_DE$mean.fdr$DESeq2.noif.tmm, 
      type='l', col=col_vector[4])
lines(DE.results.blood_2_DE$mean.discoveries$DESeq2.if.tmm, 
      DE.results.blood_2_DE$mean.fdr$DESeq2.if.tmm, 
      type='l', col=col_vector[5])
lines(DE.results.blood_2_DE$mean.discoveries$voom.tmm, 
      DE.results.blood_2_DE$mean.fdr$voom.tmm, 
      type='l', col=col_vector[6])
lines(DE.results.blood_2_DE$mean.discoveries$DSS.tmm, 
      DE.results.blood_2_DE$mean.fdr$DSS.tmm, 
      type='l', col=col_vector[7])
lines(DE.results.blood_2_DE$mean.discoveries$baySeq.tmm, 
      DE.results.blood_2_DE$mean.fdr$baySeq.tmm, 
      type='l', col=col_vector[8])
lines(DE.results.blood_2_DE$mean.discoveries$mean.MDSeq.zi.tmm, 
      DE.results.blood_2_DE$mean.fdr$mean.MDSeq.zi.tmm, 
      type='l', col=col_vector[9])
lines(DE.results.blood_2_DE$mean.discoveries$mean.MDSeq.nozi.tmm, 
      DE.results.blood_2_DE$mean.fdr$mean.MDSeq.nozi.tmm, 
      type='l', col=col_vector[10])
lines(DE.results.blood_2_DE$mean.discoveries$mean.expHM.untr.tmm, 
      DE.results.blood_2_DE$mean.fdr$mean.expHM.untr.tmm, 
      type='l', col=col_vector[11])
lines(DE.results.blood_2_DE$mean.discoveries$mean.expHM.log.tmm, 
      DE.results.blood_2_DE$mean.fdr$mean.expHM.log.tmm, 
      type='l', col=col_vector[12])
lines(DE.results.blood_2_DE$mean.discoveries$mean.lnHM.untr.tmm, 
      DE.results.blood_2_DE$mean.fdr$mean.lnHM.untr.tmm, 
      type='l', col=col_vector[13])
lines(DE.results.blood_2_DE$mean.discoveries$mean.lnHM.log.tmm, 
      DE.results.blood_2_DE$mean.fdr$mean.lnHM.log.tmm, 
      type='l', col=col_vector[14])
abline(h=0.05, col='lightgrey')
plot(DE.results.blood_5_DE$mean.discoveries$edgeR.ql.tmm, 
     DE.results.blood_5_DE$mean.fdr$edgeR.ql.tmm, 
     type='l', ylim=c(0,0.7), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("5 samples per group, TMM\n", 
                 "Differences in mean only"))
lines(DE.results.blood_5_DE$mean.discoveries$edgeR.lr.tmm, 
      DE.results.blood_5_DE$mean.fdr$edgeR.lr.tmm, 
      type='l', col=col_vector[2])
lines(DE.results.blood_5_DE$mean.discoveries$edgeR.et.tmm, 
      DE.results.blood_5_DE$mean.fdr$edgeR.et.tmm, 
      type='l', col=col_vector[3])
lines(DE.results.blood_5_DE$mean.discoveries$DESeq2.noif.tmm, 
      DE.results.blood_5_DE$mean.fdr$DESeq2.noif.tmm, 
      type='l', col=col_vector[4])
lines(DE.results.blood_5_DE$mean.discoveries$DESeq2.if.tmm, 
      DE.results.blood_5_DE$mean.fdr$DESeq2.if.tmm, 
      type='l', col=col_vector[5])
lines(DE.results.blood_5_DE$mean.discoveries$voom.tmm, 
      DE.results.blood_5_DE$mean.fdr$voom.tmm, 
      type='l', col=col_vector[6])
lines(DE.results.blood_5_DE$mean.discoveries$DSS.tmm, 
      DE.results.blood_5_DE$mean.fdr$DSS.tmm, 
      type='l', col=col_vector[7])
lines(DE.results.blood_5_DE$mean.discoveries$baySeq.tmm, 
      DE.results.blood_5_DE$mean.fdr$baySeq.tmm, 
      type='l', col=col_vector[8])
lines(DE.results.blood_5_DE$mean.discoveries$mean.MDSeq.zi.tmm, 
      DE.results.blood_5_DE$mean.fdr$mean.MDSeq.zi.tmm, 
      type='l', col=col_vector[9])
lines(DE.results.blood_5_DE$mean.discoveries$mean.MDSeq.nozi.tmm, 
      DE.results.blood_5_DE$mean.fdr$mean.MDSeq.nozi.tmm, 
      type='l', col=col_vector[10])
lines(DE.results.blood_5_DE$mean.discoveries$mean.expHM.untr.tmm, 
      DE.results.blood_5_DE$mean.fdr$mean.expHM.untr.tmm, 
      type='l', col=col_vector[11])
lines(DE.results.blood_5_DE$mean.discoveries$mean.expHM.log.tmm, 
      DE.results.blood_5_DE$mean.fdr$mean.expHM.log.tmm, 
      type='l', col=col_vector[12])
lines(DE.results.blood_5_DE$mean.discoveries$mean.lnHM.untr.tmm, 
      DE.results.blood_5_DE$mean.fdr$mean.lnHM.untr.tmm, 
      type='l', col=col_vector[13])
lines(DE.results.blood_5_DE$mean.discoveries$mean.lnHM.log.tmm, 
      DE.results.blood_5_DE$mean.fdr$mean.lnHM.log.tmm, 
      type='l', col=col_vector[14])
abline(h=0.05, col='lightgrey')
plot(DE.results.blood_10_DE$mean.discoveries$edgeR.ql.tmm, 
     DE.results.blood_10_DE$mean.fdr$edgeR.ql.tmm, 
     type='l', ylim=c(0,0.7), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("10 samples per group, TMM\n", 
                 "Differences in mean only"))
lines(DE.results.blood_10_DE$mean.discoveries$edgeR.lr.tmm, 
      DE.results.blood_10_DE$mean.fdr$edgeR.lr.tmm, 
      type='l', col=col_vector[2])
lines(DE.results.blood_10_DE$mean.discoveries$edgeR.et.tmm, 
      DE.results.blood_10_DE$mean.fdr$edgeR.et.tmm, 
      type='l', col=col_vector[3])
lines(DE.results.blood_10_DE$mean.discoveries$DESeq2.noif.tmm, 
      DE.results.blood_10_DE$mean.fdr$DESeq2.noif.tmm, 
      type='l', col=col_vector[4])
lines(DE.results.blood_10_DE$mean.discoveries$DESeq2.if.tmm, 
      DE.results.blood_10_DE$mean.fdr$DESeq2.if.tmm, 
      type='l', col=col_vector[5])
lines(DE.results.blood_10_DE$mean.discoveries$voom.tmm, 
      DE.results.blood_10_DE$mean.fdr$voom.tmm, 
      type='l', col=col_vector[6])
lines(DE.results.blood_10_DE$mean.discoveries$DSS.tmm, 
      DE.results.blood_10_DE$mean.fdr$DSS.tmm, 
      type='l', col=col_vector[7])
lines(DE.results.blood_10_DE$mean.discoveries$baySeq.tmm, 
      DE.results.blood_10_DE$mean.fdr$baySeq.tmm, 
      type='l', col=col_vector[8])
lines(DE.results.blood_10_DE$mean.discoveries$mean.MDSeq.zi.tmm, 
      DE.results.blood_10_DE$mean.fdr$mean.MDSeq.zi.tmm, 
      type='l', col=col_vector[9])
lines(DE.results.blood_10_DE$mean.discoveries$mean.MDSeq.nozi.tmm, 
      DE.results.blood_10_DE$mean.fdr$mean.MDSeq.nozi.tmm, 
      type='l', col=col_vector[10])
lines(DE.results.blood_10_DE$mean.discoveries$mean.expHM.untr.tmm, 
      DE.results.blood_10_DE$mean.fdr$mean.expHM.untr.tmm, 
      type='l', col=col_vector[11])
lines(DE.results.blood_10_DE$mean.discoveries$mean.expHM.log.tmm, 
      DE.results.blood_10_DE$mean.fdr$mean.expHM.log.tmm, 
      type='l', col=col_vector[12])
lines(DE.results.blood_10_DE$mean.discoveries$mean.lnHM.untr.tmm, 
      DE.results.blood_10_DE$mean.fdr$mean.lnHM.untr.tmm, 
      type='l', col=col_vector[13])
lines(DE.results.blood_10_DE$mean.discoveries$mean.lnHM.log.tmm, 
      DE.results.blood_10_DE$mean.fdr$mean.lnHM.log.tmm, 
      type='l', col=col_vector[14])
abline(h=0.05, col='lightgrey')
plot(DE.results.blood_20_DE$mean.discoveries$edgeR.ql.tmm, 
     DE.results.blood_20_DE$mean.fdr$edgeR.ql.tmm, 
     type='l', ylim=c(0,0.7), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("20 samples per group, TMM\n", 
                 "Differences in mean only"))
lines(DE.results.blood_20_DE$mean.discoveries$edgeR.lr.tmm, 
      DE.results.blood_20_DE$mean.fdr$edgeR.lr.tmm, 
      type='l', col=col_vector[2])
lines(DE.results.blood_20_DE$mean.discoveries$edgeR.et.tmm, 
      DE.results.blood_20_DE$mean.fdr$edgeR.et.tmm, 
      type='l', col=col_vector[3])
lines(DE.results.blood_20_DE$mean.discoveries$DESeq2.noif.tmm, 
      DE.results.blood_20_DE$mean.fdr$DESeq2.noif.tmm, 
      type='l', col=col_vector[4])
lines(DE.results.blood_20_DE$mean.discoveries$DESeq2.if.tmm, 
      DE.results.blood_20_DE$mean.fdr$DESeq2.if.tmm, 
      type='l', col=col_vector[5])
lines(DE.results.blood_20_DE$mean.discoveries$voom.tmm, 
      DE.results.blood_20_DE$mean.fdr$voom.tmm, 
      type='l', col=col_vector[6])
lines(DE.results.blood_20_DE$mean.discoveries$DSS.tmm, 
      DE.results.blood_20_DE$mean.fdr$DSS.tmm, 
      type='l', col=col_vector[7])
lines(DE.results.blood_20_DE$mean.discoveries$baySeq.tmm, 
      DE.results.blood_20_DE$mean.fdr$baySeq.tmm, 
      type='l', col=col_vector[8])
lines(DE.results.blood_20_DE$mean.discoveries$mean.MDSeq.zi.tmm, 
      DE.results.blood_20_DE$mean.fdr$mean.MDSeq.zi.tmm, 
      type='l', col=col_vector[9])
lines(DE.results.blood_20_DE$mean.discoveries$mean.MDSeq.nozi.tmm, 
      DE.results.blood_20_DE$mean.fdr$mean.MDSeq.nozi.tmm, 
      type='l', col=col_vector[10])
lines(DE.results.blood_20_DE$mean.discoveries$mean.expHM.untr.tmm, 
      DE.results.blood_20_DE$mean.fdr$mean.expHM.untr.tmm, 
      type='l', col=col_vector[11])
lines(DE.results.blood_20_DE$mean.discoveries$mean.expHM.log.tmm, 
      DE.results.blood_20_DE$mean.fdr$mean.expHM.log.tmm, 
      type='l', col=col_vector[12])
lines(DE.results.blood_20_DE$mean.discoveries$mean.lnHM.untr.tmm, 
      DE.results.blood_20_DE$mean.fdr$mean.lnHM.untr.tmm, 
      type='l', col=col_vector[13])
lines(DE.results.blood_20_DE$mean.discoveries$mean.lnHM.log.tmm, 
      DE.results.blood_20_DE$mean.fdr$mean.lnHM.log.tmm, 
      type='l', col=col_vector[14])
abline(h=0.05, col='lightgrey')
plot(DE.results.blood_50_DE$mean.discoveries$edgeR.ql.tmm, 
     DE.results.blood_50_DE$mean.fdr$edgeR.ql.tmm, 
     type='l', ylim=c(0,0.7), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("50 samples per group, TMM\n", 
                 "Differences in mean only"))
lines(DE.results.blood_50_DE$mean.discoveries$edgeR.lr.tmm, 
      DE.results.blood_50_DE$mean.fdr$edgeR.lr.tmm, 
      type='l', col=col_vector[2])
lines(DE.results.blood_50_DE$mean.discoveries$edgeR.et.tmm, 
      DE.results.blood_50_DE$mean.fdr$edgeR.et.tmm, 
      type='l', col=col_vector[3])
lines(DE.results.blood_50_DE$mean.discoveries$DESeq2.noif.tmm, 
      DE.results.blood_50_DE$mean.fdr$DESeq2.noif.tmm, 
      type='l', col=col_vector[4])
lines(DE.results.blood_50_DE$mean.discoveries$DESeq2.if.tmm, 
      DE.results.blood_50_DE$mean.fdr$DESeq2.if.tmm, 
      type='l', col=col_vector[5])
lines(DE.results.blood_50_DE$mean.discoveries$voom.tmm, 
      DE.results.blood_50_DE$mean.fdr$voom.tmm, 
      type='l', col=col_vector[6])
lines(DE.results.blood_50_DE$mean.discoveries$DSS.tmm, 
      DE.results.blood_50_DE$mean.fdr$DSS.tmm, 
      type='l', col=col_vector[7])
lines(DE.results.blood_50_DE$mean.discoveries$baySeq.tmm, 
      DE.results.blood_50_DE$mean.fdr$baySeq.tmm, 
      type='l', col=col_vector[8])
lines(DE.results.blood_50_DE$mean.discoveries$mean.MDSeq.zi.tmm, 
      DE.results.blood_50_DE$mean.fdr$mean.MDSeq.zi.tmm, 
      type='l', col=col_vector[9])
lines(DE.results.blood_50_DE$mean.discoveries$mean.MDSeq.nozi.tmm, 
      DE.results.blood_50_DE$mean.fdr$mean.MDSeq.nozi.tmm, 
      type='l', col=col_vector[10])
lines(DE.results.blood_50_DE$mean.discoveries$mean.expHM.untr.tmm, 
      DE.results.blood_50_DE$mean.fdr$mean.expHM.untr.tmm, 
      type='l', col=col_vector[11])
lines(DE.results.blood_50_DE$mean.discoveries$mean.expHM.log.tmm, 
      DE.results.blood_50_DE$mean.fdr$mean.expHM.log.tmm, 
      type='l', col=col_vector[12])
lines(DE.results.blood_50_DE$mean.discoveries$mean.lnHM.untr.tmm, 
      DE.results.blood_50_DE$mean.fdr$mean.lnHM.untr.tmm, 
      type='l', col=col_vector[13])
lines(DE.results.blood_50_DE$mean.discoveries$mean.lnHM.log.tmm, 
      DE.results.blood_50_DE$mean.fdr$mean.lnHM.log.tmm, 
      type='l', col=col_vector[14])
abline(h=0.05, col='lightgrey')
legend("topleft", bty='n', col=col_vector[1:14], lty=1, ncol=2, lwd=2, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                "voom", "DSS", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                "expHM", "expHM log", "lnHM", "lnHM log"))

plot(DE.results.blood_2_DE$mean.discoveries$edgeR.ql.rle, 
     DE.results.blood_2_DE$mean.fdr$edgeR.ql.rle, 
     type='l', ylim=c(0,0.7), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("2 samples per group, RLE\n", 
                 "Differences in mean only"))
lines(DE.results.blood_2_DE$mean.discoveries$edgeR.lr.rle, 
      DE.results.blood_2_DE$mean.fdr$edgeR.lr.rle, 
      type='l', col=col_vector[2])
lines(DE.results.blood_2_DE$mean.discoveries$edgeR.et.rle, 
      DE.results.blood_2_DE$mean.fdr$edgeR.et.rle, 
      type='l', col=col_vector[3])
lines(DE.results.blood_2_DE$mean.discoveries$DESeq2.noif.rle, 
      DE.results.blood_2_DE$mean.fdr$DESeq2.noif.rle, 
      type='l', col=col_vector[4])
lines(DE.results.blood_2_DE$mean.discoveries$DESeq2.if.rle, 
      DE.results.blood_2_DE$mean.fdr$DESeq2.if.rle, 
      type='l', col=col_vector[5])
lines(DE.results.blood_2_DE$mean.discoveries$voom.rle, 
      DE.results.blood_2_DE$mean.fdr$voom.rle, 
      type='l', col=col_vector[6])
lines(DE.results.blood_2_DE$mean.discoveries$DSS.rle, 
      DE.results.blood_2_DE$mean.fdr$DSS.rle, 
      type='l', col=col_vector[7])
lines(DE.results.blood_2_DE$mean.discoveries$baySeq.rle, 
      DE.results.blood_2_DE$mean.fdr$baySeq.rle, 
      type='l', col=col_vector[8])
lines(DE.results.blood_2_DE$mean.discoveries$mean.MDSeq.zi.rle, 
      DE.results.blood_2_DE$mean.fdr$mean.MDSeq.zi.rle, 
      type='l', col=col_vector[9])
lines(DE.results.blood_2_DE$mean.discoveries$mean.MDSeq.nozi.rle, 
      DE.results.blood_2_DE$mean.fdr$mean.MDSeq.nozi.rle, 
      type='l', col=col_vector[10])
lines(DE.results.blood_2_DE$mean.discoveries$mean.expHM.untr.rle, 
      DE.results.blood_2_DE$mean.fdr$mean.expHM.untr.rle, 
      type='l', col=col_vector[11])
lines(DE.results.blood_2_DE$mean.discoveries$mean.expHM.log.rle, 
      DE.results.blood_2_DE$mean.fdr$mean.expHM.log.rle, 
      type='l', col=col_vector[12])
lines(DE.results.blood_2_DE$mean.discoveries$mean.lnHM.untr.rle, 
      DE.results.blood_2_DE$mean.fdr$mean.lnHM.untr.rle, 
      type='l', col=col_vector[13])
lines(DE.results.blood_2_DE$mean.discoveries$mean.lnHM.log.rle, 
      DE.results.blood_2_DE$mean.fdr$mean.lnHM.log.rle, 
      type='l', col=col_vector[14])
abline(h=0.05, col='lightgrey')
plot(DE.results.blood_5_DE$mean.discoveries$edgeR.ql.rle, 
     DE.results.blood_5_DE$mean.fdr$edgeR.ql.rle, 
     type='l', ylim=c(0,0.7), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("5 samples per group, RLE\n", 
                 "Differences in mean only"))
lines(DE.results.blood_5_DE$mean.discoveries$edgeR.lr.rle, 
      DE.results.blood_5_DE$mean.fdr$edgeR.lr.rle, 
      type='l', col=col_vector[2])
lines(DE.results.blood_5_DE$mean.discoveries$edgeR.et.rle, 
      DE.results.blood_5_DE$mean.fdr$edgeR.et.rle, 
      type='l', col=col_vector[3])
lines(DE.results.blood_5_DE$mean.discoveries$DESeq2.noif.rle, 
      DE.results.blood_5_DE$mean.fdr$DESeq2.noif.rle, 
      type='l', col=col_vector[4])
lines(DE.results.blood_5_DE$mean.discoveries$DESeq2.if.rle, 
      DE.results.blood_5_DE$mean.fdr$DESeq2.if.rle, 
      type='l', col=col_vector[5])
lines(DE.results.blood_5_DE$mean.discoveries$voom.rle, 
      DE.results.blood_5_DE$mean.fdr$voom.rle, 
      type='l', col=col_vector[6])
lines(DE.results.blood_5_DE$mean.discoveries$DSS.rle, 
      DE.results.blood_5_DE$mean.fdr$DSS.rle, 
      type='l', col=col_vector[7])
lines(DE.results.blood_5_DE$mean.discoveries$baySeq.rle, 
      DE.results.blood_5_DE$mean.fdr$baySeq.rle, 
      type='l', col=col_vector[8])
lines(DE.results.blood_5_DE$mean.discoveries$mean.MDSeq.zi.rle, 
      DE.results.blood_5_DE$mean.fdr$mean.MDSeq.zi.rle, 
      type='l', col=col_vector[9])
lines(DE.results.blood_5_DE$mean.discoveries$mean.MDSeq.nozi.rle, 
      DE.results.blood_5_DE$mean.fdr$mean.MDSeq.nozi.rle, 
      type='l', col=col_vector[10])
lines(DE.results.blood_5_DE$mean.discoveries$mean.expHM.untr.rle, 
      DE.results.blood_5_DE$mean.fdr$mean.expHM.untr.rle, 
      type='l', col=col_vector[11])
lines(DE.results.blood_5_DE$mean.discoveries$mean.expHM.log.rle, 
      DE.results.blood_5_DE$mean.fdr$mean.expHM.log.rle, 
      type='l', col=col_vector[12])
lines(DE.results.blood_5_DE$mean.discoveries$mean.lnHM.untr.rle, 
      DE.results.blood_5_DE$mean.fdr$mean.lnHM.untr.rle, 
      type='l', col=col_vector[13])
lines(DE.results.blood_5_DE$mean.discoveries$mean.lnHM.log.rle, 
      DE.results.blood_5_DE$mean.fdr$mean.lnHM.log.rle, 
      type='l', col=col_vector[14])
abline(h=0.05, col='lightgrey')
plot(DE.results.blood_10_DE$mean.discoveries$edgeR.ql.rle, 
     DE.results.blood_10_DE$mean.fdr$edgeR.ql.rle, 
     type='l', ylim=c(0,0.7), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("10 samples per group, RLE\n", 
                 "Differences in mean only"))
lines(DE.results.blood_10_DE$mean.discoveries$edgeR.lr.rle, 
      DE.results.blood_10_DE$mean.fdr$edgeR.lr.rle, 
      type='l', col=col_vector[2])
lines(DE.results.blood_10_DE$mean.discoveries$edgeR.et.rle, 
      DE.results.blood_10_DE$mean.fdr$edgeR.et.rle, 
      type='l', col=col_vector[3])
lines(DE.results.blood_10_DE$mean.discoveries$DESeq2.noif.rle, 
      DE.results.blood_10_DE$mean.fdr$DESeq2.noif.rle, 
      type='l', col=col_vector[4])
lines(DE.results.blood_10_DE$mean.discoveries$DESeq2.if.rle, 
      DE.results.blood_10_DE$mean.fdr$DESeq2.if.rle, 
      type='l', col=col_vector[5])
lines(DE.results.blood_10_DE$mean.discoveries$voom.rle, 
      DE.results.blood_10_DE$mean.fdr$voom.rle, 
      type='l', col=col_vector[6])
lines(DE.results.blood_10_DE$mean.discoveries$DSS.rle, 
      DE.results.blood_10_DE$mean.fdr$DSS.rle, 
      type='l', col=col_vector[7])
lines(DE.results.blood_10_DE$mean.discoveries$baySeq.rle, 
      DE.results.blood_10_DE$mean.fdr$baySeq.rle, 
      type='l', col=col_vector[8])
lines(DE.results.blood_10_DE$mean.discoveries$mean.MDSeq.zi.rle, 
      DE.results.blood_10_DE$mean.fdr$mean.MDSeq.zi.rle, 
      type='l', col=col_vector[9])
lines(DE.results.blood_10_DE$mean.discoveries$mean.MDSeq.nozi.rle, 
      DE.results.blood_10_DE$mean.fdr$mean.MDSeq.nozi.rle, 
      type='l', col=col_vector[10])
lines(DE.results.blood_10_DE$mean.discoveries$mean.expHM.untr.rle, 
      DE.results.blood_10_DE$mean.fdr$mean.expHM.untr.rle, 
      type='l', col=col_vector[11])
lines(DE.results.blood_10_DE$mean.discoveries$mean.expHM.log.rle, 
      DE.results.blood_10_DE$mean.fdr$mean.expHM.log.rle, 
      type='l', col=col_vector[12])
lines(DE.results.blood_10_DE$mean.discoveries$mean.lnHM.untr.rle, 
      DE.results.blood_10_DE$mean.fdr$mean.lnHM.untr.rle, 
      type='l', col=col_vector[13])
lines(DE.results.blood_10_DE$mean.discoveries$mean.lnHM.log.rle, 
      DE.results.blood_10_DE$mean.fdr$mean.lnHM.log.rle, 
      type='l', col=col_vector[14])
abline(h=0.05, col='lightgrey')
plot(DE.results.blood_20_DE$mean.discoveries$edgeR.ql.rle, 
     DE.results.blood_20_DE$mean.fdr$edgeR.ql.rle, 
     type='l', ylim=c(0,0.7), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("20 samples per group, RLE\n", 
                 "Differences in mean only"))
lines(DE.results.blood_20_DE$mean.discoveries$edgeR.lr.rle, 
      DE.results.blood_20_DE$mean.fdr$edgeR.lr.rle, 
      type='l', col=col_vector[2])
lines(DE.results.blood_20_DE$mean.discoveries$edgeR.et.rle, 
      DE.results.blood_20_DE$mean.fdr$edgeR.et.rle, 
      type='l', col=col_vector[3])
lines(DE.results.blood_20_DE$mean.discoveries$DESeq2.noif.rle, 
      DE.results.blood_20_DE$mean.fdr$DESeq2.noif.rle, 
      type='l', col=col_vector[4])
lines(DE.results.blood_20_DE$mean.discoveries$DESeq2.if.rle, 
      DE.results.blood_20_DE$mean.fdr$DESeq2.if.rle, 
      type='l', col=col_vector[5])
lines(DE.results.blood_20_DE$mean.discoveries$voom.rle, 
      DE.results.blood_20_DE$mean.fdr$voom.rle, 
      type='l', col=col_vector[6])
lines(DE.results.blood_20_DE$mean.discoveries$DSS.rle, 
      DE.results.blood_20_DE$mean.fdr$DSS.rle, 
      type='l', col=col_vector[7])
lines(DE.results.blood_20_DE$mean.discoveries$baySeq.rle, 
      DE.results.blood_20_DE$mean.fdr$baySeq.rle, 
      type='l', col=col_vector[8])
lines(DE.results.blood_20_DE$mean.discoveries$mean.MDSeq.zi.rle, 
      DE.results.blood_20_DE$mean.fdr$mean.MDSeq.zi.rle, 
      type='l', col=col_vector[9])
lines(DE.results.blood_20_DE$mean.discoveries$mean.MDSeq.nozi.rle, 
      DE.results.blood_20_DE$mean.fdr$mean.MDSeq.nozi.rle, 
      type='l', col=col_vector[10])
lines(DE.results.blood_20_DE$mean.discoveries$mean.expHM.untr.rle, 
      DE.results.blood_20_DE$mean.fdr$mean.expHM.untr.rle, 
      type='l', col=col_vector[11])
lines(DE.results.blood_20_DE$mean.discoveries$mean.expHM.log.rle, 
      DE.results.blood_20_DE$mean.fdr$mean.expHM.log.rle, 
      type='l', col=col_vector[12])
lines(DE.results.blood_20_DE$mean.discoveries$mean.lnHM.untr.rle, 
      DE.results.blood_20_DE$mean.fdr$mean.lnHM.untr.rle, 
      type='l', col=col_vector[13])
lines(DE.results.blood_20_DE$mean.discoveries$mean.lnHM.log.rle, 
      DE.results.blood_20_DE$mean.fdr$mean.lnHM.log.rle, 
      type='l', col=col_vector[14])
abline(h=0.05, col='lightgrey')
plot(DE.results.blood_50_DE$mean.discoveries$edgeR.ql.rle, 
     DE.results.blood_50_DE$mean.fdr$edgeR.ql.rle, 
     type='l', ylim=c(0,0.7), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("50 samples per group, RLE\n", 
                 "Differences in mean only"))
lines(DE.results.blood_50_DE$mean.discoveries$edgeR.lr.rle, 
      DE.results.blood_50_DE$mean.fdr$edgeR.lr.rle, 
      type='l', col=col_vector[2])
lines(DE.results.blood_50_DE$mean.discoveries$edgeR.et.rle, 
      DE.results.blood_50_DE$mean.fdr$edgeR.et.rle, 
      type='l', col=col_vector[3])
lines(DE.results.blood_50_DE$mean.discoveries$DESeq2.noif.rle, 
      DE.results.blood_50_DE$mean.fdr$DESeq2.noif.rle, 
      type='l', col=col_vector[4])
lines(DE.results.blood_50_DE$mean.discoveries$DESeq2.if.rle, 
      DE.results.blood_50_DE$mean.fdr$DESeq2.if.rle, 
      type='l', col=col_vector[5])
lines(DE.results.blood_50_DE$mean.discoveries$voom.rle, 
      DE.results.blood_50_DE$mean.fdr$voom.rle, 
      type='l', col=col_vector[6])
lines(DE.results.blood_50_DE$mean.discoveries$DSS.rle, 
      DE.results.blood_50_DE$mean.fdr$DSS.rle, 
      type='l', col=col_vector[7])
lines(DE.results.blood_50_DE$mean.discoveries$baySeq.rle, 
      DE.results.blood_50_DE$mean.fdr$baySeq.rle, 
      type='l', col=col_vector[8])
lines(DE.results.blood_50_DE$mean.discoveries$mean.MDSeq.zi.rle, 
      DE.results.blood_50_DE$mean.fdr$mean.MDSeq.zi.rle, 
      type='l', col=col_vector[9])
lines(DE.results.blood_50_DE$mean.discoveries$mean.MDSeq.nozi.rle, 
      DE.results.blood_50_DE$mean.fdr$mean.MDSeq.nozi.rle, 
      type='l', col=col_vector[10])
lines(DE.results.blood_50_DE$mean.discoveries$mean.expHM.untr.rle, 
      DE.results.blood_50_DE$mean.fdr$mean.expHM.untr.rle, 
      type='l', col=col_vector[11])
lines(DE.results.blood_50_DE$mean.discoveries$mean.expHM.log.rle, 
      DE.results.blood_50_DE$mean.fdr$mean.expHM.log.rle, 
      type='l', col=col_vector[12])
lines(DE.results.blood_50_DE$mean.discoveries$mean.lnHM.untr.rle, 
      DE.results.blood_50_DE$mean.fdr$mean.lnHM.untr.rle, 
      type='l', col=col_vector[13])
lines(DE.results.blood_50_DE$mean.discoveries$mean.lnHM.log.rle, 
      DE.results.blood_50_DE$mean.fdr$mean.lnHM.log.rle, 
      type='l', col=col_vector[14])
legend("topleft", bty='n', col=col_vector[1:14], lty=1, ncol=2, lwd=2, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                "voom", "DSS", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                "expHM", "expHM log", "lnHM", "lnHM log"))
abline(h=0.05, col='lightgrey')

plot(DE.results.blood_2_DEDD$mean.discoveries$edgeR.ql.tmm, 
     DE.results.blood_2_DEDD$mean.fdr$edgeR.ql.tmm, 
     type='l', ylim=c(0,0.7), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("2 samples per group, TMM\n", 
                 "Differences in mean and dispersion"))
lines(DE.results.blood_2_DEDD$mean.discoveries$edgeR.lr.tmm, 
      DE.results.blood_2_DEDD$mean.fdr$edgeR.lr.tmm, 
      type='l', col=col_vector[2])
lines(DE.results.blood_2_DEDD$mean.discoveries$edgeR.et.tmm, 
      DE.results.blood_2_DEDD$mean.fdr$edgeR.et.tmm, 
      type='l', col=col_vector[3])
lines(DE.results.blood_2_DEDD$mean.discoveries$DESeq2.noif.tmm, 
      DE.results.blood_2_DEDD$mean.fdr$DESeq2.noif.tmm, 
      type='l', col=col_vector[4])
lines(DE.results.blood_2_DEDD$mean.discoveries$DESeq2.if.tmm, 
      DE.results.blood_2_DEDD$mean.fdr$DESeq2.if.tmm, 
      type='l', col=col_vector[5])
lines(DE.results.blood_2_DEDD$mean.discoveries$voom.tmm, 
      DE.results.blood_2_DEDD$mean.fdr$voom.tmm, 
      type='l', col=col_vector[6])
lines(DE.results.blood_2_DEDD$mean.discoveries$DSS.tmm, 
      DE.results.blood_2_DEDD$mean.fdr$DSS.tmm, 
      type='l', col=col_vector[7])
lines(DE.results.blood_2_DEDD$mean.discoveries$baySeq.tmm, 
      DE.results.blood_2_DEDD$mean.fdr$baySeq.tmm, 
      type='l', col=col_vector[8])
lines(DE.results.blood_2_DEDD$mean.discoveries$mean.MDSeq.zi.tmm, 
      DE.results.blood_2_DEDD$mean.fdr$mean.MDSeq.zi.tmm, 
      type='l', col=col_vector[9])
lines(DE.results.blood_2_DEDD$mean.discoveries$mean.MDSeq.nozi.tmm, 
      DE.results.blood_2_DEDD$mean.fdr$mean.MDSeq.nozi.tmm, 
      type='l', col=col_vector[10])
lines(DE.results.blood_2_DEDD$mean.discoveries$mean.expHM.untr.tmm, 
      DE.results.blood_2_DEDD$mean.fdr$mean.expHM.untr.tmm, 
      type='l', col=col_vector[11])
lines(DE.results.blood_2_DEDD$mean.discoveries$mean.expHM.log.tmm, 
      DE.results.blood_2_DEDD$mean.fdr$mean.expHM.log.tmm, 
      type='l', col=col_vector[12])
lines(DE.results.blood_2_DEDD$mean.discoveries$mean.lnHM.untr.tmm, 
      DE.results.blood_2_DEDD$mean.fdr$mean.lnHM.untr.tmm, 
      type='l', col=col_vector[13])
lines(DE.results.blood_2_DEDD$mean.discoveries$mean.lnHM.log.tmm, 
      DE.results.blood_2_DEDD$mean.fdr$mean.lnHM.log.tmm, 
      type='l', col=col_vector[14])
abline(h=0.05, col='lightgrey')
plot(DE.results.blood_5_DEDD$mean.discoveries$edgeR.ql.tmm, 
     DE.results.blood_5_DEDD$mean.fdr$edgeR.ql.tmm, 
     type='l', ylim=c(0,0.7), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("5 samples per group, TMM\n", 
                 "Differences in mean and dispersion"))
lines(DE.results.blood_5_DEDD$mean.discoveries$edgeR.lr.tmm, 
      DE.results.blood_5_DEDD$mean.fdr$edgeR.lr.tmm, 
      type='l', col=col_vector[2])
lines(DE.results.blood_5_DEDD$mean.discoveries$edgeR.et.tmm, 
      DE.results.blood_5_DEDD$mean.fdr$edgeR.et.tmm, 
      type='l', col=col_vector[3])
lines(DE.results.blood_5_DEDD$mean.discoveries$DESeq2.noif.tmm, 
      DE.results.blood_5_DEDD$mean.fdr$DESeq2.noif.tmm, 
      type='l', col=col_vector[4])
lines(DE.results.blood_5_DEDD$mean.discoveries$DESeq2.if.tmm, 
      DE.results.blood_5_DEDD$mean.fdr$DESeq2.if.tmm, 
      type='l', col=col_vector[5])
lines(DE.results.blood_5_DEDD$mean.discoveries$voom.tmm, 
      DE.results.blood_5_DEDD$mean.fdr$voom.tmm, 
      type='l', col=col_vector[6])
lines(DE.results.blood_5_DEDD$mean.discoveries$DSS.tmm, 
      DE.results.blood_5_DEDD$mean.fdr$DSS.tmm, 
      type='l', col=col_vector[7])
lines(DE.results.blood_5_DEDD$mean.discoveries$baySeq.tmm, 
      DE.results.blood_5_DEDD$mean.fdr$baySeq.tmm, 
      type='l', col=col_vector[8])
lines(DE.results.blood_5_DEDD$mean.discoveries$mean.MDSeq.zi.tmm, 
      DE.results.blood_5_DEDD$mean.fdr$mean.MDSeq.zi.tmm, 
      type='l', col=col_vector[9])
lines(DE.results.blood_5_DEDD$mean.discoveries$mean.MDSeq.nozi.tmm, 
      DE.results.blood_5_DEDD$mean.fdr$mean.MDSeq.nozi.tmm, 
      type='l', col=col_vector[10])
lines(DE.results.blood_5_DEDD$mean.discoveries$mean.expHM.untr.tmm, 
      DE.results.blood_5_DEDD$mean.fdr$mean.expHM.untr.tmm, 
      type='l', col=col_vector[11])
lines(DE.results.blood_5_DEDD$mean.discoveries$mean.expHM.log.tmm, 
      DE.results.blood_5_DEDD$mean.fdr$mean.expHM.log.tmm, 
      type='l', col=col_vector[12])
lines(DE.results.blood_5_DEDD$mean.discoveries$mean.lnHM.untr.tmm, 
      DE.results.blood_5_DEDD$mean.fdr$mean.lnHM.untr.tmm, 
      type='l', col=col_vector[13])
lines(DE.results.blood_5_DEDD$mean.discoveries$mean.lnHM.log.tmm, 
      DE.results.blood_5_DEDD$mean.fdr$mean.lnHM.log.tmm, 
      type='l', col=col_vector[14])
abline(h=0.05, col='lightgrey')
plot(DE.results.blood_10_DEDD$mean.discoveries$edgeR.ql.tmm, 
     DE.results.blood_10_DEDD$mean.fdr$edgeR.ql.tmm, 
     type='l', ylim=c(0,0.7), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("10 samples per group, TMM\n", 
                 "Differences in mean and dispersion"))
lines(DE.results.blood_10_DEDD$mean.discoveries$edgeR.lr.tmm, 
      DE.results.blood_10_DEDD$mean.fdr$edgeR.lr.tmm, 
      type='l', col=col_vector[2])
lines(DE.results.blood_10_DEDD$mean.discoveries$edgeR.et.tmm, 
      DE.results.blood_10_DEDD$mean.fdr$edgeR.et.tmm, 
      type='l', col=col_vector[3])
lines(DE.results.blood_10_DEDD$mean.discoveries$DESeq2.noif.tmm, 
      DE.results.blood_10_DEDD$mean.fdr$DESeq2.noif.tmm, 
      type='l', col=col_vector[4])
lines(DE.results.blood_10_DEDD$mean.discoveries$DESeq2.if.tmm, 
      DE.results.blood_10_DEDD$mean.fdr$DESeq2.if.tmm, 
      type='l', col=col_vector[5])
lines(DE.results.blood_10_DEDD$mean.discoveries$voom.tmm, 
      DE.results.blood_10_DEDD$mean.fdr$voom.tmm, 
      type='l', col=col_vector[6])
lines(DE.results.blood_10_DEDD$mean.discoveries$DSS.tmm, 
      DE.results.blood_10_DEDD$mean.fdr$DSS.tmm, 
      type='l', col=col_vector[7])
lines(DE.results.blood_10_DEDD$mean.discoveries$baySeq.tmm, 
      DE.results.blood_10_DEDD$mean.fdr$baySeq.tmm, 
      type='l', col=col_vector[8])
lines(DE.results.blood_10_DEDD$mean.discoveries$mean.MDSeq.zi.tmm, 
      DE.results.blood_10_DEDD$mean.fdr$mean.MDSeq.zi.tmm, 
      type='l', col=col_vector[9])
lines(DE.results.blood_10_DEDD$mean.discoveries$mean.MDSeq.nozi.tmm, 
      DE.results.blood_10_DEDD$mean.fdr$mean.MDSeq.nozi.tmm, 
      type='l', col=col_vector[10])
lines(DE.results.blood_10_DEDD$mean.discoveries$mean.expHM.untr.tmm, 
      DE.results.blood_10_DEDD$mean.fdr$mean.expHM.untr.tmm, 
      type='l', col=col_vector[11])
lines(DE.results.blood_10_DEDD$mean.discoveries$mean.expHM.log.tmm, 
      DE.results.blood_10_DEDD$mean.fdr$mean.expHM.log.tmm, 
      type='l', col=col_vector[12])
lines(DE.results.blood_10_DEDD$mean.discoveries$mean.lnHM.untr.tmm, 
      DE.results.blood_10_DEDD$mean.fdr$mean.lnHM.untr.tmm, 
      type='l', col=col_vector[13])
lines(DE.results.blood_10_DEDD$mean.discoveries$mean.lnHM.log.tmm, 
      DE.results.blood_10_DEDD$mean.fdr$mean.lnHM.log.tmm, 
      type='l', col=col_vector[14])
abline(h=0.05, col='lightgrey')
plot(DE.results.blood_20_DEDD$mean.discoveries$edgeR.ql.tmm, 
     DE.results.blood_20_DEDD$mean.fdr$edgeR.ql.tmm, 
     type='l', ylim=c(0,0.7), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("20 samples per group, TMM\n", 
                 "Differences in mean and dispersion"))
lines(DE.results.blood_20_DEDD$mean.discoveries$edgeR.lr.tmm, 
      DE.results.blood_20_DEDD$mean.fdr$edgeR.lr.tmm, 
      type='l', col=col_vector[2])
lines(DE.results.blood_20_DEDD$mean.discoveries$edgeR.et.tmm, 
      DE.results.blood_20_DEDD$mean.fdr$edgeR.et.tmm, 
      type='l', col=col_vector[3])
lines(DE.results.blood_20_DEDD$mean.discoveries$DESeq2.noif.tmm, 
      DE.results.blood_20_DEDD$mean.fdr$DESeq2.noif.tmm, 
      type='l', col=col_vector[4])
lines(DE.results.blood_20_DEDD$mean.discoveries$DESeq2.if.tmm, 
      DE.results.blood_20_DEDD$mean.fdr$DESeq2.if.tmm, 
      type='l', col=col_vector[5])
lines(DE.results.blood_20_DEDD$mean.discoveries$voom.tmm, 
      DE.results.blood_20_DEDD$mean.fdr$voom.tmm, 
      type='l', col=col_vector[6])
lines(DE.results.blood_20_DEDD$mean.discoveries$DSS.tmm, 
      DE.results.blood_20_DEDD$mean.fdr$DSS.tmm, 
      type='l', col=col_vector[7])
lines(DE.results.blood_20_DEDD$mean.discoveries$baySeq.tmm, 
      DE.results.blood_20_DEDD$mean.fdr$baySeq.tmm, 
      type='l', col=col_vector[8])
lines(DE.results.blood_20_DEDD$mean.discoveries$mean.MDSeq.zi.tmm, 
      DE.results.blood_20_DEDD$mean.fdr$mean.MDSeq.zi.tmm, 
      type='l', col=col_vector[9])
lines(DE.results.blood_20_DEDD$mean.discoveries$mean.MDSeq.nozi.tmm, 
      DE.results.blood_20_DEDD$mean.fdr$mean.MDSeq.nozi.tmm, 
      type='l', col=col_vector[10])
lines(DE.results.blood_20_DEDD$mean.discoveries$mean.expHM.untr.tmm, 
      DE.results.blood_20_DEDD$mean.fdr$mean.expHM.untr.tmm, 
      type='l', col=col_vector[11])
lines(DE.results.blood_20_DEDD$mean.discoveries$mean.expHM.log.tmm, 
      DE.results.blood_20_DEDD$mean.fdr$mean.expHM.log.tmm, 
      type='l', col=col_vector[12])
lines(DE.results.blood_20_DEDD$mean.discoveries$mean.lnHM.untr.tmm, 
      DE.results.blood_20_DEDD$mean.fdr$mean.lnHM.untr.tmm, 
      type='l', col=col_vector[13])
lines(DE.results.blood_20_DEDD$mean.discoveries$mean.lnHM.log.tmm, 
      DE.results.blood_20_DEDD$mean.fdr$mean.lnHM.log.tmm, 
      type='l', col=col_vector[14])
abline(h=0.05, col='lightgrey')
plot(DE.results.blood_50_DEDD$mean.discoveries$edgeR.ql.tmm, 
     DE.results.blood_50_DEDD$mean.fdr$edgeR.ql.tmm, 
     type='l', ylim=c(0,0.7), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("50 samples per group, TMM\n", 
                 "Differences in mean and dispersion"))
lines(DE.results.blood_50_DEDD$mean.discoveries$edgeR.lr.tmm, 
      DE.results.blood_50_DEDD$mean.fdr$edgeR.lr.tmm, 
      type='l', col=col_vector[2])
lines(DE.results.blood_50_DEDD$mean.discoveries$edgeR.et.tmm, 
      DE.results.blood_50_DEDD$mean.fdr$edgeR.et.tmm, 
      type='l', col=col_vector[3])
lines(DE.results.blood_50_DEDD$mean.discoveries$DESeq2.noif.tmm, 
      DE.results.blood_50_DEDD$mean.fdr$DESeq2.noif.tmm, 
      type='l', col=col_vector[4])
lines(DE.results.blood_50_DEDD$mean.discoveries$DESeq2.if.tmm, 
      DE.results.blood_50_DEDD$mean.fdr$DESeq2.if.tmm, 
      type='l', col=col_vector[5])
lines(DE.results.blood_50_DEDD$mean.discoveries$voom.tmm, 
      DE.results.blood_50_DEDD$mean.fdr$voom.tmm, 
      type='l', col=col_vector[6])
lines(DE.results.blood_50_DEDD$mean.discoveries$DSS.tmm, 
      DE.results.blood_50_DEDD$mean.fdr$DSS.tmm, 
      type='l', col=col_vector[7])
lines(DE.results.blood_50_DEDD$mean.discoveries$baySeq.tmm, 
      DE.results.blood_50_DEDD$mean.fdr$baySeq.tmm, 
      type='l', col=col_vector[8])
lines(DE.results.blood_50_DEDD$mean.discoveries$mean.MDSeq.zi.tmm, 
      DE.results.blood_50_DEDD$mean.fdr$mean.MDSeq.zi.tmm, 
      type='l', col=col_vector[9])
lines(DE.results.blood_50_DEDD$mean.discoveries$mean.MDSeq.nozi.tmm, 
      DE.results.blood_50_DEDD$mean.fdr$mean.MDSeq.nozi.tmm, 
      type='l', col=col_vector[10])
lines(DE.results.blood_50_DEDD$mean.discoveries$mean.expHM.untr.tmm, 
      DE.results.blood_50_DEDD$mean.fdr$mean.expHM.untr.tmm, 
      type='l', col=col_vector[11])
lines(DE.results.blood_50_DEDD$mean.discoveries$mean.expHM.log.tmm, 
      DE.results.blood_50_DEDD$mean.fdr$mean.expHM.log.tmm, 
      type='l', col=col_vector[12])
lines(DE.results.blood_50_DEDD$mean.discoveries$mean.lnHM.untr.tmm, 
      DE.results.blood_50_DEDD$mean.fdr$mean.lnHM.untr.tmm, 
      type='l', col=col_vector[13])
lines(DE.results.blood_50_DEDD$mean.discoveries$mean.lnHM.log.tmm, 
      DE.results.blood_50_DEDD$mean.fdr$mean.lnHM.log.tmm, 
      type='l', col=col_vector[14])
abline(h=0.05, col='lightgrey')
legend("topright", bty='n', col=col_vector[1:14], lty=1, ncol=2, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                "voom", "DSS", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                "expHM", "expHM log", "lnHM", "lnHM log"))

plot(DE.results.blood_2_DEDD$mean.discoveries$edgeR.ql.rle, 
     DE.results.blood_2_DEDD$mean.fdr$edgeR.ql.rle, 
     type='l', ylim=c(0,0.7), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("2 samples per group, RLE\n", 
                 "Differences in mean and dispersion"))
lines(DE.results.blood_2_DEDD$mean.discoveries$edgeR.lr.rle, 
      DE.results.blood_2_DEDD$mean.fdr$edgeR.lr.rle, 
      type='l', col=col_vector[2])
lines(DE.results.blood_2_DEDD$mean.discoveries$edgeR.et.rle, 
      DE.results.blood_2_DEDD$mean.fdr$edgeR.et.rle, 
      type='l', col=col_vector[3])
lines(DE.results.blood_2_DEDD$mean.discoveries$DESeq2.noif.rle, 
      DE.results.blood_2_DEDD$mean.fdr$DESeq2.noif.rle, 
      type='l', col=col_vector[4])
lines(DE.results.blood_2_DEDD$mean.discoveries$DESeq2.if.rle, 
      DE.results.blood_2_DEDD$mean.fdr$DESeq2.if.rle, 
      type='l', col=col_vector[5])
lines(DE.results.blood_2_DEDD$mean.discoveries$voom.rle, 
      DE.results.blood_2_DEDD$mean.fdr$voom.rle, 
      type='l', col=col_vector[6])
lines(DE.results.blood_2_DEDD$mean.discoveries$DSS.rle, 
      DE.results.blood_2_DEDD$mean.fdr$DSS.rle, 
      type='l', col=col_vector[7])
lines(DE.results.blood_2_DEDD$mean.discoveries$baySeq.rle, 
      DE.results.blood_2_DEDD$mean.fdr$baySeq.rle, 
      type='l', col=col_vector[8])
lines(DE.results.blood_2_DEDD$mean.discoveries$mean.MDSeq.zi.rle, 
      DE.results.blood_2_DEDD$mean.fdr$mean.MDSeq.zi.rle, 
      type='l', col=col_vector[9])
lines(DE.results.blood_2_DEDD$mean.discoveries$mean.MDSeq.nozi.rle, 
      DE.results.blood_2_DEDD$mean.fdr$mean.MDSeq.nozi.rle, 
      type='l', col=col_vector[10])
lines(DE.results.blood_2_DEDD$mean.discoveries$mean.expHM.untr.rle, 
      DE.results.blood_2_DEDD$mean.fdr$mean.expHM.untr.rle, 
      type='l', col=col_vector[11])
lines(DE.results.blood_2_DEDD$mean.discoveries$mean.expHM.log.rle, 
      DE.results.blood_2_DEDD$mean.fdr$mean.expHM.log.rle, 
      type='l', col=col_vector[12])
lines(DE.results.blood_2_DEDD$mean.discoveries$mean.lnHM.untr.rle, 
      DE.results.blood_2_DEDD$mean.fdr$mean.lnHM.untr.rle, 
      type='l', col=col_vector[13])
lines(DE.results.blood_2_DEDD$mean.discoveries$mean.lnHM.log.rle, 
      DE.results.blood_2_DEDD$mean.fdr$mean.lnHM.log.rle, 
      type='l', col=col_vector[14])
abline(h=0.05, col='lightgrey')
plot(DE.results.blood_5_DEDD$mean.discoveries$edgeR.ql.rle, 
     DE.results.blood_5_DEDD$mean.fdr$edgeR.ql.rle, 
     type='l', ylim=c(0,0.7), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("5 samples per group, RLE\n", 
                 "Differences in mean and dispersion"))
lines(DE.results.blood_5_DEDD$mean.discoveries$edgeR.lr.rle, 
      DE.results.blood_5_DEDD$mean.fdr$edgeR.lr.rle, 
      type='l', col=col_vector[2])
lines(DE.results.blood_5_DEDD$mean.discoveries$edgeR.et.rle, 
      DE.results.blood_5_DEDD$mean.fdr$edgeR.et.rle, 
      type='l', col=col_vector[3])
lines(DE.results.blood_5_DEDD$mean.discoveries$DESeq2.noif.rle, 
      DE.results.blood_5_DEDD$mean.fdr$DESeq2.noif.rle, 
      type='l', col=col_vector[4])
lines(DE.results.blood_5_DEDD$mean.discoveries$DESeq2.if.rle, 
      DE.results.blood_5_DEDD$mean.fdr$DESeq2.if.rle, 
      type='l', col=col_vector[5])
lines(DE.results.blood_5_DEDD$mean.discoveries$voom.rle, 
      DE.results.blood_5_DEDD$mean.fdr$voom.rle, 
      type='l', col=col_vector[6])
lines(DE.results.blood_5_DEDD$mean.discoveries$DSS.rle, 
      DE.results.blood_5_DEDD$mean.fdr$DSS.rle, 
      type='l', col=col_vector[7])
lines(DE.results.blood_5_DEDD$mean.discoveries$baySeq.rle, 
      DE.results.blood_5_DEDD$mean.fdr$baySeq.rle, 
      type='l', col=col_vector[8])
lines(DE.results.blood_5_DEDD$mean.discoveries$mean.MDSeq.zi.rle, 
      DE.results.blood_5_DEDD$mean.fdr$mean.MDSeq.zi.rle, 
      type='l', col=col_vector[9])
lines(DE.results.blood_5_DEDD$mean.discoveries$mean.MDSeq.nozi.rle, 
      DE.results.blood_5_DEDD$mean.fdr$mean.MDSeq.nozi.rle, 
      type='l', col=col_vector[10])
lines(DE.results.blood_5_DEDD$mean.discoveries$mean.expHM.untr.rle, 
      DE.results.blood_5_DEDD$mean.fdr$mean.expHM.untr.rle, 
      type='l', col=col_vector[11])
lines(DE.results.blood_5_DEDD$mean.discoveries$mean.expHM.log.rle, 
      DE.results.blood_5_DEDD$mean.fdr$mean.expHM.log.rle, 
      type='l', col=col_vector[12])
lines(DE.results.blood_5_DEDD$mean.discoveries$mean.lnHM.untr.rle, 
      DE.results.blood_5_DEDD$mean.fdr$mean.lnHM.untr.rle, 
      type='l', col=col_vector[13])
lines(DE.results.blood_5_DEDD$mean.discoveries$mean.lnHM.log.rle, 
      DE.results.blood_5_DEDD$mean.fdr$mean.lnHM.log.rle, 
      type='l', col=col_vector[14])
abline(h=0.05, col='lightgrey')
plot(DE.results.blood_10_DEDD$mean.discoveries$edgeR.ql.rle, 
     DE.results.blood_10_DEDD$mean.fdr$edgeR.ql.rle, 
     type='l', ylim=c(0,0.7), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("10 samples per group, RLE\n", 
                 "Differences in mean and dispersion"))
lines(DE.results.blood_10_DEDD$mean.discoveries$edgeR.lr.rle, 
      DE.results.blood_10_DEDD$mean.fdr$edgeR.lr.rle, 
      type='l', col=col_vector[2])
lines(DE.results.blood_10_DEDD$mean.discoveries$edgeR.et.rle, 
      DE.results.blood_10_DEDD$mean.fdr$edgeR.et.rle, 
      type='l', col=col_vector[3])
lines(DE.results.blood_10_DEDD$mean.discoveries$DESeq2.noif.rle, 
      DE.results.blood_10_DEDD$mean.fdr$DESeq2.noif.rle, 
      type='l', col=col_vector[4])
lines(DE.results.blood_10_DEDD$mean.discoveries$DESeq2.if.rle, 
      DE.results.blood_10_DEDD$mean.fdr$DESeq2.if.rle, 
      type='l', col=col_vector[5])
lines(DE.results.blood_10_DEDD$mean.discoveries$voom.rle, 
      DE.results.blood_10_DEDD$mean.fdr$voom.rle, 
      type='l', col=col_vector[6])
lines(DE.results.blood_10_DEDD$mean.discoveries$DSS.rle, 
      DE.results.blood_10_DEDD$mean.fdr$DSS.rle, 
      type='l', col=col_vector[7])
lines(DE.results.blood_10_DEDD$mean.discoveries$baySeq.rle, 
      DE.results.blood_10_DEDD$mean.fdr$baySeq.rle, 
      type='l', col=col_vector[8])
lines(DE.results.blood_10_DEDD$mean.discoveries$mean.MDSeq.zi.rle, 
      DE.results.blood_10_DEDD$mean.fdr$mean.MDSeq.zi.rle, 
      type='l', col=col_vector[9])
lines(DE.results.blood_10_DEDD$mean.discoveries$mean.MDSeq.nozi.rle, 
      DE.results.blood_10_DEDD$mean.fdr$mean.MDSeq.nozi.rle, 
      type='l', col=col_vector[10])
lines(DE.results.blood_10_DEDD$mean.discoveries$mean.expHM.untr.rle, 
      DE.results.blood_10_DEDD$mean.fdr$mean.expHM.untr.rle, 
      type='l', col=col_vector[11])
lines(DE.results.blood_10_DEDD$mean.discoveries$mean.expHM.log.rle, 
      DE.results.blood_10_DEDD$mean.fdr$mean.expHM.log.rle, 
      type='l', col=col_vector[12])
lines(DE.results.blood_10_DEDD$mean.discoveries$mean.lnHM.untr.rle, 
      DE.results.blood_10_DEDD$mean.fdr$mean.lnHM.untr.rle, 
      type='l', col=col_vector[13])
lines(DE.results.blood_10_DEDD$mean.discoveries$mean.lnHM.log.rle, 
      DE.results.blood_10_DEDD$mean.fdr$mean.lnHM.log.rle, 
      type='l', col=col_vector[14])
abline(h=0.05, col='lightgrey')
plot(DE.results.blood_20_DEDD$mean.discoveries$edgeR.ql.rle, 
     DE.results.blood_20_DEDD$mean.fdr$edgeR.ql.rle, 
     type='l', ylim=c(0,0.7), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("20 samples per group, RLE\n", 
                 "Differences in mean and dispersion"))
lines(DE.results.blood_20_DEDD$mean.discoveries$edgeR.lr.rle, 
      DE.results.blood_20_DEDD$mean.fdr$edgeR.lr.rle, 
      type='l', col=col_vector[2])
lines(DE.results.blood_20_DEDD$mean.discoveries$edgeR.et.rle, 
      DE.results.blood_20_DEDD$mean.fdr$edgeR.et.rle, 
      type='l', col=col_vector[3])
lines(DE.results.blood_20_DEDD$mean.discoveries$DESeq2.noif.rle, 
      DE.results.blood_20_DEDD$mean.fdr$DESeq2.noif.rle, 
      type='l', col=col_vector[4])
lines(DE.results.blood_20_DEDD$mean.discoveries$DESeq2.if.rle, 
      DE.results.blood_20_DEDD$mean.fdr$DESeq2.if.rle, 
      type='l', col=col_vector[5])
lines(DE.results.blood_20_DEDD$mean.discoveries$voom.rle, 
      DE.results.blood_20_DEDD$mean.fdr$voom.rle, 
      type='l', col=col_vector[6])
lines(DE.results.blood_20_DEDD$mean.discoveries$DSS.rle, 
      DE.results.blood_20_DEDD$mean.fdr$DSS.rle, 
      type='l', col=col_vector[7])
lines(DE.results.blood_20_DEDD$mean.discoveries$baySeq.rle, 
      DE.results.blood_20_DEDD$mean.fdr$baySeq.rle, 
      type='l', col=col_vector[8])
lines(DE.results.blood_20_DEDD$mean.discoveries$mean.MDSeq.zi.rle, 
      DE.results.blood_20_DEDD$mean.fdr$mean.MDSeq.zi.rle, 
      type='l', col=col_vector[9])
lines(DE.results.blood_20_DEDD$mean.discoveries$mean.MDSeq.nozi.rle, 
      DE.results.blood_20_DEDD$mean.fdr$mean.MDSeq.nozi.rle, 
      type='l', col=col_vector[10])
lines(DE.results.blood_20_DEDD$mean.discoveries$mean.expHM.untr.rle, 
      DE.results.blood_20_DEDD$mean.fdr$mean.expHM.untr.rle, 
      type='l', col=col_vector[11])
lines(DE.results.blood_20_DEDD$mean.discoveries$mean.expHM.log.rle, 
      DE.results.blood_20_DEDD$mean.fdr$mean.expHM.log.rle, 
      type='l', col=col_vector[12])
lines(DE.results.blood_20_DEDD$mean.discoveries$mean.lnHM.untr.rle, 
      DE.results.blood_20_DEDD$mean.fdr$mean.lnHM.untr.rle, 
      type='l', col=col_vector[13])
lines(DE.results.blood_20_DEDD$mean.discoveries$mean.lnHM.log.rle, 
      DE.results.blood_20_DEDD$mean.fdr$mean.lnHM.log.rle, 
      type='l', col=col_vector[14])
abline(h=0.05, col='lightgrey')
plot(DE.results.blood_50_DEDD$mean.discoveries$edgeR.ql.rle, 
     DE.results.blood_50_DEDD$mean.fdr$edgeR.ql.rle, 
     type='l', ylim=c(0,0.7), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("50 samples per group, RLE\n", 
                 "Differences in mean and dispersion"))
lines(DE.results.blood_50_DEDD$mean.discoveries$edgeR.lr.rle, 
      DE.results.blood_50_DEDD$mean.fdr$edgeR.lr.rle, 
      type='l', col=col_vector[2])
lines(DE.results.blood_50_DEDD$mean.discoveries$edgeR.et.rle, 
      DE.results.blood_50_DEDD$mean.fdr$edgeR.et.rle, 
      type='l', col=col_vector[3])
lines(DE.results.blood_50_DEDD$mean.discoveries$DESeq2.noif.rle, 
      DE.results.blood_50_DEDD$mean.fdr$DESeq2.noif.rle, 
      type='l', col=col_vector[4])
lines(DE.results.blood_50_DEDD$mean.discoveries$DESeq2.if.rle, 
      DE.results.blood_50_DEDD$mean.fdr$DESeq2.if.rle, 
      type='l', col=col_vector[5])
lines(DE.results.blood_50_DEDD$mean.discoveries$voom.rle, 
      DE.results.blood_50_DEDD$mean.fdr$voom.rle, 
      type='l', col=col_vector[6])
lines(DE.results.blood_50_DEDD$mean.discoveries$DSS.rle, 
      DE.results.blood_50_DEDD$mean.fdr$DSS.rle, 
      type='l', col=col_vector[7])
lines(DE.results.blood_50_DEDD$mean.discoveries$baySeq.rle, 
      DE.results.blood_50_DEDD$mean.fdr$baySeq.rle, 
      type='l', col=col_vector[8])
lines(DE.results.blood_50_DEDD$mean.discoveries$mean.MDSeq.zi.rle, 
      DE.results.blood_50_DEDD$mean.fdr$mean.MDSeq.zi.rle, 
      type='l', col=col_vector[9])
lines(DE.results.blood_50_DEDD$mean.discoveries$mean.MDSeq.nozi.rle, 
      DE.results.blood_50_DEDD$mean.fdr$mean.MDSeq.nozi.rle, 
      type='l', col=col_vector[10])
lines(DE.results.blood_50_DEDD$mean.discoveries$mean.expHM.untr.rle, 
      DE.results.blood_50_DEDD$mean.fdr$mean.expHM.untr.rle, 
      type='l', col=col_vector[11])
lines(DE.results.blood_50_DEDD$mean.discoveries$mean.expHM.log.rle, 
      DE.results.blood_50_DEDD$mean.fdr$mean.expHM.log.rle, 
      type='l', col=col_vector[12])
lines(DE.results.blood_50_DEDD$mean.discoveries$mean.lnHM.untr.rle, 
      DE.results.blood_50_DEDD$mean.fdr$mean.lnHM.untr.rle, 
      type='l', col=col_vector[13])
lines(DE.results.blood_50_DEDD$mean.discoveries$mean.lnHM.log.rle, 
      DE.results.blood_50_DEDD$mean.fdr$mean.lnHM.log.rle, 
      type='l', col=col_vector[14])
abline(h=0.05, col='lightgrey')
legend("topright", bty='n', col=col_vector[1:14], lty=1, ncol=2, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                "voom", "DSS", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                "expHM", "expHM log", "lnHM", "lnHM log"))



###################################
#### Differential distribution ####
###################################

#############
#### AUC ####
par(mfrow=c(2,3), mar=c(1,2.5,5.5,1), mgp=c(3,0.7,0))
boxplot(DEDD.results.blood_2_DD$auc, names=NA, col=col_vector[1:3], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0.35,1.05), 
        main=paste0("AUC diff dist, 2 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in dispersion only"), 
        cex.main=1.7)
lines(c(3.5,3.5), c(0.35,0.97), col="lightgrey")
legend("topleft", fill=col_vector[1:3], bty='n', cex=1.5, ncol=3, 
       legend=c("diffVar", "expHMM", "lnHMM"))
boxplot(DEDD.results.blood_5_DD$auc, names=NA, col=col_vector[1:3], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0.35,1.1), 
        main=paste0("AUC diff dist, 5 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in dispersion only"), 
        cex.main=1.7)
lines(c(3.5,3.5), c(0.35,0.97), col="lightgrey")
legend("topleft", fill=col_vector[1:3], bty='n', cex=1.5, ncol=3, 
       legend=c("diffVar", "expHMM", "lnHMM"))
boxplot(DEDD.results.blood_10_DD$auc, names=NA, col=col_vector[1:3], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0.35,1.1), 
        main=paste0("AUC diff dist, 10 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in dispersion only"), 
        cex.main=1.7)
lines(c(3.5,3.5), c(0.35,0.97), col="lightgrey")
legend("topleft", fill=col_vector[1:3], bty='n', cex=1.5, ncol=3, 
       legend=c("diffVar", "expHMM", "lnHMM"))
boxplot(DEDD.results.blood_20_DD$auc, names=NA, col=col_vector[1:3], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0.35,1.1), 
        main=paste0("AUC diff dist, 20 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in dispersion only"), 
        cex.main=1.7)
lines(c(3.5,3.5), c(0.35,0.97), col="lightgrey")
legend("topleft", fill=col_vector[1:3], bty='n', cex=1.5, ncol=3, 
       legend=c("diffVar", "expHMM", "lnHMM"))
boxplot(DEDD.results.blood_50_DD$auc, names=NA, col=col_vector[1:3], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0.35,1.1), 
        main=paste0("AUC diff dist, 50 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in dispersion only"), 
        cex.main=1.7)
lines(c(3.5,3.5), c(0.35,0.97), col="lightgrey")
legend("topleft", fill=col_vector[1:3], bty='n', cex=1.5, ncol=3, 
       legend=c("diffVar", "expHMM", "lnHMM"))
rbind(colMeans(DEDD.results.blood_2_DD$auc), 
      colMeans(DEDD.results.blood_5_DD$auc), 
      colMeans(DEDD.results.blood_10_DD$auc), 
      colMeans(DEDD.results.blood_20_DD$auc), 
      colMeans(DEDD.results.blood_50_DD$auc))

par(mfrow=c(2,3), mar=c(1,2.5,5.5,1), mgp=c(3,0.7,0))
boxplot(DEDD.results.blood_2_DE$auc, names=NA, col=col_vector[1:2], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0.35,1.05), 
        main=paste0("AUC diff dist, 2 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "differences in mean only"), 
        cex.main=1.7)
lines(c(2.5,2.5), c(0.35,0.97), col="lightgrey")
legend("topleft", fill=col_vector[1:2], bty='n',cex=1.5, ncol=3, 
       legend=c("expHMM", "lnHMM"))
boxplot(DEDD.results.blood_5_DE$auc, names=NA, col=col_vector[1:2], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0.35,1.05), 
        main=paste0("AUC diff dist, 5 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "differences in mean only"), 
        cex.main=1.7)
lines(c(2.5,2.5), c(0.35,0.97), col="lightgrey")
legend("topleft", fill=col_vector[1:2], bty='n',cex=1.5, ncol=3, 
       legend=c("expHMM", "lnHMM"))
boxplot(DEDD.results.blood_10_DE$auc, names=NA, col=col_vector[1:2], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0.35,1.05), 
        main=paste0("AUC diff dist, 10 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "differences in mean only"), 
        cex.main=1.7)
lines(c(2.5,2.5), c(0.35,0.97), col="lightgrey")
legend("topleft", fill=col_vector[1:2], bty='n',cex=1.5, ncol=3, 
       legend=c("expHMM", "lnHMM"))
boxplot(DEDD.results.blood_20_DE$auc, names=NA, col=col_vector[1:2], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0.35,1.05), 
        main=paste0("AUC diff dist, 20 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "differences in mean only"), 
        cex.main=1.7)
lines(c(2.5,2.5), c(0.35,0.97), col="lightgrey")
legend("topleft", fill=col_vector[1:2], bty='n',cex=1.5, ncol=3, 
       legend=c("expHMM", "lnHMM"))
boxplot(DEDD.results.blood_50_DE$auc, names=NA, col=col_vector[1:2], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0.35,1.05), 
        main=paste0("AUC diff dist, 50 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "differences in mean only"), 
        cex.main=1.7)
lines(c(2.5,2.5), c(0.35,0.97), col="lightgrey")
legend("topleft", fill=col_vector[1:2], bty='n',cex=1.5, ncol=3, 
       legend=c("expHMM", "lnHMM"))
rbind(colMeans(DEDD.results.blood_2_DE$auc), 
      colMeans(DEDD.results.blood_5_DE$auc), 
      colMeans(DEDD.results.blood_10_DE$auc), 
      colMeans(DEDD.results.blood_20_DE$auc), 
      colMeans(DEDD.results.blood_50_DE$auc))

par(mfrow=c(2,3), mar=c(1,2.5,5.5,1), mgp=c(3,0.7,0))
boxplot(DEDD.results.blood_2_DEDD$auc, names=NA, col=col_vector[1:9], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0.35,1.05), 
        main=paste0("AUC diff dist, 2 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "differences in mean and dispersion"), 
        cex.main=1.7)
lines(c(9.5,9.5), c(0.35,0.97), col="lightgrey")
legend("topleft", fill=col_vector[1:9], bty='n', cex=1.2, ncol=3, 
       legend=c("expHMM", "lnHMM", "diffVar", "edgeR_diffVar", "edgeR_MDSeq", 
                "edgeR_lnHM", "lnHM_diffVar", "lnHM_MDSeq", "lnHM_lnHM"))
boxplot(DEDD.results.blood_5_DEDD$auc, names=NA, col=col_vector[1:9], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0.35,1.05), 
        main=paste0("AUC diff dist, 5 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "differences in mean and dispersion"), 
        cex.main=1.7)
lines(c(9.5,9.5), c(0.35,0.97), col="lightgrey")
legend("topleft", fill=col_vector[1:9], bty='n', cex=1.2, ncol=3, 
       legend=c("expHMM", "lnHMM", "diffVar", "edgeR_diffVar", "edgeR_MDSeq", 
                "edgeR_lnHM", "lnHM_diffVar", "lnHM_MDSeq", "lnHM_lnHM"))
boxplot(DEDD.results.blood_10_DEDD$auc, names=NA, col=col_vector[1:9], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0.35,1.05), 
        main=paste0("AUC diff dist, 10 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "differences in mean and dispersion"), 
        cex.main=1.7)
lines(c(9.5,9.5), c(0.35,0.97), col="lightgrey")
legend("topleft", fill=col_vector[1:9], bty='n', cex=1.2, ncol=3, 
       legend=c("expHMM", "lnHMM", "diffVar", "edgeR_diffVar", "edgeR_MDSeq", 
                "edgeR_lnHM", "lnHM_diffVar", "lnHM_MDSeq", "lnHM_lnHM"))
boxplot(DEDD.results.blood_20_DEDD$auc, names=NA, col=col_vector[1:9], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0.35,1.05), 
        main=paste0("AUC diff dist, 20 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "differences in mean and dispersion"), 
        cex.main=1.7)
lines(c(9.5,9.5), c(0.35,0.97), col="lightgrey")
legend("topleft", fill=col_vector[1:9], bty='n', cex=1.2, ncol=3, 
       legend=c("expHMM", "lnHMM", "diffVar", "edgeR_diffVar", "edgeR_MDSeq", 
                "edgeR_lnHM", "lnHM_diffVar", "lnHM_MDSeq", "lnHM_lnHM"))
boxplot(DEDD.results.blood_50_DEDD$auc, names=NA, col=col_vector[1:9], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0.35,1.05), 
        main=paste0("AUC diff dist, 50 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "differences in mean and dispersion"), 
        cex.main=1.7)
lines(c(9.5,9.5), c(0.35,0.97), col="lightgrey")
legend("topleft", fill=col_vector[1:9], bty='n', cex=1.2, ncol=3, 
       legend=c("expHMM", "lnHMM", "diffVar", "edgeR_diffVar", "edgeR_MDSeq", 
                "edgeR_lnHM", "lnHM_diffVar", "lnHM_MDSeq", "lnHM_lnHM"))
rbind(colMeans(DEDD.results.blood_2_DEDD$auc), 
      colMeans(DEDD.results.blood_5_DEDD$auc), 
      colMeans(DEDD.results.blood_10_DEDD$auc), 
      colMeans(DEDD.results.blood_20_DEDD$auc), 
      colMeans(DEDD.results.blood_50_DEDD$auc))


##############
#### pAUC ####
par(mfrow=c(2,3), mar=c(1,2.5,5.5,1), mgp=c(3,0.7,0))
boxplot(DEDD.results.blood_2_DD$pauc, names=NA, col=col_vector[1:3], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,0.055), 
        main=paste0("Partial AUC diff dist 2 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in dispersion only"), 
        cex.main=1.7)
lines(c(3.5,3.5), c(0,0.047), col="lightgrey")
legend("topleft", fill=col_vector[1:3], bty='n', cex=1.5, ncol=3, 
       legend=c("diffVar", "expHMM", "lnHMM"))
boxplot(DEDD.results.blood_5_DD$pauc, names=NA, col=col_vector[1:3], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,0.055), 
        main=paste0("Partial AUC diff dist 5 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in dispersion only"), 
        cex.main=1.7)
lines(c(3.5,3.5), c(0,0.047), col="lightgrey")
legend("topleft", fill=col_vector[1:3], bty='n', cex=1.5, ncol=3, 
       legend=c("diffVar", "expHMM", "lnHMM"))
boxplot(DEDD.results.blood_10_DD$pauc, names=NA, col=col_vector[1:3], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,0.055), 
        main=paste0("Partial AUC diff dist 10 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in dispersion only"), 
        cex.main=1.7)
lines(c(3.5,3.5), c(0,0.047), col="lightgrey")
legend("topleft", fill=col_vector[1:3], bty='n', cex=1.5, ncol=3, 
       legend=c("diffVar", "expHMM", "lnHMM"))
boxplot(DEDD.results.blood_20_DD$pauc, names=NA, col=col_vector[1:3], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,0.055), 
        main=paste0("Partial AUC diff dist 20 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in dispersion only"), 
        cex.main=1.7)
lines(c(3.5,3.5), c(0,0.047), col="lightgrey")
legend("topleft", fill=col_vector[1:3], bty='n', cex=1.5, ncol=3, 
       legend=c("diffVar", "expHMM", "lnHMM"))
boxplot(DEDD.results.blood_50_DD$pauc, names=NA, col=col_vector[1:3], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,0.055), 
        main=paste0("Partial AUC diff dist 50 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in dispersion only"), 
        cex.main=1.7)
lines(c(3.5,3.5), c(0,0.047), col="lightgrey")
legend("topleft", fill=col_vector[1:3], bty='n', cex=1.5, ncol=3, 
       legend=c("diffVar", "expHMM", "lnHMM"))
rbind(colMeans(DEDD.results.blood_2_DD$pauc), 
      colMeans(DEDD.results.blood_5_DD$pauc), 
      colMeans(DEDD.results.blood_10_DD$pauc), 
      colMeans(DEDD.results.blood_20_DD$pauc), 
      colMeans(DEDD.results.blood_50_DD$pauc))

par(mfrow=c(2,3), mar=c(1,2.5,5.5,1), mgp=c(3,0.7,0))
boxplot(DEDD.results.blood_2_DE$pauc, names=NA, col=col_vector[1:2], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,0.055), 
        main=paste0("Partial AUC diff dist 2 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean only"), 
        cex.main=1.7)
lines(c(2.5,2.5), c(0,0.047), col="lightgrey")
legend("topleft", fill=col_vector[1:2], bty='n', cex=1.5, ncol=3, 
       legend=c("expHMM", "lnHMM"))
boxplot(DEDD.results.blood_5_DE$pauc, names=NA, col=col_vector[1:2], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,0.055), 
        main=paste0("Partial AUC diff dist 5 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean only"), 
        cex.main=1.7)
lines(c(2.5,2.5), c(0,0.047), col="lightgrey")
legend("topleft", fill=col_vector[1:2], bty='n', cex=1.5, ncol=3, 
       legend=c("expHMM", "lnHMM"))
boxplot(DEDD.results.blood_10_DE$pauc, names=NA, col=col_vector[1:2], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,0.055), 
        main=paste0("Partial AUC diff dist 10 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean only"), 
        cex.main=1.7)
lines(c(2.5,2.5), c(0,0.047), col="lightgrey")
legend("topleft", fill=col_vector[1:2], bty='n', cex=1.5, ncol=3, 
       legend=c("expHMM", "lnHMM"))
boxplot(DEDD.results.blood_20_DE$pauc, names=NA, col=col_vector[1:2], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,0.055), 
        main=paste0("Partial AUC diff dist 20 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean only"), 
        cex.main=1.7)
lines(c(2.5,2.5), c(0,0.047), col="lightgrey")
legend("topleft", fill=col_vector[1:2], bty='n', cex=1.5, ncol=3, 
       legend=c("expHMM", "lnHMM"))
boxplot(DEDD.results.blood_50_DE$pauc, names=NA, col=col_vector[1:2], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,0.055), 
        main=paste0("Partial AUC diff dist 50 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean only"), 
        cex.main=1.7)
lines(c(2.5,2.5), c(0,0.047), col="lightgrey")
legend("topleft", fill=col_vector[1:2], bty='n', cex=1.5, ncol=3, 
       legend=c("expHMM", "lnHMM"))
rbind(colMeans(DEDD.results.blood_2_DE$pauc), 
      colMeans(DEDD.results.blood_5_DE$pauc), 
      colMeans(DEDD.results.blood_10_DE$pauc), 
      colMeans(DEDD.results.blood_20_DE$pauc), 
      colMeans(DEDD.results.blood_50_DE$pauc))

par(mfrow=c(2,3), mar=c(1,2.5,5.5,1), mgp=c(3,0.7,0))
boxplot(DEDD.results.blood_2_DEDD$pauc, names=NA, col=col_vector[1:9], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,0.055), 
        main=paste0("Partial AUC diff dist 2 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean and dispersion"), 
        cex.main=1.7)
lines(c(9.5,9.5), c(0,0.047), col="lightgrey")
legend("topleft", fill=col_vector[1:9], bty='n', cex=1.2, ncol=3, 
       legend=c("diffVar", "expHMM", "lnHMM", "edgeR_diffVar", "edgeR_MDSeq", 
                "edgeR_lnHM", "lnHM_diffVar", "lnHM_MDSeq", "lnHM_lnHM"))
boxplot(DEDD.results.blood_5_DEDD$pauc, names=NA, col=col_vector[1:9], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,0.055), 
        main=paste0("Partial AUC diff dist 5 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean and dispersion"), 
        cex.main=1.7)
lines(c(9.5,9.5), c(0,0.047), col="lightgrey")
legend("topleft", fill=col_vector[1:9], bty='n', cex=1.2, ncol=3, 
       legend=c("diffVar", "expHMM", "lnHMM", "edgeR_diffVar", "edgeR_MDSeq", 
                "edgeR_lnHM", "lnHM_diffVar", "lnHM_MDSeq", "lnHM_lnHM"))
boxplot(DEDD.results.blood_10_DEDD$pauc, names=NA, col=col_vector[1:9], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,0.055), 
        main=paste0("Partial AUC diff dist 10 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean and dispersion"), 
        cex.main=1.7)
lines(c(9.5,9.5), c(0,0.047), col="lightgrey")
legend("topleft", fill=col_vector[1:9], bty='n', cex=1.2, ncol=3, 
       legend=c("diffVar", "expHMM", "lnHMM", "edgeR_diffVar", "edgeR_MDSeq", 
                "edgeR_lnHM", "lnHM_diffVar", "lnHM_MDSeq", "lnHM_lnHM"))
boxplot(DEDD.results.blood_20_DEDD$pauc, names=NA, col=col_vector[1:9], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,0.055), 
        main=paste0("Partial AUC diff dist 20 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean and dispersion"), 
        cex.main=1.7)
lines(c(9.5,9.5), c(0,0.047), col="lightgrey")
legend("topleft", fill=col_vector[1:9], bty='n', cex=1.2, ncol=3, 
       legend=c("diffVar", "expHMM", "lnHMM", "edgeR_diffVar", "edgeR_MDSeq", 
                "edgeR_lnHM", "lnHM_diffVar", "lnHM_MDSeq", "lnHM_lnHM"))
boxplot(DEDD.results.blood_50_DEDD$pauc, names=NA, col=col_vector[1:9], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,0.055), 
        main=paste0("Partial AUC diff dist 50 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean and dispersion"), 
        cex.main=1.7)
lines(c(9.5,9.5), c(0,0.047), col="lightgrey")
legend("topleft", fill=col_vector[1:9], bty='n', cex=1.2, ncol=3, 
       legend=c("diffVar", "expHMM", "lnHMM", "edgeR_diffVar", "edgeR_MDSeq", 
                "edgeR_lnHM", "lnHM_diffVar", "lnHM_MDSeq", "lnHM_lnHM"))
rbind(colMeans(DEDD.results.blood_2_DEDD$pauc), 
      colMeans(DEDD.results.blood_5_DEDD$pauc), 
      colMeans(DEDD.results.blood_10_DEDD$pauc), 
      colMeans(DEDD.results.blood_20_DEDD$pauc), 
      colMeans(DEDD.results.blood_50_DEDD$pauc))


#############
#### FDR ####
par(mfrow=c(2,3), mar=c(1,2.5,5.5,1), mgp=c(3,0.7,0))
boxplot(DEDD.results.blood_2_DD$fdr, names=NA, col=col_vector[1:7], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.5), 
        main=paste0("FDR diff dist 2 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in dispersion only"), 
        cex.main=1.7)
lines(c(7.5,7.5), c(0,1), col="lightgrey"); abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[1:7], bty='n', cex=1.4, ncol=3, 
       legend=c("diffVar", "expHMM 0.5", "expHMM thr", "expHMM bfdr", 
                "lnHMM 0.5", "lnHMM thr", "lnHMM bfdr"))
boxplot(DEDD.results.blood_5_DD$fdr, names=NA, col=col_vector[1:7], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.5), 
        main=paste0("FDR diff dist 5 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in dispersion only"), 
        cex.main=1.7)
lines(c(7.5,7.5), c(0,1), col="lightgrey"); abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[1:7], bty='n', cex=1.4, ncol=3, 
       legend=c("diffVar", "expHMM 0.5", "expHMM thr", "expHMM bfdr", 
                "lnHMM 0.5", "lnHMM thr", "lnHMM bfdr"))
boxplot(DEDD.results.blood_10_DD$fdr, names=NA, col=col_vector[1:7], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.5), 
        main=paste0("FDR diff dist 10 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in dispersion only"), 
        cex.main=1.7)
lines(c(7.5,7.5), c(0,1), col="lightgrey"); abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[1:7], bty='n', cex=1.4, ncol=3, 
       legend=c("diffVar", "expHMM 0.5", "expHMM thr", "expHMM bfdr", 
                "lnHMM 0.5", "lnHMM thr", "lnHMM bfdr"))
boxplot(DEDD.results.blood_20_DD$fdr, names=NA, col=col_vector[1:7], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.5), 
        main=paste0("FDR diff dist 20 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in dispersion only"), 
        cex.main=1.7)
lines(c(7.5,7.5), c(0,1), col="lightgrey"); abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[1:7], bty='n', cex=1.4, ncol=3, 
       legend=c("diffVar", "expHMM 0.5", "expHMM thr", "expHMM bfdr", 
                "lnHMM 0.5", "lnHMM thr", "lnHMM bfdr"))
boxplot(DEDD.results.blood_50_DD$fdr, names=NA, col=col_vector[1:7], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.5), 
        main=paste0("FDR diff dist 50 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in dispersion only"), 
        cex.main=1.7)
lines(c(7.5,7.5), c(0,1), col="lightgrey"); abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[1:7], bty='n', cex=1.4, ncol=3, 
       legend=c("diffVar", "expHMM 0.5", "expHMM thr", "expHMM bfdr", 
                "lnHMM 0.5", "lnHMM thr", "lnHMM bfdr"))
rbind(colMeans(DEDD.results.blood_2_DD$fdr), 
      colMeans(DEDD.results.blood_5_DD$fdr), 
      colMeans(DEDD.results.blood_10_DD$fdr), 
      colMeans(DEDD.results.blood_20_DD$fdr), 
      colMeans(DEDD.results.blood_50_DD$fdr))

par(mfrow=c(2,3), mar=c(1,2.5,5.5,1), mgp=c(3,0.7,0))
boxplot(DEDD.results.blood_2_DE$fdr, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.5), 
        main=paste0("FDR diff dist 2 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean only"), 
        cex.main=1.7)
lines(c(6.5,6.5), c(0,1), col="lightgrey"); abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.4, ncol=3, 
       legend=c("expHMM 0.5", "expHMM thr", "expHMM bfdr", 
                "lnHMM 0.5", "lnHMM thr", "lnHMM bfdr"))
boxplot(DEDD.results.blood_5_DE$fdr, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.5), 
        main=paste0("FDR diff dist 5 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean only"), 
        cex.main=1.7)
lines(c(6.5,6.5), c(0,1), col="lightgrey"); abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.4, ncol=3, 
       legend=c("expHMM 0.5", "expHMM thr", "expHMM bfdr", 
                "lnHMM 0.5", "lnHMM thr", "lnHMM bfdr"))
boxplot(DEDD.results.blood_10_DE$fdr, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.5), 
        main=paste0("FDR diff dist 10 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean only"), 
        cex.main=1.7)
lines(c(6.5,6.5), c(0,1), col="lightgrey"); abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.4, ncol=3, 
       legend=c("expHMM 0.5", "expHMM thr", "expHMM bfdr", 
                "lnHMM 0.5", "lnHMM thr", "lnHMM bfdr"))
boxplot(DEDD.results.blood_20_DE$fdr, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.5), 
        main=paste0("FDR diff dist 20 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean only"), 
        cex.main=1.7)
lines(c(6.5,6.5), c(0,1), col="lightgrey"); abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.4, ncol=3, 
       legend=c("expHMM 0.5", "expHMM thr", "expHMM bfdr", 
                "lnHMM 0.5", "lnHMM thr", "lnHMM bfdr"))
boxplot(DEDD.results.blood_50_DE$fdr, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.5), 
        main=paste0("FDR diff dist 50 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean only"), 
        cex.main=1.7)
lines(c(6.5,6.5), c(0,1), col="lightgrey"); abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.4, ncol=3, 
       legend=c("expHMM 0.5", "expHMM thr", "expHMM bfdr", 
                "lnHMM 0.5", "lnHMM thr", "lnHMM bfdr"))
rbind(colMeans(DEDD.results.blood_2_DE$fdr), 
      colMeans(DEDD.results.blood_5_DE$fdr), 
      colMeans(DEDD.results.blood_10_DE$fdr), 
      colMeans(DEDD.results.blood_20_DE$fdr), 
      colMeans(DEDD.results.blood_50_DE$fdr))

par(mfrow=c(2,3), mar=c(1,2.5,5.5,1), mgp=c(3,0.7,0))
boxplot(DEDD.results.blood_2_DEDD$fdr, names=NA, col=col_vector[1:13], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.5), 
        main=paste0("FDR diff dist 2 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean and dispersion"), 
        cex.main=1.7)
lines(c(13.5,13.5), c(0,1), col="lightgrey"); abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[1:13], bty='n', cex=1.2, ncol=3, 
       legend=c("diffVar", "expHMM 0.5", "expHMM thr", "expHMM bfdr", 
                "lnHMM 0.5", "lnHMM thr", "lnHMM bfdr", 
                "edgeR_diffVar", "edgeR_MDSeq", "edgeR_lnHM", 
                "lnHM_diffVar", "lnHM_MDSeq", "lnHM_lnHM"))
boxplot(DEDD.results.blood_5_DEDD$fdr, names=NA, col=col_vector[1:13], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.5), 
        main=paste0("FDR diff dist 5 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean and dispersion"), 
        cex.main=1.7)
lines(c(13.5,13.5), c(0,1), col="lightgrey"); abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[1:13], bty='n', cex=1.2, ncol=3, 
       legend=c("diffVar", "expHMM 0.5", "expHMM thr", "expHMM bfdr", 
                "lnHMM 0.5", "lnHMM thr", "lnHMM bfdr", 
                "edgeR_diffVar", "edgeR_MDSeq", "edgeR_lnHM", 
                "lnHM_diffVar", "lnHM_MDSeq", "lnHM_lnHM"))
boxplot(DEDD.results.blood_10_DEDD$fdr, names=NA, col=col_vector[1:13], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.5), 
        main=paste0("FDR diff dist 10 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean and dispersion"), 
        cex.main=1.7)
lines(c(13.5,13.5), c(0,1), col="lightgrey"); abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[1:13], bty='n', cex=1.2, ncol=3, 
       legend=c("diffVar", "expHMM 0.5", "expHMM thr", "expHMM bfdr", 
                "lnHMM 0.5", "lnHMM thr", "lnHMM bfdr", 
                "edgeR_diffVar", "edgeR_MDSeq", "edgeR_lnHM", 
                "lnHM_diffVar", "lnHM_MDSeq", "lnHM_lnHM"))
boxplot(DEDD.results.blood_20_DEDD$fdr, names=NA, col=col_vector[1:13], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.5), 
        main=paste0("FDR diff dist 20 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean and dispersion"), 
        cex.main=1.7)
lines(c(13.5,13.5), c(0,1), col="lightgrey"); abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[1:13], bty='n', cex=1.2, ncol=3, 
       legend=c("diffVar", "expHMM 0.5", "expHMM thr", "expHMM bfdr", 
                "lnHMM 0.5", "lnHMM thr", "lnHMM bfdr", 
                "edgeR_diffVar", "edgeR_MDSeq", "edgeR_lnHM", 
                "lnHM_diffVar", "lnHM_MDSeq", "lnHM_lnHM"))
boxplot(DEDD.results.blood_50_DEDD$fdr, names=NA, col=col_vector[1:13], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.5), 
        main=paste0("FDR diff dist 50 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean and dispersion"), 
        cex.main=1.7)
lines(c(13.5,13.5), c(0,1), col="lightgrey"); abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[1:13], bty='n', cex=1.2, ncol=3, 
       legend=c("diffVar", "expHMM 0.5", "expHMM thr", "expHMM bfdr", 
                "lnHMM 0.5", "lnHMM thr", "lnHMM bfdr", 
                "edgeR_diffVar", "edgeR_MDSeq", "edgeR_lnHM", 
                "lnHM_diffVar", "lnHM_MDSeq", "lnHM_lnHM"))
rbind(colMeans(DEDD.results.blood_2_DEDD$fdr), 
      colMeans(DEDD.results.blood_5_DEDD$fdr), 
      colMeans(DEDD.results.blood_10_DEDD$fdr), 
      colMeans(DEDD.results.blood_20_DEDD$fdr), 
      colMeans(DEDD.results.blood_50_DEDD$fdr))


#############
#### TPR ####
par(mfrow=c(2,3), mar=c(1,2.5,5.5,1), mgp=c(3,0.7,0))
boxplot(DEDD.results.blood_2_DD$tpr, names=NA, col=col_vector[1:7], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.2), 
        main=paste0("TPR diff dist 2 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in dispersion only"), 
        cex.main=1.7)
lines(c(7.5,7.5), c(0,0.95), col="lightgrey")
legend("topleft", fill=col_vector[1:7], bty='n', cex=1.4, ncol=3, 
       legend=c("diffVar", "expHMM 0.5", "expHMM thr", "expHMM bfdr", 
                "lnHMM 0.5", "lnHMM thr", "lnHMM bfdr"))
boxplot(DEDD.results.blood_5_DD$tpr, names=NA, col=col_vector[1:7], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.2), 
        main=paste0("TPR diff dist 5 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in dispersion only"), 
        cex.main=1.7)
lines(c(7.5,7.5), c(0,0.95), col="lightgrey")
legend("topleft", fill=col_vector[1:7], bty='n', cex=1.4, ncol=3, 
       legend=c("diffVar", "expHMM 0.5", "expHMM thr", "expHMM bfdr", 
                "lnHMM 0.5", "lnHMM thr", "lnHMM bfdr"))
boxplot(DEDD.results.blood_10_DD$tpr, names=NA, col=col_vector[1:7], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.2), 
        main=paste0("TPR diff dist 10 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in dispersion only"), 
        cex.main=1.7)
lines(c(7.5,7.5), c(0,0.95), col="lightgrey")
legend("topleft", fill=col_vector[1:7], bty='n', cex=1.4, ncol=3, 
       legend=c("diffVar", "expHMM 0.5", "expHMM thr", "expHMM bfdr", 
                "lnHMM 0.5", "lnHMM thr", "lnHMM bfdr"))
boxplot(DEDD.results.blood_20_DD$tpr, names=NA, col=col_vector[1:7], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.2), 
        main=paste0("TPR diff dist 20 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in dispersion only"), 
        cex.main=1.7)
lines(c(7.5,7.5), c(0,0.95), col="lightgrey")
legend("topleft", fill=col_vector[1:7], bty='n', cex=1.4, ncol=3, 
       legend=c("diffVar", "expHMM 0.5", "expHMM thr", "expHMM bfdr", 
                "lnHMM 0.5", "lnHMM thr", "lnHMM bfdr"))
boxplot(DEDD.results.blood_50_DD$tpr, names=NA, col=col_vector[1:7], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.2), 
        main=paste0("TPR diff dist 50 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in dispersion only"), 
        cex.main=1.7)
lines(c(7.5,7.5), c(0,0.95), col="lightgrey")
legend("topleft", fill=col_vector[1:7], bty='n', cex=1.4, ncol=3, 
       legend=c("diffVar", "expHMM 0.5", "expHMM thr", "expHMM bfdr", 
                "lnHMM 0.5", "lnHMM thr", "lnHMM bfdr"))
rbind(colMeans(DEDD.results.blood_2_DD$tpr), 
      colMeans(DEDD.results.blood_5_DD$tpr), 
      colMeans(DEDD.results.blood_10_DD$tpr), 
      colMeans(DEDD.results.blood_20_DD$tpr), 
      colMeans(DEDD.results.blood_50_DD$tpr))

par(mfrow=c(2,3), mar=c(1,2.5,5.5,1), mgp=c(3,0.7,0))
boxplot(DEDD.results.blood_2_DE$tpr, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.2), 
        main=paste0("TPR diff dist 2 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean only"), 
        cex.main=1.7)
lines(c(6.5,6.5), c(0,0.95), col="lightgrey")
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.4, ncol=3, 
       legend=c("expHMM 0.5", "expHMM thr", "expHMM bfdr", 
                "lnHMM 0.5", "lnHMM thr", "lnHMM bfdr"))
boxplot(DEDD.results.blood_5_DE$tpr, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.2), 
        main=paste0("TPR diff dist 5 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean only"), 
        cex.main=1.7)
lines(c(6.5,6.5), c(0,0.95), col="lightgrey")
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.4, ncol=3, 
       legend=c("expHMM 0.5", "expHMM thr", "expHMM bfdr", 
                "lnHMM 0.5", "lnHMM thr", "lnHMM bfdr"))
boxplot(DEDD.results.blood_10_DE$tpr, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.2), 
        main=paste0("TPR diff dist 10 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean only"), 
        cex.main=1.7)
lines(c(6.5,6.5), c(0,0.95), col="lightgrey")
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.4, ncol=3, 
       legend=c("expHMM 0.5", "expHMM thr", "expHMM bfdr", 
                "lnHMM 0.5", "lnHMM thr", "lnHMM bfdr"))
boxplot(DEDD.results.blood_20_DE$tpr, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.2), 
        main=paste0("TPR diff dist 20 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean only"), 
        cex.main=1.7)
lines(c(6.5,6.5), c(0,0.95), col="lightgrey")
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.4, ncol=3, 
       legend=c("expHMM 0.5", "expHMM thr", "expHMM bfdr", 
                "lnHMM 0.5", "lnHMM thr", "lnHMM bfdr"))
boxplot(DEDD.results.blood_50_DE$tpr, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.2), 
        main=paste0("TPR diff dist 50 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean only"), 
        cex.main=1.7)
lines(c(6.5,6.5), c(0,0.95), col="lightgrey")
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.4, ncol=3, 
       legend=c("expHMM 0.5", "expHMM thr", "expHMM bfdr", 
                "lnHMM 0.5", "lnHMM thr", "lnHMM bfdr"))
rbind(colMeans(DEDD.results.blood_2_DE$tpr), 
      colMeans(DEDD.results.blood_5_DE$tpr), 
      colMeans(DEDD.results.blood_10_DE$tpr), 
      colMeans(DEDD.results.blood_20_DE$tpr), 
      colMeans(DEDD.results.blood_50_DE$tpr))

par(mfrow=c(2,3), mar=c(1,2.5,5.5,1), mgp=c(3,0.7,0))
boxplot(DEDD.results.blood_2_DEDD$tpr, names=NA, col=col_vector[1:13], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.2), 
        main=paste0("TPR diff dist 2 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean and dispersion"), 
        cex.main=1.7)
lines(c(13.5,13.5), c(0,0.95), col="lightgrey")
legend("topleft", fill=col_vector[1:13], bty='n', cex=1.2, ncol=3, 
       legend=c("diffVar", "expHMM 0.5", "expHMM thr", "expHMM bfdr", 
                "lnHMM 0.5", "lnHMM thr", "lnHMM bfdr", 
                "edgeR_diffVar", "edgeR_MDSeq", "edgeR_lnHM", 
                "lnHM_diffVar", "lnHM_MDSeq", "lnHM_lnHM"))
boxplot(DEDD.results.blood_5_DEDD$tpr, names=NA, col=col_vector[1:13], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.2), 
        main=paste0("TPR diff dist 5 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean and dispersion"), 
        cex.main=1.7)
lines(c(13.5,13.5), c(0,0.95), col="lightgrey")
legend("topleft", fill=col_vector[1:13], bty='n', cex=1.2, ncol=3, 
       legend=c("diffVar", "expHMM 0.5", "expHMM thr", "expHMM bfdr", 
                "lnHMM 0.5", "lnHMM thr", "lnHMM bfdr", 
                "edgeR_diffVar", "edgeR_MDSeq", "edgeR_lnHM", 
                "lnHM_diffVar", "lnHM_MDSeq", "lnHM_lnHM"))
boxplot(DEDD.results.blood_10_DEDD$tpr, names=NA, col=col_vector[1:13], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.2), 
        main=paste0("TPR diff dist 10 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean and dispersion"), 
        cex.main=1.7)
lines(c(13.5,13.5), c(0,0.95), col="lightgrey")
legend("topleft", fill=col_vector[1:13], bty='n', cex=1.2, ncol=3, 
       legend=c("diffVar", "expHMM 0.5", "expHMM thr", "expHMM bfdr", 
                "lnHMM 0.5", "lnHMM thr", "lnHMM bfdr", 
                "edgeR_diffVar", "edgeR_MDSeq", "edgeR_lnHM", 
                "lnHM_diffVar", "lnHM_MDSeq", "lnHM_lnHM"))
boxplot(DEDD.results.blood_20_DEDD$tpr, names=NA, col=col_vector[1:13], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.2), 
        main=paste0("TPR diff dist 20 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean and dispersion"), 
        cex.main=1.7)
lines(c(13.5,13.5), c(0,0.95), col="lightgrey")
legend("topleft", fill=col_vector[1:13], bty='n', cex=1.2, ncol=3, 
       legend=c("diffVar", "expHMM 0.5", "expHMM thr", "expHMM bfdr", 
                "lnHMM 0.5", "lnHMM thr", "lnHMM bfdr", 
                "edgeR_diffVar", "edgeR_MDSeq", "edgeR_lnHM", 
                "lnHM_diffVar", "lnHM_MDSeq", "lnHM_lnHM"))
boxplot(DEDD.results.blood_50_DEDD$tpr, names=NA, col=col_vector[1:13], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.2), 
        main=paste0("TPR diff dist 50 samples per group\n", 
                    "TMM (left) or RLE (right) normalisation\n", 
                    "Differences in mean and dispersion"), 
        cex.main=1.7)
lines(c(13.5,13.5), c(0,0.95), col="lightgrey")
legend("topleft", fill=col_vector[1:13], bty='n', cex=1.2, ncol=3, 
       legend=c("diffVar", "expHMM 0.5", "expHMM thr", "expHMM bfdr", 
                "lnHMM 0.5", "lnHMM thr", "lnHMM bfdr", 
                "edgeR_diffVar", "edgeR_MDSeq", "edgeR_lnHM", 
                "lnHM_diffVar", "lnHM_MDSeq", "lnHM_lnHM"))
rbind(colMeans(DEDD.results.blood_2_DEDD$tpr), 
      colMeans(DEDD.results.blood_5_DEDD$tpr), 
      colMeans(DEDD.results.blood_10_DEDD$tpr), 
      colMeans(DEDD.results.blood_20_DEDD$tpr), 
      colMeans(DEDD.results.blood_50_DEDD$tpr))


###############################
#### False discovery plots ####

par(mfrow=c(6,5), mar=c(2,1.5,2.5,0.5), mgp=c(3,0.5,0))

plot(DEDD.results.blood_2_DD$mean.discoveries$diffVar.tmm_DEDD, 
     DEDD.results.blood_2_DD$mean.fdr$diffVar.tmm_DEDD, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("2 samples per group, TMM\n", 
                 "Differences in dispersion only"))
lines(DEDD.results.blood_2_DD$mean.discoveries$expHMM.tmm, 
      DEDD.results.blood_2_DD$mean.fdr$expHMM.tmm, 
      type='l', col=col_vector[2])
lines(DEDD.results.blood_2_DD$mean.discoveries$lnHMM.tmm, 
      DEDD.results.blood_2_DD$mean.fdr$lnHMM.tmm, 
      type='l', col=col_vector[3])
abline(h=0.05, col='lightgrey')
legend("right", bty='n', col=col_vector[1:3], lty=1, ncol=1, lwd=2, 
       legend=c("diffVar", "expHMM", "lnHMM"))
plot(DEDD.results.blood_5_DD$mean.discoveries$diffVar.tmm_DEDD, 
     DEDD.results.blood_5_DD$mean.fdr$diffVar.tmm_DEDD, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("5 samples per group, TMM\n", 
                 "Differences in dispersion only"))
lines(DEDD.results.blood_5_DD$mean.discoveries$expHMM.tmm, 
      DEDD.results.blood_5_DD$mean.fdr$expHMM.tmm, 
      type='l', col=col_vector[2])
lines(DEDD.results.blood_5_DD$mean.discoveries$lnHMM.tmm, 
      DEDD.results.blood_5_DD$mean.fdr$lnHMM.tmm, 
      type='l', col=col_vector[3])
abline(h=0.05, col='lightgrey')
plot(DEDD.results.blood_10_DD$mean.discoveries$diffVar.tmm_DEDD, 
     DEDD.results.blood_10_DD$mean.fdr$diffVar.tmm_DEDD, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("10 samples per group, TMM\n", 
                 "Differences in dispersion only"))
lines(DEDD.results.blood_10_DD$mean.discoveries$expHMM.tmm, 
      DEDD.results.blood_10_DD$mean.fdr$expHMM.tmm, 
      type='l', col=col_vector[2])
lines(DEDD.results.blood_10_DD$mean.discoveries$lnHMM.tmm, 
      DEDD.results.blood_10_DD$mean.fdr$lnHMM.tmm, 
      type='l', col=col_vector[3])
abline(h=0.05, col='lightgrey')
plot(DEDD.results.blood_20_DD$mean.discoveries$diffVar.tmm_DEDD, 
     DEDD.results.blood_20_DD$mean.fdr$diffVar.tmm_DEDD, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("20 samples per group, TMM\n", 
                 "Differences in dispersion only"))
lines(DEDD.results.blood_20_DD$mean.discoveries$expHMM.tmm, 
      DEDD.results.blood_20_DD$mean.fdr$expHMM.tmm, 
      type='l', col=col_vector[2])
lines(DEDD.results.blood_20_DD$mean.discoveries$lnHMM.tmm, 
      DEDD.results.blood_20_DD$mean.fdr$lnHMM.tmm, 
      type='l', col=col_vector[3])
abline(h=0.05, col='lightgrey')
plot(DEDD.results.blood_50_DD$mean.discoveries$diffVar.tmm_DEDD, 
     DEDD.results.blood_50_DD$mean.fdr$diffVar.tmm_DEDD, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("50 samples per group, TMM\n", 
                 "Differences in dispersion only"))
lines(DEDD.results.blood_50_DD$mean.discoveries$expHMM.tmm, 
      DEDD.results.blood_50_DD$mean.fdr$expHMM.tmm, 
      type='l', col=col_vector[2])
lines(DEDD.results.blood_50_DD$mean.discoveries$lnHMM.tmm, 
      DEDD.results.blood_50_DD$mean.fdr$lnHMM.tmm, 
      type='l', col=col_vector[3])
abline(h=0.05, col='lightgrey')

plot(DEDD.results.blood_2_DD$mean.discoveries$diffVar.rle_DEDD, 
     DEDD.results.blood_2_DD$mean.fdr$diffVar.rle_DEDD, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("2 samples per group, RLE\n", 
                 "Differences in dispersion only"))
lines(DEDD.results.blood_2_DD$mean.discoveries$expHMM.rle, 
      DEDD.results.blood_2_DD$mean.fdr$expHMM.rle, 
      type='l', col=col_vector[2])
lines(DEDD.results.blood_2_DD$mean.discoveries$lnHMM.rle, 
      DEDD.results.blood_2_DD$mean.fdr$lnHMM.rle, 
      type='l', col=col_vector[3])
abline(h=0.05, col='lightgrey')
legend("right", bty='n', col=col_vector[1:3], lty=1, ncol=1, lwd=2, 
       legend=c("diffVar", "expHMM", "lnHMM"))
plot(DEDD.results.blood_5_DD$mean.discoveries$diffVar.rle_DEDD, 
     DEDD.results.blood_5_DD$mean.fdr$diffVar.rle_DEDD, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("5 samples per group, RLE\n", 
                 "Differences in dispersion only"))
lines(DEDD.results.blood_5_DD$mean.discoveries$expHMM.rle, 
      DEDD.results.blood_5_DD$mean.fdr$expHMM.rle, 
      type='l', col=col_vector[2])
lines(DEDD.results.blood_5_DD$mean.discoveries$lnHMM.rle, 
      DEDD.results.blood_5_DD$mean.fdr$lnHMM.rle, 
      type='l', col=col_vector[3])
abline(h=0.05, col='lightgrey')
plot(DEDD.results.blood_10_DD$mean.discoveries$diffVar.rle_DEDD, 
     DEDD.results.blood_10_DD$mean.fdr$diffVar.rle_DEDD, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("10 samples per group, RLE\n", 
                 "Differences in dispersion only"))
lines(DEDD.results.blood_10_DD$mean.discoveries$expHMM.rle, 
      DEDD.results.blood_10_DD$mean.fdr$expHMM.rle, 
      type='l', col=col_vector[2])
lines(DEDD.results.blood_10_DD$mean.discoveries$lnHMM.rle, 
      DEDD.results.blood_10_DD$mean.fdr$lnHMM.rle, 
      type='l', col=col_vector[3])
abline(h=0.05, col='lightgrey')
plot(DEDD.results.blood_20_DD$mean.discoveries$diffVar.rle_DEDD, 
     DEDD.results.blood_20_DD$mean.fdr$diffVar.rle_DEDD, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("20 samples per group, RLE\n", 
                 "Differences in dispersion only"))
lines(DEDD.results.blood_20_DD$mean.discoveries$expHMM.rle, 
      DEDD.results.blood_20_DD$mean.fdr$expHMM.rle, 
      type='l', col=col_vector[2])
lines(DEDD.results.blood_20_DD$mean.discoveries$lnHMM.rle, 
      DEDD.results.blood_20_DD$mean.fdr$lnHMM.rle, 
      type='l', col=col_vector[3])
abline(h=0.05, col='lightgrey')
plot(DEDD.results.blood_50_DD$mean.discoveries$diffVar.rle_DEDD, 
     DEDD.results.blood_50_DD$mean.fdr$diffVar.rle_DEDD, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("50 samples per group, RLE\n", 
                 "Differences in dispersion only"))
lines(DEDD.results.blood_50_DD$mean.discoveries$expHMM.rle, 
      DEDD.results.blood_50_DD$mean.fdr$expHMM.rle, 
      type='l', col=col_vector[2])
lines(DEDD.results.blood_50_DD$mean.discoveries$lnHMM.rle, 
      DEDD.results.blood_50_DD$mean.fdr$lnHMM.rle, 
      type='l', col=col_vector[3])
abline(h=0.05, col='lightgrey')

plot(DEDD.results.blood_2_DE$mean.discoveries$expHMM.tmm, 
     DEDD.results.blood_2_DE$mean.fdr$expHMM.tmm, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("2 samples per group, TMM\n", 
                 "Differences in mean only"))
lines(DEDD.results.blood_2_DE$mean.discoveries$lnHMM.tmm, 
      DEDD.results.blood_2_DE$mean.fdr$lnHMM.tmm, 
      type='l', col=col_vector[2])
abline(h=0.05, col='lightgrey')
plot(DEDD.results.blood_5_DE$mean.discoveries$expHMM.tmm, 
     DEDD.results.blood_5_DE$mean.fdr$expHMM.tmm, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("5 samples per group, TMM\n", 
                 "Differences in mean only"))
lines(DEDD.results.blood_5_DE$mean.discoveries$lnHMM.tmm, 
      DEDD.results.blood_5_DE$mean.fdr$lnHMM.tmm, 
      type='l', col=col_vector[2])
abline(h=0.05, col='lightgrey')
plot(DEDD.results.blood_10_DE$mean.discoveries$expHMM.tmm, 
     DEDD.results.blood_10_DE$mean.fdr$expHMM.tmm, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("10 samples per group, TMM\n", 
                 "Differences in mean only"))
lines(DEDD.results.blood_10_DE$mean.discoveries$lnHMM.tmm, 
      DEDD.results.blood_10_DE$mean.fdr$lnHMM.tmm, 
      type='l', col=col_vector[2])
abline(h=0.05, col='lightgrey')
plot(DEDD.results.blood_20_DE$mean.discoveries$expHMM.tmm, 
     DEDD.results.blood_20_DE$mean.fdr$expHMM.tmm, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("20 samples per group, TMM\n", 
                 "Differences in mean only"))
lines(DEDD.results.blood_20_DE$mean.discoveries$lnHMM.tmm, 
      DEDD.results.blood_20_DE$mean.fdr$lnHMM.tmm, 
      type='l', col=col_vector[2])
abline(h=0.05, col='lightgrey')
plot(DEDD.results.blood_50_DE$mean.discoveries$expHMM.tmm, 
     DEDD.results.blood_50_DE$mean.fdr$expHMM.tmm, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("50 samples per group, TMM\n", 
                 "Differences in mean only"))
lines(DEDD.results.blood_50_DE$mean.discoveries$lnHMM.tmm, 
      DEDD.results.blood_50_DE$mean.fdr$lnHMM.tmm, 
      type='l', col=col_vector[2])
abline(h=0.05, col='lightgrey')
legend("topleft", bty='n', col=col_vector[1:2], lty=1, ncol=1, lwd=2, 
       legend=c("expHMM", "lnHMM"))

plot(DEDD.results.blood_2_DE$mean.discoveries$expHMM.rle, 
     DEDD.results.blood_2_DE$mean.fdr$expHMM.rle, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("2 samples per group, RLE\n", 
                 "Differences in mean only"))
lines(DEDD.results.blood_2_DE$mean.discoveries$lnHMM.rle, 
      DEDD.results.blood_2_DE$mean.fdr$lnHMM.rle, 
      type='l', col=col_vector[2])
abline(h=0.05, col='lightgrey')
plot(DEDD.results.blood_5_DE$mean.discoveries$expHMM.rle, 
     DEDD.results.blood_5_DE$mean.fdr$expHMM.rle, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("5 samples per group, RLE\n", 
                 "Differences in mean only"))
lines(DEDD.results.blood_5_DE$mean.discoveries$lnHMM.rle, 
      DEDD.results.blood_5_DE$mean.fdr$lnHMM.rle, 
      type='l', col=col_vector[2])
abline(h=0.05, col='lightgrey')
plot(DEDD.results.blood_10_DE$mean.discoveries$expHMM.rle, 
     DEDD.results.blood_10_DE$mean.fdr$expHMM.rle, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("10 samples per group, RLE\n", 
                 "Differences in mean only"))
lines(DEDD.results.blood_10_DE$mean.discoveries$lnHMM.rle, 
      DEDD.results.blood_10_DE$mean.fdr$lnHMM.rle, 
      type='l', col=col_vector[2])
abline(h=0.05, col='lightgrey')
plot(DEDD.results.blood_20_DE$mean.discoveries$expHMM.rle, 
     DEDD.results.blood_20_DE$mean.fdr$expHMM.rle, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("20 samples per group, RLE\n", 
                 "Differences in mean only"))
lines(DEDD.results.blood_20_DE$mean.discoveries$lnHMM.rle, 
      DEDD.results.blood_20_DE$mean.fdr$lnHMM.rle, 
      type='l', col=col_vector[2])
abline(h=0.05, col='lightgrey')
plot(DEDD.results.blood_50_DE$mean.discoveries$expHMM.rle, 
     DEDD.results.blood_50_DE$mean.fdr$expHMM.rle, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("50 samples per group, RLE\n", 
                 "Differences in mean only"))
lines(DEDD.results.blood_50_DE$mean.discoveries$lnHMM.rle, 
      DEDD.results.blood_50_DE$mean.fdr$lnHMM.rle, 
      type='l', col=col_vector[2])
abline(h=0.05, col='lightgrey')
legend("topleft", bty='n', col=col_vector[1:2], lty=1, ncol=1, lwd=2, 
       legend=c("expHMM", "lnHMM"))

plot(DEDD.results.blood_2_DEDD$mean.discoveries$diffVar.tmm_DEDD, 
     DEDD.results.blood_2_DEDD$mean.fdr$diffVar.tmm_DEDD, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("2 samples per group, TMM\n", 
                 "Differences in mean and dispersion"))
lines(DEDD.results.blood_2_DEDD$mean.discoveries$expHMM.tmm, 
      DEDD.results.blood_2_DEDD$mean.fdr$expHMM.tmm, 
      type='l', col=col_vector[2])
lines(DEDD.results.blood_2_DEDD$mean.discoveries$lnHMM.tmm, 
      DEDD.results.blood_2_DEDD$mean.fdr$lnHMM.tmm, 
      type='l', col=col_vector[3])
lines(DEDD.results.blood_2_DEDD$mean.discoveries$edgeR_dV.tmm, 
      DEDD.results.blood_2_DEDD$mean.fdr$edgeR_dV.tmm, 
      type='l', col=col_vector[4])
lines(DEDD.results.blood_2_DEDD$mean.discoveries$edgeR_MD.tmm, 
      DEDD.results.blood_2_DEDD$mean.fdr$edgeR_MD.tmm, 
      type='l', col=col_vector[5])
lines(DEDD.results.blood_2_DEDD$mean.discoveries$edgeR_HM.tmm, 
      DEDD.results.blood_2_DEDD$mean.fdr$edgeR_HM.tmm, 
      type='l', col=col_vector[6])
lines(DEDD.results.blood_2_DEDD$mean.discoveries$HM_dV.tmm, 
      DEDD.results.blood_2_DEDD$mean.fdr$HM_dV.tmm, 
      type='l', col=col_vector[7])
lines(DEDD.results.blood_2_DEDD$mean.discoveries$HM_MD.tmm, 
      DEDD.results.blood_2_DEDD$mean.fdr$HM_MD.tmm, 
      type='l', col=col_vector[8])
lines(DEDD.results.blood_2_DEDD$mean.discoveries$HM_HM.tmm, 
      DEDD.results.blood_2_DEDD$mean.fdr$HM_HM.tmm, 
      type='l', col=col_vector[9])
abline(h=0.05, col='lightgrey')
plot(DEDD.results.blood_5_DEDD$mean.discoveries$diffVar.tmm_DEDD, 
     DEDD.results.blood_5_DEDD$mean.fdr$diffVar.tmm_DEDD, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("5 samples per group, TMM\n", 
                 "Differences in mean and dispersion"))
lines(DEDD.results.blood_5_DEDD$mean.discoveries$expHMM.tmm, 
      DEDD.results.blood_5_DEDD$mean.fdr$expHMM.tmm, 
      type='l', col=col_vector[2])
lines(DEDD.results.blood_5_DEDD$mean.discoveries$lnHMM.tmm, 
      DEDD.results.blood_5_DEDD$mean.fdr$lnHMM.tmm, 
      type='l', col=col_vector[3])
lines(DEDD.results.blood_5_DEDD$mean.discoveries$edgeR_dV.tmm, 
      DEDD.results.blood_5_DEDD$mean.fdr$edgeR_dV.tmm, 
      type='l', col=col_vector[4])
lines(DEDD.results.blood_5_DEDD$mean.discoveries$edgeR_MD.tmm, 
      DEDD.results.blood_5_DEDD$mean.fdr$edgeR_MD.tmm, 
      type='l', col=col_vector[5])
lines(DEDD.results.blood_5_DEDD$mean.discoveries$edgeR_HM.tmm, 
      DEDD.results.blood_5_DEDD$mean.fdr$edgeR_HM.tmm, 
      type='l', col=col_vector[6])
lines(DEDD.results.blood_5_DEDD$mean.discoveries$HM_dV.tmm, 
      DEDD.results.blood_5_DEDD$mean.fdr$HM_dV.tmm, 
      type='l', col=col_vector[7])
lines(DEDD.results.blood_5_DEDD$mean.discoveries$HM_MD.tmm, 
      DEDD.results.blood_5_DEDD$mean.fdr$HM_MD.tmm, 
      type='l', col=col_vector[8])
lines(DEDD.results.blood_5_DEDD$mean.discoveries$HM_HM.tmm, 
      DEDD.results.blood_5_DEDD$mean.fdr$HM_HM.tmm, 
      type='l', col=col_vector[9])
abline(h=0.05, col='lightgrey')
plot(DEDD.results.blood_10_DEDD$mean.discoveries$diffVar.tmm_DEDD, 
     DEDD.results.blood_10_DEDD$mean.fdr$diffVar.tmm_DEDD, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("10 samples per group, TMM\n", 
                 "Differences in mean and dispersion"))
lines(DEDD.results.blood_10_DEDD$mean.discoveries$expHMM.tmm, 
      DEDD.results.blood_10_DEDD$mean.fdr$expHMM.tmm, 
      type='l', col=col_vector[2])
lines(DEDD.results.blood_10_DEDD$mean.discoveries$lnHMM.tmm, 
      DEDD.results.blood_10_DEDD$mean.fdr$lnHMM.tmm, 
      type='l', col=col_vector[3])
lines(DEDD.results.blood_10_DEDD$mean.discoveries$edgeR_dV.tmm, 
      DEDD.results.blood_10_DEDD$mean.fdr$edgeR_dV.tmm, 
      type='l', col=col_vector[4])
lines(DEDD.results.blood_10_DEDD$mean.discoveries$edgeR_MD.tmm, 
      DEDD.results.blood_10_DEDD$mean.fdr$edgeR_MD.tmm, 
      type='l', col=col_vector[5])
lines(DEDD.results.blood_10_DEDD$mean.discoveries$edgeR_HM.tmm, 
      DEDD.results.blood_10_DEDD$mean.fdr$edgeR_HM.tmm, 
      type='l', col=col_vector[6])
lines(DEDD.results.blood_10_DEDD$mean.discoveries$HM_dV.tmm, 
      DEDD.results.blood_10_DEDD$mean.fdr$HM_dV.tmm, 
      type='l', col=col_vector[7])
lines(DEDD.results.blood_10_DEDD$mean.discoveries$HM_MD.tmm, 
      DEDD.results.blood_10_DEDD$mean.fdr$HM_MD.tmm, 
      type='l', col=col_vector[7])
lines(DEDD.results.blood_10_DEDD$mean.discoveries$HM_HM.tmm, 
      DEDD.results.blood_10_DEDD$mean.fdr$HM_HM.tmm, 
      type='l', col=col_vector[9])
abline(h=0.05, col='lightgrey')
plot(DEDD.results.blood_20_DEDD$mean.discoveries$diffVar.tmm_DEDD, 
     DEDD.results.blood_20_DEDD$mean.fdr$diffVar.tmm_DEDD, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("20 samples per group, TMM\n", 
                 "Differences in mean and dispersion"))
lines(DEDD.results.blood_20_DEDD$mean.discoveries$expHMM.tmm, 
      DEDD.results.blood_20_DEDD$mean.fdr$expHMM.tmm, 
      type='l', col=col_vector[2])
lines(DEDD.results.blood_20_DEDD$mean.discoveries$lnHMM.tmm, 
      DEDD.results.blood_20_DEDD$mean.fdr$lnHMM.tmm, 
      type='l', col=col_vector[3])
lines(DEDD.results.blood_20_DEDD$mean.discoveries$edgeR_dV.tmm, 
      DEDD.results.blood_20_DEDD$mean.fdr$edgeR_dV.tmm, 
      type='l', col=col_vector[4])
lines(DEDD.results.blood_20_DEDD$mean.discoveries$edgeR_MD.tmm, 
      DEDD.results.blood_20_DEDD$mean.fdr$edgeR_MD.tmm, 
      type='l', col=col_vector[5])
lines(DEDD.results.blood_20_DEDD$mean.discoveries$edgeR_HM.tmm, 
      DEDD.results.blood_20_DEDD$mean.fdr$edgeR_HM.tmm, 
      type='l', col=col_vector[6])
lines(DEDD.results.blood_20_DEDD$mean.discoveries$HM_dV.tmm, 
      DEDD.results.blood_20_DEDD$mean.fdr$HM_dV.tmm, 
      type='l', col=col_vector[7])
lines(DEDD.results.blood_20_DEDD$mean.discoveries$HM_MD.tmm, 
      DEDD.results.blood_20_DEDD$mean.fdr$HM_MD.tmm, 
      type='l', col=col_vector[7])
lines(DEDD.results.blood_20_DEDD$mean.discoveries$HM_HM.tmm, 
      DEDD.results.blood_20_DEDD$mean.fdr$HM_HM.tmm, 
      type='l', col=col_vector[9])
abline(h=0.05, col='lightgrey')
plot(DEDD.results.blood_50_DEDD$mean.discoveries$diffVar.tmm_DEDD, 
     DEDD.results.blood_50_DEDD$mean.fdr$diffVar.tmm_DEDD, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("50 samples per group, TMM\n", 
                 "Differences in mean and dispersion"))
lines(DEDD.results.blood_50_DEDD$mean.discoveries$expHMM.tmm, 
      DEDD.results.blood_50_DEDD$mean.fdr$expHMM.tmm, 
      type='l', col=col_vector[2])
lines(DEDD.results.blood_50_DEDD$mean.discoveries$lnHMM.tmm, 
      DEDD.results.blood_50_DEDD$mean.fdr$lnHMM.tmm, 
      type='l', col=col_vector[3])
lines(DEDD.results.blood_50_DEDD$mean.discoveries$edgeR_dV.tmm, 
      DEDD.results.blood_50_DEDD$mean.fdr$edgeR_dV.tmm, 
      type='l', col=col_vector[4])
lines(DEDD.results.blood_50_DEDD$mean.discoveries$edgeR_MD.tmm, 
      DEDD.results.blood_50_DEDD$mean.fdr$edgeR_MD.tmm, 
      type='l', col=col_vector[5])
lines(DEDD.results.blood_50_DEDD$mean.discoveries$edgeR_HM.tmm, 
      DEDD.results.blood_50_DEDD$mean.fdr$edgeR_HM.tmm, 
      type='l', col=col_vector[6])
lines(DEDD.results.blood_50_DEDD$mean.discoveries$HM_dV.tmm, 
      DEDD.results.blood_50_DEDD$mean.fdr$HM_dV.tmm, 
      type='l', col=col_vector[7])
lines(DEDD.results.blood_50_DEDD$mean.discoveries$HM_MD.tmm, 
      DEDD.results.blood_50_DEDD$mean.fdr$HM_MD.tmm, 
      type='l', col=col_vector[7])
lines(DEDD.results.blood_50_DEDD$mean.discoveries$HM_HM.tmm, 
      DEDD.results.blood_50_DEDD$mean.fdr$HM_HM.tmm, 
      type='l', col=col_vector[9])
abline(h=0.05, col='lightgrey')
legend("topright", bty='n', col=col_vector[1:9], lty=1, ncol=2, lwd=2, 
       legend=c("diffVar", "expHMM", "lnHMM", "edgeR_diffVar", "edgeR_MDSeq", 
                "edgeR_lnHM", "lnHM_diffVar", "lnHM_MDSeq", "lnHM_lnHM"))

plot(DEDD.results.blood_2_DEDD$mean.discoveries$diffVar.rle_DEDD, 
     DEDD.results.blood_2_DEDD$mean.fdr$diffVar.rle_DEDD, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("2 samples per group, RLE\n", 
                 "Differences in mean and dispersion"))
lines(DEDD.results.blood_2_DEDD$mean.discoveries$expHMM.rle, 
      DEDD.results.blood_2_DEDD$mean.fdr$expHMM.rle, 
      type='l', col=col_vector[2])
lines(DEDD.results.blood_2_DEDD$mean.discoveries$lnHMM.rle, 
      DEDD.results.blood_2_DEDD$mean.fdr$lnHMM.rle, 
      type='l', col=col_vector[3])
lines(DEDD.results.blood_2_DEDD$mean.discoveries$edgeR_dV.rle, 
      DEDD.results.blood_2_DEDD$mean.fdr$edgeR_dV.rle, 
      type='l', col=col_vector[4])
lines(DEDD.results.blood_2_DEDD$mean.discoveries$edgeR_MD.rle, 
      DEDD.results.blood_2_DEDD$mean.fdr$edgeR_MD.rle, 
      type='l', col=col_vector[5])
lines(DEDD.results.blood_2_DEDD$mean.discoveries$edgeR_HM.rle, 
      DEDD.results.blood_2_DEDD$mean.fdr$edgeR_HM.rle, 
      type='l', col=col_vector[6])
lines(DEDD.results.blood_2_DEDD$mean.discoveries$HM_dV.rle, 
      DEDD.results.blood_2_DEDD$mean.fdr$HM_dV.rle, 
      type='l', col=col_vector[7])
lines(DEDD.results.blood_2_DEDD$mean.discoveries$HM_MD.rle, 
      DEDD.results.blood_2_DEDD$mean.fdr$HM_MD.rle, 
      type='l', col=col_vector[8])
lines(DEDD.results.blood_2_DEDD$mean.discoveries$HM_HM.rle, 
      DEDD.results.blood_2_DEDD$mean.fdr$HM_HM.rle, 
      type='l', col=col_vector[9])
abline(h=0.05, col='lightgrey')
plot(DEDD.results.blood_5_DEDD$mean.discoveries$diffVar.rle_DEDD, 
     DEDD.results.blood_5_DEDD$mean.fdr$diffVar.rle_DEDD, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("5 samples per group, RLE\n", 
                 "Differences in mean and dispersion"))
lines(DEDD.results.blood_5_DEDD$mean.discoveries$expHMM.rle, 
      DEDD.results.blood_5_DEDD$mean.fdr$expHMM.rle, 
      type='l', col=col_vector[2])
lines(DEDD.results.blood_5_DEDD$mean.discoveries$lnHMM.rle, 
      DEDD.results.blood_5_DEDD$mean.fdr$lnHMM.rle, 
      type='l', col=col_vector[3])
lines(DEDD.results.blood_5_DEDD$mean.discoveries$edgeR_dV.rle, 
      DEDD.results.blood_5_DEDD$mean.fdr$edgeR_dV.rle, 
      type='l', col=col_vector[4])
lines(DEDD.results.blood_5_DEDD$mean.discoveries$edgeR_MD.rle, 
      DEDD.results.blood_5_DEDD$mean.fdr$edgeR_MD.rle, 
      type='l', col=col_vector[5])
lines(DEDD.results.blood_5_DEDD$mean.discoveries$edgeR_HM.rle, 
      DEDD.results.blood_5_DEDD$mean.fdr$edgeR_HM.rle, 
      type='l', col=col_vector[6])
lines(DEDD.results.blood_5_DEDD$mean.discoveries$HM_dV.rle, 
      DEDD.results.blood_5_DEDD$mean.fdr$HM_dV.rle, 
      type='l', col=col_vector[7])
lines(DEDD.results.blood_5_DEDD$mean.discoveries$HM_MD.rle, 
      DEDD.results.blood_5_DEDD$mean.fdr$HM_MD.rle, 
      type='l', col=col_vector[8])
lines(DEDD.results.blood_5_DEDD$mean.discoveries$HM_HM.rle, 
      DEDD.results.blood_5_DEDD$mean.fdr$HM_HM.rle, 
      type='l', col=col_vector[9])
abline(h=0.05, col='lightgrey')
plot(DEDD.results.blood_10_DEDD$mean.discoveries$diffVar.rle_DEDD, 
     DEDD.results.blood_10_DEDD$mean.fdr$diffVar.rle_DEDD, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("10 samples per group, RLE\n", 
                 "Differences in mean and dispersion"))
lines(DEDD.results.blood_10_DEDD$mean.discoveries$expHMM.rle, 
      DEDD.results.blood_10_DEDD$mean.fdr$expHMM.rle, 
      type='l', col=col_vector[2])
lines(DEDD.results.blood_10_DEDD$mean.discoveries$lnHMM.rle, 
      DEDD.results.blood_10_DEDD$mean.fdr$lnHMM.rle, 
      type='l', col=col_vector[3])
lines(DEDD.results.blood_10_DEDD$mean.discoveries$edgeR_dV.rle, 
      DEDD.results.blood_10_DEDD$mean.fdr$edgeR_dV.rle, 
      type='l', col=col_vector[4])
lines(DEDD.results.blood_10_DEDD$mean.discoveries$edgeR_MD.rle, 
      DEDD.results.blood_10_DEDD$mean.fdr$edgeR_MD.rle, 
      type='l', col=col_vector[5])
lines(DEDD.results.blood_10_DEDD$mean.discoveries$edgeR_HM.rle, 
      DEDD.results.blood_10_DEDD$mean.fdr$edgeR_HM.rle, 
      type='l', col=col_vector[6])
lines(DEDD.results.blood_10_DEDD$mean.discoveries$HM_dV.rle, 
      DEDD.results.blood_10_DEDD$mean.fdr$HM_dV.rle, 
      type='l', col=col_vector[7])
lines(DEDD.results.blood_10_DEDD$mean.discoveries$HM_MD.rle, 
      DEDD.results.blood_10_DEDD$mean.fdr$HM_MD.rle, 
      type='l', col=col_vector[7])
lines(DEDD.results.blood_10_DEDD$mean.discoveries$HM_HM.rle, 
      DEDD.results.blood_10_DEDD$mean.fdr$HM_HM.rle, 
      type='l', col=col_vector[9])
abline(h=0.05, col='lightgrey')
plot(DEDD.results.blood_20_DEDD$mean.discoveries$diffVar.rle_DEDD, 
     DEDD.results.blood_20_DEDD$mean.fdr$diffVar.rle_DEDD, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("20 samples per group, RLE\n", 
                 "Differences in mean and dispersion"))
lines(DEDD.results.blood_20_DEDD$mean.discoveries$expHMM.rle, 
      DEDD.results.blood_20_DEDD$mean.fdr$expHMM.rle, 
      type='l', col=col_vector[2])
lines(DEDD.results.blood_20_DEDD$mean.discoveries$lnHMM.rle, 
      DEDD.results.blood_20_DEDD$mean.fdr$lnHMM.rle, 
      type='l', col=col_vector[3])
lines(DEDD.results.blood_20_DEDD$mean.discoveries$edgeR_dV.rle, 
      DEDD.results.blood_20_DEDD$mean.fdr$edgeR_dV.rle, 
      type='l', col=col_vector[4])
lines(DEDD.results.blood_20_DEDD$mean.discoveries$edgeR_MD.rle, 
      DEDD.results.blood_20_DEDD$mean.fdr$edgeR_MD.rle, 
      type='l', col=col_vector[5])
lines(DEDD.results.blood_20_DEDD$mean.discoveries$edgeR_HM.rle, 
      DEDD.results.blood_20_DEDD$mean.fdr$edgeR_HM.rle, 
      type='l', col=col_vector[6])
lines(DEDD.results.blood_20_DEDD$mean.discoveries$HM_dV.rle, 
      DEDD.results.blood_20_DEDD$mean.fdr$HM_dV.rle, 
      type='l', col=col_vector[7])
lines(DEDD.results.blood_20_DEDD$mean.discoveries$HM_MD.rle, 
      DEDD.results.blood_20_DEDD$mean.fdr$HM_MD.rle, 
      type='l', col=col_vector[7])
lines(DEDD.results.blood_20_DEDD$mean.discoveries$HM_HM.rle, 
      DEDD.results.blood_20_DEDD$mean.fdr$HM_HM.rle, 
      type='l', col=col_vector[9])
abline(h=0.05, col='lightgrey')
plot(DEDD.results.blood_50_DEDD$mean.discoveries$diffVar.rle_DEDD, 
     DEDD.results.blood_50_DEDD$mean.fdr$diffVar.rle_DEDD, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("50 samples per group, RLE\n", 
                 "Differences in mean and dispersion"))
lines(DEDD.results.blood_50_DEDD$mean.discoveries$expHMM.rle, 
      DEDD.results.blood_50_DEDD$mean.fdr$expHMM.rle, 
      type='l', col=col_vector[2])
lines(DEDD.results.blood_50_DEDD$mean.discoveries$lnHMM.rle, 
      DEDD.results.blood_50_DEDD$mean.fdr$lnHMM.rle, 
      type='l', col=col_vector[3])
lines(DEDD.results.blood_50_DEDD$mean.discoveries$edgeR_dV.rle, 
      DEDD.results.blood_50_DEDD$mean.fdr$edgeR_dV.rle, 
      type='l', col=col_vector[4])
lines(DEDD.results.blood_50_DEDD$mean.discoveries$edgeR_MD.rle, 
      DEDD.results.blood_50_DEDD$mean.fdr$edgeR_MD.rle, 
      type='l', col=col_vector[5])
lines(DEDD.results.blood_50_DEDD$mean.discoveries$edgeR_HM.rle, 
      DEDD.results.blood_50_DEDD$mean.fdr$edgeR_HM.rle, 
      type='l', col=col_vector[6])
lines(DEDD.results.blood_50_DEDD$mean.discoveries$HM_dV.rle, 
      DEDD.results.blood_50_DEDD$mean.fdr$HM_dV.rle, 
      type='l', col=col_vector[7])
lines(DEDD.results.blood_50_DEDD$mean.discoveries$HM_MD.rle, 
      DEDD.results.blood_50_DEDD$mean.fdr$HM_MD.rle, 
      type='l', col=col_vector[7])
lines(DEDD.results.blood_50_DEDD$mean.discoveries$HM_HM.rle, 
      DEDD.results.blood_50_DEDD$mean.fdr$HM_HM.rle, 
      type='l', col=col_vector[9])
abline(h=0.05, col='lightgrey')
legend("topright", bty='n', col=col_vector[1:9], lty=1, ncol=2, lwd=2, 
       legend=c("diffVar", "expHMM", "lnHMM", "edgeR_diffVar", "edgeR_MDSeq", 
                "edgeR_lnHM", "lnHM_diffVar", "lnHM_MDSeq", "lnHM_lnHM"))

