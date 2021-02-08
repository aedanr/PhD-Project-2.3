library(here)
library(RColorBrewer)

# Results initially saved by dataset, but want to assess by experiment type, 
# i.e. detecting differential expression, dispersion or distribution.

## Load data, create colour vector ####
# folder <- "Results/compcodeR DE, DD, DEDD results Feb 2020"
folder <- "Results/temp_quarantine"
for (i in c("DE", "DD", "DEDD")) {
  for (j in c("DEDD2", "DEDD5", "DEDD10", "DEDD20", "DEDD50")) {
    assign(paste0(i, ".results.", j), 
           readRDS(here(folder, paste0(i, ".results.", j, ".rds"))))
  }
}
for (i in c("DE", "DEDD")) {
  for (j in c("DE2", "DE5", "DE10", "DE20", "DE50")) {
    assign(paste0(i, ".results.", j), 
           readRDS(here(folder, paste0(i, ".results.", j, ".rds"))))
  }
}
for (i in c("DD", "DEDD")) {
  for (j in c("DD2", "DD5", "DD10", "DD20", "DD50")) {
    assign(paste0(i, ".results.", j), 
           readRDS(here(folder, paste0(i, ".results.", j, ".rds"))))
  }
}
rm(i,j)

n <- 23
qual_col_pals = brewer.pal.info[brewer.pal.info$category == "qual" & 
                                  brewer.pal.info$colorblind == T,]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# pie(rep(1,n), col=col_vector)
col_vector <- col_vector[-c(7,10,11,12,20)]


#################################
#### Differential dispersion ####
#################################

###########
### AUC ###
par(mfrow=c(2,3), mar=c(1,2,3,1), mgp=c(3,0.7,0))
boxplot(DD.results.DD2$auc, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0.5,1), 
        main=paste0("AUC for differential dispersion, 2 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in dispersion only"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0.5,0.95))
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=3, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DD.results.DD5$auc, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0.5,1), 
        main=paste0("AUC for differential dispersion, 5 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in dispersion only"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0.5,0.95))
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=3, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DD.results.DD10$auc, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0.5,1), 
        main=paste0("AUC for differential dispersion, 10 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in dispersion only"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0.5,0.95))
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=3, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DD.results.DD20$auc, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0.5,1), 
        main=paste0("AUC for differential dispersion, 20 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in dispersion only"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0.5,0.95))
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=3, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DD.results.DD50$auc, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0.5,1), 
        main=paste0("AUC for differential dispersion, 50 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in dispersion only"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0.5,0.95))
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=3, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", "lnHM", "lnHM log"))
rbind(colMeans(DD.results.DD2$auc), colMeans(DD.results.DD5$auc), 
      colMeans(DD.results.DD10$auc), colMeans(DD.results.DD20$auc), 
      colMeans(DD.results.DD50$auc))

par(mfrow=c(2,3), mar=c(1,2,3,1), mgp=c(3,0.7,0))
boxplot(DD.results.DEDD2$auc, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0.5,1), 
        main=paste0("AUC for differential dispersion, 2 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean and dispersion"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0.5,0.95))
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=3, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DD.results.DEDD5$auc, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0.5,1), 
        main=paste0("AUC for differential dispersion, 5 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean and dispersion"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0.5,0.95))
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=3, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DD.results.DEDD10$auc, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0.5,1), 
        main=paste0("AUC for differential dispersion, 10 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean and dispersion"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0.5,0.95))
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=3, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DD.results.DEDD20$auc, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0.5,1), 
        main=paste0("AUC for differential dispersion, 20 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean and dispersion"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0.5,0.95))
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=3, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DD.results.DEDD50$auc, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0.5,1), 
        main=paste0("AUC for differential dispersion, 50 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean and dispersion"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0.5,0.95))
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=3, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", "lnHM", "lnHM log"))
rbind(colMeans(DD.results.DEDD2$auc), colMeans(DD.results.DEDD5$auc), 
      colMeans(DD.results.DEDD10$auc), colMeans(DD.results.DEDD20$auc), 
      colMeans(DD.results.DEDD50$auc))


###########
### pAUC ###
par(mfrow=c(2,3), mar=c(1,2,3,1), mgp=c(3,0.7,0))
boxplot(DD.results.DD2$pauc, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,0.04), 
        main=paste0("Partial AUC for differential dispersion, 2 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in dispersion only"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0,0.035))
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=3, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DD.results.DD5$pauc, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,0.04), 
        main=paste0("Partial AUC for differential dispersion, 5 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in dispersion only"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0,0.035))
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=3, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DD.results.DD10$pauc, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,0.04), 
        main=paste0("Partial AUC for differential dispersion, 10 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in dispersion only"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0,0.035))
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=3, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DD.results.DD20$pauc, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,0.04), 
        main=paste0("Partial AUC for differential dispersion, 20 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in dispersion only"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0,0.035))
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=3, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DD.results.DD50$pauc, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,0.04), 
        main=paste0("Partial AUC for differential dispersion, 50 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in dispersion only"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0,0.035))
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=3, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", "lnHM", "lnHM log"))
rbind(colMeans(DD.results.DD2$pauc), colMeans(DD.results.DD5$pauc), 
      colMeans(DD.results.DD10$pauc), colMeans(DD.results.DD20$pauc), 
      colMeans(DD.results.DD50$pauc))

par(mfrow=c(2,3), mar=c(1,2,3,1), mgp=c(3,0.7,0))
boxplot(DD.results.DEDD2$pauc, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,0.04), 
        main=paste0("Partial AUC for differential dispersion, 2 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean and dispersion"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0,0.035))
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=3, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DD.results.DEDD5$pauc, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,0.04), 
        main=paste0("Partial AUC for differential dispersion, 5 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean and dispersion"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0,0.035))
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=3, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DD.results.DEDD10$pauc, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,0.04), 
        main=paste0("Partial AUC for differential dispersion, 10 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean and dispersion"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0,0.035))
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=3, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DD.results.DEDD20$pauc, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,0.04), 
        main=paste0("Partial AUC for differential dispersion, 20 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean and dispersion"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0,0.035))
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=3, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DD.results.DEDD50$pauc, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,0.04), 
        main=paste0("Partial AUC for differential dispersion, 50 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean and dispersion"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0,0.035))
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=3, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", "lnHM", "lnHM log"))
rbind(colMeans(DD.results.DEDD2$pauc), colMeans(DD.results.DEDD5$pauc), 
      colMeans(DD.results.DEDD10$pauc), colMeans(DD.results.DEDD20$pauc), 
      colMeans(DD.results.DEDD50$pauc))


###########
### FDR ###
par(mfrow=c(2,3), mar=c(1,2,3,1), mgp=c(3,0.7,0))
boxplot(DD.results.DD2$fdr, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,1.1), 
        main=paste0("FDR for differential dispersion, 2 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in dispersion only"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0,1)); abline(h=0.05, col='lightgrey')
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=3, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DD.results.DD5$fdr, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,1.1), 
        main=paste0("FDR for differential dispersion, 5 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in dispersion only"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0,1)); abline(h=0.05, col='lightgrey')
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=3, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DD.results.DD10$fdr, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,1.1), 
        main=paste0("FDR for differential dispersion, 10 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in dispersion only"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0,1)); abline(h=0.05, col='lightgrey')
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=3, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DD.results.DD20$fdr, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,1.1), 
        main=paste0("FDR for differential dispersion, 20 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in dispersion only"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0,1)); abline(h=0.05, col='lightgrey')
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=3, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DD.results.DD50$fdr, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,1.1), 
        main=paste0("FDR for differential dispersion, 50 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in dispersion only"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0,1)); abline(h=0.05, col='lightgrey')
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=3, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", "lnHM", "lnHM log"))
rbind(colMeans(DD.results.DD2$fdr), colMeans(DD.results.DD5$fdr), 
      colMeans(DD.results.DD10$fdr), colMeans(DD.results.DD20$fdr), 
      colMeans(DD.results.DD50$fdr))

par(mfrow=c(2,3), mar=c(1,2,3,1), mgp=c(3,0.7,0))
boxplot(DD.results.DEDD2$fdr, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,1.1), 
        main=paste0("FDR for differential dispersion, 2 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean and dispersion"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0,1)); abline(h=0.05, col='lightgrey')
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=3, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DD.results.DEDD5$fdr, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,1.1), 
        main=paste0("FDR for differential dispersion, 5 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean and dispersion"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0,1)); abline(h=0.05, col='lightgrey')
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=3, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DD.results.DEDD10$fdr, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,1.1), 
        main=paste0("FDR for differential dispersion, 10 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean and dispersion"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0,1)); abline(h=0.05, col='lightgrey')
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=3, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DD.results.DEDD20$fdr, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,1.1), 
        main=paste0("FDR for differential dispersion, 20 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean and dispersion"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0,1)); abline(h=0.05, col='lightgrey')
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=3, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DD.results.DEDD50$fdr, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,1.1), 
        main=paste0("FDR for differential dispersion, 50 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean and dispersion"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0,1)); abline(h=0.05, col='lightgrey')
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=3, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", "lnHM", "lnHM log"))
rbind(colMeans(DD.results.DEDD2$fdr), colMeans(DD.results.DEDD5$fdr), 
      colMeans(DD.results.DEDD10$fdr), colMeans(DD.results.DEDD20$fdr), 
      colMeans(DD.results.DEDD50$fdr))


###########
### TPR ###
par(mfrow=c(2,3), mar=c(1,2,3,1), mgp=c(3,0.7,0))
boxplot(DD.results.DD2$tpr, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,0.5), 
        main=paste0("TPR for differential dispersion, 2 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in dispersion only"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0,0.45))
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=3, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DD.results.DD5$tpr, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,0.5), 
        main=paste0("TPR for differential dispersion, 5 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in dispersion only"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0,0.45))
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=3, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DD.results.DD10$tpr, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,0.5), 
        main=paste0("TPR for differential dispersion, 10 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in dispersion only"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0,0.45))
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=3, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DD.results.DD20$tpr, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,0.5), 
        main=paste0("TPR for differential dispersion, 20 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in dispersion only"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0,0.45))
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=3, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DD.results.DD50$tpr, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,0.5), 
        main=paste0("TPR for differential dispersion, 50 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in dispersion only"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0,0.45))
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=3, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", "lnHM", "lnHM log"))
rbind(colMeans(DD.results.DD2$tpr), colMeans(DD.results.DD5$tpr), 
      colMeans(DD.results.DD10$tpr), colMeans(DD.results.DD20$tpr), 
      colMeans(DD.results.DD50$tpr))

par(mfrow=c(2,3), mar=c(1,2,3,1), mgp=c(3,0.7,0))
boxplot(DD.results.DEDD2$tpr, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,0.5), 
        main=paste0("TPR for differential dispersion, 2 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean and dispersion"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0,0.45))
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=3, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DD.results.DEDD5$tpr, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,0.5), 
        main=paste0("TPR for differential dispersion, 5 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean and dispersion"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0,0.45))
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=3, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DD.results.DEDD10$tpr, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,0.5), 
        main=paste0("TPR for differential dispersion, 10 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean and dispersion"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0,0.45))
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=3, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DD.results.DEDD20$tpr, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,0.5), 
        main=paste0("TPR for differential dispersion, 20 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean and dispersion"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0,0.45))
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=3, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DD.results.DEDD50$tpr, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,0.5), 
        main=paste0("TPR for differential dispersion, 50 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean and dispersion"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0,0.45))
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=3, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", "lnHM", "lnHM log"))
rbind(colMeans(DD.results.DEDD2$tpr), colMeans(DD.results.DEDD5$tpr), 
      colMeans(DD.results.DEDD10$tpr), colMeans(DD.results.DEDD20$tpr), 
      colMeans(DD.results.DEDD50$tpr))


#############################
### False discovery plots ###
par(mfrow=c(4,5), mar=c(2,2,3,1), mgp=c(3,0.7,0))

plot(DD.results.DD2$mean.discoveries$disp.MDSeq.zi.tmm, 
     DD.results.DD2$mean.fdr$disp.MDSeq.zi.tmm, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("2 samples per group\n", "Differences in dispersion only, ", 
     "TMM"))
lines(DD.results.DD2$mean.discoveries$disp.MDSeq.nozi.tmm, 
      DD.results.DD2$mean.fdr$disp.MDSeq.nozi.tmm, 
      type='l', col=col_vector[2])
lines(DD.results.DD2$mean.discoveries$disp.expHM.untr.tmm, 
      DD.results.DD2$mean.fdr$disp.expHM.untr.tmm, 
      type='l', col=col_vector[3])
lines(DD.results.DD2$mean.discoveries$disp.expHM.log.tmm, 
      DD.results.DD2$mean.fdr$disp.expHM.log.tmm, 
      type='l', col=col_vector[4])
lines(DD.results.DD2$mean.discoveries$disp.lnHM.untr.tmm, 
      DD.results.DD2$mean.fdr$disp.lnHM.untr.tmm, 
      type='l', col=col_vector[5])
lines(DD.results.DD2$mean.discoveries$disp.lnHM.log.tmm, 
      DD.results.DD2$mean.fdr$disp.lnHM.log.tmm, 
      type='l', col=col_vector[6])
legend("bottomright", bty='n', col=col_vector[1:6], lty=1, ncol=2, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", "lnHM", "lnHM log"))
plot(DD.results.DD5$mean.discoveries$disp.MDSeq.zi.tmm, 
     DD.results.DD5$mean.fdr$disp.MDSeq.zi.tmm, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("5 samples per group\n", "Differences in dispersion only, ", 
     "TMM"))
lines(DD.results.DD5$mean.discoveries$disp.MDSeq.nozi.tmm, 
      DD.results.DD5$mean.fdr$disp.MDSeq.nozi.tmm, 
      type='l', col=col_vector[2])
lines(DD.results.DD5$mean.discoveries$disp.expHM.untr.tmm, 
      DD.results.DD5$mean.fdr$disp.expHM.untr.tmm, 
      type='l', col=col_vector[3])
lines(DD.results.DD5$mean.discoveries$disp.expHM.log.tmm, 
      DD.results.DD5$mean.fdr$disp.expHM.log.tmm, 
      type='l', col=col_vector[4])
lines(DD.results.DD5$mean.discoveries$disp.lnHM.untr.tmm, 
      DD.results.DD5$mean.fdr$disp.lnHM.untr.tmm, 
      type='l', col=col_vector[5])
lines(DD.results.DD5$mean.discoveries$disp.lnHM.log.tmm, 
      DD.results.DD5$mean.fdr$disp.lnHM.log.tmm, 
      type='l', col=col_vector[6])
legend("bottomright", bty='n', col=col_vector[1:6], lty=1, ncol=2, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", "lnHM", "lnHM log"))
plot(DD.results.DD10$mean.discoveries$disp.MDSeq.zi.tmm, 
     DD.results.DD10$mean.fdr$disp.MDSeq.zi.tmm, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("10 samples per group\n", "Differences in dispersion only, ", 
     "TMM"))
lines(DD.results.DD10$mean.discoveries$disp.MDSeq.nozi.tmm, 
      DD.results.DD10$mean.fdr$disp.MDSeq.nozi.tmm, 
      type='l', col=col_vector[2])
lines(DD.results.DD10$mean.discoveries$disp.expHM.untr.tmm, 
      DD.results.DD10$mean.fdr$disp.expHM.untr.tmm, 
      type='l', col=col_vector[3])
lines(DD.results.DD10$mean.discoveries$disp.expHM.log.tmm, 
      DD.results.DD10$mean.fdr$disp.expHM.log.tmm, 
      type='l', col=col_vector[4])
lines(DD.results.DD10$mean.discoveries$disp.lnHM.untr.tmm, 
      DD.results.DD10$mean.fdr$disp.lnHM.untr.tmm, 
      type='l', col=col_vector[5])
lines(DD.results.DD10$mean.discoveries$disp.lnHM.log.tmm, 
      DD.results.DD10$mean.fdr$disp.lnHM.log.tmm, 
      type='l', col=col_vector[6])
legend("bottomright", bty='n', col=col_vector[1:6], lty=1, ncol=2, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", "lnHM", "lnHM log"))
plot(DD.results.DD20$mean.discoveries$disp.MDSeq.zi.tmm, 
     DD.results.DD20$mean.fdr$disp.MDSeq.zi.tmm, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("20 samples per group\n", "Differences in dispersion only, ", 
     "TMM"))
lines(DD.results.DD20$mean.discoveries$disp.MDSeq.nozi.tmm, 
      DD.results.DD20$mean.fdr$disp.MDSeq.nozi.tmm, 
      type='l', col=col_vector[2])
lines(DD.results.DD20$mean.discoveries$disp.expHM.untr.tmm, 
      DD.results.DD20$mean.fdr$disp.expHM.untr.tmm, 
      type='l', col=col_vector[3])
lines(DD.results.DD20$mean.discoveries$disp.expHM.log.tmm, 
      DD.results.DD20$mean.fdr$disp.expHM.log.tmm, 
      type='l', col=col_vector[4])
lines(DD.results.DD20$mean.discoveries$disp.lnHM.untr.tmm, 
      DD.results.DD20$mean.fdr$disp.lnHM.untr.tmm, 
      type='l', col=col_vector[5])
lines(DD.results.DD20$mean.discoveries$disp.lnHM.log.tmm, 
      DD.results.DD20$mean.fdr$disp.lnHM.log.tmm, 
      type='l', col=col_vector[6])
legend("bottomright", bty='n', col=col_vector[1:6], lty=1, ncol=2, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", "lnHM", "lnHM log"))
plot(DD.results.DD50$mean.discoveries$disp.MDSeq.zi.tmm, 
     DD.results.DD50$mean.fdr$disp.MDSeq.zi.tmm, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("2 samples per group\n", "Differences in dispersion only, ", 
     "TMM"))
lines(DD.results.DD50$mean.discoveries$disp.MDSeq.nozi.tmm, 
      DD.results.DD50$mean.fdr$disp.MDSeq.nozi.tmm, 
      type='l', col=col_vector[2])
lines(DD.results.DD50$mean.discoveries$disp.expHM.untr.tmm, 
      DD.results.DD50$mean.fdr$disp.expHM.untr.tmm, 
      type='l', col=col_vector[3])
lines(DD.results.DD50$mean.discoveries$disp.expHM.log.tmm, 
      DD.results.DD50$mean.fdr$disp.expHM.log.tmm, 
      type='l', col=col_vector[4])
lines(DD.results.DD50$mean.discoveries$disp.lnHM.untr.tmm, 
      DD.results.DD50$mean.fdr$disp.lnHM.untr.tmm, 
      type='l', col=col_vector[5])
lines(DD.results.DD50$mean.discoveries$disp.lnHM.log.tmm, 
      DD.results.DD50$mean.fdr$disp.lnHM.log.tmm, 
      type='l', col=col_vector[6])
legend("topleft", bty='n', col=col_vector[1:6], lty=1, ncol=2, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", "lnHM", "lnHM log"))
abline(h=0.05, col='lightgrey')

plot(DD.results.DD2$mean.discoveries$disp.MDSeq.zi.rle, 
     DD.results.DD2$mean.fdr$disp.MDSeq.zi.rle, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("2 samples per group\n", "Differences in dispersion only, ", 
     "RLE"))
lines(DD.results.DD2$mean.discoveries$disp.MDSeq.nozi.rle, 
      DD.results.DD2$mean.fdr$disp.MDSeq.nozi.rle, 
      type='l', col=col_vector[2])
lines(DD.results.DD2$mean.discoveries$disp.expHM.untr.rle, 
      DD.results.DD2$mean.fdr$disp.expHM.untr.rle, 
      type='l', col=col_vector[3])
lines(DD.results.DD2$mean.discoveries$disp.expHM.log.rle, 
      DD.results.DD2$mean.fdr$disp.expHM.log.rle, 
      type='l', col=col_vector[4])
lines(DD.results.DD2$mean.discoveries$disp.lnHM.untr.rle, 
      DD.results.DD2$mean.fdr$disp.lnHM.untr.rle, 
      type='l', col=col_vector[5])
lines(DD.results.DD2$mean.discoveries$disp.lnHM.log.rle, 
      DD.results.DD2$mean.fdr$disp.lnHM.log.rle, 
      type='l', col=col_vector[6])
legend("bottomright", bty='n', col=col_vector[1:6], lty=1, ncol=2, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", "lnHM", "lnHM log"))
plot(DD.results.DD5$mean.discoveries$disp.MDSeq.zi.rle, 
     DD.results.DD5$mean.fdr$disp.MDSeq.zi.rle, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("5 samples per group\n", "Differences in dispersion only, ", 
     "RLE"))
lines(DD.results.DD5$mean.discoveries$disp.MDSeq.nozi.rle, 
      DD.results.DD5$mean.fdr$disp.MDSeq.nozi.rle, 
      type='l', col=col_vector[2])
lines(DD.results.DD5$mean.discoveries$disp.expHM.untr.rle, 
      DD.results.DD5$mean.fdr$disp.expHM.untr.rle, 
      type='l', col=col_vector[3])
lines(DD.results.DD5$mean.discoveries$disp.expHM.log.rle, 
      DD.results.DD5$mean.fdr$disp.expHM.log.rle, 
      type='l', col=col_vector[4])
lines(DD.results.DD5$mean.discoveries$disp.lnHM.untr.rle, 
      DD.results.DD5$mean.fdr$disp.lnHM.untr.rle, 
      type='l', col=col_vector[5])
lines(DD.results.DD5$mean.discoveries$disp.lnHM.log.rle, 
      DD.results.DD5$mean.fdr$disp.lnHM.log.rle, 
      type='l', col=col_vector[6])
legend("bottomright", bty='n', col=col_vector[1:6], lty=1, ncol=2, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", "lnHM", "lnHM log"))
plot(DD.results.DD10$mean.discoveries$disp.MDSeq.zi.rle, 
     DD.results.DD10$mean.fdr$disp.MDSeq.zi.rle, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("10 samples per group\n", "Differences in dispersion only, ", 
     "RLE"))
lines(DD.results.DD10$mean.discoveries$disp.MDSeq.nozi.rle, 
      DD.results.DD10$mean.fdr$disp.MDSeq.nozi.rle, 
      type='l', col=col_vector[2])
lines(DD.results.DD10$mean.discoveries$disp.expHM.untr.rle, 
      DD.results.DD10$mean.fdr$disp.expHM.untr.rle, 
      type='l', col=col_vector[3])
lines(DD.results.DD10$mean.discoveries$disp.expHM.log.rle, 
      DD.results.DD10$mean.fdr$disp.expHM.log.rle, 
      type='l', col=col_vector[4])
lines(DD.results.DD10$mean.discoveries$disp.lnHM.untr.rle, 
      DD.results.DD10$mean.fdr$disp.lnHM.untr.rle, 
      type='l', col=col_vector[5])
lines(DD.results.DD10$mean.discoveries$disp.lnHM.log.rle, 
      DD.results.DD10$mean.fdr$disp.lnHM.log.rle, 
      type='l', col=col_vector[6])
legend("bottomright", bty='n', col=col_vector[1:6], lty=1, ncol=2, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", "lnHM", "lnHM log"))
plot(DD.results.DD20$mean.discoveries$disp.MDSeq.zi.rle, 
     DD.results.DD20$mean.fdr$disp.MDSeq.zi.rle, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("20 samples per group\n", "Differences in dispersion only, ", 
     "RLE"))
lines(DD.results.DD20$mean.discoveries$disp.MDSeq.nozi.rle, 
      DD.results.DD20$mean.fdr$disp.MDSeq.nozi.rle, 
      type='l', col=col_vector[2])
lines(DD.results.DD20$mean.discoveries$disp.expHM.untr.rle, 
      DD.results.DD20$mean.fdr$disp.expHM.untr.rle, 
      type='l', col=col_vector[3])
lines(DD.results.DD20$mean.discoveries$disp.expHM.log.rle, 
      DD.results.DD20$mean.fdr$disp.expHM.log.rle, 
      type='l', col=col_vector[4])
lines(DD.results.DD20$mean.discoveries$disp.lnHM.untr.rle, 
      DD.results.DD20$mean.fdr$disp.lnHM.untr.rle, 
      type='l', col=col_vector[5])
lines(DD.results.DD20$mean.discoveries$disp.lnHM.log.rle, 
      DD.results.DD20$mean.fdr$disp.lnHM.log.rle, 
      type='l', col=col_vector[6])
legend("bottomright", bty='n', col=col_vector[1:6], lty=1, ncol=2, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", "lnHM", "lnHM log"))
plot(DD.results.DD50$mean.discoveries$disp.MDSeq.zi.rle, 
     DD.results.DD50$mean.fdr$disp.MDSeq.zi.rle, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("2 samples per group\n", "Differences in dispersion only, ", 
     "RLE"))
lines(DD.results.DD50$mean.discoveries$disp.MDSeq.nozi.rle, 
      DD.results.DD50$mean.fdr$disp.MDSeq.nozi.rle, 
      type='l', col=col_vector[2])
lines(DD.results.DD50$mean.discoveries$disp.expHM.untr.rle, 
      DD.results.DD50$mean.fdr$disp.expHM.untr.rle, 
      type='l', col=col_vector[3])
lines(DD.results.DD50$mean.discoveries$disp.expHM.log.rle, 
      DD.results.DD50$mean.fdr$disp.expHM.log.rle, 
      type='l', col=col_vector[4])
lines(DD.results.DD50$mean.discoveries$disp.lnHM.untr.rle, 
      DD.results.DD50$mean.fdr$disp.lnHM.untr.rle, 
      type='l', col=col_vector[5])
lines(DD.results.DD50$mean.discoveries$disp.lnHM.log.rle, 
      DD.results.DD50$mean.fdr$disp.lnHM.log.rle, 
      type='l', col=col_vector[6])
legend("topleft", bty='n', col=col_vector[1:6], lty=1, ncol=2, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", "lnHM", "lnHM log"))
abline(h=0.05, col='lightgrey')

plot(DD.results.DEDD2$mean.discoveries$disp.MDSeq.zi.tmm, 
     DD.results.DEDD2$mean.fdr$disp.MDSeq.zi.tmm, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("2 samples per group\n", "Differences in mean & dispersion, ", 
     "TMM"))
lines(DD.results.DEDD2$mean.discoveries$disp.MDSeq.nozi.tmm, 
      DD.results.DEDD2$mean.fdr$disp.MDSeq.nozi.tmm, 
      type='l', col=col_vector[2])
lines(DD.results.DEDD2$mean.discoveries$disp.expHM.untr.tmm, 
      DD.results.DEDD2$mean.fdr$disp.expHM.untr.tmm, 
      type='l', col=col_vector[3])
lines(DD.results.DEDD2$mean.discoveries$disp.expHM.log.tmm, 
      DD.results.DEDD2$mean.fdr$disp.expHM.log.tmm, 
      type='l', col=col_vector[4])
lines(DD.results.DEDD2$mean.discoveries$disp.lnHM.untr.tmm, 
      DD.results.DEDD2$mean.fdr$disp.lnHM.untr.tmm, 
      type='l', col=col_vector[5])
lines(DD.results.DEDD2$mean.discoveries$disp.lnHM.log.tmm, 
      DD.results.DEDD2$mean.fdr$disp.lnHM.log.tmm, 
      type='l', col=col_vector[6])
legend("bottomright", bty='n', col=col_vector[1:6], lty=1, ncol=2, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", "lnHM", "lnHM log"))
plot(DD.results.DEDD5$mean.discoveries$disp.MDSeq.zi.tmm, 
     DD.results.DEDD5$mean.fdr$disp.MDSeq.zi.tmm, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("5 samples per group\n", "Differences in mean & dispersion, ", 
     "TMM"))
lines(DD.results.DEDD5$mean.discoveries$disp.MDSeq.nozi.tmm, 
      DD.results.DEDD5$mean.fdr$disp.MDSeq.nozi.tmm, 
      type='l', col=col_vector[2])
lines(DD.results.DEDD5$mean.discoveries$disp.expHM.untr.tmm, 
      DD.results.DEDD5$mean.fdr$disp.expHM.untr.tmm, 
      type='l', col=col_vector[3])
lines(DD.results.DEDD5$mean.discoveries$disp.expHM.log.tmm, 
      DD.results.DEDD5$mean.fdr$disp.expHM.log.tmm, 
      type='l', col=col_vector[4])
lines(DD.results.DEDD5$mean.discoveries$disp.lnHM.untr.tmm, 
      DD.results.DEDD5$mean.fdr$disp.lnHM.untr.tmm, 
      type='l', col=col_vector[5])
lines(DD.results.DEDD5$mean.discoveries$disp.lnHM.log.tmm, 
      DD.results.DEDD5$mean.fdr$disp.lnHM.log.tmm, 
      type='l', col=col_vector[6])
legend("bottomright", bty='n', col=col_vector[1:6], lty=1, ncol=2, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", "lnHM", "lnHM log"))
plot(DD.results.DEDD10$mean.discoveries$disp.MDSeq.zi.tmm, 
     DD.results.DEDD10$mean.fdr$disp.MDSeq.zi.tmm, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("10 samples per group\n", "Differences in mean & dispersion, ", 
     "TMM"))
lines(DD.results.DEDD10$mean.discoveries$disp.MDSeq.nozi.tmm, 
      DD.results.DEDD10$mean.fdr$disp.MDSeq.nozi.tmm, 
      type='l', col=col_vector[2])
lines(DD.results.DEDD10$mean.discoveries$disp.expHM.untr.tmm, 
      DD.results.DEDD10$mean.fdr$disp.expHM.untr.tmm, 
      type='l', col=col_vector[3])
lines(DD.results.DEDD10$mean.discoveries$disp.expHM.log.tmm, 
      DD.results.DEDD10$mean.fdr$disp.expHM.log.tmm, 
      type='l', col=col_vector[4])
lines(DD.results.DEDD10$mean.discoveries$disp.lnHM.untr.tmm, 
      DD.results.DEDD10$mean.fdr$disp.lnHM.untr.tmm, 
      type='l', col=col_vector[5])
lines(DD.results.DEDD10$mean.discoveries$disp.lnHM.log.tmm, 
      DD.results.DEDD10$mean.fdr$disp.lnHM.log.tmm, 
      type='l', col=col_vector[6])
legend("bottomright", bty='n', col=col_vector[1:6], lty=1, ncol=2, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", "lnHM", "lnHM log"))
plot(DD.results.DEDD20$mean.discoveries$disp.MDSeq.zi.tmm, 
     DD.results.DEDD20$mean.fdr$disp.MDSeq.zi.tmm, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("20 samples per group\n", "Differences in mean & dispersion, ", 
     "TMM"))
lines(DD.results.DEDD20$mean.discoveries$disp.MDSeq.nozi.tmm, 
      DD.results.DEDD20$mean.fdr$disp.MDSeq.nozi.tmm, 
      type='l', col=col_vector[2])
lines(DD.results.DEDD20$mean.discoveries$disp.expHM.untr.tmm, 
      DD.results.DEDD20$mean.fdr$disp.expHM.untr.tmm, 
      type='l', col=col_vector[3])
lines(DD.results.DEDD20$mean.discoveries$disp.expHM.log.tmm, 
      DD.results.DEDD20$mean.fdr$disp.expHM.log.tmm, 
      type='l', col=col_vector[4])
lines(DD.results.DEDD20$mean.discoveries$disp.lnHM.untr.tmm, 
      DD.results.DEDD20$mean.fdr$disp.lnHM.untr.tmm, 
      type='l', col=col_vector[5])
lines(DD.results.DEDD20$mean.discoveries$disp.lnHM.log.tmm, 
      DD.results.DEDD20$mean.fdr$disp.lnHM.log.tmm, 
      type='l', col=col_vector[6])
legend("topleft", bty='n', col=col_vector[1:6], lty=1, ncol=2, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", "lnHM", "lnHM log"))
abline(h=0.05, col='lightgrey')
plot(DD.results.DEDD50$mean.discoveries$disp.MDSeq.zi.tmm, 
     DD.results.DEDD50$mean.fdr$disp.MDSeq.zi.tmm, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("50 samples per group\n", "Differences in mean & dispersion, ", 
     "TMM"))
lines(DD.results.DEDD50$mean.discoveries$disp.MDSeq.nozi.tmm, 
      DD.results.DEDD50$mean.fdr$disp.MDSeq.nozi.tmm, 
      type='l', col=col_vector[2])
lines(DD.results.DEDD50$mean.discoveries$disp.expHM.untr.tmm, 
      DD.results.DEDD50$mean.fdr$disp.expHM.untr.tmm, 
      type='l', col=col_vector[3])
lines(DD.results.DEDD50$mean.discoveries$disp.expHM.log.tmm, 
      DD.results.DEDD50$mean.fdr$disp.expHM.log.tmm, 
      type='l', col=col_vector[4])
lines(DD.results.DEDD50$mean.discoveries$disp.lnHM.untr.tmm, 
      DD.results.DEDD50$mean.fdr$disp.lnHM.untr.tmm, 
      type='l', col=col_vector[5])
lines(DD.results.DEDD50$mean.discoveries$disp.lnHM.log.tmm, 
      DD.results.DEDD50$mean.fdr$disp.lnHM.log.tmm, 
      type='l', col=col_vector[6])
legend("topleft", bty='n', col=col_vector[1:6], lty=1, ncol=2, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", "lnHM", "lnHM log"))
abline(h=0.05, col='lightgrey')

plot(DD.results.DEDD2$mean.discoveries$disp.MDSeq.zi.rle, 
     DD.results.DEDD2$mean.fdr$disp.MDSeq.zi.rle, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("2 samples per group\n", "Differences in mean & dispersion, ", 
     "RLE"))
lines(DD.results.DEDD2$mean.discoveries$disp.MDSeq.nozi.rle, 
      DD.results.DEDD2$mean.fdr$disp.MDSeq.nozi.rle, 
      type='l', col=col_vector[2])
lines(DD.results.DEDD2$mean.discoveries$disp.expHM.untr.rle, 
      DD.results.DEDD2$mean.fdr$disp.expHM.untr.rle, 
      type='l', col=col_vector[3])
lines(DD.results.DEDD2$mean.discoveries$disp.expHM.log.rle, 
      DD.results.DEDD2$mean.fdr$disp.expHM.log.rle, 
      type='l', col=col_vector[4])
lines(DD.results.DEDD2$mean.discoveries$disp.lnHM.untr.rle, 
      DD.results.DEDD2$mean.fdr$disp.lnHM.untr.rle, 
      type='l', col=col_vector[5])
lines(DD.results.DEDD2$mean.discoveries$disp.lnHM.log.rle, 
      DD.results.DEDD2$mean.fdr$disp.lnHM.log.rle, 
      type='l', col=col_vector[6])
legend("bottomright", bty='n', col=col_vector[1:6], lty=1, ncol=2, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", "lnHM", "lnHM log"))
plot(DD.results.DEDD5$mean.discoveries$disp.MDSeq.zi.rle, 
     DD.results.DEDD5$mean.fdr$disp.MDSeq.zi.rle, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("5 samples per group\n", "Differences in mean & dispersion, ", 
     "RLE"))
lines(DD.results.DEDD5$mean.discoveries$disp.MDSeq.nozi.rle, 
      DD.results.DEDD5$mean.fdr$disp.MDSeq.nozi.rle, 
      type='l', col=col_vector[2])
lines(DD.results.DEDD5$mean.discoveries$disp.expHM.untr.rle, 
      DD.results.DEDD5$mean.fdr$disp.expHM.untr.rle, 
      type='l', col=col_vector[3])
lines(DD.results.DEDD5$mean.discoveries$disp.expHM.log.rle, 
      DD.results.DEDD5$mean.fdr$disp.expHM.log.rle, 
      type='l', col=col_vector[4])
lines(DD.results.DEDD5$mean.discoveries$disp.lnHM.untr.rle, 
      DD.results.DEDD5$mean.fdr$disp.lnHM.untr.rle, 
      type='l', col=col_vector[5])
lines(DD.results.DEDD5$mean.discoveries$disp.lnHM.log.rle, 
      DD.results.DEDD5$mean.fdr$disp.lnHM.log.rle, 
      type='l', col=col_vector[6])
legend("bottomright", bty='n', col=col_vector[1:6], lty=1, ncol=2, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", "lnHM", "lnHM log"))
plot(DD.results.DEDD10$mean.discoveries$disp.MDSeq.zi.rle, 
     DD.results.DEDD10$mean.fdr$disp.MDSeq.zi.rle, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("10 samples per group\n", "Differences in mean & dispersion, ", 
     "RLE"))
lines(DD.results.DEDD10$mean.discoveries$disp.MDSeq.nozi.rle, 
      DD.results.DEDD10$mean.fdr$disp.MDSeq.nozi.rle, 
      type='l', col=col_vector[2])
lines(DD.results.DEDD10$mean.discoveries$disp.expHM.untr.rle, 
      DD.results.DEDD10$mean.fdr$disp.expHM.untr.rle, 
      type='l', col=col_vector[3])
lines(DD.results.DEDD10$mean.discoveries$disp.expHM.log.rle, 
      DD.results.DEDD10$mean.fdr$disp.expHM.log.rle, 
      type='l', col=col_vector[4])
lines(DD.results.DEDD10$mean.discoveries$disp.lnHM.untr.rle, 
      DD.results.DEDD10$mean.fdr$disp.lnHM.untr.rle, 
      type='l', col=col_vector[5])
lines(DD.results.DEDD10$mean.discoveries$disp.lnHM.log.rle, 
      DD.results.DEDD10$mean.fdr$disp.lnHM.log.rle, 
      type='l', col=col_vector[6])
legend("bottomright", bty='n', col=col_vector[1:6], lty=1, ncol=2, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", "lnHM", "lnHM log"))
plot(DD.results.DEDD20$mean.discoveries$disp.MDSeq.zi.rle, 
     DD.results.DEDD20$mean.fdr$disp.MDSeq.zi.rle, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("20 samples per group\n", "Differences in mean & dispersion, ", 
     "RLE"))
lines(DD.results.DEDD20$mean.discoveries$disp.MDSeq.nozi.rle, 
      DD.results.DEDD20$mean.fdr$disp.MDSeq.nozi.rle, 
      type='l', col=col_vector[2])
lines(DD.results.DEDD20$mean.discoveries$disp.expHM.untr.rle, 
      DD.results.DEDD20$mean.fdr$disp.expHM.untr.rle, 
      type='l', col=col_vector[3])
lines(DD.results.DEDD20$mean.discoveries$disp.expHM.log.rle, 
      DD.results.DEDD20$mean.fdr$disp.expHM.log.rle, 
      type='l', col=col_vector[4])
lines(DD.results.DEDD20$mean.discoveries$disp.lnHM.untr.rle, 
      DD.results.DEDD20$mean.fdr$disp.lnHM.untr.rle, 
      type='l', col=col_vector[5])
lines(DD.results.DEDD20$mean.discoveries$disp.lnHM.log.rle, 
      DD.results.DEDD20$mean.fdr$disp.lnHM.log.rle, 
      type='l', col=col_vector[6])
legend("topleft", bty='n', col=col_vector[1:6], lty=1, ncol=2, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", "lnHM", "lnHM log"))
abline(h=0.05, col='lightgrey')
plot(DD.results.DEDD50$mean.discoveries$disp.MDSeq.zi.rle, 
     DD.results.DEDD50$mean.fdr$disp.MDSeq.zi.rle, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("50 samples per group\n", "Differences in mean & dispersion, ", 
     "RLE"))
lines(DD.results.DEDD50$mean.discoveries$disp.MDSeq.nozi.rle, 
      DD.results.DEDD50$mean.fdr$disp.MDSeq.nozi.rle, 
      type='l', col=col_vector[2])
lines(DD.results.DEDD50$mean.discoveries$disp.expHM.untr.rle, 
      DD.results.DEDD50$mean.fdr$disp.expHM.untr.rle, 
      type='l', col=col_vector[3])
lines(DD.results.DEDD50$mean.discoveries$disp.expHM.log.rle, 
      DD.results.DEDD50$mean.fdr$disp.expHM.log.rle, 
      type='l', col=col_vector[4])
lines(DD.results.DEDD50$mean.discoveries$disp.lnHM.untr.rle, 
      DD.results.DEDD50$mean.fdr$disp.lnHM.untr.rle, 
      type='l', col=col_vector[5])
lines(DD.results.DEDD50$mean.discoveries$disp.lnHM.log.rle, 
      DD.results.DEDD50$mean.fdr$disp.lnHM.log.rle, 
      type='l', col=col_vector[6])
legend("topleft", bty='n', col=col_vector[1:6], lty=1, ncol=2, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", "lnHM", "lnHM log"))
abline(h=0.05, col='lightgrey')



#################################
#### Differential expression ####
#################################

###########
### AUC ###
par(mfrow=c(2,3), mar=c(1,2,3,1), mgp=c(3,0.7,0))
boxplot(DE.results.DE2$auc, names=NA, col=col_vector[1:14], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0.75,1.05), 
        main=paste0("AUC for differential expression, 2 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean only"), 
        cex.main=1.2)
lines(c(14.5,14.5), c(0,1))
legend("topleft", fill=col_vector[1:14], bty='n', cex=1.2, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                  "voom", "DSS", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                  "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DE.results.DE5$auc, names=NA, col=col_vector[1:14], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0.75,1.05), 
        main=paste0("AUC for differential expression, 5 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean only"), 
        cex.main=1.2)
lines(c(14.5,14.5), c(0,1))
legend("topleft", fill=col_vector[1:14], bty='n', cex=1.2, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                  "voom", "DSS", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                  "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DE.results.DE10$auc, names=NA, col=col_vector[1:14], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0.75,1.05), 
        main=paste0("AUC for differential expression, 10 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean only"), 
        cex.main=1.2)
lines(c(14.5,14.5), c(0,1))
legend("topleft", fill=col_vector[1:14], bty='n', cex=1.2, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                  "voom", "DSS", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                  "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DE.results.DE20$auc, names=NA, col=col_vector[1:14], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0.75,1.05), 
        main=paste0("AUC for differential expression, 20 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean only"), 
        cex.main=1.2)
lines(c(14.5,14.5), c(0,1))
legend("topleft", fill=col_vector[1:14], bty='n', cex=1.2, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                  "voom", "DSS", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                  "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DE.results.DE50$auc, names=NA, col=col_vector[c(1:6, 8:14)], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0.75,1.05), 
        main=paste0("AUC for differential expression, 50 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean only"), 
        cex.main=1.2)
lines(c(13.5,13.5), c(0,1))
legend("topleft", fill=col_vector[c(1:6, 8:14)], bty='n', cex=1.2, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                  "voom", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                  "expHM", "expHM log", "lnHM", "lnHM log"))
rbind(colMeans(DE.results.DE2$auc), colMeans(DE.results.DE5$auc), 
      colMeans(DE.results.DE10$auc), colMeans(DE.results.DE20$auc), 
      colMeans(DE.results.DE50$auc))

par(mfrow=c(2,3), mar=c(1,2,3,1), mgp=c(3,0.7,0))
boxplot(DE.results.DEDD2$auc, names=NA, col=col_vector[1:14], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0.75,1.05), 
        main=paste0("AUC for differential expression, 2 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean and dispersion"), 
        cex.main=1.2)
lines(c(14.5,14.5), c(0,1))
legend("topleft", fill=col_vector[1:14], bty='n', cex=1.2, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                  "voom", "DSS", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                  "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DE.results.DEDD5$auc, names=NA, col=col_vector[1:14], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0.75,1.05), 
        main=paste0("AUC for differential expression, 5 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean and dispersion"), 
        cex.main=1.2)
lines(c(14.5,14.5), c(0,1))
legend("topleft", fill=col_vector[1:14], bty='n', cex=1.2, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                  "voom", "DSS", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                  "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DE.results.DEDD10$auc, names=NA, col=col_vector[1:14], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0.75,1.05), 
        main=paste0("AUC for differential expression, 10 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean and dispersion"), 
        cex.main=1.2)
lines(c(14.5,14.5), c(0,1))
legend("topleft", fill=col_vector[1:14], bty='n', cex=1.2, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                  "voom", "DSS", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                  "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DE.results.DEDD20$auc, names=NA, col=col_vector[1:14], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0.75,1.05), 
        main=paste0("AUC for differential expression, 20 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean and dispersion"), 
        cex.main=1.2)
lines(c(14.5,14.5), c(0,1))
legend("topleft", fill=col_vector[1:14], bty='n', cex=1.2, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                  "voom", "DSS", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                  "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DE.results.DEDD50$auc, names=NA, col=col_vector[c(1:6, 8:14)], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0.75,1.05), 
        main=paste0("AUC for differential expression, 50 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean and dispersion"), 
        cex.main=1.2)
lines(c(13.5,13.5), c(0,1))
legend("topleft", fill=col_vector[c(1:6, 8:14)], bty='n', cex=1.2, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                  "voom", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                  "expHM", "expHM log", "lnHM", "lnHM log"))
rbind(colMeans(DE.results.DEDD2$auc), colMeans(DE.results.DEDD5$auc), 
      colMeans(DE.results.DEDD10$auc), colMeans(DE.results.DEDD20$auc), 
      colMeans(DE.results.DEDD50$auc))


###########
### pAUC ###
par(mfrow=c(2,3), mar=c(1,2,3,1), mgp=c(3,0.7,0))
boxplot(DE.results.DE2$pauc, names=NA, col=col_vector[1:14], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,0.06), 
        main=paste0("Partial AUC for differential expression, 2 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean only"), 
        cex.main=1.2)
lines(c(14.5,14.5), c(0,0.049))
legend("topleft", fill=col_vector[1:14], bty='n', cex=1.2, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                  "voom", "DSS", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                  "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DE.results.DE5$pauc, names=NA, col=col_vector[1:14], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,0.06), 
        main=paste0("Partial AUC for differential expression, 5 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean only"), 
        cex.main=1.2)
lines(c(14.5,14.5), c(0,0.049))
legend("topleft", fill=col_vector[1:14], bty='n', cex=1.2, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                  "voom", "DSS", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                  "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DE.results.DE10$pauc, names=NA, col=col_vector[1:14], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,0.06), 
        main=paste0("Partial AUC for differential expression, 10 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean only"), 
        cex.main=1.2)
lines(c(14.5,14.5), c(0,0.049))
legend("topleft", fill=col_vector[1:14], bty='n', cex=1.2, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                  "voom", "DSS", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                  "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DE.results.DE20$pauc, names=NA, col=col_vector[1:14], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,0.06), 
        main=paste0("Partial AUC for differential expression, 20 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean only"), 
        cex.main=1.2)
lines(c(14.5,14.5), c(0,0.049))
legend("topleft", fill=col_vector[c(1:6, 8:14)], bty='n', cex=1.2, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                  "voom", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                  "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DE.results.DE50$pauc, names=NA, col=col_vector[c(1:6, 8:14)], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,0.06), 
        main=paste0("Partial AUC for differential expression, 50 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean only"), 
        cex.main=1.2)
lines(c(13.5,13.5), c(0,0.049))
legend("topleft", fill=col_vector[1:14], bty='n', cex=1.2, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                  "voom", "DSS", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                  "expHM", "expHM log", "lnHM", "lnHM log"))
rbind(colMeans(DE.results.DE2$pauc), colMeans(DE.results.DE5$pauc), 
      colMeans(DE.results.DE10$pauc), colMeans(DE.results.DE20$pauc), 
      colMeans(DE.results.DE50$pauc))

par(mfrow=c(2,3), mar=c(1,2,3,1), mgp=c(3,0.7,0))
boxplot(DE.results.DEDD2$pauc, names=NA, col=col_vector[1:14], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,0.06), 
        main=paste0("Partial AUC for differential expression, 2 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean and dispersion"), 
        cex.main=1.2)
lines(c(14.5,14.5), c(0,0.049))
legend("topleft", fill=col_vector[1:14], bty='n', cex=1.2, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                  "voom", "DSS", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                  "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DE.results.DEDD5$pauc, names=NA, col=col_vector[1:14], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,0.06), 
        main=paste0("Partial AUC for differential expression, 5 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean and dispersion"), 
        cex.main=1.2)
lines(c(14.5,14.5), c(0,0.049))
legend("topleft", fill=col_vector[1:14], bty='n', cex=1.2, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                  "voom", "DSS", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                  "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DE.results.DEDD10$pauc, names=NA, col=col_vector[1:14], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,0.06), 
        main=paste0("Partial AUC for differential expression, 10 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean and dispersion"), 
        cex.main=1.2)
lines(c(14.5,14.5), c(0,0.049))
legend("topleft", fill=col_vector[1:14], bty='n', cex=1.2, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                  "voom", "DSS", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                  "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DE.results.DEDD20$pauc, names=NA, col=col_vector[1:14], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,0.06), 
        main=paste0("Partial AUC for differential expression, 20 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean and dispersion"), 
        cex.main=1.2)
lines(c(14.5,14.5), c(0,0.049))
legend("topleft", fill=col_vector[1:14], bty='n', cex=1.2, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                  "voom", "DSS", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                  "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DE.results.DEDD50$pauc, names=NA, col=col_vector[c(1:6, 8:14)], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,0.06), 
        main=paste0("Partial AUC for differential expression, 50 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean and dispersion"), 
        cex.main=1.2)
lines(c(13.5,13.5), c(0,0.049))
legend("topleft", fill=col_vector[c(1:6, 8:14)], bty='n', cex=1.2, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                  "voom", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                  "expHM", "expHM log", "lnHM", "lnHM log"))
rbind(colMeans(DE.results.DEDD2$pauc), colMeans(DE.results.DEDD5$pauc), 
      colMeans(DE.results.DEDD10$pauc), colMeans(DE.results.DEDD20$pauc), 
      colMeans(DE.results.DEDD50$pauc))


###########
### FDR ###
par(mfrow=c(2,3), mar=c(1,2,3,1), mgp=c(3,0.7,0))
boxplot(DE.results.DE2$fdr, names=NA, col=col_vector[1:15], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,1.2), 
        main=paste0("FDR for differential expression, 2 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean only"), 
        cex.main=1.2)
lines(c(15.5,15.5), c(0,1)); abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[1:15], bty='n', cex=1.2, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                  "voom", "DSS", "DSS lfdr", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                  "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DE.results.DE5$fdr, names=NA, col=col_vector[1:15], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,1.2), 
        main=paste0("FDR for differential expression, 5 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean only"), 
        cex.main=1.2)
lines(c(15.5,15.5), c(0,1)); abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[1:15], bty='n', cex=1.2, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                  "voom", "DSS", "DSS lfdr", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                  "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DE.results.DE10$fdr, names=NA, col=col_vector[1:15], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,1.2), 
        main=paste0("FDR for differential expression, 10 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean only"), 
        cex.main=1.2)
lines(c(15.5,15.5), c(0,1)); abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[1:15], bty='n', cex=1.2, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                  "voom", "DSS", "DSS lfdr", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                  "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DE.results.DE20$fdr, names=NA, col=col_vector[1:15], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,1.2), 
        main=paste0("FDR for differential expression, 20 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean only"), 
        cex.main=1.2)
lines(c(15.5,15.5), c(0,1)); abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[1:15], bty='n', cex=1.2, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                  "voom", "DSS", "DSS lfdr", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                  "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DE.results.DE50$fdr, names=NA, col=col_vector[c(1:6, 9:15)], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,1.2), 
        main=paste0("FDR for differential expression, 50 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean only"), 
        cex.main=1.2)
lines(c(13.5,13.5), c(0,1)); abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[c(1:6, 9:15)], bty='n', cex=1.2, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                  "voom", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                  "expHM", "expHM log", "lnHM", "lnHM log"))
rbind(colMeans(DE.results.DE2$fdr), colMeans(DE.results.DE5$fdr), 
      colMeans(DE.results.DE10$fdr), colMeans(DE.results.DE20$fdr), 
      colMeans(DE.results.DE50$fdr))

par(mfrow=c(2,3), mar=c(1,2,3,1), mgp=c(3,0.7,0))
boxplot(DE.results.DEDD2$fdr, names=NA, col=col_vector[1:15], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,1.2), 
        main=paste0("FDR for differential expression, 2 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean and dispersion"), 
        cex.main=1.2)
lines(c(15.5,15.5), c(0,1)); abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[1:15], bty='n', cex=1.2, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                  "voom", "DSS", "DSS lfdr", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                  "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DE.results.DEDD5$fdr, names=NA, col=col_vector[1:15], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,1.2), 
        main=paste0("FDR for differential expression, 5 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean and dispersion"), 
        cex.main=1.2)
lines(c(15.5,15.5), c(0,1)); abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[1:15], bty='n', cex=1.2, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                  "voom", "DSS", "DSS lfdr", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                  "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DE.results.DEDD10$fdr, names=NA, col=col_vector[1:15], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,1.2), 
        main=paste0("FDR for differential expression, 10 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean and dispersion"), 
        cex.main=1.2)
lines(c(15.5,15.5), c(0,1)); abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[1:15], bty='n', cex=1.2, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                  "voom", "DSS", "DSS lfdr", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                  "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DE.results.DEDD20$fdr, names=NA, col=col_vector[1:15], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,1.2), 
        main=paste0("FDR for differential expression, 20 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean and dispersion"), 
        cex.main=1.2)
lines(c(15.5,15.5), c(0,1)); abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[1:15], bty='n', cex=1.2, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                  "voom", "DSS", "DSS lfdr", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                  "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DE.results.DEDD50$fdr, names=NA, col=col_vector[c(1:6, 9:15)], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,1.2), 
        main=paste0("FDR for differential expression, 50 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean and dispersion"), 
        cex.main=1.2)
lines(c(13.5,13.5), c(0,1)); abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[c(1:6, 9:15)], bty='n', cex=1.2, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                  "voom", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                  "expHM", "expHM log", "lnHM", "lnHM log"))
rbind(colMeans(DE.results.DEDD2$fdr), colMeans(DE.results.DEDD5$fdr), 
      colMeans(DE.results.DEDD10$fdr), colMeans(DE.results.DEDD20$fdr), 
      colMeans(DE.results.DEDD50$fdr))


###########
### TPR ###
par(mfrow=c(2,3), mar=c(1,2,3,1), mgp=c(3,0.7,0))
boxplot(DE.results.DE2$tpr, names=NA, col=col_vector[1:15], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,1.25), 
        main=paste0("TPR for differential expression, 2 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean only"), 
        cex.main=1.2)
lines(c(15.5,15.5), c(0,1))
legend("topleft", fill=col_vector[1:15], bty='n', cex=1.2, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                  "voom", "DSS", "DSS lfdr", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                  "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DE.results.DE5$tpr, names=NA, col=col_vector[1:15], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,1.25), 
        main=paste0("TPR for differential expression, 5 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean only"), 
        cex.main=1.2)
lines(c(15.5,15.5), c(0,1))
legend("topleft", fill=col_vector[1:15], bty='n', cex=1.2, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                  "voom", "DSS", "DSS lfdr", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                  "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DE.results.DE10$tpr, names=NA, col=col_vector[1:15], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,1.25), 
        main=paste0("TPR for differential expression, 10 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean only"), 
        cex.main=1.2)
lines(c(15.5,15.5), c(0,1))
legend("topleft", fill=col_vector[1:15], bty='n', cex=1.2, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                  "voom", "DSS", "DSS lfdr", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                  "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DE.results.DE20$tpr, names=NA, col=col_vector[1:15], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,1.25), 
        main=paste0("TPR for differential expression, 20 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean only"), 
        cex.main=1.2)
lines(c(15.5,15.5), c(0,1))
legend("topleft", fill=col_vector[1:15], bty='n', cex=1.2, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                  "voom", "DSS", "DSS lfdr", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                  "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DE.results.DE50$tpr, names=NA, col=col_vector[c(1:6, 9:15)], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,1.25), 
        main=paste0("TPR for differential expression, 50 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean only"), 
        cex.main=1.2)
lines(c(13.5,13.5), c(0,1))
legend("topleft", fill=col_vector[c(1:6, 9:15)], bty='n', cex=1.2, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                  "voom", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                  "expHM", "expHM log", "lnHM", "lnHM log"))
rbind(colMeans(DE.results.DE2$tpr), colMeans(DE.results.DE5$tpr), 
      colMeans(DE.results.DE10$tpr), colMeans(DE.results.DE20$tpr), 
      colMeans(DE.results.DE50$tpr))

par(mfrow=c(2,3), mar=c(1,2,3,1), mgp=c(3,0.7,0))
boxplot(DE.results.DEDD2$tpr, names=NA, col=col_vector[1:15], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,1.25), 
        main=paste0("TPR for differential expression, 2 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean and dispersion"), 
        cex.main=1.2)
lines(c(15.5,15.5), c(0,1))
legend("topleft", fill=col_vector[1:15], bty='n', cex=1.2, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                  "voom", "DSS", "DSS lfdr", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                  "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DE.results.DEDD5$tpr, names=NA, col=col_vector[1:15], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,1.25), 
        main=paste0("TPR for differential expression, 5 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean and dispersion"), 
        cex.main=1.2)
lines(c(15.5,15.5), c(0,1))
legend("topleft", fill=col_vector[1:15], bty='n', cex=1.2, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                  "voom", "DSS", "DSS lfdr", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                  "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DE.results.DEDD10$tpr, names=NA, col=col_vector[1:15], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,1.25), 
        main=paste0("TPR for differential expression, 10 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean and dispersion"), 
        cex.main=1.2)
lines(c(15.5,15.5), c(0,1))
legend("topleft", fill=col_vector[1:15], bty='n', cex=1.2, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                  "voom", "DSS", "DSS lfdr", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                  "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DE.results.DEDD20$tpr, names=NA, col=col_vector[1:15], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,1.25), 
        main=paste0("TPR for differential expression, 20 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean and dispersion"), 
        cex.main=1.2)
lines(c(15.5,15.5), c(0,1))
legend("topleft", fill=col_vector[1:15], bty='n', cex=1.2, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                  "voom", "DSS", "DSS lfdr", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                  "expHM", "expHM log", "lnHM", "lnHM log"))
boxplot(DE.results.DEDD50$tpr, names=NA, col=col_vector[c(1:6, 9:15)], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,1.25), 
        main=paste0("TPR for differential expression, 50 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean and dispersion"), 
        cex.main=1.2)
lines(c(13.5,13.5), c(0,1))
legend("topleft", fill=col_vector[c(1:6, 9:15)], bty='n', cex=1.2, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                  "voom", "baySeq", "MDSeq ZI", "MDSeq no ZI", 
                  "expHM", "expHM log", "lnHM", "lnHM log"))
rbind(colMeans(DE.results.DEDD2$tpr), colMeans(DE.results.DEDD5$tpr), 
      colMeans(DE.results.DEDD10$tpr), colMeans(DE.results.DEDD20$tpr), 
      colMeans(DE.results.DEDD50$tpr))


#############################
### False discovery plots ###
par(mfrow=c(4,5), mar=c(2,2,3,1), mgp=c(3,0.7,0))

plot(DE.results.DE2$mean.discoveries$edgeR.ql.tmm, 
     DE.results.DE2$mean.fdr$edgeR.ql.tmm, 
     type='l', ylim=c(0,0.8), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("2 samples per group\n", "Differences in mean only, ", 
     "TMM"))
lines(DE.results.DE2$mean.discoveries$edgeR.lr.tmm, 
      DE.results.DE2$mean.fdr$edgeR.lr.tmm, 
      type='l', col=col_vector[2])
lines(DE.results.DE2$mean.discoveries$edgeR.et.tmm, 
      DE.results.DE2$mean.fdr$edgeR.et.tmm, 
      type='l', col=col_vector[3])
lines(DE.results.DE2$mean.discoveries$DESeq2.noif.tmm, 
      DE.results.DE2$mean.fdr$DESeq2.noif.tmm, 
      type='l', col=col_vector[4])
lines(DE.results.DE2$mean.discoveries$DESeq2.if.tmm, 
      DE.results.DE2$mean.fdr$DESeq2.if.tmm, 
      type='l', col=col_vector[5])
lines(DE.results.DE2$mean.discoveries$voom.tmm, 
      DE.results.DE2$mean.fdr$voom.tmm, 
      type='l', col=col_vector[6])
lines(DE.results.DE2$mean.discoveries$DSS.tmm, 
      DE.results.DE2$mean.fdr$DSS.tmm, 
      type='l', col=col_vector[7])
lines(DE.results.DE2$mean.discoveries$baySeq.tmm, 
      DE.results.DE2$mean.fdr$baySeq.tmm, 
      type='l', col=col_vector[8])
lines(DE.results.DE2$mean.discoveries$mean.MDSeq.zi.tmm, 
      DE.results.DE2$mean.fdr$mean.MDSeq.zi.tmm, 
      type='l', col=col_vector[9])
lines(DE.results.DE2$mean.discoveries$mean.MDSeq.nozi.tmm, 
      DE.results.DE2$mean.fdr$mean.MDSeq.nozi.tmm, 
      type='l', col=col_vector[10])
lines(DE.results.DE2$mean.discoveries$mean.expHM.untr.tmm, 
      DE.results.DE2$mean.fdr$mean.expHM.untr.tmm, 
      type='l', col=col_vector[11])
lines(DE.results.DE2$mean.discoveries$mean.expHM.log.tmm, 
      DE.results.DE2$mean.fdr$mean.expHM.log.tmm, 
      type='l', col=col_vector[12])
lines(DE.results.DE2$mean.discoveries$mean.lnHM.untr.tmm, 
      DE.results.DE2$mean.fdr$mean.lnHM.untr.tmm, 
      type='l', col=col_vector[13])
lines(DE.results.DE2$mean.discoveries$mean.lnHM.log.tmm, 
      DE.results.DE2$mean.fdr$mean.lnHM.log.tmm, 
      type='l', col=col_vector[14])
legend("bottomright", bty='n', col=col_vector[1:14], lty=1, ncol=2, cex=0.9, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
       "voom", "DSS", "baySeq", "MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", 
       "lnHM", "lnHM log"))
plot(DE.results.DE5$mean.discoveries$edgeR.ql.tmm, 
     DE.results.DE5$mean.fdr$edgeR.ql.tmm, 
     type='l', ylim=c(0,0.8), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("5 samples per group\n", "Differences in mean only, ", 
     "TMM"))
lines(DE.results.DE5$mean.discoveries$edgeR.lr.tmm, 
      DE.results.DE5$mean.fdr$edgeR.lr.tmm, 
      type='l', col=col_vector[2])
lines(DE.results.DE5$mean.discoveries$edgeR.et.tmm, 
      DE.results.DE5$mean.fdr$edgeR.et.tmm, 
      type='l', col=col_vector[3])
lines(DE.results.DE5$mean.discoveries$DESeq2.noif.tmm, 
      DE.results.DE5$mean.fdr$DESeq2.noif.tmm, 
      type='l', col=col_vector[4])
lines(DE.results.DE5$mean.discoveries$DESeq2.if.tmm, 
      DE.results.DE5$mean.fdr$DESeq2.if.tmm, 
      type='l', col=col_vector[5])
lines(DE.results.DE5$mean.discoveries$voom.tmm, 
      DE.results.DE5$mean.fdr$voom.tmm, 
      type='l', col=col_vector[6])
lines(DE.results.DE5$mean.discoveries$DSS.tmm, 
      DE.results.DE5$mean.fdr$DSS.tmm, 
      type='l', col=col_vector[7])
lines(DE.results.DE5$mean.discoveries$baySeq.tmm, 
      DE.results.DE5$mean.fdr$baySeq.tmm, 
      type='l', col=col_vector[8])
lines(DE.results.DE5$mean.discoveries$mean.MDSeq.zi.tmm, 
      DE.results.DE5$mean.fdr$mean.MDSeq.zi.tmm, 
      type='l', col=col_vector[9])
lines(DE.results.DE5$mean.discoveries$mean.MDSeq.nozi.tmm, 
      DE.results.DE5$mean.fdr$mean.MDSeq.nozi.tmm, 
      type='l', col=col_vector[10])
lines(DE.results.DE5$mean.discoveries$mean.expHM.untr.tmm, 
      DE.results.DE5$mean.fdr$mean.expHM.untr.tmm, 
      type='l', col=col_vector[11])
lines(DE.results.DE5$mean.discoveries$mean.expHM.log.tmm, 
      DE.results.DE5$mean.fdr$mean.expHM.log.tmm, 
      type='l', col=col_vector[12])
lines(DE.results.DE5$mean.discoveries$mean.lnHM.untr.tmm, 
      DE.results.DE5$mean.fdr$mean.lnHM.untr.tmm, 
      type='l', col=col_vector[13])
lines(DE.results.DE5$mean.discoveries$mean.lnHM.log.tmm, 
      DE.results.DE5$mean.fdr$mean.lnHM.log.tmm, 
      type='l', col=col_vector[14])
abline(h=0.05, col='lightgrey')
plot(DE.results.DE10$mean.discoveries$edgeR.ql.tmm, 
     DE.results.DE10$mean.fdr$edgeR.ql.tmm, 
     type='l', ylim=c(0,0.8), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("10 samples per group\n", "Differences in mean only, ", 
     "TMM"))
lines(DE.results.DE10$mean.discoveries$edgeR.lr.tmm, 
      DE.results.DE10$mean.fdr$edgeR.lr.tmm, 
      type='l', col=col_vector[2])
lines(DE.results.DE10$mean.discoveries$edgeR.et.tmm, 
      DE.results.DE10$mean.fdr$edgeR.et.tmm, 
      type='l', col=col_vector[3])
lines(DE.results.DE10$mean.discoveries$DESeq2.noif.tmm, 
      DE.results.DE10$mean.fdr$DESeq2.noif.tmm, 
      type='l', col=col_vector[4])
lines(DE.results.DE10$mean.discoveries$DESeq2.if.tmm, 
      DE.results.DE10$mean.fdr$DESeq2.if.tmm, 
      type='l', col=col_vector[5])
lines(DE.results.DE10$mean.discoveries$voom.tmm, 
      DE.results.DE10$mean.fdr$voom.tmm, 
      type='l', col=col_vector[6])
lines(DE.results.DE10$mean.discoveries$DSS.tmm, 
      DE.results.DE10$mean.fdr$DSS.tmm, 
      type='l', col=col_vector[7])
lines(DE.results.DE10$mean.discoveries$baySeq.tmm, 
      DE.results.DE10$mean.fdr$baySeq.tmm, 
      type='l', col=col_vector[8])
lines(DE.results.DE10$mean.discoveries$mean.MDSeq.zi.tmm, 
      DE.results.DE10$mean.fdr$mean.MDSeq.zi.tmm, 
      type='l', col=col_vector[9])
lines(DE.results.DE10$mean.discoveries$mean.MDSeq.nozi.tmm, 
      DE.results.DE10$mean.fdr$mean.MDSeq.nozi.tmm, 
      type='l', col=col_vector[10])
lines(DE.results.DE10$mean.discoveries$mean.expHM.untr.tmm, 
      DE.results.DE10$mean.fdr$mean.expHM.untr.tmm, 
      type='l', col=col_vector[11])
lines(DE.results.DE10$mean.discoveries$mean.expHM.log.tmm, 
      DE.results.DE10$mean.fdr$mean.expHM.log.tmm, 
      type='l', col=col_vector[12])
lines(DE.results.DE10$mean.discoveries$mean.lnHM.untr.tmm, 
      DE.results.DE10$mean.fdr$mean.lnHM.untr.tmm, 
      type='l', col=col_vector[13])
lines(DE.results.DE10$mean.discoveries$mean.lnHM.log.tmm, 
      DE.results.DE10$mean.fdr$mean.lnHM.log.tmm, 
      type='l', col=col_vector[14])
abline(h=0.05, col='lightgrey')
plot(DE.results.DE20$mean.discoveries$edgeR.ql.tmm, 
     DE.results.DE20$mean.fdr$edgeR.ql.tmm, 
     type='l', ylim=c(0,0.8), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("20 samples per group\n", "Differences in mean only, ", 
     "TMM"))
lines(DE.results.DE20$mean.discoveries$edgeR.lr.tmm, 
      DE.results.DE20$mean.fdr$edgeR.lr.tmm, 
      type='l', col=col_vector[2])
lines(DE.results.DE20$mean.discoveries$edgeR.et.tmm, 
      DE.results.DE20$mean.fdr$edgeR.et.tmm, 
      type='l', col=col_vector[3])
lines(DE.results.DE20$mean.discoveries$DESeq2.noif.tmm, 
      DE.results.DE20$mean.fdr$DESeq2.noif.tmm, 
      type='l', col=col_vector[4])
lines(DE.results.DE20$mean.discoveries$DESeq2.if.tmm, 
      DE.results.DE20$mean.fdr$DESeq2.if.tmm, 
      type='l', col=col_vector[5])
lines(DE.results.DE20$mean.discoveries$voom.tmm, 
      DE.results.DE20$mean.fdr$voom.tmm, 
      type='l', col=col_vector[6])
lines(DE.results.DE20$mean.discoveries$DSS.tmm, 
      DE.results.DE20$mean.fdr$DSS.tmm, 
      type='l', col=col_vector[7])
lines(DE.results.DE20$mean.discoveries$baySeq.tmm, 
      DE.results.DE20$mean.fdr$baySeq.tmm, 
      type='l', col=col_vector[8])
lines(DE.results.DE20$mean.discoveries$mean.MDSeq.zi.tmm, 
      DE.results.DE20$mean.fdr$mean.MDSeq.zi.tmm, 
      type='l', col=col_vector[9])
lines(DE.results.DE20$mean.discoveries$mean.MDSeq.nozi.tmm, 
      DE.results.DE20$mean.fdr$mean.MDSeq.nozi.tmm, 
      type='l', col=col_vector[10])
lines(DE.results.DE20$mean.discoveries$mean.expHM.untr.tmm, 
      DE.results.DE20$mean.fdr$mean.expHM.untr.tmm, 
      type='l', col=col_vector[11])
lines(DE.results.DE20$mean.discoveries$mean.expHM.log.tmm, 
      DE.results.DE20$mean.fdr$mean.expHM.log.tmm, 
      type='l', col=col_vector[12])
lines(DE.results.DE20$mean.discoveries$mean.lnHM.untr.tmm, 
      DE.results.DE20$mean.fdr$mean.lnHM.untr.tmm, 
      type='l', col=col_vector[13])
lines(DE.results.DE20$mean.discoveries$mean.lnHM.log.tmm, 
      DE.results.DE20$mean.fdr$mean.lnHM.log.tmm, 
      type='l', col=col_vector[14])
abline(h=0.05, col='lightgrey')
plot(DE.results.DE50$mean.discoveries$edgeR.ql.tmm, 
     DE.results.DE50$mean.fdr$edgeR.ql.tmm, 
     type='l', ylim=c(0,0.8), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("2 samples per group\n", "Differences in mean only, ", 
     "TMM"))
lines(DE.results.DE50$mean.discoveries$edgeR.lr.tmm, 
      DE.results.DE50$mean.fdr$edgeR.lr.tmm, 
      type='l', col=col_vector[2])
lines(DE.results.DE50$mean.discoveries$edgeR.et.tmm, 
      DE.results.DE50$mean.fdr$edgeR.et.tmm, 
      type='l', col=col_vector[3])
lines(DE.results.DE50$mean.discoveries$DESeq2.noif.tmm, 
      DE.results.DE50$mean.fdr$DESeq2.noif.tmm, 
      type='l', col=col_vector[4])
lines(DE.results.DE50$mean.discoveries$DESeq2.if.tmm, 
      DE.results.DE50$mean.fdr$DESeq2.if.tmm, 
      type='l', col=col_vector[5])
lines(DE.results.DE50$mean.discoveries$voom.tmm, 
      DE.results.DE50$mean.fdr$voom.tmm, 
      type='l', col=col_vector[6])
lines(DE.results.DE50$mean.discoveries$DSS.tmm, 
      DE.results.DE50$mean.fdr$DSS.tmm, 
      type='l', col=col_vector[7])
lines(DE.results.DE50$mean.discoveries$baySeq.tmm, 
      DE.results.DE50$mean.fdr$baySeq.tmm, 
      type='l', col=col_vector[8])
lines(DE.results.DE50$mean.discoveries$mean.MDSeq.zi.tmm, 
      DE.results.DE50$mean.fdr$mean.MDSeq.zi.tmm, 
      type='l', col=col_vector[9])
lines(DE.results.DE50$mean.discoveries$mean.MDSeq.nozi.tmm, 
      DE.results.DE50$mean.fdr$mean.MDSeq.nozi.tmm, 
      type='l', col=col_vector[10])
lines(DE.results.DE50$mean.discoveries$mean.expHM.untr.tmm, 
      DE.results.DE50$mean.fdr$mean.expHM.untr.tmm, 
      type='l', col=col_vector[11])
lines(DE.results.DE50$mean.discoveries$mean.expHM.log.tmm, 
      DE.results.DE50$mean.fdr$mean.expHM.log.tmm, 
      type='l', col=col_vector[12])
lines(DE.results.DE50$mean.discoveries$mean.lnHM.untr.tmm, 
      DE.results.DE50$mean.fdr$mean.lnHM.untr.tmm, 
      type='l', col=col_vector[13])
lines(DE.results.DE50$mean.discoveries$mean.lnHM.log.tmm, 
      DE.results.DE50$mean.fdr$mean.lnHM.log.tmm, 
      type='l', col=col_vector[14])
abline(h=0.05, col='lightgrey')

plot(DE.results.DE2$mean.discoveries$edgeR.ql.rle, 
     DE.results.DE2$mean.fdr$edgeR.ql.rle, 
     type='l', ylim=c(0,0.8), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("2 samples per group\n", "Differences in mean only, ", 
     "RLE"))
lines(DE.results.DE2$mean.discoveries$edgeR.lr.rle, 
      DE.results.DE2$mean.fdr$edgeR.lr.rle, 
      type='l', col=col_vector[2])
lines(DE.results.DE2$mean.discoveries$edgeR.et.rle, 
      DE.results.DE2$mean.fdr$edgeR.et.rle, 
      type='l', col=col_vector[3])
lines(DE.results.DE2$mean.discoveries$DESeq2.noif.rle, 
      DE.results.DE2$mean.fdr$DESeq2.noif.rle, 
      type='l', col=col_vector[4])
lines(DE.results.DE2$mean.discoveries$DESeq2.if.rle, 
      DE.results.DE2$mean.fdr$DESeq2.if.rle, 
      type='l', col=col_vector[5])
lines(DE.results.DE2$mean.discoveries$voom.rle, 
      DE.results.DE2$mean.fdr$voom.rle, 
      type='l', col=col_vector[6])
lines(DE.results.DE2$mean.discoveries$DSS.rle, 
      DE.results.DE2$mean.fdr$DSS.rle, 
      type='l', col=col_vector[7])
lines(DE.results.DE2$mean.discoveries$baySeq.rle, 
      DE.results.DE2$mean.fdr$baySeq.rle, 
      type='l', col=col_vector[8])
lines(DE.results.DE2$mean.discoveries$mean.MDSeq.zi.rle, 
      DE.results.DE2$mean.fdr$mean.MDSeq.zi.rle, 
      type='l', col=col_vector[9])
lines(DE.results.DE2$mean.discoveries$mean.MDSeq.nozi.rle, 
      DE.results.DE2$mean.fdr$mean.MDSeq.nozi.rle, 
      type='l', col=col_vector[10])
lines(DE.results.DE2$mean.discoveries$mean.expHM.untr.rle, 
      DE.results.DE2$mean.fdr$mean.expHM.untr.rle, 
      type='l', col=col_vector[11])
lines(DE.results.DE2$mean.discoveries$mean.expHM.log.rle, 
      DE.results.DE2$mean.fdr$mean.expHM.log.rle, 
      type='l', col=col_vector[12])
lines(DE.results.DE2$mean.discoveries$mean.lnHM.untr.rle, 
      DE.results.DE2$mean.fdr$mean.lnHM.untr.rle, 
      type='l', col=col_vector[13])
lines(DE.results.DE2$mean.discoveries$mean.lnHM.log.rle, 
      DE.results.DE2$mean.fdr$mean.lnHM.log.rle, 
      type='l', col=col_vector[14])
legend("bottomright", bty='n', col=col_vector[1:14], lty=1, ncol=2, cex=0.9, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
       "voom", "DSS", "baySeq", "MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", 
       "lnHM", "lnHM log"))
plot(DE.results.DE5$mean.discoveries$edgeR.ql.rle, 
     DE.results.DE5$mean.fdr$edgeR.ql.rle, 
     type='l', ylim=c(0,0.8), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("5 samples per group\n", "Differences in mean only, ", 
     "RLE"))
lines(DE.results.DE5$mean.discoveries$edgeR.lr.rle, 
      DE.results.DE5$mean.fdr$edgeR.lr.rle, 
      type='l', col=col_vector[2])
lines(DE.results.DE5$mean.discoveries$edgeR.et.rle, 
      DE.results.DE5$mean.fdr$edgeR.et.rle, 
      type='l', col=col_vector[3])
lines(DE.results.DE5$mean.discoveries$DESeq2.noif.rle, 
      DE.results.DE5$mean.fdr$DESeq2.noif.rle, 
      type='l', col=col_vector[4])
lines(DE.results.DE5$mean.discoveries$DESeq2.if.rle, 
      DE.results.DE5$mean.fdr$DESeq2.if.rle, 
      type='l', col=col_vector[5])
lines(DE.results.DE5$mean.discoveries$voom.rle, 
      DE.results.DE5$mean.fdr$voom.rle, 
      type='l', col=col_vector[6])
lines(DE.results.DE5$mean.discoveries$DSS.rle, 
      DE.results.DE5$mean.fdr$DSS.rle, 
      type='l', col=col_vector[7])
lines(DE.results.DE5$mean.discoveries$baySeq.rle, 
      DE.results.DE5$mean.fdr$baySeq.rle, 
      type='l', col=col_vector[8])
lines(DE.results.DE5$mean.discoveries$mean.MDSeq.zi.rle, 
      DE.results.DE5$mean.fdr$mean.MDSeq.zi.rle, 
      type='l', col=col_vector[9])
lines(DE.results.DE5$mean.discoveries$mean.MDSeq.nozi.rle, 
      DE.results.DE5$mean.fdr$mean.MDSeq.nozi.rle, 
      type='l', col=col_vector[10])
lines(DE.results.DE5$mean.discoveries$mean.expHM.untr.rle, 
      DE.results.DE5$mean.fdr$mean.expHM.untr.rle, 
      type='l', col=col_vector[11])
lines(DE.results.DE5$mean.discoveries$mean.expHM.log.rle, 
      DE.results.DE5$mean.fdr$mean.expHM.log.rle, 
      type='l', col=col_vector[12])
lines(DE.results.DE5$mean.discoveries$mean.lnHM.untr.rle, 
      DE.results.DE5$mean.fdr$mean.lnHM.untr.rle, 
      type='l', col=col_vector[13])
lines(DE.results.DE5$mean.discoveries$mean.lnHM.log.rle, 
      DE.results.DE5$mean.fdr$mean.lnHM.log.rle, 
      type='l', col=col_vector[14])

abline(h=0.05, col='lightgrey')
plot(DE.results.DE10$mean.discoveries$edgeR.ql.rle, 
     DE.results.DE10$mean.fdr$edgeR.ql.rle, 
     type='l', ylim=c(0,0.8), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("10 samples per group\n", "Differences in mean only, ", 
     "RLE"))
lines(DE.results.DE10$mean.discoveries$edgeR.lr.rle, 
      DE.results.DE10$mean.fdr$edgeR.lr.rle, 
      type='l', col=col_vector[2])
lines(DE.results.DE10$mean.discoveries$edgeR.et.rle, 
      DE.results.DE10$mean.fdr$edgeR.et.rle, 
      type='l', col=col_vector[3])
lines(DE.results.DE10$mean.discoveries$DESeq2.noif.rle, 
      DE.results.DE10$mean.fdr$DESeq2.noif.rle, 
      type='l', col=col_vector[4])
lines(DE.results.DE10$mean.discoveries$DESeq2.if.rle, 
      DE.results.DE10$mean.fdr$DESeq2.if.rle, 
      type='l', col=col_vector[5])
lines(DE.results.DE10$mean.discoveries$voom.rle, 
      DE.results.DE10$mean.fdr$voom.rle, 
      type='l', col=col_vector[6])
lines(DE.results.DE10$mean.discoveries$DSS.rle, 
      DE.results.DE10$mean.fdr$DSS.rle, 
      type='l', col=col_vector[7])
lines(DE.results.DE10$mean.discoveries$baySeq.rle, 
      DE.results.DE10$mean.fdr$baySeq.rle, 
      type='l', col=col_vector[8])
lines(DE.results.DE10$mean.discoveries$mean.MDSeq.zi.rle, 
      DE.results.DE10$mean.fdr$mean.MDSeq.zi.rle, 
      type='l', col=col_vector[9])
lines(DE.results.DE10$mean.discoveries$mean.MDSeq.nozi.rle, 
      DE.results.DE10$mean.fdr$mean.MDSeq.nozi.rle, 
      type='l', col=col_vector[10])
lines(DE.results.DE10$mean.discoveries$mean.expHM.untr.rle, 
      DE.results.DE10$mean.fdr$mean.expHM.untr.rle, 
      type='l', col=col_vector[11])
lines(DE.results.DE10$mean.discoveries$mean.expHM.log.rle, 
      DE.results.DE10$mean.fdr$mean.expHM.log.rle, 
      type='l', col=col_vector[12])
lines(DE.results.DE10$mean.discoveries$mean.lnHM.untr.rle, 
      DE.results.DE10$mean.fdr$mean.lnHM.untr.rle, 
      type='l', col=col_vector[13])
lines(DE.results.DE10$mean.discoveries$mean.lnHM.log.rle, 
      DE.results.DE10$mean.fdr$mean.lnHM.log.rle, 
      type='l', col=col_vector[14])
abline(h=0.05, col='lightgrey')
plot(DE.results.DE20$mean.discoveries$edgeR.ql.rle, 
     DE.results.DE20$mean.fdr$edgeR.ql.rle, 
     type='l', ylim=c(0,0.8), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("20 samples per group\n", "Differences in mean only, ", 
     "RLE"))
lines(DE.results.DE20$mean.discoveries$edgeR.lr.rle, 
      DE.results.DE20$mean.fdr$edgeR.lr.rle, 
      type='l', col=col_vector[2])
lines(DE.results.DE20$mean.discoveries$edgeR.et.rle, 
      DE.results.DE20$mean.fdr$edgeR.et.rle, 
      type='l', col=col_vector[3])
lines(DE.results.DE20$mean.discoveries$DESeq2.noif.rle, 
      DE.results.DE20$mean.fdr$DESeq2.noif.rle, 
      type='l', col=col_vector[4])
lines(DE.results.DE20$mean.discoveries$DESeq2.if.rle, 
      DE.results.DE20$mean.fdr$DESeq2.if.rle, 
      type='l', col=col_vector[5])
lines(DE.results.DE20$mean.discoveries$voom.rle, 
      DE.results.DE20$mean.fdr$voom.rle, 
      type='l', col=col_vector[6])
lines(DE.results.DE20$mean.discoveries$DSS.rle, 
      DE.results.DE20$mean.fdr$DSS.rle, 
      type='l', col=col_vector[7])
lines(DE.results.DE20$mean.discoveries$baySeq.rle, 
      DE.results.DE20$mean.fdr$baySeq.rle, 
      type='l', col=col_vector[8])
lines(DE.results.DE20$mean.discoveries$mean.MDSeq.zi.rle, 
      DE.results.DE20$mean.fdr$mean.MDSeq.zi.rle, 
      type='l', col=col_vector[9])
lines(DE.results.DE20$mean.discoveries$mean.MDSeq.nozi.rle, 
      DE.results.DE20$mean.fdr$mean.MDSeq.nozi.rle, 
      type='l', col=col_vector[10])
lines(DE.results.DE20$mean.discoveries$mean.expHM.untr.rle, 
      DE.results.DE20$mean.fdr$mean.expHM.untr.rle, 
      type='l', col=col_vector[11])
lines(DE.results.DE20$mean.discoveries$mean.expHM.log.rle, 
      DE.results.DE20$mean.fdr$mean.expHM.log.rle, 
      type='l', col=col_vector[12])
lines(DE.results.DE20$mean.discoveries$mean.lnHM.untr.rle, 
      DE.results.DE20$mean.fdr$mean.lnHM.untr.rle, 
      type='l', col=col_vector[13])
lines(DE.results.DE20$mean.discoveries$mean.lnHM.log.rle, 
      DE.results.DE20$mean.fdr$mean.lnHM.log.rle, 
      type='l', col=col_vector[14])
abline(h=0.05, col='lightgrey')
plot(DE.results.DE50$mean.discoveries$edgeR.ql.rle, 
     DE.results.DE50$mean.fdr$edgeR.ql.rle, 
     type='l', ylim=c(0,0.8), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("2 samples per group\n", "Differences in mean only, ", 
     "RLE"))
lines(DE.results.DE50$mean.discoveries$edgeR.lr.rle, 
      DE.results.DE50$mean.fdr$edgeR.lr.rle, 
      type='l', col=col_vector[2])
lines(DE.results.DE50$mean.discoveries$edgeR.et.rle, 
      DE.results.DE50$mean.fdr$edgeR.et.rle, 
      type='l', col=col_vector[3])
lines(DE.results.DE50$mean.discoveries$DESeq2.noif.rle, 
      DE.results.DE50$mean.fdr$DESeq2.noif.rle, 
      type='l', col=col_vector[4])
lines(DE.results.DE50$mean.discoveries$DESeq2.if.rle, 
      DE.results.DE50$mean.fdr$DESeq2.if.rle, 
      type='l', col=col_vector[5])
lines(DE.results.DE50$mean.discoveries$voom.rle, 
      DE.results.DE50$mean.fdr$voom.rle, 
      type='l', col=col_vector[6])
lines(DE.results.DE50$mean.discoveries$DSS.rle, 
      DE.results.DE50$mean.fdr$DSS.rle, 
      type='l', col=col_vector[7])
lines(DE.results.DE50$mean.discoveries$baySeq.rle, 
      DE.results.DE50$mean.fdr$baySeq.rle, 
      type='l', col=col_vector[8])
lines(DE.results.DE50$mean.discoveries$mean.MDSeq.zi.rle, 
      DE.results.DE50$mean.fdr$mean.MDSeq.zi.rle, 
      type='l', col=col_vector[9])
lines(DE.results.DE50$mean.discoveries$mean.MDSeq.nozi.rle, 
      DE.results.DE50$mean.fdr$mean.MDSeq.nozi.rle, 
      type='l', col=col_vector[10])
lines(DE.results.DE50$mean.discoveries$mean.expHM.untr.rle, 
      DE.results.DE50$mean.fdr$mean.expHM.untr.rle, 
      type='l', col=col_vector[11])
lines(DE.results.DE50$mean.discoveries$mean.expHM.log.rle, 
      DE.results.DE50$mean.fdr$mean.expHM.log.rle, 
      type='l', col=col_vector[12])
lines(DE.results.DE50$mean.discoveries$mean.lnHM.untr.rle, 
      DE.results.DE50$mean.fdr$mean.lnHM.untr.rle, 
      type='l', col=col_vector[13])
lines(DE.results.DE50$mean.discoveries$mean.lnHM.log.rle, 
      DE.results.DE50$mean.fdr$mean.lnHM.log.rle, 
      type='l', col=col_vector[14])
abline(h=0.05, col='lightgrey')

plot(DE.results.DEDD2$mean.discoveries$edgeR.ql.tmm, 
     DE.results.DEDD2$mean.fdr$edgeR.ql.tmm, 
     type='l', ylim=c(0,0.8), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("2 samples per group\n", "Differences in mean & disp, ", 
     "TMM"))
lines(DE.results.DEDD2$mean.discoveries$edgeR.lr.tmm, 
      DE.results.DEDD2$mean.fdr$edgeR.lr.tmm, 
      type='l', col=col_vector[2])
lines(DE.results.DEDD2$mean.discoveries$edgeR.et.tmm, 
      DE.results.DEDD2$mean.fdr$edgeR.et.tmm, 
      type='l', col=col_vector[3])
lines(DE.results.DEDD2$mean.discoveries$DESeq2.noif.tmm, 
      DE.results.DEDD2$mean.fdr$DESeq2.noif.tmm, 
      type='l', col=col_vector[4])
lines(DE.results.DEDD2$mean.discoveries$DESeq2.if.tmm, 
      DE.results.DEDD2$mean.fdr$DESeq2.if.tmm, 
      type='l', col=col_vector[5])
lines(DE.results.DEDD2$mean.discoveries$voom.tmm, 
      DE.results.DEDD2$mean.fdr$voom.tmm, 
      type='l', col=col_vector[6])
lines(DE.results.DEDD2$mean.discoveries$DSS.tmm, 
      DE.results.DEDD2$mean.fdr$DSS.tmm, 
      type='l', col=col_vector[7])
lines(DE.results.DEDD2$mean.discoveries$baySeq.tmm, 
      DE.results.DEDD2$mean.fdr$baySeq.tmm, 
      type='l', col=col_vector[8])
lines(DE.results.DEDD2$mean.discoveries$mean.MDSeq.zi.tmm, 
      DE.results.DEDD2$mean.fdr$mean.MDSeq.zi.tmm, 
      type='l', col=col_vector[9])
lines(DE.results.DEDD2$mean.discoveries$mean.MDSeq.nozi.tmm, 
      DE.results.DEDD2$mean.fdr$mean.MDSeq.nozi.tmm, 
      type='l', col=col_vector[10])
lines(DE.results.DEDD2$mean.discoveries$mean.expHM.untr.tmm, 
      DE.results.DEDD2$mean.fdr$mean.expHM.untr.tmm, 
      type='l', col=col_vector[11])
lines(DE.results.DEDD2$mean.discoveries$mean.expHM.log.tmm, 
      DE.results.DEDD2$mean.fdr$mean.expHM.log.tmm, 
      type='l', col=col_vector[12])
lines(DE.results.DEDD2$mean.discoveries$mean.lnHM.untr.tmm, 
      DE.results.DEDD2$mean.fdr$mean.lnHM.untr.tmm, 
      type='l', col=col_vector[13])
lines(DE.results.DEDD2$mean.discoveries$mean.lnHM.log.tmm, 
      DE.results.DEDD2$mean.fdr$mean.lnHM.log.tmm, 
      type='l', col=col_vector[14])
plot(DE.results.DEDD5$mean.discoveries$edgeR.ql.tmm, 
     DE.results.DEDD5$mean.fdr$edgeR.ql.tmm, 
     type='l', ylim=c(0,0.8), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("5 samples per group\n", "Differences in mean & disp, ", 
     "TMM"))
lines(DE.results.DEDD5$mean.discoveries$edgeR.lr.tmm, 
      DE.results.DEDD5$mean.fdr$edgeR.lr.tmm, 
      type='l', col=col_vector[2])
lines(DE.results.DEDD5$mean.discoveries$edgeR.et.tmm, 
      DE.results.DEDD5$mean.fdr$edgeR.et.tmm, 
      type='l', col=col_vector[3])
lines(DE.results.DEDD5$mean.discoveries$DESeq2.noif.tmm, 
      DE.results.DEDD5$mean.fdr$DESeq2.noif.tmm, 
      type='l', col=col_vector[4])
lines(DE.results.DEDD5$mean.discoveries$DESeq2.if.tmm, 
      DE.results.DEDD5$mean.fdr$DESeq2.if.tmm, 
      type='l', col=col_vector[5])
lines(DE.results.DEDD5$mean.discoveries$voom.tmm, 
      DE.results.DEDD5$mean.fdr$voom.tmm, 
      type='l', col=col_vector[6])
lines(DE.results.DEDD5$mean.discoveries$DSS.tmm, 
      DE.results.DEDD5$mean.fdr$DSS.tmm, 
      type='l', col=col_vector[7])
lines(DE.results.DEDD5$mean.discoveries$baySeq.tmm, 
      DE.results.DEDD5$mean.fdr$baySeq.tmm, 
      type='l', col=col_vector[8])
lines(DE.results.DEDD5$mean.discoveries$mean.MDSeq.zi.tmm, 
      DE.results.DEDD5$mean.fdr$mean.MDSeq.zi.tmm, 
      type='l', col=col_vector[9])
lines(DE.results.DEDD5$mean.discoveries$mean.MDSeq.nozi.tmm, 
      DE.results.DEDD5$mean.fdr$mean.MDSeq.nozi.tmm, 
      type='l', col=col_vector[10])
lines(DE.results.DEDD5$mean.discoveries$mean.expHM.untr.tmm, 
      DE.results.DEDD5$mean.fdr$mean.expHM.untr.tmm, 
      type='l', col=col_vector[11])
lines(DE.results.DEDD5$mean.discoveries$mean.expHM.log.tmm, 
      DE.results.DEDD5$mean.fdr$mean.expHM.log.tmm, 
      type='l', col=col_vector[12])
lines(DE.results.DEDD5$mean.discoveries$mean.lnHM.untr.tmm, 
      DE.results.DEDD5$mean.fdr$mean.lnHM.untr.tmm, 
      type='l', col=col_vector[13])
lines(DE.results.DEDD5$mean.discoveries$mean.lnHM.log.tmm, 
      DE.results.DEDD5$mean.fdr$mean.lnHM.log.tmm, 
      type='l', col=col_vector[14])
abline(h=0.05, col='lightgrey')
plot(DE.results.DEDD10$mean.discoveries$edgeR.ql.tmm, 
     DE.results.DEDD10$mean.fdr$edgeR.ql.tmm, 
     type='l', ylim=c(0,0.8), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("10 samples per group\n", "Differences in mean & disp, ", 
     "TMM"))
lines(DE.results.DEDD10$mean.discoveries$edgeR.lr.tmm, 
      DE.results.DEDD10$mean.fdr$edgeR.lr.tmm, 
      type='l', col=col_vector[2])
lines(DE.results.DEDD10$mean.discoveries$edgeR.et.tmm, 
      DE.results.DEDD10$mean.fdr$edgeR.et.tmm, 
      type='l', col=col_vector[3])
lines(DE.results.DEDD10$mean.discoveries$DESeq2.noif.tmm, 
      DE.results.DEDD10$mean.fdr$DESeq2.noif.tmm, 
      type='l', col=col_vector[4])
lines(DE.results.DEDD10$mean.discoveries$DESeq2.if.tmm, 
      DE.results.DEDD10$mean.fdr$DESeq2.if.tmm, 
      type='l', col=col_vector[5])
lines(DE.results.DEDD10$mean.discoveries$voom.tmm, 
      DE.results.DEDD10$mean.fdr$voom.tmm, 
      type='l', col=col_vector[6])
lines(DE.results.DEDD10$mean.discoveries$DSS.tmm, 
      DE.results.DEDD10$mean.fdr$DSS.tmm, 
      type='l', col=col_vector[7])
lines(DE.results.DEDD10$mean.discoveries$baySeq.tmm, 
      DE.results.DEDD10$mean.fdr$baySeq.tmm, 
      type='l', col=col_vector[8])
lines(DE.results.DEDD10$mean.discoveries$mean.MDSeq.zi.tmm, 
      DE.results.DEDD10$mean.fdr$mean.MDSeq.zi.tmm, 
      type='l', col=col_vector[9])
lines(DE.results.DEDD10$mean.discoveries$mean.MDSeq.nozi.tmm, 
      DE.results.DEDD10$mean.fdr$mean.MDSeq.nozi.tmm, 
      type='l', col=col_vector[10])
lines(DE.results.DEDD10$mean.discoveries$mean.expHM.untr.tmm, 
      DE.results.DEDD10$mean.fdr$mean.expHM.untr.tmm, 
      type='l', col=col_vector[11])
lines(DE.results.DEDD10$mean.discoveries$mean.expHM.log.tmm, 
      DE.results.DEDD10$mean.fdr$mean.expHM.log.tmm, 
      type='l', col=col_vector[12])
lines(DE.results.DEDD10$mean.discoveries$mean.lnHM.untr.tmm, 
      DE.results.DEDD10$mean.fdr$mean.lnHM.untr.tmm, 
      type='l', col=col_vector[13])
lines(DE.results.DEDD10$mean.discoveries$mean.lnHM.log.tmm, 
      DE.results.DEDD10$mean.fdr$mean.lnHM.log.tmm, 
      type='l', col=col_vector[14])
abline(h=0.05, col='lightgrey')
plot(DE.results.DEDD20$mean.discoveries$edgeR.ql.tmm, 
     DE.results.DEDD20$mean.fdr$edgeR.ql.tmm, 
     type='l', ylim=c(0,0.8), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("20 samples per group\n", "Differences in mean & disp, ", 
     "TMM"))
lines(DE.results.DEDD20$mean.discoveries$edgeR.lr.tmm, 
      DE.results.DEDD20$mean.fdr$edgeR.lr.tmm, 
      type='l', col=col_vector[2])
lines(DE.results.DEDD20$mean.discoveries$edgeR.et.tmm, 
      DE.results.DEDD20$mean.fdr$edgeR.et.tmm, 
      type='l', col=col_vector[3])
lines(DE.results.DEDD20$mean.discoveries$DESeq2.noif.tmm, 
      DE.results.DEDD20$mean.fdr$DESeq2.noif.tmm, 
      type='l', col=col_vector[4])
lines(DE.results.DEDD20$mean.discoveries$DESeq2.if.tmm, 
      DE.results.DEDD20$mean.fdr$DESeq2.if.tmm, 
      type='l', col=col_vector[5])
lines(DE.results.DEDD20$mean.discoveries$voom.tmm, 
      DE.results.DEDD20$mean.fdr$voom.tmm, 
      type='l', col=col_vector[6])
lines(DE.results.DEDD20$mean.discoveries$DSS.tmm, 
      DE.results.DEDD20$mean.fdr$DSS.tmm, 
      type='l', col=col_vector[7])
lines(DE.results.DEDD20$mean.discoveries$baySeq.tmm, 
      DE.results.DEDD20$mean.fdr$baySeq.tmm, 
      type='l', col=col_vector[8])
lines(DE.results.DEDD20$mean.discoveries$mean.MDSeq.zi.tmm, 
      DE.results.DEDD20$mean.fdr$mean.MDSeq.zi.tmm, 
      type='l', col=col_vector[9])
lines(DE.results.DEDD20$mean.discoveries$mean.MDSeq.nozi.tmm, 
      DE.results.DEDD20$mean.fdr$mean.MDSeq.nozi.tmm, 
      type='l', col=col_vector[10])
lines(DE.results.DEDD20$mean.discoveries$mean.expHM.untr.tmm, 
      DE.results.DEDD20$mean.fdr$mean.expHM.untr.tmm, 
      type='l', col=col_vector[11])
lines(DE.results.DEDD20$mean.discoveries$mean.expHM.log.tmm, 
      DE.results.DEDD20$mean.fdr$mean.expHM.log.tmm, 
      type='l', col=col_vector[12])
lines(DE.results.DEDD20$mean.discoveries$mean.lnHM.untr.tmm, 
      DE.results.DEDD20$mean.fdr$mean.lnHM.untr.tmm, 
      type='l', col=col_vector[13])
lines(DE.results.DEDD20$mean.discoveries$mean.lnHM.log.tmm, 
      DE.results.DEDD20$mean.fdr$mean.lnHM.log.tmm, 
      type='l', col=col_vector[14])
abline(h=0.05, col='lightgrey')
plot(DE.results.DEDD50$mean.discoveries$edgeR.ql.tmm, 
     DE.results.DEDD50$mean.fdr$edgeR.ql.tmm, 
     type='l', ylim=c(0,0.8), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("2 samples per group\n", "Differences in mean & disp, ", 
     "TMM"))
lines(DE.results.DEDD50$mean.discoveries$edgeR.lr.tmm, 
      DE.results.DEDD50$mean.fdr$edgeR.lr.tmm, 
      type='l', col=col_vector[2])
lines(DE.results.DEDD50$mean.discoveries$edgeR.et.tmm, 
      DE.results.DEDD50$mean.fdr$edgeR.et.tmm, 
      type='l', col=col_vector[3])
lines(DE.results.DEDD50$mean.discoveries$DESeq2.noif.tmm, 
      DE.results.DEDD50$mean.fdr$DESeq2.noif.tmm, 
      type='l', col=col_vector[4])
lines(DE.results.DEDD50$mean.discoveries$DESeq2.if.tmm, 
      DE.results.DEDD50$mean.fdr$DESeq2.if.tmm, 
      type='l', col=col_vector[5])
lines(DE.results.DEDD50$mean.discoveries$voom.tmm, 
      DE.results.DEDD50$mean.fdr$voom.tmm, 
      type='l', col=col_vector[6])
lines(DE.results.DEDD50$mean.discoveries$DSS.tmm, 
      DE.results.DEDD50$mean.fdr$DSS.tmm, 
      type='l', col=col_vector[7])
lines(DE.results.DEDD50$mean.discoveries$baySeq.tmm, 
      DE.results.DEDD50$mean.fdr$baySeq.tmm, 
      type='l', col=col_vector[8])
lines(DE.results.DEDD50$mean.discoveries$mean.MDSeq.zi.tmm, 
      DE.results.DEDD50$mean.fdr$mean.MDSeq.zi.tmm, 
      type='l', col=col_vector[9])
lines(DE.results.DEDD50$mean.discoveries$mean.MDSeq.nozi.tmm, 
      DE.results.DEDD50$mean.fdr$mean.MDSeq.nozi.tmm, 
      type='l', col=col_vector[10])
lines(DE.results.DEDD50$mean.discoveries$mean.expHM.untr.tmm, 
      DE.results.DEDD50$mean.fdr$mean.expHM.untr.tmm, 
      type='l', col=col_vector[11])
lines(DE.results.DEDD50$mean.discoveries$mean.expHM.log.tmm, 
      DE.results.DEDD50$mean.fdr$mean.expHM.log.tmm, 
      type='l', col=col_vector[12])
lines(DE.results.DEDD50$mean.discoveries$mean.lnHM.untr.tmm, 
      DE.results.DEDD50$mean.fdr$mean.lnHM.untr.tmm, 
      type='l', col=col_vector[13])
lines(DE.results.DEDD50$mean.discoveries$mean.lnHM.log.tmm, 
      DE.results.DEDD50$mean.fdr$mean.lnHM.log.tmm, 
      type='l', col=col_vector[14])
legend("topleft", bty='n', col=col_vector[1:14], lty=1, ncol=2, cex=0.9, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                "voom", "DSS", "baySeq", "MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", 
                "lnHM", "lnHM log"))
abline(h=0.05, col='lightgrey')

plot(DE.results.DEDD2$mean.discoveries$edgeR.ql.rle, 
     DE.results.DEDD2$mean.fdr$edgeR.ql.rle, 
     type='l', ylim=c(0,0.8), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("2 samples per group\n", "Differences in mean & disp, ", 
     "RLE"))
lines(DE.results.DEDD2$mean.discoveries$edgeR.lr.rle, 
      DE.results.DEDD2$mean.fdr$edgeR.lr.rle, 
      type='l', col=col_vector[2])
lines(DE.results.DEDD2$mean.discoveries$edgeR.et.rle, 
      DE.results.DEDD2$mean.fdr$edgeR.et.rle, 
      type='l', col=col_vector[3])
lines(DE.results.DEDD2$mean.discoveries$DESeq2.noif.rle, 
      DE.results.DEDD2$mean.fdr$DESeq2.noif.rle, 
      type='l', col=col_vector[4])
lines(DE.results.DEDD2$mean.discoveries$DESeq2.if.rle, 
      DE.results.DEDD2$mean.fdr$DESeq2.if.rle, 
      type='l', col=col_vector[5])
lines(DE.results.DEDD2$mean.discoveries$voom.rle, 
      DE.results.DEDD2$mean.fdr$voom.rle, 
      type='l', col=col_vector[6])
lines(DE.results.DEDD2$mean.discoveries$DSS.rle, 
      DE.results.DEDD2$mean.fdr$DSS.rle, 
      type='l', col=col_vector[7])
lines(DE.results.DEDD2$mean.discoveries$baySeq.rle, 
      DE.results.DEDD2$mean.fdr$baySeq.rle, 
      type='l', col=col_vector[8])
lines(DE.results.DEDD2$mean.discoveries$mean.MDSeq.zi.rle, 
      DE.results.DEDD2$mean.fdr$mean.MDSeq.zi.rle, 
      type='l', col=col_vector[9])
lines(DE.results.DEDD2$mean.discoveries$mean.MDSeq.nozi.rle, 
      DE.results.DEDD2$mean.fdr$mean.MDSeq.nozi.rle, 
      type='l', col=col_vector[10])
lines(DE.results.DEDD2$mean.discoveries$mean.expHM.untr.rle, 
      DE.results.DEDD2$mean.fdr$mean.expHM.untr.rle, 
      type='l', col=col_vector[11])
lines(DE.results.DEDD2$mean.discoveries$mean.expHM.log.rle, 
      DE.results.DEDD2$mean.fdr$mean.expHM.log.rle, 
      type='l', col=col_vector[12])
lines(DE.results.DEDD2$mean.discoveries$mean.lnHM.untr.rle, 
      DE.results.DEDD2$mean.fdr$mean.lnHM.untr.rle, 
      type='l', col=col_vector[13])
lines(DE.results.DEDD2$mean.discoveries$mean.lnHM.log.rle, 
      DE.results.DEDD2$mean.fdr$mean.lnHM.log.rle, 
      type='l', col=col_vector[14])
plot(DE.results.DEDD5$mean.discoveries$edgeR.ql.rle, 
     DE.results.DEDD5$mean.fdr$edgeR.ql.rle, 
     type='l', ylim=c(0,0.8), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("5 samples per group\n", "Differences in mean & disp, ", 
     "RLE"))
lines(DE.results.DEDD5$mean.discoveries$edgeR.lr.rle, 
      DE.results.DEDD5$mean.fdr$edgeR.lr.rle, 
      type='l', col=col_vector[2])
lines(DE.results.DEDD5$mean.discoveries$edgeR.et.rle, 
      DE.results.DEDD5$mean.fdr$edgeR.et.rle, 
      type='l', col=col_vector[3])
lines(DE.results.DEDD5$mean.discoveries$DESeq2.noif.rle, 
      DE.results.DEDD5$mean.fdr$DESeq2.noif.rle, 
      type='l', col=col_vector[4])
lines(DE.results.DEDD5$mean.discoveries$DESeq2.if.rle, 
      DE.results.DEDD5$mean.fdr$DESeq2.if.rle, 
      type='l', col=col_vector[5])
lines(DE.results.DEDD5$mean.discoveries$voom.rle, 
      DE.results.DEDD5$mean.fdr$voom.rle, 
      type='l', col=col_vector[6])
lines(DE.results.DEDD5$mean.discoveries$DSS.rle, 
      DE.results.DEDD5$mean.fdr$DSS.rle, 
      type='l', col=col_vector[7])
lines(DE.results.DEDD5$mean.discoveries$baySeq.rle, 
      DE.results.DEDD5$mean.fdr$baySeq.rle, 
      type='l', col=col_vector[8])
lines(DE.results.DEDD5$mean.discoveries$mean.MDSeq.zi.rle, 
      DE.results.DEDD5$mean.fdr$mean.MDSeq.zi.rle, 
      type='l', col=col_vector[9])
lines(DE.results.DEDD5$mean.discoveries$mean.MDSeq.nozi.rle, 
      DE.results.DEDD5$mean.fdr$mean.MDSeq.nozi.rle, 
      type='l', col=col_vector[10])
lines(DE.results.DEDD5$mean.discoveries$mean.expHM.untr.rle, 
      DE.results.DEDD5$mean.fdr$mean.expHM.untr.rle, 
      type='l', col=col_vector[11])
lines(DE.results.DEDD5$mean.discoveries$mean.expHM.log.rle, 
      DE.results.DEDD5$mean.fdr$mean.expHM.log.rle, 
      type='l', col=col_vector[12])
lines(DE.results.DEDD5$mean.discoveries$mean.lnHM.untr.rle, 
      DE.results.DEDD5$mean.fdr$mean.lnHM.untr.rle, 
      type='l', col=col_vector[13])
lines(DE.results.DEDD5$mean.discoveries$mean.lnHM.log.rle, 
      DE.results.DEDD5$mean.fdr$mean.lnHM.log.rle, 
      type='l', col=col_vector[14])
abline(h=0.05, col='lightgrey')
plot(DE.results.DEDD10$mean.discoveries$edgeR.ql.rle, 
     DE.results.DEDD10$mean.fdr$edgeR.ql.rle, 
     type='l', ylim=c(0,0.8), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("10 samples per group\n", "Differences in mean & disp, ", 
     "RLE"))
lines(DE.results.DEDD10$mean.discoveries$edgeR.lr.rle, 
      DE.results.DEDD10$mean.fdr$edgeR.lr.rle, 
      type='l', col=col_vector[2])
lines(DE.results.DEDD10$mean.discoveries$edgeR.et.rle, 
      DE.results.DEDD10$mean.fdr$edgeR.et.rle, 
      type='l', col=col_vector[3])
lines(DE.results.DEDD10$mean.discoveries$DESeq2.noif.rle, 
      DE.results.DEDD10$mean.fdr$DESeq2.noif.rle, 
      type='l', col=col_vector[4])
lines(DE.results.DEDD10$mean.discoveries$DESeq2.if.rle, 
      DE.results.DEDD10$mean.fdr$DESeq2.if.rle, 
      type='l', col=col_vector[5])
lines(DE.results.DEDD10$mean.discoveries$voom.rle, 
      DE.results.DEDD10$mean.fdr$voom.rle, 
      type='l', col=col_vector[6])
lines(DE.results.DEDD10$mean.discoveries$DSS.rle, 
      DE.results.DEDD10$mean.fdr$DSS.rle, 
      type='l', col=col_vector[7])
lines(DE.results.DEDD10$mean.discoveries$baySeq.rle, 
      DE.results.DEDD10$mean.fdr$baySeq.rle, 
      type='l', col=col_vector[8])
lines(DE.results.DEDD10$mean.discoveries$mean.MDSeq.zi.rle, 
      DE.results.DEDD10$mean.fdr$mean.MDSeq.zi.rle, 
      type='l', col=col_vector[9])
lines(DE.results.DEDD10$mean.discoveries$mean.MDSeq.nozi.rle, 
      DE.results.DEDD10$mean.fdr$mean.MDSeq.nozi.rle, 
      type='l', col=col_vector[10])
lines(DE.results.DEDD10$mean.discoveries$mean.expHM.untr.rle, 
      DE.results.DEDD10$mean.fdr$mean.expHM.untr.rle, 
      type='l', col=col_vector[11])
lines(DE.results.DEDD10$mean.discoveries$mean.expHM.log.rle, 
      DE.results.DEDD10$mean.fdr$mean.expHM.log.rle, 
      type='l', col=col_vector[12])
lines(DE.results.DEDD10$mean.discoveries$mean.lnHM.untr.rle, 
      DE.results.DEDD10$mean.fdr$mean.lnHM.untr.rle, 
      type='l', col=col_vector[13])
lines(DE.results.DEDD10$mean.discoveries$mean.lnHM.log.rle, 
      DE.results.DEDD10$mean.fdr$mean.lnHM.log.rle, 
      type='l', col=col_vector[14])
abline(h=0.05, col='lightgrey')
plot(DE.results.DEDD20$mean.discoveries$edgeR.ql.rle, 
     DE.results.DEDD20$mean.fdr$edgeR.ql.rle, 
     type='l', ylim=c(0,0.8), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("20 samples per group\n", "Differences in mean & disp, ", 
     "RLE"))
lines(DE.results.DEDD20$mean.discoveries$edgeR.lr.rle, 
      DE.results.DEDD20$mean.fdr$edgeR.lr.rle, 
      type='l', col=col_vector[2])
lines(DE.results.DEDD20$mean.discoveries$edgeR.et.rle, 
      DE.results.DEDD20$mean.fdr$edgeR.et.rle, 
      type='l', col=col_vector[3])
lines(DE.results.DEDD20$mean.discoveries$DESeq2.noif.rle, 
      DE.results.DEDD20$mean.fdr$DESeq2.noif.rle, 
      type='l', col=col_vector[4])
lines(DE.results.DEDD20$mean.discoveries$DESeq2.if.rle, 
      DE.results.DEDD20$mean.fdr$DESeq2.if.rle, 
      type='l', col=col_vector[5])
lines(DE.results.DEDD20$mean.discoveries$voom.rle, 
      DE.results.DEDD20$mean.fdr$voom.rle, 
      type='l', col=col_vector[6])
lines(DE.results.DEDD20$mean.discoveries$DSS.rle, 
      DE.results.DEDD20$mean.fdr$DSS.rle, 
      type='l', col=col_vector[7])
lines(DE.results.DEDD20$mean.discoveries$baySeq.rle, 
      DE.results.DEDD20$mean.fdr$baySeq.rle, 
      type='l', col=col_vector[8])
lines(DE.results.DEDD20$mean.discoveries$mean.MDSeq.zi.rle, 
      DE.results.DEDD20$mean.fdr$mean.MDSeq.zi.rle, 
      type='l', col=col_vector[9])
lines(DE.results.DEDD20$mean.discoveries$mean.MDSeq.nozi.rle, 
      DE.results.DEDD20$mean.fdr$mean.MDSeq.nozi.rle, 
      type='l', col=col_vector[10])
lines(DE.results.DEDD20$mean.discoveries$mean.expHM.untr.rle, 
      DE.results.DEDD20$mean.fdr$mean.expHM.untr.rle, 
      type='l', col=col_vector[11])
lines(DE.results.DEDD20$mean.discoveries$mean.expHM.log.rle, 
      DE.results.DEDD20$mean.fdr$mean.expHM.log.rle, 
      type='l', col=col_vector[12])
lines(DE.results.DEDD20$mean.discoveries$mean.lnHM.untr.rle, 
      DE.results.DEDD20$mean.fdr$mean.lnHM.untr.rle, 
      type='l', col=col_vector[13])
lines(DE.results.DEDD20$mean.discoveries$mean.lnHM.log.rle, 
      DE.results.DEDD20$mean.fdr$mean.lnHM.log.rle, 
      type='l', col=col_vector[14])
abline(h=0.05, col='lightgrey')
plot(DE.results.DEDD50$mean.discoveries$edgeR.ql.rle, 
     DE.results.DEDD50$mean.fdr$edgeR.ql.rle, 
     type='l', ylim=c(0,0.8), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("2 samples per group\n", "Differences in mean & disp, ", 
     "RLE"))
lines(DE.results.DEDD50$mean.discoveries$edgeR.lr.rle, 
      DE.results.DEDD50$mean.fdr$edgeR.lr.rle, 
      type='l', col=col_vector[2])
lines(DE.results.DEDD50$mean.discoveries$edgeR.et.rle, 
      DE.results.DEDD50$mean.fdr$edgeR.et.rle, 
      type='l', col=col_vector[3])
lines(DE.results.DEDD50$mean.discoveries$DESeq2.noif.rle, 
      DE.results.DEDD50$mean.fdr$DESeq2.noif.rle, 
      type='l', col=col_vector[4])
lines(DE.results.DEDD50$mean.discoveries$DESeq2.if.rle, 
      DE.results.DEDD50$mean.fdr$DESeq2.if.rle, 
      type='l', col=col_vector[5])
lines(DE.results.DEDD50$mean.discoveries$voom.rle, 
      DE.results.DEDD50$mean.fdr$voom.rle, 
      type='l', col=col_vector[6])
lines(DE.results.DEDD50$mean.discoveries$DSS.rle, 
      DE.results.DEDD50$mean.fdr$DSS.rle, 
      type='l', col=col_vector[7])
lines(DE.results.DEDD50$mean.discoveries$baySeq.rle, 
      DE.results.DEDD50$mean.fdr$baySeq.rle, 
      type='l', col=col_vector[8])
lines(DE.results.DEDD50$mean.discoveries$mean.MDSeq.zi.rle, 
      DE.results.DEDD50$mean.fdr$mean.MDSeq.zi.rle, 
      type='l', col=col_vector[9])
lines(DE.results.DEDD50$mean.discoveries$mean.MDSeq.nozi.rle, 
      DE.results.DEDD50$mean.fdr$mean.MDSeq.nozi.rle, 
      type='l', col=col_vector[10])
lines(DE.results.DEDD50$mean.discoveries$mean.expHM.untr.rle, 
      DE.results.DEDD50$mean.fdr$mean.expHM.untr.rle, 
      type='l', col=col_vector[11])
lines(DE.results.DEDD50$mean.discoveries$mean.expHM.log.rle, 
      DE.results.DEDD50$mean.fdr$mean.expHM.log.rle, 
      type='l', col=col_vector[12])
lines(DE.results.DEDD50$mean.discoveries$mean.lnHM.untr.rle, 
      DE.results.DEDD50$mean.fdr$mean.lnHM.untr.rle, 
      type='l', col=col_vector[13])
lines(DE.results.DEDD50$mean.discoveries$mean.lnHM.log.rle, 
      DE.results.DEDD50$mean.fdr$mean.lnHM.log.rle, 
      type='l', col=col_vector[14])
legend("topleft", bty='n', col=col_vector[1:14], lty=1, ncol=2, cex=0.9, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq2 no IF", "DESeq2 IF", 
                "voom", "DSS", "baySeq", "MDSeq ZI", "MDSeq no ZI", "expHM", "expHM log", 
                "lnHM", "lnHM log"))
abline(h=0.05, col='lightgrey')



###################################
#### Differential distribution ####
###################################

###########
### AUC ###
par(mfrow=c(2,3), mar=c(1,2,3,1), mgp=c(3,0.7,0))
boxplot(DEDD.results.DD2$auc, names=NA, col=col_vector[1:2], 
        pch=1.20, cex.axis=1.2, xaxt='n', ylim=c(0.45,1.05), 
        main=paste0("AUC for differential distribution, 2 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in dispersion only"), 
        cex.main=1.2)
lines(c(2.5,2.5), c(0.45,1))
legend("topleft", fill=col_vector[1:2], bty='n', cex=1.2, 
       legend=c("expHMM", "lnHMM"))
boxplot(DEDD.results.DD5$auc, names=NA, col=col_vector[1:2], 
        pch=1.20, cex.axis=1.2, xaxt='n', ylim=c(0.45,1.05), 
        main=paste0("AUC for differential distribution, 5 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in dispersion only"), 
        cex.main=1.2)
lines(c(2.5,2.5), c(0.45,1))
legend("topleft", fill=col_vector[1:2], bty='n', cex=1.2, 
       legend=c("expHMM", "lnHMM"))
boxplot(DEDD.results.DD10$auc, names=NA, col=col_vector[1:2], 
        pch=1.20, cex.axis=1.2, xaxt='n', ylim=c(0.45,1.05), 
        main=paste0("AUC for differential distribution, 10 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in dispersion only"), 
        cex.main=1.2)
lines(c(2.5,2.5), c(0.45,1))
legend("topleft", fill=col_vector[1:2], bty='n', cex=1.2, 
       legend=c("expHMM", "lnHMM"))
boxplot(DEDD.results.DD20$auc, names=NA, col=col_vector[1:2], 
        pch=1.20, cex.axis=1.2, xaxt='n', ylim=c(0.45,1.05), 
        main=paste0("AUC for differential distribution, 20 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in dispersion only"), 
        cex.main=1.2)
lines(c(2.5,2.5), c(0.45,1))
legend("topleft", fill=col_vector[1:2], bty='n', cex=1.2, 
       legend=c("expHMM", "lnHMM"))
boxplot(DEDD.results.DD50$auc, names=NA, col=col_vector[1:2], 
        pch=1.20, cex.axis=1.2, xaxt='n', ylim=c(0.45,1.05), 
        main=paste0("AUC for differential distribution, 50 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in dispersion only"), 
        cex.main=1.2)
lines(c(2.5,2.5), c(0.45,1))
legend("topleft", fill=col_vector[1:2], bty='n', cex=1.2, 
       legend=c("expHMM", "lnHMM"))
rbind(colMeans(DEDD.results.DD2$auc), colMeans(DEDD.results.DD5$auc), 
      colMeans(DEDD.results.DD10$auc), colMeans(DEDD.results.DD20$auc), 
      colMeans(DEDD.results.DD50$auc))

par(mfrow=c(2,3), mar=c(1,2,3,1), mgp=c(3,0.7,0))
boxplot(DEDD.results.DE2$auc, names=NA, col=col_vector[1:2], 
        pch=1.20, cex.axis=1.2, xaxt='n', ylim=c(0.45,1.05), 
        main=paste0("AUC for differential distribution, 2 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean only"), 
        cex.main=1.2)
lines(c(2.5,2.5), c(0.45,1))
legend("topleft", fill=col_vector[1:2], bty='n', cex=1.2, ncol=2, 
       legend=c("expHMM", "lnHMM"))
boxplot(DEDD.results.DE5$auc, names=NA, col=col_vector[1:2], 
        pch=1.20, cex.axis=1.2, xaxt='n', ylim=c(0.45,1.05), 
        main=paste0("AUC for differential distribution, 5 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean only"), 
        cex.main=1.2)
lines(c(2.5,2.5), c(0.45,1))
legend("topleft", fill=col_vector[1:2], bty='n', cex=1.2, ncol=2, 
       legend=c("expHMM", "lnHMM"))
boxplot(DEDD.results.DE10$auc, names=NA, col=col_vector[1:2], 
        pch=1.20, cex.axis=1.2, xaxt='n', ylim=c(0.45,1.05), 
        main=paste0("AUC for differential distribution, 10 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean only"), 
        cex.main=1.2)
lines(c(2.5,2.5), c(0.45,1))
legend("topleft", fill=col_vector[1:2], bty='n', cex=1.2, ncol=2, 
       legend=c("expHMM", "lnHMM"))
boxplot(DEDD.results.DE20$auc, names=NA, col=col_vector[1:2], 
        pch=1.20, cex.axis=1.2, xaxt='n', ylim=c(0.45,1.05), 
        main=paste0("AUC for differential distribution, 20 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean only"), 
        cex.main=1.2)
lines(c(2.5,2.5), c(0.45,1))
legend("topleft", fill=col_vector[1:2], bty='n', cex=1.2, ncol=2, 
       legend=c("expHMM", "lnHMM"))
boxplot(DEDD.results.DE50$auc, names=NA, col=col_vector[1:2], 
        pch=1.20, cex.axis=1.2, xaxt='n', ylim=c(0.45,1.05), 
        main=paste0("AUC for differential distribution, 50 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean only"), 
        cex.main=1.2)
lines(c(2.5,2.5), c(0.45,1))
legend("topleft", fill=col_vector[1:2], bty='n', cex=1.2, ncol=2, 
       legend=c("expHMM", "lnHMM"))
rbind(colMeans(DEDD.results.DE2$auc), colMeans(DEDD.results.DE5$auc), 
      colMeans(DEDD.results.DE10$auc), colMeans(DEDD.results.DE20$auc), 
      colMeans(DEDD.results.DE50$auc))

par(mfrow=c(2,3), mar=c(1,2,3,1), mgp=c(3,0.7,0))
boxplot(DEDD.results.DEDD2$auc, names=NA, col=col_vector[1:6], 
        pch=1.20, cex.axis=1.2, xaxt='n', ylim=c(0.45,1.05), 
        main=paste0("AUC for differential distribution, 2 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean and dispersion"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0.45,1))
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=2, 
       legend=c("expHMM", "lnHMM", "edgeR_MDSeq", "edgeR_lnHM", "lnHM_MDSeq", "lnHM_lnHM"))
boxplot(DEDD.results.DEDD5$auc, names=NA, col=col_vector[1:6], 
        pch=1.20, cex.axis=1.2, xaxt='n', ylim=c(0.45,1.05), 
        main=paste0("AUC for differential distribution, 5 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean and dispersion"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0.45,1))
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=2, 
       legend=c("expHMM", "lnHMM", "edgeR_MDSeq", "edgeR_lnHM", "lnHM_MDSeq", "lnHM_lnHM"))
boxplot(DEDD.results.DEDD10$auc, names=NA, col=col_vector[1:6], 
        pch=1.20, cex.axis=1.2, xaxt='n', ylim=c(0.45,1.05), 
        main=paste0("AUC for differential distribution, 10 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean and dispersion"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0.45,1))
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=2, 
       legend=c("expHMM", "lnHMM", "edgeR_MDSeq", "edgeR_lnHM", "lnHM_MDSeq", "lnHM_lnHM"))
boxplot(DEDD.results.DEDD20$auc, names=NA, col=col_vector[1:6], 
        pch=1.20, cex.axis=1.2, xaxt='n', ylim=c(0.45,1.05), 
        main=paste0("AUC for differential distribution, 20 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean and dispersion"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0.45,1))
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=2, 
       legend=c("expHMM", "lnHMM", "edgeR_MDSeq", "edgeR_lnHM", "lnHM_MDSeq", "lnHM_lnHM"))
boxplot(DEDD.results.DEDD50$auc, names=NA, col=col_vector[1:6], 
        pch=1.20, cex.axis=1.2, xaxt='n', ylim=c(0.45,1.05), 
        main=paste0("AUC for differential distribution, 50 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean and dispersion"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0.45,1))
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=2, 
       legend=c("expHMM", "lnHMM", "edgeR_MDSeq", "edgeR_lnHM", "lnHM_MDSeq", "lnHM_lnHM"))
rbind(colMeans(DEDD.results.DEDD2$auc), colMeans(DEDD.results.DEDD5$auc), 
      colMeans(DEDD.results.DEDD10$auc), colMeans(DEDD.results.DEDD20$auc), 
      colMeans(DEDD.results.DEDD50$auc))


###########
### pAUC ###
par(mfrow=c(2,3), mar=c(1,2,3,1), mgp=c(3,0.7,0))
boxplot(DEDD.results.DD2$pauc, names=NA, col=col_vector[1:2], 
        pch=1.20, cex.axis=1.2, xaxt='n', ylim=c(0,0.05), 
        main=paste0("Partial AUC for differential distribution, 2 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in dispersion only"), 
        cex.main=1.2)
lines(c(2.5,2.5), c(0,0.05))
legend("topleft", fill=col_vector[1:2], bty='n', cex=1.2, ncol=2, 
       legend=c("expHMM", "lnHMM"))
boxplot(DEDD.results.DD5$pauc, names=NA, col=col_vector[1:2], 
        pch=1.20, cex.axis=1.2, xaxt='n', ylim=c(0,0.05), 
        main=paste0("Partial AUC for differential distribution, 5 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in dispersion only"), 
        cex.main=1.2)
lines(c(2.5,2.5), c(0,0.05))
legend("topleft", fill=col_vector[1:2], bty='n', cex=1.2, ncol=2, 
       legend=c("expHMM", "lnHMM"))
boxplot(DEDD.results.DD10$pauc, names=NA, col=col_vector[1:2], 
        pch=1.20, cex.axis=1.2, xaxt='n', ylim=c(0,0.05), 
        main=paste0("Partial AUC for differential distribution, 10 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in dispersion only"), 
        cex.main=1.2)
lines(c(2.5,2.5), c(0,0.05))
legend("topleft", fill=col_vector[1:2], bty='n', cex=1.2, ncol=2, 
       legend=c("expHMM", "lnHMM"))
boxplot(DEDD.results.DD20$pauc, names=NA, col=col_vector[1:2], 
        pch=1.20, cex.axis=1.2, xaxt='n', ylim=c(0,0.05), 
        main=paste0("Partial AUC for differential distribution, 20 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in dispersion only"), 
        cex.main=1.2)
lines(c(2.5,2.5), c(0,0.05))
legend("topleft", fill=col_vector[1:2], bty='n', cex=1.2, ncol=2, 
       legend=c("expHMM", "lnHMM"))
boxplot(DEDD.results.DD50$pauc, names=NA, col=col_vector[1:2], 
        pch=1.20, cex.axis=1.2, xaxt='n', ylim=c(0,0.05), 
        main=paste0("Partial AUC for differential distribution, 50 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in dispersion only"), 
        cex.main=1.2)
lines(c(2.5,2.5), c(0,0.05))
legend("topleft", fill=col_vector[1:2], bty='n', cex=1.2, ncol=2, 
       legend=c("expHMM", "lnHMM"))
rbind(colMeans(DEDD.results.DD2$pauc), colMeans(DEDD.results.DD5$pauc), 
      colMeans(DEDD.results.DD10$pauc), colMeans(DEDD.results.DD20$pauc), 
      colMeans(DEDD.results.DD50$pauc))

par(mfrow=c(2,3), mar=c(1,2,3,1), mgp=c(3,0.7,0))
boxplot(DEDD.results.DE2$pauc, names=NA, col=col_vector[1:2], 
        pch=1.20, cex.axis=1.2, xaxt='n', ylim=c(0,0.05), 
        main=paste0("Partial AUC for differential distribution, 2 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean only"), 
        cex.main=1.2)
lines(c(2.5,2.5), c(0,0.05))
legend("topleft", fill=col_vector[1:2], bty='n', cex=1.2, ncol=2, 
       legend=c("expHMM", "lnHMM"))
boxplot(DEDD.results.DE5$pauc, names=NA, col=col_vector[1:2], 
        pch=1.20, cex.axis=1.2, xaxt='n', ylim=c(0,0.05), 
        main=paste0("Partial AUC for differential distribution, 5 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean only"), 
        cex.main=1.2)
lines(c(2.5,2.5), c(0,0.05))
legend("topleft", fill=col_vector[1:2], bty='n', cex=1.2, ncol=2, 
       legend=c("expHMM", "lnHMM"))
boxplot(DEDD.results.DE10$pauc, names=NA, col=col_vector[1:2], 
        pch=1.20, cex.axis=1.2, xaxt='n', ylim=c(0,0.05), 
        main=paste0("Partial AUC for differential distribution, 10 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean only"), 
        cex.main=1.2)
lines(c(2.5,2.5), c(0,0.05))
legend("topleft", fill=col_vector[1:2], bty='n', cex=1.2, ncol=2, 
       legend=c("expHMM", "lnHMM"))
boxplot(DEDD.results.DE20$pauc, names=NA, col=col_vector[1:2], 
        pch=1.20, cex.axis=1.2, xaxt='n', ylim=c(0,0.05), 
        main=paste0("Partial AUC for differential distribution, 20 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean only"), 
        cex.main=1.2)
lines(c(2.5,2.5), c(0,0.05))
legend("topleft", fill=col_vector[1:2], bty='n', cex=1.2, ncol=2, 
       legend=c("expHMM", "lnHMM"))
boxplot(DEDD.results.DE50$pauc, names=NA, col=col_vector[1:2], 
        pch=1.20, cex.axis=1.2, xaxt='n', ylim=c(0,0.05), 
        main=paste0("Partial AUC for differential distribution, 50 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean only"), 
        cex.main=1.2)
lines(c(2.5,2.5), c(0,0.05))
legend("topleft", fill=col_vector[1:2], bty='n', cex=1.2, ncol=2, 
       legend=c("expHMM", "lnHMM"))
rbind(colMeans(DEDD.results.DE2$pauc), colMeans(DEDD.results.DE5$pauc), 
      colMeans(DEDD.results.DE10$pauc), colMeans(DEDD.results.DE20$pauc), 
      colMeans(DEDD.results.DE50$pauc))

par(mfrow=c(2,3), mar=c(1,2,3,1), mgp=c(3,0.7,0))
boxplot(DEDD.results.DEDD2$pauc, names=NA, col=col_vector[1:6], 
        pch=1.20, cex.axis=1.2, xaxt='n', ylim=c(0,0.05), 
        main=paste0("Partial AUC for differential distribution, 2 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean and dispersion"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0,0.05))
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=2, 
       legend=c("expHMM", "lnHMM", "edgeR_MDSeq", "edgeR_lnHM", "lnHM_MDSeq", "lnHM_lnHM"))
boxplot(DEDD.results.DEDD5$pauc, names=NA, col=col_vector[1:6], 
        pch=1.20, cex.axis=1.2, xaxt='n', ylim=c(0,0.05), 
        main=paste0("Partial AUC for differential distribution, 5 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean and dispersion"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0,0.05))
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=2, 
       legend=c("expHMM", "lnHMM", "edgeR_MDSeq", "edgeR_lnHM", "lnHM_MDSeq", "lnHM_lnHM"))
boxplot(DEDD.results.DEDD10$pauc, names=NA, col=col_vector[1:6], 
        pch=1.20, cex.axis=1.2, xaxt='n', ylim=c(0,0.05), 
        main=paste0("Partial AUC for differential distribution, 10 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean and dispersion"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0,0.05))
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=2, 
       legend=c("expHMM", "lnHMM", "edgeR_MDSeq", "edgeR_lnHM", "lnHM_MDSeq", "lnHM_lnHM"))
boxplot(DEDD.results.DEDD20$pauc, names=NA, col=col_vector[1:6], 
        pch=1.20, cex.axis=1.2, xaxt='n', ylim=c(0,0.05), 
        main=paste0("Partial AUC for differential distribution, 20 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean and dispersion"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0,0.05))
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=2, 
       legend=c("expHMM", "lnHMM", "edgeR_MDSeq", "edgeR_lnHM", "lnHM_MDSeq", "lnHM_lnHM"))
boxplot(DEDD.results.DEDD50$pauc, names=NA, col=col_vector[1:6], 
        pch=1.20, cex.axis=1.2, xaxt='n', ylim=c(0,0.05), 
        main=paste0("Partial AUC for differential distribution, 50 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean and dispersion"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0,0.05))
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=2, 
       legend=c("expHMM", "lnHMM", "edgeR_MDSeq", "edgeR_lnHM", "lnHM_MDSeq", "lnHM_lnHM"))
rbind(colMeans(DEDD.results.DEDD2$pauc), colMeans(DEDD.results.DEDD5$pauc), 
      colMeans(DEDD.results.DEDD10$pauc), colMeans(DEDD.results.DEDD20$pauc), 
      colMeans(DEDD.results.DEDD50$pauc))


###########
### FDR ###
par(mfrow=c(2,3), mar=c(1,2,3,1), mgp=c(3,0.7,0))
boxplot(DEDD.results.DD2$fdr, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,1.2), 
        main=paste0("FDR for differential distribution, 2 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in dispersion only"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0,1)); abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=2, 
       legend=c("expHMM 0.5", "expHMM thr", "expHMM bfdr", 
                "lnHMM 0.5", "lnHMM thr", "lnHMM bfdr"))
boxplot(DEDD.results.DD5$fdr, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,1.2), 
        main=paste0("FDR for differential distribution, 5 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in dispersion only"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0,1)); abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=2, 
       legend=c("expHMM 0.5", "expHMM thr", "expHMM bfdr", 
                "lnHMM 0.5", "lnHMM thr", "lnHMM bfdr"))
boxplot(DEDD.results.DD10$fdr, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,1.2), 
        main=paste0("FDR for differential distribution, 10 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in dispersion only"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0,1)); abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=2, 
       legend=c("expHMM 0.5", "expHMM thr", "expHMM bfdr", 
                "lnHMM 0.5", "lnHMM thr", "lnHMM bfdr"))
boxplot(DEDD.results.DD20$fdr, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,1.2), 
        main=paste0("FDR for differential distribution, 20 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in dispersion only"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0,1)); abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=2, 
       legend=c("expHMM 0.5", "expHMM thr", "expHMM bfdr", 
                "lnHMM 0.5", "lnHMM thr", "lnHMM bfdr"))
boxplot(DEDD.results.DD50$fdr, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,1.2), 
        main=paste0("FDR for differential distribution, 50 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in dispersion only"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0,1)); abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=2, 
       legend=c("expHMM 0.5", "expHMM thr", "expHMM bfdr", 
                "lnHMM 0.5", "lnHMM thr", "lnHMM bfdr"))
rbind(colMeans(DEDD.results.DD2$fdr), colMeans(DEDD.results.DD5$fdr), 
      colMeans(DEDD.results.DD10$fdr), colMeans(DEDD.results.DD20$fdr), 
      colMeans(DEDD.results.DD50$fdr))

par(mfrow=c(2,3), mar=c(1,2,3,1), mgp=c(3,0.7,0))
boxplot(DEDD.results.DE2$fdr, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,1.2), 
        main=paste0("FDR for differential distribution, 2 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean only"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0,1)); abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=2, 
       legend=c("expHMM 0.5", "expHMM thr", "expHMM bfdr", 
                "lnHMM 0.5", "lnHMM thr", "lnHMM bfdr"))
boxplot(DEDD.results.DE5$fdr, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,1.2), 
        main=paste0("FDR for differential distribution, 5 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean only"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0,1)); abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=2, 
       legend=c("expHMM 0.5", "expHMM thr", "expHMM bfdr", 
                "lnHMM 0.5", "lnHMM thr", "lnHMM bfdr"))
boxplot(DEDD.results.DE10$fdr, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,1.2), 
        main=paste0("FDR for differential distribution, 10 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean only"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0,1)); abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=2, 
       legend=c("expHMM 0.5", "expHMM thr", "expHMM bfdr", 
                "lnHMM 0.5", "lnHMM thr", "lnHMM bfdr"))
boxplot(DEDD.results.DE20$fdr, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,1.2), 
        main=paste0("FDR for differential distribution, 20 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean only"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0,1)); abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=2, 
       legend=c("expHMM 0.5", "expHMM thr", "expHMM bfdr", 
                "lnHMM 0.5", "lnHMM thr", "lnHMM bfdr"))
boxplot(DEDD.results.DE50$fdr, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,1.2), 
        main=paste0("FDR for differential distribution, 50 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean only"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0,1)); abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=2, 
       legend=c("expHMM 0.5", "expHMM thr", "expHMM bfdr", 
                "lnHMM 0.5", "lnHMM thr", "lnHMM bfdr"))
rbind(colMeans(DEDD.results.DE2$fdr), colMeans(DEDD.results.DE5$fdr), 
      colMeans(DEDD.results.DE10$fdr), colMeans(DEDD.results.DE20$fdr), 
      colMeans(DEDD.results.DE50$fdr))

par(mfrow=c(2,3), mar=c(1,2,3,1), mgp=c(3,0.7,0))
boxplot(DEDD.results.DEDD2$fdr[, c(1:3,7:9,13:16,4:6,10:12,17:20)], 
        names=NA, col=col_vector[1:10], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,1.2), 
        main=paste0("FDR for differential distribution, 2 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean and dispersion"), 
        cex.main=1.2)
lines(c(10.5,10.5), c(0,1)); abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[1:10], bty='n', cex=1.2, ncol=2, 
       legend=c("expHMM 0.5", "expHMM thr", "expHMM bfdr", 
                "lnHMM 0.5", "lnHMM thr", "lnHMM bfdr", 
                "edgeR_MDSeq", "edgeR_lnHM", "lnHM_MDSeq", "lnHM_lnHM"))
boxplot(DEDD.results.DEDD5$fdr[, c(1:3,7:9,13:16,4:6,10:12,17:20)], 
        names=NA, col=col_vector[1:10], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,1.2), 
        main=paste0("FDR for differential distribution, 5 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean and dispersion"), 
        cex.main=1.2)
lines(c(10.5,10.5), c(0,1)); abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[1:10], bty='n', cex=1.2, ncol=2, 
       legend=c("expHMM 0.5", "expHMM thr", "expHMM bfdr", 
                "lnHMM 0.5", "lnHMM thr", "lnHMM bfdr", 
                "edgeR_MDSeq", "edgeR_lnHM", "lnHM_MDSeq", "lnHM_lnHM"))
boxplot(DEDD.results.DEDD10$fdr[, c(1:3,7:9,13:16,4:6,10:12,17:20)], 
        names=NA, col=col_vector[1:10], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,1.2), 
        main=paste0("FDR for differential distribution, 10 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean and dispersion"), 
        cex.main=1.2)
lines(c(10.5,10.5), c(0,1)); abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[1:10], bty='n', cex=1.2, ncol=2, 
       legend=c("expHMM 0.5", "expHMM thr", "expHMM bfdr", 
                "lnHMM 0.5", "lnHMM thr", "lnHMM bfdr", 
                "edgeR_MDSeq", "edgeR_lnHM", "lnHM_MDSeq", "lnHM_lnHM"))
boxplot(DEDD.results.DEDD20$fdr[, c(1:3,7:9,13:16,4:6,10:12,17:20)], 
        names=NA, col=col_vector[1:10], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,1.2), 
        main=paste0("FDR for differential distribution, 20 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean and dispersion"), 
        cex.main=1.2)
lines(c(10.5,10.5), c(0,1)); abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[1:10], bty='n', cex=1.2, ncol=2, 
       legend=c("expHMM 0.5", "expHMM thr", "expHMM bfdr", 
                "lnHMM 0.5", "lnHMM thr", "lnHMM bfdr", 
                "edgeR_MDSeq", "edgeR_lnHM", "lnHM_MDSeq", "lnHM_lnHM"))
boxplot(DEDD.results.DEDD50$fdr[, c(1:3,7:9,13:16,4:6,10:12,17:20)], 
        names=NA, col=col_vector[1:10], 
        pch=20, cex.axis=1.2, xaxt='n', ylim=c(0,1.2), 
        main=paste0("FDR for differential distribution, 50 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean and dispersion"), 
        cex.main=1.2)
lines(c(10.5,10.5), c(0,1)); abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[1:10], bty='n', cex=1.2, ncol=2, 
       legend=c("expHMM 0.5", "expHMM thr", "expHMM bfdr", 
                "lnHMM 0.5", "lnHMM thr", "lnHMM bfdr", 
                "edgeR_MDSeq", "edgeR_lnHM", "lnHM_MDSeq", "lnHM_lnHM"))
rbind(colMeans(DEDD.results.DEDD2$fdr), colMeans(DEDD.results.DEDD5$fdr), 
      colMeans(DEDD.results.DEDD10$fdr), colMeans(DEDD.results.DEDD20$fdr), 
      colMeans(DEDD.results.DEDD50$fdr))[, c(1:3,7:9,13:16,4:6,10:12,17:20)]


###########
### TPR ###
par(mfrow=c(2,3), mar=c(1,2,3,1), mgp=c(3,0.7,0))
boxplot(DEDD.results.DD2$tpr, names=NA, col=col_vector[1:6], 
        pch=1.20, cex.axis=1.2, xaxt='n', ylim=c(0,1.1), 
        main=paste0("TPR for differential distribution, 2 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in dispersion only"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0,0.9))
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=2, 
       legend=c("expHMM 0.5", "expHMM thr", "expHMM bfdr", 
                "lnHMM 0.5", "lnHMM thr", "lnHMM bfdr"))
boxplot(DEDD.results.DD5$tpr, names=NA, col=col_vector[1:6], 
        pch=1.20, cex.axis=1.2, xaxt='n', ylim=c(0,1.1), 
        main=paste0("TPR for differential distribution, 5 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in dispersion only"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0,0.9))
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=2, 
       legend=c("expHMM 0.5", "expHMM thr", "expHMM bfdr", 
                "lnHMM 0.5", "lnHMM thr", "lnHMM bfdr"))
boxplot(DEDD.results.DD10$tpr, names=NA, col=col_vector[1:6], 
        pch=1.20, cex.axis=1.2, xaxt='n', ylim=c(0,1.1), 
        main=paste0("TPR for differential distribution, 10 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in dispersion only"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0,0.9))
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=2, 
       legend=c("expHMM 0.5", "expHMM thr", "expHMM bfdr", 
                "lnHMM 0.5", "lnHMM thr", "lnHMM bfdr"))
boxplot(DEDD.results.DD20$tpr, names=NA, col=col_vector[1:6], 
        pch=1.20, cex.axis=1.2, xaxt='n', ylim=c(0,1.1), 
        main=paste0("TPR for differential distribution, 20 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in dispersion only"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0,0.9))
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=2, 
       legend=c("expHMM 0.5", "expHMM thr", "expHMM bfdr", 
                "lnHMM 0.5", "lnHMM thr", "lnHMM bfdr"))
boxplot(DEDD.results.DD50$tpr, names=NA, col=col_vector[1:6], 
        pch=1.20, cex.axis=1.2, xaxt='n', ylim=c(0,1.1), 
        main=paste0("TPR for differential distribution, 50 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in dispersion only"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0,0.9))
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=2, 
       legend=c("expHMM 0.5", "expHMM thr", "expHMM bfdr", 
                "lnHMM 0.5", "lnHMM thr", "lnHMM bfdr"))
rbind(colMeans(DEDD.results.DD2$tpr), colMeans(DEDD.results.DD5$tpr), 
      colMeans(DEDD.results.DD10$tpr), colMeans(DEDD.results.DD20$tpr), 
      colMeans(DEDD.results.DD50$tpr))

par(mfrow=c(2,3), mar=c(1,2,3,1), mgp=c(3,0.7,0))
boxplot(DEDD.results.DE2$tpr, names=NA, col=col_vector[1:6], 
        pch=1.20, cex.axis=1.2, xaxt='n', ylim=c(0,1.1), 
        main=paste0("TPR for differential distribution, 2 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean only"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0,0.9))
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=2, 
       legend=c("expHMM 0.5", "expHMM thr", "expHMM bfdr", 
                "lnHMM 0.5", "lnHMM thr", "lnHMM bfdr"))
boxplot(DEDD.results.DE5$tpr, names=NA, col=col_vector[1:6], 
        pch=1.20, cex.axis=1.2, xaxt='n', ylim=c(0,1.1), 
        main=paste0("TPR for differential distribution, 5 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean only"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0,0.9))
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=2, 
       legend=c("expHMM 0.5", "expHMM thr", "expHMM bfdr", 
                "lnHMM 0.5", "lnHMM thr", "lnHMM bfdr"))
boxplot(DEDD.results.DE10$tpr, names=NA, col=col_vector[1:6], 
        pch=1.20, cex.axis=1.2, xaxt='n', ylim=c(0,1.1), 
        main=paste0("TPR for differential distribution, 10 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean only"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0,0.9))
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=2, 
       legend=c("expHMM 0.5", "expHMM thr", "expHMM bfdr", 
                "lnHMM 0.5", "lnHMM thr", "lnHMM bfdr"))
boxplot(DEDD.results.DE20$tpr, names=NA, col=col_vector[1:6], 
        pch=1.20, cex.axis=1.2, xaxt='n', ylim=c(0,1.1), 
        main=paste0("TPR for differential distribution, 20 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean only"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0,0.9))
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=2, 
       legend=c("expHMM 0.5", "expHMM thr", "expHMM bfdr", 
                "lnHMM 0.5", "lnHMM thr", "lnHMM bfdr"))
boxplot(DEDD.results.DE50$tpr, names=NA, col=col_vector[1:6], 
        pch=1.20, cex.axis=1.2, xaxt='n', ylim=c(0,1.1), 
        main=paste0("TPR for differential distribution, 50 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean only"), 
        cex.main=1.2)
lines(c(6.5,6.5), c(0,0.9))
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=2, 
       legend=c("expHMM 0.5", "expHMM thr", "expHMM bfdr", 
                "lnHMM 0.5", "lnHMM thr", "lnHMM bfdr"))
rbind(colMeans(DEDD.results.DE2$tpr), colMeans(DEDD.results.DE5$tpr), 
      colMeans(DEDD.results.DE10$tpr), colMeans(DEDD.results.DE20$tpr), 
      colMeans(DEDD.results.DE50$tpr))

par(mfrow=c(2,3), mar=c(1,2,3,1), mgp=c(3,0.7,0))
boxplot(DEDD.results.DEDD2$tpr[, c(1:3,7:9,13:16,4:6,10:12,17:20)], 
        names=NA, col=col_vector[1:10], 
        pch=1.20, cex.axis=1.2, xaxt='n', ylim=c(0,1.1), 
        main=paste0("TPR for differential distribution, 2 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean and dispersion"), 
        cex.main=1.2)
lines(c(10.5,10.5), c(0,0.9))
legend("topleft", fill=col_vector[1:10], bty='n', cex=1.2, ncol=2, 
       legend=c("expHMM 0.5", "expHMM thr", "expHMM bfdr", 
                "lnHMM 0.5", "lnHMM thr", "lnHMM bfdr", 
                "edgeR_MDSeq", "edgeR_lnHM", "lnHM_MDSeq", "lnHM_lnHM"))
boxplot(DEDD.results.DEDD5$tpr[, c(1:3,7:9,13:16,4:6,10:12,17:20)], 
        names=NA, col=col_vector[1:10], 
        pch=1.20, cex.axis=1.2, xaxt='n', ylim=c(0,1.1), 
        main=paste0("TPR for differential distribution, 5 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean and dispersion"), 
        cex.main=1.2)
lines(c(10.5,10.5), c(0,0.9))
legend("topleft", fill=col_vector[1:10], bty='n', cex=1.2, ncol=2, 
       legend=c("expHMM 0.5", "expHMM thr", "expHMM bfdr", 
                "lnHMM 0.5", "lnHMM thr", "lnHMM bfdr", 
                "edgeR_MDSeq", "edgeR_lnHM", "lnHM_MDSeq", "lnHM_lnHM"))
boxplot(DEDD.results.DEDD10$tpr[, c(1:3,7:9,13:16,4:6,10:12,17:20)], 
        names=NA, col=col_vector[1:10], 
        pch=1.20, cex.axis=1.2, xaxt='n', ylim=c(0,1.1), 
        main=paste0("TPR for differential distribution, 10 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean and dispersion"), 
        cex.main=1.2)
lines(c(10.5,10.5), c(0,0.9))
legend("topleft", fill=col_vector[1:10], bty='n', cex=1.2, ncol=2, 
       legend=c("expHMM 0.5", "expHMM thr", "expHMM bfdr", 
                "lnHMM 0.5", "lnHMM thr", "lnHMM bfdr", 
                "edgeR_MDSeq", "edgeR_lnHM", "lnHM_MDSeq", "lnHM_lnHM"))
boxplot(DEDD.results.DEDD20$tpr[, c(1:3,7:9,13:16,4:6,10:12,17:20)], 
        names=NA, col=col_vector[1:10], 
        pch=1.20, cex.axis=1.2, xaxt='n', ylim=c(0,1.1), 
        main=paste0("TPR for differential distribution, 20 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean and dispersion"), 
        cex.main=1.2)
lines(c(10.5,10.5), c(0,0.9))
legend("topleft", fill=col_vector[1:10], bty='n', cex=1.2, ncol=2, 
       legend=c("expHMM 0.5", "expHMM thr", "expHMM bfdr", 
                "lnHMM 0.5", "lnHMM thr", "lnHMM bfdr", 
                "edgeR_MDSeq", "edgeR_lnHM", "lnHM_MDSeq", "lnHM_lnHM"))
boxplot(DEDD.results.DEDD50$tpr[, c(1:3,7:9,13:16,4:6,10:12,17:20)], 
        names=NA, col=col_vector[1:10], 
        pch=1.20, cex.axis=1.2, xaxt='n', ylim=c(0,1.1), 
        main=paste0("TPR for differential distribution, 50 samples per group, ", 
                    "TMM (left)\n or RLE (right) normalisation, ", 
                    "differences in mean and dispersion"), 
        cex.main=1.2)
lines(c(10.5,10.5), c(0,0.9))
legend("topleft", fill=col_vector[1:10], bty='n', cex=1.2, ncol=2, 
       legend=c("expHMM 0.5", "expHMM thr", "expHMM bfdr", 
                "lnHMM 0.5", "lnHMM thr", "lnHMM bfdr", 
                "edgeR_MDSeq", "edgeR_lnHM", "lnHM_MDSeq", "lnHM_lnHM"))
rbind(colMeans(DEDD.results.DEDD2$tpr), colMeans(DEDD.results.DEDD5$tpr), 
      colMeans(DEDD.results.DEDD10$tpr), colMeans(DEDD.results.DEDD20$tpr), 
      colMeans(DEDD.results.DEDD50$tpr))[, c(1:3,7:9,13:16,4:6,10:12,17:20)]


#############################
### False discovery plots ###
par(mfrow=c(6,5), mar=c(2,2,3,1), mgp=c(3,0.7,0))

plot(DEDD.results.DD2$mean.discoveries$lnHMM.tmm, 
     DEDD.results.DD2$mean.fdr$lnHMM.tmm, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("2 samples per group\n", "Differences in dispersion only, ", 
     "TMM"))
lines(DEDD.results.DD2$mean.discoveries$expHMM.tmm, 
      DEDD.results.DD2$mean.fdr$expHMM.tmm, 
      type='l', col=col_vector[2])
plot(DEDD.results.DD5$mean.discoveries$lnHMM.tmm, 
     DEDD.results.DD5$mean.fdr$lnHMM.tmm, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("5 samples per group\n", "Differences in dispersion only, ", 
     "TMM"))
lines(DEDD.results.DD5$mean.discoveries$expHMM.tmm, 
      DEDD.results.DD5$mean.fdr$expHMM.tmm, 
      type='l', col=col_vector[2])
plot(DEDD.results.DD10$mean.discoveries$lnHMM.tmm, 
     DEDD.results.DD10$mean.fdr$lnHMM.tmm, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("10 samples per group\n", "Differences in dispersion only, ", 
     "TMM"))
lines(DEDD.results.DD10$mean.discoveries$expHMM.tmm, 
      DEDD.results.DD10$mean.fdr$expHMM.tmm, 
      type='l', col=col_vector[2])
plot(DEDD.results.DD20$mean.discoveries$lnHMM.tmm, 
     DEDD.results.DD20$mean.fdr$lnHMM.tmm, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("20 samples per group\n", "Differences in dispersion only, ", 
     "TMM"))
lines(DEDD.results.DD20$mean.discoveries$expHMM.tmm, 
      DEDD.results.DD20$mean.fdr$expHMM.tmm, 
      type='l', col=col_vector[2])
abline(h=0.05, col='lightgrey')
plot(DEDD.results.DD50$mean.discoveries$lnHMM.tmm, 
     DEDD.results.DD50$mean.fdr$lnHMM.tmm, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("2 samples per group\n", "Differences in dispersion only, ", 
     "TMM"))
lines(DEDD.results.DD50$mean.discoveries$expHMM.tmm, 
      DEDD.results.DD50$mean.fdr$expHMM.tmm, 
      type='l', col=col_vector[2])
abline(h=0.05, col='lightgrey')
legend("topright", bty='n', col=col_vector[1:2], lty=1, ncol=2, 
       legend=c("expHMM", "lnHMM"))

plot(DEDD.results.DD2$mean.discoveries$lnHMM.rle, 
     DEDD.results.DD2$mean.fdr$lnHMM.rle, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("2 samples per group\n", "Differences in dispersion only, ", 
     "RLE"))
lines(DEDD.results.DD2$mean.discoveries$expHMM.rle, 
      DEDD.results.DD2$mean.fdr$expHMM.rle, 
      type='l', col=col_vector[2])
plot(DEDD.results.DD5$mean.discoveries$lnHMM.rle, 
     DEDD.results.DD5$mean.fdr$lnHMM.rle, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("5 samples per group\n", "Differences in dispersion only, ", 
     "RLE"))
lines(DEDD.results.DD5$mean.discoveries$expHMM.rle, 
      DEDD.results.DD5$mean.fdr$expHMM.rle, 
      type='l', col=col_vector[2])
plot(DEDD.results.DD10$mean.discoveries$lnHMM.rle, 
     DEDD.results.DD10$mean.fdr$lnHMM.rle, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("10 samples per group\n", "Differences in dispersion only, ", 
     "RLE"))
lines(DEDD.results.DD10$mean.discoveries$expHMM.rle, 
      DEDD.results.DD10$mean.fdr$expHMM.rle, 
      type='l', col=col_vector[2])
plot(DEDD.results.DD20$mean.discoveries$lnHMM.rle, 
     DEDD.results.DD20$mean.fdr$lnHMM.rle, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("20 samples per group\n", "Differences in dispersion only, ", 
     "RLE"))
lines(DEDD.results.DD20$mean.discoveries$expHMM.rle, 
      DEDD.results.DD20$mean.fdr$expHMM.rle, 
      type='l', col=col_vector[2])
abline(h=0.05, col='lightgrey')
plot(DEDD.results.DD50$mean.discoveries$lnHMM.rle, 
     DEDD.results.DD50$mean.fdr$lnHMM.rle, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("2 samples per group\n", "Differences in dispersion only, ", 
     "RLE"))
lines(DEDD.results.DD50$mean.discoveries$expHMM.rle, 
      DEDD.results.DD50$mean.fdr$expHMM.rle, 
      type='l', col=col_vector[2])
abline(h=0.05, col='lightgrey')
legend("topright", bty='n', col=col_vector[1:2], lty=1, ncol=2, 
       legend=c("expHMM", "lnHMM"))

plot(DEDD.results.DE2$mean.discoveries$lnHMM.tmm, 
     DEDD.results.DE2$mean.fdr$lnHMM.tmm, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("2 samples per group\n", "Differences in mean only, ", 
     "TMM"))
lines(DEDD.results.DE2$mean.discoveries$expHMM.tmm, 
      DEDD.results.DE2$mean.fdr$expHMM.tmm, 
      type='l', col=col_vector[2])
plot(DEDD.results.DE5$mean.discoveries$lnHMM.tmm, 
     DEDD.results.DE5$mean.fdr$lnHMM.tmm, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("5 samples per group\n", "Differences in mean only, ", 
     "TMM"))
lines(DEDD.results.DE5$mean.discoveries$expHMM.tmm, 
      DEDD.results.DE5$mean.fdr$expHMM.tmm, 
      type='l', col=col_vector[2])
abline(h=0.05, col='lightgrey')
plot(DEDD.results.DE10$mean.discoveries$lnHMM.tmm, 
     DEDD.results.DE10$mean.fdr$lnHMM.tmm, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("10 samples per group\n", "Differences in mean only, ", 
     "TMM"))
lines(DEDD.results.DE10$mean.discoveries$expHMM.tmm, 
      DEDD.results.DE10$mean.fdr$expHMM.tmm, 
      type='l', col=col_vector[2])
abline(h=0.05, col='lightgrey')
plot(DEDD.results.DE20$mean.discoveries$lnHMM.tmm, 
     DEDD.results.DE20$mean.fdr$lnHMM.tmm, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("20 samples per group\n", "Differences in mean only, ", 
     "TMM"))
lines(DEDD.results.DE20$mean.discoveries$expHMM.tmm, 
      DEDD.results.DE20$mean.fdr$expHMM.tmm, 
      type='l', col=col_vector[2])
abline(h=0.05, col='lightgrey')
plot(DEDD.results.DE50$mean.discoveries$lnHMM.tmm, 
     DEDD.results.DE50$mean.fdr$lnHMM.tmm, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("2 samples per group\n", "Differences in mean only, ", 
     "TMM"))
lines(DEDD.results.DE50$mean.discoveries$expHMM.tmm, 
      DEDD.results.DE50$mean.fdr$expHMM.tmm, 
      type='l', col=col_vector[2])
abline(h=0.05, col='lightgrey')
legend("topleft", bty='n', col=col_vector[1:2], lty=1, ncol=2, 
       legend=c("expHMM", "lnHMM"))

plot(DEDD.results.DE2$mean.discoveries$lnHMM.rle, 
     DEDD.results.DE2$mean.fdr$lnHMM.rle, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("2 samples per group\n", "Differences in mean only, ", 
     "RLE"))
lines(DEDD.results.DE2$mean.discoveries$expHMM.rle, 
      DEDD.results.DE2$mean.fdr$expHMM.rle, 
      type='l', col=col_vector[2])
plot(DEDD.results.DE5$mean.discoveries$lnHMM.rle, 
     DEDD.results.DE5$mean.fdr$lnHMM.rle, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("5 samples per group\n", "Differences in mean only, ", 
     "RLE"))
lines(DEDD.results.DE5$mean.discoveries$expHMM.rle, 
      DEDD.results.DE5$mean.fdr$expHMM.rle, 
      type='l', col=col_vector[2])
abline(h=0.05, col='lightgrey')
plot(DEDD.results.DE10$mean.discoveries$lnHMM.rle, 
     DEDD.results.DE10$mean.fdr$lnHMM.rle, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("10 samples per group\n", "Differences in mean only, ", 
     "RLE"))
lines(DEDD.results.DE10$mean.discoveries$expHMM.rle, 
      DEDD.results.DE10$mean.fdr$expHMM.rle, 
      type='l', col=col_vector[2])
abline(h=0.05, col='lightgrey')
plot(DEDD.results.DE20$mean.discoveries$lnHMM.rle, 
     DEDD.results.DE20$mean.fdr$lnHMM.rle, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("20 samples per group\n", "Differences in mean only, ", 
     "RLE"))
lines(DEDD.results.DE20$mean.discoveries$expHMM.rle, 
      DEDD.results.DE20$mean.fdr$expHMM.rle, 
      type='l', col=col_vector[2])
abline(h=0.05, col='lightgrey')
plot(DEDD.results.DE50$mean.discoveries$lnHMM.rle, 
     DEDD.results.DE50$mean.fdr$lnHMM.rle, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("2 samples per group\n", "Differences in mean only, ", 
     "RLE"))
lines(DEDD.results.DE50$mean.discoveries$expHMM.rle, 
      DEDD.results.DE50$mean.fdr$expHMM.rle, 
      type='l', col=col_vector[2])
abline(h=0.05, col='lightgrey')
legend("topleft", bty='n', col=col_vector[1:2], lty=1, ncol=2, 
       legend=c("expHMM", "lnHMM"))

plot(DEDD.results.DEDD2$mean.discoveries$expHMM.tmm, 
     DEDD.results.DEDD2$mean.fdr$expHMM.tmm, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("2 samples per group\n", "Differences in mean & disp, ", 
     "TMM"))
lines(DEDD.results.DEDD2$mean.discoveries$lnHMM.tmm, 
      DEDD.results.DEDD2$mean.fdr$lnHMM.tmm, 
      type='l', col=col_vector[2])
lines(DEDD.results.DEDD2$mean.discoveries$edgeR_MD.tmm, 
      DEDD.results.DEDD2$mean.fdr$edgeR_MD.tmm, 
      type='l', col=col_vector[3])
lines(DEDD.results.DEDD2$mean.discoveries$edgeR_HM.tmm, 
      DEDD.results.DEDD2$mean.fdr$edgeR_HM.tmm, 
      type='l', col=col_vector[4])
lines(DEDD.results.DEDD2$mean.discoveries$HM_MD.tmm, 
      DEDD.results.DEDD2$mean.fdr$HM_MD.tmm, 
      type='l', col=col_vector[5])
lines(DEDD.results.DEDD2$mean.discoveries$HM.HM.tmm, 
      DEDD.results.DEDD2$mean.fdr$HM.HM.tmm, 
      type='l', col=col_vector[6])
plot(DEDD.results.DEDD5$mean.discoveries$expHMM.tmm, 
     DEDD.results.DEDD5$mean.fdr$expHMM.tmm, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("5 samples per group\n", "Differences in mean & disp, ", 
     "TMM"))
lines(DEDD.results.DEDD5$mean.discoveries$lnHMM.tmm, 
      DEDD.results.DEDD5$mean.fdr$lnHMM.tmm, 
      type='l', col=col_vector[2])
lines(DEDD.results.DEDD5$mean.discoveries$edgeR_MD.tmm, 
      DEDD.results.DEDD5$mean.fdr$edgeR_MD.tmm, 
      type='l', col=col_vector[3])
lines(DEDD.results.DEDD5$mean.discoveries$edgeR_HM.tmm, 
      DEDD.results.DEDD5$mean.fdr$edgeR_HM.tmm, 
      type='l', col=col_vector[4])
lines(DEDD.results.DEDD5$mean.discoveries$HM_MD.tmm, 
      DEDD.results.DEDD5$mean.fdr$HM_MD.tmm, 
      type='l', col=col_vector[5])
lines(DEDD.results.DEDD5$mean.discoveries$HM.HM.tmm, 
      DEDD.results.DEDD5$mean.fdr$HM.HM.tmm, 
      type='l', col=col_vector[6])
abline(h=0.05, col='lightgrey')
plot(DEDD.results.DEDD10$mean.discoveries$expHMM.tmm, 
     DEDD.results.DEDD10$mean.fdr$expHMM.tmm, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("10 samples per group\n", "Differences in mean & disp, ", 
     "TMM"))
lines(DEDD.results.DEDD10$mean.discoveries$lnHMM.tmm, 
      DEDD.results.DEDD10$mean.fdr$lnHMM.tmm, 
      type='l', col=col_vector[2])
lines(DEDD.results.DEDD10$mean.discoveries$edgeR_MD.tmm, 
      DEDD.results.DEDD10$mean.fdr$edgeR_MD.tmm, 
      type='l', col=col_vector[3])
lines(DEDD.results.DEDD10$mean.discoveries$edgeR_HM.tmm, 
      DEDD.results.DEDD10$mean.fdr$edgeR_HM.tmm, 
      type='l', col=col_vector[4])
lines(DEDD.results.DEDD10$mean.discoveries$HM_MD.tmm, 
      DEDD.results.DEDD10$mean.fdr$HM_MD.tmm, 
      type='l', col=col_vector[5])
lines(DEDD.results.DEDD10$mean.discoveries$HM.HM.tmm, 
      DEDD.results.DEDD10$mean.fdr$HM.HM.tmm, 
      type='l', col=col_vector[6])
abline(h=0.05, col='lightgrey')
plot(DEDD.results.DEDD20$mean.discoveries$expHMM.tmm, 
     DEDD.results.DEDD20$mean.fdr$expHMM.tmm, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("20 samples per group\n", "Differences in mean & disp, ", 
     "TMM"))
lines(DEDD.results.DEDD20$mean.discoveries$lnHMM.tmm, 
      DEDD.results.DEDD20$mean.fdr$lnHMM.tmm, 
      type='l', col=col_vector[2])
lines(DEDD.results.DEDD20$mean.discoveries$edgeR_MD.tmm, 
      DEDD.results.DEDD20$mean.fdr$edgeR_MD.tmm, 
      type='l', col=col_vector[3])
lines(DEDD.results.DEDD20$mean.discoveries$edgeR_HM.tmm, 
      DEDD.results.DEDD20$mean.fdr$edgeR_HM.tmm, 
      type='l', col=col_vector[4])
lines(DEDD.results.DEDD20$mean.discoveries$HM_MD.tmm, 
      DEDD.results.DEDD20$mean.fdr$HM_MD.tmm, 
      type='l', col=col_vector[5])
lines(DEDD.results.DEDD20$mean.discoveries$HM.HM.tmm, 
      DEDD.results.DEDD20$mean.fdr$HM.HM.tmm, 
      type='l', col=col_vector[6])
abline(h=0.05, col='lightgrey')
plot(DEDD.results.DEDD50$mean.discoveries$expHMM.tmm, 
     DEDD.results.DEDD50$mean.fdr$expHMM.tmm, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("2 samples per group\n", "Differences in mean & disp, ", 
     "TMM"))
lines(DEDD.results.DEDD50$mean.discoveries$lnHMM.tmm, 
      DEDD.results.DEDD50$mean.fdr$lnHMM.tmm, 
      type='l', col=col_vector[2])
lines(DEDD.results.DEDD50$mean.discoveries$edgeR_MD.tmm, 
      DEDD.results.DEDD50$mean.fdr$edgeR_MD.tmm, 
      type='l', col=col_vector[3])
lines(DEDD.results.DEDD50$mean.discoveries$edgeR_HM.tmm, 
      DEDD.results.DEDD50$mean.fdr$edgeR_HM.tmm, 
      type='l', col=col_vector[4])
lines(DEDD.results.DEDD50$mean.discoveries$HM_MD.tmm, 
      DEDD.results.DEDD50$mean.fdr$HM_MD.tmm, 
      type='l', col=col_vector[5])
lines(DEDD.results.DEDD50$mean.discoveries$HM.HM.tmm, 
      DEDD.results.DEDD50$mean.fdr$HM.HM.tmm, 
      type='l', col=col_vector[6])
legend("topright", bty='n', col=col_vector[1:6], lty=1, ncol=2, 
       legend=c("expHMM", "lnHMM", "edgeR_MDSeq", "edgeR_lnHM", "lnHM_MDSeq", "lnHM_lnHM"))
abline(h=0.05, col='lightgrey')

plot(DEDD.results.DEDD2$mean.discoveries$expHMM.rle, 
     DEDD.results.DEDD2$mean.fdr$expHMM.rle, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("2 samples per group\n", "Differences in mean & disp, ", 
     "RLE"))
lines(DEDD.results.DEDD2$mean.discoveries$lnHMM.rle, 
      DEDD.results.DEDD2$mean.fdr$lnHMM.rle, 
      type='l', col=col_vector[2])
lines(DEDD.results.DEDD2$mean.discoveries$edgeR_MD.rle, 
      DEDD.results.DEDD2$mean.fdr$edgeR_MD.rle, 
      type='l', col=col_vector[3])
lines(DEDD.results.DEDD2$mean.discoveries$edgeR_HM.rle, 
      DEDD.results.DEDD2$mean.fdr$edgeR_HM.rle, 
      type='l', col=col_vector[4])
lines(DEDD.results.DEDD2$mean.discoveries$HM_MD.rle, 
      DEDD.results.DEDD2$mean.fdr$HM_MD.rle, 
      type='l', col=col_vector[5])
lines(DEDD.results.DEDD2$mean.discoveries$HM.HM.rle, 
      DEDD.results.DEDD2$mean.fdr$HM.HM.rle, 
      type='l', col=col_vector[6])
plot(DEDD.results.DEDD5$mean.discoveries$expHMM.rle, 
     DEDD.results.DEDD5$mean.fdr$expHMM.rle, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("5 samples per group\n", "Differences in mean & disp, ", 
     "RLE"))
lines(DEDD.results.DEDD5$mean.discoveries$lnHMM.rle, 
      DEDD.results.DEDD5$mean.fdr$lnHMM.rle, 
      type='l', col=col_vector[2])
lines(DEDD.results.DEDD5$mean.discoveries$edgeR_MD.rle, 
      DEDD.results.DEDD5$mean.fdr$edgeR_MD.rle, 
      type='l', col=col_vector[3])
lines(DEDD.results.DEDD5$mean.discoveries$edgeR_HM.rle, 
      DEDD.results.DEDD5$mean.fdr$edgeR_HM.rle, 
      type='l', col=col_vector[4])
lines(DEDD.results.DEDD5$mean.discoveries$HM_MD.rle, 
      DEDD.results.DEDD5$mean.fdr$HM_MD.rle, 
      type='l', col=col_vector[5])
lines(DEDD.results.DEDD5$mean.discoveries$HM.HM.rle, 
      DEDD.results.DEDD5$mean.fdr$HM.HM.rle, 
      type='l', col=col_vector[6])
abline(h=0.05, col='lightgrey')
plot(DEDD.results.DEDD10$mean.discoveries$expHMM.rle, 
     DEDD.results.DEDD10$mean.fdr$expHMM.rle, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("10 samples per group\n", "Differences in mean & disp, ", 
     "RLE"))
lines(DEDD.results.DEDD10$mean.discoveries$lnHMM.rle, 
      DEDD.results.DEDD10$mean.fdr$lnHMM.rle, 
      type='l', col=col_vector[2])
lines(DEDD.results.DEDD10$mean.discoveries$edgeR_MD.rle, 
      DEDD.results.DEDD10$mean.fdr$edgeR_MD.rle, 
      type='l', col=col_vector[3])
lines(DEDD.results.DEDD10$mean.discoveries$edgeR_HM.rle, 
      DEDD.results.DEDD10$mean.fdr$edgeR_HM.rle, 
      type='l', col=col_vector[4])
lines(DEDD.results.DEDD10$mean.discoveries$HM_MD.rle, 
      DEDD.results.DEDD10$mean.fdr$HM_MD.rle, 
      type='l', col=col_vector[5])
lines(DEDD.results.DEDD10$mean.discoveries$HM.HM.rle, 
      DEDD.results.DEDD10$mean.fdr$HM.HM.rle, 
      type='l', col=col_vector[6])
abline(h=0.05, col='lightgrey')
plot(DEDD.results.DEDD20$mean.discoveries$expHMM.rle, 
     DEDD.results.DEDD20$mean.fdr$expHMM.rle, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("20 samples per group\n", "Differences in mean & disp, ", 
     "RLE"))
lines(DEDD.results.DEDD20$mean.discoveries$lnHMM.rle, 
      DEDD.results.DEDD20$mean.fdr$lnHMM.rle, 
      type='l', col=col_vector[2])
lines(DEDD.results.DEDD20$mean.discoveries$edgeR_MD.rle, 
      DEDD.results.DEDD20$mean.fdr$edgeR_MD.rle, 
      type='l', col=col_vector[3])
lines(DEDD.results.DEDD20$mean.discoveries$edgeR_HM.rle, 
      DEDD.results.DEDD20$mean.fdr$edgeR_HM.rle, 
      type='l', col=col_vector[4])
lines(DEDD.results.DEDD20$mean.discoveries$HM_MD.rle, 
      DEDD.results.DEDD20$mean.fdr$HM_MD.rle, 
      type='l', col=col_vector[5])
lines(DEDD.results.DEDD20$mean.discoveries$HM.HM.rle, 
      DEDD.results.DEDD20$mean.fdr$HM.HM.rle, 
      type='l', col=col_vector[6])
abline(h=0.05, col='lightgrey')
plot(DEDD.results.DEDD50$mean.discoveries$expHMM.rle, 
     DEDD.results.DEDD50$mean.fdr$expHMM.rle, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("2 samples per group\n", "Differences in mean & disp, ", 
     "RLE"))
lines(DEDD.results.DEDD50$mean.discoveries$lnHMM.rle, 
      DEDD.results.DEDD50$mean.fdr$lnHMM.rle, 
      type='l', col=col_vector[2])
lines(DEDD.results.DEDD50$mean.discoveries$edgeR_MD.rle, 
      DEDD.results.DEDD50$mean.fdr$edgeR_MD.rle, 
      type='l', col=col_vector[3])
lines(DEDD.results.DEDD50$mean.discoveries$edgeR_HM.rle, 
      DEDD.results.DEDD50$mean.fdr$edgeR_HM.rle, 
      type='l', col=col_vector[4])
lines(DEDD.results.DEDD50$mean.discoveries$HM_MD.rle, 
      DEDD.results.DEDD50$mean.fdr$HM_MD.rle, 
      type='l', col=col_vector[5])
lines(DEDD.results.DEDD50$mean.discoveries$HM.HM.rle, 
      DEDD.results.DEDD50$mean.fdr$HM.HM.rle, 
      type='l', col=col_vector[6])
legend("topright", bty='n', col=col_vector[1:6], lty=1, ncol=2, 
       legend=c("expHMM", "lnHMM", "edgeR_MDSeq", "edgeR_lnHM", "lnHM_MDSeq", "lnHM_lnHM"))
abline(h=0.05, col='lightgrey')

