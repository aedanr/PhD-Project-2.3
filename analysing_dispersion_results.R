library(here)
# DEDD2.TMM <- readRDS(here('Results/Dispersion estimation results Aug-Oct 2019','mse.disp.DEDD2.TMM.rds'))
# DEDD5.TMM <- readRDS(here('Results/Dispersion estimation results Aug-Oct 2019','mse.disp.DEDD5.TMM.rds'))
# DEDD10.TMM <- readRDS(here('Results/Dispersion estimation results Aug-Oct 2019','mse.disp.DEDD10.TMM.rds'))
DE2.TMM <- readRDS(here('Results/Dispersion estimation results Aug-Oct 2019','mse.disp.DE2.TMM.rds'))
DE5.TMM <- readRDS(here('Results/Dispersion estimation results Aug-Oct 2019','mse.disp.DE5.TMM.rds'))
DE10.TMM <- readRDS(here('Results/Dispersion estimation results Aug-Oct 2019','mse.disp.DE10.TMM.rds'))
DE20.TMM <- readRDS(here('Results/Dispersion estimation results Aug-Oct 2019','mse.disp.DE20.TMM.rds'))
DE2.DESeqnorm <- readRDS(here('Results/Dispersion estimation results Aug-Oct 2019','mse.disp.DE2.DESeqnorm.rds'))
DE5.DESeqnorm <- readRDS(here('Results/Dispersion estimation results Aug-Oct 2019','mse.disp.DE5.DESeqnorm.rds'))
DE10.DESeqnorm <- readRDS(here('Results/Dispersion estimation results Aug-Oct 2019','mse.disp.DE10.DESeqnorm.rds'))
DE20.DESeqnorm <- readRDS(here('Results/Dispersion estimation results Aug-Oct 2019','mse.disp.DE20.DESeqnorm.rds'))
nodiff2.TMM <- readRDS(here('Results/Dispersion estimation results Aug-Oct 2019','mse.disp.nodiff2.TMM.rds'))
nodiff5.TMM <- readRDS(here('Results/Dispersion estimation results Aug-Oct 2019','mse.disp.nodiff5.TMM.rds'))
nodiff10.TMM <- readRDS(here('Results/Dispersion estimation results Aug-Oct 2019','mse.disp.nodiff10.TMM.rds'))
nodiff20.TMM <- readRDS(here('Results/Dispersion estimation results Aug-Oct 2019','mse.disp.nodiff20.TMM.rds'))
nodiff2.DESeqnorm <- readRDS(here('Results/Dispersion estimation results Aug-Oct 2019','mse.disp.nodiff2.DESeqnorm.rds'))
nodiff5.DESeqnorm <- readRDS(here('Results/Dispersion estimation results Aug-Oct 2019','mse.disp.nodiff5.DESeqnorm.rds'))
nodiff10.DESeqnorm <- readRDS(here('Results/Dispersion estimation results Aug-Oct 2019','mse.disp.nodiff10.DESeqnorm.rds'))
nodiff20.DESeqnorm <- readRDS(here('Results/Dispersion estimation results Aug-Oct 2019','mse.disp.nodiff20.DESeqnorm.rds'))

# names(DEDD2$raw.MSEs)
# names(DEDD2$raw.MSEs$noDD.noDE)
# par(mfrow=c(2,1), mar=c(3,2,1,1))
# boxplot(DEDD2.TMM$raw.MSEs$noDD.noDE, cex.axis=0.8)
# boxplot(DEDD2.TMM$raw.MSEs$noDD.DE, cex.axis=0.8)
# colMeans(DEDD2.TMM$raw.MSEs$noDD.noDE)
# colMeans(DEDD2.TMM$raw.MSEs$noDD.DE)
# par(mfrow=c(2,1), mar=c(3,2,1,1))
# boxplot(DEDD5.TMM$raw.MSEs$noDD.noDE, cex.axis=0.8)
# boxplot(DEDD5.TMM$raw.MSEs$noDD.DE, cex.axis=0.8)
# colMeans(DEDD5.TMM$raw.MSEs$noDD.noDE)
# colMeans(DEDD5.TMM$raw.MSEs$noDD.DE)
# par(mfrow=c(2,1), mar=c(3,2,1,1))
# boxplot(DEDD10.TMM$raw.MSEs$noDD.noDE, cex.axis=0.8)
# boxplot(DEDD10.TMM$raw.MSEs$noDD.DE, cex.axis=0.8)
# colMeans(DEDD10.TMM$raw.MSEs$noDD.noDE)
# colMeans(DEDD10.TMM$raw.MSEs$noDD.DE)
# Exclude edgeR trend because apart from DEDD2 it's way worse than the rest and so makes 
# boxplots more difficult to compare for other methods.

par(mfrow=c(4,1), mar=c(3,2,1,1))
boxplot(cbind(nodiff2.TMM$raw.MSEs, nodiff2.DESeqnorm$raw.MSEs)); abline(v=11.5)
boxplot(cbind(nodiff5.TMM$raw.MSEs, nodiff5.DESeqnorm$raw.MSEs)); abline(v=11.5)
boxplot(cbind(nodiff10.TMM$raw.MSEs, nodiff10.DESeqnorm$raw.MSEs)); abline(v=11.5)
boxplot(cbind(nodiff20.TMM$raw.MSEs, nodiff20.DESeqnorm$raw.MSEs)); abline(v=11.5)
# Exclude edgeR trend because always really bad, DSS trend because usually missing for TMM 
# and similar to notrend for larger samples, and HMs on half samples
library(RColorBrewer)
colours5 <- brewer.pal(5,"Accent")
par(mfrow=c(2,2), mar=c(1,2,3,1))
boxplot(cbind(nodiff2.TMM$raw.MSEs[,-c(2,4,6,7,9,10)], nodiff2.DESeqnorm$raw.MSEs[,-c(2,4,6,7,9,10)]), 
        col=colours5, xaxt="n", pch=20, ylim=c(0,0.5), 
        main=paste0("MSE for dispersion estimation, 4 samples, 20,000 genes\n", 
                    "TMM (left) or DESeq (right) normalisation"))
abline(h=seq(0,0.45,0.05), col='lightgrey'); lines(c(5.5,5.5), c(0,0.45), col='lightgrey')
legend("top", fill=colours5, ncol=5, bty="n", cex=1.1, legend=c('edgeR', 'DESeq2', 'DSS', 'expHM', 'lnHM'))
boxplot(cbind(nodiff5.TMM$raw.MSEs[,-c(2,4,6,7,9,10)], nodiff5.DESeqnorm$raw.MSEs[,-c(2,4,6,7,9,10)]), 
        col=colours5, xaxt="n", pch=20, ylim=c(0,0.13), 
        main=paste0("MSE for dispersion estimation, 10 samples, 20,000 genes\n", 
                    "TMM (left) or DESeq (right) normalisation"))
abline(h=seq(0,0.12,0.01), col='lightgrey'); lines(c(5.5,5.5), c(0,0.12), col='lightgrey')
legend("top", fill=colours5, ncol=5, bty="n", cex=1.1, legend=c('edgeR', 'DESeq2', 'DSS', 'expHM', 'lnHM'))
boxplot(cbind(nodiff10.TMM$raw.MSEs[,-c(2,4,6,7,9,10)], nodiff10.DESeqnorm$raw.MSEs[,-c(2,4,6,7,9,10)]), 
        col=colours5, xaxt="n", pch=20, ylim=c(0,0.065), 
        main=paste0("MSE for dispersion estimation, 20 samples, 20,000 genes\n", 
                    "TMM (left) or DESeq (right) normalisation"))
abline(h=seq(0,0.06,0.005), col='lightgrey'); lines(c(5.5,5.5), c(0,0.06), col='lightgrey')
legend("top", fill=colours5, ncol=5, bty="n", cex=1.1, legend=c('edgeR', 'DESeq2', 'DSS', 'expHM', 'lnHM'))
boxplot(cbind(nodiff20.TMM$raw.MSEs[,-c(2,4,6,7,9,10)], nodiff20.DESeqnorm$raw.MSEs[,-c(2,4,6,7,9,10)]), 
        col=colours5, xaxt="n", pch=20, ylim=c(0,0.031), 
        main=paste0("MSE for dispersion estimation, 40 samples, 20,000 genes\n", 
                    "TMM (left) or DESeq (right) normalisation"))
abline(h=seq(0,0.028,0.002), col='lightgrey'); lines(c(5.5,5.5), c(0,0.028), col='lightgrey')
legend("top", fill=colours5, ncol=5, bty="n", cex=1.1, legend=c('edgeR', 'DESeq2', 'DSS', 'expHM', 'lnHM'))
round(rbind(colMeans(cbind(nodiff2.TMM$raw.MSEs, nodiff2.DESeqnorm$raw.MSEs)), 
            colMeans(cbind(nodiff5.TMM$raw.MSEs, nodiff5.DESeqnorm$raw.MSEs)), 
            colMeans(cbind(nodiff10.TMM$raw.MSEs, nodiff10.DESeqnorm$raw.MSEs)), 
            colMeans(cbind(nodiff20.TMM$raw.MSEs, nodiff20.DESeqnorm$raw.MSEs))), 3)
# edgeR.tag better with TMM, edgeR.trend very slightly better with TMM, all others better with DESeqnorm.
# lnHM best for all sample sizes, then edgeR.tag, DESeq2, expHM, DSS (but expHM < DESeq2 for 4 samples).
# Difference between lnHM and edgeR.tag gets smaller as sample size increases; nearly identical for 40.
# (Remember number of samples is double what names suggest as there are two groups.)

# With such a small improvement in dispersion estimates, is there any point in trying my estimates in 
# other DE methods? Possibly still is for very small samples, i.e. 2 per group, maybe 5.


par(mfrow=c(4,1), mar=c(3,2,1,1))
boxplot(cbind(DE2.TMM$raw.MSEs$noDE, DE2.DESeqnorm$raw.MSEs$noDE)); abline(v=11.5)
boxplot(cbind(DE5.TMM$raw.MSEs$noDE, DE5.DESeqnorm$raw.MSEs$noDE)); abline(v=11.5)
boxplot(cbind(DE10.TMM$raw.MSEs$noDE, DE10.DESeqnorm$raw.MSEs$noDE)); abline(v=11.5)
boxplot(cbind(DE20.TMM$raw.MSEs$noDE, DE20.DESeqnorm$raw.MSEs$noDE)); abline(v=11.5)
par(mfrow=c(4,1), mar=c(3,2,1,1))
boxplot(cbind(DE2.TMM$raw.MSEs$DE, DE2.DESeqnorm$raw.MSEs$DE)); abline(v=11.5)
boxplot(cbind(DE5.TMM$raw.MSEs$DE, DE5.DESeqnorm$raw.MSEs$DE)); abline(v=11.5)
boxplot(cbind(DE10.TMM$raw.MSEs$DE, DE10.DESeqnorm$raw.MSEs$DE)); abline(v=11.5)
boxplot(cbind(DE20.TMM$raw.MSEs$DE, DE20.DESeqnorm$raw.MSEs$DE)); abline(v=11.5)
# Exclude edgeR trend because always really bad and DSS trend because usually missing for TMM, 
# and per-group HMs for no DE
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual' & 
                                        brewer.pal.info$colorblind == T,]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# pie(rep(1,n), col=col_vector)
colours9 <- col_vector[c(1:6,8:10)]
par(mfrow=c(2,2), mar=c(1,2,3,1))
boxplot(cbind(DE2.TMM$raw.MSEs$noDE[,-c(2,4,6,7,9,10)], DE2.DESeqnorm$raw.MSEs$noDE[,-c(2,4,6,7,9,10)]), 
        col=colours5, xaxt="n", pch=20, ylim=c(0,0.45), 
        main=paste0("MSE for dispersion estimation, 2 samples per group\n", 
                    "Non-DE genes; TMM (left) or DESeq (right) normalisation"))
abline(h=seq(0,0.4,0.05), col='lightgrey'); lines(c(5.5,5.5), c(0,0.4), col='lightgrey')
legend("top", fill=colours5, ncol=5, bty="n", cex=1.1, legend=c('edgeR', 'DESeq2', 'DSS', 'expHM', 'lnHM'))
boxplot(cbind(DE5.TMM$raw.MSEs$noDE[,-c(2,4,6,7,9,10)], DE5.DESeqnorm$raw.MSEs$noDE[,-c(2,4,6,7,9,10)]), 
        col=colours5, xaxt="n", pch=20, ylim=c(0,0.16), 
        main=paste0("MSE for dispersion estimation, 5 samples per group\n", 
                    "Non-DE genes; TMM (left) or DESeq (right) normalisation"))
abline(h=seq(0,0.14,0.02), col='lightgrey'); lines(c(5.5,5.5), c(0,0.14), col='lightgrey')
legend("top", fill=colours5, ncol=5, bty="n", cex=1.1, legend=c('edgeR', 'DESeq2', 'DSS', 'expHM', 'lnHM'))
boxplot(cbind(DE10.TMM$raw.MSEs$noDE[,-c(2,4,6,7,9,10)], DE10.DESeqnorm$raw.MSEs$noDE[,-c(2,4,6,7,9,10)]), 
        col=colours5, xaxt="n", pch=20, ylim=c(0,0.065), 
        main=paste0("MSE for dispersion estimation, 10 samples per group\n", 
                    "Non-DE genes; TMM (left) or DESeq (right) normalisation"))
abline(h=seq(0,0.06,0.005), col='lightgrey'); lines(c(5.5,5.5), c(0,0.06), col='lightgrey')
legend("top", fill=colours5, ncol=5, bty="n", cex=1.1, legend=c('edgeR', 'DESeq2', 'DSS', 'expHM', 'lnHM'))
boxplot(cbind(DE20.TMM$raw.MSEs$noDE[,-c(2,4,6,7,9,10)], DE20.DESeqnorm$raw.MSEs$noDE[,-c(2,4,6,7,9,10)]), 
        col=colours5, xaxt="n", pch=20, ylim=c(0,0.035), 
        main=paste0("MSE for dispersion estimation, 20 samples per group\n", 
                    "Non-DE genes; TMM (left) or DESeq (right) normalisation"))
abline(h=seq(0,0.032,0.004), col='lightgrey'); lines(c(5.5,5.5), c(0,0.032), col='lightgrey')
legend("top", fill=colours5, ncol=5, bty="n", cex=1.1, legend=c('edgeR', 'DESeq2', 'DSS', 'expHM', 'lnHM'))
round(rbind(colMeans(cbind(DE2.TMM$raw.MSEs$noDE, DE2.DESeqnorm$raw.MSEs$noDE)), 
            colMeans(cbind(DE5.TMM$raw.MSEs$noDE, DE5.DESeqnorm$raw.MSEs$noDE)), 
            colMeans(cbind(DE10.TMM$raw.MSEs$noDE, DE10.DESeqnorm$raw.MSEs$noDE)), 
            colMeans(cbind(DE20.TMM$raw.MSEs$noDE, DE20.DESeqnorm$raw.MSEs$noDE))), 3)

par(mfrow=c(2,2), mar=c(1,2,3,1))
boxplot(cbind(DE2.TMM$raw.MSEs$DE[,-c(2,4)], DE2.DESeqnorm$raw.MSEs$DE[,-c(2,4)]), 
        col=colours9, xaxt="n", pch=20, ylim=c(0,3.7), 
        main=paste0("MSE for dispersion estimation, 2 samples per group\n", 
                    "DE genes; TMM (left) or DESeq (right) normalisation"))
abline(h=seq(0,3.2,0.2), col='lightgrey'); lines(c(9.5,9.5), c(0,3.2), col='lightgrey')
legend("topleft", fill=colours9, ncol=6, bty="n", cex=1.1, 
       legend=c('edgeR', 'DESeq2', 'DSS', 'expHM1', 'expHM2', 'expHM', 'lnHM1', 'lnHM2', 'lnHM'))
boxplot(cbind(DE5.TMM$raw.MSEs$DE[,-c(2,4)], DE5.DESeqnorm$raw.MSEs$DE[,-c(2,4)]), 
        col=colours9, xaxt="n", pch=20, ylim=c(0,0.26), 
        main=paste0("MSE for dispersion estimation, 5 samples per group\n", 
                    "DE genes; TMM (left) or DESeq (right) normalisation"))
abline(h=seq(0,0.22,0.02), col='lightgrey'); lines(c(9.5,9.5), c(0,0.22), col='lightgrey')
legend("topleft", fill=colours9, ncol=6, bty="n", cex=1.1, 
       legend=c('edgeR', 'DESeq2', 'DSS', 'expHM1', 'expHM2', 'expHM', 'lnHM1', 'lnHM2', 'lnHM'))
boxplot(cbind(DE10.TMM$raw.MSEs$DE[,-c(2,4)], DE10.DESeqnorm$raw.MSEs$DE[,-c(2,4)]), 
        col=colours9, xaxt="n", pch=20, ylim=c(0,0.185), 
        main=paste0("MSE for dispersion estimation, 10 samples per group\n", 
                    "DE genes; TMM (left) or DESeq (right) normalisation"))
abline(h=seq(0,0.16,0.02), col='lightgrey'); lines(c(9.5,9.5), c(0,0.16), col='lightgrey')
legend("topleft", fill=colours9, ncol=6, bty="n", cex=1.1, 
       legend=c('edgeR', 'DESeq2', 'DSS', 'expHM1', 'expHM2', 'expHM', 'lnHM1', 'lnHM2', 'lnHM'))
boxplot(cbind(DE20.TMM$raw.MSEs$DE[,-c(2,4)], DE20.DESeqnorm$raw.MSEs$DE[,-c(2,4)]), 
        col=colours9, xaxt="n", pch=20, ylim=c(0,0.13), 
        main=paste0("MSE for dispersion estimation, 20 samples per group\n", 
                    "DE genes; TMM (left) or DESeq (right) normalisation"))
abline(h=seq(0,0.11,0.01), col='lightgrey'); lines(c(9.5,9.5), c(0,0.11), col='lightgrey')
legend("topleft", fill=colours9, ncol=6, bty="n", cex=1.1, 
       legend=c('edgeR', 'DESeq2', 'DSS', 'expHM1', 'expHM2', 'expHM', 'lnHM1', 'lnHM2', 'lnHM'))
round(rbind(colMeans(cbind(DE2.TMM$raw.MSEs$DE, DE2.DESeqnorm$raw.MSEs$DE)), 
            colMeans(cbind(DE5.TMM$raw.MSEs$DE, DE5.DESeqnorm$raw.MSEs$DE)), 
            colMeans(cbind(DE10.TMM$raw.MSEs$DE, DE10.DESeqnorm$raw.MSEs$DE)), 
            colMeans(cbind(DE20.TMM$raw.MSEs$DE, DE20.DESeqnorm$raw.MSEs$DE))), 3)
# No point in discussing dispersion estimates when they're much worse than edgeR.tag and DESeq 
# for DE genes. Need to either find a way to improve estimates, or just not talk about dispersion 
# estimation as a point in favour of my model in paper. Using per-group estimates doesn't really 
# help as apart from the uncertainty in deciding which genes to use per-group estimates for, they're 
# still worse than edgeR.tag and DESeq overall estimates.
# 
# Could look at edgeR, DESeq methods to see how they manage to get such good estimates for DE genes, 
# but not worth spending much time on that until I have a paper written. For edgeR it will probably 
# be explicit shrinkage, for DESeq not sure, think they still do some sort of regression.
