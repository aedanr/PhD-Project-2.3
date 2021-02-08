library(here)
DE2 <- readRDS(here('Results/Dispersion estimation Dec 2020','mse.disp.DE2.rds'))
DE5 <- readRDS(here('Results/Dispersion estimation Dec 2020','mse.disp.DE5.rds'))
DE10 <- readRDS(here('Results/Dispersion estimation Dec 2020','mse.disp.DE10.rds'))
DE20 <- readRDS(here('Results/Dispersion estimation Dec 2020','mse.disp.DE20.rds'))
nodiff2 <- readRDS(here('Results/Dispersion estimation Dec 2020','mse.disp.nodiff2.rds'))
nodiff5 <- readRDS(here('Results/Dispersion estimation Dec 2020','mse.disp.nodiff5.rds'))
nodiff10 <- readRDS(here('Results/Dispersion estimation Dec 2020','mse.disp.nodiff10.rds'))
nodiff20 <- readRDS(here('Results/Dispersion estimation Dec 2020','mse.disp.nodiff20.rds'))

library(RColorBrewer)
colours4 <- brewer.pal(4,"Accent")
colours6 <- brewer.pal(6,"Accent")

par(mfrow=c(4,1), mar=c(3,2,1,1))
boxplot(nodiff2$raw.MSEs)
boxplot(nodiff5$raw.MSEs)
boxplot(nodiff10$raw.MSEs)
boxplot(nodiff20$raw.MSEs)
# Exclude edgeR trend because always really bad, DSS notrend because never better than trend, 
# and HMs on half samples. Only show lnHM.
# Only show raw MSEs.

nodiff2 <- nodiff2$raw.MSEs[, c(1,3,4,11)]
nodiff5 <- nodiff5$raw.MSEs[, c(1,3,4,11)]
nodiff10 <- nodiff10$raw.MSEs[, c(1,3,4,11)]
nodiff20 <- nodiff20$raw.MSEs[, c(1,3,4,11)]

par(mfrow=c(1,4), mar=c(2,4,3,2))
boxplot(nodiff2, 
        col=colours4, xaxt="n", pch=20, cex=2, 
        main="4", cex.main=3, cex.axis=3)
boxplot(nodiff5, 
        col=colours4, xaxt="n", pch=20, cex=2, 
        main="10", cex.main=3, cex.axis=3)
boxplot(nodiff10, 
        col=colours4, xaxt="n", pch=20, cex=2, 
        main="20", cex.main=3, cex.axis=3)
boxplot(nodiff20, 
        col=colours4, xaxt="n", pch=20, cex=2, 
        main="40", cex.main=3, cex.axis=3)
legend("topleft", fill=colours4, ncol=1, bty="n", cex=2.5, legend=c('edgeR', 'DESeq2', 'DSS', 'HM'))
round(rbind(colMeans(nodiff2), 
            colMeans(nodiff5), 
            colMeans(nodiff10), 
            colMeans(nodiff20)), 3)
#      edgeR.tag DESeq2 DSS.trend lnHM.overall
# [1,]     0.178  0.239     0.253        0.151
# [2,]     0.052  0.056     0.097        0.048
# [3,]     0.024  0.025     0.049        0.023
# [4,]     0.011  0.012     0.026        0.011
# lnHM better than all others, except same as edgeR for 20 samples per group; 
# much better than when using wrong normalisation adjustment

# With DE
# Keep per-group results for lnHM
DE2_all <- DE2$raw.MSEs$all_genes[, c(1,3,4,9:11)]
DE5_all <- DE5$raw.MSEs$all_genes[, c(1,3,4,9:11)]
DE10_all <- DE10$raw.MSEs$all_genes[, c(1,3,4,9:11)]
DE20_all <- DE20$raw.MSEs$all_genes[, c(1,3,4,9:11)]
DE2_DEonly <- DE2$raw.MSEs$DE[, c(1,3,4,9:11)]
DE5_DEonly <- DE5$raw.MSEs$DE[, c(1,3,4,9:11)]
DE10_DEonly <- DE10$raw.MSEs$DE[, c(1,3,4,9:11)]
DE20_DEonly <- DE20$raw.MSEs$DE[, c(1,3,4,9:11)]
DE2_noDE <- DE2$raw.MSEs$noDE[, c(1,3,4,9:11)]
DE5_noDE <- DE5$raw.MSEs$noDE[, c(1,3,4,9:11)]
DE10_noDE <- DE10$raw.MSEs$noDE[, c(1,3,4,9:11)]
DE20_noDE <- DE20$raw.MSEs$noDE[, c(1,3,4,9:11)]


par(mfrow=c(3,4), mar=c(0.5,3,3.5,1))
boxplot(DE2_all, 
        col=colours6, xaxt="n", pch=20, cex=2, ylim=c(0.1,0.7), 
        main="2 samples per group", cex.main=3, cex.axis=2.5)
boxplot(DE5_all, 
        col=colours6, xaxt="n", pch=20, cex=2, 
        main="5 samples per group", cex.main=3, cex.axis=2.5)
boxplot(DE10_all, 
        col=colours6, xaxt="n", pch=20, cex=2, 
        main="10 samples per group", cex.main=3, cex.axis=2.5)
boxplot(DE20_all, 
        col=colours6, xaxt="n", pch=20, cex=2, 
        main="20 samples per group", cex.main=3, cex.axis=2.5)

boxplot(DE2_noDE, 
        col=colours6, xaxt="n", pch=20, cex=2, ylim=c(0.1,0.71), 
        cex.main=3, cex.axis=2.5)
boxplot(DE5_noDE, 
        col=colours6, xaxt="n", pch=20, cex=2, 
        cex.main=3, cex.axis=2.5)
boxplot(DE10_noDE, 
        col=colours6, xaxt="n", pch=20, cex=2, 
        cex.main=3, cex.axis=2.5)
boxplot(DE20_noDE, 
        col=colours6, xaxt="n", pch=20, cex=2, 
        cex.main=3, cex.axis=2.5)

boxplot(DE2_DEonly, 
        col=colours6, xaxt="n", pch=20, cex=2, ylim=c(0.06,0.71), 
        cex.main=3, cex.axis=2.5)
boxplot(DE5_DEonly, 
        col=colours6, xaxt="n", pch=20, cex=2, 
        cex.main=3, cex.axis=2.5)
boxplot(DE10_DEonly, 
        col=colours6, xaxt="n", pch=20, cex=2, 
        cex.main=3, cex.axis=2.5)
boxplot(DE20_DEonly, 
        col=colours6, xaxt="n", pch=20, cex=2, #ylim=c(0.007, 0.12), 
        cex.main=3, cex.axis=2.5)
legend("topleft", fill=colours6, ncol=2, bty="n", cex=1.8, 
       legend=c('edgeR', 'DESeq2', 'DSS', 'HM (grp 1)', 'HM (grp 2)', 'HM (overall)'))

round(rbind(colMeans(DE2_all), 
            colMeans(DE5_all), 
            colMeans(DE10_all), 
            colMeans(DE20_all)), 3)
round(rbind(colMeans(DE2_noDE), 
            colMeans(DE5_noDE), 
            colMeans(DE10_noDE), 
            colMeans(DE20_noDE)), 3)
round(rbind(colMeans(DE2_DEonly), 
            colMeans(DE5_DEonly), 
            colMeans(DE10_DEonly), 
            colMeans(DE20_DEonly)), 3)


## Create and save est v true disp and density plots for one simulation per sample size ####
# Not sure whether will use or not.
for (i in c('.DE2', '.DE5', '.DE10', '.DE20')) {
        res <- readRDS(here('Results/Dispersion estimation Dec 2020', 
                            paste0('disp.results', i, '.50.rds')))
        assign(paste0('DE',i) , res$data@variable.annotations$differential.expression)
        assign(paste0('true.disps',i) , res$data@variable.annotations$truedispersions.S1)
        assign(paste0('true.disps_DE',i) , 
               res$data@variable.annotations$truedispersions.S1[which(get(paste0('DE',i)) == 1)])
        assign(paste0('true.disps_noDE',i) , 
               res$data@variable.annotations$truedispersions.S1[which(get(paste0('DE',i)) == 0)])
        assign(paste0('true.means_DE1',i) , 
               res$data@variable.annotations$truemeans.S1[which(get(paste0('DE',i)) == 1)])
        assign(paste0('true.means_DE2',i) , 
               res$data@variable.annotations$truemeans.S2[which(get(paste0('DE',i)) == 1)])
        assign(paste0('true.means_noDE',i) , 
               res$data@variable.annotations$truemeans.S1[which(get(paste0('DE',i)) == 0)])
        assign(paste0('mean.counts_DE1',i) , 
               rowMeans(res$data@count.matrix[which(get(paste0('DE',i)) == 1), 
                                              which(res$data@sample.annotations$condition == 1)]))
        assign(paste0('mean.counts_DE2',i) , 
               rowMeans(res$data@count.matrix[which(get(paste0('DE',i)) == 1), 
                                              which(res$data@sample.annotations$condition == 2)]))
        assign(paste0('mean.counts_noDE',i) , 
               rowMeans(res$data@count.matrix[which(get(paste0('DE',i)) == 0),]))
        assign(paste0('DESeq2',i), res$disps.DESeq)
        assign(paste0('DESeq2_DE',i), res$disps.DESeq[which(get(paste0('DE',i)) == 1)])
        assign(paste0('DESeq2_noDE',i), res$disps.DESeq[which(get(paste0('DE',i)) == 0)])
        assign(paste0('DSS.trend',i), res$disps.DSS.trend)
        assign(paste0('DSS.trend_DE',i), 
               res$disps.DSS.trend[which(get(paste0('DE',i)) == 1)])
        assign(paste0('DSS.trend_noDE',i), 
               res$disps.DSS.trend[which(get(paste0('DE',i)) == 0)])
        assign(paste0('DSS.trend',i), res$disps.DSS.trend)
        assign(paste0('DSS.trend_DE',i), 
               res$disps.DSS.trend[which(get(paste0('DE',i)) == 1)])
        assign(paste0('DSS.trend_noDE',i), 
               res$disps.DSS.trend[which(get(paste0('DE',i)) == 0)])
        assign(paste0('edgeR.trend',i), res$disps.edgeR.trended)
        assign(paste0('edgeR.trend_DE',i), 
               res$disps.edgeR.trended[which(get(paste0('DE',i)) == 1)])
        assign(paste0('edgeR.trend_noDE',i), 
               res$disps.edgeR.trended[which(get(paste0('DE',i)) == 0)])
        assign(paste0('edgeR.tag',i), res$disps.edgeR.tagwise)
        assign(paste0('edgeR.tag_DE',i), 
               res$disps.edgeR.tagwise[which(get(paste0('DE',i)) == 1)])
        assign(paste0('edgeR.tag_noDE',i), 
               res$disps.edgeR.tagwise[which(get(paste0('DE',i)) == 0)])
        assign(paste0('expHM.group1_est',i), res$disps.expHM.1)
        assign(paste0('expHM.group1_est_DE',i), 
               res$disps.expHM.1[which(get(paste0('DE',i)) == 1)])
        assign(paste0('expHM.group1_est_noDE',i), 
               res$disps.expHM.1[which(get(paste0('DE',i)) == 0)])
        assign(paste0('expHM.group2_est',i), res$disps.expHM.2)
        assign(paste0('expHM.group2_est_DE',i), 
               res$disps.expHM.2[which(get(paste0('DE',i)) == 1)])
        assign(paste0('expHM.group2_est_noDE',i), 
               res$disps.expHM.2[which(get(paste0('DE',i)) == 0)])
        assign(paste0('expHM.overall_est',i), res$disps.expHM)
        assign(paste0('expHM.overall_est_DE',i), 
               res$disps.expHM[which(get(paste0('DE',i)) == 1)])
        assign(paste0('expHM.overall_est_noDE',i), 
               res$disps.expHM[which(get(paste0('DE',i)) == 0)])
        assign(paste0('lnHM.group1_est',i), res$disps.lnHM.1)
        assign(paste0('lnHM.group1_est_DE',i), 
               res$disps.lnHM.1[which(get(paste0('DE',i)) == 1)])
        assign(paste0('lnHM.group1_est_noDE',i), 
               res$disps.lnHM.1[which(get(paste0('DE',i)) == 0)])
        assign(paste0('lnHM.group2_est',i), res$disps.lnHM.2)
        assign(paste0('lnHM.group2_est_DE',i), 
               res$disps.lnHM.2[which(get(paste0('DE',i)) == 1)])
        assign(paste0('lnHM.group2_est_noDE',i), 
               res$disps.lnHM.2[which(get(paste0('DE',i)) == 0)])
        assign(paste0('lnHM.overall_est',i), res$disps.lnHM)
        assign(paste0('lnHM.overall_est_DE',i), 
               res$disps.lnHM[which(get(paste0('DE',i)) == 1)])
        assign(paste0('lnHM.overall_est_noDE',i), 
               res$disps.lnHM[which(get(paste0('DE',i)) == 0)])
}

# non-DE genes
par(mfrow=c(4,4), mar=c(1.5,1.5,2,1), mgp=c(2,0.5,0))
plot(log(true.disps_noDE.DE2), log(edgeR.tag_noDE.DE2), 
     pch=20, cex=0.5, xlim=c(-6,2), ylim=c(-6,2), main='edgeR', cex.main=2)
lines(c(-6,2), c(-6,2))
legend('bottomleft', legend='2 samples per group', bty='n', cex=1.8)
plot(log(true.disps_noDE.DE2), log(DESeq2_noDE.DE2), 
     pch=20, cex=0.5, xlim=c(-6,2), ylim=c(-6,2), main='DESeq2', cex.main=2)
lines(c(-6,2), c(-6,2))
plot(log(true.disps_noDE.DE2), log(DSS.trend_noDE.DE2), 
     pch=20, cex=0.5, xlim=c(-6,2), ylim=c(-6,2), main='DSS', cex.main=2)
lines(c(-6,2), c(-6,2))
plot(log(true.disps_noDE.DE2), log(lnHM.overall_est_noDE.DE2), 
     pch=20, cex=0.5, xlim=c(-6,2), ylim=c(-6,2), main='HM', cex.main=2)
lines(c(-6,2), c(-6,2))
plot(log(true.disps_noDE.DE5), log(edgeR.tag_noDE.DE5), 
     pch=20, cex=0.5, xlim=c(-6,2), ylim=c(-6,2))
lines(c(-6,2), c(-6,2))
legend('bottomleft', legend='5 samples per group', bty='n', cex=1.8)
plot(log(true.disps_noDE.DE5), log(DESeq2_noDE.DE5), 
     pch=20, cex=0.5, xlim=c(-6,2), ylim=c(-6,2))
lines(c(-6,2), c(-6,2))
plot(log(true.disps_noDE.DE5), log(DSS.trend_noDE.DE5), 
     pch=20, cex=0.5, xlim=c(-6,2), ylim=c(-6,2))
lines(c(-6,2), c(-6,2))
plot(log(true.disps_noDE.DE5), log(lnHM.overall_est_noDE.DE5), 
     pch=20, cex=0.5, xlim=c(-6,2), ylim=c(-6,2))
lines(c(-6,2), c(-6,2))
plot(log(true.disps_noDE.DE10), log(edgeR.tag_noDE.DE10), 
     pch=20, cex=0.5, xlim=c(-6,2), ylim=c(-6,2))
lines(c(-6,2), c(-6,2))
legend('bottomleft', legend='10 samples per group', bty='n', cex=1.8)
plot(log(true.disps_noDE.DE10), log(DESeq2_noDE.DE10), 
     pch=20, cex=0.5, xlim=c(-6,2), ylim=c(-6,2))
lines(c(-6,2), c(-6,2))
plot(log(true.disps_noDE.DE10), log(DSS.trend_noDE.DE10), 
     pch=20, cex=0.5, xlim=c(-6,2), ylim=c(-6,2))
lines(c(-6,2), c(-6,2))
plot(log(true.disps_noDE.DE10), log(lnHM.overall_est_noDE.DE10), 
     pch=20, cex=0.5, xlim=c(-6,2), ylim=c(-6,2))
lines(c(-6,2), c(-6,2))
plot(log(true.disps_noDE.DE20), log(edgeR.tag_noDE.DE20), 
     pch=20, cex=0.5, xlim=c(-6,2), ylim=c(-6,2))
lines(c(-6,2), c(-6,2))
legend('bottomleft', legend='20 samples per group', bty='n', cex=1.8)
plot(log(true.disps_noDE.DE20), log(DESeq2_noDE.DE20), 
     pch=20, cex=0.5, xlim=c(-6,2), ylim=c(-6,2))
lines(c(-6,2), c(-6,2))
plot(log(true.disps_noDE.DE20), log(DSS.trend_noDE.DE20), 
     pch=20, cex=0.5, xlim=c(-6,2), ylim=c(-6,2))
lines(c(-6,2), c(-6,2))
plot(log(true.disps_noDE.DE20), log(lnHM.overall_est_noDE.DE20), 
     pch=20, cex=0.5, xlim=c(-6,2), ylim=c(-6,2))
lines(c(-6,2), c(-6,2))

# DE genes
par(mfrow=c(4,4), mar=c(1.5,1.5,2,1), mgp=c(2,0.5,0))
plot(log(true.disps_DE.DE2), log(edgeR.tag_DE.DE2), 
     pch=20, cex=0.5, xlim=c(-6,2), ylim=c(-6,2), main='edgeR', cex.main=2)
lines(c(-6,2), c(-6,2))
legend('bottomleft', legend='2 samples per group', bty='n', cex=1.8)
plot(log(true.disps_DE.DE2), log(DESeq2_DE.DE2), 
     pch=20, cex=0.5, xlim=c(-6,2), ylim=c(-6,2), main='DESeq2', cex.main=2)
lines(c(-6,2), c(-6,2))
plot(log(true.disps_DE.DE2), log(DSS.trend_DE.DE2), 
     pch=20, cex=0.5, xlim=c(-6,2), ylim=c(-6,2), main='DSS', cex.main=2)
lines(c(-6,2), c(-6,2))
plot(log(true.disps_DE.DE2), log(lnHM.overall_est_DE.DE2), 
     pch=20, cex=0.5, xlim=c(-6,2), ylim=c(-6,2), main='HM', cex.main=2)
lines(c(-6,2), c(-6,2))
plot(log(true.disps_DE.DE5), log(edgeR.tag_DE.DE5), 
     pch=20, cex=0.5, xlim=c(-6,2), ylim=c(-6,2))
lines(c(-6,2), c(-6,2))
legend('bottomleft', legend='5 samples per group', bty='n', cex=1.8)
plot(log(true.disps_DE.DE5), log(DESeq2_DE.DE5), 
     pch=20, cex=0.5, xlim=c(-6,2), ylim=c(-6,2))
lines(c(-6,2), c(-6,2))
plot(log(true.disps_DE.DE5), log(DSS.trend_DE.DE5), 
     pch=20, cex=0.5, xlim=c(-6,2), ylim=c(-6,2))
lines(c(-6,2), c(-6,2))
plot(log(true.disps_DE.DE5), log(lnHM.overall_est_DE.DE5), 
     pch=20, cex=0.5, xlim=c(-6,2), ylim=c(-6,2))
lines(c(-6,2), c(-6,2))
plot(log(true.disps_DE.DE10), log(edgeR.tag_DE.DE10), 
     pch=20, cex=0.5, xlim=c(-6,2), ylim=c(-6,2))
lines(c(-6,2), c(-6,2))
legend('bottomleft', legend='10 samples per group', bty='n', cex=1.8)
plot(log(true.disps_DE.DE10), log(DESeq2_DE.DE10), 
     pch=20, cex=0.5, xlim=c(-6,2), ylim=c(-6,2))
lines(c(-6,2), c(-6,2))
plot(log(true.disps_DE.DE10), log(DSS.trend_DE.DE10), 
     pch=20, cex=0.5, xlim=c(-6,2), ylim=c(-6,2))
lines(c(-6,2), c(-6,2))
plot(log(true.disps_DE.DE10), log(lnHM.overall_est_DE.DE10), 
     pch=20, cex=0.5, xlim=c(-6,2), ylim=c(-6,2))
lines(c(-6,2), c(-6,2))
plot(log(true.disps_DE.DE20), log(edgeR.tag_DE.DE20), 
     pch=20, cex=0.5, xlim=c(-6,2), ylim=c(-6,2))
lines(c(-6,2), c(-6,2))
legend('bottomleft', legend='20 samples per group', bty='n', cex=1.8)
plot(log(true.disps_DE.DE20), log(DESeq2_DE.DE20), 
     pch=20, cex=0.5, xlim=c(-6,2), ylim=c(-6,2))
lines(c(-6,2), c(-6,2))
plot(log(true.disps_DE.DE20), log(DSS.trend_DE.DE20), 
     pch=20, cex=0.5, xlim=c(-6,2), ylim=c(-6,2))
lines(c(-6,2), c(-6,2))
plot(log(true.disps_DE.DE20), log(lnHM.overall_est_DE.DE20), 
     pch=20, cex=0.5, xlim=c(-6,2), ylim=c(-6,2))
lines(c(-6,2), c(-6,2))


par(mfrow=c(1,2), mar=c(1.5,1.5,2,1), mgp=c(2,0.5,0))
plot(density(log(true.disps_noDE.DE20)), xlim=c(-5,2), ylim=c(0,0.6), lwd=2, 
     main='Non-DE genes', yaxt='n', cex.main=1.5)
legend("topleft", fill=c("black", "red", "orange", "blue"), bty="n", cex=1.5, border=F, 
       legend=c("True dispersions", "edgeR estimates", "DESeq2 estimates", "HM estimates"))
lines(density(log(edgeR.tag_noDE.DE20)), col='red', lwd=2)
lines(density(log(DESeq2_noDE.DE20)), col='orange', lwd=2)
lines(density(log(lnHM.overall_est_noDE.DE20)), col='blue', lwd=2)
plot(density(log(true.disps_DE.DE20)), xlim=c(-5,2), ylim=c(0,0.6), lwd=2, 
     main='DE genes', yaxt='n', cex.main=1.5)
lines(density(log(edgeR.tag_DE.DE20)), col='red', lwd=2)
lines(density(log(DESeq2_DE.DE20)), col='orange', lwd=2)
lines(density(log(lnHM.overall_est_DE.DE20)), col='blue', lwd=2)




