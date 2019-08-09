library(here)
library(coda)
library(compcodeR)
library(edgeR)
library(DESeq2)
library(DSS)
source(here('scripts','2019-04-03_exponential_hmm_adaptive_proposals_three_chains_function.R'))
source(here('scripts','2019-03-27_lognormal_hmm_adaptive_proposals_three_chains_function.R'))
source(here('scripts','2019-05-03_bfdr_function.R'))
source(here('scripts','2019-05-03_hpd_tail_prob_function.R'))
source(here('scripts','2019-05-07_symmetric_tail_prob_function.R'))
source(here('scripts','2019-05-17_compData_diff_disp_functions.R'))


# Results vectors
MSE.disp.expHM.nodiff.2 <- numeric(5)
MSE.disp.lnHM.nodiff.2 <- numeric(5)
MSE.disp.trend.edgeR.nodiff.2 <- numeric(5)
MSE.disp.tag.edgeR.nodiff.2 <- numeric(5)
MSE.disp.notrend.DSS.nodiff.2 <- numeric(5)
MSE.disp.trend.DSS.nodiff.2 <- numeric(5)
MSE.disp.DESeq.nodiff.2 <- numeric(5)

for (i in 1:5) {
  # Generate filtered data
  samples.per.cond <- 2
  counts.nodiff.2 <- generateSyntheticData(dataset='noDEDD.2', n.vars=2e4, samples.per.cond=samples.per.cond, 
                                           n.diffexp=0, filter.threshold.mediancpm=0.5)
  nf.nodiff.2 <- calcNormFactors(counts.nodiff.2@count.matrix)
  norm.nodiff.2 <- t(t(counts.nodiff.2@count.matrix) / nf.nodiff.2)
  disps.nodiff.2 <- counts.nodiff.2@variable.annotations$truedispersions.S1
  means.nodiff.2 <- counts.nodiff.2@variable.annotations$truemeans.S1
  group <- factor(c(1,1,2,2))
  
  # expHM
  expHM.nodiff.2 <- exp_hmm_adapt_3_chains(counts=t(norm.nodiff.2), groups=group)
  disps.expHM.nodiff.2 <- colMeans(as.matrix(expHM.nodiff.2$disps0))
  means.expHM.nodiff.2 <- colMeans(as.matrix(expHM.nodiff.2$means0))
  MSE.disp.expHM.nodiff.2[i] <- mean((disps.expHM.nodiff.2-disps.nodiff.2)^2)

  # lnHM
  lnHM.nodiff.2 <- ln_hmm_adapt_3_chains(counts=t(norm.nodiff.2), groups=group)
  disps.lnHM.nodiff.2 <- colMeans(as.matrix(lnHM.nodiff.2$disps0))
  means.lnHM.nodiff.2 <- colMeans(as.matrix(lnHM.nodiff.2$means0))
  MSE.disp.lnHM.nodiff.2[i] <- mean((disps.lnHM.nodiff.2-disps.nodiff.2)^2)

  # edgeR
  nodiff.2.edgeR <- DGEList(counts=counts.nodiff.2@count.matrix, norm.factors=nf.nodiff.2, group=group)
  design.edgeR <- model.matrix(~group)
  nodiff.2.edgeR <- estimateDisp(nodiff.2.edgeR, design.edgeR)
  qlfit.edgeR <- glmQLFit(nodiff.2.edgeR, design.edgeR)
  lrfit.edgeR <- glmFit(nodiff.2.edgeR, design.edgeR)
  disps.trend.edgeR.nodiff.2 <- nodiff.2.edgeR$trended.dispersion
  disps.tag.edgeR.nodiff.2 <- nodiff.2.edgeR$tagwise.dispersion
  MSE.disp.trend.edgeR.nodiff.2[i] <- mean((disps.trend.edgeR.nodiff.2 - disps.nodiff.2)^2)
  MSE.disp.tag.edgeR.nodiff.2[i] <- mean((disps.tag.edgeR.nodiff.2 - disps.nodiff.2)^2)
  
  # DESeq2
  nodiff.2.DESeq = DESeqDataSetFromMatrix(countData=counts.nodiff.2@count.matrix, colData=data.frame(group), 
                                          design=~group)
  nodiff.2.DESeq = estimateSizeFactors(nodiff.2.DESeq)
  nodiff.2.DESeq$sizeFactor <- nf.nodiff.2
  nodiff.2.DESeq <- estimateDispersions(nodiff.2.DESeq)
  disps.DESeq.nodiff.2 <- mcols(nodiff.2.DESeq)$dispersion
  MSE.disp.DESeq.nodiff.2[i] <- mean((disps.DESeq.nodiff.2 - disps.nodiff.2)^2)
  
  # DSS
  nodiff.2.DSS <- newSeqCountSet(counts=matrix(counts.nodiff.2@count.matrix, ncol=4), 
                                 designs=c(1,1,2,2), normalizationFactor=nf.nodiff.2)
  nodiff.2.notrend.DSS <- estDispersion(nodiff.2.DSS)
  nodiff.2.trend.DSS <- estDispersion(nodiff.2.DSS, trend=T)
  disps.notrend.DSS.nodiff.2 <- nodiff.2.notrend.DSS@dispersion
  disps.trend.DSS.nodiff.2 <- nodiff.2.trend.DSS@dispersion
  MSE.disp.notrend.DSS.nodiff.2[i] <- mean((disps.notrend.DSS.nodiff.2 - disps.nodiff.2)^2)
  MSE.disp.trend.DSS.nodiff.2[i] <- mean((disps.trend.DSS.nodiff.2 - disps.nodiff.2)^2)
}

MSE.disp.list.nodiff.2 <- list('expHM' = MSE.disp.expHM.nodiff.2, 
                               'lnHM' = MSE.disp.lnHM.nodiff.2, 
                               'edgeR.tag' = MSE.disp.tag.edgeR.nodiff.2,
                               'edgeR.trend' = MSE.disp.trend.edgeR.nodiff.2, 
                               'DESeq2' = MSE.disp.DESeq.nodiff.2, 
                               'DSS.notrend' = MSE.disp.notrend.DSS.nodiff.2, 
                               'DSS.trend' = MSE.disp.trend.DSS.nodiff.2)
filename <- paste0('MSE.disp.nodiff.2.', format(Sys.time(), "%Y-%m-%d_%H%M%S"), '.rds')
saveRDS(MSE.disp.list.nodiff.2, file=here(filename))


# Import results from multiple runs from files
MSElist1 <- readRDS(here('MSE.disp.nodiff.2.2019-05-27_090228.rds'))
MSElist2 <- readRDS(here('MSE.disp.nodiff.2.2019-05-27_212928.rds'))
MSElist3 <- readRDS(here('MSE.disp.nodiff.2.2019-05-28_001239.rds'))
MSElist4 <- readRDS(here('MSE.disp.nodiff.2.2019-05-28_005421.rds'))
MSElist5 <- readRDS(here('MSE.disp.nodiff.2.2019-05-29_174609.rds'))
MSEresults <- mapply(c, MSElist1, MSElist2, MSElist3, MSElist4)
boxplot(MSEresults)
colMeans(MSEresults)
#     expHM        lnHM   edgeR.tag edgeR.trend      DESeq2 DSS.notrend   DSS.trend 
# 0.1906832   0.1716861   0.1841471   0.2504281   0.2700229   0.2895091   0.2551108 

saveRDS(MSEresults, file=here('Results','2019-05-30_MSE_disp_comparison_4_samples_25_runs.rds'))

MSEdispresults <- readRDS(here('Results','2019-05-30_MSE_disp_comparison_4_samples_25_runs.rds'))
library(RColorBrewer)
par(mfrow=c(1,2), mgp=c(0.8,0.5,0), mar=c(1,3,2,1))
colours5 <- brewer.pal(5,"Accent")
colours7 <- brewer.pal(7,"Accent")
boxplot(MSEdispresults, xaxt='n', yaxt='n', main="Dispersion",cex.main=2, col=colours7, 
        ylab='MSE', cex.lab=2)
legend(x='topleft', legend=c('Hierarchical model v1', 'Hierarchical model v2', 
                             'edgeR tag', 'edgeR trend', 'DESeq2', 
                             'DSS no trend', 'DSS trend'), fill=colours7)
colMeans(MSEdispresults)
apply(MSEdispresults,2,median)


