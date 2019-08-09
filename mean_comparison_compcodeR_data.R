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
MSE.mean.expHM.nodiff.2 <- numeric(5)
MSE.mean.lnHM.nodiff.2 <- numeric(5)
MSE.mean.ql.edgeR.nodiff.2 <- numeric(5)
MSE.mean.lr.edgeR.nodiff.2 <- numeric(5)
MSE.mean.DESeq.nodiff.2 <- numeric(5)

for (i in 2:5) {
  # Generate filtered data
  samples.per.cond <- 2
  counts.nodiff.2 <- generateSyntheticData(dataset='noDEDD.2', n.vars=2e4, samples.per.cond=samples.per.cond, 
                                           n.diffexp=0, filter.threshold.mediancpm=0.5)
  nf.nodiff.2 <- calcNormFactors(counts.nodiff.2@count.matrix)
  norm.nodiff.2 <- t(t(counts.nodiff.2@count.matrix) / nf.nodiff.2)
  disps.nodiff.2 <- counts.nodiff.2@variable.annotations$truedispersions.S1
  means.nodiff.2 <- counts.nodiff.2@variable.annotations$truemeans.S1
  group <- factor(c(1,1,2,2))
  sample.means.nodiff.2 <- rowMeans(counts.nodiff.2@count.matrix)
  sample.means.nodiff.2.norm <- rowMeans(norm.nodiff.2)
  
  # expHM
  expHM.nodiff.2 <- exp_hmm_adapt_3_chains(counts=t(norm.nodiff.2), groups=group)
  means.expHM.nodiff.2 <- colMeans(as.matrix(expHM.nodiff.2$means0))
  MSE.mean.expHM.nodiff.2[i] <- mean((means.expHM.nodiff.2-means.nodiff.2)^2)

  # lnHM
  lnHM.nodiff.2 <- ln_hmm_adapt_3_chains(counts=t(norm.nodiff.2), groups=group)
  disps.lnHM.nodiff.2 <- colMeans(as.matrix(lnHM.nodiff.2$disps0))
  means.lnHM.nodiff.2 <- colMeans(as.matrix(lnHM.nodiff.2$means0))
  MSE.mean.lnHM.nodiff.2[i] <- mean((means.lnHM.nodiff.2-means.nodiff.2)^2)

  # edgeR
  nodiff.2.edgeR <- DGEList(counts=counts.nodiff.2@count.matrix, norm.factors=nf.nodiff.2, group=group)
  design.edgeR <- model.matrix(~group)
  nodiff.2.edgeR <- estimateDisp(nodiff.2.edgeR, design.edgeR)
  qlfit.edgeR <- glmQLFit(nodiff.2.edgeR, design.edgeR)
  lrfit.edgeR <- glmFit(nodiff.2.edgeR, design.edgeR)
  means.ql.edgeR <- (rowMeans(qlfit.edgeR$fitted.values))
  means.lr.edgeR <- (rowMeans(lrfit.edgeR$fitted.values))
  MSE.mean.ql.edgeR.nodiff.2[i] <- mean((means.ql.edgeR - means.nodiff.2)^2)
  MSE.mean.lr.edgeR.nodiff.2[i] <- mean((means.lr.edgeR - means.nodiff.2)^2)
  
  # DESeq2
  nodiff.2.DESeq = DESeqDataSetFromMatrix(countData=counts.nodiff.2@count.matrix, colData=data.frame(group), 
                                          design=~group)
  nodiff.2.DESeq = estimateSizeFactors(nodiff.2.DESeq)
  nodiff.2.DESeq$sizeFactor <- nf.nodiff.2
  nodiff.2.DESeq <- estimateDispersions(nodiff.2.DESeq)
  nodiff.2.DESeq <- DESeq(nodiff.2.DESeq)
  means.DESeq <- rowMeans(assays(nodiff.2.DESeq)$mu)
  MSE.mean.DESeq.nodiff.2[i] <- mean((means.DESeq - means.nodiff.2)^2)
}

MSE.mean.list.nodiff.2 <- list('expHM' = MSE.mean.expHM.nodiff.2, 
                               'lnHM' = MSE.mean.lnHM.nodiff.2, 
                               'edgeR.ql' = MSE.mean.ql.edgeR.nodiff.2,
                               'edgeR.lr' = MSE.mean.lr.edgeR.nodiff.2, 
                               'DESeq2' = MSE.mean.DESeq.nodiff.2)
filename <- paste0('MSE.mean.nodiff.2.', format(Sys.time(), "%Y-%m-%d_%H%M%S"), '.rds')
saveRDS(MSE.mean.list.nodiff.2, file=here(filename))


# Import results from multiple runs from files
MSElist1 <- readRDS(here('MSE.mean.nodiff.2.2019-05-29_053425.rds'))
MSElist2 <- readRDS(here('MSE.mean.nodiff.2.2019-05-29_062229.rds'))
MSElist3 <- readRDS(here('MSE.mean.nodiff.2.2019-05-29_073344.rds'))
MSElist4 <- readRDS(here('MSE.mean.nodiff.2.2019-05-29_074834.rds'))
MSElist5 <- readRDS(here('MSE.mean.nodiff.2.2019-05-29_082231.rds'))
MSEmeanresults <- mapply(c, MSElist1, MSElist2, MSElist3, MSElist4, MSElist5)

saveRDS(MSEresults, file=here('Results','2019-05-29_MSE_mean_comparison_4_samples_25_runs.rds'))

MSEmeanresults <- readRDS(here('Results','2019-05-29_MSE_mean_comparison_4_samples_25_runs.rds'))
library(RColorBrewer)
par(mfrow=c(1,2), mgp=c(0.8,0.5,0), mar=c(1,3,2,1))
colours5 <- brewer.pal(5,"Accent")
colours7 <- brewer.pal(7,"Accent")
boxplot(MSEmeanresults, xaxt='n', yaxt='n', main="Mean", cex.main=2, col=colours5, 
        ylab='MSE',  cex.lab=2)
legend(x='topleft', legend=c('Hierarchical model v1', 'Hierarchical model v2', 
                             'edgeR QL', 'edgeR LR', 'DESeq2'), fill=colours5)
colMeans(MSEmeanresults)/1e6
apply(MSEmeanresults,2,median)/1e6



