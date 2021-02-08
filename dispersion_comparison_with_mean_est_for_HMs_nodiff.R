library(here)
library(coda)
library(compcodeR)
library(edgeR)
library(DESeq2)
library(DSS)
source(here("scripts", "2019-05-03_bfdr_function.R"))
source(here("scripts", "2020-07-23_hpd_tail_prob_function.R"))

samples.per.cond <- 20
group <- factor(c(rep(1,samples.per.cond), rep(2,samples.per.cond)))
design <- model.matrix(~group)

for (i in 29:50) {
  # Load data
  filename <- paste0('nodiff', samples.per.cond, '.', i)
  counts <- readRDS(here('Simulated data', paste0(filename, '.rds')))

  # Normalise and create DGEList object
  libsizes <- colSums(counts@count.matrix)
  nf <- calcNormFactors(counts@count.matrix, method="TMM")
  sf <- nf * libsizes / exp(mean(log(libsizes)))
  norm.counts <- t(t(counts@count.matrix) / sf)
  dge <- DGEList(counts=counts@count.matrix, norm.factors=nf, group=group)
  
  ## edgeR
  dat.edgeR <- estimateDisp(dge, design)
  disps.edgeR.trended <- dat.edgeR$trended.dispersion
  disps.edgeR.tagwise <- dat.edgeR$tagwise.dispersion
  rm(dat.edgeR)

  ## DESeq2
  dat.DESeq2 <- DESeqDataSetFromMatrix(countData=counts@count.matrix, colData=data.frame(group), 
                                      design=~group)
  sizeFactors(dat.DESeq2) <- sf
  dat.DESeq2 <- estimateDispersions(dat.DESeq2)
  disps.DESeq2 <- mcols(dat.DESeq2)$dispersion
  rm(dat.DESeq2)
  
  ## DSS
  dat.DSS <- newSeqCountSet(counts=matrix(counts@count.matrix, ncol=length(group)), 
                            designs=as.numeric(group), normalizationFactor=sf)
  notrend.DSS <- estDispersion(dat.DSS)
  trend.DSS <- estDispersion(dat.DSS, trend=T)
  disps.DSS.notrend <- notrend.DSS@dispersion
  disps.DSS.trend <- trend.DSS@dispersion
  rm(dat.DSS)

  ## expHM
  source(here("scripts", "2020-07-17_conditional_posterior_functions_exponential_hmm.R"))
  source(here("scripts", "2020-07-23_exponential_hmm_one_chain_function.R"))
  source(here("scripts", "2020-07-23_exponential_hmm_three_chains_function.R"))
  source(here("scripts", "2020-07-25_exponential_hmm_adaptive_proposals_three_chains_function.R"))
  source(here("scripts", "2020-09-02_run_expHMM.R"))
  expHM <- expHMM(counts=t(norm.counts), groups=group, return.raw.results=T)
  disps.expHM <- colMeans(as.matrix(expHM$disps0))
  disps.expHM.1 <- colMeans(as.matrix(expHM$disps1))
  disps.expHM.2 <- colMeans(as.matrix(expHM$disps2))
  means.expHM <- colMeans(as.matrix(expHM$means0))
  means.expHM.1 <- colMeans(as.matrix(expHM$means1))
  means.expHM.2 <- colMeans(as.matrix(expHM$means2))
  rm(expHM)
  
  # lnHM
  source(here("scripts", "2020-07-17_conditional_posterior_functions_lognormal_hmm.R"))
  source(here("scripts", "2020-07-23_lognormal_hmm_one_chain_function.R"))
  source(here("scripts", "2020-07-23_lognormal_hmm_three_chains_function.R"))
  source(here("scripts", "2020-07-25_lognormal_hmm_adaptive_proposals_three_chains_function.R"))
  source(here("scripts", "2020-09-02_run_lnHMM.R"))
  lnHM <- lnHMM(counts=t(norm.counts), groups=group, return.raw.results=T)
  disps.lnHM <- colMeans(as.matrix(lnHM$disps0))
  disps.lnHM.1 <- colMeans(as.matrix(lnHM$disps1))
  disps.lnHM.2 <- colMeans(as.matrix(lnHM$disps2))
  means.lnHM <- colMeans(as.matrix(lnHM$means0))
  means.lnHM.1 <- colMeans(as.matrix(lnHM$means1))
  means.lnHM.2 <- colMeans(as.matrix(lnHM$means2))
  rm(lnHM)
  
  results <- list(data = counts, 
                  disps.edgeR.trended = disps.edgeR.trended, 
                  disps.edgeR.tagwise = disps.edgeR.tagwise, 
                  disps.DESeq2 = disps.DESeq2, 
                  disps.DSS.notrend = disps.DSS.notrend, 
                  disps.DSS.trend = disps.DSS.trend, 
                  disps.expHM = disps.expHM, 
                  disps.expHM.1 = disps.expHM.1, 
                  disps.expHM.2 = disps.expHM.2, 
                  means.expHM = means.expHM, 
                  means.expHM.1 = means.expHM.1, 
                  means.expHM.2 = means.expHM.2, 
                  disps.lnHM = disps.lnHM, 
                  disps.lnHM.1 = disps.lnHM.1, 
                  disps.lnHM.2 = disps.lnHM.2, 
                  means.lnHM = means.lnHM, 
                  means.lnHM.1 = means.lnHM.1, 
                  means.lnHM.2 = means.lnHM.2)
  
  filename <- paste0('disp.results.', filename, '.rds')
  saveRDS(results, file=here('Results/Dispersion estimation Dec 2020', filename))
  
  rm(list=c('norm.counts', 'nf', 'sf', 'dge', 'notrend.DSS', 'trend.DSS', 'filename', 'results'))
}

filename <- paste0('sessionInfo.dispersion.nodiff', samples.per.cond, '.rds')
saveRDS(sessionInfo(), file=here(filename))

