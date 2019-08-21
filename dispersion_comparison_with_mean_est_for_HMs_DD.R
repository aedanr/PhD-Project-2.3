library(here)
library(coda)
library(compcodeR)
library(edgeR)
library(DESeq2)
library(DSS)
source(here('scripts','2019-04-03_exponential_hmm_adaptive_proposals_three_chains_function.R'))
source(here('scripts','2019-03-27_lognormal_hmm_adaptive_proposals_three_chains_function.R'))

samples.per.cond <- 2
group <- factor(c(rep(1,samples.per.cond), rep(2,samples.per.cond)))
design <- model.matrix(~group)

for (i in 1:50) {
  # Load data
  filename <- paste0('DD', samples.per.cond, '.', i)
  #counts <- readRDS(paste0(filename,'.rds'))
  counts <- readRDS(here('Data',paste0(filename,'.rds')))
  DD <- counts@variable.annotations$differential.expression

  # Normalise and create DGEList object
  nf <- calcNormFactors(counts@count.matrix)
  norm <- t(t(counts@count.matrix) / nf)
  dge <- DGEList(counts=counts@count.matrix, norm.factors=nf, group=group)

  ## edgeR
  dat.edgeR <- estimateDisp(dge, design)
  disps.edgeR.trended <- dat.edgeR$trended.dispersion
  disps.edgeR.tagwise <- dat.edgeR$tagwise.dispersion

  ## DESeq2
  dat.DESeq <- DESeqDataSetFromMatrix(countData=counts@count.matrix, colData=data.frame(group), 
                                      design=~group)
  dat.DESeq <- estimateSizeFactors(dat.DESeq)
  dat.DESeq$sizeFactor <- nf
  dat.DESeq <- estimateDispersions(dat.DESeq)
  disps.DESeq <- mcols(dat.DESeq)$dispersion
  
  ## DSS
  dat.DSS <- newSeqCountSet(counts=matrix(counts@count.matrix, ncol=length(group)), 
                            designs=as.numeric(group), normalizationFactor=nf)
  notrend.DSS <- estDispersion(dat.DSS)
  trend.DSS <- estDispersion(dat.DSS, trend=T)
  disps.DSS.notrend <- notrend.DSS@dispersion
  disps.DSS.trend <- trend.DSS@dispersion

  ## expHM
  expHM <- exp_hmm_adapt_3_chains(counts=t(norm), groups=group)
  disps.expHM <- colMeans(as.matrix(expHM$disps0))
  disps.expHM.1 <- colMeans(as.matrix(expHM$disps1))
  disps.expHM.2 <- colMeans(as.matrix(expHM$disps2))
  means.expHM <- colMeans(as.matrix(expHM$means0))
  means.expHM.1 <- colMeans(as.matrix(expHM$means1))
  means.expHM.2 <- colMeans(as.matrix(expHM$means2))
  
  # lnHM
  lnHM <- ln_hmm_adapt_3_chains(counts=t(norm), groups=group)
  disps.lnHM <- colMeans(as.matrix(lnHM$disps0))
  disps.lnHM.1 <- colMeans(as.matrix(lnHM$disps1))
  disps.lnHM.2 <- colMeans(as.matrix(lnHM$disps2))
  means.lnHM <- colMeans(as.matrix(lnHM$means0))
  means.lnHM.1 <- colMeans(as.matrix(lnHM$means1))
  means.lnHM.2 <- colMeans(as.matrix(lnHM$means2))
  
  results <- list(data = counts, 
                  DD = DD, 
                  disps.edgeR.trended = disps.edgeR.trended, 
                  disps.edgeR.tagwise = disps.edgeR.tagwise, 
                  disps.DESeq = disps.DESeq, 
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
  saveRDS(results, file=here(filename))
  
  rm(list=c('counts', 'DD', 'nf', 'norm', 'dge', 
            'dat.edgeR', 'dat.DESeq', 'notrend.DSS', 'trend.DSS',  'dat.DSS', 
            'filename', 'results'))
}

filename <- paste0('sessionInfo.dispersion.DD', samples.per.cond, '.rds')
saveRDS(sessionInfo(), file=here(filename))

