library(here)
library(coda)
library(HDInterval)
library(compcodeR)
library(edgeR)
library(MDSeq)
source(here('scripts','2019-04-03_exponential_hmm_adaptive_proposals_three_chains_function.R'))
source(here('scripts','2019-03-27_lognormal_hmm_adaptive_proposals_three_chains_function.R'))
source(here('scripts','2019-06-26_hpd_tail_prob_function.R'))
source(here('scripts','2019-05-17_compData_diff_disp_functions.R'))

samples.per.cond <- 2
group <- factor(c(rep(1,samples.per.cond), rep(2,samples.per.cond)))
cores <- detectCores()

for (i in 1:50) {
  # Generate filtered data
  #counts <- simulate.DD.only.data(dataset='DDonly', n.vars=20000, samples.per.cond=samples.per.cond)
  filename <- paste0('DD', samples.per.cond, '.', i)
  #counts <- readRDS(paste0(filename,'.rds))
  counts <- readRDS(here('Data',paste0(filename,'.rds')))
  DD <- counts@variable.annotations$differential.dispersion
  lfcd1 <- abs(counts@variable.annotations$truelog2fcdispersion) > 1
  lfcd2 <- abs(counts@variable.annotations$truelog2fcdispersion) > 2
  
  # Normalise and create DGEList object
  nf <- calcNormFactors(counts@count.matrix)
  norm <- t(t(counts@count.matrix) / nf)

  ## MDSeq
  contrasts <- get.model.matrix(group)
  fit.zi.MDSeq <- MDSeq(norm, contrast=contrasts, mc.cores=cores-1)
  fit.nozi.MDSeq <- MDSeq(norm, contrast=contrasts, test.ZI=F, mc.cores=cores-1)
  res.zi.MDSeq <- extract.ZIMD(fit.zi.MDSeq, compare=list(A="1",B="2"))
  res.zi.lfc1.MDSeq <- extract.ZIMD(fit.zi.MDSeq, compare=list(A="1",B="2"), log2FC.threshold=1)
  res.zi.lfc2.MDSeq <- extract.ZIMD(fit.zi.MDSeq, compare=list(A="1",B="2"), log2FC.threshold=2)
  res.nozi.MDSeq <- extract.ZIMD(fit.nozi.MDSeq, compare=list(A="1",B="2"))
  res.nozi.lfc1.MDSeq <- extract.ZIMD(fit.nozi.MDSeq, compare=list(A="1",B="2"), log2FC.threshold=1)
  res.nozi.lfc2.MDSeq <- extract.ZIMD(fit.nozi.MDSeq, compare=list(A="1",B="2"), log2FC.threshold=2)
  p.disp.zi.MDSeq <- res.zi.MDSeq$Pvalue.dispersion
  p.disp.zi.lfc1.MDSeq <- res.zi.lfc1.MDSeq$Pvalue.dispersion
  p.disp.zi.lfc2.MDSeq <- res.zi.lfc2.MDSeq$Pvalue.dispersion
  p.disp.nozi.MDSeq <- res.nozi.MDSeq$Pvalue.dispersion
  p.disp.nozi.lfc1.MDSeq <- res.nozi.lfc1.MDSeq$Pvalue.dispersion
  p.disp.nozi.lfc2.MDSeq <- res.nozi.lfc2.MDSeq$Pvalue.dispersion
  rm(list=c('contrasts', 'fit.zi.MDSeq', 'fit.nozi.MDSeq', 'res.zi.MDSeq', 
            'res.zi.lfc1.MDSeq', 'res.zi.lfc2.MDSeq', 'res.nozi.MDSeq', 'res.nozi.lfc1.MDSeq', 
            'res.nozi.lfc2.MDSeq'))
  
  ## expHM
  expHM <- exp_hmm_adapt_3_chains(counts=t(norm), groups=group)
  prob.expHM <- colMeans(as.matrix(expHM$indicators))
  post.prop.expHM <- mean(as.matrix(expHM$proportion))
  mean.diff.expHM <- as.matrix(expHM$means1) - as.matrix(expHM$means2)
  log.mean.diff.expHM <- log(as.matrix(expHM$means1)) - log(as.matrix(expHM$means2))
  disp.diff.expHM <- as.matrix(expHM$disps1) - as.matrix(expHM$disps2)
  log.disp.diff.expHM <- log(as.matrix(expHM$disps1)) - log(as.matrix(expHM$disps2))
  p.disp.expHM <- apply(disp.diff.expHM,2,hpd.pval, m=0)
  p.ldisp.expHM <- apply(log.disp.diff.expHM,2,hpd.pval, m=0)
  p.disp.lfc1.expHM <- apply(log.disp.diff.expHM,2,hpd.pval, m=1)
  p.disp.lfc2.expHM <- apply(log.disp.diff.expHM,2,hpd.pval, m=2)
  rm(list=c('expHM', 'mean.diff.expHM', 'log.mean.diff.expHM', 
            'disp.diff.expHM', 'log.disp.diff.expHM'))
  
  ## lnHM
  lnHM <- ln_hmm_adapt_3_chains(counts=t(norm), groups=group)
  prob.lnHM <- colMeans(as.matrix(lnHM$indicators))
  post.prop.lnHM <- mean(as.matrix(lnHM$proportion))
  mean.diff.lnHM <- as.matrix(lnHM$means1) - as.matrix(lnHM$means2)
  log.mean.diff.lnHM <- log(as.matrix(lnHM$means1)) - log(as.matrix(lnHM$means2))
  disp.diff.lnHM <- as.matrix(lnHM$disps1) - as.matrix(lnHM$disps2)
  log.disp.diff.lnHM <- log(as.matrix(lnHM$disps1)) - log(as.matrix(lnHM$disps2))
  p.disp.lnHM <- apply(disp.diff.lnHM,2,hpd.pval)
  p.ldisp.lnHM <- apply(log.disp.diff.lnHM,2,hpd.pval)
  p.disp.lfc1.lnHM <- apply(log.disp.diff.lnHM,2,hpd.pval, m=1)
  p.disp.lfc2.lnHM <- apply(log.disp.diff.lnHM,2,hpd.pval, m=2)
  rm(list=c('lnHM', 'mean.diff.lnHM', 'log.mean.diff.lnHM', 
            'disp.diff.lnHM', 'log.disp.diff.lnHM'))
  
  results <- list(data = counts, 
                  DD = DD, 
                  lfcd1 = lfcd1, 
                  lfcd2 = lfcd2, 
                  p.disp.zi.MDSeq = p.disp.zi.MDSeq, 
                  p.disp.zi.lfc1.MDSeq = p.disp.zi.lfc1.MDSeq, 
                  p.disp.zi.lfc2.MDSeq = p.disp.zi.lfc2.MDSeq, 
                  p.disp.nozi.MDSeq = p.disp.nozi.MDSeq, 
                  p.disp.nozi.lfc1.MDSeq = p.disp.nozi.lfc1.MDSeq, 
                  p.disp.nozi.lfc2.MDSeq = p.disp.nozi.lfc2.MDSeq, 
                  prob.expHM = prob.expHM, 
                  prop.expHM = post.prop.expHM, 
                  p.disp.expHM = p.disp.expHM, 
                  p.ldisp.expHM = p.ldisp.expHM, 
                  p.disp.lfc1.expHM = p.disp.lfc1.expHM, 
                  p.disp.lfc2.expHM = p.disp.lfc2.expHM, 
                  prob.lnHM = prob.lnHM, 
                  prop.lnHM = post.prop.lnHM, 
                  p.disp.lnHM = p.disp.lnHM, 
                  p.ldisp.lnHM = p.ldisp.lnHM, 
                  p.disp.lfc1.lnHM = p.disp.lfc1.lnHM, 
                  p.disp.lfc2.lnHM = p.disp.lfc2.lnHM)
  
  filename <- paste0('results.', filename, '.rds')
  saveRDS(results, file=here(filename))

  rm(list=c('counts', 'DD', 'lfcd1', 'lfcd2', 'nf', 'norm', 'filename', 'results'))
}

filename <- paste0('sessionInfo.DD', samples.per.cond, '.rds')
saveRDS(sessionInfo(), file=here(filename))
