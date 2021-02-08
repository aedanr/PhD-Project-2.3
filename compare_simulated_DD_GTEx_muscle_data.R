library(here)
library(coda)
library(HDInterval)
library(compcodeR)
library(limma)
library(edgeR)
library(DESeq2)
library(DSS)
library(baySeq)
library(MDSeq)
library(missMethyl)

cluster <- T

if (cluster) {
  source(here('Data/scripts', '2019-04-03_exponential_hmm_adaptive_proposals_three_chains_function.R'))
  source(here('Data/scripts', '2019-03-27_lognormal_hmm_adaptive_proposals_three_chains_function.R'))
  source(here('Data/scripts', '2019-09-05_hpd_tail_prob_function.R'))
  source(here('Data/scripts', '2019-05-03_bfdr_function.R'))
} else {
  source(here('scripts', '2019-04-03_exponential_hmm_adaptive_proposals_three_chains_function.R'))
  source(here('scripts', '2019-03-27_lognormal_hmm_adaptive_proposals_three_chains_function.R'))
  source(here('scripts', '2019-09-05_hpd_tail_prob_function.R'))
  source(here('scripts', '2019-05-03_bfdr_function.R'))
}

samples.per.cond <- 2
group <- factor(c(rep(1, samples.per.cond), rep(2, samples.per.cond)))
design <- model.matrix(~group)
contrasts <- get.model.matrix(group)
# cores <- detectCores() - 1

for (i in 1:10) {
  # Data
  if (cluster) {
    folder <- paste0("Data/recount data/GTEx/muscle_", samples.per.cond, "_samples_per_group")
  } else {
    folder <- paste0("recount data/GTEx/muscle_", samples.per.cond, "_samples_per_group")
  }
  filename <- paste0("muscle_", samples.per.cond, "_set", i, "_DD")
  data <- readRDS(here(folder, paste0(filename, ".rds")))
  counts <- data$counts
  DD <- data$DDindex
  DEDD <- DD
  FC <- data$FC
  rm(data)

  # Normalise
  libsizes <- colSums(counts)
  nf.tmm <- calcNormFactors(counts, method="TMM")
  els.tmm <- nf.tmm * libsizes
  sf.tmm <- els.tmm / exp(mean(log(libsizes)))
  norm.counts.tmm <- t(t(counts) / sf.tmm)
  nf.rle <- calcNormFactors(counts, method="RLE")
  els.rle <- nf.rle * libsizes
  sf.rle <- els.rle / exp(mean(log(libsizes)))
  norm.counts.rle <- t(t(counts) / sf.rle)
  
  ## MDSeq
  # fit.MDSeq.zi.tmm <- MDSeq(counts, offsets=sf.tmm, contrast=contrasts, mc.cores=cores-1)
  # hangs on Phoenix with mc.cores
  fit.MDSeq.zi.tmm <- MDSeq(counts, offsets=sf.tmm, contrast=contrasts)
  res.MDSeq.zi.tmm <- extract.ZIMD(fit.MDSeq.zi.tmm, compare=list(A="1",B="2"))
  # fit.MDSeq.nozi.tmm <- MDSeq(counts, offsets=sf.tmm, contrast=contrasts, test.ZI=F, mc.cores=cores-1)
  # hangs on Phoenix with mc.cores
  fit.MDSeq.nozi.tmm <- MDSeq(counts, offsets=sf.tmm, contrast=contrasts, test.ZI=F)
  res.MDSeq.nozi.tmm <- extract.ZIMD(fit.MDSeq.nozi.tmm, compare=list(A="1",B="2"))
  p.disp.MDSeq.zi.tmm <- res.MDSeq.zi.tmm$Pvalue.disp
  q.disp.MDSeq.zi.tmm <- res.MDSeq.zi.tmm$FDR.disp
  p.disp.MDSeq.nozi.tmm <- res.MDSeq.nozi.tmm$Pvalue.disp
  q.disp.MDSeq.nozi.tmm <- res.MDSeq.nozi.tmm$FDR.disp
  rm(fit.MDSeq.zi.tmm, res.MDSeq.zi.tmm, fit.MDSeq.nozi.tmm, res.MDSeq.nozi.tmm)
  # fit.MDSeq.zi.rle <- MDSeq(counts, offsets=sf.rle, contrast=contrasts, mc.cores=cores-1)
  # hangs on Phoenix with mc.cores
  fit.MDSeq.zi.rle <- MDSeq(counts, offsets=sf.rle, contrast=contrasts)
  res.MDSeq.zi.rle <- extract.ZIMD(fit.MDSeq.zi.rle, compare=list(A="1",B="2"))
  # fit.MDSeq.nozi.rle <- MDSeq(counts, offsets=sf.rle, contrast=contrasts, test.ZI=F, mc.cores=cores-1)
  # hangs on Phoenix with mc.cores
  fit.MDSeq.nozi.rle <- MDSeq(counts, offsets=sf.rle, contrast=contrasts, test.ZI=F)
  res.MDSeq.nozi.rle <- extract.ZIMD(fit.MDSeq.nozi.rle, compare=list(A="1",B="2"))
  p.disp.MDSeq.zi.rle <- res.MDSeq.zi.rle$Pvalue.disp
  q.disp.MDSeq.zi.rle <- res.MDSeq.zi.rle$FDR.disp
  p.disp.MDSeq.nozi.rle <- res.MDSeq.nozi.rle$Pvalue.disp
  q.disp.MDSeq.nozi.rle <- res.MDSeq.nozi.rle$FDR.disp
  rm(fit.MDSeq.zi.rle, res.MDSeq.zi.rle, fit.MDSeq.nozi.rle, res.MDSeq.nozi.rle)
  
  ## diffVar
  fit.diffVar.tmm <- varFit(norm.counts.tmm, design=design, coef=c(1,2))
  res.diffVar.tmm <- topVar(fit.diffVar.tmm, coef=2, number=nrow(norm.counts.tmm), sort=F)
  p.diffVar.tmm <- res.diffVar.tmm$P.Value
  q.diffVar.tmm <- res.diffVar.tmm$Adj.P.Value
  rm(fit.diffVar.tmm, res.diffVar.tmm)
  fit.diffVar.rle <- varFit(norm.counts.rle, design=design, coef=c(1,2))
  res.diffVar.rle <- topVar(fit.diffVar.rle, coef=2, number=nrow(norm.counts.rle), sort=F)
  p.diffVar.rle <- res.diffVar.rle$P.Value
  q.diffVar.rle <- res.diffVar.rle$Adj.P.Value
  rm(fit.diffVar.rle, res.diffVar.rle)
  
  ## expHM
  expHM.tmm <- exp_hmm_adapt_3_chains(counts=t(norm.counts.tmm), groups=group)
  prob.expHMM.tmm <- unname(colMeans(as.matrix(expHM.tmm$indicators)))
  post.prop.expHMM.tmm <- mean(as.matrix(expHM.tmm$proportion))
  thr.expHMM.tmm <- sort(prob.expHMM.tmm, decreasing=T)[round(nrow(counts) * post.prop.expHMM.tmm)]
  bfdr.expHMM.tmm <- bfdr(prob.expHMM.tmm)
  disp.diff.expHM.untr.tmm <- unname(as.matrix(expHM.tmm$disps1) - as.matrix(expHM.tmm$disps2))
  p.disp.expHM.untr.tmm <- apply(disp.diff.expHM.untr.tmm,2,hpd.pval)
  rm(disp.diff.expHM.untr.tmm); gc()
  disp.diff.expHM.log.tmm <- unname(log(as.matrix(expHM.tmm$disps1)) - log(as.matrix(expHM.tmm$disps2)))
  rm(expHM.tmm); gc()
  p.disp.expHM.log.tmm <- apply(disp.diff.expHM.log.tmm,2,hpd.pval)
  rm(disp.diff.expHM.log.tmm); gc()
  q.disp.expHM.untr.tmm <- p.adjust(p.disp.expHM.untr.tmm, method='BH')
  q.disp.expHM.log.tmm <- p.adjust(p.disp.expHM.log.tmm, method='BH')
  expHM.rle <- exp_hmm_adapt_3_chains(counts=t(norm.counts.rle), groups=group)
  prob.expHMM.rle <- unname(colMeans(as.matrix(expHM.rle$indicators)))
  post.prop.expHMM.rle <- mean(as.matrix(expHM.rle$proportion))
  thr.expHMM.rle <- sort(prob.expHMM.rle, decreasing=T)[round(nrow(counts) * post.prop.expHMM.rle)]
  bfdr.expHMM.rle <- bfdr(prob.expHMM.rle)
  disp.diff.expHM.untr.rle <- unname(as.matrix(expHM.rle$disps1) - as.matrix(expHM.rle$disps2))
  p.disp.expHM.untr.rle <- apply(disp.diff.expHM.untr.rle,2,hpd.pval)
  rm(disp.diff.expHM.untr.rle); gc()
  disp.diff.expHM.log.rle <- unname(log(as.matrix(expHM.rle$disps1)) - log(as.matrix(expHM.rle$disps2)))
  rm(expHM.rle); gc()
  p.disp.expHM.log.rle <- apply(disp.diff.expHM.log.rle,2,hpd.pval)
  rm(disp.diff.expHM.log.rle); gc()
  q.disp.expHM.untr.rle <- p.adjust(p.disp.expHM.untr.rle, method='BH')
  q.disp.expHM.log.rle <- p.adjust(p.disp.expHM.log.rle, method='BH')
  
  ## lnHM
  lnHM.tmm <- ln_hmm_adapt_3_chains(counts=t(norm.counts.tmm), groups=group)
  prob.lnHMM.tmm <- unname(colMeans(as.matrix(lnHM.tmm$indicators)))
  post.prop.lnHMM.tmm <- mean(as.matrix(lnHM.tmm$proportion))
  thr.lnHMM.tmm <- sort(prob.lnHMM.tmm, decreasing=T)[round(nrow(counts) * post.prop.lnHMM.tmm)]
  bfdr.lnHMM.tmm <- bfdr(prob.lnHMM.tmm)
  disp.diff.lnHM.untr.tmm <- unname(as.matrix(lnHM.tmm$disps1) - as.matrix(lnHM.tmm$disps2))
  p.disp.lnHM.untr.tmm <- apply(disp.diff.lnHM.untr.tmm,2,hpd.pval)
  rm(disp.diff.lnHM.untr.tmm); gc()
  disp.diff.lnHM.log.tmm <- unname(log(as.matrix(lnHM.tmm$disps1)) - log(as.matrix(lnHM.tmm$disps2)))
  rm(lnHM.tmm); gc()
  p.disp.lnHM.log.tmm <- apply(disp.diff.lnHM.log.tmm,2,hpd.pval)
  rm(disp.diff.lnHM.log.tmm); gc()
  q.disp.lnHM.untr.tmm <- p.adjust(p.disp.lnHM.untr.tmm, method='BH')
  q.disp.lnHM.log.tmm <- p.adjust(p.disp.lnHM.log.tmm, method='BH')
  lnHM.rle <- ln_hmm_adapt_3_chains(counts=t(norm.counts.rle), groups=group)
  prob.lnHMM.rle <- unname(colMeans(as.matrix(lnHM.rle$indicators)))
  post.prop.lnHMM.rle <- mean(as.matrix(lnHM.rle$proportion))
  thr.lnHMM.rle <- sort(prob.lnHMM.rle, decreasing=T)[round(nrow(counts) * post.prop.lnHMM.rle)]
  bfdr.lnHMM.rle <- bfdr(prob.lnHMM.rle)
  disp.diff.lnHM.untr.rle <- unname(as.matrix(lnHM.rle$disps1) - as.matrix(lnHM.rle$disps2))
  p.disp.lnHM.untr.rle <- apply(disp.diff.lnHM.untr.rle,2,hpd.pval)
  rm(disp.diff.lnHM.untr.rle); gc()
  disp.diff.lnHM.log.rle <- unname(log(as.matrix(lnHM.rle$disps1)) - log(as.matrix(lnHM.rle$disps2)))
  rm(lnHM.rle); gc()
  p.disp.lnHM.log.rle <- apply(disp.diff.lnHM.log.rle,2,hpd.pval)
  rm(disp.diff.lnHM.log.rle); gc()
  q.disp.lnHM.untr.rle <- p.adjust(p.disp.lnHM.untr.rle, method='BH')
  q.disp.lnHM.log.rle <- p.adjust(p.disp.lnHM.log.rle, method='BH')
  
  results <- list(counts = counts, 
                  DD = DD, 
                  DEDD = DEDD, 
                  p.disp.MDSeq.zi.tmm = p.disp.MDSeq.zi.tmm, 
                  p.disp.MDSeq.nozi.tmm = p.disp.MDSeq.nozi.tmm, 
                  q.disp.MDSeq.zi.tmm = q.disp.MDSeq.zi.tmm, 
                  q.disp.MDSeq.nozi.tmm = q.disp.MDSeq.nozi.tmm, 
                  p.disp.MDSeq.zi.rle = p.disp.MDSeq.zi.rle, 
                  p.disp.MDSeq.nozi.rle = p.disp.MDSeq.nozi.rle, 
                  q.disp.MDSeq.zi.rle = q.disp.MDSeq.zi.rle, 
                  q.disp.MDSeq.nozi.rle = q.disp.MDSeq.nozi.rle, 
                  p.diffVar.tmm = p.diffVar.tmm, 
                  p.diffVar.rle = p.diffVar.rle, 
                  q.diffVar.tmm = q.diffVar.tmm, 
                  q.diffVar.rle = q.diffVar.rle, 
                  prob.expHMM.tmm = prob.expHMM.tmm, 
                  prop.expHMM.tmm = post.prop.expHMM.tmm, 
                  thr.expHMM.tmm = thr.expHMM.tmm, 
                  bfdr.expHMM.tmm = bfdr.expHMM.tmm, 
                  p.disp.expHM.untr.tmm = p.disp.expHM.untr.tmm, 
                  p.disp.expHM.log.tmm = p.disp.expHM.log.tmm, 
                  q.disp.expHM.untr.tmm = q.disp.expHM.untr.tmm, 
                  q.disp.expHM.log.tmm = q.disp.expHM.log.tmm, 
                  prob.expHMM.rle = prob.expHMM.rle, 
                  prop.expHMM.rle = post.prop.expHMM.rle, 
                  thr.expHMM.rle = thr.expHMM.rle, 
                  bfdr.expHMM.rle = bfdr.expHMM.rle, 
                  p.disp.expHM.untr.rle = p.disp.expHM.untr.rle, 
                  p.disp.expHM.log.rle = p.disp.expHM.log.rle, 
                  q.disp.expHM.untr.rle = q.disp.expHM.untr.rle, 
                  q.disp.expHM.log.rle = q.disp.expHM.log.rle, 
                  prob.lnHMM.tmm = prob.lnHMM.tmm, 
                  prop.lnHMM.tmm = post.prop.lnHMM.tmm, 
                  thr.lnHMM.tmm = thr.lnHMM.tmm, 
                  bfdr.lnHMM.tmm = bfdr.lnHMM.tmm, 
                  p.disp.lnHM.untr.tmm = p.disp.lnHM.untr.tmm, 
                  p.disp.lnHM.log.tmm = p.disp.lnHM.log.tmm, 
                  q.disp.lnHM.untr.tmm = q.disp.lnHM.untr.tmm, 
                  q.disp.lnHM.log.tmm = q.disp.lnHM.log.tmm, 
                  prob.lnHMM.rle = prob.lnHMM.rle, 
                  prop.lnHMM.rle = post.prop.lnHMM.rle, 
                  thr.lnHMM.rle = thr.lnHMM.rle, 
                  bfdr.lnHMM.rle = bfdr.lnHMM.rle, 
                  p.disp.lnHM.untr.rle = p.disp.lnHM.untr.rle, 
                  p.disp.lnHM.log.rle = p.disp.lnHM.log.rle, 
                  q.disp.lnHM.untr.rle = q.disp.lnHM.untr.rle, 
                  q.disp.lnHM.log.rle = q.disp.lnHM.log.rle)
  
  if (cluster) {
    folder <- paste0("Data/GTEx muscle simulated DE, DD, DEDD results March 2020")
  } else {
    folder <- paste0("Results/GTEx muscle simulated DE, DD, DEDD results March 2020")
  }
  filename <- paste0("results.", filename, ".rds")
  saveRDS(results, file=here(folder, filename))
  
  rm(counts, DD, DEDD, FC, libsizes, 
     nf.tmm, els.tmm, sf.tmm, norm.counts.tmm, 
     nf.rle, els.rle, sf.rle, norm.counts.rle)
  
}

filename <- paste0("sessionInfo.GTEx_muscle_", samples.per.cond, "_set", i, "_DD.rds")
saveRDS(sessionInfo, file=here(folder, filename))
