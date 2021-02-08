# MDSeq ZI only, lnHM log only since just want lfc1 results as quickly as possible for paper and thesis.
# Include results without lfc threshold to allow comparison with previous results for QC.

library(here)
library(coda)
library(HDInterval)
library(compcodeR)
library(edgeR)
library(MDSeq)

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

for (i in 1:10) {
  # Data
  if (cluster) {
    folder <- paste0("Data/recount data/GTEx/blood_", samples.per.cond, "_samples_per_group")
  } else {
    folder <- paste0("recount data/GTEx/blood_", samples.per.cond, "_samples_per_group")
  }
  filename <- paste0("blood_", samples.per.cond, "_set", i, "_DEDD")
  data <- readRDS(here(folder, paste0(filename, ".rds")))
  counts <- data$counts
  DD <- data$DDindex
  DE <- data$DEindex
  DEDD <- as.numeric(DE == 1 | DD == 1)
  FC.mean <- data$FC.mean
  FC.disp <- data$FC.disp
  lfcd1 <- abs(log2(data$FC.disp)) > 1
  rm(data)

  # Normalise
  libsizes <- colSums(counts)
  nf <- calcNormFactors(counts, method="TMM")
  els <- nf * libsizes
  sf <- els / exp(mean(log(libsizes)))
  norm.counts <- t(t(counts) / sf)

  ## MDSeq
  fit.MDSeq.zi <- MDSeq(counts, offsets=sf, contrast=contrasts)
  res.MDSeq.zi <- extract.ZIMD(fit.MDSeq.zi, compare=list(A="1",B="2"))
  lfc1.MDSeq.zi <- extract.ZIMD(fit.MDSeq.zi, compare=list(A="1",B="2"), log2FC.threshold=1)
  p.disp.MDSeq.zi <- res.MDSeq.zi$Pvalue.disp
  p.disp.lfc1.MDSeq.zi <- lfc1.MDSeq.zi$Pvalue.disp
  q.disp.MDSeq.zi <- res.MDSeq.zi$FDR.disp
  q.disp.lfc1.MDSeq.zi <- lfc1.MDSeq.zi$FDR.disp
  rm(fit.MDSeq.zi, res.MDSeq.zi, lfc1.MDSeq.zi)

  ## lnHM
  lnHM <- ln_hmm_adapt_3_chains(counts=t(norm.counts), groups=group)
  prob.lnHMM <- unname(colMeans(as.matrix(lnHM$indicators)))
  post.prop.lnHMM <- mean(as.matrix(lnHM$proportion))
  thr.lnHMM <- sort(prob.lnHMM, decreasing=T)[round(nrow(counts) * post.prop.lnHMM)]
  bfdr.lnHMM <- bfdr(prob.lnHMM)
  disp.diff.lnHM <- unname(log(as.matrix(lnHM$disps1)) - log(as.matrix(lnHM$disps2)))
  rm(lnHM); gc()
  p.disp.lnHM <- apply(disp.diff.lnHM,2,hpd.pval)
  p.disp.lfc1.lnHM <- apply(disp.diff.lnHM,2,hpd.pval, m=1)
  rm(disp.diff.lnHM); gc()
  q.disp.lnHM <- p.adjust(p.disp.lnHM, method='BH')
  q.disp.lfc1.lnHM <- p.adjust(p.disp.lfc1.lnHM, method='BH')
  
  results <- list(counts = counts, 
                  DE = DE, 
                  DD = DD, 
                  DEDD = DEDD, 
                  FC.mean = FC.mean, 
                  FC.disp = FC.disp, 
                  lfcd1 = lfcd1, 
                  p.disp.MDSeq.zi = p.disp.MDSeq.zi, 
                  p.disp.lfc1.MDSeq.zi = p.disp.lfc1.MDSeq.zi, 
                  q.disp.MDSeq.zi = q.disp.MDSeq.zi, 
                  q.disp.lfc1.MDSeq.zi = q.disp.lfc1.MDSeq.zi, 
                  prob.lnHMM = prob.lnHMM, 
                  prop.lnHMM = post.prop.lnHMM, 
                  thr.lnHMM = thr.lnHMM, 
                  bfdr.lnHMM = bfdr.lnHMM, 
                  p.disp.lnHM = p.disp.lnHM, 
                  p.disp.lfc1.lnHM = p.disp.lfc1.lnHM, 
                  q.disp.lnHM = q.disp.lnHM, 
                  q.disp.lfc1.lnHM = q.disp.lfc1.lnHM)
  
  if (cluster) {
    folder <- paste0("Data/GTEx DD lfc1 results Jan 2021")
  } else {
    folder <- paste0("Results/GTEx DD lfc1 results Jan 2021")
  }
  filename <- paste0("results.", filename, ".rds")
  saveRDS(results, file=here(folder, filename))
  
  rm(counts, DE, DD, DEDD, FC.mean, FC.disp, libsizes, nf, els, sf, norm.counts)
  
}

filename <- paste0("sessionInfo.GTEx_blood_", samples.per.cond, "_set", i, "_DD.lfc1.rds")
saveRDS(sessionInfo, file=here(folder, filename))
