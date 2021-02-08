library(here)
library(compcodeR)
library(limma)
library(edgeR)
library(missMethyl)

samples.per.cond <- 50
group <- factor(c(rep(1,samples.per.cond), rep(2,samples.per.cond)))
design <- model.matrix(~group)

for (i in 1:50) {
  # Import filtered data
  filename <- paste0('DD', samples.per.cond, '.', i)
  # data <- readRDS(here('Data/Simulated data', paste0(filename,'.rds')))
  data <- readRDS(here('Simulated data', paste0(filename,'.rds')))
  counts <- data@count.matrix
  DD <- data@variable.annotations$differential.dispersion
  DEDD <- DD
  lfcm1 <- abs(data@variable.annotations$truelog2foldchanges) > 1
  lfcm2 <- abs(data@variable.annotations$truelog2foldchanges) > 2
  lfcd1 <- abs(data@variable.annotations$truelog2fcdispersion) > 1
  lfcd2 <- abs(data@variable.annotations$truelog2fcdispersion) > 2
  
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
  
  results <- list(data = data, 
                  DD = DD, 
                  DEDD = DEDD, 
                  lfcm1 = lfcm1, 
                  lfcm2 = lfcm2, 
                  lfcd1 = lfcd1, 
                  lfcd2 = lfcd2, 
                  p.diffVar.tmm = p.diffVar.tmm, 
                  p.diffVar.rle = p.diffVar.rle, 
                  q.diffVar.tmm = q.diffVar.tmm, 
                  q.diffVar.rle = q.diffVar.rle)
  
  filename <- paste0('results.', filename, '.rds')
  # saveRDS(results, file=here('Data', filename))
  saveRDS(results, file=here('Results/diffVar compcodeR results Dec 2020', filename))
  
  rm(data, counts, DD, DEDD, lfcm1, lfcm2, lfcd1, lfcd2, libsizes, 
     nf.tmm, els.tmm, sf.tmm, norm.counts.tmm, nf.rle, els.rle, sf.rle, norm.counts.rle)
  
}

filename <- paste0('sessionInfo.DD', samples.per.cond, '.rds')
# saveRDS(sessionInfo(), file=here('Data', filename))
saveRDS(sessionInfo(), file=here('Results/diffVar compcodeR results Dec 2020', filename))

