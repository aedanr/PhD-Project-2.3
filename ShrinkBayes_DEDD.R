library(here)
library(compcodeR)
library(edgeR)
library(ShrinkBayes)

samples.per.cond <- 2
group <- factor(c(rep(1,samples.per.cond), rep(2,samples.per.cond)))
cores <- detectCores()

for (i in 1:50) {
  # Generate filtered data
  filename <- paste0('DEDD', samples.per.cond, '.', i)
  counts <- readRDS(here('Simulated data',paste0(filename,'.rds')))
  DE <- counts@variable.annotations$differential.expression
  DD <- counts@variable.annotations$differential.dispersion
  lfcm1 <- abs(counts@variable.annotations$truelog2foldchanges) > 1
  lfcm2 <- abs(counts@variable.annotations$truelog2foldchanges) > 2
  lfcd1 <- abs(counts@variable.annotations$truelog2fcdispersion) > 1
  lfcd2 <- abs(counts@variable.annotations$truelog2fcdispersion) > 2
  
  # Normalise and create DGEList object
  nf <- calcNormFactors(counts@count.matrix)
  norm <- t(t(counts@count.matrix) / nf)

  ## ShrinkBayes
  form <- y ~ 1 + group
  subset <- sample(nrow(norm), 500)
  norm.subset <- norm[subset,]
  shrink <- ShrinkSeq(form=form, dat=norm.subset, fams='zinb', shrinkfixed='group', ncpus=cores)
  fit.subset.ShrinkBayes <- FitAllShrink(form, norm.subset, fams='zinb', shrinksimul=shrink, ncpus=cores)
  npprior.mix <- MixtureUpdatePrior(fitall=fit.subset.ShrinkBayes, shrinkpara='group', ncpus=cores)
  npprior.np <- NonParaUpdatePrior(fitall=fit.subset.ShrinkBayes, shrinkpara='group', ncpus=cores)
  npprior.np.lfc <- NonParaUpdatePrior(fitall=fit.subset.ShrinkBayes, shrinkpara='group', ncpus=cores, 
                                       includeP0=F)
  fit.ShrinkBayes <- FitAllShrink(form, norm, fams='zinb', shrinksimul=shrink, ncpus=cores)
  nppost.mix <- MixtureUpdatePosterior(fit.ShrinkBayes, npprior.mix, ncpus=cores)
  nppost.np <- NonParaUpdatePosterior(fit.ShrinkBayes, npprior.np, ncpus=cores)
  nppost.np.lfc <- NonParaUpdatePosterior(fit.ShrinkBayes, npprior.np.lfc, ncpus=cores)
  lfdr.mix.ShrinkBayes <- SummaryWrap(nppost.mix)
  lfdr.np.ShrinkBayes <- SummaryWrap(nppost.np)
  lfdr.np.lfc1.ShrinkBayes <- SummaryWrap(nppost.np.lfc, thr=1, direction='two-sided')
  lfdr.np.lfc2.ShrinkBayes <- SummaryWrap(nppost.np.lfc, thr=2, direction='two-sided')
  bfdr.mix.ShrinkBayes <- BFDR(lfdr.mix.ShrinkBayes)
  bfdr.np.ShrinkBayes <- BFDR(lfdr.np.ShrinkBayes)
  bfdr.np.lfc1.ShrinkBayes <- BFDR(lfdr.np.lfc1.ShrinkBayes)
  bfdr.np.lfc2.ShrinkBayes <- BFDR(lfdr.np.lfc2.ShrinkBayes)
  rm(list=c('form', 'subset', 'norm.subset', 'shrink', 'fit.subset.ShrinkBayes', 'npprior.mix', 
            'npprior.np', 'npprior.np.lfc', 'fit.ShrinkBayes', 'nppost.mix', 'nppost.np', 
            'nppost.np.lfc'))

  results <- list(data = counts, 
                  DE = DE, 
                  DD = DD, 
                  lfcm1 = lfcm1, 
                  lfcm2 = lfcm2, 
                  lfcd1 = lfcd1, 
                  lfcd2 = lfcd2, 
                  lfdr.mix.ShrinkBayes = lfdr.mix.ShrinkBayes, 
                  lfdr.np.ShrinkBayes = lfdr.np.ShrinkBayes, 
                  lfdr.np.lfc1.ShrinkBayes = lfdr.np.lfc1.ShrinkBayes, 
                  lfdr.np.lfc2.ShrinkBayes = lfdr.np.lfc2.ShrinkBayes, 
                  bfdr.mix.ShrinkBayes = bfdr.mix.ShrinkBayes, 
                  bfdr.np.ShrinkBayes = bfdr.np.ShrinkBayes, 
                  bfdr.np.lfc1.ShrinkBayes = bfdr.np.lfc1.ShrinkBayes, 
                  bfdr.np.lfc2.ShrinkBayes = bfdr.np.lfc2.ShrinkBayes)
  
  filename <- paste0('results.', filename, '.ShrinkBayes.rds')
  saveRDS(results, file=here(filename))

  rm(list=c('counts', 'DE', 'DD', 'lfcm1', 'lfcm2', 'lfcd1', 'lfcd2', 
            'nf', 'norm', 'filename', 'results'))
}

filename <- paste0('sessionInfo.DEDD', samples.per.cond, '.ShrinkBayes.rds')
saveRDS(sessionInfo(), file=here(filename))
