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
source(here('scripts','2019-04-03_exponential_hmm_adaptive_proposals_three_chains_function.R'))
source(here('scripts','2019-03-27_lognormal_hmm_adaptive_proposals_three_chains_function.R'))
source(here('scripts','2019-06-26_hpd_tail_prob_function.R'))
source(here('scripts','2019-05-17_compData_diff_disp_functions.R'))

samples.per.cond <- 2
group <- factor(c(rep(1,samples.per.cond), rep(2,samples.per.cond)))
design <- model.matrix(~group)
cores <- detectCores()
if(require("parallel")) cl <- makeCluster(cores) else cl <- NULL

for (i in 1:50) {
  # Generate filtered data
  #counts <- simulate.DE.DD.data(dataset='DEDD2', n.vars=20000, samples.per.cond=samples.per.cond)
  filename <- paste0('DEDD', samples.per.cond, '.', i)
  #counts <- readRDS(paste0(filename,'.rds'))
  counts <- readRDS(here('Data',paste0(filename,'.rds')))
  DE <- counts@variable.annotations$differential.expression
  DD <- counts@variable.annotations$differential.dispersion
  lfcm1 <- abs(counts@variable.annotations$truelog2foldchanges) > 1
  lfcm2 <- abs(counts@variable.annotations$truelog2foldchanges) > 2
  lfcd1 <- abs(counts@variable.annotations$truelog2fcdispersion) > 1
  lfcd2 <- abs(counts@variable.annotations$truelog2fcdispersion) > 2

  # Normalise and create DGEList object
  # nf <- calcNormFactors(counts@count.matrix)
  dat.DESeq <- DESeqDataSetFromMatrix(countData=counts@count.matrix, colData=data.frame(group), 
                                      design=~group)
  dat.DESeq <- estimateSizeFactors(dat.DESeq)
  nf <- dat.DESeq$sizeFactor
  norm <- t(t(counts@count.matrix) / nf)
  dge <- DGEList(counts=counts@count.matrix, norm.factors=nf, group=group)
  
  ## edgeR
  dat.edgeR <- estimateDisp(dge, design)
  qlfit.edgeR <- glmQLFit(dat.edgeR, design)
  qltest.edgeR <- glmQLFTest(qlfit.edgeR)
  ql.lfc1.edgeR <- glmTreat(qlfit.edgeR, lfc=1)
  ql.lfc2.edgeR <- glmTreat(qlfit.edgeR, lfc=2)
  lrfit.edgeR <- glmFit(dat.edgeR, design)
  lrtest.edgeR <- glmLRT(lrfit.edgeR)
  lr.lfc1.edgeR <- glmTreat(lrfit.edgeR, lfc=1)
  lr.lfc2.edgeR <- glmTreat(lrfit.edgeR, lfc=2)
  et.edgeR <- exactTest(dat.edgeR)
  p.ql.edgeR <- qltest.edgeR$table$PValue
  p.ql.lfc1.edgeR <- ql.lfc1.edgeR$table$PValue
  p.ql.lfc2.edgeR <- ql.lfc2.edgeR$table$PValue
  p.lr.edgeR <- lrtest.edgeR$table$PValue
  p.lr.lfc1.edgeR <- lr.lfc1.edgeR$table$PValue
  p.lr.lfc2.edgeR <- lr.lfc2.edgeR$table$PValue
  p.et.edgeR <- et.edgeR$table$PValue
  q.ql.edgeR <- topTags(qltest.edgeR, n=nrow(dat.edgeR$counts), sort='none')$table$FDR
  q.ql.lfc1.edgeR <- topTags(ql.lfc1.edgeR, n=nrow(dat.edgeR$counts), sort='none')$table$FDR
  q.ql.lfc2.edgeR <- topTags(ql.lfc2.edgeR, n=nrow(dat.edgeR$counts), sort='none')$table$FDR
  q.lr.edgeR <- topTags(lrtest.edgeR, n=nrow(dat.edgeR$counts), sort='none')$table$FDR
  q.lr.lfc1.edgeR <- topTags(lr.lfc1.edgeR, n=nrow(dat.edgeR$counts), sort='none')$table$FDR
  q.lr.lfc2.edgeR <- topTags(lr.lfc2.edgeR, n=nrow(dat.edgeR$counts), sort='none')$table$FDR
  q.et.edgeR <- topTags(et.edgeR, n=nrow(dat.edgeR$counts), sort='none')$table$FDR
  rm(list=c('dat.edgeR', 'qlfit.edgeR', 'qltest.edgeR', 'ql.lfc1.edgeR', 'ql.lfc2.edgeR', 'lrfit.edgeR', 
            'lrtest.edgeR', 'lr.lfc1.edgeR', 'lr.lfc2.edgeR', 'et.edgeR'))
  
  ## DESeq2
  dat.DESeq <- DESeq(dat.DESeq, minReplicatesForReplace=Inf)
  res.noif.DESeq <- results(dat.DESeq, independentFiltering=F, cooksCutoff=F)
  res.if.DESeq <- results(dat.DESeq, cooksCutoff=F, alpha=0.05)
  res.lfc1.DESeq <- results(dat.DESeq, lfcThreshold=1, independentFiltering=F, cooksCutoff=F)
  res.lfc2.DESeq <- results(dat.DESeq, lfcThreshold=2, independentFiltering=F, cooksCutoff=F)
  p.noif.DESeq <- res.noif.DESeq$pvalue
  p.if.DESeq <- res.noif.DESeq$pvalue
  p.lfc1.DESeq <- res.lfc1.DESeq$pvalue
  p.lfc2.DESeq <- res.lfc2.DESeq$pvalue
  q.noif.DESeq <- res.noif.DESeq$padj
  q.if.DESeq <- res.noif.DESeq$padj
  q.lfc1.DESeq <- res.lfc1.DESeq$padj
  q.lfc2.DESeq <- res.lfc2.DESeq$padj
  rm(list=c('dat.DESeq', 'res.noif.DESeq', 'res.if.DESeq', 'res.lfc1.DESeq', 
            'res.lfc2.DESeq'))
  
  ## limma-voom
  dat.voom <- voom(dge)
  fit.voom <- lmFit(dat.voom, design)
  res.voom <- eBayes(fit.voom)
  lfc1.voom <- treat(fit.voom, lfc=1)
  lfc2.voom <- treat(fit.voom, lfc=2)
  p.voom <- topTable(res.voom, number=Inf, sort.by='none')$P.Value
  p.lfc1.voom <- topTreat(lfc1.voom, number=Inf, sort.by='none', coef=2)$P.Value
  p.lfc2.voom <- topTreat(lfc2.voom, number=Inf, sort.by='none', coef=2)$P.Value
  q.voom <- topTable(res.voom, number=Inf, sort.by='none')$adj.P.Val
  q.lfc1.voom <- topTreat(lfc1.voom, number=Inf, sort.by='none', coef=2)$adj.P.Val
  q.lfc2.voom <- topTreat(lfc2.voom, number=Inf, sort.by='none', coef=2)$adj.P.Val
  rm(list=c('dat.voom', 'fit.voom', 'res.voom', 'lfc1.voom', 'lfc2.voom'))
  
  ## DSS
  dat.DSS <- newSeqCountSet(counts=matrix(counts@count.matrix, ncol=length(group)), 
                            designs=as.numeric(group), normalizationFactor=nf)
  notrend.DSS <- estDispersion(dat.DSS)
  trend.DSS <- estDispersion(dat.DSS, trend=T)
  res.notrend.DSS <- waldTest(notrend.DSS, 1, 2)[order(waldTest(notrend.DSS, 1, 2)$geneIndex),]
  res.trend.DSS <- waldTest(trend.DSS, 1, 2)[order(waldTest(trend.DSS, 1, 2)$geneIndex),]
  p.notrend.DSS <- res.notrend.DSS$pval
  q.notrend.DSS <- res.notrend.DSS$fdr
  lfdr.notrend.DSS <- res.notrend.DSS$local.fdr
  p.trend.DSS <- res.trend.DSS$pval
  q.trend.DSS <- res.trend.DSS$fdr
  lfdr.trend.DSS <- res.trend.DSS$local.fdr
  rm(list=c('dat.DSS', 'notrend.DSS', 'trend.DSS', 'res.notrend.DSS', 'res.trend.DSS'))
  
  ## baySeq
  dat.baySeq <- new('countData', data=counts@count.matrix, replicates=group, 
                    groups=list(NDE=rep(1,length(group)), DE=as.numeric(group)))
  dat.baySeq@annotation <- data.frame(name = 1:nrow(dat.baySeq@data))
  libsizes(dat.baySeq) <- nf
  dat.baySeq <- getPriors.NB(dat.baySeq, samplesize=10000, cl=cl)
  dat.baySeq <- getLikelihoods(dat.baySeq, cl=cl)
  res.baySeq <- topCounts(dat.baySeq, group='DE', 
                          number=Inf)[order(topCounts(dat.baySeq, group='DE', number=Inf)$annotation),]
  # $annotation on cluster (presumably different in older version), $name on my laptop
  prob.baySeq <- res.baySeq$Likelihood
  # $Likelihood on cluster (presumably different in older version), $likes on my laptop
  q.baySeq <- res.baySeq$FDR.DE
  rm(list=c('dat.baySeq', 'res.baySeq'))
  
  ## ShrinkBayes
#  library(ShrinkBayes)
#  form <- y ~ 1 + group
#  subset <- sample(nrow(norm), 500)
#  norm.subset <- norm[subset,]
#  shrink <- ShrinkSeq(form=form, dat=norm.subset, fams='zinb', shrinkfixed='group', ncpus=cores)
#  fit.subset.ShrinkBayes <- FitAllShrink(form, norm.subset, fams='zinb', shrinksimul=shrink, ncpus=cores)
#  npprior.mix <- MixtureUpdatePrior(fitall=fit.subset.ShrinkBayes, shrinkpara='group', ncpus=cores)
#  npprior.np <- NonParaUpdatePrior(fitall=fit.subset.ShrinkBayes, shrinkpara='group', ncpus=cores)
#  npprior.np.lfc <- NonParaUpdatePrior(fitall=fit.subset.ShrinkBayes, shrinkpara='group', ncpus=cores, 
#                                       includeP0=F)
#  fit.ShrinkBayes <- FitAllShrink(form, norm, fams='zinb', shrinksimul=shrink, ncpus=cores)
#  nppost.mix <- MixtureUpdatePosterior(fit.ShrinkBayes, npprior.mix, ncpus=cores)
#  nppost.np <- NonParaUpdatePosterior(fit.ShrinkBayes, npprior.np, ncpus=cores)
#  nppost.np.lfc <- NonParaUpdatePosterior(fit.ShrinkBayes, npprior.np.lfc, ncpus=cores)
#  lfdr.mix.ShrinkBayes <- SummaryWrap(nppost.mix)
#  lfdr.np.ShrinkBayes <- SummaryWrap(nppost.np)
#  lfdr.np.lfc1.ShrinkBayes <- SummaryWrap(nppost.np.lfc, thr=1, direction='two-sided')
#  lfdr.np.lfc2.ShrinkBayes <- SummaryWrap(nppost.np.lfc, thr=2, direction='two-sided')
#  bfdr.mix.ShrinkBayes <- BFDR(lfdr.mix.ShrinkBayes)
#  bfdr.np.ShrinkBayes <- BFDR(lfdr.np.ShrinkBayes)
#  bfdr.np.lfc1.ShrinkBayes <- BFDR(lfdr.np.lfc1.ShrinkBayes)
#  bfdr.np.lfc2.ShrinkBayes <- BFDR(lfdr.np.lfc2.ShrinkBayes)
#  rm(list=c('form', 'subset', 'norm.subset', 'shrink', 'fit.subset.ShrinkBayes', 'npprior.mix', 
#            'npprior.np', 'npprior.np.lfc', 'fit.ShrinkBayes', 'nppost.mix', 'nppost.np', 
#            'nppost.np.lfc'))
#  detach("package:ShrinkBayes", unload=TRUE)
#  detach("package:snowfall", unload=TRUE)
#  detach("package:snow", unload=TRUE)
  
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
  p.mean.zi.MDSeq <- res.zi.MDSeq$Pvalue.mean
  p.mean.zi.lfc1.MDSeq <- res.zi.lfc1.MDSeq$Pvalue.mean
  p.mean.zi.lfc2.MDSeq <- res.zi.lfc2.MDSeq$Pvalue.mean
  p.mean.nozi.MDSeq <- res.nozi.MDSeq$Pvalue.mean
  p.mean.nozi.lfc1.MDSeq <- res.nozi.lfc1.MDSeq$Pvalue.mean
  p.mean.nozi.lfc2.MDSeq <- res.nozi.lfc2.MDSeq$Pvalue.mean
  p.disp.zi.MDSeq <- res.zi.MDSeq$Pvalue.dispersion
  p.disp.zi.lfc1.MDSeq <- res.zi.lfc1.MDSeq$Pvalue.dispersion
  p.disp.zi.lfc2.MDSeq <- res.zi.lfc2.MDSeq$Pvalue.dispersion
  p.disp.nozi.MDSeq <- res.nozi.MDSeq$Pvalue.dispersion
  p.disp.nozi.lfc1.MDSeq <- res.nozi.lfc1.MDSeq$Pvalue.dispersion
  p.disp.nozi.lfc2.MDSeq <- res.nozi.lfc2.MDSeq$Pvalue.dispersion
  q.mean.zi.MDSeq <- res.zi.MDSeq$FDR.mean
  q.mean.zi.lfc1.MDSeq <- res.zi.lfc1.MDSeq$FDR.mean
  q.mean.zi.lfc2.MDSeq <- res.zi.lfc2.MDSeq$FDR.mean
  q.mean.nozi.MDSeq <- res.nozi.MDSeq$FDR.mean
  q.mean.nozi.lfc1.MDSeq <- res.nozi.lfc1.MDSeq$FDR.mean
  q.mean.nozi.lfc2.MDSeq <- res.nozi.lfc2.MDSeq$FDR.mean
  q.disp.zi.MDSeq <- res.zi.MDSeq$FDR.dispersion
  q.disp.zi.lfc1.MDSeq <- res.zi.lfc1.MDSeq$FDR.dispersion
  q.disp.zi.lfc2.MDSeq <- res.zi.lfc2.MDSeq$FDR.dispersion
  q.disp.nozi.MDSeq <- res.nozi.MDSeq$FDR.dispersion
  q.disp.nozi.lfc1.MDSeq <- res.nozi.lfc1.MDSeq$FDR.dispersion
  q.disp.nozi.lfc2.MDSeq <- res.nozi.lfc2.MDSeq$FDR.dispersion
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
  p.mean.expHM <- apply(mean.diff.expHM,2,hpd.pval, m=0)
  p.lmean.expHM <- apply(log.mean.diff.expHM,2,hpd.pval, m=0)
  p.mean.lfc1.expHM <- apply(log.mean.diff.expHM,2,hpd.pval, m=1)
  p.mean.lfc2.expHM <- apply(log.mean.diff.expHM,2,hpd.pval, m=2)
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
  p.mean.lnHM <- apply(mean.diff.lnHM,2,hpd.pval)
  p.lmean.lnHM <- apply(log.mean.diff.lnHM,2,hpd.pval)
  p.mean.lfc1.lnHM <- apply(log.mean.diff.lnHM,2,hpd.pval, m=1)
  p.mean.lfc2.lnHM <- apply(log.mean.diff.lnHM,2,hpd.pval, m=2)
  p.disp.lnHM <- apply(disp.diff.lnHM,2,hpd.pval)
  p.ldisp.lnHM <- apply(log.disp.diff.lnHM,2,hpd.pval)
  p.disp.lfc1.lnHM <- apply(log.disp.diff.lnHM,2,hpd.pval, m=1)
  p.disp.lfc2.lnHM <- apply(log.disp.diff.lnHM,2,hpd.pval, m=2)
  rm(list=c('lnHM', 'mean.diff.lnHM', 'log.mean.diff.lnHM', 
            'disp.diff.lnHM', 'log.disp.diff.lnHM'))
  
  results <- list(data = counts, 
                  DE = DE, 
                  DD = DD, 
                  lfcm1 = lfcm1, 
                  lfcm2 = lfcm2, 
                  lfcd1 = lfcd1, 
                  lfcd2 = lfcd2, 
                  p.ql.edgeR = p.ql.edgeR, 
                  p.ql.lfc1.edgeR = p.ql.lfc1.edgeR, 
                  p.ql.lfc2.edgeR = p.ql.lfc2.edgeR, 
                  p.lr.edgeR = p.lr.edgeR, 
                  p.lr.lfc1.edgeR = p.lr.lfc1.edgeR, 
                  p.lr.lfc2.edgeR = p.lr.lfc2.edgeR, 
                  p.et.edgeR = p.et.edgeR, 
                  q.ql.edgeR = q.ql.edgeR, 
                  q.ql.lfc1.edgeR = q.ql.lfc1.edgeR, 
                  q.ql.lfc2.edgeR = q.ql.lfc2.edgeR, 
                  q.lr.edgeR = q.lr.edgeR, 
                  q.lr.lfc1.edgeR = q.lr.lfc1.edgeR, 
                  q.lr.lfc2.edgeR = q.lr.lfc2.edgeR, 
                  q.et.edgeR = q.et.edgeR, 
                  p.noif.DESeq = p.noif.DESeq, 
                  p.if.DESeq = p.if.DESeq, 
                  p.lfc1.DESeq = p.lfc1.DESeq, 
                  p.lfc2.DESeq = p.lfc2.DESeq, 
                  q.noif.DESeq = q.noif.DESeq, 
                  q.if.DESeq = q.if.DESeq, 
                  q.lfc1.DESeq = q.lfc1.DESeq, 
                  q.lfc2.DESeq = q.lfc2.DESeq, 
                  p.voom = p.voom, 
                  p.lfc1.voom = p.lfc1.voom, 
                  p.lfc2.voom = p.lfc2.voom, 
                  q.voom = q.voom, 
                  q.lfc1.voom = q.lfc1.voom, 
                  q.lfc2.voom = q.lfc2.voom, 
                  p.notrend.DSS = p.notrend.DSS, 
                  q.notrend.DSS = q.notrend.DSS, 
                  lfdr.notrend.DSS = lfdr.notrend.DSS, 
                  p.trend.DSS = p.trend.DSS, 
                  q.trend.DSS = q.trend.DSS, 
                  lfdr.trend.DSS = lfdr.trend.DSS, 
                  prob.baySeq = prob.baySeq, 
                  q.baySeq = q.baySeq, 
                  #lfdr.mix.ShrinkBayes = lfdr.mix.ShrinkBayes, 
                  #lfdr.np.ShrinkBayes = lfdr.np.ShrinkBayes, 
                  #lfdr.np.lfc1.ShrinkBayes = lfdr.np.lfc1.ShrinkBayes, 
                  #lfdr.np.lfc2.ShrinkBayes = lfdr.np.lfc2.ShrinkBayes, 
                  #bfdr.mix.ShrinkBayes = bfdr.mix.ShrinkBayes, 
                  #bfdr.np.ShrinkBayes = bfdr.np.ShrinkBayes, 
                  #bfdr.np.lfc1.ShrinkBayes = bfdr.np.lfc1.ShrinkBayes, 
                  #bfdr.np.lfc2.ShrinkBayes = bfdr.np.lfc2.ShrinkBayes, 
                  p.mean.zi.MDSeq = p.mean.zi.MDSeq, 
                  p.mean.zi.lfc1.MDSeq = p.mean.zi.lfc1.MDSeq, 
                  p.mean.zi.lfc2.MDSeq = p.mean.zi.lfc2.MDSeq, 
                  p.mean.nozi.MDSeq = p.mean.nozi.MDSeq, 
                  p.mean.nozi.lfc1.MDSeq = p.mean.nozi.lfc1.MDSeq, 
                  p.mean.nozi.lfc2.MDSeq = p.mean.nozi.lfc2.MDSeq, 
                  p.disp.zi.MDSeq = p.disp.zi.MDSeq, 
                  p.disp.zi.lfc1.MDSeq = p.disp.zi.lfc1.MDSeq, 
                  p.disp.zi.lfc2.MDSeq = p.disp.zi.lfc2.MDSeq, 
                  p.disp.nozi.MDSeq = p.disp.nozi.MDSeq, 
                  p.disp.nozi.lfc1.MDSeq = p.disp.nozi.lfc1.MDSeq, 
                  p.disp.nozi.lfc2.MDSeq = p.disp.nozi.lfc2.MDSeq, 
                  q.mean.zi.MDSeq = q.mean.zi.MDSeq, 
                  q.mean.zi.lfc1.MDSeq = q.mean.zi.lfc1.MDSeq, 
                  q.mean.zi.lfc2.MDSeq = q.mean.zi.lfc2.MDSeq, 
                  q.mean.nozi.MDSeq = q.mean.nozi.MDSeq, 
                  q.mean.nozi.lfc1.MDSeq = q.mean.nozi.lfc1.MDSeq, 
                  q.mean.nozi.lfc2.MDSeq = q.mean.nozi.lfc2.MDSeq, 
                  q.disp.zi.MDSeq = q.disp.zi.MDSeq, 
                  q.disp.zi.lfc1.MDSeq = q.disp.zi.lfc1.MDSeq, 
                  q.disp.zi.lfc2.MDSeq = q.disp.zi.lfc2.MDSeq, 
                  q.disp.nozi.MDSeq = q.disp.nozi.MDSeq, 
                  q.disp.nozi.lfc1.MDSeq = q.disp.nozi.lfc1.MDSeq, 
                  q.disp.nozi.lfc2.MDSeq = q.disp.nozi.lfc2.MDSeq, 
                  prob.expHM = prob.expHM, 
                  prop.expHM = post.prop.expHM, 
                  p.mean.expHM = p.mean.expHM, 
                  p.lmean.expHM = p.lmean.expHM, 
                  p.mean.lfc1.expHM = p.mean.lfc1.expHM, 
                  p.mean.lfc2.expHM = p.mean.lfc2.expHM, 
                  p.disp.expHM = p.disp.expHM, 
                  p.ldisp.expHM = p.ldisp.expHM, 
                  p.disp.lfc1.expHM = p.disp.lfc1.expHM, 
                  p.disp.lfc2.expHM = p.disp.lfc2.expHM, 
                  prob.lnHM = prob.lnHM, 
                  prop.lnHM = post.prop.lnHM, 
                  p.mean.lnHM = p.mean.lnHM, 
                  p.lmean.lnHM = p.lmean.lnHM, 
                  p.mean.lfc1.lnHM = p.mean.lfc1.lnHM, 
                  p.mean.lfc2.lnHM = p.mean.lfc2.lnHM, 
                  p.disp.lnHM = p.disp.lnHM, 
                  p.ldisp.lnHM = p.ldisp.lnHM, 
                  p.disp.lfc1.lnHM = p.disp.lfc1.lnHM, 
                  p.disp.lfc2.lnHM = p.disp.lfc2.lnHM)
  
  filename <- paste0('results.', filename, '.DESeqnorm.rds')
  saveRDS(results, file=here(filename))
  
  rm(list=c('counts', 'DE', 'DD', 'lfcm1', 'lfcm2', 'lfcd1', 'lfcd2', 
            'nf', 'norm', 'dge', 'filename', 'results'))
}

filename <- paste0('sessionInfo.DEDD', samples.per.cond, '.DESeqnorm.rds')
saveRDS(sessionInfo(), file=here(filename))


