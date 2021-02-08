library(here)
library(coda)
library(HDInterval)
library(edgeR)
library(limma)
library(DESeq2)
library(DSS)
library(baySeq)
library(MDSeq)
source(here('scripts','2019-04-03_exponential_hmm_adaptive_proposals_three_chains_function.R'))
source(here('scripts','2019-03-27_lognormal_hmm_adaptive_proposals_three_chains_function.R'))
source(here('scripts','2019-06-26_hpd_tail_prob_function.R'))

cores <- detectCores() - 1
if(require("parallel")) cl <- makeCluster(cores) else cl <- NULL

# DE ####
for (i in c(4, 10, 20, 100)) {
  # Data
  filename <- paste0("blood_", i, "_DE")
  data <- readRDS(here("recount data/GTEx", paste0(filename, ".rds")))
  counts <- data$counts
  DE <- data$DEindex
  FC <- data$FC
  group <- factor(c(rep(1,i/2), rep(2,i/2)))
  design <- model.matrix(~group)
  
  # Normalisation
  nf.TMM <- calcNormFactors(counts)
  norm.TMM <- t(t(counts) / nf.TMM)
  dat.edgeR <- DGEList(counts=counts, 
                       norm.factors=nf.TMM, 
                       group=group)
  dat.DESeq2 <- DESeqDataSetFromMatrix(countData=counts, 
                                       colData=data.frame(group), 
                                       design=~group)
  dat.DESeq2 <- estimateSizeFactors(dat.DESeq2)
  nf.DESeq2 <- dat.DESeq2$sizeFactor
  norm.DESeq2 <- t(t(counts) / nf.DESeq2)
  
  # edgeR
  dat.edgeR <- estimateDisp(dat.edgeR, design)
  qlfit.edgeR <- glmQLFit(dat.edgeR, design)
  qltest.edgeR <- glmQLFTest(qlfit.edgeR)
  lrfit.edgeR <- glmFit(dat.edgeR, design)
  lrtest.edgeR <- glmLRT(lrfit.edgeR)
  et.edgeR <- exactTest(dat.edgeR)
  p.ql.edgeR <- qltest.edgeR$table$PValue
  p.lr.edgeR <- lrtest.edgeR$table$PValue
  p.et.edgeR <- et.edgeR$table$PValue
  q.ql.edgeR <- topTags(qltest.edgeR, n=nrow(dat.edgeR$counts), sort='none')$table$FDR
  q.lr.edgeR <- topTags(lrtest.edgeR, n=nrow(dat.edgeR$counts), sort='none')$table$FDR
  q.et.edgeR <- topTags(et.edgeR, n=nrow(dat.edgeR$counts), sort='none')$table$FDR
  rm(qlfit.edgeR, qltest.edgeR, lrfit.edgeR, lrtest.edgeR, et.edgeR)
  
  # DESeq2
  dat.DESeq2 <- DESeq(dat.DESeq2, minReplicatesForReplace=Inf)
  res.if.DESeq2 <- results(dat.DESeq2, cooksCutoff=F, alpha=0.05)
  res.noif.DESeq2 <- results(dat.DESeq2, independentFiltering=F, cooksCutoff=F)
  p.noif.DESeq2 <- res.noif.DESeq2$pvalue
  p.if.DESeq2 <- res.if.DESeq2$pvalue
  q.noif.DESeq2 <- res.noif.DESeq2$padj
  q.if.DESeq2 <- res.if.DESeq2$padj
  rm(list=c('dat.DESeq2', 'res.noif.DESeq2', 'res.if.DESeq2'))
  
  # limma-voom
  dat.voom <- voom(dat.edgeR)
  fit.voom <- lmFit(dat.voom, design)
  res.voom <- eBayes(fit.voom)
  p.voom <- topTable(res.voom, number=Inf, sort.by='none')$P.Value
  q.voom <- topTable(res.voom, number=Inf, sort.by='none')$adj.P.Val
  rm(list=c('dat.edgeR', 'dat.voom', 'fit.voom', 'res.voom'))
  
  # DSS
  dat.DSS <- newSeqCountSet(counts=matrix(counts, ncol=length(group)), 
                            designs=as.numeric(group), 
                            normalizationFactor=nf.DESeq2)
  notrend.DSS <- estDispersion(dat.DSS)
  # trend.DSS <- estDispersion(dat.DSS, trend=T)
  res.notrend.DSS <- waldTest(notrend.DSS, 1, 2)[order(waldTest(notrend.DSS, 1, 2)$geneIndex),]
  # res.trend.DSS <- waldTest(trend.DSS, 1, 2)[order(waldTest(trend.DSS, 1, 2)$geneIndex),]
  p.notrend.DSS <- res.notrend.DSS$pval
  q.notrend.DSS <- res.notrend.DSS$fdr
  lfdr.notrend.DSS <- res.notrend.DSS$local.fdr
  # p.trend.DSS <- res.trend.DSS$pval
  # q.trend.DSS <- res.trend.DSS$fdr
  # lfdr.trend.DSS <- res.trend.DSS$local.fdr
  # rm(list=c('dat.DSS', 'notrend.DSS', 'trend.DSS', 'res.notrend.DSS', 'res.trend.DSS'))
  rm(list=c('dat.DSS', 'notrend.DSS', 'res.notrend.DSS'))
  
  # baySeq
  dat.baySeq <- new("countData",
                    data=counts, 
                    replicates=group, 
                    groups=list(NDE=rep(1, length(group)), DE=as.numeric(group)))
  dat.baySeq@annotation <- data.frame(name = 1:nrow(dat.baySeq@data))
  libsizes(dat.baySeq) <- nf.DESeq2
  dat.baySeq <- getPriors.NB(dat.baySeq, samplesize=10000, cl=cl)
  dat.baySeq <- getLikelihoods(dat.baySeq, cl=cl)
  res.baySeq <- topCounts(dat.baySeq, group='DE', 
                          number=Inf)[order(topCounts(dat.baySeq, group='DE', number=Inf)$name),]
  # $annotation on cluster (presumably different in older version), $name on my laptop
  prob.baySeq <- res.baySeq$likes
  # $Likelihood on cluster (presumably different in older version), $likes on my laptop
  q.baySeq <- res.baySeq$FDR.DE
  rm(list=c('dat.baySeq', 'res.baySeq'))
  
  # MDSeq
  contrasts <- get.model.matrix(group)
  fit.zi.MDSeq <- MDSeq(norm.DESeq2, contrast=contrasts, mc.cores=cores-1)
  fit.nozi.MDSeq <- MDSeq(norm.DESeq2, contrast=contrasts, test.ZI=F, mc.cores=cores-1)
  res.zi.MDSeq <- extract.ZIMD(fit.zi.MDSeq, compare=list(A="1",B="2"))
  res.nozi.MDSeq <- extract.ZIMD(fit.nozi.MDSeq, compare=list(A="1",B="2"))
  p.mean.zi.MDSeq <- res.zi.MDSeq$Pvalue.mean
  p.mean.nozi.MDSeq <- res.nozi.MDSeq$Pvalue.mean
  q.mean.zi.MDSeq <- res.zi.MDSeq$FDR.mean
  q.mean.nozi.MDSeq <- res.nozi.MDSeq$FDR.mean
  rm(list=c('contrasts', 'fit.zi.MDSeq', 'fit.nozi.MDSeq', 'res.zi.MDSeq', 'res.nozi.MDSeq'))
  
  # expHM
  expHM <- exp_hmm_adapt_3_chains(counts=t(norm.DESeq2), groups=group)
  prob.expHM <- colMeans(as.matrix(expHM$indicators))
  post.prop.expHM <- mean(as.matrix(expHM$proportion))
  mean.diff.expHM <- as.matrix(expHM$means1) - as.matrix(expHM$means2)
  log.mean.diff.expHM <- log(as.matrix(expHM$means1)) - log(as.matrix(expHM$means2))
  p.mean.expHM <- apply(mean.diff.expHM,2,hpd.pval, m=0)
  p.lmean.expHM <- apply(log.mean.diff.expHM,2,hpd.pval, m=0)
  q.mean.expHM <- p.adjust(p.mean.expHM, method="BH")
  q.lmean.expHM <- p.adjust(p.lmean.expHM, method="BH")
  rm(list=c('expHM', 'mean.diff.expHM', 'log.mean.diff.expHM'))
  
  # lnHM
  lnHM <- ln_hmm_adapt_3_chains(counts=t(norm.DESeq2), groups=group)
  prob.lnHM <- colMeans(as.matrix(lnHM$indicators))
  post.prop.lnHM <- mean(as.matrix(lnHM$proportion))
  mean.diff.lnHM <- as.matrix(lnHM$means1) - as.matrix(lnHM$means2)
  log.mean.diff.lnHM <- log(as.matrix(lnHM$means1)) - log(as.matrix(lnHM$means2))
  p.mean.lnHM <- apply(mean.diff.lnHM,2,hpd.pval)
  p.lmean.lnHM <- apply(log.mean.diff.lnHM,2,hpd.pval)
  q.mean.lnHM <- p.adjust(p.mean.lnHM, method="BH")
  q.lmean.lnHM <- p.adjust(p.lmean.lnHM, method="BH")
  rm(list=c('lnHM', 'mean.diff.lnHM', 'log.mean.diff.lnHM'))
  
  # Results
  results <- list(counts = counts, 
                  DE = DE, 
                  FC = FC, 
                  p.ql.edgeR = p.ql.edgeR, 
                  p.lr.edgeR = p.lr.edgeR, 
                  p.et.edgeR = p.et.edgeR, 
                  q.ql.edgeR = q.ql.edgeR, 
                  q.lr.edgeR = q.lr.edgeR, 
                  q.et.edgeR = q.et.edgeR, 
                  p.noif.DESeq2 = p.noif.DESeq2, 
                  p.if.DESeq2 = p.if.DESeq2, 
                  q.noif.DESeq2 = q.noif.DESeq2, 
                  q.if.DESeq2 = q.if.DESeq2, 
                  p.voom = p.voom, 
                  q.voom = q.voom, 
                  p.notrend.DSS = p.notrend.DSS, 
                  q.notrend.DSS = q.notrend.DSS, 
                  lfdr.notrend.DSS = lfdr.notrend.DSS, 
                  # p.trend.DSS = p.trend.DSS, 
                  # q.trend.DSS = q.trend.DSS, 
                  # lfdr.trend.DSS = lfdr.trend.DSS, 
                  prob.baySeq = prob.baySeq, 
                  q.baySeq = q.baySeq, 
                  p.mean.zi.MDSeq = p.mean.zi.MDSeq, 
                  p.mean.nozi.MDSeq = p.mean.nozi.MDSeq, 
                  q.mean.zi.MDSeq = q.mean.zi.MDSeq, 
                  q.mean.nozi.MDSeq = q.mean.nozi.MDSeq, 
                  prob.expHM = prob.expHM, 
                  prop.expHM = post.prop.expHM, 
                  p.mean.expHM = p.mean.expHM, 
                  p.lmean.expHM = p.lmean.expHM, 
                  q.mean.expHM = q.mean.expHM, 
                  q.lmean.expHM = q.lmean.expHM, 
                  prob.lnHM = prob.lnHM, 
                  prop.lnHM = post.prop.lnHM, 
                  p.mean.lnHM = p.mean.lnHM, 
                  p.lmean.lnHM = p.lmean.lnHM, 
                  q.mean.lnHM = q.mean.lnHM, 
                  q.lmean.lnHM = q.lmean.lnHM
  )
  
  filename <- paste0("results.", filename, ".rds")
  saveRDS(results, file=here("recount data/GTEx/quarantine", filename))
  rm("counts", "DE", "nf.DESeq2", "nf.TMM", "norm.DESeq2", "norm.TMM", 
     "filename", "results")
}


# DD ####
for (i in c(4, 10, 20, 100)) {
  # Data
  filename <- paste0("blood_", i, "_DD")
  data <- readRDS(here("recount data/GTEx", paste0(filename, ".rds")))
  counts <- data$counts
  DD <- data$DDindex
  FC <- data$FC
  group <- factor(c(rep(1,i/2), rep(2,i/2)))

  # Normalisation
  dat.DESeq2 <- DESeqDataSetFromMatrix(countData=counts, 
                                       colData=data.frame(group), 
                                       design=~group)
  dat.DESeq2 <- estimateSizeFactors(dat.DESeq2)
  nf.DESeq2 <- dat.DESeq2$sizeFactor
  norm.DESeq2 <- t(t(counts) / nf.DESeq2)
  
  # MDSeq
  contrasts <- get.model.matrix(group)
  fit.zi.MDSeq <- MDSeq(norm.DESeq2, contrast=contrasts, mc.cores=cores-1)
  fit.nozi.MDSeq <- MDSeq(norm.DESeq2, contrast=contrasts, test.ZI=F, mc.cores=cores-1)
  res.zi.MDSeq <- extract.ZIMD(fit.zi.MDSeq, compare=list(A="1",B="2"))
  res.nozi.MDSeq <- extract.ZIMD(fit.nozi.MDSeq, compare=list(A="1",B="2"))
  p.disp.zi.MDSeq <- res.zi.MDSeq$Pvalue.disp
  p.disp.nozi.MDSeq <- res.nozi.MDSeq$Pvalue.disp
  q.disp.zi.MDSeq <- res.zi.MDSeq$FDR.disp
  q.disp.nozi.MDSeq <- res.nozi.MDSeq$FDR.disp
  rm(list=c('contrasts', 'fit.zi.MDSeq', 'fit.nozi.MDSeq', 'res.zi.MDSeq', 'res.nozi.MDSeq'))
  
  # expHM
  expHM <- exp_hmm_adapt_3_chains(counts=t(norm.DESeq2), groups=group)
  prob.expHM <- colMeans(as.matrix(expHM$indicators))
  post.prop.expHM <- mean(as.matrix(expHM$proportion))
  disp.diff.expHM <- as.matrix(expHM$disps1) - as.matrix(expHM$disps2)
  log.disp.diff.expHM <- log(as.matrix(expHM$disps1)) - log(as.matrix(expHM$disps2))
  p.disp.expHM <- apply(disp.diff.expHM,2,hpd.pval, m=0)
  p.ldisp.expHM <- apply(log.disp.diff.expHM,2,hpd.pval, m=0)
  q.disp.expHM <- p.adjust(p.disp.expHM, method="BH")
  q.ldisp.expHM <- p.adjust(p.ldisp.expHM, method="BH")
  rm(list=c('expHM', 'disp.diff.expHM', 'log.disp.diff.expHM'))
  
  # lnHM
  lnHM <- ln_hmm_adapt_3_chains(counts=t(norm.DESeq2), groups=group)
  prob.lnHM <- colMeans(as.matrix(lnHM$indicators))
  post.prop.lnHM <- mean(as.matrix(lnHM$proportion))
  disp.diff.lnHM <- as.matrix(lnHM$disps1) - as.matrix(lnHM$disps2)
  log.disp.diff.lnHM <- log(as.matrix(lnHM$disps1)) - log(as.matrix(lnHM$disps2))
  p.disp.lnHM <- apply(disp.diff.lnHM,2,hpd.pval)
  p.ldisp.lnHM <- apply(log.disp.diff.lnHM,2,hpd.pval)
  q.disp.lnHM <- p.adjust(p.disp.lnHM, method="BH")
  q.ldisp.lnHM <- p.adjust(p.ldisp.lnHM, method="BH")
  rm(list=c('lnHM', 'disp.diff.lnHM', 'log.disp.diff.lnHM'))
  
  # Results
  results <- list(counts = counts, 
                  DD = DD, 
                  FC = FC, 
                  p.disp.zi.MDSeq = p.disp.zi.MDSeq, 
                  p.disp.nozi.MDSeq = p.disp.nozi.MDSeq, 
                  q.disp.zi.MDSeq = q.disp.zi.MDSeq, 
                  q.disp.nozi.MDSeq = q.disp.nozi.MDSeq, 
                  prob.expHM = prob.expHM, 
                  prop.expHM = post.prop.expHM, 
                  p.disp.expHM = p.disp.expHM, 
                  p.ldisp.expHM = p.ldisp.expHM, 
                  q.disp.expHM = q.disp.expHM, 
                  q.ldisp.expHM = q.ldisp.expHM, 
                  prob.lnHM = prob.lnHM, 
                  prop.lnHM = post.prop.lnHM, 
                  p.disp.lnHM = p.disp.lnHM, 
                  p.ldisp.lnHM = p.ldisp.lnHM, 
                  q.disp.lnHM = q.disp.lnHM, 
                  q.ldisp.lnHM = q.ldisp.lnHM
  )
  
  filename <- paste0("results.", filename, ".rds")
  saveRDS(results, file=here("recount data/GTEx/quarantine", filename))
  rm("counts", "DD", "nf.DESeq2", "norm.DESeq2", "filename", "results")
}


# DEDD ####
for (i in c(4, 10, 20, 100)) {
  # Data
  filename <- paste0("blood_", i, "_DEDD")
  data <- readRDS(here("recount data/GTEx", paste0(filename, ".rds")))
  counts <- data$counts
  DE <- data$DEindex
  DD <- data$DDindex
  DEDD <- data$DEDDindex
  DEonly <- data$DEonlyindex
  DDonly <- data$DDonlyindex
  DEplusDD <- data$DEplusDDindex
  FC.mean <- data$FC.mean
  FC.disp <- data$FC.disp
  group <- factor(c(rep(1,i/2), rep(2,i/2)))
  design <- model.matrix(~group)
  
  # Normalisation
  nf.TMM <- calcNormFactors(counts)
  norm.TMM <- t(t(counts) / nf.TMM)
  dat.edgeR <- DGEList(counts=counts, 
                       norm.factors=nf.TMM, 
                       group=group)
  dat.DESeq2 <- DESeqDataSetFromMatrix(countData=counts, 
                                       colData=data.frame(group), 
                                       design=~group)
  dat.DESeq2 <- estimateSizeFactors(dat.DESeq2)
  nf.DESeq2 <- dat.DESeq2$sizeFactor
  norm.DESeq2 <- t(t(counts) / nf.DESeq2)
  
  # edgeR
  dat.edgeR <- estimateDisp(dat.edgeR, design)
  qlfit.edgeR <- glmQLFit(dat.edgeR, design)
  qltest.edgeR <- glmQLFTest(qlfit.edgeR)
  lrfit.edgeR <- glmFit(dat.edgeR, design)
  lrtest.edgeR <- glmLRT(lrfit.edgeR)
  et.edgeR <- exactTest(dat.edgeR)
  p.ql.edgeR <- qltest.edgeR$table$PValue
  p.lr.edgeR <- lrtest.edgeR$table$PValue
  p.et.edgeR <- et.edgeR$table$PValue
  q.ql.edgeR <- topTags(qltest.edgeR, n=nrow(dat.edgeR$counts), sort='none')$table$FDR
  q.lr.edgeR <- topTags(lrtest.edgeR, n=nrow(dat.edgeR$counts), sort='none')$table$FDR
  q.et.edgeR <- topTags(et.edgeR, n=nrow(dat.edgeR$counts), sort='none')$table$FDR
  rm(qlfit.edgeR, qltest.edgeR, lrfit.edgeR, lrtest.edgeR, et.edgeR)
  
  # DESeq2
  dat.DESeq2 <- DESeq(dat.DESeq2, minReplicatesForReplace=Inf)
  res.if.DESeq2 <- results(dat.DESeq2, cooksCutoff=F, alpha=0.05)
  res.noif.DESeq2 <- results(dat.DESeq2, independentFiltering=F, cooksCutoff=F)
  p.noif.DESeq2 <- res.noif.DESeq2$pvalue
  p.if.DESeq2 <- res.if.DESeq2$pvalue
  q.noif.DESeq2 <- res.noif.DESeq2$padj
  q.if.DESeq2 <- res.if.DESeq2$padj
  rm(list=c('dat.DESeq2', 'res.noif.DESeq2', 'res.if.DESeq2'))
  
  # limma-voom
  dat.voom <- voom(dat.edgeR)
  fit.voom <- lmFit(dat.voom, design)
  res.voom <- eBayes(fit.voom)
  p.voom <- topTable(res.voom, number=Inf, sort.by='none')$P.Value
  q.voom <- topTable(res.voom, number=Inf, sort.by='none')$adj.P.Val
  rm(list=c('dat.edgeR', 'dat.voom', 'fit.voom', 'res.voom'))
  
  # DSS
  dat.DSS <- newSeqCountSet(counts=matrix(counts, ncol=length(group)), 
                            designs=as.numeric(group), 
                            normalizationFactor=nf.DESeq2)
  notrend.DSS <- estDispersion(dat.DSS)
  # trend.DSS <- estDispersion(dat.DSS, trend=T)
  res.notrend.DSS <- waldTest(notrend.DSS, 1, 2)[order(waldTest(notrend.DSS, 1, 2)$geneIndex),]
  # res.trend.DSS <- waldTest(trend.DSS, 1, 2)[order(waldTest(trend.DSS, 1, 2)$geneIndex),]
  p.notrend.DSS <- res.notrend.DSS$pval
  q.notrend.DSS <- res.notrend.DSS$fdr
  lfdr.notrend.DSS <- res.notrend.DSS$local.fdr
  # p.trend.DSS <- res.trend.DSS$pval
  # q.trend.DSS <- res.trend.DSS$fdr
  # lfdr.trend.DSS <- res.trend.DSS$local.fdr
  # rm(list=c('dat.DSS', 'notrend.DSS', 'trend.DSS', 'res.notrend.DSS', 'res.trend.DSS'))
  rm(list=c('dat.DSS', 'notrend.DSS', 'res.notrend.DSS'))
  
  # baySeq
  dat.baySeq <- new("countData",
                    data=counts, 
                    replicates=group, 
                    groups=list(NDE=rep(1, length(group)), DE=as.numeric(group)))
  dat.baySeq@annotation <- data.frame(name = 1:nrow(dat.baySeq@data))
  libsizes(dat.baySeq) <- nf.DESeq2
  dat.baySeq <- getPriors.NB(dat.baySeq, samplesize=10000, cl=cl)
  dat.baySeq <- getLikelihoods(dat.baySeq, cl=cl)
  res.baySeq <- topCounts(dat.baySeq, group='DE', 
                          number=Inf)[order(topCounts(dat.baySeq, group='DE', number=Inf)$name),]
  # $annotation on cluster (presumably different in older version), $name on my laptop
  prob.baySeq <- res.baySeq$likes
  # $Likelihood on cluster (presumably different in older version), $likes on my laptop
  q.baySeq <- res.baySeq$FDR.DE
  rm(list=c('dat.baySeq', 'res.baySeq'))
  
  # MDSeq
  contrasts <- get.model.matrix(group)
  fit.zi.MDSeq <- MDSeq(norm.DESeq2, contrast=contrasts, mc.cores=cores-1)
  fit.nozi.MDSeq <- MDSeq(norm.DESeq2, contrast=contrasts, test.ZI=F, mc.cores=cores-1)
  res.zi.MDSeq <- extract.ZIMD(fit.zi.MDSeq, compare=list(A="1",B="2"))
  res.nozi.MDSeq <- extract.ZIMD(fit.nozi.MDSeq, compare=list(A="1",B="2"))
  p.mean.zi.MDSeq <- res.zi.MDSeq$Pvalue.mean
  p.mean.nozi.MDSeq <- res.nozi.MDSeq$Pvalue.mean
  p.disp.zi.MDSeq <- res.zi.MDSeq$Pvalue.disp
  p.disp.nozi.MDSeq <- res.nozi.MDSeq$Pvalue.disp
  q.mean.zi.MDSeq <- res.zi.MDSeq$FDR.mean
  q.mean.nozi.MDSeq <- res.nozi.MDSeq$FDR.mean
  q.disp.zi.MDSeq <- res.zi.MDSeq$FDR.disp
  q.disp.nozi.MDSeq <- res.nozi.MDSeq$FDR.disp
  rm(list=c('contrasts', 'fit.zi.MDSeq', 'fit.nozi.MDSeq', 'res.zi.MDSeq', 'res.nozi.MDSeq'))
  
  # expHM
  expHM <- exp_hmm_adapt_3_chains(counts=t(norm.DESeq2), groups=group)
  prob.expHM <- colMeans(as.matrix(expHM$indicators))
  post.prop.expHM <- mean(as.matrix(expHM$proportion))
  mean.diff.expHM <- as.matrix(expHM$means1) - as.matrix(expHM$means2)
  log.mean.diff.expHM <- log(as.matrix(expHM$means1)) - log(as.matrix(expHM$means2))
  disp.diff.expHM <- as.matrix(expHM$disps1) - as.matrix(expHM$disps2)
  log.disp.diff.expHM <- log(as.matrix(expHM$disps1)) - log(as.matrix(expHM$disps2))
  p.mean.expHM <- apply(mean.diff.expHM,2,hpd.pval, m=0)
  p.lmean.expHM <- apply(log.mean.diff.expHM,2,hpd.pval, m=0)
  p.disp.expHM <- apply(disp.diff.expHM,2,hpd.pval, m=0)
  p.ldisp.expHM <- apply(log.disp.diff.expHM,2,hpd.pval, m=0)
  q.mean.expHM <- p.adjust(p.mean.expHM, method="BH")
  q.lmean.expHM <- p.adjust(p.lmean.expHM, method="BH")
  q.disp.expHM <- p.adjust(p.disp.expHM, method="BH")
  q.ldisp.expHM <- p.adjust(p.ldisp.expHM, method="BH")
  rm(list=c('expHM', 'mean.diff.expHM', 'log.mean.diff.expHM', 
            'disp.diff.expHM', 'log.disp.diff.expHM'))
  
  # lnHM
  lnHM <- ln_hmm_adapt_3_chains(counts=t(norm.DESeq2), groups=group)
  prob.lnHM <- colMeans(as.matrix(lnHM$indicators))
  post.prop.lnHM <- mean(as.matrix(lnHM$proportion))
  mean.diff.lnHM <- as.matrix(lnHM$means1) - as.matrix(lnHM$means2)
  log.mean.diff.lnHM <- log(as.matrix(lnHM$means1)) - log(as.matrix(lnHM$means2))
  disp.diff.lnHM <- as.matrix(lnHM$disps1) - as.matrix(lnHM$disps2)
  log.disp.diff.lnHM <- log(as.matrix(lnHM$disps1)) - log(as.matrix(lnHM$disps2))
  p.mean.lnHM <- apply(mean.diff.lnHM,2,hpd.pval)
  p.lmean.lnHM <- apply(log.mean.diff.lnHM,2,hpd.pval)
  p.disp.lnHM <- apply(disp.diff.lnHM,2,hpd.pval)
  p.ldisp.lnHM <- apply(log.disp.diff.lnHM,2,hpd.pval)
  q.mean.lnHM <- p.adjust(p.mean.lnHM, method="BH")
  q.lmean.lnHM <- p.adjust(p.lmean.lnHM, method="BH")
  q.disp.lnHM <- p.adjust(p.disp.lnHM, method="BH")
  q.ldisp.lnHM <- p.adjust(p.ldisp.lnHM, method="BH")
  rm(list=c('lnHM', 'mean.diff.lnHM', 'log.mean.diff.lnHM', 
            'disp.diff.lnHM', 'log.disp.diff.lnHM'))
  
  # Results
  results <- list(counts = counts, 
                  DE = DE, 
                  DD = DD, 
                  DEDD = DEDD, 
                  DEonly = DEonly, 
                  DDonly = DDonly, 
                  DEplusDD = DEplusDD, 
                  FC.mean = FC.mean, 
                  FC.disp = FC.disp, 
                  p.ql.edgeR = p.ql.edgeR, 
                  p.lr.edgeR = p.lr.edgeR, 
                  p.et.edgeR = p.et.edgeR, 
                  q.ql.edgeR = q.ql.edgeR, 
                  q.lr.edgeR = q.lr.edgeR, 
                  q.et.edgeR = q.et.edgeR, 
                  p.noif.DESeq2 = p.noif.DESeq2, 
                  p.if.DESeq2 = p.if.DESeq2, 
                  q.noif.DESeq2 = q.noif.DESeq2, 
                  q.if.DESeq2 = q.if.DESeq2, 
                  p.voom = p.voom, 
                  q.voom = q.voom, 
                  p.notrend.DSS = p.notrend.DSS, 
                  q.notrend.DSS = q.notrend.DSS, 
                  lfdr.notrend.DSS = lfdr.notrend.DSS, 
                  # p.trend.DSS = p.trend.DSS, 
                  # q.trend.DSS = q.trend.DSS, 
                  # lfdr.trend.DSS = lfdr.trend.DSS, 
                  prob.baySeq = prob.baySeq, 
                  q.baySeq = q.baySeq, 
                  p.mean.zi.MDSeq = p.mean.zi.MDSeq, 
                  p.mean.nozi.MDSeq = p.mean.nozi.MDSeq, 
                  p.disp.zi.MDSeq = p.disp.zi.MDSeq, 
                  p.disp.nozi.MDSeq = p.disp.nozi.MDSeq, 
                  q.mean.zi.MDSeq = q.mean.zi.MDSeq, 
                  q.mean.nozi.MDSeq = q.mean.nozi.MDSeq, 
                  q.disp.zi.MDSeq = q.disp.zi.MDSeq, 
                  q.disp.nozi.MDSeq = q.disp.nozi.MDSeq, 
                  prob.expHM = prob.expHM, 
                  prop.expHM = post.prop.expHM, 
                  p.mean.expHM = p.mean.expHM, 
                  p.lmean.expHM = p.lmean.expHM, 
                  p.disp.expHM = p.disp.expHM, 
                  p.ldisp.expHM = p.ldisp.expHM, 
                  q.mean.expHM = q.mean.expHM, 
                  q.lmean.expHM = q.lmean.expHM, 
                  q.disp.expHM = q.disp.expHM, 
                  q.ldisp.expHM = q.ldisp.expHM, 
                  prob.lnHM = prob.lnHM, 
                  prop.lnHM = post.prop.lnHM, 
                  p.mean.lnHM = p.mean.lnHM, 
                  p.lmean.lnHM = p.lmean.lnHM, 
                  p.disp.lnHM = p.disp.lnHM, 
                  p.ldisp.lnHM = p.ldisp.lnHM, 
                  q.mean.lnHM = q.mean.lnHM, 
                  q.lmean.lnHM = q.lmean.lnHM, 
                  q.disp.lnHM = q.disp.lnHM, 
                  q.ldisp.lnHM = q.ldisp.lnHM
  )
  
  filename <- paste0("results.", filename, ".rds")
  saveRDS(results, file=here("recount data/GTEx/quarantine", filename))
  rm("counts", "DE", "DD", "DEDD", "DEonly", "DDonly", "DEplusDD",  
     "nf.DESeq2", "nf.TMM", "norm.DESeq2", "norm.TMM", 
     "filename", "results")
}


saveRDS(sessionInfo(), file=here("recount data/GTEx/quarantine", "blood_sessionInfo.rds"))




