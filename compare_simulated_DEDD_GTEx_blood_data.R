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

samples.per.cond <- 50
group <- factor(c(rep(1, samples.per.cond), rep(2, samples.per.cond)))
design <- model.matrix(~group)
contrasts <- get.model.matrix(group)
# cores <- detectCores() - 1
if(require("parallel")) cl <- makeCluster(4) else cl <- NULL

for (i in 10) {
  # Data
  if (cluster) {
    folder <- paste0("Data/recount data/GTEx/blood_", samples.per.cond, "_samples_per_group")
  } else {
    folder <- paste0("recount data/GTEx/blood_", samples.per.cond, "_samples_per_group")
  }
  filename <- paste0("blood_", samples.per.cond, "_set", i, "_DEDD")
  data <- readRDS(here(folder, paste0(filename, ".rds")))
  counts <- data$counts
  DE <- data$DEindex
  DD <- data$DDindex
  DEDD <- as.numeric(DE == 1 | DD == 1)
  FC.mean <- data$FC.mean
  FC.disp <- data$FC.disp
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
  
  # Create data objects
  dat.edgeR.tmm <- DGEList(counts=counts, norm.factors=nf.tmm, group=group)
  dat.edgeR.tmm <- estimateDisp(dat.edgeR.tmm, design)
  dat.edgeR.rle <- DGEList(counts=counts, norm.factors=nf.rle, group=group)
  dat.edgeR.rle <- estimateDisp(dat.edgeR.rle, design)
  dat.DESeq2.tmm <- DESeqDataSetFromMatrix(countData=counts, colData=data.frame(group), design=~group)
  sizeFactors(dat.DESeq2.tmm) <- sf.tmm
  dat.DESeq2.rle <- DESeqDataSetFromMatrix(countData=counts, colData=data.frame(group), design=~group)
  dat.DESeq2.rle <- estimateSizeFactors(dat.DESeq2.rle)
  dat.DSS.tmm <- newSeqCountSet(counts=matrix(counts, ncol=length(group)), 
                                designs=as.numeric(group), normalizationFactor=sf.tmm)
  dat.DSS.tmm <- estDispersion(dat.DSS.tmm)
  dat.DSS.rle <- newSeqCountSet(counts=matrix(counts, ncol=length(group)), 
                                designs=as.numeric(group), normalizationFactor=sf.rle)
  dat.DSS.rle <- estDispersion(dat.DSS.rle)
  dat.baySeq.tmm <- new('countData', data=counts, replicates=group, 
                        groups=list(NDE=rep(1,length(group)), DE=as.numeric(group)))
  dat.baySeq.tmm@annotation <- data.frame(name = 1:nrow(dat.baySeq.tmm@data))
  libsizes(dat.baySeq.tmm) <- els.tmm
  dat.baySeq.rle <- new('countData', data=counts, replicates=group, 
                        groups=list(NDE=rep(1,length(group)), DE=as.numeric(group)))
  dat.baySeq.rle@annotation <- data.frame(name = 1:nrow(dat.baySeq.rle@data))
  libsizes(dat.baySeq.rle) <- els.rle
  
  ## edgeR
  fit.edgeR.ql.tmm <- glmQLFit(dat.edgeR.tmm, design)
  test.edgeR.ql.tmm <- glmQLFTest(fit.edgeR.ql.tmm)
  fit.edgeR.lr.tmm <- glmFit(dat.edgeR.tmm, design)
  test.edgeR.lr.tmm <- glmLRT(fit.edgeR.lr.tmm)
  edgeR.et.tmm <- exactTest(dat.edgeR.tmm)
  p.edgeR.ql.tmm <- test.edgeR.ql.tmm$table$PValue
  p.edgeR.lr.tmm <- test.edgeR.lr.tmm$table$PValue
  p.edgeR.et.tmm <- edgeR.et.tmm$table$PValue
  q.edgeR.ql.tmm <- topTags(test.edgeR.ql.tmm, n=nrow(dat.edgeR.tmm$counts), sort='none')$table$FDR
  q.edgeR.lr.tmm <- topTags(test.edgeR.lr.tmm, n=nrow(dat.edgeR.tmm$counts), sort='none')$table$FDR
  q.edgeR.et.tmm <- topTags(edgeR.et.tmm, n=nrow(dat.edgeR.tmm$counts), sort='none')$table$FDR
  rm(fit.edgeR.ql.tmm, test.edgeR.ql.tmm, fit.edgeR.lr.tmm, test.edgeR.lr.tmmedgeR.et.tmm)
  fit.edgeR.ql.rle <- glmQLFit(dat.edgeR.rle, design)
  test.edgeR.ql.rle <- glmQLFTest(fit.edgeR.ql.rle)
  fit.edgeR.lr.rle <- glmFit(dat.edgeR.rle, design)
  test.edgeR.lr.rle <- glmLRT(fit.edgeR.lr.rle)
  edgeR.et.rle <- exactTest(dat.edgeR.rle)
  p.edgeR.ql.rle <- test.edgeR.ql.rle$table$PValue
  p.edgeR.lr.rle <- test.edgeR.lr.rle$table$PValue
  p.edgeR.et.rle <- edgeR.et.rle$table$PValue
  q.edgeR.ql.rle <- topTags(test.edgeR.ql.rle, n=nrow(dat.edgeR.rle$counts), sort='none')$table$FDR
  q.edgeR.lr.rle <- topTags(test.edgeR.lr.rle, n=nrow(dat.edgeR.rle$counts), sort='none')$table$FDR
  q.edgeR.et.rle <- topTags(edgeR.et.rle, n=nrow(dat.edgeR.rle$counts), sort='none')$table$FDR
  rm(fit.edgeR.ql.rle, test.edgeR.ql.rle, fit.edgeR.lr.rle, test.edgeR.lr.rle, edgeR.et.rle)
  
  ## DESeq2
  fit.DESeq2.tmm <- DESeq(dat.DESeq2.tmm, minReplicatesForReplace=Inf)
  res.DESeq2.noif.tmm <- results(fit.DESeq2.tmm, independentFiltering=F, cooksCutoff=F)
  res.DESeq2.if.tmm <- results(fit.DESeq2.tmm, cooksCutoff=F, alpha=0.05)
  p.DESeq2.noif.tmm <- res.DESeq2.noif.tmm$pvalue
  p.DESeq2.if.tmm <- res.DESeq2.if.tmm$pvalue
  q.DESeq2.noif.tmm <- res.DESeq2.noif.tmm$padj
  q.DESeq2.if.tmm <- res.DESeq2.if.tmm$padj
  rm(dat.DESeq2.tmm, fit.DESeq2.tmm, res.DESeq2.noif.tmm, res.DESeq2.if.tmm)
  fit.DESeq2.rle <- DESeq(dat.DESeq2.rle, minReplicatesForReplace=Inf)
  res.DESeq2.noif.rle <- results(fit.DESeq2.rle, independentFiltering=F, cooksCutoff=F)
  res.DESeq2.if.rle <- results(fit.DESeq2.rle, cooksCutoff=F, alpha=0.05)
  p.DESeq2.noif.rle <- res.DESeq2.noif.rle$pvalue
  p.DESeq2.if.rle <- res.DESeq2.if.rle$pvalue
  q.DESeq2.noif.rle <- res.DESeq2.noif.rle$padj
  q.DESeq2.if.rle <- res.DESeq2.if.rle$padj
  rm(dat.DESeq2.rle, fit.DESeq2.rle, res.DESeq2.noif.rle, res.DESeq2.if.rle)
  
  ## limma-voom
  dat.voom.tmm <- voom(dat.edgeR.tmm)
  fit.voom.tmm <- lmFit(dat.voom.tmm, design)
  res.voom.tmm <- eBayes(fit.voom.tmm)
  p.voom.tmm <- topTable(res.voom.tmm, number=Inf, sort.by='none')$P.Value
  q.voom.tmm <- topTable(res.voom.tmm, number=Inf, sort.by='none')$adj.P.Val
  rm(dat.edgeR.tmm, dat.voom.tmm, fit.voom.tmm, res.voom.tmm)
  dat.voom.rle <- voom(dat.edgeR.rle)
  fit.voom.rle <- lmFit(dat.voom.rle, design)
  res.voom.rle <- eBayes(fit.voom.rle)
  p.voom.rle <- topTable(res.voom.rle, number=Inf, sort.by='none')$P.Value
  q.voom.rle <- topTable(res.voom.rle, number=Inf, sort.by='none')$adj.P.Val
  rm(dat.edgeR.rle, dat.voom.rle, fit.voom.rle, res.voom.rle)
  
  ## DSS
  res.DSS.tmm <- waldTest(dat.DSS.tmm, 1, 2)[order(waldTest(dat.DSS.tmm, 1, 2)$geneIndex),]
  p.DSS.tmm <- res.DSS.tmm$pval
  q.DSS.tmm <- res.DSS.tmm$fdr
  lfdr.DSS.tmm <- res.DSS.tmm$local.fdr
  rm(dat.DSS.tmm, res.DSS.tmm)
  res.DSS.rle <- waldTest(dat.DSS.rle, 1, 2)[order(waldTest(dat.DSS.rle, 1, 2)$geneIndex),]
  p.DSS.rle <- res.DSS.rle$pval
  q.DSS.rle <- res.DSS.rle$fdr
  lfdr.DSS.rle <- res.DSS.rle$local.fdr
  rm(dat.DSS.rle, res.DSS.rle)
  
  ## baySeq
  dat.baySeq.tmm <- getPriors.NB(dat.baySeq.tmm, samplesize=10000, cl=cl)
  # dat.baySeq.tmm <- getLikelihoods(dat.baySeq.tmm, cl=cl) # error on cluster using cl here
  dat.baySeq.tmm <- getLikelihoods(dat.baySeq.tmm)
  res.baySeq.tmm <- topCounts(dat.baySeq.tmm, group='DE', number=Inf)[
    order(topCounts(dat.baySeq.tmm, group='DE', number=Inf)$name), ]
  prob.baySeq.tmm <- res.baySeq.tmm$likes
  q.baySeq.tmm <- res.baySeq.tmm$FDR.DE
  rm(dat.baySeq.tmm, res.baySeq.tmm)
  dat.baySeq.rle <- getPriors.NB(dat.baySeq.rle, samplesize=10000, cl=cl)
  # dat.baySeq.rle <- getLikelihoods(dat.baySeq.rle, cl=cl) error on cluster using cl here
  dat.baySeq.rle <- getLikelihoods(dat.baySeq.rle)
  res.baySeq.rle <- topCounts(dat.baySeq.rle, group='DE', number=Inf)[
    order(topCounts(dat.baySeq.rle, group='DE', number=Inf)$name), ]
  prob.baySeq.rle <- res.baySeq.rle$likes
  q.baySeq.rle <- res.baySeq.rle$FDR.DE
  # $name may be $annotation in older versions
  # $likes may be $Likelihood in older versions
  rm(dat.baySeq.rle, res.baySeq.rle)
  
  ## MDSeq
  # fit.MDSeq.zi.tmm <- MDSeq(counts, offsets=sf.tmm, contrast=contrasts, mc.cores=cores-1)
  # hangs on Phoenix with mc.cores
  fit.MDSeq.zi.tmm <- MDSeq(counts, offsets=sf.tmm, contrast=contrasts)
  res.MDSeq.zi.tmm <- extract.ZIMD(fit.MDSeq.zi.tmm, compare=list(A="1",B="2"))
  # fit.MDSeq.nozi.tmm <- MDSeq(counts, offsets=sf.tmm, contrast=contrasts, test.ZI=F, mc.cores=cores-1)
  # hangs on Phoenix with mc.cores
  fit.MDSeq.nozi.tmm <- MDSeq(counts, offsets=sf.tmm, contrast=contrasts, test.ZI=F)
  res.MDSeq.nozi.tmm <- extract.ZIMD(fit.MDSeq.nozi.tmm, compare=list(A="1",B="2"))
  p.mean.MDSeq.zi.tmm <- res.MDSeq.zi.tmm$Pvalue.mean
  q.mean.MDSeq.zi.tmm <- res.MDSeq.zi.tmm$FDR.mean
  p.mean.MDSeq.nozi.tmm <- res.MDSeq.nozi.tmm$Pvalue.mean
  q.mean.MDSeq.nozi.tmm <- res.MDSeq.nozi.tmm$FDR.mean
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
  p.mean.MDSeq.zi.rle <- res.MDSeq.zi.rle$Pvalue.mean
  q.mean.MDSeq.zi.rle <- res.MDSeq.zi.rle$FDR.mean
  p.mean.MDSeq.nozi.rle <- res.MDSeq.nozi.rle$Pvalue.mean
  q.mean.MDSeq.nozi.rle <- res.MDSeq.nozi.rle$FDR.mean
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
  mean.diff.expHM.untr.tmm <- unname(as.matrix(expHM.tmm$means1) - as.matrix(expHM.tmm$means2))
  p.mean.expHM.untr.tmm <- apply(mean.diff.expHM.untr.tmm,2,hpd.pval)
  rm(mean.diff.expHM.untr.tmm); gc()
  mean.diff.expHM.log.tmm <- unname(log(as.matrix(expHM.tmm$means1)) - log(as.matrix(expHM.tmm$means2)))
  p.mean.expHM.log.tmm <- apply(mean.diff.expHM.log.tmm,2,hpd.pval)
  rm(mean.diff.expHM.log.tmm); gc()
  q.mean.expHM.untr.tmm <- p.adjust(p.mean.expHM.untr.tmm, method='BH')
  q.mean.expHM.log.tmm <- p.adjust(p.mean.expHM.log.tmm, method='BH')
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
  mean.diff.expHM.untr.rle <- unname(as.matrix(expHM.rle$means1) - as.matrix(expHM.rle$means2))
  p.mean.expHM.untr.rle <- apply(mean.diff.expHM.untr.rle,2,hpd.pval)
  rm(mean.diff.expHM.untr.rle); gc()
  mean.diff.expHM.log.rle <- unname(log(as.matrix(expHM.rle$means1)) - log(as.matrix(expHM.rle$means2)))
  p.mean.expHM.log.rle <- apply(mean.diff.expHM.log.rle,2,hpd.pval)
  rm(mean.diff.expHM.log.rle); gc()
  q.mean.expHM.untr.rle <- p.adjust(p.mean.expHM.untr.rle, method='BH')
  q.mean.expHM.log.rle <- p.adjust(p.mean.expHM.log.rle, method='BH')
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
  mean.diff.lnHM.untr.tmm <- unname(as.matrix(lnHM.tmm$means1) - as.matrix(lnHM.tmm$means2))
  p.mean.lnHM.untr.tmm <- apply(mean.diff.lnHM.untr.tmm,2,hpd.pval)
  rm(mean.diff.lnHM.untr.tmm); gc()
  mean.diff.lnHM.log.tmm <- unname(log(as.matrix(lnHM.tmm$means1)) - log(as.matrix(lnHM.tmm$means2)))
  p.mean.lnHM.log.tmm <- apply(mean.diff.lnHM.log.tmm,2,hpd.pval)
  rm(mean.diff.lnHM.log.tmm); gc()
  q.mean.lnHM.untr.tmm <- p.adjust(p.mean.lnHM.untr.tmm, method='BH')
  q.mean.lnHM.log.tmm <- p.adjust(p.mean.lnHM.log.tmm, method='BH')
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
  mean.diff.lnHM.untr.rle <- unname(as.matrix(lnHM.rle$means1) - as.matrix(lnHM.rle$means2))
  p.mean.lnHM.untr.rle <- apply(mean.diff.lnHM.untr.rle,2,hpd.pval)
  rm(mean.diff.lnHM.untr.rle); gc()
  mean.diff.lnHM.log.rle <- unname(log(as.matrix(lnHM.rle$means1)) - log(as.matrix(lnHM.rle$means2)))
  p.mean.lnHM.log.rle <- apply(mean.diff.lnHM.log.rle,2,hpd.pval)
  rm(mean.diff.lnHM.log.rle); gc()
  q.mean.lnHM.untr.rle <- p.adjust(p.mean.lnHM.untr.rle, method='BH')
  q.mean.lnHM.log.rle <- p.adjust(p.mean.lnHM.log.rle, method='BH')
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
                  DE = DE, 
                  DD = DD, 
                  DEDD = DEDD, 
                  p.edgeR.ql.tmm = p.edgeR.ql.tmm, 
                  p.edgeR.lr.tmm = p.edgeR.lr.tmm, 
                  p.edgeR.et.tmm = p.edgeR.et.tmm, 
                  q.edgeR.ql.tmm = q.edgeR.ql.tmm, 
                  q.edgeR.lr.tmm = q.edgeR.lr.tmm, 
                  q.edgeR.et.tmm = q.edgeR.et.tmm, 
                  p.edgeR.ql.rle = p.edgeR.ql.rle, 
                  p.edgeR.lr.rle = p.edgeR.lr.rle, 
                  p.edgeR.et.rle = p.edgeR.et.rle, 
                  q.edgeR.ql.rle = q.edgeR.ql.rle, 
                  q.edgeR.lr.rle = q.edgeR.lr.rle, 
                  q.edgeR.et.rle = q.edgeR.et.rle, 
                  p.DESeq2.noif.tmm = p.DESeq2.noif.tmm, 
                  p.DESeq2.if.tmm = p.DESeq2.if.tmm, 
                  q.DESeq2.noif.tmm = q.DESeq2.noif.tmm, 
                  q.DESeq2.if.tmm = q.DESeq2.if.tmm, 
                  p.DESeq2.noif.rle = p.DESeq2.noif.rle, 
                  p.DESeq2.if.rle = p.DESeq2.if.rle, 
                  q.DESeq2.noif.rle = q.DESeq2.noif.rle, 
                  q.DESeq2.if.rle = q.DESeq2.if.rle, 
                  p.voom.tmm = p.voom.tmm, 
                  q.voom.tmm = q.voom.tmm, 
                  p.voom.rle = p.voom.rle, 
                  q.voom.rle = q.voom.rle, 
                  p.DSS.tmm = p.DSS.tmm, 
                  q.DSS.tmm = q.DSS.tmm, 
                  lfdr.DSS.tmm = lfdr.DSS.tmm, 
                  p.DSS.rle = p.DSS.rle, 
                  q.DSS.rle = q.DSS.rle, 
                  lfdr.DSS.rle = lfdr.DSS.rle, 
                  prob.baySeq.tmm = prob.baySeq.tmm, 
                  q.baySeq.tmm = q.baySeq.tmm, 
                  prob.baySeq.rle = prob.baySeq.rle, 
                  q.baySeq.rle = q.baySeq.rle, 
                  p.mean.MDSeq.zi.tmm = p.mean.MDSeq.zi.tmm, 
                  p.mean.MDSeq.nozi.tmm = p.mean.MDSeq.nozi.tmm, 
                  q.mean.MDSeq.zi.tmm = q.mean.MDSeq.zi.tmm, 
                  q.mean.MDSeq.nozi.tmm = q.mean.MDSeq.nozi.tmm, 
                  p.disp.MDSeq.zi.tmm = p.disp.MDSeq.zi.tmm, 
                  p.disp.MDSeq.nozi.tmm = p.disp.MDSeq.nozi.tmm, 
                  q.disp.MDSeq.zi.tmm = q.disp.MDSeq.zi.tmm, 
                  q.disp.MDSeq.nozi.tmm = q.disp.MDSeq.nozi.tmm, 
                  p.mean.MDSeq.zi.rle = p.mean.MDSeq.zi.rle, 
                  p.mean.MDSeq.nozi.rle = p.mean.MDSeq.nozi.rle, 
                  q.mean.MDSeq.zi.rle = q.mean.MDSeq.zi.rle, 
                  q.mean.MDSeq.nozi.rle = q.mean.MDSeq.nozi.rle, 
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
                  p.mean.expHM.untr.tmm = p.mean.expHM.untr.tmm, 
                  p.mean.expHM.log.tmm = p.mean.expHM.log.tmm, 
                  p.disp.expHM.untr.tmm = p.disp.expHM.untr.tmm, 
                  p.disp.expHM.log.tmm = p.disp.expHM.log.tmm, 
                  q.mean.expHM.untr.tmm = q.mean.expHM.untr.tmm, 
                  q.mean.expHM.log.tmm = q.mean.expHM.log.tmm, 
                  q.disp.expHM.untr.tmm = q.disp.expHM.untr.tmm, 
                  q.disp.expHM.log.tmm = q.disp.expHM.log.tmm, 
                  prob.expHMM.rle = prob.expHMM.rle, 
                  prop.expHMM.rle = post.prop.expHMM.rle, 
                  thr.expHMM.rle = thr.expHMM.rle, 
                  bfdr.expHMM.rle = bfdr.expHMM.rle, 
                  p.mean.expHM.untr.rle = p.mean.expHM.untr.rle, 
                  p.mean.expHM.log.rle = p.mean.expHM.log.rle, 
                  p.disp.expHM.untr.rle = p.disp.expHM.untr.rle, 
                  p.disp.expHM.log.rle = p.disp.expHM.log.rle, 
                  q.mean.expHM.untr.rle = q.mean.expHM.untr.rle, 
                  q.mean.expHM.log.rle = q.mean.expHM.log.rle, 
                  q.disp.expHM.untr.rle = q.disp.expHM.untr.rle, 
                  q.disp.expHM.log.rle = q.disp.expHM.log.rle, 
                  prob.lnHMM.tmm = prob.lnHMM.tmm, 
                  prop.lnHMM.tmm = post.prop.lnHMM.tmm, 
                  thr.lnHMM.tmm = thr.lnHMM.tmm, 
                  bfdr.lnHMM.tmm = bfdr.lnHMM.tmm, 
                  p.mean.lnHM.untr.tmm = p.mean.lnHM.untr.tmm, 
                  p.mean.lnHM.log.tmm = p.mean.lnHM.log.tmm, 
                  p.disp.lnHM.untr.tmm = p.disp.lnHM.untr.tmm, 
                  p.disp.lnHM.log.tmm = p.disp.lnHM.log.tmm, 
                  q.mean.lnHM.untr.tmm = q.mean.lnHM.untr.tmm, 
                  q.mean.lnHM.log.tmm = q.mean.lnHM.log.tmm, 
                  q.disp.lnHM.untr.tmm = q.disp.lnHM.untr.tmm, 
                  q.disp.lnHM.log.tmm = q.disp.lnHM.log.tmm, 
                  prob.lnHMM.rle = prob.lnHMM.rle, 
                  prop.lnHMM.rle = post.prop.lnHMM.rle, 
                  thr.lnHMM.rle = thr.lnHMM.rle, 
                  bfdr.lnHMM.rle = bfdr.lnHMM.rle, 
                  p.mean.lnHM.untr.rle = p.mean.lnHM.untr.rle, 
                  p.mean.lnHM.log.rle = p.mean.lnHM.log.rle, 
                  p.disp.lnHM.untr.rle = p.disp.lnHM.untr.rle, 
                  p.disp.lnHM.log.rle = p.disp.lnHM.log.rle, 
                  q.mean.lnHM.untr.rle = q.mean.lnHM.untr.rle, 
                  q.mean.lnHM.log.rle = q.mean.lnHM.log.rle, 
                  q.disp.lnHM.untr.rle = q.disp.lnHM.untr.rle, 
                  q.disp.lnHM.log.rle = q.disp.lnHM.log.rle)
  
  if (cluster) {
    folder <- paste0("Data/GTEx blood simulated DE, DD, DEDD results Feb 2020")
  } else {
    folder <- paste0("Results/GTEx blood simulated DE, DD, DEDD results Feb 2020")
  }
  filename <- paste0("results.", filename, ".rds")
  saveRDS(results, file=here(folder, filename))
  
  rm(counts, DE, DD, DEDD, FC.mean, FC.disp, libsizes, 
     nf.tmm, els.tmm, sf.tmm, norm.counts.tmm, 
	   nf.rle, els.rle, sf.rle, norm.counts.rle)
  
}

filename <- paste0("sessionInfo.GTEx_blood_", samples.per.cond, "_set", i, "_DEDD.rds")
saveRDS(sessionInfo, file=here(folder, filename))
