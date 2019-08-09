library(here)
library(coda)
library(HDInterval)
library(ROCR)
library(qvalue)
library(compcodeR)
library(limma)
library(edgeR)
library(DESeq2)
library(DSS)
#library(baySeq)
#library(ShrinkBayes)
#library(MDSeq)
source(here('scripts','2019-04-03_exponential_hmm_adaptive_proposals_three_chains_function.R'))
source(here('scripts','2019-03-27_lognormal_hmm_adaptive_proposals_three_chains_function.R'))
source(here('scripts','2019-05-03_bfdr_function.R'))
source(here('scripts','2019-05-03_hpd_tail_prob_function.R'))
source(here('scripts','2019-05-07_symmetric_tail_prob_function.R'))
source(here('scripts','2019-05-17_compData_diff_disp_functions.R'))

samples.per.cond <- 5
group <- factor(c(rep(1,samples.per.cond), rep(2,samples.per.cond)))
FPR.DEDD2.ql.edgeR <- numeric(5); FPR.DEDD2.lr.edgeR <- numeric(5); FPR.DEDD2.et.edgeR <- numeric(5)
FDR.DEDD2.ql.edgeR <- numeric(5); FDR.DEDD2.lr.edgeR <- numeric(5); FDR.DEDD2.et.edgeR <- numeric(5)
TPR.DEDD2.ql.edgeR <- numeric(5); TPR.DEDD2.lr.edgeR <- numeric(5); TPR.DEDD2.et.edgeR <- numeric(5)
AUC.DEDD.2.ql.edgeR <- numeric(5); AUC.DEDD.2.lr.edgeR <- numeric(5); AUC.DEDD.2.et.edgeR <- numeric(5)
FPR.DEDD2.hmm.5.expHM <- numeric(5); FPR.DEDD2.hmm.prop.expHM <- numeric(5); FPR.DEDD.2.mean.expHM <- numeric(5)
FPR.DEDD.2.log.mean.expHM <- numeric(5); FPR.DEDD.2.disp.expHM <- numeric(5); FPR.DEDD.2.log.disp.expHM <- numeric(5)
FDR.DEDD2.hmm.5.expHM <- numeric(5); FDR.DEDD2.hmm.prop.expHM <- numeric(5); FDR.DEDD2.hmm.q <- numeric(5)
FDR.DEDD.2.mean.expHM <- numeric(5); FDR.DEDD.2.log.mean.expHM <- numeric(5); FDR.DEDD.2.disp.expHM <- numeric(5)
FDR.DEDD.2.log.disp.expHM <- numeric(5); TPR.DEDD2.hmm.5.expHM <- numeric(5); TPR.DEDD2.hmm.prop.expHM <- numeric(5)
TPR.DEDD2.hmm.q <- numeric(5); TPR.DEDD.2.mean.expHM <- numeric(5); TPR.DEDD.2.log.mean.expHM <- numeric(5)
TPR.DEDD.2.disp.expHM <- numeric(5); TPR.DEDD.2.log.disp.expHM <- numeric(5); AUC.DEDD.2.hmm.expHM <- numeric(5)
AUC.DEDD.2.mean.expHM <- numeric(5); AUC.DEDD.2.log.mean.expHM <- numeric(5); AUC.DEDD.2.disp.expHM <- numeric(5)
AUC.DEDD.2.log.disp.expHM <- numeric(5); FPR.DEDD2.hmm.5.lnHM <- numeric(5); FPR.DEDD2.hmm.prop.lnHM <- numeric(5)
FPR.DEDD.2.mean.lnHM <- numeric(5); FPR.DEDD.2.log.mean.lnHM <- numeric(5); FPR.DEDD.2.disp.lnHM <- numeric(5)
FPR.DEDD.2.log.disp.lnHM <- numeric(5); FDR.DEDD2.hmm.5.lnHM <- numeric(5); FDR.DEDD2.hmm.prop.lnHM <- numeric(5)
FDR.DEDD2.hmm.q <- numeric(5); FDR.DEDD.2.mean.lnHM <- numeric(5); FDR.DEDD.2.log.mean.lnHM <- numeric(5)
FDR.DEDD.2.disp.lnHM <- numeric(5); FDR.DEDD.2.log.disp.lnHM <- numeric(5); TPR.DEDD2.hmm.5.lnHM <- numeric(5)
TPR.DEDD2.hmm.prop.lnHM <- numeric(5); TPR.DEDD2.hmm.q <- numeric(5); TPR.DEDD.2.mean.lnHM <- numeric(5)
TPR.DEDD.2.log.mean.lnHM <- numeric(5); TPR.DEDD.2.disp.lnHM <- numeric(5); TPR.DEDD.2.log.disp.lnHM <- numeric(5)
AUC.DEDD.2.hmm.lnHM <- numeric(5); AUC.DEDD.2.mean.lnHM <- numeric(5); AUC.DEDD.2.log.mean.lnHM <- numeric(5)
AUC.DEDD.2.disp.lnHM <- numeric(5); AUC.DEDD.2.log.disp.lnHM <- numeric(5)

for (i in 1:5) {
  # Generate filtered data
  counts.DEDD.2 <- simulate.DE.DD.data(dataset='DEDD.2', n.vars=2e4, samples.per.cond=samples.per.cond)
  #dim(counts.DEDD.2@count.matrix)
  #c(sum(counts.DEDD.2@variable.annotations$upregulation), sum(counts.DEDD.2@variable.annotations$downregulation), 
  #  sum(counts.DEDD.2@variable.annotations$differential.expression), 
  #  nrow(counts.DEDD.2@count.matrix)-sum(counts.DEDD.2@variable.annotations$differential.expression)) / 
  #  nrow(counts.DEDD.2@count.matrix)
  # [1] 0.05406397 0.04672497 0.10078894 0.89921106 (DE up, DE down, DE total, not DE)
  #c(sum(counts.DEDD.2@variable.annotations$updispersion), sum(counts.DEDD.2@variable.annotations$downdispersion), 
  #  sum(counts.DEDD.2@variable.annotations$differential.dispersion), 
  #  nrow(counts.DEDD.2@count.matrix)-sum(counts.DEDD.2@variable.annotations$differential.dispersion)) / 
  #  nrow(counts.DEDD.2@count.matrix)
  # [1] 0.05002752 0.04990520 0.09993273 0.90006727 (DD up, DD down, DD total, not DD)
  #c(mean(counts.DEDD.2@variable.annotations$differential.expression==0 & 
  #         counts.DEDD.2@variable.annotations$differential.dispersion==0), 
  #  mean(counts.DEDD.2@variable.annotations$differential.expression==1 & 
  #         counts.DEDD.2@variable.annotations$differential.dispersion==0), 
  #  mean(counts.DEDD.2@variable.annotations$differential.expression==0 & 
  #         counts.DEDD.2@variable.annotations$differential.dispersion==1), 
  #  mean(counts.DEDD.2@variable.annotations$differential.expression==1 & 
  #         counts.DEDD.2@variable.annotations$differential.dispersion==1))
  # [1] 0.84918354 0.05088374 0.05002752 0.04990520 (no DE/DD, DE only, DD only, both)
  # Looks like filtering hasn't affected proportions of DE/DD.
  DE.DEDD.2 <- counts.DEDD.2@variable.annotations$differential.expression
  DD.DEDD.2 <- counts.DEDD.2@variable.annotations$differential.dispersion
  #lfcm1.DEDD.2 <- abs(counts.DEDD.2@variable.annotations$truelog2foldchanges) > 1
  #lfcm2.DEDD.2 <- abs(counts.DEDD.2@variable.annotations$truelog2foldchanges) > 2
  #lfcd1.DEDD.2 <- abs(counts.DEDD.2@variable.annotations$truelog2fcdispersion) > 1
  #lfcd2.DEDD.2 <- abs(counts.DEDD.2@variable.annotations$truelog2fcdispersion) > 2
  truemeans.DEDD.2.1 <- counts.DEDD.2@variable.annotations$truemeans.S1
  truemeans.DEDD.2.2 <- counts.DEDD.2@variable.annotations$truemeans.S2
  truedisps.DEDD.2.1 <- counts.DEDD.2@variable.annotations$truedispersions.S1
  truedisps.DEDD.2.2 <- counts.DEDD.2@variable.annotations$truedispersions.S2
  
  # Normalise
  nf.DEDD.2 <- calcNormFactors(counts.DEDD.2@count.matrix)
  norm.DEDD.2 <- t(t(counts.DEDD.2@count.matrix) / nf.DEDD.2)
  
  ## edgeR
  DEDD.2.edgeR <- DGEList(counts=counts.DEDD.2@count.matrix, norm.factors=nf.DEDD.2, group=group)
  design.edgeR <- model.matrix(~group)
  DEDD.2.edgeR <- estimateDisp(DEDD.2.edgeR, design.edgeR)
  qlfit.edgeR <- glmQLFit(DEDD.2.edgeR, design.edgeR)
  qltest.edgeR <- glmQLFTest(qlfit.edgeR)
#  ql.lfc1.edgeR <- glmTreat(qlfit.edgeR, lfc=1)
#  ql.lfc2.edgeR <- glmTreat(qlfit.edgeR, lfc=2)
  lrfit.edgeR <- glmFit(DEDD.2.edgeR, design.edgeR)
  lrtest.edgeR <- glmLRT(lrfit.edgeR)
#  lr.lfc1.edgeR <- glmTreat(lrfit.edgeR, lfc=1)
#  lr.lfc2.edgeR <- glmTreat(lrfit.edgeR, lfc=2)
  et.edgeR <- exactTest(DEDD.2.edgeR)
  pval.DEDD.2.ql.edgeR <- qltest.edgeR$table$PValue
  #pval.DEDD.2.ql.lfc1.edgeR <- ql.lfc1.edgeR$table$PValue
  #pval.DEDD.2.ql.lfc2.edgeR <- ql.lfc2.edgeR$table$PValue
  pval.DEDD.2.lr.edgeR <- lrtest.edgeR$table$PValue
  #pval.DEDD.2.lr.lfc1.edgeR <- lr.lfc1.edgeR$table$PValue
  #pval.DEDD.2.lr.lfc2.edgeR <- lr.lfc2.edgeR$table$PValue
  pval.DEDD.2.et.edgeR <- et.edgeR$table$PValue
  qval.DEDD.2.ql.edgeR <- topTags(qltest.edgeR, n=nrow(DEDD.2.edgeR$counts), sort='none')$table$FDR
  #qval.DEDD.2.ql.lfc1.edgeR <- topTags(ql.lfc1.edgeR, n=nrow(DEDD.2.edgeR$counts), sort='none')$table$FDR
  #qval.DEDD.2.ql.lfc2.edgeR <- topTags(ql.lfc2.edgeR, n=nrow(DEDD.2.edgeR$counts), sort='none')$table$FDR
  qval.DEDD.2.lr.edgeR <- topTags(lrtest.edgeR, n=nrow(DEDD.2.edgeR$counts), sort='none')$table$FDR
  #qval.DEDD.2.lr.lfc1.edgeR <- topTags(lr.lfc1.edgeR, n=nrow(DEDD.2.edgeR$counts), sort='none')$table$FDR
  #qval.DEDD.2.lr.lfc2.edgeR <- topTags(lr.lfc2.edgeR, n=nrow(DEDD.2.edgeR$counts), sort='none')$table$FDR
  qval.DEDD.2.et.edgeR <- topTags(et.edgeR, n=nrow(DEDD.2.edgeR$counts), sort='none')$table$FDR
  pr.DEDD.2.ql.edgeR <- prediction(1-pval.DEDD.2.ql.edgeR, DE.DEDD.2)
  #pr.DEDD.2.ql.lfc1.edgeR <- prediction(1-pval.DEDD.2.ql.lfc1.edgeR, lfcm1.DEDD.2)
  #pr.DEDD.2.ql.lfc2.edgeR <- prediction(1-pval.DEDD.2.ql.lfc2.edgeR, lfcm2.DEDD.2)
  pr.DEDD.2.lr.edgeR <- prediction(1-pval.DEDD.2.lr.edgeR, DE.DEDD.2)
  #pr.DEDD.2.lr.lfc1.edgeR <- prediction(1-pval.DEDD.2.lr.lfc1.edgeR, lfcm1.DEDD.2)
  #pr.DEDD.2.lr.lfc2.edgeR <- prediction(1-pval.DEDD.2.lr.lfc2.edgeR, lfcm2.DEDD.2)
  pr.DEDD.2.et.edgeR <- prediction(1-pval.DEDD.2.et.edgeR, DE.DEDD.2)
  
  FPR.DEDD2.ql.edgeR[i] <- sum(pval.DEDD.2.ql.edgeR<0.05 & DE.DEDD.2==0) / sum(DE.DEDD.2==0)
  #FPR.DEDD2.ql.lfc1.edgeR <- sum(pval.DEDD.2.ql.lfc1.edgeR<0.05 & lfcm1.DEDD.2==0) / sum(lfcm1.DEDD.2==0)
  #FPR.DEDD2.ql.lfc2.edgeR <- sum(pval.DEDD.2.ql.lfc2.edgeR<0.05 & lfcm2.DEDD.2==0) / sum(lfcm2.DEDD.2==0)
  FPR.DEDD2.lr.edgeR[i] <- sum(pval.DEDD.2.lr.edgeR<0.05 & DE.DEDD.2==0) / sum(DE.DEDD.2==0)
  #FPR.DEDD2.lr.lfc1.edgeR <- sum(pval.DEDD.2.lr.lfc1.edgeR<0.05 & lfcm1.DEDD.2==0) / sum(lfcm1.DEDD.2==0)
  #FPR.DEDD2.lr.lfc2.edgeR <- sum(pval.DEDD.2.lr.lfc2.edgeR<0.05 & lfcm2.DEDD.2==0) / sum(lfcm2.DEDD.2==0)
  FPR.DEDD2.et.edgeR[i] <- sum(pval.DEDD.2.et.edgeR<0.05 & DE.DEDD.2==0) / sum(DE.DEDD.2==0)
  FDR.DEDD2.ql.edgeR[i] <- sum(qval.DEDD.2.ql.edgeR<0.05 & DE.DEDD.2==0) / sum(qval.DEDD.2.ql.edgeR<0.05)
  #FDR.DEDD2.ql.lfc1.edgeR <- sum(qval.DEDD.2.ql.lfc1.edgeR<0.05 & lfcm1.DEDD.2==0) / sum(qval.DEDD.2.ql.lfc1.edgeR<0.05)
  #FDR.DEDD2.ql.lfc2.edgeR <- sum(qval.DEDD.2.ql.lfc2.edgeR<0.05 & lfcm2.DEDD.2==0) / sum(qval.DEDD.2.ql.lfc2.edgeR<0.05)
  FDR.DEDD2.lr.edgeR[i] <- sum(qval.DEDD.2.lr.edgeR<0.05 & DE.DEDD.2==0) / sum(qval.DEDD.2.lr.edgeR<0.05)
  #FDR.DEDD2.lr.lfc1.edgeR <- sum(qval.DEDD.2.lr.lfc1.edgeR<0.05 & lfcm1.DEDD.2==0) / sum(qval.DEDD.2.lr.lfc1.edgeR<0.05)
  #FDR.DEDD2.lr.lfc2.edgeR <- sum(qval.DEDD.2.lr.lfc2.edgeR<0.05 & lfcm2.DEDD.2==0) / sum(qval.DEDD.2.lr.lfc2.edgeR<0.05)
  FDR.DEDD2.et.edgeR[i] <- sum(qval.DEDD.2.et.edgeR<0.05 & DE.DEDD.2==0) / sum(qval.DEDD.2.et.edgeR<0.05)
  TPR.DEDD2.ql.edgeR[i] <- sum(qval.DEDD.2.ql.edgeR<0.05 & DE.DEDD.2==1) / sum(DE.DEDD.2==1)
  #TPR.DEDD2.ql.lfc1.edgeR <- sum(qval.DEDD.2.ql.lfc1.edgeR<0.05 & lfcm1.DEDD.2==1) / sum(lfcm1.DEDD.2==1)
  #TPR.DEDD2.ql.lfc2.edgeR <- sum(qval.DEDD.2.ql.lfc2.edgeR<0.05 & lfcm2.DEDD.2==1) / sum(lfcm2.DEDD.2==1)
  TPR.DEDD2.lr.edgeR[i] <- sum(qval.DEDD.2.lr.edgeR<0.05 & DE.DEDD.2==1) / sum(DE.DEDD.2==1)
  #TPR.DEDD2.lr.lfc1.edgeR <- sum(qval.DEDD.2.lr.lfc1.edgeR<0.05 & lfcm1.DEDD.2==1) / sum(lfcm1.DEDD.2==1)
  #TPR.DEDD2.lr.lfc2.edgeR <- sum(qval.DEDD.2.lr.lfc2.edgeR<0.05 & lfcm2.DEDD.2==1) / sum(lfcm2.DEDD.2==1)
  TPR.DEDD2.et.edgeR[i] <- sum(qval.DEDD.2.et.edgeR<0.05 & DE.DEDD.2==1) / sum(DE.DEDD.2==1)
  AUC.DEDD.2.ql.edgeR[i] <- performance(pr.DEDD.2.ql.edgeR, measure='auc')@y.values[[1]]
  #AUC.DEDD.2.ql.lfc1.edgeR <- performance(pr.DEDD.2.ql.lfc1.edgeR, measure='auc')@y.values[[1]]
  #AUC.DEDD.2.ql.lfc2.edgeR <- performance(pr.DEDD.2.ql.lfc2.edgeR, measure='auc')@y.values[[1]]
  AUC.DEDD.2.lr.edgeR[i] <- performance(pr.DEDD.2.lr.edgeR, measure='auc')@y.values[[1]]
  #AUC.DEDD.2.lr.lfc1.edgeR <- performance(pr.DEDD.2.lr.lfc1.edgeR, measure='auc')@y.values[[1]]
  #AUC.DEDD.2.lr.lfc2.edgeR <- performance(pr.DEDD.2.lr.lfc2.edgeR, measure='auc')@y.values[[1]]
  AUC.DEDD.2.et.edgeR[i] <- performance(pr.DEDD.2.et.edgeR, measure='auc')@y.values[[1]]
#  MSE.m.lr.edgeR <- mean((rowMeans(lrfit.edgeR$fitted.values) - truemeans.DEDD.2.1)[which(DE.DEDD.2==0)]^2)
#  MSE.m.ql.edgeR <- mean((rowMeans(qlfit.edgeR$fitted.values) - truemeans.DEDD.2.1)[which(DE.DEDD.2==0)]^2)
#  MSE.d.trend.edgeR <- mean((DEDD.2.edgeR$trended.dispersion - truedisps.DEDD.2.1)[which(DD.DEDD.2==0)]^2)
#  MSE.d.tag.edgeR <- mean((DEDD.2.edgeR$tagwise.dispersion - truedisps.DEDD.2.1)[which(DD.DEDD.2==0)]^2)
  
  
  ## DESeq2
  
  
  ## limma-voom
  
  
  ## baySeq
  
  
  ## ShrinkSeq
  
  
  ## MDSeq
  
  
  ## expHM
  DEDD.2.expHM <- exp_hmm_adapt_3_chains(counts=t(norm.DEDD.2), groups=group)
  prob.DEDD.2.expHM <- colMeans(as.matrix(DEDD.2.expHM$indicators))
  post.prop.DEDD.2.expHM <- mean(as.matrix(DEDD.2.expHM$proportion))
  thr.DEDD.2.expHM <- sort(prob.DEDD.2.expHM, decreasing=T)[round(nrow(counts.DEDD.2@count.matrix)*post.prop.DEDD.2.expHM)]
  mean.diff.DEDD.2.expHM <- as.matrix(DEDD.2.expHM$means1) - as.matrix(DEDD.2.expHM$means2)
  log.mean.diff.DEDD.2.expHM <- log(as.matrix(DEDD.2.expHM$means1)) - log(as.matrix(DEDD.2.expHM$means2))
  disp.diff.DEDD.2.expHM <- as.matrix(DEDD.2.expHM$disps1) - as.matrix(DEDD.2.expHM$disps2)
  log.disp.diff.DEDD.2.expHM <- log(as.matrix(DEDD.2.expHM$disps1)) - log(as.matrix(DEDD.2.expHM$disps2))
  p.mean.DEDD.2.expHM <- apply(mean.diff.DEDD.2.expHM,2,hpd.pval)
  p.log.mean.DEDD.2.expHM <- apply(log.mean.diff.DEDD.2.expHM,2,hpd.pval)
  p.disp.DEDD.2.expHM <- apply(disp.diff.DEDD.2.expHM,2,hpd.pval)
  p.log.disp.DEDD.2.expHM <- apply(log.disp.diff.DEDD.2.expHM,2,hpd.pval)
  q.hmm.DEDD.2.expHM <- bfdr(prob.DEDD.2.expHM)
  q.mean.DEDD.2.expHM <- qvalue(p.mean.DEDD.2.expHM)$qvalues
  q.log.mean.DEDD.2.expHM <- qvalue(p.log.mean.DEDD.2.expHM)$qvalues
  q.disp.DEDD.2.expHM <- qvalue(p.disp.DEDD.2.expHM)$qvalues
  q.log.disp.DEDD.2.expHM <- qvalue(p.log.disp.DEDD.2.expHM)$qvalues
  pr.DEDD.2.hmm.expHM <- prediction(prob.DEDD.2.expHM, (DE.DEDD.2 | DD.DEDD.2))
  pr.DEDD.2.mean.expHM <- prediction(1-p.mean.DEDD.2.expHM, DE.DEDD.2)
  pr.DEDD.2.log.mean.expHM <- prediction(1-p.log.mean.DEDD.2.expHM, DE.DEDD.2)
  pr.DEDD.2.disp.expHM <- prediction(1-p.disp.DEDD.2.expHM, DD.DEDD.2)
  pr.DEDD.2.log.disp.expHM <- prediction(1-p.log.disp.DEDD.2.expHM, DD.DEDD.2)
  
  FPR.DEDD2.hmm.5.expHM[i] <- sum(prob.DEDD.2.expHM>0.5 & DE.DEDD.2==0 & DD.DEDD.2==0) / 
    sum(DE.DEDD.2==0 & DD.DEDD.2==0)
  FPR.DEDD2.hmm.prop.expHM[i] <- sum(prob.DEDD.2.expHM>thr.DEDD.2.expHM & DE.DEDD.2==0 & DD.DEDD.2==0) / 
    sum(DE.DEDD.2==0 & DD.DEDD.2==0)
  FPR.DEDD.2.mean.expHM[i] <- sum(p.mean.DEDD.2.expHM<0.05 & DE.DEDD.2==0) / sum(DE.DEDD.2==0)
  FPR.DEDD.2.log.mean.expHM[i] <- sum(p.log.mean.DEDD.2.expHM<0.05 & DE.DEDD.2==0) / sum(DE.DEDD.2==0)
  FPR.DEDD.2.disp.expHM[i] <- sum(p.disp.DEDD.2.expHM<0.05 & DD.DEDD.2==0) / sum(DD.DEDD.2==0)
  FPR.DEDD.2.log.disp.expHM[i] <- sum(p.log.disp.DEDD.2.expHM<0.05 & DD.DEDD.2==0) / sum(DD.DEDD.2==0)
  FDR.DEDD2.hmm.5.expHM[i] <- sum(prob.DEDD.2.expHM>0.5 & DE.DEDD.2==0 & DD.DEDD.2==0) / sum(prob.DEDD.2.expHM>0.5)
  FDR.DEDD2.hmm.prop.expHM[i] <- sum(prob.DEDD.2.expHM>thr.DEDD.2.expHM & DE.DEDD.2==0 & DD.DEDD.2==0) / 
    sum(prob.DEDD.2.expHM>thr.DEDD.2.expHM)
  FDR.DEDD2.hmm.q[i] <- sum(q.hmm.DEDD.2.expHM<0.05 & DE.DEDD.2==0 & DD.DEDD.2==0) / sum(q.hmm.DEDD.2.expHM<0.05)
  FDR.DEDD.2.mean.expHM[i] <- sum(q.mean.DEDD.2.expHM<0.05 & DE.DEDD.2==0) / sum(q.mean.DEDD.2.expHM<0.05)
  FDR.DEDD.2.log.mean.expHM[i] <- sum(q.log.mean.DEDD.2.expHM<0.05 & DE.DEDD.2==0) / sum(q.log.mean.DEDD.2.expHM<0.05)
  FDR.DEDD.2.disp.expHM[i] <- sum(q.disp.DEDD.2.expHM<0.05 & DD.DEDD.2==0) / sum(q.disp.DEDD.2.expHM<0.05)
  FDR.DEDD.2.log.disp.expHM[i] <- sum(q.log.disp.DEDD.2.expHM<0.05 & DD.DEDD.2==0) / sum(q.log.disp.DEDD.2.expHM<0.05)
  TPR.DEDD2.hmm.5.expHM[i] <- sum(prob.DEDD.2.expHM>0.5 & (DE.DEDD.2==1 | DD.DEDD.2==1)) / 
    sum(DE.DEDD.2==1 | DD.DEDD.2==1)
  TPR.DEDD2.hmm.prop.expHM[i] <- sum(prob.DEDD.2.expHM>thr.DEDD.2.expHM & (DE.DEDD.2==1 | DD.DEDD.2==1)) / 
    sum(DE.DEDD.2==1 | DD.DEDD.2==1)
  TPR.DEDD2.hmm.q[i] <- sum(q.hmm.DEDD.2.expHM<0.05 & (DE.DEDD.2==1 | DD.DEDD.2==1)) / sum(DE.DEDD.2==1 | DD.DEDD.2==1)
  TPR.DEDD.2.mean.expHM[i] <- sum(q.mean.DEDD.2.expHM<0.05 & DE.DEDD.2==1) / sum(DE.DEDD.2==1)
  TPR.DEDD.2.log.mean.expHM[i] <- sum(q.log.mean.DEDD.2.expHM<0.05 & DE.DEDD.2==1) / sum(DE.DEDD.2==1)
  TPR.DEDD.2.disp.expHM[i] <- sum(q.disp.DEDD.2.expHM<0.05 & DD.DEDD.2==1) / sum(DD.DEDD.2==1)
  TPR.DEDD.2.log.disp.expHM[i] <- sum(q.log.disp.DEDD.2.expHM<0.05 & DD.DEDD.2==1) / sum(DD.DEDD.2==1)
  AUC.DEDD.2.hmm.expHM[i] <- performance(pr.DEDD.2.hmm.expHM, measure='auc')@y.values[[1]]
  AUC.DEDD.2.mean.expHM[i] <- performance(pr.DEDD.2.mean.expHM, measure='auc')@y.values[[1]]
  AUC.DEDD.2.log.mean.expHM[i] <- performance(pr.DEDD.2.log.mean.expHM, measure='auc')@y.values[[1]]
  AUC.DEDD.2.disp.expHM[i] <- performance(pr.DEDD.2.disp.expHM, measure='auc')@y.values[[1]]
  AUC.DEDD.2.log.disp.expHM[i] <- performance(pr.DEDD.2.log.disp.expHM, measure='auc')@y.values[[1]]
#  MSE.m.expHM <- mean((colMeans(as.matrix(DEDD.2.expHM$means0)) - truemeans.DEDD.2.1)[which(DE.DEDD.2==0)]^2)
#  MSE.d.expHM <- mean((colMeans(as.matrix(DEDD.2.expHM$disps0)) - truedisps.DEDD.2.1)[which(DD.DEDD.2==0)]^2)
  
  
  ## lnHM
  DEDD.2.lnHM <- ln_hmm_adapt_3_chains(counts=t(norm.DEDD.2), groups=group)
  prob.DEDD.2.lnHM <- colMeans(as.matrix(DEDD.2.lnHM$indicators))
  post.prop.DEDD.2.lnHM <- mean(as.matrix(DEDD.2.lnHM$proportion))
  thr.DEDD.2.lnHM <- sort(prob.DEDD.2.lnHM, decreasing=T)[round(nrow(counts.DEDD.2@count.matrix)*post.prop.DEDD.2.lnHM)]
  mean.diff.DEDD.2.lnHM <- as.matrix(DEDD.2.lnHM$means1) - as.matrix(DEDD.2.lnHM$means2)
  log.mean.diff.DEDD.2.lnHM <- log(as.matrix(DEDD.2.lnHM$means1)) - log(as.matrix(DEDD.2.lnHM$means2))
  disp.diff.DEDD.2.lnHM <- as.matrix(DEDD.2.lnHM$disps1) - as.matrix(DEDD.2.lnHM$disps2)
  log.disp.diff.DEDD.2.lnHM <- log(as.matrix(DEDD.2.lnHM$disps1)) - log(as.matrix(DEDD.2.lnHM$disps2))
  p.mean.DEDD.2.lnHM <- apply(mean.diff.DEDD.2.lnHM,2,hpd.pval)
  p.log.mean.DEDD.2.lnHM <- apply(log.mean.diff.DEDD.2.lnHM,2,hpd.pval)
  p.disp.DEDD.2.lnHM <- apply(disp.diff.DEDD.2.lnHM,2,hpd.pval)
  p.log.disp.DEDD.2.lnHM <- apply(log.disp.diff.DEDD.2.lnHM,2,hpd.pval)
  q.hmm.DEDD.2.lnHM <- bfdr(prob.DEDD.2.lnHM)
  q.mean.DEDD.2.lnHM <- qvalue(p.mean.DEDD.2.lnHM)$qvalues
  q.log.mean.DEDD.2.lnHM <- qvalue(p.log.mean.DEDD.2.lnHM)$qvalues
  q.disp.DEDD.2.lnHM <- qvalue(p.disp.DEDD.2.lnHM)$qvalues
  q.log.disp.DEDD.2.lnHM <- qvalue(p.log.disp.DEDD.2.lnHM)$qvalues
  pr.DEDD.2.hmm.lnHM <- prediction(prob.DEDD.2.lnHM, (DE.DEDD.2 | DD.DEDD.2))
  pr.DEDD.2.mean.lnHM <- prediction(1-p.mean.DEDD.2.lnHM, DE.DEDD.2)
  pr.DEDD.2.log.mean.lnHM <- prediction(1-p.log.mean.DEDD.2.lnHM, DE.DEDD.2)
  pr.DEDD.2.disp.lnHM <- prediction(1-p.disp.DEDD.2.lnHM, DD.DEDD.2)
  pr.DEDD.2.log.disp.lnHM <- prediction(1-p.log.disp.DEDD.2.lnHM, DD.DEDD.2)
  
  FPR.DEDD2.hmm.5.lnHM[i] <- sum(prob.DEDD.2.lnHM>0.5 & DE.DEDD.2==0 & DD.DEDD.2==0) / sum(DE.DEDD.2==0 & DD.DEDD.2==0)
  FPR.DEDD2.hmm.prop.lnHM[i] <- sum(prob.DEDD.2.lnHM>thr.DEDD.2.lnHM & DE.DEDD.2==0 & DD.DEDD.2==0) / 
    sum(DE.DEDD.2==0 & DD.DEDD.2==0)
  FPR.DEDD.2.mean.lnHM[i] <- sum(p.mean.DEDD.2.lnHM<0.05 & DE.DEDD.2==0) / sum(DE.DEDD.2==0)
  FPR.DEDD.2.log.mean.lnHM[i] <- sum(p.log.mean.DEDD.2.lnHM<0.05 & DE.DEDD.2==0) / sum(DE.DEDD.2==0)
  FPR.DEDD.2.disp.lnHM[i] <- sum(p.disp.DEDD.2.lnHM<0.05 & DD.DEDD.2==0) / sum(DD.DEDD.2==0)
  FPR.DEDD.2.log.disp.lnHM[i] <- sum(p.log.disp.DEDD.2.lnHM<0.05 & DD.DEDD.2==0) / sum(DD.DEDD.2==0)
  FDR.DEDD2.hmm.5.lnHM[i] <- sum(prob.DEDD.2.lnHM>0.5 & DE.DEDD.2==0 & DD.DEDD.2==0) / sum(prob.DEDD.2.lnHM>0.5)
  FDR.DEDD2.hmm.prop.lnHM[i] <- sum(prob.DEDD.2.lnHM>thr.DEDD.2.lnHM & DE.DEDD.2==0 & DD.DEDD.2==0) / 
    sum(prob.DEDD.2.lnHM>thr.DEDD.2.lnHM)
  FDR.DEDD2.hmm.q[i] <- sum(q.hmm.DEDD.2.lnHM<0.05 & DE.DEDD.2==0 & DD.DEDD.2==0) / sum(q.hmm.DEDD.2.lnHM<0.05)
  FDR.DEDD.2.mean.lnHM[i] <- sum(q.mean.DEDD.2.lnHM<0.05 & DE.DEDD.2==0) / sum(q.mean.DEDD.2.lnHM<0.05)
  FDR.DEDD.2.log.mean.lnHM[i] <- sum(q.log.mean.DEDD.2.lnHM<0.05 & DE.DEDD.2==0) / sum(q.log.mean.DEDD.2.lnHM<0.05)
  FDR.DEDD.2.disp.lnHM[i] <- sum(q.disp.DEDD.2.lnHM<0.05 & DD.DEDD.2==0) / sum(q.disp.DEDD.2.lnHM<0.05)
  FDR.DEDD.2.log.disp.lnHM[i] <- sum(q.log.disp.DEDD.2.lnHM<0.05 & DD.DEDD.2==0) / sum(q.log.disp.DEDD.2.lnHM<0.05)
  TPR.DEDD2.hmm.5.lnHM[i] <- sum(prob.DEDD.2.lnHM>0.5 & (DE.DEDD.2==1 | DD.DEDD.2==1)) / sum(DE.DEDD.2==1 | DD.DEDD.2==1)
  TPR.DEDD2.hmm.prop.lnHM[i] <- sum(prob.DEDD.2.lnHM>thr.DEDD.2.lnHM & (DE.DEDD.2==1 | DD.DEDD.2==1)) / 
    sum(DE.DEDD.2==1 | DD.DEDD.2==1)
  TPR.DEDD2.hmm.q[i] <- sum(q.hmm.DEDD.2.lnHM<0.05 & (DE.DEDD.2==1 | DD.DEDD.2==1)) / sum(DE.DEDD.2==1 | DD.DEDD.2==1)
  TPR.DEDD.2.mean.lnHM[i] <- sum(q.mean.DEDD.2.lnHM<0.05 & DE.DEDD.2==1) / sum(DE.DEDD.2==1)
  TPR.DEDD.2.log.mean.lnHM[i] <- sum(q.log.mean.DEDD.2.lnHM<0.05 & DE.DEDD.2==1) / sum(DE.DEDD.2==1)
  TPR.DEDD.2.disp.lnHM[i] <- sum(q.disp.DEDD.2.lnHM<0.05 & DD.DEDD.2==1) / sum(DD.DEDD.2==1)
  TPR.DEDD.2.log.disp.lnHM[i] <- sum(q.log.disp.DEDD.2.lnHM<0.05 & DD.DEDD.2==1) / sum(DD.DEDD.2==1)
  AUC.DEDD.2.hmm.lnHM[i] <- performance(pr.DEDD.2.hmm.lnHM, measure='auc')@y.values[[1]]
  AUC.DEDD.2.mean.lnHM[i] <- performance(pr.DEDD.2.mean.lnHM, measure='auc')@y.values[[1]]
  AUC.DEDD.2.log.mean.lnHM[i] <- performance(pr.DEDD.2.log.mean.lnHM, measure='auc')@y.values[[1]]
  AUC.DEDD.2.disp.lnHM[i] <- performance(pr.DEDD.2.disp.lnHM, measure='auc')@y.values[[1]]
  AUC.DEDD.2.log.disp.lnHM[i] <- performance(pr.DEDD.2.log.disp.lnHM, measure='auc')@y.values[[1]]
#  MSE.m.lnHM <- mean((colMeans(as.matrix(DEDD.2.lnHM$means0)) - truemeans.DEDD.2.1)[which(DE.DEDD.2==0)]^2)
#  MSE.d.lnHM <- mean((colMeans(as.matrix(DEDD.2.lnHM$disps0)) - truedisps.DEDD.2.1)[which(DD.DEDD.2==0)]^2)
}

results <- list('FPR.DEDD2.ql.edgeR' = FPR.DEDD2.ql.edgeR, 
                'FPR.DEDD2.lr.edgeR' = FPR.DEDD2.lr.edgeR, 
                'FPR.DEDD2.et.edgeR' = FPR.DEDD2.et.edgeR, 
                'FDR.DEDD2.ql.edgeR' = FDR.DEDD2.ql.edgeR, 
                'FDR.DEDD2.lr.edgeR' = FDR.DEDD2.lr.edgeR, 
                'FDR.DEDD2.et.edgeR' = FDR.DEDD2.et.edgeR, 
                'TPR.DEDD2.ql.edgeR' = TPR.DEDD2.ql.edgeR, 
                'TPR.DEDD2.lr.edgeR' = TPR.DEDD2.lr.edgeR, 
                'TPR.DEDD2.et.edgeR' = TPR.DEDD2.et.edgeR, 
                'AUC.DEDD.2.ql.edgeR' = AUC.DEDD.2.ql.edgeR, 
                'AUC.DEDD.2.lr.edgeR' = AUC.DEDD.2.lr.edgeR, 
                'AUC.DEDD.2.et.edgeR' = AUC.DEDD.2.et.edgeR, 
                'FPR.DEDD2.hmm.5.expHM' = FPR.DEDD2.hmm.5.expHM, 
                'FPR.DEDD2.hmm.prop.expHM' = FPR.DEDD2.hmm.prop.expHM, 
                'FPR.DEDD.2.mean.expHM' = FPR.DEDD.2.mean.expHM, 
                'FPR.DEDD.2.log.mean.expHM' = FPR.DEDD.2.log.mean.expHM, 
                'FPR.DEDD.2.disp.expHM' = FPR.DEDD.2.disp.expHM, 
                'FPR.DEDD.2.log.disp.expHM' = FPR.DEDD.2.log.disp.expHM, 
                'FDR.DEDD2.hmm.5.expHM' = FDR.DEDD2.hmm.5.expHM, 
                'FDR.DEDD2.hmm.prop.expHM' = FDR.DEDD2.hmm.prop.expHM, 
                'FDR.DEDD2.hmm.q' = FDR.DEDD2.hmm.q, 
                'FDR.DEDD.2.mean.expHM' = FDR.DEDD.2.mean.expHM, 
                'FDR.DEDD.2.log.mean.expHM' = FDR.DEDD.2.log.mean.expHM, 
                'FDR.DEDD.2.disp.expHM' = FDR.DEDD.2.disp.expHM, 
                'FDR.DEDD.2.log.disp.expHM' = FDR.DEDD.2.log.disp.expHM, 
                'TPR.DEDD2.hmm.5.expHM' = TPR.DEDD2.hmm.5.expHM, 
                'TPR.DEDD2.hmm.prop.expHM' = TPR.DEDD2.hmm.prop.expHM, 
                'TPR.DEDD2.hmm.q' = TPR.DEDD2.hmm.q, 
                'TPR.DEDD.2.mean.expHM' = TPR.DEDD.2.mean.expHM, 
                'TPR.DEDD.2.log.mean.expHM' = TPR.DEDD.2.log.mean.expHM, 
                'TPR.DEDD.2.disp.expHM' = TPR.DEDD.2.disp.expHM, 
                'TPR.DEDD.2.log.disp.expHM' = TPR.DEDD.2.log.disp.expHM, 
                'AUC.DEDD.2.hmm.expHM' = AUC.DEDD.2.hmm.expHM, 
                'AUC.DEDD.2.mean.expHM' = AUC.DEDD.2.mean.expHM, 
                'AUC.DEDD.2.log.mean.expHM' = AUC.DEDD.2.log.mean.expHM, 
                'AUC.DEDD.2.disp.expHM' = AUC.DEDD.2.disp.expHM, 
                'AUC.DEDD.2.log.disp.expHM' = AUC.DEDD.2.log.disp.expHM, 
                'FPR.DEDD2.hmm.5.lnHM' = FPR.DEDD2.hmm.5.lnHM, 
                'FPR.DEDD2.hmm.prop.lnHM' = FPR.DEDD2.hmm.prop.lnHM, 
                'FPR.DEDD.2.mean.lnHM' = FPR.DEDD.2.mean.lnHM, 
                'FPR.DEDD.2.log.mean.lnHM' = FPR.DEDD.2.log.mean.lnHM, 
                'FPR.DEDD.2.disp.lnHM' = FPR.DEDD.2.disp.lnHM, 
                'FPR.DEDD.2.log.disp.lnHM' = FPR.DEDD.2.log.disp.lnHM, 
                'FDR.DEDD2.hmm.5.lnHM' = FDR.DEDD2.hmm.5.lnHM, 
                'FDR.DEDD2.hmm.prop.lnHM' = FDR.DEDD2.hmm.prop.lnHM, 
                'FDR.DEDD2.hmm.q' = FDR.DEDD2.hmm.q, 
                'FDR.DEDD.2.mean.lnHM' = FDR.DEDD.2.mean.lnHM, 
                'FDR.DEDD.2.log.mean.lnHM' = FDR.DEDD.2.log.mean.lnHM, 
                'FDR.DEDD.2.disp.lnHM' = FDR.DEDD.2.disp.lnHM, 
                'FDR.DEDD.2.log.disp.lnHM' = FDR.DEDD.2.log.disp.lnHM, 
                'TPR.DEDD2.hmm.5.lnHM' = TPR.DEDD2.hmm.5.lnHM, 
                'TPR.DEDD2.hmm.prop.lnHM' = TPR.DEDD2.hmm.prop.lnHM, 
                'TPR.DEDD2.hmm.q' = TPR.DEDD2.hmm.q, 
                'TPR.DEDD.2.mean.lnHM' = TPR.DEDD.2.mean.lnHM, 
                'TPR.DEDD.2.log.mean.lnHM' = TPR.DEDD.2.log.mean.lnHM, 
                'TPR.DEDD.2.disp.lnHM' = TPR.DEDD.2.disp.lnHM, 
                'TPR.DEDD.2.log.disp.lnHM' = TPR.DEDD.2.log.disp.lnHM, 
                'AUC.DEDD.2.hmm.lnHM' = AUC.DEDD.2.hmm.lnHM, 
                'AUC.DEDD.2.mean.lnHM' = AUC.DEDD.2.mean.lnHM, 
                'AUC.DEDD.2.log.mean.lnHM' = AUC.DEDD.2.log.mean.lnHM, 
                'AUC.DEDD.2.disp.lnHM' = AUC.DEDD.2.disp.lnHM, 
                'AUC.DEDD.2.log.disp.lnHM' = AUC.DEDD.2.log.disp.lnHM)
filename <- paste0('results.DEDD.5.', format(Sys.time(), "%Y-%m-%d_%H%M%S"), '.rds')
saveRDS(results, file=here(filename))








#bq.mean.DEDD.2.lnHM <- bfdr(1-p.mean.DEDD.2.lnHM)
#bq.log.mean.DEDD.2.lnHM <- bfdr(1-p.log.mean.DEDD.2.lnHM)
#bq.disp.DEDD.2.lnHM <- bfdr(1-p.disp.DEDD.2.lnHM)
#bq.log.disp.DEDD.2.lnHM <- bfdr(1-p.log.disp.DEDD.2.lnHM)
#q.bq.mean.DEDD.2.lnHM <- qvalue(bq.mean.DEDD.2.lnHM, pi0=1)$qvalues
#q.bq.log.mean.DEDD.2.lnHM <- qvalue(bq.log.mean.DEDD.2.lnHM, pi0=1)$qvalues
#q.bq.disp.DEDD.2.lnHM <- qvalue(bq.disp.DEDD.2.lnHM, pi0=1)$qvalues
#q.bq.log.disp.DEDD.2.lnHM <- qvalue(bq.log.disp.DEDD.2.lnHM, pi0=1)$qvalues
#qbFD.DEDD.2.mean.lnHM <- sum(q.bq.mean.DEDD.2.lnHM<0.05 & DE.DEDD.2==0) / sum(q.bq.mean.DEDD.2.lnHM<0.05)
#qbFD.DEDD.2.log.mean.lnHM <- sum(q.bq.log.mean.DEDD.2.lnHM<0.05 & DE.DEDD.2==0) / sum(q.bq.log.mean.DEDD.2.lnHM<0.05)
#qbFD.DEDD.2.disp.lnHM <- sum(q.bq.disp.DEDD.2.lnHM<0.05 & DD.DEDD.2==0) / sum(q.bq.disp.DEDD.2.lnHM<0.05)
#qbFD.DEDD.2.log.disp.lnHM <- sum(q.bq.log.disp.DEDD.2.lnHM<0.05 & DD.DEDD.2==0) / sum(q.bq.log.disp.DEDD.2.lnHM<0.05)
#TPR.qb.DEDD.2.mean.lnHM <- sum(q.bq.mean.DEDD.2.lnHM<0.05 & DE.DEDD.2==1) / sum(DE.DEDD.2==1)
#TPR.qb.DEDD.2.log.mean.lnHM <- sum(q.bq.log.mean.DEDD.2.lnHM<0.05 & DE.DEDD.2==1) / sum(DE.DEDD.2==1)
#TPR.qb.DEDD.2.disp.lnHM <- sum(q.bq.disp.DEDD.2.lnHM<0.05 & DD.DEDD.2==1) / sum(DD.DEDD.2==1)
#TPR.qb.DEDD.2.log.disp.lnHM <- sum(q.bq.log.disp.DEDD.2.lnHM<0.05 & DD.DEDD.2==1) / sum(DD.DEDD.2==1)

#FDR.p.DEDD.2.mean.lnHM <- sum(p.mean.DEDD.2.lnHM<0.05 & DE.DEDD.2==0) / sum(p.mean.DEDD.2.lnHM<0.05)
#FDR.p.DEDD.2.log.mean.lnHM <- sum(p.log.mean.DEDD.2.lnHM<0.05 & DE.DEDD.2==0) / sum(p.log.mean.DEDD.2.lnHM<0.05)
#FDR.p.DEDD.2.disp.lnHM <- sum(p.disp.DEDD.2.lnHM<0.05 & DD.DEDD.2==0) / sum(p.disp.DEDD.2.lnHM<0.05)
#FDR.p.DEDD.2.log.disp.lnHM <- sum(p.log.disp.DEDD.2.lnHM<0.05 & DD.DEDD.2==0) / sum(p.log.disp.DEDD.2.lnHM<0.05)
#TPR.p.DEDD.2.mean.lnHM <- sum(p.mean.DEDD.2.lnHM<0.05 & DE.DEDD.2==1) / sum(DE.DEDD.2==1)
#TPR.p.DEDD.2.log.mean.lnHM <- sum(p.log.mean.DEDD.2.lnHM<0.05 & DE.DEDD.2==1) / sum(DE.DEDD.2==1)
#TPR.p.DEDD.2.disp.lnHM <- sum(p.disp.DEDD.2.lnHM<0.05 & DD.DEDD.2==1) / sum(DD.DEDD.2==1)
#TPR.p.DEDD.2.log.disp.lnHM <- sum(p.log.disp.DEDD.2.lnHM<0.05 & DD.DEDD.2==1) / sum(DD.DEDD.2==1)

#rbind(c(FPR.DEDD.2.mean.lnHM, FPR.DEDD.2.log.mean.lnHM, FPR.DEDD.2.disp.lnHM, FPR.DEDD.2.log.disp.lnHM), 
#      c(FDR.DEDD.2.mean.lnHM, FDR.DEDD.2.log.mean.lnHM, FDR.DEDD.2.disp.lnHM, FDR.DEDD.2.log.disp.lnHM), 
#      c(FDR.p.DEDD.2.mean.lnHM, FDR.p.DEDD.2.log.mean.lnHM, FDR.p.DEDD.2.disp.lnHM, FDR.p.DEDD.2.log.disp.lnHM), 
#      c(qbFD.DEDD.2.mean.lnHM, qbFD.DEDD.2.log.mean.lnHM, qbFD.DEDD.2.disp.lnHM, qbFD.DEDD.2.log.disp.lnHM), 
#      c(TPR.DEDD.2.mean.lnHM, TPR.DEDD.2.log.mean.lnHM, TPR.DEDD.2.disp.lnHM, TPR.DEDD.2.log.disp.lnHM), 
#      c(TPR.p.DEDD.2.mean.lnHM, TPR.p.DEDD.2.log.mean.lnHM, TPR.p.DEDD.2.disp.lnHM, TPR.p.DEDD.2.log.disp.lnHM), 
#      c(TPR.qb.DEDD.2.mean.lnHM, TPR.qb.DEDD.2.log.mean.lnHM, TPR.qb.DEDD.2.disp.lnHM, TPR.qb.DEDD.2.log.disp.lnHM), 
#      c(AUC.DEDD.2.mean.lnHM, AUC.DEDD.2.log.mean.lnHM, AUC.DEDD.2.disp.lnHM, AUC.DEDD.2.log.disp.lnHM))




#disc.mean.DEDD.2.expHM <- numeric(1000); disc.log.mean.DEDD.2.expHM <- numeric(1000)
#disc.disp.DEDD.2.expHM <- numeric(1000); disc.log.disp.DEDD.2.expHM <- numeric(1000)
#disc.mean.DEDD.2.lnHM <- numeric(1000); disc.log.mean.DEDD.2.lnHM <- numeric(1000)
#disc.disp.DEDD.2.lnHM <- numeric(1000); disc.log.disp.DEDD.2.lnHM <- numeric(1000)
#fd.mean.DEDD.2.expHM <- numeric(1000); fd.log.mean.DEDD.2.expHM <- numeric(1000)
#fd.disp.DEDD.2.expHM <- numeric(1000); fd.log.disp.DEDD.2.expHM <- numeric(1000)
#fd.mean.DEDD.2.lnHM <- numeric(1000); fd.log.mean.DEDD.2.lnHM <- numeric(1000)
#fd.disp.DEDD.2.lnHM <- numeric(1000); fd.log.disp.DEDD.2.lnHM <- numeric(1000)
#for (t in 1:1000) {
#  disc.mean.DEDD.2.expHM[t] <- sum(p.mean.DEDD.2.expHM < t/1000)
#  disc.log.mean.DEDD.2.expHM[t] <- sum(p.log.mean.DEDD.2.expHM < t/1000)
#  disc.disp.DEDD.2.expHM[t] <- sum(p.disp.DEDD.2.expHM < t/1000)
#  disc.log.disp.DEDD.2.expHM[t] <- sum(p.log.disp.DEDD.2.expHM < t/1000)
#  disc.mean.DEDD.2.lnHM[t] <- sum(p.mean.DEDD.2.lnHM < t/1000)
#  disc.log.mean.DEDD.2.lnHM[t] <- sum(p.log.mean.DEDD.2.lnHM < t/1000)
#  disc.disp.DEDD.2.lnHM[t] <- sum(p.disp.DEDD.2.lnHM < t/1000)
#  disc.log.disp.DEDD.2.lnHM[t] <- sum(p.log.disp.DEDD.2.lnHM < t/1000)
#  fd.mean.DEDD.2.expHM[t] <- sum(p.mean.DEDD.2.expHM < t/1000 & DE.DEDD.2==0)
#  fd.log.mean.DEDD.2.expHM[t] <- sum(p.log.mean.DEDD.2.expHM < t/1000 & DE.DEDD.2==0)
#  fd.disp.DEDD.2.expHM[t] <- sum(p.disp.DEDD.2.expHM < t/1000 & DD.DEDD.2==0)
#  fd.log.disp.DEDD.2.expHM[t] <- sum(p.log.disp.DEDD.2.expHM < t/1000 & DD.DEDD.2==0)
#  fd.mean.DEDD.2.lnHM[t] <- sum(p.mean.DEDD.2.lnHM < t/1000 & DE.DEDD.2==0)
#  fd.log.mean.DEDD.2.lnHM[t] <- sum(p.log.mean.DEDD.2.lnHM < t/1000 & DE.DEDD.2==0)
#  fd.disp.DEDD.2.lnHM[t] <- sum(p.disp.DEDD.2.lnHM < t/1000 & DD.DEDD.2==0)
#  fd.log.disp.DEDD.2.lnHM[t] <- sum(p.log.disp.DEDD.2.lnHM < t/1000 & DD.DEDD.2==0)
#}
#fdr.mean.DEDD.2.expHM <- fd.mean.DEDD.2.expHM / disc.mean.DEDD.2.expHM  
#fdr.log.mean.DEDD.2.expHM <- fd.log.mean.DEDD.2.expHM / disc.log.mean.DEDD.2.expHM  
#fdr.disp.DEDD.2.expHM <- fd.disp.DEDD.2.expHM / disc.disp.DEDD.2.expHM  
#fdr.log.disp.DEDD.2.expHM <- fd.log.disp.DEDD.2.expHM / disc.log.disp.DEDD.2.expHM  
#fdr.mean.DEDD.2.lnHM <- fd.mean.DEDD.2.lnHM / disc.mean.DEDD.2.lnHM  
#fdr.log.mean.DEDD.2.lnHM <- fd.log.mean.DEDD.2.lnHM / disc.log.mean.DEDD.2.lnHM  
#fdr.disp.DEDD.2.lnHM <- fd.disp.DEDD.2.lnHM / disc.disp.DEDD.2.lnHM  
#fdr.log.disp.DEDD.2.lnHM <- fd.log.disp.DEDD.2.lnHM / disc.log.disp.DEDD.2.lnHM  
#par(mfrow=c(4,2))
#plot(disc.mean.DEDD.2.expHM, fd.mean.DEDD.2.expHM, type='l')
#plot(disc.log.mean.DEDD.2.expHM, fd.log.mean.DEDD.2.expHM, type='l')
#plot(disc.disp.DEDD.2.expHM, fd.disp.DEDD.2.expHM, type='l')
#plot(disc.log.disp.DEDD.2.expHM, fd.log.disp.DEDD.2.expHM, type='l')
#plot(disc.mean.DEDD.2.lnHM, fd.mean.DEDD.2.lnHM, type='l')
#plot(disc.log.mean.DEDD.2.lnHM, fd.log.mean.DEDD.2.lnHM, type='l')
#plot(disc.disp.DEDD.2.lnHM, fd.disp.DEDD.2.lnHM, type='l')
#plot(disc.log.disp.DEDD.2.lnHM, fd.log.disp.DEDD.2.lnHM, type='l')
#plot(disc.mean.DEDD.2.expHM, fdr.mean.DEDD.2.expHM, type='l')
#plot(disc.log.mean.DEDD.2.expHM, fdr.log.mean.DEDD.2.expHM, type='l')
#plot(disc.disp.DEDD.2.expHM, fdr.disp.DEDD.2.expHM, type='l')
#plot(disc.log.disp.DEDD.2.expHM, fdr.log.disp.DEDD.2.expHM, type='l')
#plot(disc.mean.DEDD.2.lnHM, fdr.mean.DEDD.2.lnHM, type='l')
#plot(disc.log.mean.DEDD.2.lnHM, fdr.log.mean.DEDD.2.lnHM, type='l')
#plot(disc.disp.DEDD.2.lnHM, fdr.disp.DEDD.2.lnHM, type='l')
#plot(disc.log.disp.DEDD.2.lnHM, fdr.log.disp.DEDD.2.lnHM, type='l')





