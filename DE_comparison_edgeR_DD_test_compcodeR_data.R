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
source(here('scripts','2019-04-03_exponential_hmm_adaptive_proposals_three_chains_function.R'))
source(here('scripts','2019-03-27_lognormal_hmm_adaptive_proposals_three_chains_function.R'))
source(here('scripts','2019-05-03_bfdr_function.R'))
source(here('scripts','2019-05-03_hpd_tail_prob_function.R'))
source(here('scripts','2019-05-07_symmetric_tail_prob_function.R'))
source(here('scripts','2019-05-17_compData_diff_disp_functions.R'))

samples.per.cond <- 5
group <- factor(c(rep(1,samples.per.cond), rep(2,samples.per.cond)))
FPR.DEDD.5.ql.edgeR <- numeric(5); FPR.DEDD.5.lr.edgeR <- numeric(5); FPR.DEDD.5.et.edgeR <- numeric(5)
FDR.DEDD.5.ql.edgeR <- numeric(5); FDR.DEDD.5.lr.edgeR <- numeric(5); FDR.DEDD.5.et.edgeR <- numeric(5)
TPR.DEDD.5.ql.edgeR <- numeric(5); TPR.DEDD.5.lr.edgeR <- numeric(5); TPR.DEDD.5.et.edgeR <- numeric(5)
AUC.DEDD.5.ql.edgeR <- numeric(5); AUC.DEDD.5.lr.edgeR <- numeric(5); AUC.DEDD.5.et.edgeR <- numeric(5)
FPR.DEDD.5.hmm.5.expHM <- numeric(5); FPR.DEDD.5.hmm.prop.expHM <- numeric(5); FPR.DEDD.5.mean.expHM <- numeric(5)
FPR.DEDD.5.log.mean.expHM <- numeric(5); FPR.DEDD.5.disp.expHM <- numeric(5); FPR.DEDD.5.log.disp.expHM <- numeric(5)
FDR.DEDD.5.hmm.5.expHM <- numeric(5); FDR.DEDD.5.hmm.prop.expHM <- numeric(5); FDR.DEDD.5.hmm.q.expHM <- numeric(5)
FDR.DEDD.5.mean.expHM <- numeric(5); FDR.DEDD.5.log.mean.expHM <- numeric(5); FDR.DEDD.5.disp.expHM <- numeric(5)
FDR.DEDD.5.log.disp.expHM <- numeric(5); TPR.DEDD.5.hmm.5.expHM <- numeric(5); TPR.DEDD.5.hmm.prop.expHM <- numeric(5)
TPR.DEDD.5.hmm.q.expHM <- numeric(5); TPR.DEDD.5.mean.expHM <- numeric(5); TPR.DEDD.5.log.mean.expHM <- numeric(5)
TPR.DEDD.5.disp.expHM <- numeric(5); TPR.DEDD.5.log.disp.expHM <- numeric(5); AUC.DEDD.5.hmm.expHM <- numeric(5)
AUC.DEDD.5.mean.expHM <- numeric(5); AUC.DEDD.5.log.mean.expHM <- numeric(5); AUC.DEDD.5.disp.expHM <- numeric(5)
AUC.DEDD.5.log.disp.expHM <- numeric(5); FPR.DEDD.5.hmm.5.lnHM <- numeric(5); FPR.DEDD.5.hmm.prop.lnHM <- numeric(5)
FPR.DEDD.5.mean.lnHM <- numeric(5); FPR.DEDD.5.log.mean.lnHM <- numeric(5); FPR.DEDD.5.disp.lnHM <- numeric(5)
FPR.DEDD.5.log.disp.lnHM <- numeric(5); FDR.DEDD.5.hmm.5.lnHM <- numeric(5); FDR.DEDD.5.hmm.prop.lnHM <- numeric(5)
FDR.DEDD.5.hmm.q.lnHM <- numeric(5); FDR.DEDD.5.mean.lnHM <- numeric(5); FDR.DEDD.5.log.mean.lnHM <- numeric(5)
FDR.DEDD.5.disp.lnHM <- numeric(5); FDR.DEDD.5.log.disp.lnHM <- numeric(5); TPR.DEDD.5.hmm.5.lnHM <- numeric(5)
TPR.DEDD.5.hmm.prop.lnHM <- numeric(5); TPR.DEDD.5.hmm.q.lnHM <- numeric(5); TPR.DEDD.5.mean.lnHM <- numeric(5)
TPR.DEDD.5.log.mean.lnHM <- numeric(5); TPR.DEDD.5.disp.lnHM <- numeric(5); TPR.DEDD.5.log.disp.lnHM <- numeric(5)
AUC.DEDD.5.hmm.lnHM <- numeric(5); AUC.DEDD.5.mean.lnHM <- numeric(5); AUC.DEDD.5.log.mean.lnHM <- numeric(5)
AUC.DEDD.5.disp.lnHM <- numeric(5); AUC.DEDD.5.log.disp.lnHM <- numeric(5)

for (i in 1:5) {
  # Generate filtered data
  counts.DEDD.5 <- simulate.DE.DD.data(dataset='DEDD.5', n.vars=2e4, samples.per.cond=samples.per.cond)
  DE.DEDD.5 <- counts.DEDD.5@variable.annotations$differential.expression
  DD.DEDD.5 <- counts.DEDD.5@variable.annotations$differential.dispersion
  truemeans.DEDD.5.1 <- counts.DEDD.5@variable.annotations$truemeans.S1
  truemeans.DEDD.5.2 <- counts.DEDD.5@variable.annotations$truemeans.S2
  truedisps.DEDD.5.1 <- counts.DEDD.5@variable.annotations$truedispersions.S1
  truedisps.DEDD.5.2 <- counts.DEDD.5@variable.annotations$truedispersions.S2
  
  # Normalise
  nf.DEDD.5 <- calcNormFactors(counts.DEDD.5@count.matrix)
  norm.DEDD.5 <- t(t(counts.DEDD.5@count.matrix) / nf.DEDD.5)
  
  ## edgeR
  DEDD.5.edgeR <- DGEList(counts=counts.DEDD.5@count.matrix, norm.factors=nf.DEDD.5, group=group)
  design.edgeR <- model.matrix(~group)
  DEDD.5.edgeR <- estimateDisp(DEDD.5.edgeR, design.edgeR)
  qlfit.edgeR <- glmQLFit(DEDD.5.edgeR, design.edgeR)
  qltest.edgeR <- glmQLFTest(qlfit.edgeR)
  lrfit.edgeR <- glmFit(DEDD.5.edgeR, design.edgeR)
  lrtest.edgeR <- glmLRT(lrfit.edgeR)
  et.edgeR <- exactTest(DEDD.5.edgeR)
  pval.DEDD.5.ql.edgeR <- qltest.edgeR$table$PValue
  pval.DEDD.5.lr.edgeR <- lrtest.edgeR$table$PValue
  pval.DEDD.5.et.edgeR <- et.edgeR$table$PValue
  qval.DEDD.5.ql.edgeR <- topTags(qltest.edgeR, n=nrow(DEDD.5.edgeR$counts), sort='none')$table$FDR
  qval.DEDD.5.lr.edgeR <- topTags(lrtest.edgeR, n=nrow(DEDD.5.edgeR$counts), sort='none')$table$FDR
  qval.DEDD.5.et.edgeR <- topTags(et.edgeR, n=nrow(DEDD.5.edgeR$counts), sort='none')$table$FDR
  pr.DEDD.5.ql.edgeR <- prediction(1-pval.DEDD.5.ql.edgeR, DE.DEDD.5)
  pr.DEDD.5.lr.edgeR <- prediction(1-pval.DEDD.5.lr.edgeR, DE.DEDD.5)
  pr.DEDD.5.et.edgeR <- prediction(1-pval.DEDD.5.et.edgeR, DE.DEDD.5)
  
  FPR.DEDD.5.ql.edgeR[i] <- sum(pval.DEDD.5.ql.edgeR<0.05 & DE.DEDD.5==0) / sum(DE.DEDD.5==0)
  FPR.DEDD.5.lr.edgeR[i] <- sum(pval.DEDD.5.lr.edgeR<0.05 & DE.DEDD.5==0) / sum(DE.DEDD.5==0)
  FPR.DEDD.5.et.edgeR[i] <- sum(pval.DEDD.5.et.edgeR<0.05 & DE.DEDD.5==0) / sum(DE.DEDD.5==0)
  FDR.DEDD.5.ql.edgeR[i] <- sum(qval.DEDD.5.ql.edgeR<0.05 & DE.DEDD.5==0) / sum(qval.DEDD.5.ql.edgeR<0.05)
  FDR.DEDD.5.lr.edgeR[i] <- sum(qval.DEDD.5.lr.edgeR<0.05 & DE.DEDD.5==0) / sum(qval.DEDD.5.lr.edgeR<0.05)
  FDR.DEDD.5.et.edgeR[i] <- sum(qval.DEDD.5.et.edgeR<0.05 & DE.DEDD.5==0) / sum(qval.DEDD.5.et.edgeR<0.05)
  TPR.DEDD.5.ql.edgeR[i] <- sum(qval.DEDD.5.ql.edgeR<0.05 & DE.DEDD.5==1) / sum(DE.DEDD.5==1)
  TPR.DEDD.5.lr.edgeR[i] <- sum(qval.DEDD.5.lr.edgeR<0.05 & DE.DEDD.5==1) / sum(DE.DEDD.5==1)
  TPR.DEDD.5.et.edgeR[i] <- sum(qval.DEDD.5.et.edgeR<0.05 & DE.DEDD.5==1) / sum(DE.DEDD.5==1)
  AUC.DEDD.5.ql.edgeR[i] <- performance(pr.DEDD.5.ql.edgeR, measure='auc')@y.values[[1]]
  AUC.DEDD.5.lr.edgeR[i] <- performance(pr.DEDD.5.lr.edgeR, measure='auc')@y.values[[1]]
  AUC.DEDD.5.et.edgeR[i] <- performance(pr.DEDD.5.et.edgeR, measure='auc')@y.values[[1]]

  ## expHM
  DEDD.5.expHM <- exp_hmm_adapt_3_chains(counts=t(norm.DEDD.5), groups=group)
  prob.DEDD.5.expHM <- colMeans(as.matrix(DEDD.5.expHM$indicators))
  post.prop.DEDD.5.expHM <- mean(as.matrix(DEDD.5.expHM$proportion))
  thr.DEDD.5.expHM <- sort(prob.DEDD.5.expHM, decreasing=T)[round(nrow(counts.DEDD.5@count.matrix) * 
                                                                    post.prop.DEDD.5.expHM)]
  mean.diff.DEDD.5.expHM <- as.matrix(DEDD.5.expHM$means1) - as.matrix(DEDD.5.expHM$means2)
  log.mean.diff.DEDD.5.expHM <- log(as.matrix(DEDD.5.expHM$means1)) - log(as.matrix(DEDD.5.expHM$means2))
  disp.diff.DEDD.5.expHM <- as.matrix(DEDD.5.expHM$disps1) - as.matrix(DEDD.5.expHM$disps2)
  log.disp.diff.DEDD.5.expHM <- log(as.matrix(DEDD.5.expHM$disps1)) - log(as.matrix(DEDD.5.expHM$disps2))
  p.mean.DEDD.5.expHM <- apply(mean.diff.DEDD.5.expHM,2,hpd.pval)
  p.log.mean.DEDD.5.expHM <- apply(log.mean.diff.DEDD.5.expHM,2,hpd.pval)
  p.disp.DEDD.5.expHM <- apply(disp.diff.DEDD.5.expHM,2,hpd.pval)
  p.log.disp.DEDD.5.expHM <- apply(log.disp.diff.DEDD.5.expHM,2,hpd.pval)
  q.hmm.DEDD.5.expHM <- bfdr(prob.DEDD.5.expHM)
  q.mean.DEDD.5.expHM <- qvalue(p.mean.DEDD.5.expHM)$qvalues
  q.log.mean.DEDD.5.expHM <- qvalue(p.log.mean.DEDD.5.expHM)$qvalues
  q.disp.DEDD.5.expHM <- qvalue(p.disp.DEDD.5.expHM)$qvalues
  q.log.disp.DEDD.5.expHM <- qvalue(p.log.disp.DEDD.5.expHM)$qvalues
  pr.DEDD.5.hmm.expHM <- prediction(prob.DEDD.5.expHM, (DE.DEDD.5 | DD.DEDD.5))
  pr.DEDD.5.mean.expHM <- prediction(1-p.mean.DEDD.5.expHM, DE.DEDD.5)
  pr.DEDD.5.log.mean.expHM <- prediction(1-p.log.mean.DEDD.5.expHM, DE.DEDD.5)
  pr.DEDD.5.disp.expHM <- prediction(1-p.disp.DEDD.5.expHM, DD.DEDD.5)
  pr.DEDD.5.log.disp.expHM <- prediction(1-p.log.disp.DEDD.5.expHM, DD.DEDD.5)
  
  FPR.DEDD.5.hmm.5.expHM[i] <- sum(prob.DEDD.5.expHM>0.5 & DE.DEDD.5==0 & DD.DEDD.5==0) / 
    sum(DE.DEDD.5==0 & DD.DEDD.5==0)
  FPR.DEDD.5.hmm.prop.expHM[i] <- sum(prob.DEDD.5.expHM>thr.DEDD.5.expHM & DE.DEDD.5==0 & DD.DEDD.5==0) / 
    sum(DE.DEDD.5==0 & DD.DEDD.5==0)
  FPR.DEDD.5.mean.expHM[i] <- sum(p.mean.DEDD.5.expHM<0.05 & DE.DEDD.5==0) / sum(DE.DEDD.5==0)
  FPR.DEDD.5.log.mean.expHM[i] <- sum(p.log.mean.DEDD.5.expHM<0.05 & DE.DEDD.5==0) / sum(DE.DEDD.5==0)
  FPR.DEDD.5.disp.expHM[i] <- sum(p.disp.DEDD.5.expHM<0.05 & DD.DEDD.5==0) / sum(DD.DEDD.5==0)
  FPR.DEDD.5.log.disp.expHM[i] <- sum(p.log.disp.DEDD.5.expHM<0.05 & DD.DEDD.5==0) / sum(DD.DEDD.5==0)
  FDR.DEDD.5.hmm.5.expHM[i] <- sum(prob.DEDD.5.expHM>0.5 & DE.DEDD.5==0 & DD.DEDD.5==0) / sum(prob.DEDD.5.expHM>0.5)
  FDR.DEDD.5.hmm.prop.expHM[i] <- sum(prob.DEDD.5.expHM>thr.DEDD.5.expHM & DE.DEDD.5==0 & DD.DEDD.5==0) / 
    sum(prob.DEDD.5.expHM>thr.DEDD.5.expHM)
  FDR.DEDD.5.hmm.q.expHM[i] <- sum(q.hmm.DEDD.5.expHM<0.05 & DE.DEDD.5==0 & DD.DEDD.5==0) / sum(q.hmm.DEDD.5.expHM<0.05)
  FDR.DEDD.5.mean.expHM[i] <- sum(q.mean.DEDD.5.expHM<0.05 & DE.DEDD.5==0) / sum(q.mean.DEDD.5.expHM<0.05)
  FDR.DEDD.5.log.mean.expHM[i] <- sum(q.log.mean.DEDD.5.expHM<0.05 & DE.DEDD.5==0) / sum(q.log.mean.DEDD.5.expHM<0.05)
  FDR.DEDD.5.disp.expHM[i] <- sum(q.disp.DEDD.5.expHM<0.05 & DD.DEDD.5==0) / sum(q.disp.DEDD.5.expHM<0.05)
  FDR.DEDD.5.log.disp.expHM[i] <- sum(q.log.disp.DEDD.5.expHM<0.05 & DD.DEDD.5==0) / sum(q.log.disp.DEDD.5.expHM<0.05)
  TPR.DEDD.5.hmm.5.expHM[i] <- sum(prob.DEDD.5.expHM>0.5 & (DE.DEDD.5==1 | DD.DEDD.5==1)) / 
    sum(DE.DEDD.5==1 | DD.DEDD.5==1)
  TPR.DEDD.5.hmm.prop.expHM[i] <- sum(prob.DEDD.5.expHM>thr.DEDD.5.expHM & (DE.DEDD.5==1 | DD.DEDD.5==1)) / 
    sum(DE.DEDD.5==1 | DD.DEDD.5==1)
  TPR.DEDD.5.hmm.q.expHM[i] <- sum(q.hmm.DEDD.5.expHM<0.05 & (DE.DEDD.5==1 | DD.DEDD.5==1)) / 
    sum(DE.DEDD.5==1 | DD.DEDD.5==1)
  TPR.DEDD.5.mean.expHM[i] <- sum(q.mean.DEDD.5.expHM<0.05 & DE.DEDD.5==1) / sum(DE.DEDD.5==1)
  TPR.DEDD.5.log.mean.expHM[i] <- sum(q.log.mean.DEDD.5.expHM<0.05 & DE.DEDD.5==1) / sum(DE.DEDD.5==1)
  TPR.DEDD.5.disp.expHM[i] <- sum(q.disp.DEDD.5.expHM<0.05 & DD.DEDD.5==1) / sum(DD.DEDD.5==1)
  TPR.DEDD.5.log.disp.expHM[i] <- sum(q.log.disp.DEDD.5.expHM<0.05 & DD.DEDD.5==1) / sum(DD.DEDD.5==1)
  AUC.DEDD.5.hmm.expHM[i] <- performance(pr.DEDD.5.hmm.expHM, measure='auc')@y.values[[1]]
  AUC.DEDD.5.mean.expHM[i] <- performance(pr.DEDD.5.mean.expHM, measure='auc')@y.values[[1]]
  AUC.DEDD.5.log.mean.expHM[i] <- performance(pr.DEDD.5.log.mean.expHM, measure='auc')@y.values[[1]]
  AUC.DEDD.5.disp.expHM[i] <- performance(pr.DEDD.5.disp.expHM, measure='auc')@y.values[[1]]
  AUC.DEDD.5.log.disp.expHM[i] <- performance(pr.DEDD.5.log.disp.expHM, measure='auc')@y.values[[1]]

  ## lnHM
  DEDD.5.lnHM <- ln_hmm_adapt_3_chains(counts=t(norm.DEDD.5), groups=group)
  prob.DEDD.5.lnHM <- colMeans(as.matrix(DEDD.5.lnHM$indicators))
  post.prop.DEDD.5.lnHM <- mean(as.matrix(DEDD.5.lnHM$proportion))
  thr.DEDD.5.lnHM <- sort(prob.DEDD.5.lnHM, decreasing=T)[round(nrow(counts.DEDD.5@count.matrix)*post.prop.DEDD.5.lnHM)]
  mean.diff.DEDD.5.lnHM <- as.matrix(DEDD.5.lnHM$means1) - as.matrix(DEDD.5.lnHM$means2)
  log.mean.diff.DEDD.5.lnHM <- log(as.matrix(DEDD.5.lnHM$means1)) - log(as.matrix(DEDD.5.lnHM$means2))
  disp.diff.DEDD.5.lnHM <- as.matrix(DEDD.5.lnHM$disps1) - as.matrix(DEDD.5.lnHM$disps2)
  log.disp.diff.DEDD.5.lnHM <- log(as.matrix(DEDD.5.lnHM$disps1)) - log(as.matrix(DEDD.5.lnHM$disps2))
  p.mean.DEDD.5.lnHM <- apply(mean.diff.DEDD.5.lnHM,2,hpd.pval)
  p.log.mean.DEDD.5.lnHM <- apply(log.mean.diff.DEDD.5.lnHM,2,hpd.pval)
  p.disp.DEDD.5.lnHM <- apply(disp.diff.DEDD.5.lnHM,2,hpd.pval)
  p.log.disp.DEDD.5.lnHM <- apply(log.disp.diff.DEDD.5.lnHM,2,hpd.pval)
  q.hmm.DEDD.5.lnHM <- bfdr(prob.DEDD.5.lnHM)
  q.mean.DEDD.5.lnHM <- qvalue(p.mean.DEDD.5.lnHM)$qvalues
  q.log.mean.DEDD.5.lnHM <- qvalue(p.log.mean.DEDD.5.lnHM)$qvalues
  q.disp.DEDD.5.lnHM <- qvalue(p.disp.DEDD.5.lnHM)$qvalues
  q.log.disp.DEDD.5.lnHM <- qvalue(p.log.disp.DEDD.5.lnHM)$qvalues
  pr.DEDD.5.hmm.lnHM <- prediction(prob.DEDD.5.lnHM, (DE.DEDD.5 | DD.DEDD.5))
  pr.DEDD.5.mean.lnHM <- prediction(1-p.mean.DEDD.5.lnHM, DE.DEDD.5)
  pr.DEDD.5.log.mean.lnHM <- prediction(1-p.log.mean.DEDD.5.lnHM, DE.DEDD.5)
  pr.DEDD.5.disp.lnHM <- prediction(1-p.disp.DEDD.5.lnHM, DD.DEDD.5)
  pr.DEDD.5.log.disp.lnHM <- prediction(1-p.log.disp.DEDD.5.lnHM, DD.DEDD.5)
  
  FPR.DEDD.5.hmm.5.lnHM[i] <- sum(prob.DEDD.5.lnHM>0.5 & DE.DEDD.5==0 & DD.DEDD.5==0) / sum(DE.DEDD.5==0 & DD.DEDD.5==0)
  FPR.DEDD.5.hmm.prop.lnHM[i] <- sum(prob.DEDD.5.lnHM>thr.DEDD.5.lnHM & DE.DEDD.5==0 & DD.DEDD.5==0) / 
    sum(DE.DEDD.5==0 & DD.DEDD.5==0)
  FPR.DEDD.5.mean.lnHM[i] <- sum(p.mean.DEDD.5.lnHM<0.05 & DE.DEDD.5==0) / sum(DE.DEDD.5==0)
  FPR.DEDD.5.log.mean.lnHM[i] <- sum(p.log.mean.DEDD.5.lnHM<0.05 & DE.DEDD.5==0) / sum(DE.DEDD.5==0)
  FPR.DEDD.5.disp.lnHM[i] <- sum(p.disp.DEDD.5.lnHM<0.05 & DD.DEDD.5==0) / sum(DD.DEDD.5==0)
  FPR.DEDD.5.log.disp.lnHM[i] <- sum(p.log.disp.DEDD.5.lnHM<0.05 & DD.DEDD.5==0) / sum(DD.DEDD.5==0)
  FDR.DEDD.5.hmm.5.lnHM[i] <- sum(prob.DEDD.5.lnHM>0.5 & DE.DEDD.5==0 & DD.DEDD.5==0) / sum(prob.DEDD.5.lnHM>0.5)
  FDR.DEDD.5.hmm.prop.lnHM[i] <- sum(prob.DEDD.5.lnHM>thr.DEDD.5.lnHM & DE.DEDD.5==0 & DD.DEDD.5==0) / 
    sum(prob.DEDD.5.lnHM>thr.DEDD.5.lnHM)
  FDR.DEDD.5.hmm.q.lnHM[i] <- sum(q.hmm.DEDD.5.lnHM<0.05 & DE.DEDD.5==0 & DD.DEDD.5==0) / sum(q.hmm.DEDD.5.lnHM<0.05)
  FDR.DEDD.5.mean.lnHM[i] <- sum(q.mean.DEDD.5.lnHM<0.05 & DE.DEDD.5==0) / sum(q.mean.DEDD.5.lnHM<0.05)
  FDR.DEDD.5.log.mean.lnHM[i] <- sum(q.log.mean.DEDD.5.lnHM<0.05 & DE.DEDD.5==0) / sum(q.log.mean.DEDD.5.lnHM<0.05)
  FDR.DEDD.5.disp.lnHM[i] <- sum(q.disp.DEDD.5.lnHM<0.05 & DD.DEDD.5==0) / sum(q.disp.DEDD.5.lnHM<0.05)
  FDR.DEDD.5.log.disp.lnHM[i] <- sum(q.log.disp.DEDD.5.lnHM<0.05 & DD.DEDD.5==0) / sum(q.log.disp.DEDD.5.lnHM<0.05)
  TPR.DEDD.5.hmm.5.lnHM[i] <- sum(prob.DEDD.5.lnHM>0.5 & (DE.DEDD.5==1 | DD.DEDD.5==1)) / 
    sum(DE.DEDD.5==1 | DD.DEDD.5==1)
  TPR.DEDD.5.hmm.prop.lnHM[i] <- sum(prob.DEDD.5.lnHM>thr.DEDD.5.lnHM & (DE.DEDD.5==1 | DD.DEDD.5==1)) / 
    sum(DE.DEDD.5==1 | DD.DEDD.5==1)
  TPR.DEDD.5.hmm.q.lnHM[i] <- sum(q.hmm.DEDD.5.lnHM<0.05 & (DE.DEDD.5==1 | DD.DEDD.5==1)) / 
    sum(DE.DEDD.5==1 | DD.DEDD.5==1)
  TPR.DEDD.5.mean.lnHM[i] <- sum(q.mean.DEDD.5.lnHM<0.05 & DE.DEDD.5==1) / sum(DE.DEDD.5==1)
  TPR.DEDD.5.log.mean.lnHM[i] <- sum(q.log.mean.DEDD.5.lnHM<0.05 & DE.DEDD.5==1) / sum(DE.DEDD.5==1)
  TPR.DEDD.5.disp.lnHM[i] <- sum(q.disp.DEDD.5.lnHM<0.05 & DD.DEDD.5==1) / sum(DD.DEDD.5==1)
  TPR.DEDD.5.log.disp.lnHM[i] <- sum(q.log.disp.DEDD.5.lnHM<0.05 & DD.DEDD.5==1) / sum(DD.DEDD.5==1)
  AUC.DEDD.5.hmm.lnHM[i] <- performance(pr.DEDD.5.hmm.lnHM, measure='auc')@y.values[[1]]
  AUC.DEDD.5.mean.lnHM[i] <- performance(pr.DEDD.5.mean.lnHM, measure='auc')@y.values[[1]]
  AUC.DEDD.5.log.mean.lnHM[i] <- performance(pr.DEDD.5.log.mean.lnHM, measure='auc')@y.values[[1]]
  AUC.DEDD.5.disp.lnHM[i] <- performance(pr.DEDD.5.disp.lnHM, measure='auc')@y.values[[1]]
  AUC.DEDD.5.log.disp.lnHM[i] <- performance(pr.DEDD.5.log.disp.lnHM, measure='auc')@y.values[[1]]
}

results <- list('FPR.DEDD.5.ql.edgeR' = FPR.DEDD.5.ql.edgeR, 
                'FPR.DEDD.5.lr.edgeR' = FPR.DEDD.5.lr.edgeR, 
                'FPR.DEDD.5.et.edgeR' = FPR.DEDD.5.et.edgeR, 
                'FDR.DEDD.5.ql.edgeR' = FDR.DEDD.5.ql.edgeR, 
                'FDR.DEDD.5.lr.edgeR' = FDR.DEDD.5.lr.edgeR, 
                'FDR.DEDD.5.et.edgeR' = FDR.DEDD.5.et.edgeR, 
                'TPR.DEDD.5.ql.edgeR' = TPR.DEDD.5.ql.edgeR, 
                'TPR.DEDD.5.lr.edgeR' = TPR.DEDD.5.lr.edgeR, 
                'TPR.DEDD.5.et.edgeR' = TPR.DEDD.5.et.edgeR, 
                'AUC.DEDD.5.ql.edgeR' = AUC.DEDD.5.ql.edgeR, 
                'AUC.DEDD.5.lr.edgeR' = AUC.DEDD.5.lr.edgeR, 
                'AUC.DEDD.5.et.edgeR' = AUC.DEDD.5.et.edgeR, 
                'FPR.DEDD.5.hmm.5.expHM' = FPR.DEDD.5.hmm.5.expHM, 
                'FPR.DEDD.5.hmm.prop.expHM' = FPR.DEDD.5.hmm.prop.expHM, 
                'FPR.DEDD.5.mean.expHM' = FPR.DEDD.5.mean.expHM, 
                'FPR.DEDD.5.log.mean.expHM' = FPR.DEDD.5.log.mean.expHM, 
                'FPR.DEDD.5.disp.expHM' = FPR.DEDD.5.disp.expHM, 
                'FPR.DEDD.5.log.disp.expHM' = FPR.DEDD.5.log.disp.expHM, 
                'FDR.DEDD.5.hmm.5.expHM' = FDR.DEDD.5.hmm.5.expHM, 
                'FDR.DEDD.5.hmm.prop.expHM' = FDR.DEDD.5.hmm.prop.expHM, 
                'FDR.DEDD.5.hmm.q.expHM' = FDR.DEDD.5.hmm.q.expHM, 
                'FDR.DEDD.5.mean.expHM' = FDR.DEDD.5.mean.expHM, 
                'FDR.DEDD.5.log.mean.expHM' = FDR.DEDD.5.log.mean.expHM, 
                'FDR.DEDD.5.disp.expHM' = FDR.DEDD.5.disp.expHM, 
                'FDR.DEDD.5.log.disp.expHM' = FDR.DEDD.5.log.disp.expHM, 
                'TPR.DEDD.5.hmm.5.expHM' = TPR.DEDD.5.hmm.5.expHM, 
                'TPR.DEDD.5.hmm.prop.expHM' = TPR.DEDD.5.hmm.prop.expHM, 
                'TPR.DEDD.5.hmm.q.expHM' = TPR.DEDD.5.hmm.q.expHM, 
                'TPR.DEDD.5.mean.expHM' = TPR.DEDD.5.mean.expHM, 
                'TPR.DEDD.5.log.mean.expHM' = TPR.DEDD.5.log.mean.expHM, 
                'TPR.DEDD.5.disp.expHM' = TPR.DEDD.5.disp.expHM, 
                'TPR.DEDD.5.log.disp.expHM' = TPR.DEDD.5.log.disp.expHM, 
                'AUC.DEDD.5.hmm.expHM' = AUC.DEDD.5.hmm.expHM, 
                'AUC.DEDD.5.mean.expHM' = AUC.DEDD.5.mean.expHM, 
                'AUC.DEDD.5.log.mean.expHM' = AUC.DEDD.5.log.mean.expHM, 
                'AUC.DEDD.5.disp.expHM' = AUC.DEDD.5.disp.expHM, 
                'AUC.DEDD.5.log.disp.expHM' = AUC.DEDD.5.log.disp.expHM, 
                'FPR.DEDD.5.hmm.5.lnHM' = FPR.DEDD.5.hmm.5.lnHM, 
                'FPR.DEDD.5.hmm.prop.lnHM' = FPR.DEDD.5.hmm.prop.lnHM, 
                'FPR.DEDD.5.mean.lnHM' = FPR.DEDD.5.mean.lnHM, 
                'FPR.DEDD.5.log.mean.lnHM' = FPR.DEDD.5.log.mean.lnHM, 
                'FPR.DEDD.5.disp.lnHM' = FPR.DEDD.5.disp.lnHM, 
                'FPR.DEDD.5.log.disp.lnHM' = FPR.DEDD.5.log.disp.lnHM, 
                'FDR.DEDD.5.hmm.5.lnHM' = FDR.DEDD.5.hmm.5.lnHM, 
                'FDR.DEDD.5.hmm.prop.lnHM' = FDR.DEDD.5.hmm.prop.lnHM, 
                'FDR.DEDD.5.hmm.q.lnHM' = FDR.DEDD.5.hmm.q.lnHM, 
                'FDR.DEDD.5.mean.lnHM' = FDR.DEDD.5.mean.lnHM, 
                'FDR.DEDD.5.log.mean.lnHM' = FDR.DEDD.5.log.mean.lnHM, 
                'FDR.DEDD.5.disp.lnHM' = FDR.DEDD.5.disp.lnHM, 
                'FDR.DEDD.5.log.disp.lnHM' = FDR.DEDD.5.log.disp.lnHM, 
                'TPR.DEDD.5.hmm.5.lnHM' = TPR.DEDD.5.hmm.5.lnHM, 
                'TPR.DEDD.5.hmm.prop.lnHM' = TPR.DEDD.5.hmm.prop.lnHM, 
                'TPR.DEDD.5.hmm.q.lnHM' = TPR.DEDD.5.hmm.q.lnHM, 
                'TPR.DEDD.5.mean.lnHM' = TPR.DEDD.5.mean.lnHM, 
                'TPR.DEDD.5.log.mean.lnHM' = TPR.DEDD.5.log.mean.lnHM, 
                'TPR.DEDD.5.disp.lnHM' = TPR.DEDD.5.disp.lnHM, 
                'TPR.DEDD.5.log.disp.lnHM' = TPR.DEDD.5.log.disp.lnHM, 
                'AUC.DEDD.5.hmm.lnHM' = AUC.DEDD.5.hmm.lnHM, 
                'AUC.DEDD.5.mean.lnHM' = AUC.DEDD.5.mean.lnHM, 
                'AUC.DEDD.5.log.mean.lnHM' = AUC.DEDD.5.log.mean.lnHM, 
                'AUC.DEDD.5.disp.lnHM' = AUC.DEDD.5.disp.lnHM, 
                'AUC.DEDD.5.log.disp.lnHM' = AUC.DEDD.5.log.disp.lnHM)
filename <- paste0('results.DEDD.5.', format(Sys.time(), "%Y-%m-%d_%H%M%S"), '.rds')
saveRDS(results, file=here(filename))

res1 <- readRDS(here('results.DEDD.5.2019-05-29_190501.rds')) # 2 runs
res2 <- readRDS(here('results.DEDD.5.2019-05-29_203110.rds')) # 2 runs
res3 <- readRDS(here('results.DEDD.5.2019-05-29_233554.rds')) # 1 run
res4 <- readRDS(here('results.DEDD.5.2019-05-30_020313.rds')) # 2 runs
res5 <- readRDS(here('results.DEDD.5.2019-05-30_022216.rds')) # 2 runs
res6 <- readRDS(here('results.DEDD.5.2019-05-30_042507.rds')) # 5 runs
results <- mapply(c,res1,res2,res3,res4,res5,res6)
empty <- which(results[,62]==0)
results <- results[-empty,]
colnames(results)
FPR <- as.data.frame(results[,c(1:3,13:18,38:43)])
FDR <- as.data.frame(results[,c(4:6,19:25,44:50)])
TPR <- as.data.frame(results[,c(7:9,26:32,51:57)])
AUC <- as.data.frame(results[,c(10:12,33:37,58:62)])
names(FPR) <- c("FPR.DEDD.5.ql.edgeR", "FPR.DEDD.5.lr.edgeR", "FPR.DEDD.5.et.edgeR", "FPR.DEDD.5.hmm.5.expHM", 
                "FPR.DEDD.5.hmm.prop.expHM", "FPR.DEDD.5.mean.expHM", "FPR.DEDD.5.log.mean.expHM", 
                "FPR.DEDD.5.disp.expHM", "FPR.DEDD.5.log.disp.expHM", "FPR.DEDD.5.hmm.5.lnHM", 
                "FPR.DEDD.5.hmm.prop.lnHM", "FPR.DEDD.5.mean.lnHM", "FPR.DEDD.5.log.mean.lnHM", "FPR.DEDD.5.disp.lnHM", 
                "FPR.DEDD.5.log.disp.lnHM")
names(FDR) <- c("FDR.DEDD.5.ql.edgeR", "FDR.DEDD.5.lr.edgeR", "FDR.DEDD.5.et.edgeR", "FDR.DEDD.5.hmm.5.expHM", 
                "FDR.DEDD.5.hmm.prop.expHM", "FDR.DEDD.5.hmm.q.expHM", "FDR.DEDD.5.mean.expHM", 
                "FDR.DEDD.5.log.mean.expHM", "FDR.DEDD.5.disp.expHM", "FDR.DEDD.5.log.disp.expHM", 
                "FDR.DEDD.5.hmm.5.lnHM", "FDR.DEDD.5.hmm.prop.lnHM", "FDR.DEDD.5.hmm.q.lnHM", "FDR.DEDD.5.mean.lnHM", 
                "FDR.DEDD.5.log.mean.lnHM", "FDR.DEDD.5.disp.lnHM", "FDR.DEDD.5.log.disp.lnHM")
names(TPR) <- c("TPR.DEDD.5.ql.edgeR", "TPR.DEDD.5.lr.edgeR", "TPR.DEDD.5.et.edgeR", "TPR.DEDD.5.hmm.5.expHM", 
                "TPR.DEDD.5.hmm.prop.expHM", "TPR.DEDD.5.hmm.q.expHM", "TPR.DEDD.5.mean.expHM", 
                "TPR.DEDD.5.log.mean.expHM", "TPR.DEDD.5.disp.expHM", "TPR.DEDD.5.log.disp.expHM", 
                "TPR.DEDD.5.hmm.5.lnHM", "TPR.DEDD.5.hmm.prop.lnHM", "TPR.DEDD.5.hmm.q.lnHM", "TPR.DEDD.5.mean.lnHM", 
                "TPR.DEDD.5.log.mean.lnHM", "TPR.DEDD.5.disp.lnHM", "TPR.DEDD.5.log.disp.lnHM")
names(AUC) <- c("AUC.DEDD.5.ql.edgeR", "AUC.DEDD.5.lr.edgeR", "AUC.DEDD.5.et.edgeR", "AUC.DEDD.5.hmm.expHM", 
                "AUC.DEDD.5.mean.expHM", "AUC.DEDD.5.log.mean.expHM", "AUC.DEDD.5.disp.expHM", 
                "AUC.DEDD.5.log.disp.expHM", "AUC.DEDD.5.hmm.lnHM", "AUC.DEDD.5.mean.lnHM", "AUC.DEDD.5.log.mean.lnHM", 
                "AUC.DEDD.5.disp.lnHM", "AUC.DEDD.5.log.disp.lnHM")
results <- list('FPR'=FPR, 'FDR'=FDR, 'TPR'=TPR, 'AUC'=AUC)
saveRDS(results, file=here('Results','2019-05-30_DEDD.5_comparison_14_runs.rds'))

results <- readRDS(here('Results','2019-05-30_DEDD.5_comparison_14_runs.rds'))
colours1 <- c(rep('orange',3),rep('blue',6),rep('red',6))
colours2 <- c(rep('orange',3),rep('blue',7),rep('red',7))
colours3 <- c(rep('orange',3),rep('blue',5),rep('red',5))
par(mfrow=c(4,1), mar=c(3,3,1,1), mgp=c(2,0.5,0))
boxplot(results$FPR, col=colours1, names=NULL); abline(h=0.05, v=c(3.5,5.5,7.5,9.5,11.5,13.5))
colMeans(results$FPR*100)
# All edgeR methods very good. HMM all zero. Mean conservative but log lnHM reasonable. 
# Dispersion very conservative, log expHM best.
boxplot(results$FDR, col=colours2, names=NULL); abline(h=0.05, v=c(3.5,6.5,8.5,10.5,13.5,15.5))
colMeans(results$FDR*100)
# edgeR QL good, LR and ET very liberal. HMM conservative but highly variable; 50% threshold particularly 
# variable but usually close to zero for expHM; posterior proportion threshold also very variable but not 
# too bad on average; BFDR always 0 or NaN, i.e. always no false discoveries and sometimes no true discoveries. 
# Mean reasonable but very variable; log expHM and untransformed lnHM closest to 0.05 on average. Dispersion 
# always NaN, i.e. no false or true discoveries.
boxplot(results$TPR, col=colours2, names=NULL); abline(v=c(3.5,6.5,8.5,10.5,13.5,15.5))
colMeans(results$TPR)
# edgeR clearly better than HM; LR and ET better than QL but since they have poor FDR control, QL is best to compare. 
# HMM all very low; posterior proportion highest, followed by 50%; BFDR often zero and if not then close to it; 
# since posterior proportion has highest power and decent FDR control generally, this is the best to use. Highest 
# power for mean is log transformed lnHM, but this has slightly liberal FDR control; untransformed lnHM has next 
# highest power and reasonable FDR control, so this is probably best to use. Dispersion always zero, i.e. no true
# positives (from FDR results, also no false positives). Could consider using raw p-values for dispersion.
boxplot(results$AUC, col=colours3, names=NULL); abline(v=c(3.5,4.5,6.5,8.5,9.5,11.5))
colMeans(results$AUC)
# edgeR very good and all variants similar, average 0.90. lnHM better than expHM for HMM, average 0.72. 
# Means good and all similar but slightly worse than edgeR; lnHM slightly better than expHM, average 0.89. 
# Dispersions poor and all similar; expHM slightly better than lnHM, average 0.56.


