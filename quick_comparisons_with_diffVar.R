library(here)
library(missMethyl)
library(limma)
library(ROCR)
library(caret)
library(compcodeR)

## DD, simulated data, 10 samples per group
DD10sim <- readRDS(here("Results/DD compcodeR data results DESeq norm Sept 2019", 
                     "results.DD10.1.DESeqnorm.rds"))
counts <- DD10sim$data@count.matrix
DD <- DD10sim$DD
group <- factor(c(rep("1",10), rep("2",10)))
design <- model.matrix(~group)
fit <- varFit(counts, design=design, coef=c(1,2))
DV <- topVar(fit, coef=2, number=nrow(counts), sort=F)
p <- DV$P.Value
q <- DV$Adj.P.Value
bh.disp.lnHM <- p.adjust(DD10sim$p.disp.lnHM, method="BH")
by.disp.lnHM <- p.adjust(DD10sim$p.disp.lnHM, method="BY")
holm.disp.lnHM <- p.adjust(DD10sim$p.disp.lnHM, method="holm")

fdr.MDSeq.DD10sim <- 1-precision(factor(DD10sim$q.disp.zi.MDSeq < 0.05, levels=c("TRUE", "FALSE")), 
                         factor(DD==1, levels=c("TRUE", "FALSE")))
fdr.bh.DD10sim <- 1-precision(factor(bh.disp.lnHM < 0.05, levels=c("TRUE", "FALSE")), 
                      factor(DD==1, levels=c("TRUE", "FALSE")))
fdr.by.DD10sim <- 1-precision(factor(by.disp.lnHM < 0.05, levels=c("TRUE", "FALSE")), 
                      factor(DD==1, levels=c("TRUE", "FALSE")))
fdr.holm.DD10sim <- 1-precision(factor(holm.disp.lnHM < 0.05, levels=c("TRUE", "FALSE")), 
                        factor(DD==1, levels=c("TRUE", "FALSE")))
fdr.diffvar.DD10sim <- 1-precision(factor(q < 0.05, levels=c("TRUE", "FALSE")), 
                           factor(DD==1, levels=c("TRUE", "FALSE")))

tpr.MDSeq.DD10sim <- sensitivity(factor(DD10sim$q.disp.zi.MDSeq < 0.05, levels=c("TRUE", "FALSE")), 
                         factor(DD==1, levels=c("TRUE", "FALSE")))
tpr.bh.DD10sim <- sensitivity(factor(bh.disp.lnHM < 0.05, levels=c("TRUE", "FALSE")), 
                      factor(DD==1, levels=c("TRUE", "FALSE")))
tpr.by.DD10sim <- sensitivity(factor(by.disp.lnHM < 0.05, levels=c("TRUE", "FALSE")), 
                      factor(DD==1, levels=c("TRUE", "FALSE")))
tpr.holm.DD10sim <- sensitivity(factor(holm.disp.lnHM < 0.05, levels=c("TRUE", "FALSE")), 
                        factor(DD==1, levels=c("TRUE", "FALSE")))
tpr.diffvar.DD10sim <- sensitivity(factor(q < 0.05, levels=c("TRUE", "FALSE")), 
                           factor(DD==1, levels=c("TRUE", "FALSE")))

pred.MDSeq.DD10sim <- prediction(1-DD10sim$p.disp.zi.MDSeq, DD)
auc.MDSeq.DD10sim <- performance(pred.MDSeq.DD10sim, measure='auc')@y.values[[1]]
pauc.MDSeq.DD10sim <- performance(pred.MDSeq.DD10sim, measure='auc', fpr.stop=0.02)@y.values[[1]]
pred.lnHM.DD10sim <- prediction(1-DD10sim$p.disp.lnHM, DD)
auc.lnHM.DD10sim <- performance(pred.lnHM.DD10sim, measure='auc')@y.values[[1]]
pauc.lnHM.DD10sim <- performance(pred.lnHM.DD10sim, measure='auc', fpr.stop=0.02)@y.values[[1]]
pred.diffvar.DD10sim <- prediction(1-p, DD)
auc.diffvar.DD10sim <- performance(pred.diffvar.DD10sim, measure='auc')@y.values[[1]]
pauc.diffvar.DD10sim <- performance(pred.diffvar.DD10sim, measure='auc', 
                                    fpr.stop=0.02)@y.values[[1]]

c(fdr.MDSeq.DD10sim, fdr.bh.DD10sim, fdr.by.DD10sim, fdr.holm.DD10sim, fdr.diffvar.DD10sim)
c(tpr.MDSeq.DD10sim, tpr.bh.DD10sim, tpr.by.DD10sim, tpr.holm.DD10sim, tpr.diffvar.DD10sim)
c(auc.MDSeq.DD10sim, auc.lnHM.DD10sim, auc.diffvar.DD10sim)
c(pauc.MDSeq.DD10sim, pauc.lnHM.DD10sim, pauc.diffvar.DD10sim)
# Highest AUC and pAUC for HM but no positive calls

roc.MDSeq.DD10sim <- performance(pred.MDSeq.DD10sim, "tpr", "fpr")
roc.lnHM.DD10sim <- performance(pred.lnHM.DD10sim, "tpr", "fpr")
roc.diffvar.DD10sim <- performance(pred.diffvar.DD10sim, "tpr", "fpr")
fp.zi.MDSeq.DD10sim <- pred.MDSeq.DD10sim@fp[[1]]
fp.lnHM.DD10sim <- pred.lnHM.DD10sim@fp[[1]]
fp.diffvar.DD10sim <- pred.diffvar.DD10sim@fp[[1]]
pos.pred.MDSeq.DD10sim <- pred.MDSeq.DD10sim@n.pos.pred[[1]]
pos.pred.lnHM.DD10sim <- pred.lnHM.DD10sim@n.pos.pred[[1]]
pos.pred.diffvar.DD10sim <- pred.diffvar.DD10sim@n.pos.pred[[1]]

par(mfrow=c(2,2), mar=c(2,2,1,1), mgp=c(1.5,1,0))
plot(roc.MDSeq.DD10sim, col='red')
plot(roc.lnHM.DD10sim, col='blue', add=T)
plot(roc.diffvar.DD10sim, col='orange', add=T)
lines(c(0,1), c(0,1), col='grey')
plot(roc.MDSeq.DD10sim, col='red', xlim=c(0,0.05), ylim=c(0,0.3))
plot(roc.lnHM.DD10sim, col='blue', add=T)
plot(roc.diffvar.DD10sim, col='orange', add=T)
lines(c(0,1), c(0,1), col='grey')
plot(pos.pred.MDSeq.DD10sim, fp.zi.MDSeq.DD10sim / pos.pred.MDSeq.DD10sim,
     type='l', col='red', ylim=c(0,1), xlim=c(0,1000))
lines(pos.pred.lnHM.DD10sim, fp.lnHM.DD10sim / pos.pred.lnHM.DD10sim, 
      type='l', col='blue')
lines(pos.pred.diffvar.DD10sim, fp.diffvar.DD10sim / pos.pred.diffvar.DD10sim, 
      type='l', col='orange')
abline(h=0.05, col='grey')
plot(pos.pred.MDSeq.DD10sim, fp.zi.MDSeq.DD10sim / pos.pred.MDSeq.DD10sim,
     type='l', col='red', ylim=c(0,1), xlim=c(0,100))
lines(pos.pred.lnHM.DD10sim, fp.lnHM.DD10sim / pos.pred.lnHM.DD10sim, 
      type='l', col='blue')
lines(pos.pred.diffvar.DD10sim, fp.diffvar.DD10sim / pos.pred.diffvar.DD10sim, 
      type='l', col='orange')
abline(h=0.05, col='grey')
# HM always above MDSeq, both way above diffVar.
# diffVar better than HM only in that it has FDR below 0.05 where HM doesn't 
# make any predictions. Smallest number of discoveries for HM is about 12, and 
# FDR is always far lower than diffVar after that, and lower than MDSeq.


## DD, simulated data, 20 samples per group
DD20sim <- readRDS(here("Results/DD compcodeR data results DESeq norm Sept 2019", 
                        "results.DD20.1.DESeqnorm.rds"))
counts <- DD20sim$data@count.matrix
DD <- DD20sim$DD
group <- factor(c(rep("1",20), rep("2",20)))
design <- model.matrix(~group)
fit <- varFit(counts, design=design, coef=c(1,2))
DV <- topVar(fit, coef=2, number=nrow(counts), sort=F)
p <- DV$P.Value
q <- DV$Adj.P.Value
bh.disp.lnHM <- p.adjust(DD20sim$p.disp.lnHM, method="BH")
by.disp.lnHM <- p.adjust(DD20sim$p.disp.lnHM, method="BY")
holm.disp.lnHM <- p.adjust(DD20sim$p.disp.lnHM, method="holm")

fdr.MDSeq.DD20sim <- 1-precision(factor(DD20sim$q.disp.zi.MDSeq < 0.05, levels=c("TRUE", "FALSE")), 
                                 factor(DD==1, levels=c("TRUE", "FALSE")))
fdr.bh.DD20sim <- 1-precision(factor(bh.disp.lnHM < 0.05, levels=c("TRUE", "FALSE")), 
                              factor(DD==1, levels=c("TRUE", "FALSE")))
fdr.by.DD20sim <- 1-precision(factor(by.disp.lnHM < 0.05, levels=c("TRUE", "FALSE")), 
                              factor(DD==1, levels=c("TRUE", "FALSE")))
fdr.holm.DD20sim <- 1-precision(factor(holm.disp.lnHM < 0.05, levels=c("TRUE", "FALSE")), 
                                factor(DD==1, levels=c("TRUE", "FALSE")))
fdr.diffvar.DD20sim <- 1-precision(factor(q < 0.05, levels=c("TRUE", "FALSE")), 
                                   factor(DD==1, levels=c("TRUE", "FALSE")))

tpr.MDSeq.DD20sim <- sensitivity(factor(DD20sim$q.disp.zi.MDSeq < 0.05, levels=c("TRUE", "FALSE")), 
                                 factor(DD==1, levels=c("TRUE", "FALSE")))
tpr.bh.DD20sim <- sensitivity(factor(bh.disp.lnHM < 0.05, levels=c("TRUE", "FALSE")), 
                              factor(DD==1, levels=c("TRUE", "FALSE")))
tpr.by.DD20sim <- sensitivity(factor(by.disp.lnHM < 0.05, levels=c("TRUE", "FALSE")), 
                              factor(DD==1, levels=c("TRUE", "FALSE")))
tpr.holm.DD20sim <- sensitivity(factor(holm.disp.lnHM < 0.05, levels=c("TRUE", "FALSE")), 
                                factor(DD==1, levels=c("TRUE", "FALSE")))
tpr.diffvar.DD20sim <- sensitivity(factor(q < 0.05, levels=c("TRUE", "FALSE")), 
                                   factor(DD==1, levels=c("TRUE", "FALSE")))

pred.MDSeq.DD20sim <- prediction(1-DD20sim$p.disp.zi.MDSeq, DD)
auc.MDSeq.DD20sim <- performance(pred.MDSeq.DD20sim, measure='auc')@y.values[[1]]
pauc.MDSeq.DD20sim <- performance(pred.MDSeq.DD20sim, measure='auc', fpr.stop=0.02)@y.values[[1]]
pred.lnHM.DD20sim <- prediction(1-DD20sim$p.disp.lnHM, DD)
auc.lnHM.DD20sim <- performance(pred.lnHM.DD20sim, measure='auc')@y.values[[1]]
pauc.lnHM.DD20sim <- performance(pred.lnHM.DD20sim, measure='auc', fpr.stop=0.02)@y.values[[1]]
pred.diffvar.DD20sim <- prediction(1-p, DD)
auc.diffvar.DD20sim <- performance(pred.diffvar.DD20sim, measure='auc')@y.values[[1]]
pauc.diffvar.DD20sim <- performance(pred.diffvar.DD20sim, 
                                    measure='auc', fpr.stop=0.02)@y.values[[1]]

c(fdr.MDSeq.DD20sim, fdr.bh.DD20sim, fdr.by.DD20sim, fdr.holm.DD20sim, fdr.diffvar.DD20sim)
c(tpr.MDSeq.DD20sim, tpr.bh.DD20sim, tpr.by.DD20sim, tpr.holm.DD20sim, tpr.diffvar.DD20sim)
c(auc.MDSeq.DD20sim, auc.lnHM.DD20sim, auc.diffvar.DD20sim)
c(pauc.MDSeq.DD20sim, pauc.lnHM.DD20sim, pauc.diffvar.DD20sim)
# Reasonable FDR for BH, no positive calls for others or diffVar
# Highest AUC and pAUC for HM

roc.MDSeq.DD20sim <- performance(pred.MDSeq.DD20sim, "tpr", "fpr")
roc.lnHM.DD20sim <- performance(pred.lnHM.DD20sim, "tpr", "fpr")
roc.diffvar.DD20sim <- performance(pred.diffvar.DD20sim, "tpr", "fpr")
fp.zi.MDSeq.DD20sim <- pred.MDSeq.DD20sim@fp[[1]]
fp.lnHM.DD20sim <- pred.lnHM.DD20sim@fp[[1]]
fp.diffvar.DD20sim <- pred.diffvar.DD20sim@fp[[1]]
pos.pred.MDSeq.DD20sim <- pred.MDSeq.DD20sim@n.pos.pred[[1]]
pos.pred.lnHM.DD20sim <- pred.lnHM.DD20sim@n.pos.pred[[1]]
pos.pred.diffvar.DD20sim <- pred.diffvar.DD20sim@n.pos.pred[[1]]

par(mfrow=c(2,2), mar=c(2,2,1,1), mgp=c(1.5,1,0))
plot(roc.MDSeq.DD20sim, col='red')
plot(roc.lnHM.DD20sim, col='blue', add=T)
plot(roc.diffvar.DD20sim, col='orange', add=T)
lines(c(0,1), c(0,1), col='grey')
plot(roc.MDSeq.DD20sim, col='red', xlim=c(0,0.05), ylim=c(0,0.4))
plot(roc.lnHM.DD20sim, col='blue', add=T)
plot(roc.diffvar.DD20sim, col='orange', add=T)
lines(c(0,1), c(0,1), col='grey')
plot(pos.pred.MDSeq.DD20sim, fp.zi.MDSeq.DD20sim / pos.pred.MDSeq.DD20sim,
     type='l', col='red', ylim=c(0,1), xlim=c(0,1000))
lines(pos.pred.lnHM.DD20sim, fp.lnHM.DD20sim / pos.pred.lnHM.DD20sim, 
      type='l', col='blue')
lines(pos.pred.diffvar.DD20sim, fp.diffvar.DD20sim / pos.pred.diffvar.DD20sim, 
      type='l', col='orange')
abline(h=0.05, col='grey')
plot(pos.pred.MDSeq.DD20sim, fp.zi.MDSeq.DD20sim / pos.pred.MDSeq.DD20sim,
     type='l', col='red', ylim=c(0,1), xlim=c(0,100))
lines(pos.pred.lnHM.DD20sim, fp.lnHM.DD20sim / pos.pred.lnHM.DD20sim, 
      type='l', col='blue')
lines(pos.pred.diffvar.DD20sim, fp.diffvar.DD20sim / pos.pred.diffvar.DD20sim, 
      type='l', col='orange')
abline(h=0.05, col='grey')
# HM alway has highest curve except a tiny tiny bit at about 0.0001 where MDSeq 
# is higher. Both way above diffVar.
# diffVar has FDR below 0.05 up to about 4 discoveries, others never do, but 
# both way lower than diffVar, and HM lower than MDSeq after about 60 
# discoveries.


## DD, simulated data, 50 samples per group
DD50sim <- readRDS(here("Results/DD compcodeR data results DESeq norm Sept 2019", 
                        "results.DD50.1.DESeqnorm.rds"))
counts <- DD50sim$data@count.matrix
DD <- DD50sim$DD
group <- factor(c(rep("1",50), rep("2",50)))
design <- model.matrix(~group)
fit <- varFit(counts, design=design, coef=c(1,2))
DV <- topVar(fit, coef=2, number=nrow(counts), sort=F)
p <- DV$P.Value
q <- DV$Adj.P.Value
bh.disp.lnHM <- p.adjust(DD50sim$p.disp.lnHM, method="BH")
by.disp.lnHM <- p.adjust(DD50sim$p.disp.lnHM, method="BY")
holm.disp.lnHM <- p.adjust(DD50sim$p.disp.lnHM, method="holm")

fdr.MDSeq.DD50sim <- 1-precision(factor(DD50sim$q.disp.zi.MDSeq < 0.05, levels=c("TRUE", "FALSE")), 
                                 factor(DD==1, levels=c("TRUE", "FALSE")))
fdr.bh.DD50sim <- 1-precision(factor(bh.disp.lnHM < 0.05, levels=c("TRUE", "FALSE")), 
                              factor(DD==1, levels=c("TRUE", "FALSE")))
fdr.by.DD50sim <- 1-precision(factor(by.disp.lnHM < 0.05, levels=c("TRUE", "FALSE")), 
                              factor(DD==1, levels=c("TRUE", "FALSE")))
fdr.holm.DD50sim <- 1-precision(factor(holm.disp.lnHM < 0.05, levels=c("TRUE", "FALSE")), 
                                factor(DD==1, levels=c("TRUE", "FALSE")))
fdr.diffvar.DD50sim <- 1-precision(factor(q < 0.05, levels=c("TRUE", "FALSE")), 
                                   factor(DD==1, levels=c("TRUE", "FALSE")))

tpr.MDSeq.DD50sim <- sensitivity(factor(DD50sim$q.disp.zi.MDSeq < 0.05, levels=c("TRUE", "FALSE")), 
                                 factor(DD==1, levels=c("TRUE", "FALSE")))
tpr.bh.DD50sim <- sensitivity(factor(bh.disp.lnHM < 0.05, levels=c("TRUE", "FALSE")), 
                              factor(DD==1, levels=c("TRUE", "FALSE")))
tpr.by.DD50sim <- sensitivity(factor(by.disp.lnHM < 0.05, levels=c("TRUE", "FALSE")), 
                              factor(DD==1, levels=c("TRUE", "FALSE")))
tpr.holm.DD50sim <- sensitivity(factor(holm.disp.lnHM < 0.05, levels=c("TRUE", "FALSE")), 
                                factor(DD==1, levels=c("TRUE", "FALSE")))
tpr.diffvar.DD50sim <- sensitivity(factor(q < 0.05, levels=c("TRUE", "FALSE")), 
                                   factor(DD==1, levels=c("TRUE", "FALSE")))

pred.MDSeq.DD50sim <- prediction(1-DD50sim$p.disp.zi.MDSeq, DD)
auc.MDSeq.DD50sim <- performance(pred.MDSeq.DD50sim, measure='auc')@y.values[[1]]
pauc.MDSeq.DD50sim <- performance(pred.MDSeq.DD50sim, measure='auc', fpr.stop=0.02)@y.values[[1]]
pred.lnHM.DD50sim <- prediction(1-DD50sim$p.disp.lnHM, DD)
auc.lnHM.DD50sim <- performance(pred.lnHM.DD50sim, measure='auc')@y.values[[1]]
pauc.lnHM.DD50sim <- performance(pred.lnHM.DD50sim, measure='auc', fpr.stop=0.02)@y.values[[1]]
pred.diffvar.DD50sim <- prediction(1-p, DD)
auc.diffvar.DD50sim <- performance(pred.diffvar.DD50sim, measure='auc')@y.values[[1]]
pauc.diffvar.DD50sim <- performance(pred.diffvar.DD50sim, 
                                    measure='auc', fpr.stop=0.02)@y.values[[1]]

c(fdr.MDSeq.DD50sim, fdr.bh.DD50sim, fdr.by.DD50sim, fdr.holm.DD50sim, fdr.diffvar.DD50sim)
c(tpr.MDSeq.DD50sim, tpr.bh.DD50sim, tpr.by.DD50sim, tpr.holm.DD50sim, tpr.diffvar.DD50sim)
c(auc.MDSeq.DD50sim, auc.lnHM.DD50sim, auc.diffvar.DD50sim)
c(pauc.MDSeq.DD50sim, pauc.lnHM.DD50sim, pauc.diffvar.DD50sim)
# FDR < 0.05 for MDSeq, slightly over for BH, no positive calls for others, 
# high for diffVar. BH way higher power than MDSeq.
# Highest AUC and pAUC for HM, but AUC only just higher than MDSeq.

roc.MDSeq.DD50sim <- performance(pred.MDSeq.DD50sim, "tpr", "fpr")
roc.lnHM.DD50sim <- performance(pred.lnHM.DD50sim, "tpr", "fpr")
roc.diffvar.DD50sim <- performance(pred.diffvar.DD50sim, "tpr", "fpr")
fp.zi.MDSeq.DD50sim <- pred.MDSeq.DD50sim@fp[[1]]
fp.lnHM.DD50sim <- pred.lnHM.DD50sim@fp[[1]]
fp.diffvar.DD50sim <- pred.diffvar.DD50sim@fp[[1]]
pos.pred.MDSeq.DD50sim <- pred.MDSeq.DD50sim@n.pos.pred[[1]]
pos.pred.lnHM.DD50sim <- pred.lnHM.DD50sim@n.pos.pred[[1]]
pos.pred.diffvar.DD50sim <- pred.diffvar.DD50sim@n.pos.pred[[1]]

par(mfrow=c(2,2), mar=c(2,2,1,1), mgp=c(1.5,1,0))
plot(roc.MDSeq.DD50sim, col='red')
plot(roc.lnHM.DD50sim, col='blue', add=T)
plot(roc.diffvar.DD50sim, col='orange', add=T)
lines(c(0,1), c(0,1), col='grey')
plot(roc.MDSeq.DD50sim, col='red', xlim=c(0,0.05), ylim=c(0,0.6))
plot(roc.lnHM.DD50sim, col='blue', add=T)
plot(roc.diffvar.DD50sim, col='orange', add=T)
lines(c(0,1), c(0,1), col='grey')
plot(pos.pred.MDSeq.DD50sim, fp.zi.MDSeq.DD50sim / pos.pred.MDSeq.DD50sim,
     type='l', col='red', ylim=c(0,0.8), xlim=c(0,2000))
lines(pos.pred.lnHM.DD50sim, fp.lnHM.DD50sim / pos.pred.lnHM.DD50sim, 
      type='l', col='blue')
lines(pos.pred.diffvar.DD50sim, fp.diffvar.DD50sim / pos.pred.diffvar.DD50sim, 
      type='l', col='orange')
abline(h=0.05, col='grey')
plot(pos.pred.MDSeq.DD50sim, fp.zi.MDSeq.DD50sim / pos.pred.MDSeq.DD50sim,
     type='l', col='red', ylim=c(0,0.2), xlim=c(0,400))
lines(pos.pred.lnHM.DD50sim, fp.lnHM.DD50sim / pos.pred.lnHM.DD50sim, 
      type='l', col='blue')
lines(pos.pred.diffvar.DD50sim, fp.diffvar.DD50sim / pos.pred.diffvar.DD50sim, 
      type='l', col='orange')
abline(h=0.05, col='grey')
# HM always above or matching MDSeq, both always above diffVar.
# diffVar has FDR below 0.05 for very very small number of discoveries then 
# shoots up way above others, MDSeq below 0.05 up to about 220 discoveries, 
# HM up to about 230 but doesn't make any discoveries below about 220, and 
# nearly always below MDSeq after that.


## DD, GTEx blood data, 5 samples per group
DD5blood <- readRDS(here("Results/GTEx blood artificial DD, DE, DEDD results Dec 2019", 
                        "results.blood_10_DD.rds"))
counts <- DD5blood$counts
DD <- DD5blood$DD
group <- factor(c(rep("1",5), rep("2",5)))
design <- model.matrix(~group)
fit <- varFit(counts, design=design, coef=c(1,2))
DV <- topVar(fit, coef=2, number=nrow(counts), sort=F)
p <- DV$P.Value
q <- DV$Adj.P.Value
bh.disp.lnHM <- p.adjust(DD5blood$p.disp.lnHM, method="BH")
by.disp.lnHM <- p.adjust(DD5blood$p.disp.lnHM, method="BY")
holm.disp.lnHM <- p.adjust(DD5blood$p.disp.lnHM, method="holm")

fdr.MDSeq.DD5blood <- 1-precision(factor(DD5blood$q.disp.zi.MDSeq < 0.05, 
                                         levels=c("TRUE", "FALSE")), 
                                 factor(DD==1, levels=c("TRUE", "FALSE")))
fdr.bh.DD5blood <- 1-precision(factor(bh.disp.lnHM < 0.05, levels=c("TRUE", "FALSE")), 
                              factor(DD==1, levels=c("TRUE", "FALSE")))
fdr.by.DD5blood <- 1-precision(factor(by.disp.lnHM < 0.05, levels=c("TRUE", "FALSE")), 
                              factor(DD==1, levels=c("TRUE", "FALSE")))
fdr.holm.DD5blood <- 1-precision(factor(holm.disp.lnHM < 0.05, levels=c("TRUE", "FALSE")), 
                                factor(DD==1, levels=c("TRUE", "FALSE")))
fdr.diffvar.DD5blood <- 1-precision(factor(q < 0.05, levels=c("TRUE", "FALSE")), 
                                   factor(DD==1, levels=c("TRUE", "FALSE")))

tpr.MDSeq.DD5blood <- sensitivity(factor(DD5blood$q.disp.zi.MDSeq < 0.05, 
                                         levels=c("TRUE", "FALSE")), 
                                 factor(DD==1, levels=c("TRUE", "FALSE")))
tpr.bh.DD5blood <- sensitivity(factor(bh.disp.lnHM < 0.05, levels=c("TRUE", "FALSE")), 
                              factor(DD==1, levels=c("TRUE", "FALSE")))
tpr.by.DD5blood <- sensitivity(factor(by.disp.lnHM < 0.05, levels=c("TRUE", "FALSE")), 
                              factor(DD==1, levels=c("TRUE", "FALSE")))
tpr.holm.DD5blood <- sensitivity(factor(holm.disp.lnHM < 0.05, levels=c("TRUE", "FALSE")), 
                                factor(DD==1, levels=c("TRUE", "FALSE")))
tpr.diffvar.DD5blood <- sensitivity(factor(q < 0.05, levels=c("TRUE", "FALSE")), 
                                   factor(DD==1, levels=c("TRUE", "FALSE")))

pred.MDSeq.DD5blood <- prediction(1-DD5blood$p.disp.zi.MDSeq, DD)
auc.MDSeq.DD5blood <- performance(pred.MDSeq.DD5blood, measure='auc')@y.values[[1]]
pauc.MDSeq.DD5blood <- performance(pred.MDSeq.DD5blood, measure='auc', fpr.stop=0.02)@y.values[[1]]
pred.lnHM.DD5blood <- prediction(1-DD5blood$p.disp.lnHM, DD)
auc.lnHM.DD5blood <- performance(pred.lnHM.DD5blood, measure='auc')@y.values[[1]]
pauc.lnHM.DD5blood <- performance(pred.lnHM.DD5blood, measure='auc', fpr.stop=0.02)@y.values[[1]]
pred.diffvar.DD5blood <- prediction(1-p, DD)
auc.diffvar.DD5blood <- performance(pred.diffvar.DD5blood, measure='auc')@y.values[[1]]
pauc.diffvar.DD5blood <- performance(pred.diffvar.DD5blood, measure='auc', 
                                    fpr.stop=0.02)@y.values[[1]]

c(fdr.MDSeq.DD5blood, fdr.bh.DD5blood, fdr.by.DD5blood, fdr.holm.DD5blood, fdr.diffvar.DD5blood)
c(tpr.MDSeq.DD5blood, tpr.bh.DD5blood, tpr.by.DD5blood, tpr.holm.DD5blood, tpr.diffvar.DD5blood)
c(auc.MDSeq.DD5blood, auc.lnHM.DD5blood, auc.diffvar.DD5blood)
c(pauc.MDSeq.DD5blood, pauc.lnHM.DD5blood, pauc.diffvar.DD5blood)
# Only MDSeq makes any positive calls, and has extremely high FDR.
# Highest AUC and pAUC for HM; AUCs quite close, diffVar higher than 
# MDSeq. diffVar pAUC also higher than MDSeq, HM way higher than both.

roc.MDSeq.DD5blood <- performance(pred.MDSeq.DD5blood, "tpr", "fpr")
roc.lnHM.DD5blood <- performance(pred.lnHM.DD5blood, "tpr", "fpr")
roc.diffvar.DD5blood <- performance(pred.diffvar.DD5blood, "tpr", "fpr")
fp.zi.MDSeq.DD5blood <- pred.MDSeq.DD5blood@fp[[1]]
fp.lnHM.DD5blood <- pred.lnHM.DD5blood@fp[[1]]
fp.diffvar.DD5blood <- pred.diffvar.DD5blood@fp[[1]]
pos.pred.MDSeq.DD5blood <- pred.MDSeq.DD5blood@n.pos.pred[[1]]
pos.pred.lnHM.DD5blood <- pred.lnHM.DD5blood@n.pos.pred[[1]]
pos.pred.diffvar.DD5blood <- pred.diffvar.DD5blood@n.pos.pred[[1]]

par(mfrow=c(2,2), mar=c(2,2,1,1), mgp=c(1.5,1,0))
plot(roc.MDSeq.DD5blood, col='red')
plot(roc.lnHM.DD5blood, col='blue', add=T)
plot(roc.diffvar.DD5blood, col='orange', add=T)
lines(c(0,1), c(0,1), col='grey')
plot(roc.MDSeq.DD5blood, col='red', xlim=c(0,0.05), ylim=c(0,0.16))
plot(roc.lnHM.DD5blood, col='blue', add=T)
plot(roc.diffvar.DD5blood, col='orange', add=T)
lines(c(0,1), c(0,1), col='grey')
plot(pos.pred.MDSeq.DD5blood, fp.zi.MDSeq.DD5blood / pos.pred.MDSeq.DD5blood,
     type='l', col='red', ylim=c(0,1), xlim=c(0,1000))
lines(pos.pred.lnHM.DD5blood, fp.lnHM.DD5blood / pos.pred.lnHM.DD5blood, 
      type='l', col='blue')
lines(pos.pred.diffvar.DD5blood, fp.diffvar.DD5blood / pos.pred.diffvar.DD5blood, 
      type='l', col='orange')
abline(h=0.05, col='grey')
plot(pos.pred.MDSeq.DD5blood, fp.zi.MDSeq.DD5blood / pos.pred.MDSeq.DD5blood,
     type='l', col='red', ylim=c(0,1), xlim=c(0,100))
lines(pos.pred.lnHM.DD5blood, fp.lnHM.DD5blood / pos.pred.lnHM.DD5blood, 
      type='l', col='blue')
lines(pos.pred.diffvar.DD5blood, fp.diffvar.DD5blood / pos.pred.diffvar.DD5blood, 
      type='l', col='orange')
abline(h=0.05, col='grey')
# HM always highest, others similar; HM way higher at low end.
# Only HM ever has FDR below 0.05, but only up to 20 discoveries, but 
# always way lower FDR than others, which are both similar.

## DD, GTEx blood data, 10 samples per group
DD10blood <- readRDS(here("Results/GTEx blood artificial DD, DE, DEDD results Dec 2019", 
                         "results.blood_20_DD.rds"))
counts <- DD10blood$counts
DD <- DD10blood$DD
group <- factor(c(rep("1",10), rep("2",10)))
design <- model.matrix(~group)
fit <- varFit(counts, design=design, coef=c(1,2))
DV <- topVar(fit, coef=2, number=nrow(counts), sort=F)
p <- DV$P.Value
q <- DV$Adj.P.Value
bh.disp.lnHM <- p.adjust(DD10blood$p.disp.lnHM, method="BH")
by.disp.lnHM <- p.adjust(DD10blood$p.disp.lnHM, method="BY")
holm.disp.lnHM <- p.adjust(DD10blood$p.disp.lnHM, method="holm")

fdr.MDSeq.DD10blood <- 1-precision(factor(DD10blood$q.disp.zi.MDSeq < 0.05, 
                                         levels=c("TRUE", "FALSE")), 
                                  factor(DD==1, levels=c("TRUE", "FALSE")))
fdr.bh.DD10blood <- 1-precision(factor(bh.disp.lnHM < 0.05, levels=c("TRUE", "FALSE")), 
                               factor(DD==1, levels=c("TRUE", "FALSE")))
fdr.by.DD10blood <- 1-precision(factor(by.disp.lnHM < 0.05, levels=c("TRUE", "FALSE")), 
                               factor(DD==1, levels=c("TRUE", "FALSE")))
fdr.holm.DD10blood <- 1-precision(factor(holm.disp.lnHM < 0.05, levels=c("TRUE", "FALSE")), 
                                 factor(DD==1, levels=c("TRUE", "FALSE")))
fdr.diffvar.DD10blood <- 1-precision(factor(q < 0.05, levels=c("TRUE", "FALSE")), 
                                    factor(DD==1, levels=c("TRUE", "FALSE")))

tpr.MDSeq.DD10blood <- sensitivity(factor(DD10blood$q.disp.zi.MDSeq < 0.05, 
                                         levels=c("TRUE", "FALSE")), 
                                  factor(DD==1, levels=c("TRUE", "FALSE")))
tpr.bh.DD10blood <- sensitivity(factor(bh.disp.lnHM < 0.05, levels=c("TRUE", "FALSE")), 
                               factor(DD==1, levels=c("TRUE", "FALSE")))
tpr.by.DD10blood <- sensitivity(factor(by.disp.lnHM < 0.05, levels=c("TRUE", "FALSE")), 
                               factor(DD==1, levels=c("TRUE", "FALSE")))
tpr.holm.DD10blood <- sensitivity(factor(holm.disp.lnHM < 0.05, levels=c("TRUE", "FALSE")), 
                                 factor(DD==1, levels=c("TRUE", "FALSE")))
tpr.diffvar.DD10blood <- sensitivity(factor(q < 0.05, levels=c("TRUE", "FALSE")), 
                                    factor(DD==1, levels=c("TRUE", "FALSE")))

pred.MDSeq.DD10blood <- prediction(1-DD10blood$p.disp.zi.MDSeq, DD)
auc.MDSeq.DD10blood <- performance(pred.MDSeq.DD10blood, measure='auc')@y.values[[1]]
pauc.MDSeq.DD10blood <- performance(pred.MDSeq.DD10blood, measure='auc', fpr.stop=0.02)@y.values[[1]]
pred.lnHM.DD10blood <- prediction(1-DD10blood$p.disp.lnHM, DD)
auc.lnHM.DD10blood <- performance(pred.lnHM.DD10blood, measure='auc')@y.values[[1]]
pauc.lnHM.DD10blood <- performance(pred.lnHM.DD10blood, measure='auc', fpr.stop=0.02)@y.values[[1]]
pred.diffvar.DD10blood <- prediction(1-p, DD)
auc.diffvar.DD10blood <- performance(pred.diffvar.DD10blood, measure='auc')@y.values[[1]]
pauc.diffvar.DD10blood <- performance(pred.diffvar.DD10blood, measure='auc', 
                                     fpr.stop=0.02)@y.values[[1]]

c(fdr.MDSeq.DD10blood, fdr.bh.DD10blood, fdr.by.DD10blood, fdr.holm.DD10blood, fdr.diffvar.DD10blood)
c(tpr.MDSeq.DD10blood, tpr.bh.DD10blood, tpr.by.DD10blood, tpr.holm.DD10blood, tpr.diffvar.DD10blood)
c(auc.MDSeq.DD10blood, auc.lnHM.DD10blood, auc.diffvar.DD10blood)
c(pauc.MDSeq.DD10blood, pauc.lnHM.DD10blood, pauc.diffvar.DD10blood)
# Only MDSeq, BH, BY make any positive calls, all with extremely high FDR.
# HM has highest AUC but diffVar highest pAUC.

roc.MDSeq.DD10blood <- performance(pred.MDSeq.DD10blood, "tpr", "fpr")
roc.lnHM.DD10blood <- performance(pred.lnHM.DD10blood, "tpr", "fpr")
roc.diffvar.DD10blood <- performance(pred.diffvar.DD10blood, "tpr", "fpr")
fp.zi.MDSeq.DD10blood <- pred.MDSeq.DD10blood@fp[[1]]
fp.lnHM.DD10blood <- pred.lnHM.DD10blood@fp[[1]]
fp.diffvar.DD10blood <- pred.diffvar.DD10blood@fp[[1]]
pos.pred.MDSeq.DD10blood <- pred.MDSeq.DD10blood@n.pos.pred[[1]]
pos.pred.lnHM.DD10blood <- pred.lnHM.DD10blood@n.pos.pred[[1]]
pos.pred.diffvar.DD10blood <- pred.diffvar.DD10blood@n.pos.pred[[1]]

par(mfrow=c(2,2), mar=c(2,2,1,1), mgp=c(1.5,1,0))
plot(roc.MDSeq.DD10blood, col='red')
plot(roc.lnHM.DD10blood, col='blue', add=T)
plot(roc.diffvar.DD10blood, col='orange', add=T)
lines(c(0,1), c(0,1), col='grey')
plot(roc.MDSeq.DD10blood, col='red', xlim=c(0,0.05), ylim=c(0,0.15))
plot(roc.lnHM.DD10blood, col='blue', add=T)
plot(roc.diffvar.DD10blood, col='orange', add=T)
lines(c(0,1), c(0,1), col='grey')
plot(pos.pred.MDSeq.DD10blood, fp.zi.MDSeq.DD10blood / pos.pred.MDSeq.DD10blood,
     type='l', col='red', ylim=c(0,0.9), xlim=c(0,1000))
lines(pos.pred.lnHM.DD10blood, fp.lnHM.DD10blood / pos.pred.lnHM.DD10blood, 
      type='l', col='blue')
lines(pos.pred.diffvar.DD10blood, fp.diffvar.DD10blood / pos.pred.diffvar.DD10blood, 
      type='l', col='orange')
abline(h=0.05, col='grey')
plot(pos.pred.MDSeq.DD10blood, fp.zi.MDSeq.DD10blood / pos.pred.MDSeq.DD10blood,
     type='l', col='red', ylim=c(0,0.9), xlim=c(0,100))
lines(pos.pred.lnHM.DD10blood, fp.lnHM.DD10blood / pos.pred.lnHM.DD10blood, 
      type='l', col='blue')
lines(pos.pred.diffvar.DD10blood, fp.diffvar.DD10blood / pos.pred.diffvar.DD10blood, 
      type='l', col='orange')
abline(h=0.05, col='grey')
# HM always highest except at very low end where diffVar is higher.
# Only diffVar ever has FDR below 0.05, but only for about 2 discoveries. 
# HM starts at over 400 discoveries and is always slightly above diffVar, and 
# both well below MDSeq.


## DD, GTEx blood data, 50 samples per group
DD50blood <- readRDS(here("Results/GTEx blood artificial DD, DE, DEDD results Dec 2019", 
                          "results.blood_100_DD.rds"))
counts <- DD50blood$counts
DD <- DD50blood$DD
group <- factor(c(rep("1",50), rep("2",50)))
design <- model.matrix(~group)
fit <- varFit(counts, design=design, coef=c(1,2))
DV <- topVar(fit, coef=2, number=nrow(counts), sort=F)
p <- DV$P.Value
q <- DV$Adj.P.Value
bh.disp.lnHM <- p.adjust(DD50blood$p.disp.lnHM, method="BH")
by.disp.lnHM <- p.adjust(DD50blood$p.disp.lnHM, method="BY")
holm.disp.lnHM <- p.adjust(DD50blood$p.disp.lnHM, method="holm")

fdr.MDSeq.DD50blood <- 1-precision(factor(DD50blood$q.disp.zi.MDSeq < 0.05, 
                                          levels=c("TRUE", "FALSE")), 
                                   factor(DD==1, levels=c("TRUE", "FALSE")))
fdr.bh.DD50blood <- 1-precision(factor(bh.disp.lnHM < 0.05, levels=c("TRUE", "FALSE")), 
                                factor(DD==1, levels=c("TRUE", "FALSE")))
fdr.by.DD50blood <- 1-precision(factor(by.disp.lnHM < 0.05, levels=c("TRUE", "FALSE")), 
                                factor(DD==1, levels=c("TRUE", "FALSE")))
fdr.holm.DD50blood <- 1-precision(factor(holm.disp.lnHM < 0.05, levels=c("TRUE", "FALSE")), 
                                  factor(DD==1, levels=c("TRUE", "FALSE")))
fdr.diffvar.DD50blood <- 1-precision(factor(q < 0.05, levels=c("TRUE", "FALSE")), 
                                     factor(DD==1, levels=c("TRUE", "FALSE")))

tpr.MDSeq.DD50blood <- sensitivity(factor(DD50blood$q.disp.zi.MDSeq < 0.05, 
                                          levels=c("TRUE", "FALSE")), 
                                   factor(DD==1, levels=c("TRUE", "FALSE")))
tpr.bh.DD50blood <- sensitivity(factor(bh.disp.lnHM < 0.05, levels=c("TRUE", "FALSE")), 
                                factor(DD==1, levels=c("TRUE", "FALSE")))
tpr.by.DD50blood <- sensitivity(factor(by.disp.lnHM < 0.05, levels=c("TRUE", "FALSE")), 
                                factor(DD==1, levels=c("TRUE", "FALSE")))
tpr.holm.DD50blood <- sensitivity(factor(holm.disp.lnHM < 0.05, levels=c("TRUE", "FALSE")), 
                                  factor(DD==1, levels=c("TRUE", "FALSE")))
tpr.diffvar.DD50blood <- sensitivity(factor(q < 0.05, levels=c("TRUE", "FALSE")), 
                                     factor(DD==1, levels=c("TRUE", "FALSE")))

pred.MDSeq.DD50blood <- prediction(1-DD50blood$p.disp.zi.MDSeq, DD)
auc.MDSeq.DD50blood <- performance(pred.MDSeq.DD50blood, measure='auc')@y.values[[1]]
pauc.MDSeq.DD50blood <- performance(pred.MDSeq.DD50blood, measure='auc', fpr.stop=0.02)@y.values[[1]]
pred.lnHM.DD50blood <- prediction(1-DD50blood$p.disp.lnHM, DD)
auc.lnHM.DD50blood <- performance(pred.lnHM.DD50blood, measure='auc')@y.values[[1]]
pauc.lnHM.DD50blood <- performance(pred.lnHM.DD50blood, measure='auc', fpr.stop=0.02)@y.values[[1]]
pred.diffvar.DD50blood <- prediction(1-p, DD)
auc.diffvar.DD50blood <- performance(pred.diffvar.DD50blood, measure='auc')@y.values[[1]]
pauc.diffvar.DD50blood <- performance(pred.diffvar.DD50blood, measure='auc', 
                                      fpr.stop=0.02)@y.values[[1]]

c(fdr.MDSeq.DD50blood, fdr.bh.DD50blood, fdr.by.DD50blood, fdr.holm.DD50blood, fdr.diffvar.DD50blood)
c(tpr.MDSeq.DD50blood, tpr.bh.DD50blood, tpr.by.DD50blood, tpr.holm.DD50blood, tpr.diffvar.DD50blood)
c(auc.MDSeq.DD50blood, auc.lnHM.DD50blood, auc.diffvar.DD50blood)
c(pauc.MDSeq.DD50blood, pauc.lnHM.DD50blood, pauc.diffvar.DD50blood)
# diffVar close to reasonable FDR, others very high; HM better than MDSeq, 
# BY better than BH, no positive calls for Holm. BH, BY higher power than 
# MDSeq despite lower FDRs.
# HM highest AUC, then MDSeq. diffVar highest pAUC, MDSeq way worse than HM.

roc.MDSeq.DD50blood <- performance(pred.MDSeq.DD50blood, "tpr", "fpr")
roc.lnHM.DD50blood <- performance(pred.lnHM.DD50blood, "tpr", "fpr")
roc.diffvar.DD50blood <- performance(pred.diffvar.DD50blood, "tpr", "fpr")
fp.zi.MDSeq.DD50blood <- pred.MDSeq.DD50blood@fp[[1]]
fp.lnHM.DD50blood <- pred.lnHM.DD50blood@fp[[1]]
fp.diffvar.DD50blood <- pred.diffvar.DD50blood@fp[[1]]
pos.pred.MDSeq.DD50blood <- pred.MDSeq.DD50blood@n.pos.pred[[1]]
pos.pred.lnHM.DD50blood <- pred.lnHM.DD50blood@n.pos.pred[[1]]
pos.pred.diffvar.DD50blood <- pred.diffvar.DD50blood@n.pos.pred[[1]]

par(mfrow=c(2,2), mar=c(2,2,1,1), mgp=c(1.5,1,0))
plot(roc.MDSeq.DD50blood, col='red')
plot(roc.lnHM.DD50blood, col='blue', add=T)
plot(roc.diffvar.DD50blood, col='orange', add=T)
lines(c(0,1), c(0,1), col='grey')
plot(roc.MDSeq.DD50blood, col='red', xlim=c(0,0.05), ylim=c(0,0.5))
plot(roc.lnHM.DD50blood, col='blue', add=T)
plot(roc.diffvar.DD50blood, col='orange', add=T)
lines(c(0,1), c(0,1), col='grey')
plot(pos.pred.MDSeq.DD50blood, fp.zi.MDSeq.DD50blood / pos.pred.MDSeq.DD50blood,
     type='l', col='red', ylim=c(0,0.9), xlim=c(0,5000))
lines(pos.pred.lnHM.DD50blood, fp.lnHM.DD50blood / pos.pred.lnHM.DD50blood, 
      type='l', col='blue')
lines(pos.pred.diffvar.DD50blood, fp.diffvar.DD50blood / pos.pred.diffvar.DD50blood, 
      type='l', col='orange')
abline(h=0.05, col='grey')
plot(pos.pred.MDSeq.DD50blood, fp.zi.MDSeq.DD50blood / pos.pred.MDSeq.DD50blood,
     type='l', col='red', ylim=c(0,0.9), xlim=c(0,500))
lines(pos.pred.lnHM.DD50blood, fp.lnHM.DD50blood / pos.pred.lnHM.DD50blood, 
      type='l', col='blue')
lines(pos.pred.diffvar.DD50blood, fp.diffvar.DD50blood / pos.pred.diffvar.DD50blood, 
      type='l', col='orange')
abline(h=0.05, col='grey')
# HM always highest except at low end where diffVar is higher. HM looks to go up in
# a straight line from (0,0) to about (0.03,0.5), possibly because of inability to 
# separate genes at low end, suggesting that longer runs might help.
# Only diffVar ever has FDR below 0.05, which it has up to nearly 300. diffVar is 
# below MDSeq up to about 2500. HM doesn't make any positive calls until about 
# 1500, but after that is always below diffVar.


