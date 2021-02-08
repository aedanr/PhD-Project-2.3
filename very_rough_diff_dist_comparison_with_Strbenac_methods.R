# Strbenac 2016 looked at differential distribution with a focus on feature selection for 
# classification. They tested four different feature selection methods, but chose "DMD" - 
# differences in medians and deviations - as the best. The method can be implemented via the 
# getLocationsAndScales() function from their ClassifyR package. It's not really a method 
# for identifying differential distribution in its own right, being focused purely on 
# feature selection, but best to compare to my method anyway to be sure that there isn't 
# something existing that's better. Would be surprising if it was given it's a very empirical, 
# non-parametric method.
# Below is a very quick comparison which shows that DMD is far inferior to HMMs for one 
# dataset with 20 samples per group and differences in means and dispersions. The differences 
# in AUC and partial AUC are so big that this is enough to convince me that there's no way 
# this method could compete with mine. I don't see a need to do a more robus comparison unless 
# asked by a reviewer.

library(compcodeR)
library(edgeR)
library(DESeq2)
library(ClassifyR)
library(ROCR)
library(caret)
data <- readRDS('Simulated data/DEDD20.1.rds')
groups <- as.factor(data@sample.annotations$condition)
dat.DESeq <- DESeqDataSetFromMatrix(countData=data@count.matrix, colData=data.frame(groups), 
                                    design=~groups)
dat.DESeq <- estimateSizeFactors(dat.DESeq)
nf <- dat.DESeq$sizeFactor
norm <- t(t(data@count.matrix) / nf)
dim(norm)
# 16357, 40
MD1 <- getLocationsAndScales(measurements=norm[, groups == 1], location="median", scale="Qn")
MD2 <- getLocationsAndScales(measurements=norm[, groups == 2], location="median", scale="Qn")
DMD <- abs(MD1$median - MD2$median) + abs(MD1$Qn - MD2$Qn)
length(DMD)
plot(density(DMD))
plot(density(log(DMD)))
range(DMD)
DE <- as.numeric(data@variable.annotations$differential.expression == 1)
DD <- as.numeric(data@variable.annotations$differential.dispersion == 1)
DEDD <- pmax(DE, DD)
mean(DMD[which(DEDD == 0)])
mean(DMD[which(DEDD == 1)])
mean(DMD[which(DE == 0 & DD == 1)])
mean(DMD[which(DE == 1 & DD == 0)])

pred <- prediction(DMD, DEDD)
auc <- performance(pred, measure='auc')@y.values[[1]]
pauc <- performance(pred, measure='auc', fpr.stop=0.05)@y.values[[1]]
c(auc, pauc) # 0.70, 0.0075

pred <- prediction(DMD, DE)
auc <- performance(pred, measure='auc')@y.values[[1]]
pauc <- performance(pred, measure='auc', fpr.stop=0.05)@y.values[[1]]
c(auc, pauc) # 0.77, 0.0103

pred <- prediction(DMD, DD)
auc <- performance(pred, measure='auc')@y.values[[1]]
pauc <- performance(pred, measure='auc', fpr.stop=0.05)@y.values[[1]]
c(auc, pauc) # 0.64, 0.0047

# AUC for differential distribution using DMD 0.70 for 20 samples per group
# Mean for HMMs around 0.86
# pAUC using DMD 0.008, HMMs 0.031

