library(here)
library(compcodeR)
library(limma)
library(edgeR)
library(DESeq2)

#### Number of positive calls by sample set with, without DE ####
npos.DD <- list("2" = data.frame(edgeR.ql = numeric(10), edgeR.lr = numeric(10), edgeR.et = numeric(10),
                                 DESeq2 = numeric(10), voom = numeric(10), n = numeric(10)), 
                "5" = data.frame(edgeR.ql = numeric(10), edgeR.lr = numeric(10), edgeR.et = numeric(10),
                                 DESeq2 = numeric(10), voom = numeric(10), n = numeric(10)), 
                "10" = data.frame(edgeR.ql = numeric(10), edgeR.lr = numeric(10), edgeR.et = numeric(10),
                                  DESeq2 = numeric(10), voom = numeric(10), n = numeric(10)), 
                "20" = data.frame(edgeR.ql = numeric(10), edgeR.lr = numeric(10), edgeR.et = numeric(10),
                                  DESeq2 = numeric(10), voom = numeric(10), n = numeric(10)), 
                "50" = data.frame(edgeR.ql = numeric(10), edgeR.lr = numeric(10), edgeR.et = numeric(10),
                                  DESeq2 = numeric(10), voom = numeric(10), n = numeric(10)))

for (samples.per.cond in c(2,5,10,20,50)) {
  group <- factor(c(rep(1, samples.per.cond), rep(2, samples.per.cond)))
  design <- model.matrix(~group)
  folder <- paste0("recount data/GTEx/blood_", samples.per.cond, "_samples_per_group")
  
  for (i in 1:10) {
    filename <- paste0("blood_", samples.per.cond, "_set", i, "_DD")
    data <- readRDS(here(folder, paste0(filename, ".rds")))
    counts <- data$counts
    npos.DD[[paste0(samples.per.cond)]]$n[i] <- nrow(counts)
    rm(data)
    
    libsizes <- colSums(counts)
    nf <- calcNormFactors(counts, method="TMM")
    els <- nf * libsizes
    sf <- els / exp(mean(log(libsizes)))
    norm.counts <- t(t(counts) / sf)
    
    dat.edgeR <- DGEList(counts=counts, norm.factors=nf, group=group)
    dat.edgeR <- estimateDisp(dat.edgeR, design)
    dat.DESeq2 <- DESeqDataSetFromMatrix(countData=counts, colData=data.frame(group), design=~group)
    sizeFactors(dat.DESeq2) <- sf
    
    fit.edgeR.ql <- glmQLFit(dat.edgeR, design)
    test.edgeR.ql <- glmQLFTest(fit.edgeR.ql)
    fit.edgeR.lr <- glmFit(dat.edgeR, design)
    test.edgeR.lr <- glmLRT(fit.edgeR.lr)
    edgeR.et <- exactTest(dat.edgeR)
    q.edgeR.ql <- topTags(test.edgeR.ql, n=nrow(dat.edgeR$counts), sort='none')$table$FDR
    q.edgeR.lr <- topTags(test.edgeR.lr, n=nrow(dat.edgeR$counts), sort='none')$table$FDR
    q.edgeR.et <- topTags(edgeR.et, n=nrow(dat.edgeR$counts), sort='none')$table$FDR
    npos.DD[[paste0(samples.per.cond)]]$edgeR.ql[i] <- sum(q.edgeR.ql < 0.05)
    npos.DD[[paste0(samples.per.cond)]]$edgeR.lr[i] <- sum(q.edgeR.lr < 0.05)
    npos.DD[[paste0(samples.per.cond)]]$edgeR.et[i] <- sum(q.edgeR.et < 0.05)
    rm(fit.edgeR.ql, test.edgeR.ql, q.edgeR.ql, 
       fit.edgeR.lr, test.edgeR.lr, q.edgeR.lr, 
       edgeR.et, q.edgeR.et)
    
    ## DESeq2
    fit.DESeq2 <- DESeq(dat.DESeq2, minReplicatesForReplace=Inf)
    res.DESeq2 <- results(fit.DESeq2, cooksCutoff=F, alpha=0.05)
    q.DESeq2 <- res.DESeq2$padj
    npos.DD[[paste0(samples.per.cond)]]$DESeq2[i] <- sum(q.DESeq2 < 0.05, na.rm=T)
    rm(dat.DESeq2, fit.DESeq2, res.DESeq2, q.DESeq2)
    
    ## limma-voom
    dat.voom <- voom(dat.edgeR)
    fit.voom <- lmFit(dat.voom, design)
    res.voom <- eBayes(fit.voom)
    q.voom <- topTable(res.voom, number=Inf, sort.by='none')$adj.P.Val
    npos.DD[[paste0(samples.per.cond)]]$voom[i] <- sum(q.voom < 0.05)
    rm(dat.edgeR, dat.voom, fit.voom, res.voom, q.voom)
  }
}

npos # Positive calls for filtered data with no DE
npos.DE # Positive calls for DE data
npos.DE_noDE # Positive calls for DE data without DE genes
npos.DD # Positive calls for DD data

npos$'2'
npos$'5'
npos$'10'
npos$'20'
npos$'50'
npos.DE$'50'
npos.DE_noDE$'50'
# npos.DD$'50'
DE.results.blood_50_DE$fdr[c(1,2,3,5,6)]
DE.results.blood_50_DE$auc[c(1,2,3,5,6)]
DE.results.blood_50_DE$tpr[c(1,2,3,5,6)]


# Some sets give very high false positives - very clearly distinct from others, so must
# be something about those sets that means there are differences when I've assumed 
# there wouldn't be. Nothing to do with introducing DE/DD, as pattern is same in data 
# with no changes. Seems likely there are some differences between samples, possibly 
# some aren't 'normal', or something has gone wrong in the data collection or processing. 
# Can start by looking at which samples are included in each set and see if there are any 
# obviously suspect ones. Indices were sampled without replacement for 2, 5, 10 and 20 
# samples per group, so can only really compare between different sample sizes for those. 
# Indices were sampled with replacement for 50 samples per group, so may be able to 
# compare between sets within 50 samples per group.

indices_2 <- readRDS(file=here("recount data/GTEx/blood_2_samples_per_group", 
                               "indices_from_whole_blood_PAXGene_RIN_at_least_6point9.rds"))
indices_5 <- readRDS(file=here("recount data/GTEx/blood_5_samples_per_group", 
                               "indices_from_whole_blood_PAXGene_RIN_at_least_6point9.rds"))
indices_10 <- readRDS(file=here("recount data/GTEx/blood_10_samples_per_group", 
                               "indices_from_whole_blood_PAXGene_RIN_at_least_6point9.rds"))
indices_20 <- readRDS(file=here("recount data/GTEx/blood_20_samples_per_group", 
                               "indices_from_whole_blood_PAXGene_RIN_at_least_6point9.rds"))
indices_50 <- readRDS(file=here("recount data/GTEx/blood_50_samples_per_group", 
                               "indices_from_whole_blood_PAXGene_RIN_at_least_6point9.rds"))
# From above, excess false positives in:
# 2 samples per group: 2, 5, 6, 10
# 5 samples per group: 1, 3, 6, 7, 8
# 10 samples per group: 5, 8
# 20 samples per group: 1, 3, 8
# 50 samples per group: 1, 4, 5, 7, 10
n2 <- c(indices_2$set2, 
        indices_2$set5, 
        indices_2$set6, 
        indices_2$set10)
n5 <- c(indices_5$set1, 
        indices_5$set3, 
        indices_5$set6, 
        indices_5$set7, 
        indices_5$set8)
n10 <- c(indices_10$set5, 
         indices_10$set8)
n20 <- c(indices_20$set1, 
         indices_20$set3, 
         indices_20$set8)
n50 <- c(indices_50$set1, 
         indices_50$set4, 
         indices_50$set5, 
         indices_50$set7, 
         indices_50$set10)
sum(n2 %in% n5) # 2
sum(n2 %in% n10) # 1
sum(n2 %in% n20) # 7
sum(n5 %in% n10) # 4
sum(n5 %in% n20) # 8
sum(n10 %in% n20) # 9
table(table(n50)) # 105 x 2, 38 x 3, 4 x 4
n2[which(n2 %in% n5)] # 259 91
n2[which(n2 %in% n10)] # 291
n2[which(n2 %in% n20)] # 7 222 217 259 226 273 244
n5[which(n5 %in% n10)] # 236 149 242 206
n5[which(n5 %in% n20)] # 73 114 131 259 203 161 109 59
n10[which(n10 %in% n20)] # 56 72 324 256 147 167 175 379 45
as.numeric(names(table(n50)[which(table(n50) == 4)])) # 17 100 257 401
table(c(n2[which(n2 %in% n5)], 
        n2[which(n2 %in% n10)], 
        n2[which(n2 %in% n20)], 
        n5[which(n5 %in% n10)], 
        n5[which(n5 %in% n20)], 
        n10[which(n10 %in% n20)], 
        as.numeric(names(table(n50)[which(table(n50) == 4)]))))
# Only 259 appears more than once, and it appears 3 times.
table(c(n2[which(n2 %in% n5)], 
        n2[which(n2 %in% n10)], 
        n2[which(n2 %in% n20)], 
        n5[which(n5 %in% n10)], 
        n5[which(n5 %in% n20)], 
        n10[which(n10 %in% n20)], 
        as.numeric(names(table(n50)[which(table(n50) == 4 | table(n50) == 3)]))))
# 161, 206, 217 and 324 appear twice.

259 %in% names(table(unlist(indices_10)))
# 259  in 2, 5 and 20 but not 10
161 %in% names(table(unlist(indices_2)))
161 %in% names(table(unlist(indices_10)))
# 161 in 5 and 20 but not 2 or 10
206 %in% names(table(unlist(indices_2)))
206 %in% names(table(unlist(indices_20)))
# 206 in 5 and 10 but not 2
217 %in% names(table(unlist(indices_2)))
217 %in% names(table(unlist(indices_10)))
# 217 in 5 and 20, not in 2 but in 10
which(unlist(indices_10) == 217)
# 217 in set 4 for 10 samples per group, which doesn't have excess false positives
324 %in% names(table(unlist(indices_2)))
324 %in% names(table(unlist(indices_5)))
# 324 in 10 and 20, not in 2 but in 5
which(unlist(indices_5) == 324)
# 324 in set 10 for 5 samples per group, which doesn't have excess false positives


# Nothing really stands out to suggest that certain samples are causing excess 
# false positives in some sample sets. AUCs and TPRs don't seem unusual for the 
# sets with high FDRs, which also suggests that some sets having true differences 
# between groups isn't the issue. False discovery plots also clearly improve with 
# sample size for all methods and for DD, DE and DEDD detection, so it looks like 
# the issue is in FDR control, or in the absolute but not relative distribution of 
# p-values.


#### Different tissue sources ####
# See 2020-02-26_auc_pauc_fdr_tpr_false_discovery_plots_GTEx_muscle_data_edgeR_DESeq2_voom.R
# edgeR, DESeq2 and voom only, DE only, skeletal muscle
# Trend in FDRs isn't really there for muscle as for blood. FDRs generally decrease with 
# increasing sample size.
# Should do a full analysis on skeletal muscle data.


