library(here)
library(RColorBrewer)

## Load data, create colour vector ####

# Blood
folder <- "Results/GTEx blood simulated DE, DD, DEDD results Feb 2020"
for (i in c("DE", "DD", "DEDD")) {
  for (j in c("2_DEDD", "5_DEDD", "10_DEDD", "20_DEDD", "50_DEDD")) {
    assign(paste0(i, ".results.blood_", j), 
           readRDS(here(folder, paste0(i, ".results.blood_", j, ".rds"))))
  }
}
for (i in c("DE", "DEDD")) {
  for (j in c("2_DE", "5_DE", "10_DE", "20_DE", "50_DE")) {
    assign(paste0(i, ".results.blood_", j), 
           readRDS(here(folder, paste0(i, ".results.blood_", j, ".rds"))))
  }
}
for (i in c("DD", "DEDD")) {
  for (j in c("2_DD", "5_DD", "10_DD", "20_DD", "50_DD")) {
    assign(paste0(i, ".results.blood_", j), 
           readRDS(here(folder, paste0(i, ".results.blood_", j, ".rds"))))
  }
}

# Muscle
folder <- "Results/GTEx muscle simulated DE, DD, DEDD results March 2020"
for (i in c("DE", "DD", "DEDD")) {
  for (j in c("2_DEDD", "5_DEDD", "10_DEDD", "20_DEDD", "50_DEDD")) {
    assign(paste0(i, ".results.muscle_", j), 
           readRDS(here(folder, paste0(i, ".results.muscle_", j, ".rds"))))
  }
}
for (i in c("DE", "DEDD")) {
  for (j in c("2_DE", "5_DE", "10_DE", "20_DE", "50_DE")) {
    assign(paste0(i, ".results.muscle_", j), 
           readRDS(here(folder, paste0(i, ".results.muscle_", j, ".rds"))))
  }
}
for (i in c("DD", "DEDD")) {
  for (j in c("2_DD", "5_DD", "10_DD", "20_DD", "50_DD")) {
    assign(paste0(i, ".results.muscle_", j), 
           readRDS(here(folder, paste0(i, ".results.muscle_", j, ".rds"))))
  }
}
rm(i,j,folder)

n <- 23
qual_col_pals = brewer.pal.info[brewer.pal.info$category == "qual" & 
                                  brewer.pal.info$colorblind == T,]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, 
                           rownames(qual_col_pals)))
# pie(rep(1,n), col=col_vector)
col_vector <- col_vector[c(2,4,9,18)]

## Compile into single list for blood and tissue, TMM only ####
DD.results.2_DEDD <- list()
DD.results.5_DEDD <- list()
DD.results.10_DEDD <- list()
DD.results.20_DEDD <- list()
DD.results.50_DEDD <- list()
DE.results.2_DEDD <- list()
DE.results.5_DEDD <- list()
DE.results.10_DEDD <- list()
DE.results.20_DEDD <- list()
DE.results.50_DEDD <- list()
DEDD.results.2_DEDD <- list()
DEDD.results.5_DEDD <- list()
DEDD.results.10_DEDD <- list()
DEDD.results.20_DEDD <- list()
DEDD.results.50_DEDD <- list()
for (i in 1:7) {
  DD.results.2_DEDD[[i]] <- cbind(
    DD.results.blood_2_DEDD[[i]][
      , grep("tmm", names(DD.results.blood_2_DEDD[[i]]))], 
    DD.results.muscle_2_DEDD[[i]][
      , grep("tmm", names(DD.results.muscle_2_DEDD[[i]]))])
  methods <- c("diffVar_blood", "MDSeq.zi_blood", "MDSeq.nozi_blood", 
               "expHM_blood", "expHM.log_blood", "lnHM_blood", "lnHM.log_blood", 
               "diffVar_muscle", "MDSeq.zi_muscle", "MDSeq.nozi_muscle", 
               "expHM_muscle", "expHM.log_muscle", "lnHM_muscle", "lnHM.log_muscle")
  names(DD.results.2_DEDD[[i]]) <- methods
  DD.results.5_DEDD[[i]] <- cbind(
    DD.results.blood_5_DEDD[[i]][
      , grep("tmm", names(DD.results.blood_5_DEDD[[i]]))], 
    DD.results.muscle_5_DEDD[[i]][
      , grep("tmm", names(DD.results.muscle_5_DEDD[[i]]))])
  names(DD.results.5_DEDD[[i]]) <- methods
  DD.results.10_DEDD[[i]] <- cbind(
    DD.results.blood_10_DEDD[[i]][
      , grep("tmm", names(DD.results.blood_10_DEDD[[i]]))], 
    DD.results.muscle_10_DEDD[[i]][
      , grep("tmm", names(DD.results.muscle_10_DEDD[[i]]))])
  names(DD.results.10_DEDD[[i]]) <- methods
  DD.results.20_DEDD[[i]] <- cbind(
    DD.results.blood_20_DEDD[[i]][
      , grep("tmm", names(DD.results.blood_20_DEDD[[i]]))], 
    DD.results.muscle_20_DEDD[[i]][
      , grep("tmm", names(DD.results.muscle_20_DEDD[[i]]))])
  names(DD.results.20_DEDD[[i]]) <- methods
  DD.results.50_DEDD[[i]] <- cbind(
    DD.results.blood_50_DEDD[[i]][
      , grep("tmm", names(DD.results.blood_50_DEDD[[i]]))], 
    DD.results.muscle_50_DEDD[[i]][
      , grep("tmm", names(DD.results.muscle_50_DEDD[[i]]))])
  names(DD.results.50_DEDD[[i]]) <- methods
  DE.results.2_DEDD[[i]] <- cbind(
    DE.results.blood_2_DEDD[[i]][
      , grep("tmm", names(DE.results.blood_2_DEDD[[i]]))], 
    DE.results.muscle_2_DEDD[[i]][
      , grep("tmm", names(DE.results.muscle_2_DEDD[[i]]))])
  if (i %in% c(1,2,6,7)) {
    methods <- c("edgeR.ql_blood", "edgeR.lr_blood", "edgeR.et_blood", 
                 "DESeq2.noif_blood", "DESeq2.if_blood", "voom_blood", 
                 "DSS_blood", "baySeq_blood", 
                 "MDSeq.zi_blood", "MDSeq.nozi_blood", 
                 "expHM_blood", "expHM.log_blood", 
                 "lnHM_blood", "lnHM.log_blood", 
                 "edgeR.ql_muscle", "edgeR.lr_muscle", "edgeR.et_muscle", 
                 "DESeq2.noif_muscle", "DESeq2.if_muscle", "voom_muscle", 
                 "DSS_muscle", "baySeq_muscle", 
                 "MDSeq.zi_muscle", "MDSeq.nozi_muscle", 
                 "expHM_muscle", "expHM.log_muscle", 
                 "lnHM_muscle", "lnHM.log_muscle")
  } else {
    methods <- c("edgeR.ql_blood", "edgeR.lr_blood", "edgeR.et_blood", 
                 "DESeq2.noif_blood", "DESeq2.if_blood", "voom_blood", 
                 "DSS_blood", "DSS.lfdr_blood", "baySeq_blood", 
                 "MDSeq.zi_blood", "MDSeq.nozi_blood", 
                 "expHM_blood", "expHM.log_blood", 
                 "lnHM_blood", "lnHM.log_blood", 
                 "edgeR.ql_muscle", "edgeR.lr_muscle", "edgeR.et_muscle", 
                 "DESeq2.noif_muscle", "DESeq2.if_muscle", "voom_muscle", 
                 "DSS_muscle", "DSS.lfdr_muscle", "baySeq_muscle", 
                 "MDSeq.zi_muscle", "MDSeq.nozi_muscle", 
                 "expHM_muscle", "expHM.log_muscle", 
                 "lnHM_muscle", "lnHM.log_muscle")
  }
  names(DE.results.2_DEDD[[i]]) <- methods
  DE.results.5_DEDD[[i]] <- cbind(
    DE.results.blood_5_DEDD[[i]][
      , grep("tmm", names(DE.results.blood_5_DEDD[[i]]))], 
    DE.results.muscle_5_DEDD[[i]][
      , grep("tmm", names(DE.results.muscle_5_DEDD[[i]]))])
  names(DE.results.5_DEDD[[i]]) <- methods
  DE.results.10_DEDD[[i]] <- cbind(
    DE.results.blood_10_DEDD[[i]][
      , grep("tmm", names(DE.results.blood_10_DEDD[[i]]))], 
    DE.results.muscle_10_DEDD[[i]][
      , grep("tmm", names(DE.results.muscle_10_DEDD[[i]]))])
  names(DE.results.10_DEDD[[i]]) <- methods
  DE.results.20_DEDD[[i]] <- cbind(
    DE.results.blood_20_DEDD[[i]][
      , grep("tmm", names(DE.results.blood_20_DEDD[[i]]))], 
    DE.results.muscle_20_DEDD[[i]][
      , grep("tmm", names(DE.results.muscle_20_DEDD[[i]]))])
  names(DE.results.20_DEDD[[i]]) <- methods
  DE.results.50_DEDD[[i]] <- cbind(
    DE.results.blood_50_DEDD[[i]][
      , grep("tmm", names(DE.results.blood_50_DEDD[[i]]))], 
    DE.results.muscle_50_DEDD[[i]][
      , grep("tmm", names(DE.results.muscle_50_DEDD[[i]]))])
  names(DE.results.50_DEDD[[i]]) <- methods
  DEDD.results.2_DEDD[[i]] <- cbind(
    DEDD.results.blood_2_DEDD[[i]][
      , grep("tmm", names(DEDD.results.blood_2_DEDD[[i]]))], 
    DEDD.results.muscle_2_DEDD[[i]][
      , grep("tmm", names(DEDD.results.muscle_2_DEDD[[i]]))])
  if (i %in% c(1,2,6,7)) {
    methods <- c("diffVar_blood", "expHMM_blood", "lnHMM_blood", 
                 "edgeR_diffVar_blood", "edgeR_MDSeq_blood", "edgeR_lnHM_blood", 
                 "lnHM_diffVar_blood", "lnHM_MDSeq_blood", "lnHM_lnHM_blood", 
                 "diffVar_muscle", "expHMM_muscle", "lnHMM_muscle", 
                 "edgeR_diffVar_muscle", "edgeR_MDSeq_muscle", "edgeR_lnHM_muscle", 
                 "lnHM_diffVar_muscle", "lnHM_MDSeq_muscle", "lnHM_lnHM_muscle")
  } else {
    methods <- c("diffVar_blood", 
                 "expHMM.point5_blood", "expHMM.thr_blood", "expHMM.bfdr_blood", 
                 "lnHMM.point5_blood", "lnHMM.thr_blood", "lnHMM.bfdr_blood", 
                 "edgeR_diffVar_blood", "edgeR_MDSeq_blood", "edgeR_lnHM_blood", 
                 "lnHM_diffVar_blood", "lnHM_MDSeq_blood", "lnHM_lnHM_blood", 
                 "diffVar_muscle", 
                 "expHMM_muscle.point5", "expHMM.thr_muscle", "expHMM.bfdr_muscle", 
                 "lnHMM_muscle.point5", "lnHMM.thr_muscle", "lnHMM.bfdr_muscle", 
                 "edgeR_diffVar_muscle", "edgeR_MDSeq_muscle", "edgeR_lnHM_muscle", 
                 "lnHM_diffVar_muscle", "lnHM_MDSeq_muscle", "lnHM_lnHM_muscle")
    
  }
  names(DEDD.results.2_DEDD[[i]]) <- methods
  DEDD.results.5_DEDD[[i]] <- cbind(
    DEDD.results.blood_5_DEDD[[i]][
      , grep("tmm", names(DEDD.results.blood_5_DEDD[[i]]))], 
    DEDD.results.muscle_5_DEDD[[i]][
      , grep("tmm", names(DEDD.results.muscle_5_DEDD[[i]]))])
  names(DEDD.results.5_DEDD[[i]]) <- methods
  DEDD.results.10_DEDD[[i]] <- cbind(
    DEDD.results.blood_10_DEDD[[i]][
      , grep("tmm", names(DEDD.results.blood_10_DEDD[[i]]))], 
    DEDD.results.muscle_10_DEDD[[i]][
      , grep("tmm", names(DEDD.results.muscle_10_DEDD[[i]]))])
  names(DEDD.results.10_DEDD[[i]]) <- methods
  DEDD.results.20_DEDD[[i]] <- cbind(
    DEDD.results.blood_20_DEDD[[i]][
      , grep("tmm", names(DEDD.results.blood_20_DEDD[[i]]))], 
    DEDD.results.muscle_20_DEDD[[i]][
      , grep("tmm", names(DEDD.results.muscle_20_DEDD[[i]]))])
  names(DEDD.results.20_DEDD[[i]]) <- methods
  DEDD.results.50_DEDD[[i]] <- cbind(
    DEDD.results.blood_50_DEDD[[i]][
      , grep("tmm", names(DEDD.results.blood_50_DEDD[[i]]))], 
    DEDD.results.muscle_50_DEDD[[i]][
      , grep("tmm", names(DEDD.results.muscle_50_DEDD[[i]]))])
  names(DEDD.results.50_DEDD[[i]]) <- methods
}
names(DD.results.2_DEDD) <- c("pauc", "auc", "fpr", "fdr", "tpr", "mean.fdr", 
                              "mean.discoveries")
names(DD.results.5_DEDD) <- c("pauc", "auc", "fpr", "fdr", "tpr", "mean.fdr", 
                              "mean.discoveries")
names(DD.results.10_DEDD) <- c("pauc", "auc", "fpr", "fdr", "tpr", "mean.fdr", 
                               "mean.discoveries")
names(DD.results.20_DEDD) <- c("pauc", "auc", "fpr", "fdr", "tpr", "mean.fdr", 
                               "mean.discoveries")
names(DD.results.50_DEDD) <- c("pauc", "auc", "fpr", "fdr", "tpr", "mean.fdr", 
                               "mean.discoveries")
names(DE.results.2_DEDD) <- c("pauc", "auc", "fpr", "fdr", "tpr", "mean.fdr", 
                              "mean.discoveries")
names(DE.results.5_DEDD) <- c("pauc", "auc", "fpr", "fdr", "tpr", "mean.fdr", 
                              "mean.discoveries")
names(DE.results.10_DEDD) <- c("pauc", "auc", "fpr", "fdr", "tpr", "mean.fdr", 
                               "mean.discoveries")
names(DE.results.20_DEDD) <- c("pauc", "auc", "fpr", "fdr", "tpr", "mean.fdr", 
                               "mean.discoveries")
names(DE.results.50_DEDD) <- c("pauc", "auc", "fpr", "fdr", "tpr", "mean.fdr", 
                               "mean.discoveries")
names(DEDD.results.2_DEDD) <- c("pauc", "auc", "fpr", "fdr", "tpr", "mean.fdr", 
                              "mean.discoveries")
names(DEDD.results.5_DEDD) <- c("pauc", "auc", "fpr", "fdr", "tpr", "mean.fdr", 
                              "mean.discoveries")
names(DEDD.results.10_DEDD) <- c("pauc", "auc", "fpr", "fdr", "tpr", "mean.fdr", 
                               "mean.discoveries")
names(DEDD.results.20_DEDD) <- c("pauc", "auc", "fpr", "fdr", "tpr", "mean.fdr", 
                               "mean.discoveries")
names(DEDD.results.50_DEDD) <- c("pauc", "auc", "fpr", "fdr", "tpr", "mean.fdr", 
                               "mean.discoveries")


############################################
#### Differences in mean and dispersion ####
############################################

################################
### Differential dispersion ####
################################

# Compare diffVar, MDSeq with ZI, log-transformed lnHM
# 5 and 50 samples per group

###########
## AUC ####
methods <- c("diffVar", "MDSeq", "HM")
par(mfrow=c(1,2), mar=c(1,2.5,5.5,1), mgp=c(3,0.7,0))
boxplot(DD.results.5_DEDD$auc[, c(1,2,7,8,9,14)], names=NA, col=col_vector[1:3], 
        pch=20, cex.axis=2, xaxt='n', ylim=c(0.4,1), cex.main=2, 
        main=paste0("5 samples per group\n", 
                    "Blood (left), muscle (right)"))
lines(c(3.5,3.5), c(0.4,0.95), col="lightgrey")
legend("topleft", fill=col_vector[1:3], bty='n', cex=2, ncol=3, legend=methods)
boxplot(DD.results.50_DEDD$auc[, c(1,2,7,8,9,14)], names=NA, col=col_vector[1:3], 
        pch=20, cex.axis=2, xaxt='n', ylim=c(0.4,1), cex.main=2, 
        main=paste0("50 samples per group\n", 
                    "Blood (left), muscle (right)"))
lines(c(3.5,3.5), c(0.4,0.95), col="lightgrey")
legend("topleft", fill=col_vector[1:3], bty='n', cex=2, ncol=3, legend=methods)
rbind(colMeans(DD.results.2_DEDD$auc[, c(1,2,7,8,9,14)]), 
      colMeans(DD.results.5_DEDD$auc[, c(1,2,7,8,9,14)]), 
      colMeans(DD.results.10_DEDD$auc[, c(1,2,7,8,9,14)]), 
      colMeans(DD.results.20_DEDD$auc[, c(1,2,7,8,9,14)]), 
      colMeans(DD.results.50_DEDD$auc[, c(1,2,7,8,9,14)]))

############
## pAUC ####
methods <- c("diffVar", "MDSeq", "HM")
par(mfrow=c(1,2), mar=c(1,2.5,5.5,1), mgp=c(3,0.7,0))
boxplot(DD.results.5_DEDD$pauc[, c(1,2,7,8,9,14)], names=NA, col=col_vector[1:3], 
        pch=20, cex.axis=2, xaxt='n', ylim=c(0,0.04), cex.main=2, 
        main=paste0("5 samples per group\n", 
                    "Blood (left), muscle (right)"))
lines(c(3.5,3.5), c(0,0.035), col="lightgrey")
legend("topleft", fill=col_vector[1:3], bty='n', cex=2, ncol=3, legend=methods)
boxplot(DD.results.50_DEDD$pauc[, c(1,2,7,8,9,14)], names=NA, col=col_vector[1:3], 
        pch=20, cex.axis=2, xaxt='n', ylim=c(0,0.04), cex.main=2, 
        main=paste0("50 samples per group\n", 
                    "Blood (left), muscle (right)"))
lines(c(3.5,3.5), c(0,0.035), col="lightgrey")
legend("topleft", fill=col_vector[1:3], bty='n', cex=2, ncol=3, legend=methods)
rbind(colMeans(DD.results.2_DEDD$pauc[, c(1,2,7,8,9,14)]), 
      colMeans(DD.results.5_DEDD$pauc[, c(1,2,7,8,9,14)]), 
      colMeans(DD.results.10_DEDD$pauc[, c(1,2,7,8,9,14)]), 
      colMeans(DD.results.20_DEDD$pauc[, c(1,2,7,8,9,14)]), 
      colMeans(DD.results.50_DEDD$pauc[, c(1,2,7,8,9,14)]))

###########
## FDR ####
methods <- c("diffVar", "MDSeq", "HM")
par(mfrow=c(1,2), mar=c(1,2.5,5.5,1), mgp=c(3,0.7,0))
boxplot(DD.results.5_DEDD$fdr[, c(1,2,7,8,9,14)], names=NA, col=col_vector[1:3], 
        pch=20, cex.axis=2, xaxt='n', ylim=c(0,1.1), cex.main=2, 
        main=paste0("5 samples per group\n", 
                    "Blood (left), muscle (right)"))
lines(c(3.5,3.5), c(0,1), col="lightgrey"); abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[1:3], bty='n', cex=2, ncol=3, legend=methods)
boxplot(DD.results.50_DEDD$fdr[, c(1,2,7,8,9,14)], names=NA, col=col_vector[1:3], 
        pch=20, cex.axis=2, xaxt='n', ylim=c(0,1.1), cex.main=2, 
        main=paste0("50 samples per group\n", 
                    "Blood (left), muscle (right)"))
lines(c(3.5,3.5), c(0,1), col="lightgrey"); abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[1:3], bty='n', cex=2, ncol=3, legend=methods)
rbind(colMeans(DD.results.2_DEDD$fdr[, c(1,2,7,8,9,14)]), 
      colMeans(DD.results.5_DEDD$fdr[, c(1,2,7,8,9,14)]), 
      colMeans(DD.results.10_DEDD$fdr[, c(1,2,7,8,9,14)]), 
      colMeans(DD.results.20_DEDD$fdr[, c(1,2,7,8,9,14)]), 
      colMeans(DD.results.50_DEDD$fdr[, c(1,2,7,8,9,14)]))

###########
## TPR ####
methods <- c("diffVar", "MDSeq", "HM")
par(mfrow=c(1,2), mar=c(1,2.5,5.5,1), mgp=c(3,0.7,0))
boxplot(DD.results.5_DEDD$tpr[, c(1,2,7,8,9,14)], names=NA, col=col_vector[1:3], 
        pch=20, cex.axis=2, xaxt='n', ylim=c(0,0.7), cex.main=2, 
        main=paste0("TPR, 5 samples per group\n", 
                    "Blood (left), muscle (right)"))
lines(c(3.5,3.5), c(0,0.65), col="lightgrey")
legend("topleft", fill=col_vector[1:3], bty='n', cex=2, ncol=3, legend=methods)
boxplot(DD.results.50_DEDD$tpr[, c(1,2,7,8,9,14)], names=NA, col=col_vector[1:3], 
        pch=20, cex.axis=2, xaxt='n', ylim=c(0,0.7), cex.main=2, 
        main=paste0("TPR, 50 samples per group\n", 
                    "Blood (left), muscle (right)"))
lines(c(3.5,3.5), c(0,0.65), col="lightgrey")
legend("topleft", fill=col_vector[1:3], bty='n', cex=2, ncol=3, legend=methods)
rbind(colMeans(DD.results.2_DEDD$tpr[, c(1,2,7,8,9,14)]), 
      colMeans(DD.results.5_DEDD$tpr[, c(1,2,7,8,9,14)]), 
      colMeans(DD.results.10_DEDD$tpr[, c(1,2,7,8,9,14)]), 
      colMeans(DD.results.20_DEDD$tpr[, c(1,2,7,8,9,14)]), 
      colMeans(DD.results.50_DEDD$tpr[, c(1,2,7,8,9,14)]))

#############################
## False discovery plots ####
methods <- c("diffVar", "MDSeq", "HM")
par(mfrow=c(1,2), mar=c(4.5,5,4.5,1), mgp=c(3,1,0))
plot(DD.results.5_DEDD$mean.discoveries$diffVar_blood, 
     DD.results.5_DEDD$mean.fdr$diffVar_blood, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], lty=1, lwd=2, 
     main=paste0("5 samples per group\n", 
                 "Solid: blood; dashed: muscle"), 
     cex.main=2, xlab="Discoveries", ylab="FDR", cex.lab=1.8, cex.axis=1.5)
lines(DD.results.5_DEDD$mean.discoveries$MDSeq.zi_blood, 
      DD.results.5_DEDD$mean.fdr$MDSeq.zi_blood, 
      lty=1, lwd=2, col=col_vector[2])
lines(DD.results.5_DEDD$mean.discoveries$lnHM.log_blood, 
      DD.results.5_DEDD$mean.fdr$lnHM.log_blood, 
      lty=1, lwd=2, col=col_vector[3])
lines(DD.results.5_DEDD$mean.discoveries$diffVar_muscle, 
      DD.results.5_DEDD$mean.fdr$diffVar_muscle, 
      lty=2, lwd=2, col=col_vector[1])
lines(DD.results.5_DEDD$mean.discoveries$MDSeq.zi_muscle, 
      DD.results.5_DEDD$mean.fdr$MDSeq.zi_muscle, 
      lty=2, lwd=2, col=col_vector[2])
lines(DD.results.5_DEDD$mean.discoveries$lnHM.log_muscle, 
      DD.results.5_DEDD$mean.fdr$lnHM.log_muscle, 
      lty=2, lwd=2, col=col_vector[3])
lines(c(0,2000), c(0.05,0.05), col="lightgrey")
plot(DD.results.50_DEDD$mean.discoveries$diffVar_blood, 
     DD.results.50_DEDD$mean.fdr$diffVar_blood, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], lty=1, lwd=2, 
     main=paste0("50 samples per group\n", 
                 "Solid: blood; dashed: muscle"), 
     cex.main=2, xlab="Discoveries", ylab="FDR", cex.lab=1.8, cex.axis=1.5)
lines(DD.results.50_DEDD$mean.discoveries$MDSeq.zi_blood, 
      DD.results.50_DEDD$mean.fdr$MDSeq.zi_blood, 
      lty=1, lwd=2, col=col_vector[2])
lines(DD.results.50_DEDD$mean.discoveries$lnHM.log_blood, 
      DD.results.50_DEDD$mean.fdr$lnHM.log_blood, 
      lty=1, lwd=2, col=col_vector[3])
lines(DD.results.50_DEDD$mean.discoveries$diffVar_muscle, 
      DD.results.50_DEDD$mean.fdr$diffVar_muscle, 
      lty=2, lwd=2, col=col_vector[1])
lines(DD.results.50_DEDD$mean.discoveries$MDSeq.zi_muscle, 
      DD.results.50_DEDD$mean.fdr$MDSeq.zi_muscle, 
      lty=2, lwd=2, col=col_vector[2])
lines(DD.results.50_DEDD$mean.discoveries$lnHM.log_muscle, 
      DD.results.50_DEDD$mean.fdr$lnHM.log_muscle, 
      lty=2, lwd=2, col=col_vector[3])
lines(c(0,2000), c(0.05,0.05), col="lightgrey")
legend("topright", bty='n', col=col_vector[1:3], lty=1, ncol=1, lwd=2, cex=1.5, 
       legend=methods)


################################
### Differential expression ####
################################

# Compare edgeR QL, DESeq2 with IF, limma-voom, log-transformed lnHM
# 5, 20 and 50 samples per group

###########
## AUC ####
methods <- c("edgeR", "DESeq2", "voom", "HM")
par(mfrow=c(1,3), mar=c(1,2.5,5.5,1), mgp=c(3,0.7,0))
boxplot(DE.results.5_DEDD$auc[, c(1,5,6,14,15,19,20,28)], 
        names=NA, col=col_vector[1:4], 
        pch=20, cex.axis=2, xaxt='n', ylim=c(0.77,1.039), cex.main=2.5, 
        main=paste0("5 samples per group\n", "Blood (left), muscle (right)"))
lines(c(4.5,4.5), c(0.77,1), col="lightgrey")
legend("topleft", fill=col_vector[1:4], bty='n', cex=2.5, ncol=2, legend=methods)
boxplot(DE.results.20_DEDD$auc[, c(1,5,6,14,15,19,20,28)], 
        names=NA, col=col_vector[1:4], 
        pch=20, cex.axis=2, xaxt='n', ylim=c(0.77,1.039), cex.main=2.5, 
        main=paste0("20 samples per group\n", "Blood (left), muscle (right)"))
lines(c(4.5,4.5), c(0.77,1), col="lightgrey")
legend("topleft", fill=col_vector[1:4], bty='n', cex=2.5, ncol=2, legend=methods)
boxplot(DE.results.50_DEDD$auc[, c(1,5,6,14,15,19,20,28)], 
        names=NA, col=col_vector[1:4], 
        pch=20, cex.axis=2, xaxt='n', ylim=c(0.77,1.039), cex.main=2.5, 
        main=paste0("50 samples per group\n", "Blood (left), muscle (right)"))
lines(c(4.5,4.5), c(0.77,1), col="lightgrey")
legend("topleft", fill=col_vector[1:4], bty='n', cex=2.5, ncol=2, legend=methods)
rbind(colMeans(DE.results.2_DEDD$auc), 
      colMeans(DE.results.5_DEDD$auc), 
      colMeans(DE.results.10_DEDD$auc), 
      colMeans(DE.results.20_DEDD$auc), 
      colMeans(DE.results.50_DEDD$auc))

############
## pAUC ####
methods <- c("edgeR", "DESeq2", "voom", "HM")
par(mfrow=c(1,3), mar=c(1,2.5,5.5,1), mgp=c(3,0.7,0))
boxplot(DE.results.5_DEDD$pauc[, c(1,5,6,14,15,19,20,28)], 
        names=NA, col=col_vector[1:4], 
        pch=20, cex.axis=2, xaxt='n', ylim=c(0,0.05), cex.main=2.5, 
        main=paste0("5 samples per group\n", "Blood (left), muscle (right)"))
lines(c(4.5,4.5), c(0.005,0.05), col="lightgrey")
legend("bottomleft", fill=col_vector[1:4], bty='n', cex=2.5, ncol=2, legend=methods)
boxplot(DE.results.20_DEDD$pauc[, c(1,5,6,14,15,19,20,28)], 
        names=NA, col=col_vector[1:4], 
        pch=20, cex.axis=2, xaxt='n', ylim=c(0,0.05), cex.main=2.5, 
        main=paste0("20 samples per group\n", "Blood (left), muscle (right)"))
lines(c(4.5,4.5), c(0.005,0.05), col="lightgrey")
legend("bottomleft", fill=col_vector[1:4], bty='n', cex=2.5, ncol=2, legend=methods)
boxplot(DE.results.50_DEDD$pauc[, c(1,5,6,14,15,19,20,28)], 
        names=NA, col=col_vector[1:4], 
        pch=20, cex.axis=2, xaxt='n', ylim=c(0,0.05), cex.main=2.5, 
        main=paste0("50 samples per group\n", "Blood (left), muscle (right)"))
lines(c(4.5,4.5), c(0.005,0.05), col="lightgrey")
legend("bottomleft", fill=col_vector[1:4], bty='n', cex=2.5, ncol=2, legend=methods)
rbind(colMeans(DE.results.2_DEDD$pauc), 
      colMeans(DE.results.5_DEDD$pauc), 
      colMeans(DE.results.10_DEDD$pauc), 
      colMeans(DE.results.20_DEDD$pauc), 
      colMeans(DE.results.50_DEDD$pauc))

###########
## FDR ####
methods <- c("edgeR", "DESeq2", "voom", "HM")
par(mfrow=c(1,3), mar=c(1,2.5,5.5,1), mgp=c(3,0.7,0))
boxplot(DE.results.5_DEDD$fdr[, c(1,5,6,15,16,20,21,30)], 
        names=NA, col=col_vector[1:4], 
        pch=20, cex.axis=2, xaxt='n', ylim=c(0,0.85), cex.main=2.5, 
        main=paste0("5 samples per group\n", "Blood (left), muscle (right)"))
lines(c(4.5,4.5), c(0,0.75), col="lightgrey"); abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[1:4], bty='n', cex=2.5, ncol=2, legend=methods)
boxplot(DE.results.20_DEDD$fdr[, c(1,5,6,15,16,20,21,30)], 
        names=NA, col=col_vector[1:4], 
        pch=20, cex.axis=2, xaxt='n', ylim=c(0,0.85), cex.main=2.5, 
        main=paste0("20 samples per group\n", "Blood (left), muscle (right)"))
lines(c(4.5,4.5), c(0,0.75), col="lightgrey"); abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[1:4], bty='n', cex=2.5, ncol=2, legend=methods)
boxplot(DE.results.50_DEDD$fdr[, c(1,5,6,15,16,20,21,30)], 
        names=NA, col=col_vector[1:4], 
        pch=20, cex.axis=2, xaxt='n', ylim=c(0,0.85), cex.main=2.5, 
        main=paste0("50 samples per group\n", "Blood (left), muscle (right)"))
lines(c(4.5,4.5), c(0,0.75), col="lightgrey"); abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[1:4], bty='n', cex=2.5, ncol=2, legend=methods)
rbind(colMeans(DE.results.2_DEDD$fdr), 
      colMeans(DE.results.5_DEDD$fdr), 
      colMeans(DE.results.10_DEDD$fdr), 
      colMeans(DE.results.20_DEDD$fdr), 
      colMeans(DE.results.50_DEDD$fdr))

###########
## TPR ####
methods <- c("edgeR", "DESeq2", "voom", "HM")
par(mfrow=c(1,3), mar=c(1,2.5,5.5,1), mgp=c(3,0.7,0))
boxplot(DE.results.5_DEDD$tpr[, c(1,5,6,15,16,20,21,30)], 
        names=NA, col=col_vector[1:4], 
        pch=20, cex.axis=2.5, xaxt='n', ylim=c(0,1.1), cex.main=2.5, 
        main=paste0("5 samples per group\n", "Blood (left), muscle (right)"))
lines(c(4.5,4.5), c(0,1), col="lightgrey")
legend("topleft", fill=col_vector[1:4], bty='n', cex=2, ncol=2, legend=methods)
boxplot(DE.results.20_DEDD$tpr[, c(1,5,6,15,16,20,21,30)], 
        names=NA, col=col_vector[1:4], 
        pch=20, cex.axis=2.5, xaxt='n', ylim=c(0,1.1), cex.main=2.5, 
        main=paste0("20 samples per group\n", "Blood (left), muscle (right)"))
lines(c(4.5,4.5), c(0,1), col="lightgrey")
legend("topleft", fill=col_vector[1:4], bty='n', cex=2, ncol=2, legend=methods)
boxplot(DE.results.50_DEDD$tpr[, c(1,5,6,15,16,20,21,30)], 
        names=NA, col=col_vector[1:4], 
        pch=20, cex.axis=2.5, xaxt='n', ylim=c(0,1.1), cex.main=2.5, 
        main=paste0("50 samples per group\n", "Blood (left), muscle (right)"))
lines(c(4.5,4.5), c(0,1), col="lightgrey")
legend("topleft", fill=col_vector[1:4], bty='n', cex=2, ncol=2, legend=methods)
rbind(colMeans(DE.results.2_DEDD$tpr), 
      colMeans(DE.results.5_DEDD$tpr), 
      colMeans(DE.results.10_DEDD$tpr), 
      colMeans(DE.results.20_DEDD$tpr), 
      colMeans(DE.results.50_DEDD$tpr))

#############################
## False discovery plots ####
methods <- c("edgeR", "DESeq2", "voom", "HM")
par(mfrow=c(1,3), mar=c(4.5,5,5,1), mgp=c(2.7,1,0))
plot(DE.results.5_DEDD$mean.discoveries$edgeR.ql_blood, 
     DE.results.5_DEDD$mean.fdr$edgeR.ql_blood, 
     type='l', ylim=c(0,0.5), xlim=c(0,2000), col=col_vector[1], lty=1, lwd=2, 
     main=paste0("5 samples per group\n", 
                 "Solid: blood; dashed: muscle"), 
     cex.main=2.5, xlab="Discoveries", ylab="FDR", cex.lab=2, cex.axis=1.8)
lines(DE.results.5_DEDD$mean.discoveries$DESeq2.if_blood, 
      DE.results.5_DEDD$mean.fdr$DESeq2.if_blood, lty=1, lwd=2, col=col_vector[2])
lines(DE.results.5_DEDD$mean.discoveries$voom_blood, 
      DE.results.5_DEDD$mean.fdr$voom_blood, lty=1, lwd=2, col=col_vector[3])
lines(DE.results.5_DEDD$mean.discoveries$lnHM.log_blood, 
      DE.results.5_DEDD$mean.fdr$lnHM.log_blood, lty=1, lwd=2, col=col_vector[4])
lines(DE.results.5_DEDD$mean.discoveries$edgeR.ql_muscle, 
      DE.results.5_DEDD$mean.fdr$edgeR.ql_muscle, lty=2, lwd=2, col=col_vector[1])
lines(DE.results.5_DEDD$mean.discoveries$DESeq2.if_muscle, 
      DE.results.5_DEDD$mean.fdr$DESeq2.if_muscle, lty=2, lwd=2, col=col_vector[2])
lines(DE.results.5_DEDD$mean.discoveries$voom_muscle, 
      DE.results.5_DEDD$mean.fdr$voom_muscle, lty=2, lwd=2, col=col_vector[3])
lines(DE.results.5_DEDD$mean.discoveries$lnHM.log_muscle, 
      DE.results.5_DEDD$mean.fdr$lnHM.log_muscle, lty=2, lwd=2, col=col_vector[4])
abline(h=0.05, col="lightgrey")
plot(DE.results.20_DEDD$mean.discoveries$edgeR.ql_blood, 
     DE.results.20_DEDD$mean.fdr$edgeR.ql_blood, 
     type='l', ylim=c(0,0.5), xlim=c(0,2000), col=col_vector[1], lty=1, lwd=2, 
     main=paste0("20 samples per group\n", 
                 "Solid: blood; dashed: muscle"), 
     cex.main=2.5, xlab="Discoveries", ylab="FDR", cex.lab=2, cex.axis=1.8)
lines(DE.results.20_DEDD$mean.discoveries$DESeq2.if_blood, 
      DE.results.20_DEDD$mean.fdr$DESeq2.if_blood, lty=1, lwd=2, col=col_vector[2])
lines(DE.results.20_DEDD$mean.discoveries$voom_blood, 
      DE.results.20_DEDD$mean.fdr$voom_blood, lty=1, lwd=2, col=col_vector[3])
lines(DE.results.20_DEDD$mean.discoveries$lnHM.log_blood, 
      DE.results.20_DEDD$mean.fdr$lnHM.log_blood, lty=1, lwd=2, col=col_vector[4])
lines(DE.results.20_DEDD$mean.discoveries$edgeR.ql_muscle, 
      DE.results.20_DEDD$mean.fdr$edgeR.ql_muscle, lty=2, lwd=2, col=col_vector[1])
lines(DE.results.20_DEDD$mean.discoveries$DESeq2.if_muscle, 
      DE.results.20_DEDD$mean.fdr$DESeq2.if_muscle, lty=2, lwd=2, col=col_vector[2])
lines(DE.results.20_DEDD$mean.discoveries$voom_muscle, 
      DE.results.20_DEDD$mean.fdr$voom_muscle, lty=2, lwd=2, col=col_vector[3])
lines(DE.results.20_DEDD$mean.discoveries$lnHM.log_muscle, 
      DE.results.20_DEDD$mean.fdr$lnHM.log_muscle, lty=2, lwd=2, col=col_vector[4])
abline(h=0.05, col="lightgrey")
plot(DE.results.50_DEDD$mean.discoveries$edgeR.ql_blood, 
     DE.results.50_DEDD$mean.fdr$edgeR.ql_blood, 
     type='l', ylim=c(0,0.5), xlim=c(0,2000), col=col_vector[1], lty=1, lwd=2, 
     main=paste0("50 samples per group\n", 
                 "Solid: blood; dashed: muscle"), 
     cex.main=2.5, xlab="Discoveries", ylab="FDR", cex.lab=2, cex.axis=1.8)
lines(DE.results.50_DEDD$mean.discoveries$DESeq2.if_blood, 
      DE.results.50_DEDD$mean.fdr$DESeq2.if_blood, lty=1, lwd=2, col=col_vector[2])
lines(DE.results.50_DEDD$mean.discoveries$voom_blood, 
      DE.results.50_DEDD$mean.fdr$voom_blood, lty=1, lwd=2, col=col_vector[3])
lines(DE.results.50_DEDD$mean.discoveries$lnHM.log_blood, 
      DE.results.50_DEDD$mean.fdr$lnHM.log_blood, lty=1, lwd=2, col=col_vector[4])
lines(DE.results.50_DEDD$mean.discoveries$edgeR.ql_muscle, 
      DE.results.50_DEDD$mean.fdr$edgeR.ql_muscle, lty=2, lwd=2, col=col_vector[1])
lines(DE.results.50_DEDD$mean.discoveries$DESeq2.if_muscle, 
      DE.results.50_DEDD$mean.fdr$DESeq2.if_muscle, lty=2, lwd=2, col=col_vector[2])
lines(DE.results.50_DEDD$mean.discoveries$voom_muscle, 
      DE.results.50_DEDD$mean.fdr$voom_muscle, lty=2, lwd=2, col=col_vector[3])
lines(DE.results.50_DEDD$mean.discoveries$lnHM.log_muscle, 
      DE.results.50_DEDD$mean.fdr$lnHM.log_muscle, lty=2, lwd=2, col=col_vector[4])
abline(h=0.05, col="lightgrey")
legend("topleft", bty='n', col=col_vector[1:4], lty=1, ncol=1, lwd=2, cex=2, 
       legend=methods)


##################################
### Differential distribution ####
##################################

# Compare diffVar, edgeR combined with log-transformed lnHM, lnHMM with BFDR
# (since posterior threshold doesn't attempt to control FDR so isn't directly 
# comparable with other methods for FDR control)
# 5 and 50 samples per group

###########
## AUC ####
methods <- c("diffVar", "edgeR/HM", "HMM")
par(mfrow=c(1,2), mar=c(1,2.5,5.5,1), mgp=c(3,0.7,0))
boxplot(DEDD.results.5_DEDD$auc[, c(1,6,3,10,15,12)], 
        names=NA, col=col_vector[1:3], 
        pch=20, cex.axis=2, xaxt='n', ylim=c(0.3,1), cex.main=2, 
        main=paste0("5 samples per group\n", "Blood (left), muscle (right)"))
lines(c(3.5,3.5), c(0.4,1), col="lightgrey")
legend("bottomleft", fill=col_vector[1:3], bty='n', cex=1.8, ncol=3, legend=methods)
boxplot(DEDD.results.50_DEDD$auc[, c(1,6,3,10,15,12)], 
        names=NA, col=col_vector[1:3], 
        pch=20, cex.axis=2, xaxt='n', ylim=c(0.3,1), cex.main=2, 
        main=paste0("50 samples per group\n", "Blood (left), muscle (right)"))
lines(c(3.5,3.5), c(0.4,1), col="lightgrey")
legend("bottomleft", fill=col_vector[1:3], bty='n', cex=1.8, ncol=3, legend=methods)
rbind(colMeans(DEDD.results.2_DEDD$auc), 
      colMeans(DEDD.results.5_DEDD$auc), 
      colMeans(DEDD.results.10_DEDD$auc), 
      colMeans(DEDD.results.20_DEDD$auc), 
      colMeans(DEDD.results.50_DEDD$auc))

############
## pAUC ####
methods <- c("diffVar", "edgeR/HM", "HMM")
par(mfrow=c(1,2), mar=c(1,2.5,5.5,1), mgp=c(3,0.7,0))
boxplot(DEDD.results.5_DEDD$pauc[, c(1,6,3,10,15,12)], 
        names=NA, col=col_vector[1:3], 
        pch=20, cex.axis=2, xaxt='n', ylim=c(0,0.048), cex.main=2, 
        main=paste0("5 samples per group\n", "Blood (left), muscle (right)"))
lines(c(3.5,3.5), c(0,0.045), col="lightgrey")
legend("topleft", fill=col_vector[1:3], bty='n', cex=1.8, ncol=3, legend=methods)
boxplot(DEDD.results.50_DEDD$pauc[, c(1,6,3,10,15,12)], 
        names=NA, col=col_vector[1:3], 
        pch=20, cex.axis=2, xaxt='n', ylim=c(0,0.048), cex.main=2, 
        main=paste0("50 samples per group\n", "Blood (left), muscle (right)"))
lines(c(3.5,3.5), c(0,0.045), col="lightgrey")
legend("topleft", fill=col_vector[1:3], bty='n', cex=1.8, ncol=3, legend=methods)
rbind(colMeans(DEDD.results.2_DEDD$pauc), 
      colMeans(DEDD.results.5_DEDD$pauc), 
      colMeans(DEDD.results.10_DEDD$pauc), 
      colMeans(DEDD.results.20_DEDD$pauc), 
      colMeans(DEDD.results.50_DEDD$pauc))

###########
## FDR ####
methods <- c("diffVar", "edgeR/HM", "HMM")
par(mfrow=c(1,2), mar=c(1,2.5,5.5,1), mgp=c(3,0.7,0))
boxplot(DEDD.results.5_DEDD$fdr[, c(1,10,7,14,23,20)], 
        names=NA, col=col_vector[1:3], 
        pch=20, cex.axis=2, xaxt='n', ylim=c(0,1.05), cex.main=2, 
        main=paste0("5 samples per group\n", "Blood (left), muscle (right)"))
lines(c(3.5,3.5), c(0,0.95), col="lightgrey"); abline(h=c(0,0.05), col="lightgrey")
legend("topleft", fill=col_vector[1:3], bty='n', cex=1.8, ncol=3, legend=methods)
boxplot(DEDD.results.50_DEDD$fdr[, c(1,10,7,14,23,20)], 
        names=NA, col=col_vector[1:3], 
        pch=20, cex.axis=2, xaxt='n', ylim=c(0,1.05), cex.main=2, 
        main=paste0("50 samples per group\n", "Blood (left), muscle (right)"))
lines(c(3.5,3.5), c(0,0.95), col="lightgrey"); abline(h=c(0,0.05), col="lightgrey")
legend("topleft", fill=col_vector[1:3], bty='n', cex=1.8, ncol=3, legend=methods)
rbind(colMeans(DEDD.results.2_DEDD$fdr), 
      colMeans(DEDD.results.5_DEDD$fdr), 
      colMeans(DEDD.results.10_DEDD$fdr), 
      colMeans(DEDD.results.20_DEDD$fdr), 
      colMeans(DEDD.results.50_DEDD$fdr))

###########
## TPR ####
methods <- c("diffVar", "edgeR/HM", "HMM")
par(mfrow=c(1,2), mar=c(1,2.5,5.5,1), mgp=c(3,0.7,0))
boxplot(DEDD.results.5_DEDD$tpr[, c(1,10,7,14,23,20)], 
        names=NA, col=col_vector[1:3], 
        pch=20, cex.axis=2, xaxt='n', ylim=c(0,1), cex.main=2, 
        main=paste0("5 samples per group\n", "Blood (left), muscle (right)"))
lines(c(3.5,3.5), c(0,0.9), col="lightgrey")
legend("topleft", fill=col_vector[1:3], bty='n', cex=1.8, ncol=3, legend=methods)
boxplot(DEDD.results.50_DEDD$tpr[, c(1,10,7,14,23,20)], 
        names=NA, col=col_vector[1:3], 
        pch=20, cex.axis=2, xaxt='n', ylim=c(0,1), cex.main=2, 
        main=paste0("50 samples per group\n", "Blood (left), muscle (right)"))
lines(c(3.5,3.5), c(0,0.9), col="lightgrey")
legend("topleft", fill=col_vector[1:3], bty='n', cex=1.8, ncol=3, legend=methods)
rbind(colMeans(DEDD.results.2_DEDD$tpr), 
      colMeans(DEDD.results.5_DEDD$tpr), 
      colMeans(DEDD.results.10_DEDD$tpr), 
      colMeans(DEDD.results.20_DEDD$tpr), 
      colMeans(DEDD.results.50_DEDD$tpr))

#############################
## False discovery plots ####
methods <- c("diffVar", "edgeR/HM", "HMM")
par(mfrow=c(1,2), mar=c(4.5,5,4.5,1), mgp=c(3,1,0))
plot(DEDD.results.5_DEDD$mean.discoveries$diffVar_blood, 
     DEDD.results.5_DEDD$mean.fdr$diffVar_blood, 
     type='l', ylim=c(0,0.7), xlim=c(0,2000), col=col_vector[1], lty=1, lwd=2, 
     main=paste0("5 samples per group\n", 
                 "Solid: blood; dashed: muscle"), 
     cex.main=2, xlab="Discoveries", ylab="FDR", cex.lab=1.8, cex.axis=1.5)
lines(DEDD.results.5_DEDD$mean.discoveries$edgeR_lnHM_blood, 
      DEDD.results.5_DEDD$mean.fdr$edgeR_lnHM_blood, 
      lty=1, lwd=2, col=col_vector[2])
lines(DEDD.results.5_DEDD$mean.discoveries$lnHMM_blood, 
      DEDD.results.5_DEDD$mean.fdr$lnHMM_blood, 
      lty=1, lwd=2, col=col_vector[3])
lines(DEDD.results.5_DEDD$mean.discoveries$diffVar_muscle, 
      DEDD.results.5_DEDD$mean.fdr$diffVar_muscle, 
      lty=2, lwd=2, col=col_vector[1])
lines(DEDD.results.5_DEDD$mean.discoveries$edgeR_lnHM_muscle, 
      DEDD.results.5_DEDD$mean.fdr$edgeR_lnHM_muscle, 
      lty=2, lwd=2, col=col_vector[2])
lines(DEDD.results.5_DEDD$mean.discoveries$lnHMM_muscle, 
      DEDD.results.5_DEDD$mean.fdr$lnHMM_muscle, 
      lty=2, lwd=2, col=col_vector[3])
abline(h=0.05, col="lightgrey")
plot(DEDD.results.50_DEDD$mean.discoveries$diffVar_blood, 
     DEDD.results.50_DEDD$mean.fdr$diffVar_blood, 
     type='l', ylim=c(0,0.7), xlim=c(0,2000), col=col_vector[1], lty=1, lwd=2, 
     main=paste0("50 samples per group\n", 
                 "Solid: blood; dashed: muscle"), 
     cex.main=2, xlab="Discoveries", ylab="FDR", cex.lab=1.8, cex.axis=1.5)
lines(DEDD.results.50_DEDD$mean.discoveries$edgeR_lnHM_blood, 
      DEDD.results.50_DEDD$mean.fdr$edgeR_lnHM_blood, 
      lty=1, lwd=2, col=col_vector[2])
lines(DEDD.results.50_DEDD$mean.discoveries$lnHMM_blood, 
      DEDD.results.50_DEDD$mean.fdr$lnHMM_blood, 
      lty=1, lwd=2, col=col_vector[3])
lines(DEDD.results.50_DEDD$mean.discoveries$diffVar_muscle, 
      DEDD.results.50_DEDD$mean.fdr$diffVar_muscle, 
      lty=2, lwd=2, col=col_vector[1])
lines(DEDD.results.50_DEDD$mean.discoveries$edgeR_lnHM_muscle, 
      DEDD.results.50_DEDD$mean.fdr$edgeR_lnHM_muscle, 
      lty=2, lwd=2, col=col_vector[2])
lines(DEDD.results.50_DEDD$mean.discoveries$lnHMM_muscle, 
      DEDD.results.50_DEDD$mean.fdr$lnHMM_muscle, 
      lty=2, lwd=2, col=col_vector[3])
abline(h=0.05, col="lightgrey")
legend("topright", bty='n', col=col_vector[1:3], lty=1, ncol=1, lwd=2, cex=1.5, 
       legend=methods)


