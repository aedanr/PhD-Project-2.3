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
col_vector <- col_vector[c(3,4,9,19)]
legend=c("expHM", "expHM log", "lnHM", "lnHM log")

## Compile into single list for blood and tissue, TMM only ####
{
  DD.results.2_DD <- list()
  DD.results.5_DD <- list()
  DD.results.10_DD <- list()
  DD.results.20_DD <- list()
  DD.results.50_DD <- list()
  DEDD.results.2_DD <- list()
  DEDD.results.5_DD <- list()
  DEDD.results.10_DD <- list()
  DEDD.results.20_DD <- list()
  DEDD.results.50_DD <- list()
  for (i in 1:7) {
    DD.results.2_DD[[i]] <- cbind(
      DD.results.blood_2_DD[[i]][
        , grep("tmm", names(DD.results.blood_2_DD[[i]]))], 
      DD.results.muscle_2_DD[[i]][
        , grep("tmm", names(DD.results.muscle_2_DD[[i]]))])
    methods <- c("diffVar_blood", "MDSeq.zi_blood", "MDSeq.nozi_blood", 
                 "expHM_blood", "expHM.log_blood", "lnHM_blood", "lnHM.log_blood", 
                 "diffVar_muscle", "MDSeq.zi_muscle", "MDSeq.nozi_muscle", 
                 "expHM_muscle", "expHM.log_muscle", "lnHM_muscle", "lnHM.log_muscle")
    names(DD.results.2_DD[[i]]) <- methods
    DD.results.5_DD[[i]] <- cbind(
      DD.results.blood_5_DD[[i]][
        , grep("tmm", names(DD.results.blood_5_DD[[i]]))], 
      DD.results.muscle_5_DD[[i]][
        , grep("tmm", names(DD.results.muscle_5_DD[[i]]))])
    names(DD.results.5_DD[[i]]) <- methods
    DD.results.10_DD[[i]] <- cbind(
      DD.results.blood_10_DD[[i]][
        , grep("tmm", names(DD.results.blood_10_DD[[i]]))], 
      DD.results.muscle_10_DD[[i]][
        , grep("tmm", names(DD.results.muscle_10_DD[[i]]))])
    names(DD.results.10_DD[[i]]) <- methods
    DD.results.20_DD[[i]] <- cbind(
      DD.results.blood_20_DD[[i]][
        , grep("tmm", names(DD.results.blood_20_DD[[i]]))], 
      DD.results.muscle_20_DD[[i]][
        , grep("tmm", names(DD.results.muscle_20_DD[[i]]))])
    names(DD.results.20_DD[[i]]) <- methods
    DD.results.50_DD[[i]] <- cbind(
      DD.results.blood_50_DD[[i]][
        , grep("tmm", names(DD.results.blood_50_DD[[i]]))], 
      DD.results.muscle_50_DD[[i]][
        , grep("tmm", names(DD.results.muscle_50_DD[[i]]))])
    names(DD.results.50_DD[[i]]) <- methods
    DEDD.results.2_DD[[i]] <- cbind(
      DEDD.results.blood_2_DD[[i]][
        , grep("tmm", names(DEDD.results.blood_2_DD[[i]]))], 
      DEDD.results.muscle_2_DD[[i]][
        , grep("tmm", names(DEDD.results.muscle_2_DD[[i]]))])
    if (i %in% c(1,2,6,7)) {
      methods <- c("diffVar_blood", "expHMM_blood", "lnHMM_blood", 
                   "diffVar_muscle", "expHMM_muscle", "lnHMM_muscle")
    } else {
      methods <- c("diffVar_blood", 
                   "expHMM.point5_blood", "expHMM.thr_blood", "expHMM.bfdr_blood", 
                   "lnHMM.point5_blood", "lnHMM.thr_blood", "lnHMM.bfdr_blood", 
                   "diffVar_muscle", 
                   "expHMM_muscle.point5", "expHMM.thr_muscle", "expHMM.bfdr_muscle", 
                   "lnHMM_muscle.point5", "lnHMM.thr_muscle", "lnHMM.bfdr_muscle")
      
    }
    names(DEDD.results.2_DD[[i]]) <- methods
    DEDD.results.5_DD[[i]] <- cbind(
      DEDD.results.blood_5_DD[[i]][
        , grep("tmm", names(DEDD.results.blood_5_DD[[i]]))], 
      DEDD.results.muscle_5_DD[[i]][
        , grep("tmm", names(DEDD.results.muscle_5_DD[[i]]))])
    names(DEDD.results.5_DD[[i]]) <- methods
    DEDD.results.10_DD[[i]] <- cbind(
      DEDD.results.blood_10_DD[[i]][
        , grep("tmm", names(DEDD.results.blood_10_DD[[i]]))], 
      DEDD.results.muscle_10_DD[[i]][
        , grep("tmm", names(DEDD.results.muscle_10_DD[[i]]))])
    names(DEDD.results.10_DD[[i]]) <- methods
    DEDD.results.20_DD[[i]] <- cbind(
      DEDD.results.blood_20_DD[[i]][
        , grep("tmm", names(DEDD.results.blood_20_DD[[i]]))], 
      DEDD.results.muscle_20_DD[[i]][
        , grep("tmm", names(DEDD.results.muscle_20_DD[[i]]))])
    names(DEDD.results.20_DD[[i]]) <- methods
    DEDD.results.50_DD[[i]] <- cbind(
      DEDD.results.blood_50_DD[[i]][
        , grep("tmm", names(DEDD.results.blood_50_DD[[i]]))], 
      DEDD.results.muscle_50_DD[[i]][
        , grep("tmm", names(DEDD.results.muscle_50_DD[[i]]))])
    names(DEDD.results.50_DD[[i]]) <- methods
  }
  names(DD.results.2_DD) <- c("pauc", "auc", "fpr", "fdr", "tpr", "mean.fdr", 
                              "mean.discoveries")
  names(DD.results.5_DD) <- c("pauc", "auc", "fpr", "fdr", "tpr", "mean.fdr", 
                              "mean.discoveries")
  names(DD.results.10_DD) <- c("pauc", "auc", "fpr", "fdr", "tpr", "mean.fdr", 
                               "mean.discoveries")
  names(DD.results.20_DD) <- c("pauc", "auc", "fpr", "fdr", "tpr", "mean.fdr", 
                               "mean.discoveries")
  names(DD.results.50_DD) <- c("pauc", "auc", "fpr", "fdr", "tpr", "mean.fdr", 
                               "mean.discoveries")
  names(DEDD.results.2_DD) <- c("pauc", "auc", "fpr", "fdr", "tpr", "mean.fdr", 
                                "mean.discoveries")
  names(DEDD.results.5_DD) <- c("pauc", "auc", "fpr", "fdr", "tpr", "mean.fdr", 
                                "mean.discoveries")
  names(DEDD.results.10_DD) <- c("pauc", "auc", "fpr", "fdr", "tpr", "mean.fdr", 
                                 "mean.discoveries")
  names(DEDD.results.20_DD) <- c("pauc", "auc", "fpr", "fdr", "tpr", "mean.fdr", 
                                 "mean.discoveries")
  names(DEDD.results.50_DD) <- c("pauc", "auc", "fpr", "fdr", "tpr", "mean.fdr", 
                                 "mean.discoveries")
}

{
  DE.results.2_DE <- list()
  DE.results.5_DE <- list()
  DE.results.10_DE <- list()
  DE.results.20_DE <- list()
  DE.results.50_DE <- list()
  DEDD.results.2_DE <- list()
  DEDD.results.5_DE <- list()
  DEDD.results.10_DE <- list()
  DEDD.results.20_DE <- list()
  DEDD.results.50_DE <- list()
  for (i in 1:7) {
    DE.results.2_DE[[i]] <- cbind(
      DE.results.blood_2_DE[[i]][
        , grep("tmm", names(DE.results.blood_2_DE[[i]]))], 
      DE.results.muscle_2_DE[[i]][
        , grep("tmm", names(DE.results.muscle_2_DE[[i]]))])
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
    names(DE.results.2_DE[[i]]) <- methods
    DE.results.5_DE[[i]] <- cbind(
      DE.results.blood_5_DE[[i]][
        , grep("tmm", names(DE.results.blood_5_DE[[i]]))], 
      DE.results.muscle_5_DE[[i]][
        , grep("tmm", names(DE.results.muscle_5_DE[[i]]))])
    names(DE.results.5_DE[[i]]) <- methods
    DE.results.10_DE[[i]] <- cbind(
      DE.results.blood_10_DE[[i]][
        , grep("tmm", names(DE.results.blood_10_DE[[i]]))], 
      DE.results.muscle_10_DE[[i]][
        , grep("tmm", names(DE.results.muscle_10_DE[[i]]))])
    names(DE.results.10_DE[[i]]) <- methods
    DE.results.20_DE[[i]] <- cbind(
      DE.results.blood_20_DE[[i]][
        , grep("tmm", names(DE.results.blood_20_DE[[i]]))], 
      DE.results.muscle_20_DE[[i]][
        , grep("tmm", names(DE.results.muscle_20_DE[[i]]))])
    names(DE.results.20_DE[[i]]) <- methods
    DE.results.50_DE[[i]] <- cbind(
      DE.results.blood_50_DE[[i]][
        , grep("tmm", names(DE.results.blood_50_DE[[i]]))], 
      DE.results.muscle_50_DE[[i]][
        , grep("tmm", names(DE.results.muscle_50_DE[[i]]))])
    names(DE.results.50_DE[[i]]) <- methods
    DEDD.results.2_DE[[i]] <- cbind(
      DEDD.results.blood_2_DE[[i]][
        , grep("tmm", names(DEDD.results.blood_2_DE[[i]]))], 
      DEDD.results.muscle_2_DE[[i]][
        , grep("tmm", names(DEDD.results.muscle_2_DE[[i]]))])
    if (i %in% c(1,2,6,7)) {
      methods <- c("expHMM_blood", "lnHMM_blood", 
                   "expHMM_muscle", "lnHMM_muscle")
    } else {
      methods <- c("expHMM.point5_blood", "expHMM.thr_blood", "expHMM.bfdr_blood", 
                   "lnHMM.point5_blood", "lnHMM.thr_blood", "lnHMM.bfdr_blood", 
                   "expHMM_muscle.point5", "expHMM.thr_muscle", "expHMM.bfdr_muscle", 
                   "lnHMM_muscle.point5", "lnHMM.thr_muscle", "lnHMM.bfdr_muscle")
      
    }
    names(DEDD.results.2_DE[[i]]) <- methods
    DEDD.results.5_DE[[i]] <- cbind(
      DEDD.results.blood_5_DE[[i]][
        , grep("tmm", names(DEDD.results.blood_5_DE[[i]]))], 
      DEDD.results.muscle_5_DE[[i]][
        , grep("tmm", names(DEDD.results.muscle_5_DE[[i]]))])
    names(DEDD.results.5_DE[[i]]) <- methods
    DEDD.results.10_DE[[i]] <- cbind(
      DEDD.results.blood_10_DE[[i]][
        , grep("tmm", names(DEDD.results.blood_10_DE[[i]]))], 
      DEDD.results.muscle_10_DE[[i]][
        , grep("tmm", names(DEDD.results.muscle_10_DE[[i]]))])
    names(DEDD.results.10_DE[[i]]) <- methods
    DEDD.results.20_DE[[i]] <- cbind(
      DEDD.results.blood_20_DE[[i]][
        , grep("tmm", names(DEDD.results.blood_20_DE[[i]]))], 
      DEDD.results.muscle_20_DE[[i]][
        , grep("tmm", names(DEDD.results.muscle_20_DE[[i]]))])
    names(DEDD.results.20_DE[[i]]) <- methods
    DEDD.results.50_DE[[i]] <- cbind(
      DEDD.results.blood_50_DE[[i]][
        , grep("tmm", names(DEDD.results.blood_50_DE[[i]]))], 
      DEDD.results.muscle_50_DE[[i]][
        , grep("tmm", names(DEDD.results.muscle_50_DE[[i]]))])
    names(DEDD.results.50_DE[[i]]) <- methods
  }
  names(DE.results.2_DE) <- c("pauc", "auc", "fpr", "fdr", "tpr", "mean.fdr", 
                              "mean.discoveries")
  names(DE.results.5_DE) <- c("pauc", "auc", "fpr", "fdr", "tpr", "mean.fdr", 
                              "mean.discoveries")
  names(DE.results.10_DE) <- c("pauc", "auc", "fpr", "fdr", "tpr", "mean.fdr", 
                               "mean.discoveries")
  names(DE.results.20_DE) <- c("pauc", "auc", "fpr", "fdr", "tpr", "mean.fdr", 
                               "mean.discoveries")
  names(DE.results.50_DE) <- c("pauc", "auc", "fpr", "fdr", "tpr", "mean.fdr", 
                               "mean.discoveries")
  names(DEDD.results.2_DE) <- c("pauc", "auc", "fpr", "fdr", "tpr", "mean.fdr", 
                                "mean.discoveries")
  names(DEDD.results.5_DE) <- c("pauc", "auc", "fpr", "fdr", "tpr", "mean.fdr", 
                                "mean.discoveries")
  names(DEDD.results.10_DE) <- c("pauc", "auc", "fpr", "fdr", "tpr", "mean.fdr", 
                                 "mean.discoveries")
  names(DEDD.results.20_DE) <- c("pauc", "auc", "fpr", "fdr", "tpr", "mean.fdr", 
                                 "mean.discoveries")
  names(DEDD.results.50_DE) <- c("pauc", "auc", "fpr", "fdr", "tpr", "mean.fdr", 
                                 "mean.discoveries")
}

{
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
}

rm(i, methods, n, qual_col_pals)
rm(list=ls()[c(grep("blood", ls()), grep("muscle", ls()))])

## Remove all except HM/HMM results ####
{
  {
    DD.results.2_DD$auc <- DD.results.2_DD$auc[, grep("HM", names(DD.results.2_DD$auc))]
    DD.results.2_DD$pauc <- DD.results.2_DD$pauc[, grep("HM", names(DD.results.2_DD$pauc))]
    DD.results.2_DD$fdr <- DD.results.2_DD$fdr[, grep("HM", names(DD.results.2_DD$fdr))]
    DD.results.2_DD$tpr <- DD.results.2_DD$tpr[, grep("HM", names(DD.results.2_DD$tpr))]
    DD.results.2_DD$mean.fdr <- DD.results.2_DD$mean.fdr[, grep("HM", names(DD.results.2_DD$mean.fdr))]
    DD.results.2_DD$mean.discoveries <- DD.results.2_DD$mean.discoveries[, grep("HM", names(DD.results.2_DD$mean.discoveries))]
    DD.results.5_DD$auc <- DD.results.5_DD$auc[, grep("HM", names(DD.results.5_DD$auc))]
    DD.results.5_DD$pauc <- DD.results.5_DD$pauc[, grep("HM", names(DD.results.5_DD$pauc))]
    DD.results.5_DD$fdr <- DD.results.5_DD$fdr[, grep("HM", names(DD.results.5_DD$fdr))]
    DD.results.5_DD$tpr <- DD.results.5_DD$tpr[, grep("HM", names(DD.results.5_DD$tpr))]
    DD.results.5_DD$mean.fdr <- DD.results.5_DD$mean.fdr[, grep("HM", names(DD.results.5_DD$mean.fdr))]
    DD.results.5_DD$mean.discoveries <- DD.results.5_DD$mean.discoveries[, grep("HM", names(DD.results.5_DD$mean.discoveries))]
    DD.results.10_DD$auc <- DD.results.10_DD$auc[, grep("HM", names(DD.results.10_DD$auc))]
    DD.results.10_DD$pauc <- DD.results.10_DD$pauc[, grep("HM", names(DD.results.10_DD$pauc))]
    DD.results.10_DD$fdr <- DD.results.10_DD$fdr[, grep("HM", names(DD.results.10_DD$fdr))]
    DD.results.10_DD$tpr <- DD.results.10_DD$tpr[, grep("HM", names(DD.results.10_DD$tpr))]
    DD.results.10_DD$mean.fdr <- DD.results.10_DD$mean.fdr[, grep("HM", names(DD.results.10_DD$mean.fdr))]
    DD.results.10_DD$mean.discoveries <- DD.results.10_DD$mean.discoveries[, grep("HM", names(DD.results.10_DD$mean.discoveries))]
    DD.results.20_DD$auc <- DD.results.20_DD$auc[, grep("HM", names(DD.results.20_DD$auc))]
    DD.results.20_DD$pauc <- DD.results.20_DD$pauc[, grep("HM", names(DD.results.20_DD$pauc))]
    DD.results.20_DD$fdr <- DD.results.20_DD$fdr[, grep("HM", names(DD.results.20_DD$fdr))]
    DD.results.20_DD$tpr <- DD.results.20_DD$tpr[, grep("HM", names(DD.results.20_DD$tpr))]
    DD.results.20_DD$mean.fdr <- DD.results.20_DD$mean.fdr[, grep("HM", names(DD.results.20_DD$mean.fdr))]
    DD.results.20_DD$mean.discoveries <- DD.results.20_DD$mean.discoveries[, grep("HM", names(DD.results.20_DD$mean.discoveries))]
    DD.results.50_DD$auc <- DD.results.50_DD$auc[, grep("HM", names(DD.results.50_DD$auc))]
    DD.results.50_DD$pauc <- DD.results.50_DD$pauc[, grep("HM", names(DD.results.50_DD$pauc))]
    DD.results.50_DD$fdr <- DD.results.50_DD$fdr[, grep("HM", names(DD.results.50_DD$fdr))]
    DD.results.50_DD$tpr <- DD.results.50_DD$tpr[, grep("HM", names(DD.results.50_DD$tpr))]
    DD.results.50_DD$mean.fdr <- DD.results.50_DD$mean.fdr[, grep("HM", names(DD.results.50_DD$mean.fdr))]
    DD.results.50_DD$mean.discoveries <- DD.results.50_DD$mean.discoveries[, grep("HM", names(DD.results.50_DD$mean.discoveries))]
    DEDD.results.2_DD$auc <- DEDD.results.2_DD$auc[, grep("HM", names(DEDD.results.2_DD$auc))]
    DEDD.results.2_DD$pauc <- DEDD.results.2_DD$pauc[, grep("HM", names(DEDD.results.2_DD$pauc))]
    DEDD.results.2_DD$fdr <- DEDD.results.2_DD$fdr[, grep("HM", names(DEDD.results.2_DD$fdr))]
    DEDD.results.2_DD$tpr <- DEDD.results.2_DD$tpr[, grep("HM", names(DEDD.results.2_DD$tpr))]
    DEDD.results.2_DD$mean.fdr <- DEDD.results.2_DD$mean.fdr[, grep("HM", names(DEDD.results.2_DD$mean.fdr))]
    DEDD.results.2_DD$mean.discoveries <- DEDD.results.2_DD$mean.discoveries[
      , grep("HM", names(DEDD.results.2_DD$mean.discoveries))
    ]
    DEDD.results.5_DD$auc <- DEDD.results.5_DD$auc[, grep("HM", names(DEDD.results.5_DD$auc))]
    DEDD.results.5_DD$pauc <- DEDD.results.5_DD$pauc[, grep("HM", names(DEDD.results.5_DD$pauc))]
    DEDD.results.5_DD$fdr <- DEDD.results.5_DD$fdr[, grep("HM", names(DEDD.results.5_DD$fdr))]
    DEDD.results.5_DD$tpr <- DEDD.results.5_DD$tpr[, grep("HM", names(DEDD.results.5_DD$tpr))]
    DEDD.results.5_DD$mean.fdr <- DEDD.results.5_DD$mean.fdr[, grep("HM", names(DEDD.results.5_DD$mean.fdr))]
    DEDD.results.5_DD$mean.discoveries <- DEDD.results.5_DD$mean.discoveries[
      , grep("HM", names(DEDD.results.5_DD$mean.discoveries))]
    DEDD.results.10_DD$auc <- DEDD.results.10_DD$auc[, grep("HM", names(DEDD.results.10_DD$auc))
                                                     ]
    DEDD.results.10_DD$pauc <- DEDD.results.10_DD$pauc[, grep("HM", names(DEDD.results.10_DD$pauc))]
    DEDD.results.10_DD$fdr <- DEDD.results.10_DD$fdr[, grep("HM", names(DEDD.results.10_DD$fdr))]
    DEDD.results.10_DD$tpr <- DEDD.results.10_DD$tpr[, grep("HM", names(DEDD.results.10_DD$tpr))]
    DEDD.results.10_DD$mean.fdr <- DEDD.results.10_DD$mean.fdr[, grep("HM", names(DEDD.results.10_DD$mean.fdr))]
    DEDD.results.10_DD$mean.discoveries <- DEDD.results.10_DD$mean.discoveries[
      , grep("HM", names(DEDD.results.10_DD$mean.discoveries))
      ]
    DEDD.results.20_DD$auc <- DEDD.results.20_DD$auc[, grep("HM", names(DEDD.results.20_DD$auc))]
    DEDD.results.20_DD$pauc <- DEDD.results.20_DD$pauc[, grep("HM", names(DEDD.results.20_DD$pauc))]
    DEDD.results.20_DD$fdr <- DEDD.results.20_DD$fdr[, grep("HM", names(DEDD.results.20_DD$fdr))]
    DEDD.results.20_DD$tpr <- DEDD.results.20_DD$tpr[, grep("HM", names(DEDD.results.20_DD$tpr))]
    DEDD.results.20_DD$mean.fdr <- DEDD.results.20_DD$mean.fdr[, grep("HM", names(DEDD.results.20_DD$mean.fdr))]
    DEDD.results.20_DD$mean.discoveries <- DEDD.results.20_DD$mean.discoveries[
      , grep("HM", names(DEDD.results.20_DD$mean.discoveries))
      ]
    DEDD.results.50_DD$auc <- DEDD.results.50_DD$auc[, grep("HM", names(DEDD.results.50_DD$auc))]
    DEDD.results.50_DD$pauc <- DEDD.results.50_DD$pauc[, grep("HM", names(DEDD.results.50_DD$pauc))]
    DEDD.results.50_DD$fdr <- DEDD.results.50_DD$fdr[, grep("HM", names(DEDD.results.50_DD$fdr))]
    DEDD.results.50_DD$tpr <- DEDD.results.50_DD$tpr[, grep("HM", names(DEDD.results.50_DD$tpr))]
    DEDD.results.50_DD$mean.fdr <- DEDD.results.50_DD$mean.fdr[, grep("HM", names(DEDD.results.50_DD$mean.fdr))]
    DEDD.results.50_DD$mean.discoveries <- DEDD.results.50_DD$mean.discoveries[
      , grep("HM", names(DEDD.results.50_DD$mean.discoveries))
      ]
  }
  
  {
    DE.results.2_DE$auc <- DE.results.2_DE$auc[, grep("HM", names(DE.results.2_DE$auc))]
    DE.results.2_DE$pauc <- DE.results.2_DE$pauc[, grep("HM", names(DE.results.2_DE$pauc))]
    DE.results.2_DE$fdr <- DE.results.2_DE$fdr[, grep("HM", names(DE.results.2_DE$fdr))]
    DE.results.2_DE$tpr <- DE.results.2_DE$tpr[, grep("HM", names(DE.results.2_DE$tpr))]
    DE.results.2_DE$mean.fdr <- DE.results.2_DE$mean.fdr[, grep("HM", names(DE.results.2_DE$mean.fdr))]
    DE.results.2_DE$mean.discoveries <- DE.results.2_DE$mean.discoveries[, grep("HM", names(DE.results.2_DE$mean.discoveries))]
    DE.results.5_DE$auc <- DE.results.5_DE$auc[, grep("HM", names(DE.results.5_DE$auc))]
    DE.results.5_DE$pauc <- DE.results.5_DE$pauc[, grep("HM", names(DE.results.5_DE$pauc))]
    DE.results.5_DE$fdr <- DE.results.5_DE$fdr[, grep("HM", names(DE.results.5_DE$fdr))]
    DE.results.5_DE$tpr <- DE.results.5_DE$tpr[, grep("HM", names(DE.results.5_DE$tpr))]
    DE.results.5_DE$mean.fdr <- DE.results.5_DE$mean.fdr[, grep("HM", names(DE.results.5_DE$mean.fdr))]
    DE.results.5_DE$mean.discoveries <- DE.results.5_DE$mean.discoveries[, grep("HM", names(DE.results.5_DE$mean.discoveries))]
    DE.results.10_DE$auc <- DE.results.10_DE$auc[, grep("HM", names(DE.results.10_DE$auc))]
    DE.results.10_DE$pauc <- DE.results.10_DE$pauc[, grep("HM", names(DE.results.10_DE$pauc))]
    DE.results.10_DE$fdr <- DE.results.10_DE$fdr[, grep("HM", names(DE.results.10_DE$fdr))]
    DE.results.10_DE$tpr <- DE.results.10_DE$tpr[, grep("HM", names(DE.results.10_DE$tpr))]
    DE.results.10_DE$mean.fdr <- DE.results.10_DE$mean.fdr[, grep("HM", names(DE.results.10_DE$mean.fdr))]
    DE.results.10_DE$mean.discoveries <- DE.results.10_DE$mean.discoveries[, grep("HM", names(DE.results.10_DE$mean.discoveries))]
    DE.results.20_DE$auc <- DE.results.20_DE$auc[, grep("HM", names(DE.results.20_DE$auc))]
    DE.results.20_DE$pauc <- DE.results.20_DE$pauc[, grep("HM", names(DE.results.20_DE$pauc))]
    DE.results.20_DE$fdr <- DE.results.20_DE$fdr[, grep("HM", names(DE.results.20_DE$fdr))]
    DE.results.20_DE$tpr <- DE.results.20_DE$tpr[, grep("HM", names(DE.results.20_DE$tpr))]
    DE.results.20_DE$mean.fdr <- DE.results.20_DE$mean.fdr[, grep("HM", names(DE.results.20_DE$mean.fdr))]
    DE.results.20_DE$mean.discoveries <- DE.results.20_DE$mean.discoveries[, grep("HM", names(DE.results.20_DE$mean.discoveries))]
    DE.results.50_DE$auc <- DE.results.50_DE$auc[, grep("HM", names(DE.results.50_DE$auc))]
    DE.results.50_DE$pauc <- DE.results.50_DE$pauc[, grep("HM", names(DE.results.50_DE$pauc))]
    DE.results.50_DE$fdr <- DE.results.50_DE$fdr[, grep("HM", names(DE.results.50_DE$fdr))]
    DE.results.50_DE$tpr <- DE.results.50_DE$tpr[, grep("HM", names(DE.results.50_DE$tpr))]
    DE.results.50_DE$mean.fdr <- DE.results.50_DE$mean.fdr[, grep("HM", names(DE.results.50_DE$mean.fdr))]
    DE.results.50_DE$mean.discoveries <- DE.results.50_DE$mean.discoveries[, grep("HM", names(DE.results.50_DE$mean.discoveries))]
    DEDD.results.2_DE$auc <- DEDD.results.2_DE$auc[, grep("HM", names(DEDD.results.2_DE$auc))]
    DEDD.results.2_DE$pauc <- DEDD.results.2_DE$pauc[, grep("HM", names(DEDD.results.2_DE$pauc))]
    DEDD.results.2_DE$fdr <- DEDD.results.2_DE$fdr[, grep("HM", names(DEDD.results.2_DE$fdr))]
    DEDD.results.2_DE$tpr <- DEDD.results.2_DE$tpr[, grep("HM", names(DEDD.results.2_DE$tpr))]
    DEDD.results.2_DE$mean.fdr <- DEDD.results.2_DE$mean.fdr[, grep("HM", names(DEDD.results.2_DE$mean.fdr))]
    DEDD.results.2_DE$mean.discoveries <- DEDD.results.2_DE$mean.discoveries[
      , grep("HM", names(DEDD.results.2_DE$mean.discoveries))
      ]
    DEDD.results.5_DE$auc <- DEDD.results.5_DE$auc[, grep("HM", names(DEDD.results.5_DE$auc))]
    DEDD.results.5_DE$pauc <- DEDD.results.5_DE$pauc[, grep("HM", names(DEDD.results.5_DE$pauc))]
    DEDD.results.5_DE$fdr <- DEDD.results.5_DE$fdr[, grep("HM", names(DEDD.results.5_DE$fdr))]
    DEDD.results.5_DE$tpr <- DEDD.results.5_DE$tpr[, grep("HM", names(DEDD.results.5_DE$tpr))]
    DEDD.results.5_DE$mean.fdr <- DEDD.results.5_DE$mean.fdr[, grep("HM", names(DEDD.results.5_DE$mean.fdr))]
    DEDD.results.5_DE$mean.discoveries <- DEDD.results.5_DE$mean.discoveries[
      , grep("HM", names(DEDD.results.5_DE$mean.discoveries))]
    DEDD.results.10_DE$auc <- DEDD.results.10_DE$auc[, grep("HM", names(DEDD.results.10_DE$auc))
                                                     ]
    DEDD.results.10_DE$pauc <- DEDD.results.10_DE$pauc[, grep("HM", names(DEDD.results.10_DE$pauc))]
    DEDD.results.10_DE$fdr <- DEDD.results.10_DE$fdr[, grep("HM", names(DEDD.results.10_DE$fdr))]
    DEDD.results.10_DE$tpr <- DEDD.results.10_DE$tpr[, grep("HM", names(DEDD.results.10_DE$tpr))]
    DEDD.results.10_DE$mean.fdr <- DEDD.results.10_DE$mean.fdr[, grep("HM", names(DEDD.results.10_DE$mean.fdr))]
    DEDD.results.10_DE$mean.discoveries <- DEDD.results.10_DE$mean.discoveries[
      , grep("HM", names(DEDD.results.10_DE$mean.discoveries))
      ]
    DEDD.results.20_DE$auc <- DEDD.results.20_DE$auc[, grep("HM", names(DEDD.results.20_DE$auc))]
    DEDD.results.20_DE$pauc <- DEDD.results.20_DE$pauc[, grep("HM", names(DEDD.results.20_DE$pauc))]
    DEDD.results.20_DE$fdr <- DEDD.results.20_DE$fdr[, grep("HM", names(DEDD.results.20_DE$fdr))]
    DEDD.results.20_DE$tpr <- DEDD.results.20_DE$tpr[, grep("HM", names(DEDD.results.20_DE$tpr))]
    DEDD.results.20_DE$mean.fdr <- DEDD.results.20_DE$mean.fdr[, grep("HM", names(DEDD.results.20_DE$mean.fdr))]
    DEDD.results.20_DE$mean.discoveries <- DEDD.results.20_DE$mean.discoveries[
      , grep("HM", names(DEDD.results.20_DE$mean.discoveries))
      ]
    DEDD.results.50_DE$auc <- DEDD.results.50_DE$auc[, grep("HM", names(DEDD.results.50_DE$auc))]
    DEDD.results.50_DE$pauc <- DEDD.results.50_DE$pauc[, grep("HM", names(DEDD.results.50_DE$pauc))]
    DEDD.results.50_DE$fdr <- DEDD.results.50_DE$fdr[, grep("HM", names(DEDD.results.50_DE$fdr))]
    DEDD.results.50_DE$tpr <- DEDD.results.50_DE$tpr[, grep("HM", names(DEDD.results.50_DE$tpr))]
    DEDD.results.50_DE$mean.fdr <- DEDD.results.50_DE$mean.fdr[, grep("HM", names(DEDD.results.50_DE$mean.fdr))]
    DEDD.results.50_DE$mean.discoveries <- DEDD.results.50_DE$mean.discoveries[
      , grep("HM", names(DEDD.results.50_DE$mean.discoveries))
      ]
  }
  
  {
    DD.results.2_DEDD$auc <- DD.results.2_DEDD$auc[, grep("HM", names(DD.results.2_DEDD$auc))]
    DD.results.2_DEDD$pauc <- DD.results.2_DEDD$pauc[, grep("HM", names(DD.results.2_DEDD$pauc))]
    DD.results.2_DEDD$fdr <- DD.results.2_DEDD$fdr[, grep("HM", names(DD.results.2_DEDD$fdr))]
    DD.results.2_DEDD$tpr <- DD.results.2_DEDD$tpr[, grep("HM", names(DD.results.2_DEDD$tpr))]
    DD.results.2_DEDD$mean.fdr <- DD.results.2_DEDD$mean.fdr[, grep("HM", names(DD.results.2_DEDD$mean.fdr))]
    DD.results.2_DEDD$mean.discoveries <- DD.results.2_DEDD$mean.discoveries[
      , grep("HM", names(DD.results.2_DEDD$mean.discoveries))
      ]
    DD.results.5_DEDD$auc <- DD.results.5_DEDD$auc[, grep("HM", names(DD.results.5_DEDD$auc))]
    DD.results.5_DEDD$pauc <- DD.results.5_DEDD$pauc[, grep("HM", names(DD.results.5_DEDD$pauc))]
    DD.results.5_DEDD$fdr <- DD.results.5_DEDD$fdr[, grep("HM", names(DD.results.5_DEDD$fdr))]
    DD.results.5_DEDD$tpr <- DD.results.5_DEDD$tpr[, grep("HM", names(DD.results.5_DEDD$tpr))]
    DD.results.5_DEDD$mean.fdr <- DD.results.5_DEDD$mean.fdr[, grep("HM", names(DD.results.5_DEDD$mean.fdr))]
    DD.results.5_DEDD$mean.discoveries <- DD.results.5_DEDD$mean.discoveries[
      , grep("HM", names(DD.results.5_DEDD$mean.discoveries))
      ]
    DD.results.10_DEDD$auc <- DD.results.10_DEDD$auc[, grep("HM", names(DD.results.10_DEDD$auc))]
    DD.results.10_DEDD$pauc <- DD.results.10_DEDD$pauc[, grep("HM", names(DD.results.10_DEDD$pauc))]
    DD.results.10_DEDD$fdr <- DD.results.10_DEDD$fdr[, grep("HM", names(DD.results.10_DEDD$fdr))]
    DD.results.10_DEDD$tpr <- DD.results.10_DEDD$tpr[, grep("HM", names(DD.results.10_DEDD$tpr))]
    DD.results.10_DEDD$mean.fdr <- DD.results.10_DEDD$mean.fdr[, grep("HM", names(DD.results.10_DEDD$mean.fdr))]
    DD.results.10_DEDD$mean.discoveries <- DD.results.10_DEDD$mean.discoveries[
      , grep("HM", names(DD.results.10_DEDD$mean.discoveries))
      ]
    DD.results.20_DEDD$auc <- DD.results.20_DEDD$auc[, grep("HM", names(DD.results.20_DEDD$auc))]
    DD.results.20_DEDD$pauc <- DD.results.20_DEDD$pauc[, grep("HM", names(DD.results.20_DEDD$pauc))]
    DD.results.20_DEDD$fdr <- DD.results.20_DEDD$fdr[, grep("HM", names(DD.results.20_DEDD$fdr))]
    DD.results.20_DEDD$tpr <- DD.results.20_DEDD$tpr[, grep("HM", names(DD.results.20_DEDD$tpr))]
    DD.results.20_DEDD$mean.fdr <- DD.results.20_DEDD$mean.fdr[, grep("HM", names(DD.results.20_DEDD$mean.fdr))]
    DD.results.20_DEDD$mean.discoveries <- DD.results.20_DEDD$mean.discoveries[
      , grep("HM", names(DD.results.20_DEDD$mean.discoveries))]
    DD.results.50_DEDD$auc <- DD.results.50_DEDD$auc[, grep("HM", names(DD.results.50_DEDD$auc))
                                                     ]
    DD.results.50_DEDD$pauc <- DD.results.50_DEDD$pauc[, grep("HM", names(DD.results.50_DEDD$pauc))]
    DD.results.50_DEDD$fdr <- DD.results.50_DEDD$fdr[, grep("HM", names(DD.results.50_DEDD$fdr))]
    DD.results.50_DEDD$tpr <- DD.results.50_DEDD$tpr[, grep("HM", names(DD.results.50_DEDD$tpr))]
    DD.results.50_DEDD$mean.fdr <- DD.results.50_DEDD$mean.fdr[, grep("HM", names(DD.results.50_DEDD$mean.fdr))]
    DD.results.50_DEDD$mean.discoveries <- DD.results.50_DEDD$mean.discoveries[
      , grep("HM", names(DD.results.50_DEDD$mean.discoveries))]
    DE.results.2_DEDD$auc <- DE.results.2_DEDD$auc[, grep("HM", names(DE.results.2_DEDD$auc))
                                                   ]
    DE.results.2_DEDD$pauc <- DE.results.2_DEDD$pauc[, grep("HM", names(DE.results.2_DEDD$pauc))]
    DE.results.2_DEDD$fdr <- DE.results.2_DEDD$fdr[, grep("HM", names(DE.results.2_DEDD$fdr))]
    DE.results.2_DEDD$tpr <- DE.results.2_DEDD$tpr[, grep("HM", names(DE.results.2_DEDD$tpr))]
    DE.results.2_DEDD$mean.fdr <- DE.results.2_DEDD$mean.fdr[, grep("HM", names(DE.results.2_DEDD$mean.fdr))]
    DE.results.2_DEDD$mean.discoveries <- DE.results.2_DEDD$mean.discoveries[
      , grep("HM", names(DE.results.2_DEDD$mean.discoveries))
      ]
    DE.results.5_DEDD$auc <- DE.results.5_DEDD$auc[, grep("HM", names(DE.results.5_DEDD$auc))]
    DE.results.5_DEDD$pauc <- DE.results.5_DEDD$pauc[, grep("HM", names(DE.results.5_DEDD$pauc))]
    DE.results.5_DEDD$fdr <- DE.results.5_DEDD$fdr[, grep("HM", names(DE.results.5_DEDD$fdr))]
    DE.results.5_DEDD$tpr <- DE.results.5_DEDD$tpr[, grep("HM", names(DE.results.5_DEDD$tpr))]
    DE.results.5_DEDD$mean.fdr <- DE.results.5_DEDD$mean.fdr[, grep("HM", names(DE.results.5_DEDD$mean.fdr))]
    DE.results.5_DEDD$mean.discoveries <- DE.results.5_DEDD$mean.discoveries[
      , grep("HM", names(DE.results.5_DEDD$mean.discoveries))
      ]
    DE.results.10_DEDD$auc <- DE.results.10_DEDD$auc[, grep("HM", names(DE.results.10_DEDD$auc))]
    DE.results.10_DEDD$pauc <- DE.results.10_DEDD$pauc[, grep("HM", names(DE.results.10_DEDD$pauc))]
    DE.results.10_DEDD$fdr <- DE.results.10_DEDD$fdr[, grep("HM", names(DE.results.10_DEDD$fdr))]
    DE.results.10_DEDD$tpr <- DE.results.10_DEDD$tpr[, grep("HM", names(DE.results.10_DEDD$tpr))]
    DE.results.10_DEDD$mean.fdr <- DE.results.10_DEDD$mean.fdr[, grep("HM", names(DE.results.10_DEDD$mean.fdr))]
    DE.results.10_DEDD$mean.discoveries <- DE.results.10_DEDD$mean.discoveries[
      , grep("HM", names(DE.results.10_DEDD$mean.discoveries))
      ]
    DE.results.20_DEDD$auc <- DE.results.20_DEDD$auc[, grep("HM", names(DE.results.20_DEDD$auc))]
    DE.results.20_DEDD$pauc <- DE.results.20_DEDD$pauc[, grep("HM", names(DE.results.20_DEDD$pauc))]
    DE.results.20_DEDD$fdr <- DE.results.20_DEDD$fdr[, grep("HM", names(DE.results.20_DEDD$fdr))]
    DE.results.20_DEDD$tpr <- DE.results.20_DEDD$tpr[, grep("HM", names(DE.results.20_DEDD$tpr))]
    DE.results.20_DEDD$mean.fdr <- DE.results.20_DEDD$mean.fdr[, grep("HM", names(DE.results.20_DEDD$mean.fdr))]
    DE.results.20_DEDD$mean.discoveries <- DE.results.20_DEDD$mean.discoveries[
      , grep("HM", names(DE.results.20_DEDD$mean.discoveries))
      ]
    DE.results.50_DEDD$auc <- DE.results.50_DEDD$auc[, grep("HM", names(DE.results.50_DEDD$auc))]
    DE.results.50_DEDD$pauc <- DE.results.50_DEDD$pauc[, grep("HM", names(DE.results.50_DEDD$pauc))]
    DE.results.50_DEDD$fdr <- DE.results.50_DEDD$fdr[, grep("HM", names(DE.results.50_DEDD$fdr))]
    DE.results.50_DEDD$tpr <- DE.results.50_DEDD$tpr[, grep("HM", names(DE.results.50_DEDD$tpr))]
    DE.results.50_DEDD$mean.fdr <- DE.results.50_DEDD$mean.fdr[, grep("HM", names(DE.results.50_DEDD$mean.fdr))]
    DE.results.50_DEDD$mean.discoveries <- DE.results.50_DEDD$mean.discoveries[
      , grep("HM", names(DE.results.50_DEDD$mean.discoveries))
      ]
    DEDD.results.2_DEDD$auc <- DEDD.results.2_DEDD$auc[, grep("HMM", names(DEDD.results.2_DEDD$auc))]
    DEDD.results.2_DEDD$pauc <- DEDD.results.2_DEDD$pauc[, grep("HMM", names(DEDD.results.2_DEDD$pauc))]
    DEDD.results.2_DEDD$fdr <- DEDD.results.2_DEDD$fdr[, grep("HMM", names(DEDD.results.2_DEDD$fdr))]
    DEDD.results.2_DEDD$tpr <- DEDD.results.2_DEDD$tpr[, grep("HMM", names(DEDD.results.2_DEDD$tpr))]
    DEDD.results.2_DEDD$mean.fdr <- DEDD.results.2_DEDD$mean.fdr[, grep("HMM", names(DEDD.results.2_DEDD$mean.fdr))]
    DEDD.results.2_DEDD$mean.discoveries <- DEDD.results.2_DEDD$mean.discoveries[
      , grep("HMM", names(DEDD.results.2_DEDD$mean.discoveries))
      ]
    DEDD.results.5_DEDD$auc <- DEDD.results.5_DEDD$auc[, grep("HMM", names(DEDD.results.5_DEDD$auc))]
    DEDD.results.5_DEDD$pauc <- DEDD.results.5_DEDD$pauc[, grep("HMM", names(DEDD.results.5_DEDD$pauc))]
    DEDD.results.5_DEDD$fdr <- DEDD.results.5_DEDD$fdr[, grep("HMM", names(DEDD.results.5_DEDD$fdr))]
    DEDD.results.5_DEDD$tpr <- DEDD.results.5_DEDD$tpr[, grep("HMM", names(DEDD.results.5_DEDD$tpr))]
    DEDD.results.5_DEDD$mean.fdr <- DEDD.results.5_DEDD$mean.fdr[, grep("HMM", names(DEDD.results.5_DEDD$mean.fdr))]
    DEDD.results.5_DEDD$mean.discoveries <- DEDD.results.5_DEDD$mean.discoveries[
      , grep("HMM", names(DEDD.results.5_DEDD$mean.discoveries))
      ]
    DEDD.results.10_DEDD$auc <- DEDD.results.10_DEDD$auc[, grep("HMM", names(DEDD.results.10_DEDD$auc))]
    DEDD.results.10_DEDD$pauc <- DEDD.results.10_DEDD$pauc[, grep("HMM", names(DEDD.results.10_DEDD$pauc))]
    DEDD.results.10_DEDD$fdr <- DEDD.results.10_DEDD$fdr[, grep("HMM", names(DEDD.results.10_DEDD$fdr))]
    DEDD.results.10_DEDD$tpr <- DEDD.results.10_DEDD$tpr[, grep("HMM", names(DEDD.results.10_DEDD$tpr))]
    DEDD.results.10_DEDD$mean.fdr <- DEDD.results.10_DEDD$mean.fdr[, grep("HMM", names(DEDD.results.10_DEDD$mean.fdr))]
    DEDD.results.10_DEDD$mean.discoveries <- DEDD.results.10_DEDD$mean.discoveries[
      , grep("HMM", names(DEDD.results.10_DEDD$mean.discoveries))
      ]
    DEDD.results.20_DEDD$auc <- DEDD.results.20_DEDD$auc[, grep("HMM", names(DEDD.results.20_DEDD$auc))]
    DEDD.results.20_DEDD$pauc <- DEDD.results.20_DEDD$pauc[, grep("HMM", names(DEDD.results.20_DEDD$pauc))]
    DEDD.results.20_DEDD$fdr <- DEDD.results.20_DEDD$fdr[, grep("HMM", names(DEDD.results.20_DEDD$fdr))]
    DEDD.results.20_DEDD$tpr <- DEDD.results.20_DEDD$tpr[, grep("HMM", names(DEDD.results.20_DEDD$tpr))]
    DEDD.results.20_DEDD$mean.fdr <- DEDD.results.20_DEDD$mean.fdr[, grep("HMM", names(DEDD.results.20_DEDD$mean.fdr))]
    DEDD.results.20_DEDD$mean.discoveries <- DEDD.results.20_DEDD$mean.discoveries[
      , grep("HMM", names(DEDD.results.20_DEDD$mean.discoveries))
      ]
    DEDD.results.50_DEDD$auc <- DEDD.results.50_DEDD$auc[, grep("HMM", names(DEDD.results.50_DEDD$auc))]
    DEDD.results.50_DEDD$pauc <- DEDD.results.50_DEDD$pauc[, grep("HMM", names(DEDD.results.50_DEDD$pauc))]
    DEDD.results.50_DEDD$fdr <- DEDD.results.50_DEDD$fdr[, grep("HMM", names(DEDD.results.50_DEDD$fdr))]
    DEDD.results.50_DEDD$tpr <- DEDD.results.50_DEDD$tpr[, grep("HMM", names(DEDD.results.50_DEDD$tpr))]
    DEDD.results.50_DEDD$mean.fdr <- DEDD.results.50_DEDD$mean.fdr[, grep("HMM", names(DEDD.results.50_DEDD$mean.fdr))]
    DEDD.results.50_DEDD$mean.discoveries <- DEDD.results.50_DEDD$mean.discoveries[
      , grep("HMM", names(DEDD.results.50_DEDD$mean.discoveries))
      ]
  }
}


############################################
#### Differences in mean and dispersion ####
############################################

################################
### Differential dispersion ####
################################

###########
## AUC ####
par(mfrow=c(2,5), mar=c(1,2,2,1), mgp=c(3,0.7,0))
boxplot(DD.results.2_DD$auc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.5,0.9), main=paste0("AUC diff disp, DD2"), cex.main=1.5)
abline(v=4.5, col="lightgrey")
legend("topleft", fill=col_vector[1:4], bty='n', cex=1.5, ncol=2, legend=legend)
boxplot(DD.results.5_DD$auc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.45,0.95), main=paste0("AUC diff disp, DD5"), cex.main=1.5)
abline(v=4.5, col="lightgrey")
boxplot(DD.results.10_DD$auc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.45,0.95), main=paste0("AUC diff disp, DD10"), cex.main=1.5)
abline(v=4.5, col="lightgrey")
boxplot(DD.results.20_DD$auc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.45,0.95), main=paste0("AUC diff disp, DD20"), cex.main=1.5)
abline(v=4.5, col="lightgrey")
boxplot(DD.results.50_DD$auc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.45,0.95), main=paste0("AUC diff disp, DD50"), cex.main=1.5)
abline(v=4.5, col="lightgrey")
boxplot(DD.results.2_DEDD$auc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.45,0.95), main=paste0("AUC diff disp, DEDD2"), cex.main=1.5)
abline(v=4.5, col="lightgrey")
boxplot(DD.results.5_DEDD$auc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.45,0.95), main=paste0("AUC diff disp, DEDD5"), cex.main=1.5)
abline(v=4.5, col="lightgrey")
boxplot(DD.results.10_DEDD$auc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.45,0.95), main=paste0("AUC diff disp, DEDD10"), cex.main=1.5)
abline(v=4.5, col="lightgrey")
boxplot(DD.results.20_DEDD$auc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.45,0.95), main=paste0("AUC diff disp, DEDD50"), cex.main=1.5)
abline(v=4.5, col="lightgrey")
boxplot(DD.results.50_DEDD$auc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.45,0.95), main=paste0("AUC diff disp, DEDD50"), cex.main=1.5)
abline(v=4.5, col="lightgrey")
rbind(colMeans(DD.results.2_DD$auc), colMeans(DD.results.5_DD$auc), colMeans(DD.results.10_DD$auc), 
      colMeans(DD.results.20_DD$auc), colMeans(DD.results.50_DD$auc))
rbind(colMeans(DD.results.2_DEDD$auc), colMeans(DD.results.5_DEDD$auc), colMeans(DD.results.10_DEDD$auc), 
      colMeans(DD.results.20_DEDD$auc), colMeans(DD.results.50_DEDD$auc))
# Generally very little difference between methods, but nearly always very slightly higher AUCs for log than untransformed and for 
# expHM than lnHM.

############
## pAUC ####
par(mfrow=c(2,5), mar=c(1,2,2,1), mgp=c(3,0.7,0))
boxplot(DD.results.2_DD$pauc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.035), main=paste0("Partial AUC diff disp, DD2"), cex.main=1.5)
abline(v=4.5, col="lightgrey")
legend("topleft", fill=col_vector[1:4], bty='n', cex=1.5, ncol=2, legend=legend)
boxplot(DD.results.5_DD$pauc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.035), main=paste0("Partial AUC diff disp, DD5"), cex.main=1.5)
abline(v=4.5, col="lightgrey")
boxplot(DD.results.10_DD$pauc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.035), main=paste0("Partial AUC diff disp, DD10"), cex.main=1.5)
abline(v=4.5, col="lightgrey")
boxplot(DD.results.20_DD$pauc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.035), main=paste0("Partial AUC diff disp, DD20"), cex.main=1.5)
abline(v=4.5, col="lightgrey")
boxplot(DD.results.50_DD$pauc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.035), main=paste0("Partial AUC diff disp, DD50"), cex.main=1.5)
abline(v=4.5, col="lightgrey")
boxplot(DD.results.2_DEDD$pauc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.035), main=paste0("Partial AUC diff disp, DEDD2"), cex.main=1.5)
abline(v=4.5, col="lightgrey")
boxplot(DD.results.5_DEDD$pauc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.035), main=paste0("Partial AUC diff disp, DEDD5"), cex.main=1.5)
abline(v=4.5, col="lightgrey")
boxplot(DD.results.10_DEDD$pauc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.035), main=paste0("Partial AUC diff disp, DEDD10"), cex.main=1.5)
abline(v=4.5, col="lightgrey")
boxplot(DD.results.20_DEDD$pauc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.035), main=paste0("Partial AUC diff disp, DEDD20"), cex.main=1.5)
abline(v=4.5, col="lightgrey")
boxplot(DD.results.50_DEDD$pauc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.035), main=paste0("Partial AUC diff disp, DEDD50"), cex.main=1.5)
abline(v=4.5, col="lightgrey")
rbind(colMeans(DD.results.2_DD$pauc), colMeans(DD.results.5_DD$pauc), colMeans(DD.results.10_DD$pauc), 
      colMeans(DD.results.20_DD$pauc), colMeans(DD.results.50_DD$pauc))
rbind(colMeans(DD.results.2_DEDD$pauc), colMeans(DD.results.5_DEDD$pauc), colMeans(DD.results.10_DEDD$pauc), 
      colMeans(DD.results.20_DEDD$pauc), colMeans(DD.results.50_DEDD$pauc))
# Consistently slightly higher pAUCs for log than untransformed. No consistent pattern for expHM v lnHM - lnHM better for blood 
# data, expHM better for muscle.

###########
## FDR ####
par(mfrow=c(2,5), mar=c(1,2,2,1), mgp=c(3,0.7,0))
boxplot(DD.results.2_DD$fdr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.4), main=paste0("FDR diff disp, DD2"), cex.main=1.5)
abline(h=0.05, v=4.5, col='lightgrey')
legend("topleft", fill=col_vector[1:4], bty='n', cex=1.5, ncol=2, legend=legend)
boxplot(DD.results.5_DD$fdr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.4), main=paste0("FDR diff disp, DD5"), cex.main=1.5)
abline(h=0.05, v=4.5, col='lightgrey')
boxplot(DD.results.10_DD$fdr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.4), main=paste0("FDR diff disp, DD10"), cex.main=1.5)
abline(h=0.05, v=4.5, col='lightgrey')
boxplot(DD.results.20_DD$fdr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.4), main=paste0("FDR diff disp, DD20"), cex.main=1.5)
abline(h=0.05, v=4.5, col='lightgrey')
boxplot(DD.results.50_DD$fdr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.4), main=paste0("FDR diff disp, DD50"), cex.main=1.5)
abline(h=0.05, v=4.5, col='lightgrey')
boxplot(DD.results.2_DEDD$fdr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.4), main=paste0("FDR diff disp, DEDD2"), cex.main=1.5)
abline(h=0.05, v=4.5, col='lightgrey')
boxplot(DD.results.5_DEDD$fdr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.4), main=paste0("FDR diff disp, DEDD5"), cex.main=1.5)
abline(h=0.05, v=4.5, col='lightgrey')
boxplot(DD.results.10_DEDD$fdr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.4), main=paste0("FDR diff disp, DEDD10"), cex.main=1.5)
abline(h=0.05, v=4.5, col='lightgrey')
boxplot(DD.results.20_DEDD$fdr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.4), main=paste0("FDR diff disp, DEDD20"), cex.main=1.5)
abline(h=0.05, v=4.5, col='lightgrey')
boxplot(DD.results.50_DEDD$fdr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.4), main=paste0("FDR diff disp, DEDD50"), cex.main=1.5)
abline(h=0.05, v=4.5, col='lightgrey')
# Mean FDRs:
rbind(colMeans(DD.results.2_DD$fdr, na.rm=T), colMeans(DD.results.5_DD$fdr, na.rm=T), colMeans(DD.results.10_DD$fdr, na.rm=T), 
      colMeans(DD.results.20_DD$fdr, na.rm=T), colMeans(DD.results.50_DD$fdr, na.rm=T))
rbind(colMeans(DD.results.2_DEDD$fdr, na.rm=T), colMeans(DD.results.5_DEDD$fdr, na.rm=T), colMeans(DD.results.10_DEDD$fdr, na.rm=T), 
      colMeans(DD.results.20_DEDD$fdr, na.rm=T), colMeans(DD.results.50_DEDD$fdr, na.rm=T))
# Mean squared distances from 0.05:
rbind(colMeans((DD.results.2_DD$fdr - 0.05)^2, na.rm=T), colMeans((DD.results.5_DD$fdr - 0.05)^2, na.rm=T), 
      colMeans((DD.results.10_DD$fdr - 0.05)^2, na.rm=T), colMeans((DD.results.20_DD$fdr - 0.05)^2, na.rm=T), 
      colMeans((DD.results.50_DD$fdr - 0.05)^2, na.rm=T))
rbind(colMeans((DD.results.2_DEDD$fdr - 0.05)^2, na.rm=T), colMeans((DD.results.5_DEDD$fdr - 0.05)^2, na.rm=T), 
      colMeans((DD.results.10_DEDD$fdr - 0.05)^2, na.rm=T), colMeans((DD.results.20_DEDD$fdr - 0.05)^2, na.rm=T), 
      colMeans((DD.results.50_DEDD$fdr - 0.05)^2, na.rm=T))
# Generally untransformed closer to 0.05 than log, and lnHM closer to 0.05 than expHM, but all really bad.

###########
## TPR ####
par(mfrow=c(2,5), mar=c(1,2,2,1), mgp=c(3,0.7,0))
boxplot(DD.results.2_DD$tpr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.65), main=paste0("TPR diff disp, DD2"), cex.main=1.5)
abline(v=4.5, col="lightgrey")
legend("topleft", fill=col_vector[1:4], bty='n', cex=1.5, ncol=2, legend=legend)
boxplot(DD.results.5_DD$tpr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.65), main=paste0("TPR diff disp, DD5"), cex.main=1.5)
abline(v=4.5, col="lightgrey")
boxplot(DD.results.10_DD$tpr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.65), main=paste0("TPR diff disp, DD10"), cex.main=1.5)
abline(v=4.5, col="lightgrey")
boxplot(DD.results.20_DD$tpr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.65), main=paste0("TPR diff disp, DD20"), cex.main=1.5)
abline(v=4.5, col="lightgrey")
boxplot(DD.results.50_DD$tpr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.65), main=paste0("TPR diff disp, DD50"), cex.main=1.5)
abline(v=4.5, col="lightgrey")
boxplot(DD.results.2_DEDD$tpr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.65), main=paste0("TPR diff disp, DEDD2"), cex.main=1.5)
abline(v=4.5, col="lightgrey")
boxplot(DD.results.5_DEDD$tpr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.65), main=paste0("TPR diff disp, DEDD5"), cex.main=1.5)
abline(v=4.5, col="lightgrey")
boxplot(DD.results.10_DEDD$tpr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.65), main=paste0("TPR diff disp, DEDD10"), cex.main=1.5)
abline(v=4.5, col="lightgrey")
boxplot(DD.results.20_DEDD$tpr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.65), main=paste0("TPR diff disp, DEDD20"), cex.main=1.5)
abline(v=4.5, col="lightgrey")
boxplot(DD.results.50_DEDD$tpr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.65), main=paste0("TPR diff disp, DEDD50"), cex.main=1.5)
abline(v=4.5, col="lightgrey")
rbind(colMeans(DD.results.2_DD$tpr), colMeans(DD.results.5_DD$tpr), colMeans(DD.results.10_DD$tpr), 
      colMeans(DD.results.20_DD$tpr), colMeans(DD.results.50_DD$tpr))
rbind(colMeans(DD.results.2_DEDD$tpr), colMeans(DD.results.5_DEDD$tpr), colMeans(DD.results.10_DEDD$tpr), 
      colMeans(DD.results.20_DEDD$tpr), colMeans(DD.results.50_DEDD$tpr))
# TPR consistently higher for log than untransformed. lnHM generally better than expHM for small samples, and expHM generally 
# slightly better than lnHM for large samples (20 or 50 per group).

#############################
## False discovery plots ####
par(mfrow=c(2,5), mar=c(2,2,2,1), mgp=c(3,0.7,0))

plot(DD.results.2_DD$mean.discoveries$expHM_blood, DD.results.2_DD$mean.fdr$expHM_blood, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], main=paste0("DD2"), cex.main=1.5)
lines(DD.results.2_DD$mean.discoveries$expHM.log_blood, DD.results.2_DD$mean.fdr$expHM.log_blood, 
      col=col_vector[2])
lines(DD.results.2_DD$mean.discoveries$lnHM_blood, DD.results.2_DD$mean.fdr$lnHM_blood, 
      col=col_vector[3])
lines(DD.results.2_DD$mean.discoveries$lnHM.log_blood, DD.results.2_DD$mean.fdr$lnHM.log_blood, 
      col=col_vector[4])
lines(DD.results.2_DD$mean.discoveries$expHM_muscle, DD.results.2_DD$mean.fdr$expHM_muscle, 
      lty=2, col=col_vector[1])
lines(DD.results.2_DD$mean.discoveries$expH.logM_muscle, DD.results.2_DD$mean.fdr$expHM.log_muscle, 
      lty=2, col=col_vector[2])
lines(DD.results.2_DD$mean.discoveries$lnHM_muscle, DD.results.2_DD$mean.fdr$lnHM_muscle, 
      lty=2, col=col_vector[3])
lines(DD.results.2_DD$mean.discoveries$lnHM.log_muscle, DD.results.2_DD$mean.fdr$lnHM.log_muscle, 
      lty=2, col=col_vector[4])
abline(h=0.05, col='lightgrey')
plot(DD.results.5_DD$mean.discoveries$expHM_blood, DD.results.5_DD$mean.fdr$expHM_blood, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], main=paste0("DD5"), cex.main=1.5)
lines(DD.results.5_DD$mean.discoveries$expHM.log_blood, DD.results.5_DD$mean.fdr$expHM.log_blood, 
      col=col_vector[2])
lines(DD.results.5_DD$mean.discoveries$lnHM_blood, DD.results.5_DD$mean.fdr$lnHM_blood, 
      col=col_vector[3])
lines(DD.results.5_DD$mean.discoveries$lnHM.log_blood, DD.results.5_DD$mean.fdr$lnHM.log_blood, 
      col=col_vector[4])
lines(DD.results.5_DD$mean.discoveries$expHM_muscle, DD.results.5_DD$mean.fdr$expHM_muscle, 
      lty=2, col=col_vector[1])
lines(DD.results.5_DD$mean.discoveries$expHM.log_muscle, DD.results.5_DD$mean.fdr$expHM.log_muscle, 
      lty=2, col=col_vector[2])
lines(DD.results.5_DD$mean.discoveries$lnHM_muscle, DD.results.5_DD$mean.fdr$lnHM_muscle, 
      lty=2, col=col_vector[3])
lines(DD.results.5_DD$mean.discoveries$lnHM.log_muscle, DD.results.5_DD$mean.fdr$lnHM.log_muscle, 
      lty=2, col=col_vector[4])
abline(h=0.05, col='lightgrey')
plot(DD.results.10_DD$mean.discoveries$expHM_blood, DD.results.10_DD$mean.fdr$expHM_blood, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], main=paste0("DD10"), cex.main=1.5)
lines(DD.results.10_DD$mean.discoveries$expHM.log_blood, DD.results.10_DD$mean.fdr$expHM.log_blood, 
      col=col_vector[2])
lines(DD.results.10_DD$mean.discoveries$lnHM_blood, DD.results.10_DD$mean.fdr$lnHM_blood, 
      col=col_vector[3])
lines(DD.results.10_DD$mean.discoveries$lnHM.log_blood, DD.results.10_DD$mean.fdr$lnHM.log_blood, 
      col=col_vector[4])
lines(DD.results.10_DD$mean.discoveries$expHM_muscle, DD.results.10_DD$mean.fdr$expHM_muscle, 
      lty=2, col=col_vector[1])
lines(DD.results.10_DD$mean.discoveries$expHM.log_muscle, DD.results.10_DD$mean.fdr$expHM.log_muscle, 
      lty=2, col=col_vector[2])
lines(DD.results.10_DD$mean.discoveries$lnHM_muscle, DD.results.10_DD$mean.fdr$lnHM_muscle, 
      lty=2, col=col_vector[3])
lines(DD.results.10_DD$mean.discoveries$lnHM.log_muscle, DD.results.10_DD$mean.fdr$lnHM.log_muscle, 
      lty=2, col=col_vector[4])
abline(h=0.05, col='lightgrey')
plot(DD.results.20_DD$mean.discoveries$expHM_blood, DD.results.20_DD$mean.fdr$expHM_blood, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], main=paste0("DD20"), cex.main=1.5)
lines(DD.results.20_DD$mean.discoveries$expHM.log_blood, DD.results.20_DD$mean.fdr$expHM.log_blood, 
      col=col_vector[2])
lines(DD.results.20_DD$mean.discoveries$lnHM_blood, DD.results.20_DD$mean.fdr$lnHM_blood, 
      col=col_vector[3])
lines(DD.results.20_DD$mean.discoveries$lnHM.log_blood, DD.results.20_DD$mean.fdr$lnHM.log_blood, 
      col=col_vector[4])
lines(DD.results.20_DD$mean.discoveries$expHM_muscle, DD.results.20_DD$mean.fdr$expHM_muscle, 
      lty=2, col=col_vector[1])
lines(DD.results.20_DD$mean.discoveries$expHM.log_muscle, DD.results.20_DD$mean.fdr$expHM.log_muscle, 
      lty=2, col=col_vector[2])
lines(DD.results.20_DD$mean.discoveries$lnHM_muscle, DD.results.20_DD$mean.fdr$lnHM_muscle, 
      lty=2, col=col_vector[3])
lines(DD.results.20_DD$mean.discoveries$lnHM.log_muscle, DD.results.20_DD$mean.fdr$lnHM.log_muscle, 
      lty=2, col=col_vector[4])
abline(h=0.05, col='lightgrey')
plot(DD.results.50_DD$mean.discoveries$expHM_blood, DD.results.50_DD$mean.fdr$expHM_blood, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], main=paste0("DD50"), cex.main=1.5)
lines(DD.results.50_DD$mean.discoveries$expHM.log_blood, DD.results.50_DD$mean.fdr$expHM.log_blood, 
      col=col_vector[2])
lines(DD.results.50_DD$mean.discoveries$lnHM_blood, DD.results.50_DD$mean.fdr$lnHM_blood, 
      col=col_vector[3])
lines(DD.results.50_DD$mean.discoveries$lnHM.log_blood, DD.results.50_DD$mean.fdr$lnHM.log_blood, 
      col=col_vector[4])
lines(DD.results.50_DD$mean.discoveries$expHM_muscle, DD.results.50_DD$mean.fdr$expHM_muscle, 
      lty=2, col=col_vector[1])
lines(DD.results.50_DD$mean.discoveries$expHM.log_muscle, DD.results.50_DD$mean.fdr$expHM.log_muscle, 
      lty=2, col=col_vector[2])
lines(DD.results.50_DD$mean.discoveries$lnHM_muscle, DD.results.50_DD$mean.fdr$lnHM_muscle, 
      lty=2, col=col_vector[3])
lines(DD.results.50_DD$mean.discoveries$lnHM.log_muscle, DD.results.50_DD$mean.fdr$lnHM.log_muscle, 
      lty=2, col=col_vector[4])
abline(h=0.05, col='lightgrey')
legend("topleft", bty='n', col=col_vector[1:4], lty=1, ncol=1, cex=1.5, legend=legend)

plot(DD.results.2_DEDD$mean.discoveries$expHM_blood, DD.results.2_DEDD$mean.fdr$expHM_blood, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], main=paste0("DEDD2"), cex.main=1.5)
lines(DD.results.2_DEDD$mean.discoveries$expHM.log_blood, DD.results.2_DEDD$mean.fdr$expHM.log_blood, 
      col=col_vector[2])
lines(DD.results.2_DEDD$mean.discoveries$lnHM_blood, DD.results.2_DEDD$mean.fdr$lnHM_blood, 
      col=col_vector[3])
lines(DD.results.2_DEDD$mean.discoveries$lnHM.log_blood, DD.results.2_DEDD$mean.fdr$lnHM.log_blood, 
      col=col_vector[4])
lines(DD.results.2_DEDD$mean.discoveries$expHM_muscle, DD.results.2_DEDD$mean.fdr$expHM_muscle, 
      lty=2, col=col_vector[1])
lines(DD.results.2_DEDD$mean.discoveries$expHM.log_muscle, DD.results.2_DEDD$mean.fdr$expHM.log_muscle, 
      lty=2, col=col_vector[2])
lines(DD.results.2_DEDD$mean.discoveries$lnHM_muscle, DD.results.2_DEDD$mean.fdr$lnHM_muscle, 
      lty=2, col=col_vector[3])
lines(DD.results.2_DEDD$mean.discoveries$lnHM.log_muscle, DD.results.2_DEDD$mean.fdr$lnHM.log_muscle, 
      lty=2, col=col_vector[4])
abline(h=0.05, col='lightgrey')
plot(DD.results.5_DEDD$mean.discoveries$expHM_blood, DD.results.5_DEDD$mean.fdr$expHM_blood, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], main=paste0("DEDD5"), cex.main=1.5)
lines(DD.results.5_DEDD$mean.discoveries$expHM.log_blood, DD.results.5_DEDD$mean.fdr$expHM.log_blood, 
      col=col_vector[2])
lines(DD.results.5_DEDD$mean.discoveries$lnHM_blood, DD.results.5_DEDD$mean.fdr$lnHM_blood, 
      col=col_vector[3])
lines(DD.results.5_DEDD$mean.discoveries$lnHM.log_blood, DD.results.5_DEDD$mean.fdr$lnHM.log_blood, 
      col=col_vector[4])
lines(DD.results.5_DEDD$mean.discoveries$expHM_muscle, DD.results.5_DEDD$mean.fdr$expHM_muscle, 
      lty=2, col=col_vector[1])
lines(DD.results.5_DEDD$mean.discoveries$expHM.log_muscle, DD.results.5_DEDD$mean.fdr$expHM.log_muscle, 
      lty=2, col=col_vector[2])
lines(DD.results.5_DEDD$mean.discoveries$lnHM_muscle, DD.results.5_DEDD$mean.fdr$lnHM_muscle, 
      lty=2, col=col_vector[3])
lines(DD.results.5_DEDD$mean.discoveries$lnHM.log_muscle, DD.results.5_DEDD$mean.fdr$lnHM.log_muscle, 
      lty=2, col=col_vector[4])
abline(h=0.05, col='lightgrey')
plot(DD.results.10_DEDD$mean.discoveries$expHM_blood, DD.results.10_DEDD$mean.fdr$expHM_blood, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], main=paste0("DEDD10"), cex.main=1.5)
lines(DD.results.10_DEDD$mean.discoveries$expHM.log_blood, DD.results.10_DEDD$mean.fdr$expHM.log_blood, 
      col=col_vector[2])
lines(DD.results.10_DEDD$mean.discoveries$lnHM_blood, DD.results.10_DEDD$mean.fdr$lnHM_blood, 
      col=col_vector[3])
lines(DD.results.10_DEDD$mean.discoveries$lnHM.log_blood, DD.results.10_DEDD$mean.fdr$lnHM.log_blood, 
      col=col_vector[4])
lines(DD.results.10_DEDD$mean.discoveries$expHM_muscle, DD.results.10_DEDD$mean.fdr$expHM_muscle, 
      lty=2, col=col_vector[1])
lines(DD.results.10_DEDD$mean.discoveries$expHM.log_muscle, DD.results.10_DEDD$mean.fdr$expHM.log_muscle, 
      lty=2, col=col_vector[2])
lines(DD.results.10_DEDD$mean.discoveries$lnHM_muscle, DD.results.10_DEDD$mean.fdr$lnHM_muscle, 
      lty=2, col=col_vector[3])
lines(DD.results.10_DEDD$mean.discoveries$lnHM.log_muscle, DD.results.10_DEDD$mean.fdr$lnHM.log_muscle, 
      lty=2, col=col_vector[4])
abline(h=0.05, col='lightgrey')
plot(DD.results.20_DEDD$mean.discoveries$expHM_blood, DD.results.20_DEDD$mean.fdr$expHM_blood, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], main=paste0("DEDD20"), cex.main=1.5)
lines(DD.results.20_DEDD$mean.discoveries$expHM.log_blood, DD.results.20_DEDD$mean.fdr$expHM.log_blood, 
      col=col_vector[2])
lines(DD.results.20_DEDD$mean.discoveries$lnHM_blood, DD.results.20_DEDD$mean.fdr$lnHM_blood, 
      col=col_vector[3])
lines(DD.results.20_DEDD$mean.discoveries$lnHM.log_blood, DD.results.20_DEDD$mean.fdr$lnHM.log_blood, 
      col=col_vector[4])
lines(DD.results.20_DEDD$mean.discoveries$expHM_muscle, DD.results.20_DEDD$mean.fdr$expHM_muscle, 
      lty=2, col=col_vector[1])
lines(DD.results.20_DEDD$mean.discoveries$expHM.log_muscle, DD.results.20_DEDD$mean.fdr$expHM.log_muscle, 
      lty=2, col=col_vector[2])
lines(DD.results.20_DEDD$mean.discoveries$lnHM_muscle, DD.results.20_DEDD$mean.fdr$lnHM_muscle, 
      lty=2, col=col_vector[3])
lines(DD.results.20_DEDD$mean.discoveries$lnHM.log_muscle, DD.results.20_DEDD$mean.fdr$lnHM.log_muscle, 
      lty=2, col=col_vector[4])
abline(h=0.05, col='lightgrey')
plot(DD.results.50_DEDD$mean.discoveries$expHM_blood, DD.results.50_DEDD$mean.fdr$expHM_blood, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], main=paste0("DEDD50"), cex.main=1.5)
lines(DD.results.50_DEDD$mean.discoveries$expHM.log_blood, DD.results.50_DEDD$mean.fdr$expHM.log_blood, 
      col=col_vector[2])
lines(DD.results.50_DEDD$mean.discoveries$lnHM_blood, DD.results.50_DEDD$mean.fdr$lnHM_blood, 
      col=col_vector[3])
lines(DD.results.50_DEDD$mean.discoveries$lnHM.log_blood, DD.results.50_DEDD$mean.fdr$lnHM.log_blood, 
      col=col_vector[4])
lines(DD.results.50_DEDD$mean.discoveries$expHM_muscle, DD.results.50_DEDD$mean.fdr$expHM_muscle, 
      lty=2, col=col_vector[1])
lines(DD.results.50_DEDD$mean.discoveries$expHM.log_muscle, DD.results.50_DEDD$mean.fdr$expHM.log_muscle, 
      lty=2, col=col_vector[2])
lines(DD.results.50_DEDD$mean.discoveries$lnHM_muscle, DD.results.50_DEDD$mean.fdr$lnHM_muscle, 
      lty=2, col=col_vector[3])
lines(DD.results.50_DEDD$mean.discoveries$lnHM.log_muscle, DD.results.50_DEDD$mean.fdr$lnHM.log_muscle, 
      lty=2, col=col_vector[4])
abline(h=0.05, col='lightgrey')
legend("topleft", bty='n', col=col_vector[1:4], lty=1, ncol=1, cex=1.5, legend=legend)
# FDR curves virtually indistinguishable, just very slightly better for lnHM with 5 or 10 samples per group for blood, very 
# slightly better for expHM for some sample sizes for muscle, and generally very very better for log than untransformed.


################################
### Differential expression ####
################################

###########
## AUC ####
par(mfrow=c(2,5), mar=c(1,2,2,1), mgp=c(3,0.7,0))
boxplot(DE.results.2_DE$auc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.65,1), main=paste0("AUC diff exp, DE2"), cex.main=1.5)
abline(v=4.5, col="lightgrey")
legend("topleft", fill=col_vector[1:4], bty='n', cex=1.5, ncol=2, legend=legend)
boxplot(DE.results.5_DE$auc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.65,1), main=paste0("AUC diff exp, DE5"), cex.main=1.5)
abline(v=4.5, col="lightgrey")
boxplot(DE.results.10_DE$auc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.65,1), main=paste0("AUC diff exp, DE10"), cex.main=1.5)
abline(v=4.5, col="lightgrey")
boxplot(DE.results.20_DE$auc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.65,1), main=paste0("AUC diff exp, DE20"), cex.main=1.5)
abline(v=4.5, col="lightgrey")
boxplot(DE.results.50_DE$auc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.65,1), main=paste0("AUC diff exp, DE50"), cex.main=1.5)
abline(v=4.5, col="lightgrey")
boxplot(DE.results.2_DEDD$auc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.65,1), main=paste0("AUC diff exp, DEDD2"), cex.main=1.5)
abline(v=4.5, col="lightgrey")
boxplot(DE.results.5_DEDD$auc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.65,1), main=paste0("AUC diff exp, DEDD5"), cex.main=1.5)
abline(v=4.5, col="lightgrey")
boxplot(DE.results.10_DEDD$auc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.65,1), main=paste0("AUC diff exp, DEDD10"), cex.main=1.5)
abline(v=4.5, col="lightgrey")
boxplot(DE.results.20_DEDD$auc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.65,1), main=paste0("AUC diff exp, DEDD50"), cex.main=1.5)
abline(v=4.5, col="lightgrey")
boxplot(DE.results.50_DEDD$auc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.65,1), main=paste0("AUC diff exp, DEDD50"), cex.main=1.5)
abline(v=4.5, col="lightgrey")
rbind(colMeans(DE.results.2_DE$auc), colMeans(DE.results.5_DE$auc), colMeans(DE.results.10_DE$auc), 
      colMeans(DE.results.20_DE$auc), colMeans(DE.results.50_DE$auc))
rbind(colMeans(DE.results.2_DEDD$auc), colMeans(DE.results.5_DEDD$auc), colMeans(DE.results.10_DEDD$auc), 
      colMeans(DE.results.20_DEDD$auc), colMeans(DE.results.50_DEDD$auc))
# AUCs higher for untranslated than log for 2 and 5 samples per group, higher for log for 20 and 50. Consistently higher for lnHM 
# than expHM.

############
## pAUC ####
par(mfrow=c(2,5), mar=c(1,2,2,1), mgp=c(3,0.7,0))
boxplot(DE.results.2_DE$pauc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.005,0.047), main=paste0("Partial AUC diff exp, DE2"), cex.main=1.5)
abline(v=4.5, col="lightgrey")
legend("topleft", fill=col_vector[1:4], bty='n', cex=1.5, ncol=2, legend=legend)
boxplot(DE.results.5_DE$pauc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.005,0.047), main=paste0("Partial AUC diff exp, DE5"), cex.main=1.5)
abline(v=4.5, col="lightgrey")
boxplot(DE.results.10_DE$pauc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.005,0.047), main=paste0("Partial AUC diff exp, DE10"), cex.main=1.5)
abline(v=4.5, col="lightgrey")
boxplot(DE.results.20_DE$pauc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.005,0.047), main=paste0("Partial AUC diff exp, DE20"), cex.main=1.5)
abline(v=4.5, col="lightgrey")
boxplot(DE.results.50_DE$pauc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.005,0.047), main=paste0("Partial AUC diff exp, DE50"), cex.main=1.5)
abline(v=4.5, col="lightgrey")
boxplot(DE.results.2_DEDD$pauc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.005,0.047), main=paste0("Partial AUC diff exp, DEDD2"), cex.main=1.5)
abline(v=4.5, col="lightgrey")
boxplot(DE.results.5_DEDD$pauc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.005,0.047), main=paste0("Partial AUC diff exp, DEDD5"), cex.main=1.5)
abline(v=4.5, col="lightgrey")
boxplot(DE.results.10_DEDD$pauc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.005,0.047), main=paste0("Partial AUC diff exp, DEDD10"), cex.main=1.5)
abline(v=4.5, col="lightgrey")
boxplot(DE.results.20_DEDD$pauc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.005,0.047), main=paste0("Partial AUC diff exp, DEDD20"), cex.main=1.5)
abline(v=4.5, col="lightgrey")
boxplot(DE.results.50_DEDD$pauc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.005,0.047), main=paste0("Partial AUC diff exp, DEDD50"), cex.main=1.5)
abline(v=4.5, col="lightgrey")
rbind(colMeans(DE.results.2_DE$pauc), colMeans(DE.results.5_DE$pauc), colMeans(DE.results.10_DE$pauc), 
      colMeans(DE.results.20_DE$pauc), colMeans(DE.results.50_DE$pauc))
rbind(colMeans(DE.results.2_DEDD$pauc), colMeans(DE.results.5_DEDD$pauc), colMeans(DE.results.10_DEDD$pauc), 
      colMeans(DE.results.20_DEDD$pauc), colMeans(DE.results.50_DEDD$pauc))
# Partial AUCs higher for untranslated than log, except for 50 samples per group, and consistently higher for lnHM than expHM.

###########
## FDR ####
par(mfrow=c(2,5), mar=c(1,2,2,1), mgp=c(3,0.7,0))
boxplot(DE.results.2_DE$fdr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.4), main=paste0("FDR diff exp, DE2"), cex.main=1.5)
abline(h=0.05, v=4.5, col='lightgrey')
legend("topleft", fill=col_vector[1:4], bty='n', cex=1.5, ncol=2, legend=legend)
boxplot(DE.results.5_DE$fdr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.4), main=paste0("FDR diff exp, DE5"), cex.main=1.5)
abline(h=0.05, v=4.5, col='lightgrey')
boxplot(DE.results.10_DE$fdr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.4), main=paste0("FDR diff exp, DE10"), cex.main=1.5)
abline(h=0.05, v=4.5, col='lightgrey')
boxplot(DE.results.20_DE$fdr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.4), main=paste0("FDR diff exp, DE20"), cex.main=1.5)
abline(h=0.05, v=4.5, col='lightgrey')
boxplot(DE.results.50_DE$fdr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.4), main=paste0("FDR diff exp, DE50"), cex.main=1.5)
abline(h=0.05, v=4.5, col='lightgrey')
boxplot(DE.results.2_DEDD$fdr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.4), main=paste0("FDR diff exp, DEDD2"), cex.main=1.5)
abline(h=0.05, v=4.5, col='lightgrey')
boxplot(DE.results.5_DEDD$fdr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.4), main=paste0("FDR diff exp, DEDD5"), cex.main=1.5)
abline(h=0.05, v=4.5, col='lightgrey')
boxplot(DE.results.10_DEDD$fdr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.4), main=paste0("FDR diff exp, DEDD10"), cex.main=1.5)
abline(h=0.05, v=4.5, col='lightgrey')
boxplot(DE.results.20_DEDD$fdr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.4), main=paste0("FDR diff exp, DEDD20"), cex.main=1.5)
abline(h=0.05, v=4.5, col='lightgrey')
boxplot(DE.results.50_DEDD$fdr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.4), main=paste0("FDR diff exp, DEDD50"), cex.main=1.5)
abline(h=0.05, v=4.5, col='lightgrey')
# Mean FDRs:
rbind(colMeans(DE.results.2_DE$fdr, na.rm=T), colMeans(DE.results.5_DE$fdr, na.rm=T), colMeans(DE.results.10_DE$fdr, na.rm=T), 
      colMeans(DE.results.20_DE$fdr, na.rm=T), colMeans(DE.results.50_DE$fdr, na.rm=T))
rbind(colMeans(DE.results.2_DEDD$fdr, na.rm=T), colMeans(DE.results.5_DEDD$fdr, na.rm=T), colMeans(DE.results.10_DEDD$fdr, na.rm=T), 
      colMeans(DE.results.20_DEDD$fdr, na.rm=T), colMeans(DE.results.50_DEDD$fdr, na.rm=T))
# Mean squared distances from 0.05:
rbind(colMeans((DE.results.2_DE$fdr - 0.05)^2, na.rm=T), colMeans((DE.results.5_DE$fdr - 0.05)^2, na.rm=T), 
      colMeans((DE.results.10_DE$fdr - 0.05)^2, na.rm=T), colMeans((DE.results.20_DE$fdr - 0.05)^2, na.rm=T), 
      colMeans((DE.results.50_DE$fdr - 0.05)^2, na.rm=T))
rbind(colMeans((DE.results.2_DEDD$fdr - 0.05)^2, na.rm=T), colMeans((DE.results.5_DEDD$fdr - 0.05)^2, na.rm=T), 
      colMeans((DE.results.10_DEDD$fdr - 0.05)^2, na.rm=T), colMeans((DE.results.20_DEDD$fdr - 0.05)^2, na.rm=T), 
      colMeans((DE.results.50_DEDD$fdr - 0.05)^2, na.rm=T))
# FDRs always closer to 0.05 for untranslated than log, and generally closer to 0.05 for lnHM than expHM, but always really bad.

###########
## TPR ####
par(mfrow=c(2,5), mar=c(1,2,2,1), mgp=c(3,0.7,0))
boxplot(DE.results.2_DE$tpr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.96), main=paste0("TPR diff exp, DE2"), cex.main=1.5)
legend("topleft", fill=col_vector[1:4], bty='n', cex=1.5, ncol=2, legend=legend)
abline(v=4.5, col="lightgrey")
boxplot(DE.results.5_DE$tpr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.96), main=paste0("TPR diff exp, DE5"), cex.main=1.5)
abline(v=4.5, col="lightgrey")
boxplot(DE.results.10_DE$tpr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.96), main=paste0("TPR diff exp, DE10"), cex.main=1.5)
abline(v=4.5, col="lightgrey")
boxplot(DE.results.20_DE$tpr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.96), main=paste0("TPR diff exp, DE20"), cex.main=1.5)
abline(v=4.5, col="lightgrey")
boxplot(DE.results.50_DE$tpr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.96), main=paste0("TPR diff exp, DE50"), cex.main=1.5)
abline(v=4.5, col="lightgrey")
boxplot(DE.results.2_DEDD$tpr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.96), main=paste0("TPR diff exp, DEDD2"), cex.main=1.5)
abline(v=4.5, col="lightgrey")
boxplot(DE.results.5_DEDD$tpr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.96), main=paste0("TPR diff exp, DEDD5"), cex.main=1.5)
abline(v=4.5, col="lightgrey")
boxplot(DE.results.10_DEDD$tpr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.96), main=paste0("TPR diff exp, DEDD10"), cex.main=1.5)
abline(v=4.5, col="lightgrey")
boxplot(DE.results.20_DEDD$tpr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.96), main=paste0("TPR diff exp, DEDD20"), cex.main=1.5)
abline(v=4.5, col="lightgrey")
boxplot(DE.results.50_DEDD$tpr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.96), main=paste0("TPR diff exp, DEDD50"), cex.main=1.5)
abline(v=4.5, col="lightgrey")
rbind(colMeans(DE.results.2_DE$tpr), colMeans(DE.results.5_DE$tpr), colMeans(DE.results.10_DE$tpr), 
      colMeans(DE.results.20_DE$tpr), colMeans(DE.results.50_DE$tpr))
rbind(colMeans(DE.results.2_DEDD$tpr), colMeans(DE.results.5_DEDD$tpr), colMeans(DE.results.10_DEDD$tpr), 
      colMeans(DE.results.20_DEDD$tpr), colMeans(DE.results.50_DEDD$tpr))
# TPRs always higher for log than untransformed, and generally higher for lnHM than expHM.

#############################
## False discovery plots ####
par(mfrow=c(2,5), mar=c(2,2,2,1), mgp=c(3,0.7,0))

plot(DE.results.2_DE$mean.discoveries$expHM_blood, DE.results.2_DE$mean.fdr$expHM_blood, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], main=paste0("DE2"), cex.main=1.5)
lines(DE.results.2_DE$mean.discoveries$expHM.log_blood, DE.results.2_DE$mean.fdr$expHM.log_blood, 
      col=col_vector[2])
lines(DE.results.2_DE$mean.discoveries$lnHM_blood, DE.results.2_DE$mean.fdr$lnHM_blood, 
      col=col_vector[3])
lines(DE.results.2_DE$mean.discoveries$lnHM.log_blood, DE.results.2_DE$mean.fdr$lnHM.log_blood, 
      col=col_vector[4])
lines(DE.results.2_DE$mean.discoveries$expHM_muscle, DE.results.2_DE$mean.fdr$expHM_muscle, 
      lty=2, col=col_vector[1])
lines(DE.results.2_DE$mean.discoveries$expHM.log_muscle, DE.results.2_DE$mean.fdr$expHM.log_muscle, 
      lty=2, col=col_vector[2])
lines(DE.results.2_DE$mean.discoveries$lnHM_muscle, DE.results.2_DE$mean.fdr$lnHM_muscle, 
      lty=2, col=col_vector[3])
lines(DE.results.2_DE$mean.discoveries$lnHM.log_muscle, DE.results.2_DE$mean.fdr$lnHM.log_muscle, 
      lty=2, col=col_vector[4])
abline(h=0.05, col='lightgrey')
plot(DE.results.5_DE$mean.discoveries$expHM_blood, DE.results.5_DE$mean.fdr$expHM_blood, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], main=paste0("DE5"), cex.main=1.5)
lines(DE.results.5_DE$mean.discoveries$expHM.log_blood, DE.results.5_DE$mean.fdr$expHM.log_blood, 
      col=col_vector[2])
lines(DE.results.5_DE$mean.discoveries$lnHM_blood, DE.results.5_DE$mean.fdr$lnHM_blood, 
      col=col_vector[3])
lines(DE.results.5_DE$mean.discoveries$lnHM.log_blood, DE.results.5_DE$mean.fdr$lnHM.log_blood, 
      col=col_vector[4])
lines(DE.results.5_DE$mean.discoveries$expHM_muscle, DE.results.5_DE$mean.fdr$expHM_muscle, 
      lty=2, col=col_vector[1])
lines(DE.results.5_DE$mean.discoveries$expHM.log_muscle, DE.results.5_DE$mean.fdr$expHM.log_muscle, 
      lty=2, col=col_vector[2])
lines(DE.results.5_DE$mean.discoveries$lnHM_muscle, DE.results.5_DE$mean.fdr$lnHM_muscle, 
      lty=2, col=col_vector[3])
lines(DE.results.5_DE$mean.discoveries$lnHM.log_muscle, DE.results.5_DE$mean.fdr$lnHM.log_muscle, 
      lty=2, col=col_vector[4])
abline(h=0.05, col='lightgrey')
plot(DE.results.10_DE$mean.discoveries$expHM_blood, DE.results.10_DE$mean.fdr$expHM_blood, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], main=paste0("DE10"), cex.main=1.5)
lines(DE.results.10_DE$mean.discoveries$expHM.log_blood, DE.results.10_DE$mean.fdr$expHM.log_blood, 
      col=col_vector[2])
lines(DE.results.10_DE$mean.discoveries$lnHM_blood, DE.results.10_DE$mean.fdr$lnHM_blood, 
      col=col_vector[3])
lines(DE.results.10_DE$mean.discoveries$lnHM.log_blood, DE.results.10_DE$mean.fdr$lnHM.log_blood, 
      col=col_vector[4])
lines(DE.results.10_DE$mean.discoveries$expHM_muscle, DE.results.10_DE$mean.fdr$expHM_muscle, 
      lty=2, col=col_vector[1])
lines(DE.results.10_DE$mean.discoveries$expHM.log_muscle, DE.results.10_DE$mean.fdr$expHM.log_muscle, 
      lty=2, col=col_vector[2])
lines(DE.results.10_DE$mean.discoveries$lnHM_muscle, DE.results.10_DE$mean.fdr$lnHM_muscle, 
      lty=2, col=col_vector[3])
lines(DE.results.10_DE$mean.discoveries$lnHM.log_muscle, DE.results.10_DE$mean.fdr$lnHM.log_muscle, 
      lty=2, col=col_vector[4])
abline(h=0.05, col='lightgrey')
plot(DE.results.20_DE$mean.discoveries$expHM_blood, DE.results.20_DE$mean.fdr$expHM_blood, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], main=paste0("DE20"), cex.main=1.5)
lines(DE.results.20_DE$mean.discoveries$expHM.log_blood, DE.results.20_DE$mean.fdr$expHM.log_blood, 
      col=col_vector[2])
lines(DE.results.20_DE$mean.discoveries$lnHM_blood, DE.results.20_DE$mean.fdr$lnHM_blood, 
      col=col_vector[3])
lines(DE.results.20_DE$mean.discoveries$lnHM.log_blood, DE.results.20_DE$mean.fdr$lnHM.log_blood, 
      col=col_vector[4])
lines(DE.results.20_DE$mean.discoveries$expHM_muscle, DE.results.20_DE$mean.fdr$expHM_muscle, 
      lty=2, col=col_vector[1])
lines(DE.results.20_DE$mean.discoveries$expHM.log_muscle, DE.results.20_DE$mean.fdr$expHM.log_muscle, 
      lty=2, col=col_vector[2])
lines(DE.results.20_DE$mean.discoveries$lnHM_muscle, DE.results.20_DE$mean.fdr$lnHM_muscle, 
      lty=2, col=col_vector[3])
lines(DE.results.20_DE$mean.discoveries$lnHM.log_muscle, DE.results.20_DE$mean.fdr$lnHM.log_muscle, 
      lty=2, col=col_vector[4])
abline(h=0.05, col='lightgrey')
plot(DE.results.50_DE$mean.discoveries$expHM_blood, DE.results.50_DE$mean.fdr$expHM_blood, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], main=paste0("DE50"), cex.main=1.5)
lines(DE.results.50_DE$mean.discoveries$expHM.log_blood, DE.results.50_DE$mean.fdr$expHM.log_blood, 
      col=col_vector[2])
lines(DE.results.50_DE$mean.discoveries$lnHM_blood, DE.results.50_DE$mean.fdr$lnHM_blood, 
      col=col_vector[3])
lines(DE.results.50_DE$mean.discoveries$lnHM.log_blood, DE.results.50_DE$mean.fdr$lnHM.log_blood, 
      col=col_vector[4])
lines(DE.results.50_DE$mean.discoveries$expHM_muscle, DE.results.50_DE$mean.fdr$expHM_muscle, 
      lty=2, col=col_vector[1])
lines(DE.results.50_DE$mean.discoveries$expHM.log_muscle, DE.results.50_DE$mean.fdr$expHM.log_muscle, 
      lty=2, col=col_vector[2])
lines(DE.results.50_DE$mean.discoveries$lnHM_muscle, DE.results.50_DE$mean.fdr$lnHM_muscle, 
      lty=2, col=col_vector[3])
lines(DE.results.50_DE$mean.discoveries$lnHM.log_muscle, DE.results.50_DE$mean.fdr$lnHM.log_muscle, 
      lty=2, col=col_vector[4])
abline(h=0.05, col='lightgrey')
legend("topleft", bty='n', col=col_vector[1:4], lty=1, ncol=1, cex=1.5, legend=legend)

plot(DE.results.2_DEDD$mean.discoveries$expHM_blood, DE.results.2_DEDD$mean.fdr$expHM_blood, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], main=paste0("DEDD2"), cex.main=1.5)
lines(DE.results.2_DEDD$mean.discoveries$expHM.log_blood, DE.results.2_DEDD$mean.fdr$expHM.log_blood, 
      col=col_vector[2])
lines(DE.results.2_DEDD$mean.discoveries$lnHM_blood, DE.results.2_DEDD$mean.fdr$lnHM_blood, 
      col=col_vector[3])
lines(DE.results.2_DEDD$mean.discoveries$lnHM.log_blood, DE.results.2_DEDD$mean.fdr$lnHM.log_blood, 
      col=col_vector[4])
lines(DE.results.2_DEDD$mean.discoveries$expHM_muscle, DE.results.2_DEDD$mean.fdr$expHM_muscle, 
      lty=2, col=col_vector[1])
lines(DE.results.2_DEDD$mean.discoveries$expHM.log_muscle, DE.results.2_DEDD$mean.fdr$expHM.log_muscle, 
      lty=2, col=col_vector[2])
lines(DE.results.2_DEDD$mean.discoveries$lnHM_muscle, DE.results.2_DEDD$mean.fdr$lnHM_muscle, 
      lty=2, col=col_vector[3])
lines(DE.results.2_DEDD$mean.discoveries$lnHM.log_muscle, DE.results.2_DEDD$mean.fdr$lnHM.log_muscle, 
      lty=2, col=col_vector[4])
abline(h=0.05, col='lightgrey')
plot(DE.results.5_DEDD$mean.discoveries$expHM_blood, DE.results.5_DEDD$mean.fdr$expHM_blood, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], main=paste0("DEDD5"), cex.main=1.5)
lines(DE.results.5_DEDD$mean.discoveries$expHM.log_blood, DE.results.5_DEDD$mean.fdr$expHM.log_blood, 
      col=col_vector[2])
lines(DE.results.5_DEDD$mean.discoveries$lnHM_blood, DE.results.5_DEDD$mean.fdr$lnHM_blood, 
      col=col_vector[3])
lines(DE.results.5_DEDD$mean.discoveries$lnHM.log_blood, DE.results.5_DEDD$mean.fdr$lnHM.log_blood, 
      col=col_vector[4])
lines(DE.results.5_DEDD$mean.discoveries$expHM_muscle, DE.results.5_DEDD$mean.fdr$expHM_muscle, 
      lty=2, col=col_vector[1])
lines(DE.results.5_DEDD$mean.discoveries$expHM.log_muscle, DE.results.5_DEDD$mean.fdr$expHM.log_muscle, 
      lty=2, col=col_vector[2])
lines(DE.results.5_DEDD$mean.discoveries$lnHM_muscle, DE.results.5_DEDD$mean.fdr$lnHM_muscle, 
      lty=2, col=col_vector[3])
lines(DE.results.5_DEDD$mean.discoveries$lnHM.log_muscle, DE.results.5_DEDD$mean.fdr$lnHM.log_muscle, 
      lty=2, col=col_vector[4])
abline(h=0.05, col='lightgrey')
plot(DE.results.10_DEDD$mean.discoveries$expHM_blood, DE.results.10_DEDD$mean.fdr$expHM_blood, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], main=paste0("DEDD10"), cex.main=1.5)
lines(DE.results.10_DEDD$mean.discoveries$expHM.log_blood, DE.results.10_DEDD$mean.fdr$expHM.log_blood, 
      col=col_vector[2])
lines(DE.results.10_DEDD$mean.discoveries$lnHM_blood, DE.results.10_DEDD$mean.fdr$lnHM_blood, 
      col=col_vector[3])
lines(DE.results.10_DEDD$mean.discoveries$lnHM.log_blood, DE.results.10_DEDD$mean.fdr$lnHM.log_blood, 
      col=col_vector[4])
lines(DE.results.10_DEDD$mean.discoveries$expHM_muscle, DE.results.10_DEDD$mean.fdr$expHM_muscle, 
      lty=2, col=col_vector[1])
lines(DE.results.10_DEDD$mean.discoveries$expHM.log_muscle, DE.results.10_DEDD$mean.fdr$expHM.log_muscle, 
      lty=2, col=col_vector[2])
lines(DE.results.10_DEDD$mean.discoveries$lnHM_muscle, DE.results.10_DEDD$mean.fdr$lnHM_muscle, 
      lty=2, col=col_vector[3])
lines(DE.results.10_DEDD$mean.discoveries$lnHM.log_muscle, DE.results.10_DEDD$mean.fdr$lnHM.log_muscle, 
      lty=2, col=col_vector[4])
abline(h=0.05, col='lightgrey')
plot(DE.results.20_DEDD$mean.discoveries$expHM_blood, DE.results.20_DEDD$mean.fdr$expHM_blood, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], main=paste0("DEDD20"), cex.main=1.5)
lines(DE.results.20_DEDD$mean.discoveries$expHM.log_blood, DE.results.20_DEDD$mean.fdr$expHM.log_blood, 
      col=col_vector[2])
lines(DE.results.20_DEDD$mean.discoveries$lnHM_blood, DE.results.20_DEDD$mean.fdr$lnHM_blood, 
      col=col_vector[3])
lines(DE.results.20_DEDD$mean.discoveries$lnHM.log_blood, DE.results.20_DEDD$mean.fdr$lnHM.log_blood, 
      col=col_vector[4])
lines(DE.results.20_DEDD$mean.discoveries$expHM_muscle, DE.results.20_DEDD$mean.fdr$expHM_muscle, 
      lty=2, col=col_vector[1])
lines(DE.results.20_DEDD$mean.discoveries$expHM.log_muscle, DE.results.20_DEDD$mean.fdr$expHM.log_muscle, 
      lty=2, col=col_vector[2])
lines(DE.results.20_DEDD$mean.discoveries$lnHM_muscle, DE.results.20_DEDD$mean.fdr$lnHM_muscle, 
      lty=2, col=col_vector[3])
lines(DE.results.20_DEDD$mean.discoveries$lnHM.log_muscle, DE.results.20_DEDD$mean.fdr$lnHM.log_muscle, 
      lty=2, col=col_vector[4])
abline(h=0.05, col='lightgrey')
plot(DE.results.50_DEDD$mean.discoveries$expHM_blood, DE.results.50_DEDD$mean.fdr$expHM_blood, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], main=paste0("DED50"), cex.main=1.5)
lines(DE.results.50_DEDD$mean.discoveries$expHM.log_blood, DE.results.50_DEDD$mean.fdr$expHM.log_blood, 
      col=col_vector[2])
lines(DE.results.50_DEDD$mean.discoveries$lnHM_blood, DE.results.50_DEDD$mean.fdr$lnHM_blood, 
      col=col_vector[3])
lines(DE.results.50_DEDD$mean.discoveries$lnHM.log_blood, DE.results.50_DEDD$mean.fdr$lnHM.log_blood, 
      col=col_vector[4])
lines(DE.results.50_DEDD$mean.discoveries$expHM_muscle, DE.results.50_DEDD$mean.fdr$expHM_muscle, 
      lty=2, col=col_vector[1])
lines(DE.results.50_DEDD$mean.discoveries$expHM.log_muscle, DE.results.50_DEDD$mean.fdr$expHM.log_muscle, 
      lty=2, col=col_vector[2])
lines(DE.results.50_DEDD$mean.discoveries$lnHM_muscle, DE.results.50_DEDD$mean.fdr$lnHM_muscle, 
      lty=2, col=col_vector[3])
lines(DE.results.50_DEDD$mean.discoveries$lnHM.log_muscle, DE.results.50_DEDD$mean.fdr$lnHM.log_muscle, 
      lty=2, col=col_vector[4])
abline(h=0.05, col='lightgrey')
legend("topleft", bty='n', col=col_vector[1:4], lty=1, ncol=1, cex=1.5, legend=legend)
# FDR curves lower for untransformed than log, and lower for lnHM than expHM. Difference is bigger for untransformed v log than for 
# expHM v lnHM for small samples, and all differences are nearly indistinguishable for 20 and 50 samples per group.


##################################
### Differential distribution ####
##################################
legend <- c("expHMM", "lnHMM")

###########
## AUC ####
par(mfrow=c(3,5), mar=c(1,2,2,1), mgp=c(3,0.7,0))
boxplot(DEDD.results.2_DD$auc, names=NA, col=col_vector[1:2], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.43,0.87), main=paste0("AUC diff exp, DD2"), cex.main=1.5)
legend("topleft", fill=col_vector[1:2], bty='n', cex=1.5, ncol=2, legend=legend)
abline(v=2.5, col="lightgrey")
boxplot(DEDD.results.5_DD$auc, names=NA, col=col_vector[1:2], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.43,0.87), main=paste0("AUC diff exp, DD5"), cex.main=1.5)
abline(v=2.5, col="lightgrey")
boxplot(DEDD.results.10_DD$auc, names=NA, col=col_vector[1:2], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.43,0.87), main=paste0("AUC diff exp, DD10"), cex.main=1.5)
abline(v=2.5, col="lightgrey")
boxplot(DEDD.results.20_DD$auc, names=NA, col=col_vector[1:2], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.43,0.87), main=paste0("AUC diff exp, DD20"), cex.main=1.5)
abline(v=2.5, col="lightgrey")
boxplot(DEDD.results.50_DD$auc, names=NA, col=col_vector[1:2], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.43,0.87), main=paste0("AUC diff exp, DD50"), cex.main=1.5)
abline(v=2.5, col="lightgrey")
boxplot(DEDD.results.2_DE$auc, names=NA, col=col_vector[1:2], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.5,1), main=paste0("AUC diff exp, DE2"), cex.main=1.5)
abline(v=2.5, col="lightgrey")
boxplot(DEDD.results.5_DE$auc, names=NA, col=col_vector[1:2], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.5,1), main=paste0("AUC diff exp, DE5"), cex.main=1.5)
abline(v=2.5, col="lightgrey")
boxplot(DEDD.results.10_DE$auc, names=NA, col=col_vector[1:2], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.5,1), main=paste0("AUC diff exp, DE10"), cex.main=1.5)
abline(v=2.5, col="lightgrey")
boxplot(DEDD.results.20_DE$auc, names=NA, col=col_vector[1:2], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.5,1), main=paste0("AUC diff exp, DE20"), cex.main=1.5)
abline(v=2.5, col="lightgrey")
boxplot(DEDD.results.50_DE$auc, names=NA, col=col_vector[1:2], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.5,1), main=paste0("AUC diff exp, DE50"), cex.main=1.5)
abline(v=2.5, col="lightgrey")
boxplot(DEDD.results.2_DEDD$auc, names=NA, col=col_vector[1:2], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.5,0.94), main=paste0("AUC diff exp, DEDD2"), cex.main=1.5)
abline(v=2.5, col="lightgrey")
boxplot(DEDD.results.5_DEDD$auc, names=NA, col=col_vector[1:2], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.5,0.94), main=paste0("AUC diff exp, DEDD5"), cex.main=1.5)
abline(v=2.5, col="lightgrey")
boxplot(DEDD.results.10_DEDD$auc, names=NA, col=col_vector[1:2], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.5,0.94), main=paste0("AUC diff exp, DEDD10"), cex.main=1.5)
abline(v=2.5, col="lightgrey")
boxplot(DEDD.results.20_DEDD$auc, names=NA, col=col_vector[1:2], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.5,0.94), main=paste0("AUC diff exp, DEDD50"), cex.main=1.5)
abline(v=2.5, col="lightgrey")
boxplot(DEDD.results.50_DEDD$auc, names=NA, col=col_vector[1:2], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.5,0.94), main=paste0("AUC diff exp, DEDD50"), cex.main=1.5)
abline(v=2.5, col="lightgrey")
rbind(colMeans(DEDD.results.2_DD$auc), colMeans(DEDD.results.5_DD$auc), colMeans(DEDD.results.10_DD$auc), 
      colMeans(DEDD.results.20_DD$auc), colMeans(DEDD.results.50_DD$auc))
rbind(colMeans(DEDD.results.2_DE$auc), colMeans(DEDD.results.5_DE$auc), colMeans(DEDD.results.10_DE$auc), 
      colMeans(DEDD.results.20_DE$auc), colMeans(DEDD.results.50_DE$auc))
rbind(colMeans(DEDD.results.2_DEDD$auc), colMeans(DEDD.results.5_DEDD$auc), colMeans(DEDD.results.10_DEDD$auc), 
      colMeans(DEDD.results.20_DEDD$auc), colMeans(DEDD.results.50_DEDD$auc))
# AUCs always higher for expHMM than lnHMM with differences in dispersion only, and nearly always with differences in mean and 
# dispersion.

############
## pAUC ####
par(mfrow=c(3,5), mar=c(1,2,2,1), mgp=c(3,0.7,0))
boxplot(DEDD.results.2_DD$pauc, names=NA, col=col_vector[1:2], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.03), main=paste0("Partial AUC diff exp, DE2"), cex.main=1.5)
legend("topleft", fill=col_vector[1:2], bty='n', cex=1.5, ncol=2, legend=legend)
abline(v=2.5, col="lightgrey")
boxplot(DEDD.results.5_DD$pauc, names=NA, col=col_vector[1:2], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.03), main=paste0("Partial AUC diff exp, DE5"), cex.main=1.5)
abline(v=2.5, col="lightgrey")
boxplot(DEDD.results.10_DD$pauc, names=NA, col=col_vector[1:2], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.03), main=paste0("Partial AUC diff exp, DE10"), cex.main=1.5)
abline(v=2.5, col="lightgrey")
boxplot(DEDD.results.20_DD$pauc, names=NA, col=col_vector[1:2], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.03), main=paste0("Partial AUC diff exp, DE20"), cex.main=1.5)
abline(v=2.5, col="lightgrey")
boxplot(DEDD.results.50_DD$pauc, names=NA, col=col_vector[1:2], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.03), main=paste0("Partial AUC diff exp, DE50"), cex.main=1.5)
abline(v=2.5, col="lightgrey")
boxplot(DEDD.results.2_DE$pauc, names=NA, col=col_vector[1:2], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.048), main=paste0("Partial AUC diff exp, DE2"), cex.main=1.5)
abline(v=2.5, col="lightgrey")
boxplot(DEDD.results.5_DE$pauc, names=NA, col=col_vector[1:2], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.048), main=paste0("Partial AUC diff exp, DE5"), cex.main=1.5)
abline(v=2.5, col="lightgrey")
boxplot(DEDD.results.10_DE$pauc, names=NA, col=col_vector[1:2], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.048), main=paste0("Partial AUC diff exp, DE10"), cex.main=1.5)
abline(v=2.5, col="lightgrey")
boxplot(DEDD.results.20_DE$pauc, names=NA, col=col_vector[1:2], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.048), main=paste0("Partial AUC diff exp, DE20"), cex.main=1.5)
abline(v=2.5, col="lightgrey")
boxplot(DEDD.results.50_DE$pauc, names=NA, col=col_vector[1:2], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.048), main=paste0("Partial AUC diff exp, DE50"), cex.main=1.5)
abline(v=2.5, col="lightgrey")
boxplot(DEDD.results.2_DEDD$pauc, names=NA, col=col_vector[1:2], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.043), main=paste0("Partial AUC diff exp, DEDD2"), cex.main=1.5)
abline(v=2.5, col="lightgrey")
boxplot(DEDD.results.5_DEDD$pauc, names=NA, col=col_vector[1:2], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.043), main=paste0("Partial AUC diff exp, DEDD5"), cex.main=1.5)
abline(v=2.5, col="lightgrey")
boxplot(DEDD.results.10_DEDD$pauc, names=NA, col=col_vector[1:2], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.043), main=paste0("Partial AUC diff exp, DEDD10"), cex.main=1.5)
abline(v=2.5, col="lightgrey")
boxplot(DEDD.results.20_DEDD$pauc, names=NA, col=col_vector[1:2], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.043), main=paste0("Partial AUC diff exp, DEDD20"), cex.main=1.5)
abline(v=2.5, col="lightgrey")
boxplot(DEDD.results.50_DEDD$pauc, names=NA, col=col_vector[1:2], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.043), main=paste0("Partial AUC diff exp, DEDD50"), cex.main=1.5)
abline(v=2.5, col="lightgrey")
rbind(colMeans(DEDD.results.2_DD$pauc), colMeans(DEDD.results.5_DD$pauc), colMeans(DEDD.results.10_DD$pauc), 
      colMeans(DEDD.results.20_DD$pauc), colMeans(DEDD.results.50_DD$pauc))
rbind(colMeans(DEDD.results.2_DE$pauc), colMeans(DEDD.results.5_DE$pauc), colMeans(DEDD.results.10_DE$pauc), 
      colMeans(DEDD.results.20_DE$pauc), colMeans(DEDD.results.50_DE$pauc))
rbind(colMeans(DEDD.results.2_DEDD$pauc), colMeans(DEDD.results.5_DEDD$pauc), colMeans(DEDD.results.10_DEDD$pauc), 
      colMeans(DEDD.results.20_DEDD$pauc), colMeans(DEDD.results.50_DEDD$pauc))
# Partial AUCs higher for expHMM than lnHMM except  for small sample sizes with differences in mean only or mean and dispersion.

###########
## FDR ####
par(mfrow=c(3,5), mar=c(1,2,2,1), mgp=c(3,0.7,0))
boxplot(DEDD.results.2_DD$fdr, names=NA, col=col_vector[1:3], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.4), main=paste0("FDR diff exp, DD2"), cex.main=1.5)
abline(h=0.05, v=6.5, col='lightgrey')
legend("topleft", fill=col_vector[1:3], bty='n', cex=1.5, ncol=1, 
       legend=c("p = 0.5 threshold", "Posterior threshold", "bfdr"))
boxplot(DEDD.results.5_DD$fdr, names=NA, col=col_vector[1:3], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.4), main=paste0("FDR diff exp, DD5"), cex.main=1.5)
abline(h=0.05, v=6.5, col='lightgrey')
boxplot(DEDD.results.10_DD$fdr, names=NA, col=col_vector[1:3], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.4), main=paste0("FDR diff exp, DD10"), cex.main=1.5)
abline(h=0.05, v=6.5, col='lightgrey')
boxplot(DEDD.results.20_DD$fdr, names=NA, col=col_vector[1:3], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.4), main=paste0("FDR diff exp, DD20"), cex.main=1.5)
abline(h=0.05, v=6.5, col='lightgrey')
boxplot(DEDD.results.50_DD$fdr, names=NA, col=col_vector[1:3], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.4), main=paste0("FDR diff exp, DD50"), cex.main=1.5)
abline(h=0.05, v=6.5, col='lightgrey')
boxplot(DEDD.results.2_DE$fdr, names=NA, col=col_vector[1:3], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.4), main=paste0("FDR diff exp, DE2"), cex.main=1.5)
abline(h=0.05, v=6.5, col='lightgrey')
boxplot(DEDD.results.5_DE$fdr, names=NA, col=col_vector[1:3], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.4), main=paste0("FDR diff exp, DE5"), cex.main=1.5)
abline(h=0.05, v=6.5, col='lightgrey')
boxplot(DEDD.results.10_DE$fdr, names=NA, col=col_vector[1:3], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.4), main=paste0("FDR diff exp, DE10"), cex.main=1.5)
abline(h=0.05, v=6.5, col='lightgrey')
boxplot(DEDD.results.20_DE$fdr, names=NA, col=col_vector[1:3], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.4), main=paste0("FDR diff exp, DE20"), cex.main=1.5)
abline(h=0.05, v=6.5, col='lightgrey')
boxplot(DEDD.results.50_DE$fdr, names=NA, col=col_vector[1:3], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.4), main=paste0("FDR diff exp, DE50"), cex.main=1.5)
abline(h=0.05, v=6.5, col='lightgrey')
boxplot(DEDD.results.2_DEDD$fdr, names=NA, col=col_vector[1:3], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.4), main=paste0("FDR diff exp, DEDD2"), cex.main=1.5)
abline(h=0.05, v=6.5, col='lightgrey')
boxplot(DEDD.results.5_DEDD$fdr, names=NA, col=col_vector[1:3], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.4), main=paste0("FDR diff exp, DEDD5"), cex.main=1.5)
abline(h=0.05, v=6.5, col='lightgrey')
boxplot(DEDD.results.10_DEDD$fdr, names=NA, col=col_vector[1:3], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.4), main=paste0("FDR diff exp, DEDD10"), cex.main=1.5)
abline(h=0.05, v=6.5, col='lightgrey')
boxplot(DEDD.results.20_DEDD$fdr, names=NA, col=col_vector[1:3], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.4), main=paste0("FDR diff exp, DEDD20"), cex.main=1.5)
abline(h=0.05, v=6.5, col='lightgrey')
boxplot(DEDD.results.50_DEDD$fdr, names=NA, col=col_vector[1:3], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.4), main=paste0("FDR diff exp, DEDD50"), cex.main=1.5)
abline(h=0.05, v=6.5, col='lightgrey')
# Mean FDRs:
rbind(colMeans(DEDD.results.2_DD$fdr, na.rm=T), colMeans(DEDD.results.5_DD$fdr, na.rm=T), 
      colMeans(DEDD.results.10_DD$fdr, na.rm=T), colMeans(DEDD.results.20_DD$fdr, na.rm=T), 
      colMeans(DEDD.results.50_DD$fdr, na.rm=T))
rbind(colMeans(DEDD.results.2_DE$fdr, na.rm=T), colMeans(DEDD.results.5_DE$fdr, na.rm=T), 
      colMeans(DEDD.results.10_DE$fdr, na.rm=T), colMeans(DEDD.results.20_DE$fdr, na.rm=T), 
      colMeans(DEDD.results.50_DE$fdr, na.rm=T))
rbind(colMeans(DEDD.results.2_DEDD$fdr, na.rm=T), colMeans(DEDD.results.5_DEDD$fdr, na.rm=T), 
      colMeans(DEDD.results.10_DEDD$fdr, na.rm=T), colMeans(DEDD.results.20_DEDD$fdr, na.rm=T), 
      colMeans(DEDD.results.50_DEDD$fdr, na.rm=T))
# Mean squared distances from 0.05:
rbind(colMeans((DEDD.results.2_DD$fdr - 0.05)^2, na.rm=T), colMeans((DEDD.results.5_DD$fdr - 0.05)^2, na.rm=T), 
      colMeans((DEDD.results.10_DD$fdr - 0.05)^2, na.rm=T), colMeans((DEDD.results.20_DD$fdr - 0.05)^2, na.rm=T), 
      colMeans((DEDD.results.50_DD$fdr - 0.05)^2, na.rm=T))
rbind(colMeans((DEDD.results.2_DE$fdr - 0.05)^2, na.rm=T), colMeans((DEDD.results.5_DE$fdr - 0.05)^2, na.rm=T), 
      colMeans((DEDD.results.10_DE$fdr - 0.05)^2, na.rm=T), colMeans((DEDD.results.20_DE$fdr - 0.05)^2, na.rm=T), 
      colMeans((DEDD.results.50_DE$fdr - 0.05)^2, na.rm=T))
rbind(colMeans((DEDD.results.2_DEDD$fdr - 0.05)^2, na.rm=T), colMeans((DEDD.results.5_DEDD$fdr - 0.05)^2, na.rm=T), 
      colMeans((DEDD.results.10_DEDD$fdr - 0.05)^2, na.rm=T), colMeans((DEDD.results.20_DEDD$fdr - 0.05)^2, na.rm=T), 
      colMeans((DEDD.results.50_DEDD$fdr - 0.05)^2, na.rm=T))
# BFDR generally closer to 0.05 than posterior threshold method, although posterior threshold method doesn't aim for any particular 
# FDR. FDRs are pretty poor generally and don't improve with sample size. FDRs nearly always closer to 0.05 for expHMM than lnHMM.

###########
## TPR ####
par(mfrow=c(3,5), mar=c(1,2,2,1), mgp=c(3,0.7,0))
boxplot(DEDD.results.2_DD$tpr, names=NA, col=col_vector[1:3], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.45), main=paste0("TPR diff exp, DE2"), cex.main=1.5)
legend("topleft", fill=col_vector[1:3], bty='n', cex=1.5, ncol=1, 
       legend=c("p = 0.5 threshold", "Posterior threshold", "bfdr"))
abline(v=6.5, col="lightgrey")
boxplot(DEDD.results.5_DD$tpr, names=NA, col=col_vector[1:3], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.45), main=paste0("TPR diff exp, DE5"), cex.main=1.5)
abline(v=6.5, col="lightgrey")
boxplot(DEDD.results.10_DD$tpr, names=NA, col=col_vector[1:3], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.45), main=paste0("TPR diff exp, DE10"), cex.main=1.5)
abline(v=6.5, col="lightgrey")
boxplot(DEDD.results.20_DD$tpr, names=NA, col=col_vector[1:3], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.45), main=paste0("TPR diff exp, DE20"), cex.main=1.5)
abline(v=6.5, col="lightgrey")
boxplot(DEDD.results.50_DD$tpr, names=NA, col=col_vector[1:3], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.45), main=paste0("TPR diff exp, DE50"), cex.main=1.5)
abline(v=6.5, col="lightgrey")
boxplot(DEDD.results.2_DE$tpr, names=NA, col=col_vector[1:3], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.9), main=paste0("TPR diff exp, DE2"), cex.main=1.5)
abline(v=6.5, col="lightgrey")
boxplot(DEDD.results.5_DE$tpr, names=NA, col=col_vector[1:3], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.9), main=paste0("TPR diff exp, DE5"), cex.main=1.5)
abline(v=6.5, col="lightgrey")
boxplot(DEDD.results.10_DE$tpr, names=NA, col=col_vector[1:3], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.9), main=paste0("TPR diff exp, DE10"), cex.main=1.5)
abline(v=6.5, col="lightgrey")
boxplot(DEDD.results.20_DE$tpr, names=NA, col=col_vector[1:3], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.9), main=paste0("TPR diff exp, DE20"), cex.main=1.5)
abline(v=6.5, col="lightgrey")
boxplot(DEDD.results.50_DE$tpr, names=NA, col=col_vector[1:3], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.9), main=paste0("TPR diff exp, DE50"), cex.main=1.5)
abline(v=6.5, col="lightgrey")
boxplot(DEDD.results.2_DEDD$tpr, names=NA, col=col_vector[1:3], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.75), main=paste0("TPR diff exp, DEDD2"), cex.main=1.5)
abline(v=6.5, col="lightgrey")
boxplot(DEDD.results.5_DEDD$tpr, names=NA, col=col_vector[1:3], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.75), main=paste0("TPR diff exp, DEDD5"), cex.main=1.5)
abline(v=6.5, col="lightgrey")
boxplot(DEDD.results.10_DEDD$tpr, names=NA, col=col_vector[1:3], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.75), main=paste0("TPR diff exp, DEDD10"), cex.main=1.5)
abline(v=6.5, col="lightgrey")
boxplot(DEDD.results.20_DEDD$tpr, names=NA, col=col_vector[1:3], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.75), main=paste0("TPR diff exp, DEDD20"), cex.main=1.5)
abline(v=6.5, col="lightgrey")
boxplot(DEDD.results.50_DEDD$tpr, names=NA, col=col_vector[1:3], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.75), main=paste0("TPR diff exp, DEDD50"), cex.main=1.5)
abline(v=6.5, col="lightgrey")
rbind(colMeans(DEDD.results.2_DD$tpr), colMeans(DEDD.results.5_DD$tpr), colMeans(DEDD.results.10_DD$tpr), 
      colMeans(DEDD.results.20_DD$tpr), colMeans(DEDD.results.50_DD$tpr))
rbind(colMeans(DEDD.results.2_DE$tpr), colMeans(DEDD.results.5_DE$tpr), colMeans(DEDD.results.10_DE$tpr), 
      colMeans(DEDD.results.20_DE$tpr), colMeans(DEDD.results.50_DE$tpr))
rbind(colMeans(DEDD.results.2_DEDD$tpr), colMeans(DEDD.results.5_DEDD$tpr), colMeans(DEDD.results.10_DEDD$tpr), 
      colMeans(DEDD.results.20_DEDD$tpr), colMeans(DEDD.results.50_DEDD$tpr))
# Consistently highest TPRs using posterior threshold method except sometimes for 50 samples per group, and higher for lnHMM than 
# expHMM, but at expense of higher FDRs.

#############################
## False discovery plots ####
par(mfrow=c(3,5), mar=c(2,2,2,1), mgp=c(3,0.7,0))

plot(DEDD.results.2_DD$mean.discoveries$expHMM_blood, DEDD.results.2_DD$mean.fdr$expHMM_blood, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], main=paste0("DD2"), cex.main=1.5)
lines(DEDD.results.2_DD$mean.discoveries$lnHMM_blood, DEDD.results.2_DD$mean.fdr$lnHMM_blood, 
      col=col_vector[4])
lines(DEDD.results.2_DD$mean.discoveries$expHMM_muscle, DEDD.results.2_DD$mean.fdr$expHMM_muscle, 
      lty=2, col=col_vector[1])
lines(DEDD.results.2_DD$mean.discoveries$lnHMM_muscle, DEDD.results.2_DD$mean.fdr$lnHMM_muscle, 
      lty=2, col=col_vector[4])
abline(h=0.05, col='lightgrey')
plot(DEDD.results.5_DD$mean.discoveries$expHMM_blood, DEDD.results.5_DD$mean.fdr$expHMM_blood, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], main=paste0("DD5"), cex.main=1.5)
lines(DEDD.results.5_DD$mean.discoveries$lnHMM_blood, DEDD.results.5_DD$mean.fdr$lnHMM_blood, 
      col=col_vector[4])
lines(DEDD.results.5_DD$mean.discoveries$expHMM_muscle, DEDD.results.5_DD$mean.fdr$expHMM_muscle, 
      lty=2, col=col_vector[1])
lines(DEDD.results.5_DD$mean.discoveries$lnHMM_muscle, DEDD.results.5_DD$mean.fdr$lnHMM_muscle, 
      lty=2, col=col_vector[4])
abline(h=0.05, col='lightgrey')
plot(DEDD.results.10_DD$mean.discoveries$expHMM_blood, DEDD.results.10_DD$mean.fdr$expHMM_blood, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], main=paste0("DD10"), cex.main=1.5)
lines(DEDD.results.10_DD$mean.discoveries$lnHMM_blood, DEDD.results.10_DD$mean.fdr$lnHMM_blood, 
      col=col_vector[4])
lines(DEDD.results.10_DD$mean.discoveries$expHMM_muscle, DEDD.results.10_DD$mean.fdr$expHMM_muscle, 
      lty=2, col=col_vector[1])
lines(DEDD.results.10_DD$mean.discoveries$lnHMM_muscle, DEDD.results.10_DD$mean.fdr$lnHMM_muscle, 
      lty=2, col=col_vector[4])
abline(h=0.05, col='lightgrey')
plot(DEDD.results.20_DD$mean.discoveries$expHMM_blood, DEDD.results.20_DD$mean.fdr$expHMM_blood, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], main=paste0("DD20"), cex.main=1.5)
lines(DEDD.results.20_DD$mean.discoveries$lnHMM_blood, DEDD.results.20_DD$mean.fdr$lnHMM_blood, 
      col=col_vector[4])
lines(DEDD.results.20_DD$mean.discoveries$expHMM_muscle, DEDD.results.20_DD$mean.fdr$expHMM_muscle, 
      lty=2, col=col_vector[1])
lines(DEDD.results.20_DD$mean.discoveries$lnHMM_muscle, DEDD.results.20_DD$mean.fdr$lnHMM_muscle, 
      lty=2, col=col_vector[4])
abline(h=0.05, col='lightgrey')
plot(DEDD.results.50_DD$mean.discoveries$expHMM_blood, DEDD.results.50_DD$mean.fdr$expHMM_blood, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], main=paste0("DD50"), cex.main=1.5)
lines(DEDD.results.50_DD$mean.discoveries$lnHMM_blood, DEDD.results.50_DD$mean.fdr$lnHMM_blood, 
      col=col_vector[4])
lines(DEDD.results.50_DD$mean.discoveries$expHMM_muscle, DEDD.results.50_DD$mean.fdr$expHMM_muscle, 
      lty=2, col=col_vector[1])
lines(DEDD.results.50_DD$mean.discoveries$lnHMM_muscle, DEDD.results.50_DD$mean.fdr$lnHMM_muscle, 
      lty=2, col=col_vector[4])
abline(h=0.05, col='lightgrey')
legend("topleft", bty='n', col=col_vector[c(1,4)], lty=1, ncol=1, cex=1.5, legend=legend)

plot(DEDD.results.2_DE$mean.discoveries$expHMM_blood, DEDD.results.2_DE$mean.fdr$expHMM_blood, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], main=paste0("DE2"), cex.main=1.5)
lines(DEDD.results.2_DE$mean.discoveries$lnHMM_blood, DEDD.results.2_DE$mean.fdr$lnHMM_blood, 
      col=col_vector[4])
lines(DEDD.results.2_DE$mean.discoveries$expHMM_muscle, DEDD.results.2_DE$mean.fdr$expHMM_muscle, 
      lty=2, col=col_vector[1])
lines(DEDD.results.2_DE$mean.discoveries$lnHMM_muscle, DEDD.results.2_DE$mean.fdr$lnHMM_muscle, 
      lty=2, col=col_vector[4])
abline(h=0.05, col='lightgrey')
plot(DEDD.results.5_DE$mean.discoveries$expHMM_blood, DEDD.results.5_DE$mean.fdr$expHMM_blood, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], main=paste0("DE5"), cex.main=1.5)
lines(DEDD.results.5_DE$mean.discoveries$lnHMM_blood, DEDD.results.5_DE$mean.fdr$lnHMM_blood, 
      col=col_vector[4])
lines(DEDD.results.5_DE$mean.discoveries$expHMM_muscle, DEDD.results.5_DE$mean.fdr$expHMM_muscle, 
      lty=2, col=col_vector[1])
lines(DEDD.results.5_DE$mean.discoveries$lnHMM_muscle, DEDD.results.5_DE$mean.fdr$lnHMM_muscle, 
      lty=2, col=col_vector[4])
abline(h=0.05, col='lightgrey')
plot(DEDD.results.10_DE$mean.discoveries$expHMM_blood, DEDD.results.10_DE$mean.fdr$expHMM_blood, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], main=paste0("DE10"), cex.main=1.5)
lines(DEDD.results.10_DE$mean.discoveries$lnHMM_blood, DEDD.results.10_DE$mean.fdr$lnHMM_blood, 
      col=col_vector[4])
lines(DEDD.results.10_DE$mean.discoveries$expHMM_muscle, DEDD.results.10_DE$mean.fdr$expHMM_muscle, 
      lty=2, col=col_vector[1])
lines(DEDD.results.10_DE$mean.discoveries$lnHMM_muscle, DEDD.results.10_DE$mean.fdr$lnHMM_muscle, 
      lty=2, col=col_vector[4])
abline(h=0.05, col='lightgrey')
plot(DEDD.results.20_DE$mean.discoveries$expHMM_blood, DEDD.results.20_DE$mean.fdr$expHMM_blood, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], main=paste0("DE20"), cex.main=1.5)
lines(DEDD.results.20_DE$mean.discoveries$lnHMM_blood, DEDD.results.20_DE$mean.fdr$lnHMM_blood, 
      col=col_vector[4])
lines(DEDD.results.20_DE$mean.discoveries$expHMM_muscle, DEDD.results.20_DE$mean.fdr$expHMM_muscle, 
      lty=2, col=col_vector[1])
lines(DEDD.results.20_DE$mean.discoveries$lnHMM_muscle, DEDD.results.20_DE$mean.fdr$lnHMM_muscle, 
      lty=2, col=col_vector[4])
abline(h=0.05, col='lightgrey')
plot(DEDD.results.50_DE$mean.discoveries$expHMM_blood, DEDD.results.50_DE$mean.fdr$expHMM_blood, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], main=paste0("DD50"), cex.main=1.5)
lines(DEDD.results.50_DE$mean.discoveries$lnHMM_blood, DEDD.results.50_DE$mean.fdr$lnHMM_blood, 
      col=col_vector[4])
lines(DEDD.results.50_DE$mean.discoveries$expHMM_muscle, DEDD.results.50_DE$mean.fdr$expHMM_muscle, 
      lty=2, col=col_vector[1])
lines(DEDD.results.50_DE$mean.discoveries$lnHMM_muscle, DEDD.results.50_DE$mean.fdr$lnHMM_muscle, 
      lty=2, col=col_vector[4])
abline(h=0.05, col='lightgrey')
legend("topleft", bty='n', col=col_vector[c(1,4)], lty=1, ncol=1, cex=1.5, legend=legend)

plot(DEDD.results.2_DEDD$mean.discoveries$expHMM_blood, DEDD.results.2_DEDD$mean.fdr$expHMM_blood, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], main=paste0("DEDD2"), cex.main=1.5)
lines(DEDD.results.2_DEDD$mean.discoveries$lnHMM_blood, DEDD.results.2_DEDD$mean.fdr$lnHMM_blood, 
      col=col_vector[4])
lines(DEDD.results.2_DEDD$mean.discoveries$expHMM_muscle, DEDD.results.2_DEDD$mean.fdr$expHMM_muscle, 
      lty=2, col=col_vector[1])
lines(DEDD.results.2_DEDD$mean.discoveries$lnHMM_muscle, DEDD.results.2_DEDD$mean.fdr$lnHMM_muscle, 
      lty=2, col=col_vector[4])
abline(h=0.05, col='lightgrey')
plot(DEDD.results.5_DEDD$mean.discoveries$expHMM_blood, DEDD.results.5_DEDD$mean.fdr$expHMM_blood, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], main=paste0("DEDD5"), cex.main=1.5)
lines(DEDD.results.5_DEDD$mean.discoveries$lnHMM_blood, DEDD.results.5_DEDD$mean.fdr$lnHMM_blood, 
      col=col_vector[4])
lines(DEDD.results.5_DEDD$mean.discoveries$expHMM_muscle, DEDD.results.5_DEDD$mean.fdr$expHMM_muscle, 
      lty=2, col=col_vector[1])
lines(DEDD.results.5_DEDD$mean.discoveries$lnHMM_muscle, DEDD.results.5_DEDD$mean.fdr$lnHMM_muscle, 
      lty=2, col=col_vector[4])
abline(h=0.05, col='lightgrey')
plot(DEDD.results.10_DEDD$mean.discoveries$expHMM_blood, DEDD.results.10_DEDD$mean.fdr$expHMM_blood, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], main=paste0("DEDD10"), cex.main=1.5)
lines(DEDD.results.10_DEDD$mean.discoveries$lnHMM_blood, DEDD.results.10_DEDD$mean.fdr$lnHMM_blood, 
      col=col_vector[4])
lines(DEDD.results.10_DEDD$mean.discoveries$expHMM_muscle, DEDD.results.10_DEDD$mean.fdr$expHMM_muscle, 
      lty=2, col=col_vector[1])
lines(DEDD.results.10_DEDD$mean.discoveries$lnHMM_muscle, DEDD.results.10_DEDD$mean.fdr$lnHMM_muscle, 
      lty=2, col=col_vector[4])
abline(h=0.05, col='lightgrey')
plot(DEDD.results.20_DEDD$mean.discoveries$expHMM_blood, DEDD.results.20_DEDD$mean.fdr$expHMM_blood, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], main=paste0("DEDD20"), cex.main=1.5)
lines(DEDD.results.20_DEDD$mean.discoveries$lnHMM_blood, DEDD.results.20_DEDD$mean.fdr$lnHMM_blood, 
      col=col_vector[4])
lines(DEDD.results.20_DEDD$mean.discoveries$expHMM_muscle, DEDD.results.20_DEDD$mean.fdr$expHMM_muscle, 
      lty=2, col=col_vector[1])
lines(DEDD.results.20_DEDD$mean.discoveries$lnHMM_muscle, DEDD.results.20_DEDD$mean.fdr$lnHMM_muscle, 
      lty=2, col=col_vector[4])
abline(h=0.05, col='lightgrey')
plot(DEDD.results.50_DEDD$mean.discoveries$expHMM_blood, DEDD.results.50_DEDD$mean.fdr$expHMM_blood, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], main=paste0("DEDD50"), cex.main=1.5)
lines(DEDD.results.50_DEDD$mean.discoveries$lnHMM_blood, DEDD.results.50_DEDD$mean.fdr$lnHMM_blood, 
      col=col_vector[4])
lines(DEDD.results.50_DEDD$mean.discoveries$expHMM_muscle, DEDD.results.50_DEDD$mean.fdr$expHMM_muscle, 
      lty=2, col=col_vector[1])
lines(DEDD.results.50_DEDD$mean.discoveries$lnHMM_muscle, DEDD.results.50_DEDD$mean.fdr$lnHMM_muscle, 
      lty=2, col=col_vector[4])
abline(h=0.05, col='lightgrey')
legend("topleft", bty='n', col=col_vector[c(1,4)], lty=1, ncol=1, cex=1.5, legend=legend)
# FDR curves consistently lower for expHMM than lnHMM for differences in dispersion only, but very little difference between them 
# for differences in mean only or mean and dispersion.
