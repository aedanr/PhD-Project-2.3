# Create overall combined results which can be subsetted for different comparisons (e.g. 
# blood and muscle together or separately, short v long chain HMs), for DEDD data only.
# TMM only, hybrids using log-transformed HM results and voom instead of edgeR.

library(here)

# Import blood results
folder <- "Results/GTEx blood updated results Sept 2020"
for (i in c("2", "5", "10", "20", "50")) {
  assign(paste0("DD.results.blood", i), 
         readRDS(here(folder, paste0("DD.results.blood_", i, "_DEDD.rds"))))
  assign(paste0("DE.results.blood", i), 
         readRDS(here(folder, paste0("DE.results.blood_", i, "_DEDD.rds"))))
  assign(paste0("DEDD.results.blood", i), 
         readRDS(here(folder, paste0("DEDD.results.blood_", i, "_DEDD.rds"))))
}

# Import muscle results - original and long chain lnHM
folder <- "Results/GTEx muscle updated results Sept 2020"
for (i in c("2", "5", "10", "20", "50")) {
  assign(paste0("DD.results.muscle", i), 
         readRDS(here(folder, paste0("DD.results.muscle_", i, "_DEDD.rds"))))
  assign(paste0("DE.results.muscle", i), 
         readRDS(here(folder, paste0("DE.results.muscle_", i, "_DEDD.rds"))))
  assign(paste0("DEDD.results.muscle", i), 
         readRDS(here(folder, paste0("DEDD.results.muscle_", i, "_DEDD.rds"))))
}
folder <- "Results/GTEx muscle long chain results Aug 2020"
for (i in c("2", "5", "10", "20", "50")) {
  assign(paste0("DD.results.muscle_long", i), 
         readRDS(here(folder, paste0("DD.results.lnHM.muscle_", i, "_DEDD.rds"))))
  assign(paste0("DE.results.muscle_long", i), 
         readRDS(here(folder, paste0("DE.results.lnHM.muscle_", i, "_DEDD.rds"))))
  assign(paste0("DEDD.results.muscle_long", i), 
         readRDS(here(folder, paste0("DEDD.results.lnHM.muscle_", i, "_DEDD.rds"))))
}
rm(i,folder)

# Combine all results for each sample size and analysis type
for (size in c("2", "5", "10", "20", "50")) {
  for (i in c("pauc", "auc", "fpr", "fdr", "tpr", "mean.fdr", "mean.discoveries")) {
    temp <- data.frame(
      cbind(
        get(paste0("DD.results.blood", size))[[i]], 
        get(paste0("DD.results.muscle", size))[[i]], 
        get(paste0("DD.results.muscle_long", size))[[i]]
      )
    )
    names(temp) <- c(
      'blood_dV', 'blood_MD.zi', 'blood_MD.nozi', 'blood_expHM.untr', 'blood_expHM.log', 
      'blood_lnHM.untr', 'blood_lnHM.log', 
      'muscle_dV', 'muscle_MD.zi', 'muscle_MD.nozi', 'muscle_expHM.untr', 'muscle_expHM.log', 
      'muscle_lnHM.untr', 'muscle_lnHM.log', 'muscle_lnHM.untr.lc', 'muscle_lnHM.log.lc'
    )
    assign(i, temp)
  }
  assign(
    paste0(
      "DD.results.DEDD", size), 
    list(pauc = pauc, 
         auc = auc, 
         fpr = fpr, 
         fdr = fdr, 
         tpr = tpr, 
         mean.fdr = mean.fdr, 
         mean.discoveries = mean.discoveries
    )
  )
}

for (size in c("2", "5", "10", "20", "50")) {
  for (i in c("pauc", "auc", "mean.fdr", "mean.discoveries")) {
    temp <- data.frame(
      cbind(
        get(paste0("DE.results.blood", size))[[i]], 
        get(paste0("DE.results.muscle", size))[[i]], 
        get(paste0("DE.results.muscle_long", size))[[i]]
      )
    )
    names(temp) <- c(
      'blood_eR.ql', 'blood_eR.lr', 'blood_eR.et', 'blood_DES.noif', 'blood_DES.if', 
      'blood_voom', 'blood_DSS', 'blood_baySeq', 'blood_MD.zi', 'blood_MD.nozi', 
      'blood_expHM.untr', 'blood_expHM.log', 'blood_lnHM.untr', 'blood_lnHM.log', 
      'muscle_eR.ql', 'muscle_eR.lr', 'muscle_eR.et', 'muscle_DES.noif', 'muscle_DES.if', 
      'muscle_voom', 'muscle_DSS', 'muscle_baySeq', 'muscle_MD.zi', 'muscle_MD.nozi', 
      'muscle_expHM.untr', 'muscle_expHM.log', 'muscle_lnHM.untr', 'muscle_lnHM.log', 
      'muscle_lnHM.untr.lc', 'muscle_lnHM.log.lc'
    )
    assign(i, temp)
  }
  for (i in c("fpr", "fdr", "tpr")) {
    temp <- data.frame(
      cbind(
        get(paste0("DE.results.blood", size))[[i]], 
        get(paste0("DE.results.muscle", size))[[i]], 
        get(paste0("DE.results.muscle_long", size))[[i]]
      )
    )
    names(temp) <- c(
      'blood_eR.ql', 'blood_eR.lr', 'blood_eR.et', 'blood_DES.noif', 'blood_DES.if', 
      'blood_voom', 'blood_DSS', 'blood_DSS.lfdr', 'blood_baySeq', 
      'blood_MD.zi', 'blood_MD.nozi', 
      'blood_expHM.untr', 'blood_expHM.log', 'blood_lnHM.untr', 'blood_lnHM.log', 
      'muscle_eR.ql', 'muscle_eR.lr', 'muscle_eR.et', 'muscle_DES.noif', 'muscle_DES.if', 
      'muscle_voom', 'muscle_DSS', 'muscle_DSS.lfdr', 'muscle_baySeq', 
      'muscle_MD.zi', 'muscle_MD.nozi', 
      'muscle_expHM.untr', 'muscle_expHM.log', 'muscle_lnHM.untr', 'muscle_lnHM.log', 
      'muscle_lnHM.untr.lc', 'muscle_lnHM.log.lc'
    )
    assign(i, temp)
  }
  assign(
    paste0(
      "DE.results.DEDD", size), 
    list(pauc = pauc, 
         auc = auc, 
         fpr = fpr, 
         fdr = fdr, 
         tpr = tpr, 
         mean.fdr = mean.fdr, 
         mean.discoveries = mean.discoveries
    )
  )
}

for (size in c("2", "5", "10", "20", "50")) {
  for (i in c("pauc", "auc", "mean.fdr", "mean.discoveries")) {
    temp <- data.frame(
      cbind(
        get(paste0("DEDD.results.blood", size))[[i]], 
        get(paste0("DEDD.results.muscle", size))[[i]], 
        get(paste0("DEDD.results.muscle_long", size))[[i]]
      )
    )
    names(temp) <- c(
      'blood_dV', 'blood_expHMM', 'blood_lnHMM', 
      'blood_voom_dV', 'blood_voom_MD', 'blood_voom_HM', 
      'blood_HM_dV', 'blood_HM_MD', 'blood_HM_HM', 
      'muscle_dV', 'muscle_expHMM', 'muscle_lnHMM', 
      'muscle_voom_dV', 'muscle_voom_MD', 'muscle_voom_HM', 
      'muscle_HM_dV', 'muscle_HM_MD', 'muscle_HM_HM', 
      'muscle_lnHMM.lc', 'muscle_voom_HM.lc', 
      'muscle_HM.lc_dV', 'muscle_HM.lc_MD', 'muscle_HM_HM.lc'
    )
    assign(i, temp)
  }
  for (i in c("fpr", "fdr", "tpr")) {
    temp <- data.frame(
      cbind(
        get(paste0("DEDD.results.blood", size))[[i]], 
        get(paste0("DEDD.results.muscle", size))[[i]], 
        get(paste0("DEDD.results.muscle_long", size))[[i]]
      )
    )
    names(temp) <- c(
      'blood_dV', 'blood_expHMM.p5', 'blood_expHMM.thr', 'blood_expHMM.bfdr', 
      'blood_lnHMM.p5', 'blood_lnHMM.thr', 'blood_lnHMM.bfdr', 
      'blood_voom_dV', 'blood_voom_MD', 'blood_voom_HM', 
      'blood_HM_dV', 'blood_HM_MD', 'blood_HM_HM', 
      'muscle_dV', 'muscle_expHMM.p5', 'muscle_expHMM.thr', 'muscle_expHMM.bfdr', 
      'muscle_lnHMM.p5', 'muscle_lnHMM.thr', 'muscle_lnHMM.bfdr', 
      'muscle_voom_dV', 'muscle_voom_MD', 'muscle_voom_HM', 
      'muscle_HM_dV', 'muscle_HM_MD', 'muscle_HM_HM', 
      'muscle_lnHMM.p5.lc', 'muscle_lnHMM.thr.lc', 'muscle_lnHMM.bfdr.lc', 
      'muscle_voom_HM.lc', 'muscle_HM.lc_dV', 'muscle_HM.lc_MD', 'muscle_HM_HM.lc'
    )
    assign(i, temp)
  }
  assign(
    paste0(
      "DEDD.results.DEDD", size), 
    list(pauc = pauc, 
         auc = auc, 
         fpr = fpr, 
         fdr = fdr, 
         tpr = tpr, 
         mean.fdr = mean.fdr, 
         mean.discoveries = mean.discoveries
    )
  )
}

# Save results
for (i in c("DE", "DD", "DEDD")) {
  for (j in c("2", "5", "10", "20", "50")) {
    saveRDS(get(paste0(i, ".results.DEDD", j)), 
            here("Results/GTEx combined results Sept 2020", 
                 paste0(i, ".results.DEDD", j, ".rds")))
  }
}
