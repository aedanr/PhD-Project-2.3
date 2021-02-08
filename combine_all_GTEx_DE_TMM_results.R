# Create overall combined results which can be subsetted for different comparisons (e.g. 
# blood and muscle together or separately), for DE data and TMM only.

library(here)

# Import blood results
folder <- "Results/GTEx blood simulated DE, DD, DEDD results Feb 2020"
for (i in c("2", "5", "10", "20", "50")) {
  assign(paste0("DE.results.blood", i), 
         readRDS(here(folder, paste0("DE.results.blood_", i, "_DE.rds"))))
  assign(paste0("DEDD.results.blood", i), 
         readRDS(here(folder, paste0("DEDD.results.blood_", i, "_DE.rds"))))
}

# Import muscle results
folder <- "Results/GTEx muscle simulated DE, DD, DEDD results March 2020"
for (i in c("2", "5", "10", "20", "50")) {
  assign(paste0("DE.results.muscle", i), 
         readRDS(here(folder, paste0("DE.results.muscle_", i, "_DE.rds"))))
  assign(paste0("DEDD.results.muscle", i), 
         readRDS(here(folder, paste0("DEDD.results.muscle_", i, "_DE.rds"))))
}
rm(i,folder)

# Combine all results for each sample size and analysis type, keeping only TMM
for (size in c("2", "5", "10", "20", "50")) {
  for (i in c("pauc", "auc", "mean.fdr", "mean.discoveries")) {
    temp <- data.frame(
      cbind(
        get(paste0("DE.results.blood", size))[[i]], 
        get(paste0("DE.results.muscle", size))[[i]]
      )
    )
    temp <- temp[, c(1:14, c(1:14) + 28)]
    names(temp) <- c(
      'blood_eR.ql', 'blood_eR.lr', 'blood_eR.et', 'blood_DES.noif', 'blood_DES.if', 
      'blood_voom', 'blood_DSS', 'blood_baySeq', 'blood_MD.zi', 'blood_MD.nozi', 
      'blood_expHM.untr', 'blood_expHM.log', 'blood_lnHM.untr', 'blood_lnHM.log', 
      'muscle_eR.ql', 'muscle_eR.lr', 'muscle_eR.et', 'muscle_DES.noif', 'muscle_DES.if', 
      'muscle_voom', 'muscle_DSS', 'muscle_baySeq', 'muscle_MD.zi', 'muscle_MD.nozi', 
      'muscle_expHM.untr', 'muscle_expHM.log', 'muscle_lnHM.untr', 'muscle_lnHM.log'
    )
    assign(i, temp)
  }
  for (i in c("fpr", "fdr", "tpr")) {
    temp <- data.frame(
      cbind(
        get(paste0("DE.results.blood", size))[[i]], 
        get(paste0("DE.results.muscle", size))[[i]]
      )
    )
    temp <- temp[, c(1:15, c(1:15) + 30)]
    names(temp) <- c(
      'blood_eR.ql', 'blood_eR.lr', 'blood_eR.et', 'blood_DES.noif', 'blood_DES.if', 
      'blood_voom', 'blood_DSS', 'blood_DSS.lfdr', 'blood_baySeq', 
      'blood_MD.zi', 'blood_MD.nozi', 
      'blood_expHM.untr', 'blood_expHM.log', 'blood_lnHM.untr', 'blood_lnHM.log', 
      'muscle_eR.ql', 'muscle_eR.lr', 'muscle_eR.et', 'muscle_DES.noif', 'muscle_DES.if', 
      'muscle_voom', 'muscle_DSS', 'muscle_DSS.lfdr', 'muscle_baySeq', 
      'muscle_MD.zi', 'muscle_MD.nozi', 
      'muscle_expHM.untr', 'muscle_expHM.log', 'muscle_lnHM.untr', 'muscle_lnHM.log'
    )
    assign(i, temp)
  }
  assign(
    paste0(
      "DE.results.DE", size), 
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
        get(paste0("DEDD.results.muscle", size))[[i]]
      )
    )
    temp <- temp[, c(1,3,5,7)]
    names(temp) <- c(
      'blood_expHMM', 'blood_lnHMM', 'muscle_expHMM', 'muscle_lnHMM'
    )
    assign(i, temp)
  }
  for (i in c("fpr", "fdr", "tpr")) {
    temp <- data.frame(
      cbind(
        get(paste0("DEDD.results.blood", size))[[i]], 
        get(paste0("DEDD.results.muscle", size))[[i]]
      )
    )
    temp <- temp[, c(1:3, 7:9, 13:15, 19:21)]
    names(temp) <- c(
      'blood_expHMM.p5', 'blood_expHMM.thr', 'blood_expHMM.bfdr', 
      'blood_lnHMM.p5', 'blood_lnHMM.thr', 'blood_lnHMM.bfdr', 
      'muscle_expHMM.p5', 'muscle_expHMM.thr', 'muscle_expHMM.bfdr', 
      'muscle_lnHMM.p5', 'muscle_lnHMM.thr', 'muscle_lnHMM.bfdr'
    )
    assign(i, temp)
  }
  assign(
    paste0(
      "DEDD.results.DE", size), 
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
for (i in c("DE", "DEDD")) {
  for (j in c("2", "5", "10", "20", "50")) {
    saveRDS(get(paste0(i, ".results.DE", j)), 
            here("Results/GTEx combined results Sept 2020", 
                 paste0(i, ".results.DE", j, ".rds")))
  }
}
