# Create overall combined results which can be subsetted for different comparisons (e.g. 
# blood and muscle together or separately), for DD data and TMM only.

library(here)

# Import blood results
folder <- "Results/GTEx blood simulated DE, DD, DEDD results Feb 2020"
for (i in c("2", "5", "10", "20", "50")) {
  assign(paste0("DD.results.blood", i), 
         readRDS(here(folder, paste0("DD.results.blood_", i, "_DD.rds"))))
  assign(paste0("DEDD.results.blood", i), 
         readRDS(here(folder, paste0("DEDD.results.blood_", i, "_DD.rds"))))
}

# Import muscle results - original and long chain lnHM
folder <- "Results/GTEx muscle simulated DE, DD, DEDD results March 2020"
for (i in c("2", "5", "10", "20", "50")) {
  assign(paste0("DD.results.muscle", i), 
         readRDS(here(folder, paste0("DD.results.muscle_", i, "_DD.rds"))))
  assign(paste0("DEDD.results.muscle", i), 
         readRDS(here(folder, paste0("DEDD.results.muscle_", i, "_DD.rds"))))
}
rm(i,folder)

# Combine all results for each sample size and analysis type
for (size in c("2", "5", "10", "20", "50")) {
  for (i in c("pauc", "auc", "fpr", "fdr", "tpr", "mean.fdr", "mean.discoveries")) {
    temp <- data.frame(
      cbind(
        get(paste0("DD.results.blood", size))[[i]], 
        get(paste0("DD.results.muscle", size))[[i]]
      )
    )
    temp <- temp[, c(1:7, 15:21)]
    names(temp) <- c(
      'blood_dV', 'blood_MD.zi', 'blood_MD.nozi', 'blood_expHM.untr', 'blood_expHM.log', 
      'blood_lnHM.untr', 'blood_lnHM.log', 
      'muscle_dV', 'muscle_MD.zi', 'muscle_MD.nozi', 'muscle_expHM.untr', 'muscle_expHM.log', 
      'muscle_lnHM.untr', 'muscle_lnHM.log'
    )
    assign(i, temp)
  }
  assign(
    paste0(
      "DD.results.DD", size), 
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
    temp <- temp[, c(1:3, 7:9)]
    names(temp) <- c(
      'blood_dV', 'blood_expHMM', 'blood_lnHMM', 
      'muscle_dV', 'muscle_expHMM', 'muscle_lnHMM'
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
    temp <- temp[, c(1:7, 15:21)]
    names(temp) <- c(
      'blood_dV', 'blood_expHMM.p5', 'blood_expHMM.thr', 'blood_expHMM.bfdr', 
      'blood_lnHMM.p5', 'blood_lnHMM.thr', 'blood_lnHMM.bfdr', 
      'muscle_dV', 'muscle_expHMM.p5', 'muscle_expHMM.thr', 'muscle_expHMM.bfdr', 
      'muscle_lnHMM.p5', 'muscle_lnHMM.thr', 'muscle_lnHMM.bfdr'
    )
    assign(i, temp)
  }
  assign(
    paste0(
      "DEDD.results.DD", size), 
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
for (i in c("DD", "DEDD")) {
  for (j in c("2", "5", "10", "20", "50")) {
    saveRDS(get(paste0(i, ".results.DD", j)), 
            here("Results/GTEx combined results Sept 2020", 
                 paste0(i, ".results.DD", j, ".rds")))
  }
}
