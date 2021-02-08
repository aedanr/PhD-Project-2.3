library(here)

## GTEx, DEDD data only since don't have combined results for DD data (yet)
folder <- "Results/GTEx combined results Sept 2020"
for (i in c("2", "5", "10", "20", "50")) {
  assign(paste0("DD.results.GTEx.DD", i), 
         readRDS(here(folder, paste0("DD.results.DD", i, ".rds"))))
}

folder <- "Results/ExpVarQuant GTEx results Dec 2020"
for (i in c("2", "5", "10", "20", "50")) {
  assign(paste0("DD.results.ExpVarQuant.GTEx.blood.DD", i), 
         readRDS(here(folder, paste0("DD.results.blood_", i, "_DD.rds"))))
  assign(paste0("DD.results.ExpVarQuant.GTEx.muscle.DD", i), 
         readRDS(here(folder, paste0("DD.results.muscle_", i, "_DD.rds"))))
}
rm(i,folder)

for (i in c("2", "5", "10", "20", "50")) {
  target <- get(paste0("DD.results.GTEx.DD", i))
  blood <- get(paste0("DD.results.ExpVarQuant.GTEx.blood.DD", i))
  muscle <- get(paste0("DD.results.ExpVarQuant.GTEx.muscle.DD", i))
  for (j in c("auc", "pauc", "fpr", "fdr", "tpr", "mean.fdr", "mean.discoveries")) {
    target[[j]] <- cbind(target[[j]], blood[[j]], muscle[[j]])
    names(target[[j]])[15:16] <- c("blood_EVQ", "muscle_EVQ")
    target[[j]] <- target[[j]][, c(1:3, 15, 4:10, 16, 11:14)]
  }
  saveRDS(target, here("Results/GTEx combined results Dec 2020", 
                       paste0("DD.results.DD", i, ".rds")))
}

