library(here)

## GTEx, DEDD data only since don't have combined results for DD data (yet)
folder <- "Results/GTEx combined results Sept 2020"
for (i in c("2", "5", "10", "20", "50")) {
  assign(paste0("DD.results.GTEx.DEDD", i), 
         readRDS(here(folder, paste0("DD.results.DEDD", i, ".rds"))))
}

folder <- "Results/ExpVarQuant GTEx results Dec 2020"
for (i in c("2", "5", "10", "20", "50")) {
  assign(paste0("DD.results.ExpVarQuant.GTEx.blood.DEDD", i), 
         readRDS(here(folder, paste0("DD.results.blood_", i, "_DEDD.rds"))))
  assign(paste0("DD.results.ExpVarQuant.GTEx.muscle.DEDD", i), 
         readRDS(here(folder, paste0("DD.results.muscle_", i, "_DEDD.rds"))))
}
rm(i,folder)

for (i in c("2", "5", "10", "20", "50")) {
  target <- get(paste0("DD.results.GTEx.DEDD", i))
  blood <- get(paste0("DD.results.ExpVarQuant.GTEx.blood.DEDD", i))
  muscle <- get(paste0("DD.results.ExpVarQuant.GTEx.muscle.DEDD", i))
  for (j in c("auc", "pauc", "fpr", "fdr", "tpr", "mean.fdr", "mean.discoveries")) {
    target[[j]] <- cbind(target[[j]], blood[[j]], muscle[[j]])
    names(target[[j]])[17:18] <- c("blood_EVQ", "muscle_EVQ")
    target[[j]] <- target[[j]][, c(1:3, 17, 4:10, 18, 11:16)]
  }
  saveRDS(target, here("Results/GTEx combined results Dec 2020", 
                       paste0("DD.results.DEDD", i, ".rds")))
}


## compcodeR, DD and DEDD data
folder <- "Results/compcodeR DE, DD, DEDD results Feb 2020"
for (i in c("DEDD2", "DEDD5", "DEDD10", "DEDD20", "DEDD50")) {
  assign(paste0("DD.results.compcodeR.", i), 
         readRDS(here(folder, paste0("DD.results.", i, ".rds"))))
}
for (i in c("DD2", "DD5", "DD10", "DD20", "DD50")) {
  assign(paste0("DD.results.compcodeR.", i), 
         readRDS(here(folder, paste0("DD.results.", i, ".rds"))))
}
folder <- "Results/ExpVarQuant compcodeR results Dec 2020"
for (i in c("2", "5", "10", "20", "50")) {
  assign(paste0("DD.results.ExpVarQuant.compcodeR.DD", i), 
         readRDS(here(folder, paste0("DD.results.DD", i, ".rds"))))
  assign(paste0("DD.results.ExpVarQuant.compcodeR.DEDD", i), 
         readRDS(here(folder, paste0("DD.results.DEDD", i, ".rds"))))
}
rm(i,folder)

for (i in c("2", "5", "10", "20", "50")) {
  for (k in c("DD", "DEDD")) {
    target <- get(paste0("DD.results.compcodeR.", k, i))
    source <- get(paste0("DD.results.ExpVarQuant.compcodeR.", k, i))
    for (j in c("auc", "pauc", "fpr", "fdr", "tpr", "mean.fdr", "mean.discoveries")) {
      target[[j]] <- cbind(target[[j]], source[[j]])
      names(target[[j]])[13] <- c("EVQ.tmm")
      target[[j]] <- target[[j]][, c(1:2, 13, 3:12)]
    }
    saveRDS(target, here("Results/compcodeR combined results Dec 2020",
                         paste0("DD.results.", k, i, ".rds")))
  }
}
