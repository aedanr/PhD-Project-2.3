library(here)

folder <- "Results/compcodeR combined results Dec 2020"
for (i in c("DEDD2", "DEDD5", "DEDD10", "DEDD20", "DEDD50")) {
  assign(paste0("DD.results.compcodeR.", i), 
         readRDS(here(folder, paste0("DD.results.", i, ".rds"))))
}
for (i in c("DD2", "DD5", "DD10", "DD20", "DD50")) {
  assign(paste0("DD.results.compcodeR.", i), 
         readRDS(here(folder, paste0("DD.results.", i, ".rds"))))
}
folder <- "Results/diffVar compcodeR results Dec 2020"
for (i in c("2", "5", "10", "20", "50")) {
  assign(paste0("DD.results.diffVar.compcodeR.DD", i), 
         readRDS(here(folder, paste0("DD.results.DD", i, ".rds"))))
  assign(paste0("DD.results.diffVar.compcodeR.DEDD", i), 
         readRDS(here(folder, paste0("DD.results.DEDD", i, ".rds"))))
}
rm(i,folder)

for (i in c("2", "5", "10", "20", "50")) {
  for (k in c("DD", "DEDD")) {
    target <- get(paste0("DD.results.compcodeR.", k, i))
    source <- get(paste0("DD.results.diffVar.compcodeR.", k, i))
    for (j in c("auc", "pauc", "fpr", "fdr", "tpr", "mean.fdr", "mean.discoveries")) {
      target[[j]] <- cbind(target[[j]], source[[j]])
      names(target[[j]])[14] <- c("dV.tmm")
      target[[j]] <- target[[j]][, c(1:3, 14, 4:13)]
    }
    saveRDS(target, here("Results/compcodeR combined results Dec 2020 including diffVar",
                         paste0("DD.results.", k, i, ".rds")))
  }
}
