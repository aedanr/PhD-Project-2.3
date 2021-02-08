library(here)

folder <- "Results/compcodeR DE, DD, DEDD results Feb 2020"
for (i in c("DEDD2", "DEDD5", "DEDD10", "DEDD20", "DEDD50")) {
  assign(paste0("DEDD.results.compcodeR.", i), 
         readRDS(here(folder, paste0("DEDD.results.", i, ".rds"))))
}
folder <- "Results/diffVar compcodeR results Dec 2020"
for (i in c("2", "5", "10", "20", "50")) {
  assign(paste0("DEDD.results.diffVar.compcodeR.DD", i), 
         readRDS(here(folder, paste0("DEDD.results.DD", i, ".rds"))))
  assign(paste0("DEDD.results.diffVar.compcodeR.DEDD", i), 
         readRDS(here(folder, paste0("DEDD.results.DEDD", i, ".rds"))))
}
rm(i,folder)

for (i in c("2", "5", "10", "20", "50")) {
  target <- get(paste0("DEDD.results.compcodeR.DEDD", i))
  source <- get(paste0("DEDD.results.diffVar.compcodeR.DEDD", i))
  for (j in c("fpr", "fdr", "tpr")) {
    target[[j]] <- cbind(target[[j]], source[[j]])
    names(target[[j]])[21] <- c("dV.tmm")
    target[[j]] <- target[[j]][, c(1:12, 21, 13:20)]
  }
  for (j in c("auc", "pauc", "mean.fdr", "mean.discoveries")) {
    target[[j]] <- cbind(target[[j]], source[[j]])
    names(target[[j]])[13] <- c("dV.tmm")
    target[[j]] <- target[[j]][, c(1:4, 13, 5:12)]
  }
  saveRDS(target, here("Results/compcodeR combined results Dec 2020 including diffVar",
                       paste0("DEDD.results.DEDD", i, ".rds")))
}
