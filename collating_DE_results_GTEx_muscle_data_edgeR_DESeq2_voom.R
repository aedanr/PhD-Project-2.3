library(here)
library(compcodeR)
library(ROCR)
library(caret)

# First create variables to store results
results.list <- c("edgeR.ql.tmm", "edgeR.lr.tmm", "edgeR.et.tmm", 
                  "DESeq2.noif.tmm", "DESeq2.if.tmm", "voom.tmm")

for (samples_per_group in c(2, 5, 10, 20, 50)) {
  for (i in results.list) {
    for (j in c("fpr.", "fdr.", "tpr.")) {
      assign(paste0(j, i), numeric(10))
    }
    assign(paste0("auc.", i), numeric(10))
    assign(paste0("pauc.", i), numeric(10))
    assign(paste0("discoveries.", i), matrix(nrow=10, ncol=20000))
    assign(paste0("false.discoveries.", i), matrix(nrow=10, ncol=20000))
  }
  
  # Then import each run in turn and save results
  folder <- "Results/GTEx muscle simulated DE, DD, DEDD results Feb 2020/edgeR, DESeq2, voom only"
  for (i in 1:10) {
    filename <- paste0("results.muscle_", samples_per_group, "_set", i, "_DE.rds")
    res <- readRDS(here(folder, filename))
    DE <- factor(res$DE, levels=c("1", "0"))
    for (j in c(results.list)) {
      assign(paste0("call.", j), factor(as.numeric(get(paste0("q.", j), res) < 0.05), 
                                        levels=c("1", "0")))
      assign(paste0("pred.", j), prediction(1 - get(paste0("p.",j), res), DE))
      assign(paste0("fpr.", j), `[<-`(get(paste0("fpr.", j)), i,
                                      value=1 - specificity(get(paste0("call.", j)), DE)))
      assign(paste0("fdr.", j), `[<-`(get(paste0("fdr.", j)), i,
                                      value=1 - precision(get(paste0("call.", j)), DE)))
      assign(paste0("tpr.", j), `[<-`(get(paste0("tpr.", j)), i,
                                      value=sensitivity(get(paste0("call.", j)), DE)))
      assign(paste0("auc.", j), `[<-`(get(paste0("auc.", j)), i, 
                                      value=performance(get(paste0("pred.", j)), 
                                                        measure="auc")@y.values[[1]]))
      assign(paste0("pauc.", j), `[<-`(get(paste0("pauc.", j)), i, 
                                       value=performance(get(paste0("pred.", j)), measure="auc", 
                                                         fpr.stop=0.05)@y.values[[1]]))
      assign(paste0("false.discoveries.", j), `[<-`(get(paste0("false.discoveries.", j)), i,, 
                                                    value=c(get(paste0("pred.", j))@fp[[1]], 
                                                            rep(NA, 20000 - 
                                                                  length(get(paste0("pred.", j))@fp[[1]])))))
      assign(paste0("discoveries.", j), `[<-`(get(paste0("discoveries.", j)), i,, 
                                              value=c(get(paste0("pred.", j))@n.pos.pred[[1]], 
                                                      rep(NA, 20000 - 
                                                            length(get(paste0("pred.", j))@n.pos.pred[[1]])))))
    }
  }
  
  for (i in results.list) {
    assign(paste0("mean.discoveries.", i), colMeans(get(paste0("discoveries.", i))))
    assign(paste0("mean.fdr.", i), colMeans(get(paste0("false.discoveries.", i)) 
                                            / get(paste0("discoveries.", i))))
  }
  
  # Clean up temporary variables no longer needed
  {
    rm(list=ls()[grep("^pred", ls())])
    rm(list=ls()[grep("^q", ls())])
    rm(list=ls()[grep("^thr", ls())])
    rm(list=ls()[grep("^call", ls())])
    rm(list=ls()[grep("^discoveries", ls())])
    rm(list=ls()[grep("^false.discoveries", ls())])
    rm(list=c("i", "j", "filename", "res", "DE"))
  }
  
  # Re-format results to save by experiment type ####
  
  # Set up variables
  names.DE <- c("edgeR.ql.tmm", "edgeR.lr.tmm", "edgeR.et.tmm", 
                "DESeq2.noif.tmm", "DESeq2.if.tmm", "voom.tmm")
  
  results.DE <- list()
  for (j in c("pauc", "auc", "fpr", "fdr", "tpr")) {
    results.DE[[j]] <- data.frame(matrix(nrow=10, ncol=0))
  }
  for (j in c("mean.fdr", "mean.discoveries")) {
    results.DE[[j]] <- data.frame(matrix(nrow=20000, ncol=0))
  }
  
  # Put results into variables
  for (i in names.DE) {
    results.DE$auc[[i]] <- get(paste0("auc.", i))
    results.DE$pauc[[i]] <- get(paste0("pauc.", i))
    results.DE$fpr[[i]] <- get(paste0("fpr.", i))
    results.DE$fdr[[i]] <- get(paste0("fdr.", i))
    results.DE$tpr[[i]] <- get(paste0("tpr.", i))
    results.DE$mean.fdr[[i]] <- get(paste0("mean.fdr.", i))
    results.DE$mean.discoveries[[i]] <- get(paste0("mean.discoveries.", i))
  }
  
  # Clean up temporary variables no longer needed and save results
  {
    rm(list=ls()[grep("^auc", ls())])
    rm(list=ls()[grep("^pauc", ls())])
    rm(list=ls()[grep("^fdr", ls())])
    rm(list=ls()[grep("^fpr", ls())])
    rm(list=ls()[grep("^mean", ls())])
    rm(list=ls()[grep("^tpr", ls())])
    rm("i", "j")
  }
  
  saveRDS(results.DE, file=here(folder, paste0("DE.results.muscle_", samples_per_group, "_DE.rds")))

}
