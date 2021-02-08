library(here)
library(ROCR)
library(caret)

# First create variables to store results
{
  names.MDSeq <- c("disp.MDSeq.zi", "disp.lfc1.MDSeq.zi")
  names.lnHM <- c("disp.lnHM", "disp.lfc1.lnHM")
  results.list <- c(names.MDSeq, names.lnHM)
}

for (samples_per_group in c(2,5,10,20,50)) {
  
  for (i in results.list) {
    for (j in c("auc.", "pauc.", "fpr.", "fdr.", "tpr.")) {
      assign(paste0(j, i), numeric(10))
    }
    assign(paste0("discoveries.", i), matrix(nrow=10, ncol=20000))
    assign(paste0("false.discoveries.", i), matrix(nrow=10, ncol=20000))
  }
  
  # Then import each run in turn and save results
  folder <- "Results/GTEx DD lfc1 results Jan 2021"
  for (i in 1:10) {
    filename <- paste0("results.muscle_", samples_per_group, "_set", i, "_DD.rds")
    res <- readRDS(here(folder, filename))
    DD <- factor(res$DD, levels=c("1", "0"))
    lfcd1 <- factor(as.numeric(res$lfcd1), levels=c("1", "0"))
    for (j in c(results.list)) {
      assign(paste0("call.", j), factor(as.numeric(get(paste0("q.", j), res) < 0.05), 
                                        levels=c("1", "0")))
      if (grepl("lfc1", j)) {
        assign(paste0("pred.", j), prediction(1 - get(paste0("p.",j), res)[which(!is.na(get(paste0("p.", j), res)))], 
                                              lfcd1[which(!is.na(get(paste0("p.", j), res)))]))
        assign(paste0("fpr.", j), `[<-`(get(paste0("fpr.", j)), i,
                                        value=1 - specificity(get(paste0("call.", j)), lfcd1)))
        assign(paste0("fdr.", j), `[<-`(get(paste0("fdr.", j)), i,
                                        value=1 - precision(get(paste0("call.", j)), lfcd1)))
        assign(paste0("tpr.", j), `[<-`(get(paste0("tpr.", j)), i,
                                        value=sensitivity(get(paste0("call.", j)), lfcd1)))
      }
      else {
        assign(paste0("pred.", j), prediction(1 - get(paste0("p.", j), res)[which(!is.na(get(paste0("p.", j), res)))], 
                                              DD[which(!is.na(get(paste0("p.", j), res)))]))
        assign(paste0("fpr.", j), `[<-`(get(paste0("fpr.", j)), i,
                                        value=1 - specificity(get(paste0("call.", j)), DD)))
        assign(paste0("fdr.", j), `[<-`(get(paste0("fdr.", j)), i,
                                        value=1 - precision(get(paste0("call.", j)), DD)))
        assign(paste0("tpr.", j), `[<-`(get(paste0("tpr.", j)), i,
                                        value=sensitivity(get(paste0("call.", j)), DD)))
      }
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
    rm(list=c("i", "j", "filename", "res", "DD", "lfcd1"))
  }
  
  # Re-format results to save by experiment type ####
  
  # Set up variables
  {
    names.DD <- c("disp.MDSeq.zi", "disp.lnHM")
    names.DD.lfc1 <- c("disp.lfc1.MDSeq.zi", "disp.lfc1.lnHM")
  }
  
  for (i in c("DD", "DD.lfc1")) {
    assign(paste0("results.", i), list())
    for (j in c("pauc", "auc", "fpr", "fdr", "tpr")) {
      assign(paste0("results.", i), `[[<-`(get(paste0("results.", i)), j, data.frame(matrix(nrow=10, ncol=0))))
    }
    for (j in c("mean.fdr", "mean.discoveries")) {
      assign(paste0("results.", i), `[[<-`(get(paste0("results.", i)), j, data.frame(matrix(nrow=20000, ncol=0))))
    }
  }
  
  {
    names.fpr.fdr.tpr <- vector()
    for (i in results.list) {
      names.fpr.fdr.tpr <- c(names.fpr.fdr.tpr, i)
    }
  }
  
  # Put results into variables
  {
    for (i in names.DD) {
      results.DD$auc[[i]] <- get(paste0("auc.", i))
      results.DD$pauc[[i]] <- get(paste0("pauc.", i))
      results.DD$fpr[[i]] <- get(paste0("fpr.", i))
      results.DD$fdr[[i]] <- get(paste0("fdr.", i))
      results.DD$tpr[[i]] <- get(paste0("tpr.", i))
      results.DD$mean.fdr[[i]] <- get(paste0("mean.fdr.", i))
      results.DD$mean.discoveries[[i]] <- get(paste0("mean.discoveries.", i))
    }
    for (i in names.DD.lfc1) {
      results.DD.lfc1$auc[[i]] <- get(paste0("auc.", i))
      results.DD.lfc1$pauc[[i]] <- get(paste0("pauc.", i))
      results.DD.lfc1$fpr[[i]] <- get(paste0("fpr.", i))
      results.DD.lfc1$fdr[[i]] <- get(paste0("fdr.", i))
      results.DD.lfc1$tpr[[i]] <- get(paste0("tpr.", i))
      results.DD.lfc1$mean.fdr[[i]] <- get(paste0("mean.fdr.", i))
      results.DD.lfc1$mean.discoveries[[i]] <- get(paste0("mean.discoveries.", i))
    }
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
  
  saveRDS(results.DD, file=here(folder, paste0("DD.results.muscle_", samples_per_group, "_DD.rds")))
  saveRDS(results.DD.lfc1, file=here(folder, paste0("DD.lfc1.results.muscle_", samples_per_group, "_DD.rds")))
  
}
