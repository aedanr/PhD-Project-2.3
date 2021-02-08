library(here)
library(compcodeR)
library(ROCR)
library(caret)

# First create variables to store results
{
  names.MDSeq <- c("disp.MDSeq.zi.tmm", "disp.lfc1.MDSeq.zi.tmm", "disp.lfc2.MDSeq.zi.tmm", 
                   "disp.MDSeq.nozi.tmm", "disp.lfc1.MDSeq.nozi.tmm", "disp.lfc2.MDSeq.nozi.tmm", 
                   "disp.MDSeq.zi.rle", "disp.lfc1.MDSeq.zi.rle", "disp.lfc2.MDSeq.zi.rle", 
                   "disp.MDSeq.nozi.rle", "disp.lfc1.MDSeq.nozi.rle", "disp.lfc2.MDSeq.nozi.rle")
  names.expHMM <- c("expHMM.tmm", "expHMM.rle")
  names.expHM <- c("disp.expHM.untr.tmm", "disp.expHM.log.tmm", 
                   "disp.lfc1.expHM.tmm", "disp.lfc2.expHM.tmm", 
                   "disp.expHM.untr.rle", "disp.expHM.log.rle", 
                   "disp.lfc1.expHM.rle", "disp.lfc2.expHM.rle")
  names.lnHMM <- c("lnHMM.tmm", "lnHMM.rle")
  names.lnHM <- c("disp.lnHM.untr.tmm", "disp.lnHM.log.tmm", 
                  "disp.lfc1.lnHM.tmm", "disp.lfc2.lnHM.tmm", 
                  "disp.lnHM.untr.rle", "disp.lnHM.log.rle", 
                  "disp.lfc1.lnHM.rle", "disp.lfc2.lnHM.rle")
  results.list <- c(names.MDSeq, names.expHMM, names.expHM, names.lnHMM, names.lnHM)
}


for (samples_per_group in c(2,5,10,20,50)) {
  {
    for (i in c(names.MDSeq, names.expHM, names.lnHM)) {
      for (j in c("fpr.", "fdr.", "tpr.")) {
        assign(paste0(j, i), numeric(50))
      }
    }
    for (i in c(names.expHMM, names.lnHMM)) {
      for (j in c("fpr.point5.", "fdr.point5.", "tpr.point5.", 
                  "fpr.thr.", "fdr.thr.", "tpr.thr.", 
                  "fpr.bfdr.", "fdr.bfdr.", "tpr.bfdr.")) {
        assign(paste0(j, i), numeric(50))
      }
    }
    for (i in results.list) {
      assign(paste0("auc.", i), numeric(50))
      assign(paste0("pauc.", i), numeric(50))
      assign(paste0("discoveries.", i), matrix(nrow=50, ncol=20000))
      assign(paste0("false.discoveries.", i), matrix(nrow=50, ncol=20000))
    }
  }
  
  # Then import each run in turn and save results
  # folder <- "Results/compcodeR DE, DD, DEDD results Feb 2020"
  folder <- "Results/temp_quarantine"
  for (i in 1:50) {
    filename <- paste0("results.DD", samples_per_group, ".", i, ".rds")
    res <- readRDS(here(folder, filename))
    DD <- factor(res$DD, levels=c("1", "0"))
    DEDD <- DD
    lfcd1 <- factor(as.numeric(res$lfcd1), levels=c("1", "0"))
    lfcd2 <- factor(as.numeric(res$lfcd2), levels=c("1", "0"))
    for (j in c(names.MDSeq, names.expHM, names.lnHM)) {
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
      else if (grepl("lfc2", j)) {
        assign(paste0("pred.", j), prediction(1 - get(paste0("p.",j), res)[which(!is.na(get(paste0("p.", j), res)))], 
                                              lfcd2[which(!is.na(get(paste0("p.", j), res)))]))
        assign(paste0("fpr.", j), `[<-`(get(paste0("fpr.",j)), i,
                                        value=1 - specificity(get(paste0("call.",j)), lfcd2)))
        assign(paste0("fdr.", j), `[<-`(get(paste0("fdr.",j)), i,
                                        value=1 - precision(get(paste0("call.",j)), lfcd2)))
        assign(paste0("tpr.", j), `[<-`(get(paste0("tpr.",j)), i,
                                        value=sensitivity(get(paste0("call.",j)), lfcd2)))
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
    for (j in c(names.expHMM, names.lnHMM)) {
      assign(paste0("call.point5.", j), factor(as.numeric(get(paste0("prob.", j), res) > 0.5), 
                                               levels=c("1", "0")))
      assign(paste0("call.thr.", j), 
             factor(as.numeric(get(paste0("prob.", j), res) > get(paste0("thr.", j), res)), levels=c("1", "0")))
      assign(paste0("call.bfdr.", j), factor(as.numeric(get(paste0("bfdr.", j), res) < 0.05), 
                                             levels=c("1", "0")))
      assign(paste0("pred.", j), prediction(get(paste0("prob.", j), res), DEDD))
      assign(paste0("fpr.point5.", j), `[<-`(get(paste0("fpr.point5.", j)), i, 
                                             value=1-specificity(get(paste0("call.point5.",j)), DEDD)))
      assign(paste0("fpr.thr.", j), `[<-`(get(paste0("fpr.thr.", j)), i, 
                                          value=1-specificity(get(paste0("call.thr.",j)), DEDD)))
      assign(paste0("fpr.bfdr.", j), `[<-`(get(paste0("fpr.bfdr.", j)), i, 
                                           value=1-specificity(get(paste0("call.bfdr.",j)), DEDD)))
      assign(paste0("fdr.point5.", j), `[<-`(get(paste0("fdr.point5.", j)), i, 
                                             value=1-precision(get(paste0("call.point5.",j)), DEDD)))
      assign(paste0("fdr.thr.", j), `[<-`(get(paste0("fdr.thr.", j)), i, 
                                          value=1-precision(get(paste0("call.thr.",j)), DEDD)))
      assign(paste0("fdr.bfdr.", j), `[<-`(get(paste0("fdr.bfdr.", j)), i, 
                                           value=1-precision(get(paste0("call.bfdr.",j)), DEDD)))
      assign(paste0("tpr.point5.", j), `[<-`(get(paste0("tpr.point5.", j)), i, 
                                             value=sensitivity(get(paste0("call.point5.",j)), DEDD)))
      assign(paste0("tpr.thr.", j), `[<-`(get(paste0("tpr.thr.", j)), i, 
                                          value=sensitivity(get(paste0("call.thr.",j)), DEDD)))
      assign(paste0("tpr.bfdr.", j), `[<-`(get(paste0("tpr.bfdr.", j)), i, 
                                           value=sensitivity(get(paste0("call.bfdr.",j)), DEDD)))
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
                                                            length(get(paste0("pred." ,j))@n.pos.pred[[1]])))))
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
    rm(list=c("i", "j", "filename", "res", "DD", "DEDD", "lfcd1", "lfcd2"))
  }
  
  # Re-format results to save by experiment type ####
  
  # Set up variables
  {
    names.DD <- c("disp.MDSeq.zi.tmm", "disp.MDSeq.nozi.tmm", 
                  "disp.expHM.untr.tmm", "disp.expHM.log.tmm", 
                  "disp.lnHM.untr.tmm", "disp.lnHM.log.tmm", 
                  "disp.MDSeq.zi.rle", "disp.MDSeq.nozi.rle", 
                  "disp.expHM.untr.rle", "disp.expHM.log.rle", 
                  "disp.lnHM.untr.rle", "disp.lnHM.log.rle")
    names.DD.lfc1 <- c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.MDSeq.nozi.tmm", 
                       "disp.lfc1.expHM.tmm", "disp.lfc1.lnHM.tmm", 
                       "disp.lfc1.MDSeq.zi.rle", "disp.lfc1.MDSeq.nozi.rle", 
                       "disp.lfc1.expHM.rle", "disp.lfc1.lnHM.rle")
    names.DD.lfc2 <- c("disp.lfc2.MDSeq.zi.tmm", "disp.lfc2.MDSeq.nozi.tmm", 
                       "disp.lfc2.expHM.tmm", "disp.lfc2.lnHM.tmm", 
                       "disp.lfc2.MDSeq.zi.rle", "disp.lfc2.MDSeq.nozi.rle", 
                       "disp.lfc2.expHM.rle", "disp.lfc2.lnHM.rle")
    names.DEDD <- c("expHMM.tmm", "expHMM.rle", "lnHMM.tmm", "lnHMM.rle")
  }
  
  for (i in c("DD", "DD.lfc1", "DD.lfc2", "DEDD")) {
    assign(paste0("results.", i), list())
    for (j in c("pauc", "auc", "fpr", "fdr", "tpr")) {
      assign(paste0("results.", i), `[[<-`(get(paste0("results.", i)), j, data.frame(matrix(nrow=50, ncol=0))))
    }
    for (j in c("mean.fdr", "mean.discoveries")) {
      assign(paste0("results.", i), `[[<-`(get(paste0("results.", i)), j, data.frame(matrix(nrow=20000, ncol=0))))
    }
  }
  
  {
    names.fpr.fdr.tpr <- vector()
    for (i in c(names.MDSeq, names.expHM, names.lnHM)) {
      names.fpr.fdr.tpr <- c(names.fpr.fdr.tpr, i)
    }
    for (i in c(names.expHMM, names.lnHMM)) {
      for (j in c("point5.", "thr.", "bfdr.")) {
        names.fpr.fdr.tpr <- c(names.fpr.fdr.tpr, paste0(j, i))
      }
    }
  }
  
  # Put results into variables
  {
    for (i in names.DD) {
      results.DD$auc[[i]] <- get(paste0("auc.", i))
      results.DD$pauc[[i]] <- get(paste0("pauc.", i))
      results.DD$mean.fdr[[i]] <- get(paste0("mean.fdr.", i))
      results.DD$mean.discoveries[[i]] <- get(paste0("mean.discoveries.", i))
      for (j in grep(i, names.fpr.fdr.tpr, value=T)) {
        results.DD$fpr[[j]] <- get(paste0("fpr.", j))
        results.DD$fdr[[j]] <- get(paste0("fdr.", j))
        results.DD$tpr[[j]] <- get(paste0("tpr.", j))
      }
    }
    for (i in names.DD.lfc1) {
      results.DD.lfc1$auc[[i]] <- get(paste0("auc.", i))
      results.DD.lfc1$pauc[[i]] <- get(paste0("pauc.", i))
      results.DD.lfc1$mean.fdr[[i]] <- get(paste0("mean.fdr.", i))
      results.DD.lfc1$mean.discoveries[[i]] <- get(paste0("mean.discoveries.", i))
      for (j in grep(i, names.fpr.fdr.tpr, value=T)) {
        results.DD.lfc1$fpr[[j]] <- get(paste0("fpr.", j))
        results.DD.lfc1$fdr[[j]] <- get(paste0("fdr.", j))
        results.DD.lfc1$tpr[[j]] <- get(paste0("tpr.", j))
      }
    }
    for (i in names.DD.lfc2) {
      results.DD.lfc2$auc[[i]] <- get(paste0("auc.", i))
      results.DD.lfc2$pauc[[i]] <- get(paste0("pauc.", i))
      results.DD.lfc2$mean.fdr[[i]] <- get(paste0("mean.fdr.", i))
      results.DD.lfc2$mean.discoveries[[i]] <- get(paste0("mean.discoveries.", i))
      for (j in grep(i, names.fpr.fdr.tpr, value=T)) {
        results.DD.lfc2$fpr[[j]] <- get(paste0("fpr.", j))
        results.DD.lfc2$fdr[[j]] <- get(paste0("fdr.", j))
        results.DD.lfc2$tpr[[j]] <- get(paste0("tpr.", j))
      }
    }
    for (i in names.DEDD) {
      results.DEDD$auc[[i]] <- get(paste0("auc.", i))
      results.DEDD$pauc[[i]] <- get(paste0("pauc.", i))
      results.DEDD$mean.fdr[[i]] <- get(paste0("mean.fdr.", i))
      results.DEDD$mean.discoveries[[i]] <- get(paste0("mean.discoveries.", i))
      for (j in grep(i, names.fpr.fdr.tpr, value=T)) {
        results.DEDD$fpr[[j]] <- get(paste0("fpr.", j))
        results.DEDD$fdr[[j]] <- get(paste0("fdr.", j))
        results.DEDD$tpr[[j]] <- get(paste0("tpr.", j))
      }
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
  
  saveRDS(results.DD, file=here(folder, paste0("DD.results.DD", samples_per_group, ".rds")))
  saveRDS(results.DD.lfc1, file=here(folder, paste0("DD.lfc1.results.DD", samples_per_group, ".rds")))
  saveRDS(results.DD.lfc2, file=here(folder, paste0("DD.lfc2.results.DD", samples_per_group, ".rds")))
  saveRDS(results.DEDD, file=here(folder, paste0("DEDD.results.DD", samples_per_group, ".rds")))
  
}
