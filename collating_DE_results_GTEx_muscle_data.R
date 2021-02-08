library(here)
library(compcodeR)
library(ROCR)
library(caret)

# First create variables to store results
{
  names.edgeR <- c("edgeR.ql.tmm", "edgeR.lr.tmm", "edgeR.et.tmm", 
                   "edgeR.ql.rle", "edgeR.lr.rle", "edgeR.et.rle")
  names.DESeq2 <- c("DESeq2.noif.tmm", "DESeq2.if.tmm", 
                    "DESeq2.noif.rle", "DESeq2.if.rle")
  names.voom <- c("voom.tmm", "voom.rle")
  names.DSS <- c("DSS.tmm", "DSS.rle")
  names.baySeq <- c("baySeq.tmm", "baySeq.rle")
  names.MDSeq <- c("mean.MDSeq.zi.tmm", "mean.MDSeq.nozi.tmm", 
                   "mean.MDSeq.zi.rle", "mean.MDSeq.nozi.rle")
  names.expHMM <- c("expHMM.tmm", "expHMM.rle")
  names.expHM <- c("mean.expHM.untr.tmm", "mean.expHM.log.tmm", 
                   "mean.expHM.untr.rle", "mean.expHM.log.rle")
  names.lnHMM <- c("lnHMM.tmm", "lnHMM.rle")
  names.lnHM <- c("mean.lnHM.untr.tmm", "mean.lnHM.log.tmm", 
                  "mean.lnHM.untr.rle", "mean.lnHM.log.rle")
  results.list <- c(names.edgeR, names.DESeq2, names.voom, 
                    names.DSS,
                    names.baySeq, 
                    names.MDSeq, names.expHMM, names.expHM, names.lnHMM, names.lnHM)
}

for (samples_per_group in c(2, 5, 10, 20, 50)) {
  {
    for (i in c(names.edgeR, names.DESeq2, names.voom, names.DSS, names.baySeq, 
                names.MDSeq, names.expHM, names.lnHM)) {
      for (j in c("fpr.", "fdr.", "tpr.")) {
        assign(paste0(j, i), rep(NA, 10))
      }
    }
    for (i in names.DSS) {
      for (j in c("fpr.lfdr.", "fdr.lfdr.", "tpr.lfdr.")) {
        assign(paste0(j, i), rep(NA, 10))
      }
    }
    for (i in c(names.expHMM, names.lnHMM)) {
      for (j in c("fpr.point5.", "fdr.point5.", "tpr.point5.", 
                  "fpr.thr.", "fdr.thr.", "tpr.thr.", 
                  "fpr.bfdr.", "fdr.bfdr.", "tpr.bfdr.")) {
        assign(paste0(j, i), rep(NA, 10))
      }
    }
    for (i in results.list) {
      assign(paste0("auc.", i), rep(NA, 10))
      assign(paste0("pauc.", i), rep(NA, 10))
      assign(paste0("discoveries.", i), matrix(nrow=10, ncol=20000))
      assign(paste0("false.discoveries.", i), matrix(nrow=10, ncol=20000))
    }
  }
  
  # Then import each run in turn and save results
  folder <- "Results/GTEx muscle simulated DE, DD, DEDD results March 2020"
  for (i in 1:10) {
    filename <- paste0("results.muscle_", samples_per_group, "_set", i, "_DE.rds")
    res <- readRDS(here(folder, filename))
    DE <- factor(res$DE, levels=c("1", "0"))
    DEDD <- DE
    for (j in c(names.edgeR, names.DESeq2, names.voom, names.DSS, 
                names.MDSeq, names.expHM, names.lnHM)) {
      if (mean(is.na(get(paste0("p.",j), res))) < 1) { # added because no data for DSS for some runs
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
    for (j in names.DSS) {
      if (mean(is.na(get(paste0("p.",j), res))) < 1) {
        assign(paste0("call.lfdr.", j), factor(as.numeric(get(paste0("lfdr.", j), res) < 0.05), 
                                               levels=c("1", "0")))
        assign(paste0("fpr.lfdr.", j), `[<-`(get(paste0("fpr.lfdr.",j)), i,
                                             value=1 - specificity(get(paste0("call.lfdr.", j)), DE)))
        assign(paste0("fdr.lfdr.", j), `[<-`(get(paste0("fdr.lfdr.",j)), i,
                                             value=1 - precision(get(paste0("call.lfdr.", j)), DE)))
        assign(paste0("tpr.lfdr.", j), `[<-`(get(paste0("tpr.lfdr.",j)), i,
                                             value=sensitivity(get(paste0("call.lfdr.", j)), DE)))
      }
    }
    for (j in names.baySeq) {
      assign(paste0("call.", j), factor(as.numeric(get(paste0("q.", j), res) < 0.05), 
                                        levels=c("1", "0")))
      assign(paste0("pred.", j), prediction(get(paste0("prob.",j), res), DE))
      assign(paste0("fpr.", j), `[<-`(get(paste0("fpr.",j)), i,
                                      value=1 - specificity(get(paste0("call.", j)), DE)))
      assign(paste0("fdr.", j), `[<-`(get(paste0("fdr.",j)), i,
                                      value=1 - precision(get(paste0("call.", j)), DE)))
      assign(paste0("tpr.", j), `[<-`(get(paste0("tpr.",j)), i,
                                      value=sensitivity(get(paste0("call.", j)), DE)))
      assign(paste0("auc.", j), `[<-`(get(paste0("auc.",j)), i,
                                      value=performance(get(paste0("pred.", j)), 
                                                        measure="auc")@y.values[[1]]))
      assign(paste0("pauc.", j), `[<-`(get(paste0("pauc.",j)), i,
                                       value=performance(get(paste0("pred.", j)), measure="auc", 
                                                         fpr.stop=0.05)@y.values[[1]]))
      assign(paste0("false.discoveries.", j), `[<-`(get(paste0("false.discoveries.", j)), i,,
                                                    value=c(get(paste0("pred.",j))@fp[[1]],
                                                            rep(NA, 20000 - 
                                                                  length(get(paste0("pred.", j))@fp[[1]])))))
      assign(paste0("discoveries.", j), `[<-`(get(paste0("discoveries.",j)), i,,
                                              value=c(get(paste0("pred.",j))@n.pos.pred[[1]],
                                                      rep(NA, 20000 - 
                                                            length(get(paste0("pred.",j))@n.pos.pred[[1]])))))
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
    rm(list=c("i", "j", "filename", "res", "DE", "DEDD"))
  }
  
  # Re-format results to save by experiment type ####
  
  # Set up variables
  {
    names.DE <- c("edgeR.ql.tmm", "edgeR.lr.tmm", "edgeR.et.tmm", "DESeq2.noif.tmm", "DESeq2.if.tmm", 
                  "voom.tmm", 
                  "DSS.tmm", 
                  "baySeq.tmm", "mean.MDSeq.zi.tmm", "mean.MDSeq.nozi.tmm", 
                  "mean.expHM.untr.tmm", "mean.expHM.log.tmm", "mean.lnHM.untr.tmm", "mean.lnHM.log.tmm", 
                  "edgeR.ql.rle", "edgeR.lr.rle", "edgeR.et.rle", "DESeq2.noif.rle", "DESeq2.if.rle", 
                  "voom.rle", "DSS.rle", "baySeq.rle", "mean.MDSeq.zi.rle", "mean.MDSeq.nozi.rle", 
                  "mean.expHM.untr.rle", "mean.expHM.log.rle", "mean.lnHM.untr.rle", "mean.lnHM.log.rle")
    names.DEDD <- c("expHMM.tmm", "expHMM.rle", "lnHMM.tmm", "lnHMM.rle")
  }
  
  for (i in c("DE", "DEDD")) {
    assign(paste0("results.", i), list())
    for (j in c("pauc", "auc", "fpr", "fdr", "tpr")) {
      assign(paste0("results.", i), `[[<-`(get(paste0("results.", i)), j, 
                                           data.frame(matrix(nrow=10, ncol=0))))
    }
    for (j in c("mean.fdr", "mean.discoveries")) {
      assign(paste0("results.", i), `[[<-`(get(paste0("results.", i)), j, 
                                           data.frame(matrix(nrow=20000, ncol=0))))
    }
  }
  
  {
    names.fpr.fdr.tpr <- vector()
    for (i in c(names.edgeR, names.DESeq2, names.voom, names.DSS, names.baySeq, 
                names.MDSeq, names.expHM, names.lnHM)) {
      names.fpr.fdr.tpr <- c(names.fpr.fdr.tpr, i)
    }
    for (i in names.DSS) {
      names.fpr.fdr.tpr <- c(names.fpr.fdr.tpr, paste0("lfdr.", i))
    }
    for (i in c(names.expHMM, names.lnHMM)) {
      for (j in c("point5.", "thr.", "bfdr.")) {
        names.fpr.fdr.tpr <- c(names.fpr.fdr.tpr, paste0(j, i))
      }
    }
  }
  
  # Put results into variables
  {
    for (i in names.DE) {
      results.DE$auc[[i]] <- get(paste0("auc.", i))
      results.DE$pauc[[i]] <- get(paste0("pauc.", i))
      results.DE$mean.fdr[[i]] <- get(paste0("mean.fdr.", i))
      results.DE$mean.discoveries[[i]] <- get(paste0("mean.discoveries.", i))
      for (j in grep(i, names.fpr.fdr.tpr, value=T)) {
        results.DE$fpr[[j]] <- get(paste0("fpr.", j))
        results.DE$fdr[[j]] <- get(paste0("fdr.", j))
        results.DE$tpr[[j]] <- get(paste0("tpr.", j))
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
  
  saveRDS(results.DE, file=here(folder, paste0("DE.results.muscle_", samples_per_group, "_DE.rds")))
  saveRDS(results.DEDD, file=here(folder, paste0("DEDD.results.muscle_", samples_per_group, "_DE.rds")))
  
}
