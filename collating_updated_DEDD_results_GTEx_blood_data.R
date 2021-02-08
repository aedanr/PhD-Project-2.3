# A few changes from previous results. These should be the final ones, and used for comparisons with 
# other results, e.g. long chain.
# Using log-transformed HM results to make hybrid DEDD predictions instead of untransformed. Not 
# sure why I used untransformed in the first place, but I've since decided it's better to consistently 
# use log-transformed, at least partly because they allow for LFC-based comparisons, and there doesn't 
# seem to be much difference between results with and without transformation.
# Using voom instead of edgeR for hybrid models since voom is better on GTEx data.
# Only TMM, no RLE.

library(here)
library(compcodeR)
library(ROCR)
library(caret)

# First create variables to store results
{
  names.edgeR <- c("edgeR.ql", "edgeR.lr", "edgeR.et")
  names.DESeq2 <- c("DESeq2.noif", "DESeq2.if")
  names.voom <- c("voom")
  names.DSS <- c("DSS")
  names.baySeq <- c("baySeq")
  names.diffVar <- c("diffVar")
  names.MDSeq <- c("mean.MDSeq.zi", "mean.MDSeq.nozi", 
                   "disp.MDSeq.zi", "disp.MDSeq.nozi")
  names.expHMM <- c("expHMM")
  names.expHM <- c("mean.expHM.untr", "mean.expHM.log", 
                   "disp.expHM.untr", "disp.expHM.log")
  names.lnHMM <- c("lnHMM")
  names.lnHM <- c("mean.lnHM.untr", "mean.lnHM.log", 
                  "disp.lnHM.untr", "disp.lnHM.log")
  names.combined_DE <- c("voom", "mean.lnHM.log")
  names.combined_DD <- c("diffVar", "disp.MDSeq.zi", "disp.lnHM.log")
  names.combined_long <- c("voom_diffVar", "voom_disp.MDSeq.zi", 
                           "voom_disp.lnHM.log", 
                           "mean.lnHM.log_diffVar", "mean.lnHM.log_disp.MDSeq.zi", 
                           "mean.lnHM.log_disp.lnHM.log")
  names.combined <- c("voom_dV", "voom_MD", "voom_HM", 
                      "HM_dV", "HM_MD", "HM_HM")
  names.diffVar_DEDD <- c("diffVar_DEDD")
  results.list <- c(names.edgeR, names.DESeq2, names.voom, 
                    names.DSS,
                    names.baySeq, 
                    names.diffVar, names.MDSeq, names.expHMM, names.expHM, 
                    names.lnHMM, names.lnHM, names.combined, names.diffVar_DEDD)
}


for (samples_per_group in c(2, 5, 10, 20, 50)) {
  {
    for (i in c(names.edgeR, names.DESeq2, names.voom, names.DSS, names.baySeq, names.diffVar, 
                names.MDSeq, names.expHM, names.lnHM, names.combined, names.diffVar_DEDD)) {
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
  folder <- "Results/GTEx blood simulated DE, DD, DEDD results Feb 2020"
  for (i in 1:10) {
    filename <- paste0("results.blood_", samples_per_group, "_set", i, "_DEDD.rds")
    res <- readRDS(here(folder, filename))
    DE <- factor(res$DE, levels=c("1", "0"))
    DD <- factor(res$DD, levels=c("1", "0"))
    DEDD <- factor(as.numeric(res$DE==1 | res$DD==1), levels=c("1", "0"))
    for (k in names.combined_DE) {
      for (l in names.combined_DD) {
        assign(paste0("p.", k, "_", l), pmin(get(paste0("p.", k, ".tmm"), res), get(paste0("p.", l, ".tmm"), res)))
        assign(paste0("q.", k, "_", l), p.adjust(get(paste0("p.", k, "_", l)), method="BH"))
      }
    }
    for (j in seq_len(length(names.combined))) {
      assign(paste0("p.", names.combined[j]), get(paste0("p.", names.combined_long[j])))
      assign(paste0("q.", names.combined[j]), get(paste0("q.", names.combined_long[j])))
    }
    for (j in c(names.edgeR, names.DESeq2, names.voom, names.DSS,
                names.diffVar, names.MDSeq, names.expHM, names.lnHM)) {
      if (mean(is.na(get(paste0("p.", j, ".tmm"), res))) < 1) { # added because no data for DSS for some runs
        assign(paste0("call.", j), factor(as.numeric(get(paste0("q.", j, ".tmm"), res) < 0.05), 
                                          levels=c("1", "0")))
        if (grepl("disp", j) | grepl("diffVar", j)) {
          assign(paste0("pred.", j), 
                 prediction(1 - get(paste0("p.", j, ".tmm"), res)[which(!is.na(get(paste0("p.",j, ".tmm"), res)))], 
                            DD[which(!is.na(get(paste0("p.",j, ".tmm"), res)))]))
          assign(paste0("fpr.", j), `[<-`(get(paste0("fpr.", j)), i,
                                          value=1 - specificity(get(paste0("call.", j)), DD)))
          assign(paste0("fdr.", j), `[<-`(get(paste0("fdr.", j)), i,
                                          value=1 - precision(get(paste0("call.", j)), DD)))
          assign(paste0("tpr.", j), `[<-`(get(paste0("tpr.", j)), i,
                                          value=sensitivity(get(paste0("call.", j)), DD)))
        }
        else {
          assign(paste0("pred.", j), 
                 prediction(1 - get(paste0("p.",j, ".tmm"), res)[which(!is.na(get(paste0("p.",j, ".tmm"), res)))], 
                            DE[which(!is.na(get(paste0("p.",j, ".tmm"), res)))]))
          assign(paste0("fpr.", j), `[<-`(get(paste0("fpr.", j)), i,
                                          value=1 - specificity(get(paste0("call.", j)), DE)))
          assign(paste0("fdr.", j), `[<-`(get(paste0("fdr.", j)), i,
                                          value=1 - precision(get(paste0("call.", j)), DE)))
          assign(paste0("tpr.", j), `[<-`(get(paste0("tpr.", j)), i,
                                          value=sensitivity(get(paste0("call.", j)), DE)))
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
    for (j in names.DSS) {
      if (mean(is.na(get(paste0("p.", j, ".tmm"), res))) < 1) {
        assign(paste0("call.lfdr.", j), factor(as.numeric(get(paste0("lfdr.", j, ".tmm"), res) < 0.05),
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
      assign(paste0("call.", j), factor(as.numeric(get(paste0("q.", j, ".tmm"), res) < 0.05), 
                                        levels=c("1", "0")))
      assign(paste0("pred.", j), prediction(get(paste0("prob.",j, ".tmm"), res), DE))
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
      assign(paste0("call.point5.", j), factor(as.numeric(get(paste0("prob.", j, ".tmm"), res) > 0.5), 
                                               levels=c("1", "0")))
      assign(paste0("call.thr.", j), 
             factor(as.numeric(get(paste0("prob.", j, ".tmm"), res) > 
                                 get(paste0("thr.", j, ".tmm"), res)), levels=c("1", "0")))
      assign(paste0("call.bfdr.", j), factor(as.numeric(get(paste0("bfdr.", j, ".tmm"), res) < 0.05), 
                                             levels=c("1", "0")))
      assign(paste0("pred.", j), prediction(get(paste0("prob.", j, ".tmm"), res), DEDD))
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
    for (j in names.combined) {
      assign(paste0("call.", j), factor(as.numeric(get(paste0("q.", j)) < 0.05), 
                                        levels=c("1", "0")))
      assign(paste0("pred.", j), 
             prediction(1 - get(paste0("p.", j))[which(!is.na(get(paste0("p.",j))))], 
                        DEDD[which(!is.na(get(paste0("p.",j))))]))
      assign(paste0("fpr.", j), `[<-`(get(paste0("fpr.", j)), i, 
                                      value=1 - specificity(get(paste0("call.",j)), DEDD)))
      assign(paste0("fdr.", j), `[<-`(get(paste0("fdr.", j)), i, 
                                      value=1 - precision(get(paste0("call.",j)), DEDD)))
      assign(paste0("tpr.", j), `[<-`(get(paste0("tpr.", j)), i, 
                                      value=sensitivity(get(paste0("call.",j)), DEDD)))
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
    for (j in names.diffVar) {
      assign(paste0("pred.", j, "_DEDD"), prediction(1 - get(paste0("p.", j, ".tmm"), res), DEDD))
      assign(paste0("fpr.", j, "_DEDD"), `[<-`(get(paste0("fpr.", j, "_DEDD")), i, 
                                               value=1 - specificity(get(paste0("call.",j)), DEDD)))
      assign(paste0("fdr.", j, "_DEDD"), `[<-`(get(paste0("fdr.", j, "_DEDD")), i, 
                                               value=1 - precision(get(paste0("call.",j)), DEDD)))
      assign(paste0("tpr.", j, "_DEDD"), `[<-`(get(paste0("tpr.", j, "_DEDD")), i, 
                                               value=sensitivity(get(paste0("call.",j)), DEDD)))
      assign(paste0("auc.", j, "_DEDD"), `[<-`(get(paste0("auc.", j, "_DEDD")), i, 
                                               value=performance(get(paste0("pred.", j, "_DEDD")), 
                                                                 measure="auc")@y.values[[1]]))
      assign(paste0("pauc.", j, "_DEDD"), `[<-`(get(paste0("pauc.", j, "_DEDD")), i,
                                                value=performance(get(paste0("pred.", j, "_DEDD")), measure="auc", 
                                                                  fpr.stop=0.05)@y.values[[1]]))
      assign(paste0("false.discoveries.", j, "_DEDD"), `[<-`(get(paste0("false.discoveries.", j, "_DEDD")), i,, 
                                                             value=c(get(paste0("pred.", j, "_DEDD"))@fp[[1]], 
                                                                     rep(NA, 20000 - 
                                                                           length(get(paste0("pred.", j, 
                                                                                             "_DEDD"))@fp[[1]])))))
      assign(paste0("discoveries.", j, "_DEDD"), `[<-`(get(paste0("discoveries.", j, "_DEDD")), i,, 
                                                       value=c(get(paste0("pred.", j, "_DEDD"))@n.pos.pred[[1]], 
                                                               rep(NA, 20000 - 
                                                                     length(get(paste0("pred." ,j, 
                                                                                       "_DEDD"))@n.pos.pred[[1]])))))
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
    rm(list=c("i", "j", "filename", "res", "DE", "DD", "DEDD"))
  }
  
  # Re-format results to save by experiment type ####
  
  # Set up variables
  {
    names.DE <- c("edgeR.ql", "edgeR.lr", "edgeR.et", "DESeq2.noif", "DESeq2.if", 
                  "voom", "DSS", "baySeq", "mean.MDSeq.zi", "mean.MDSeq.nozi", 
                  "mean.expHM.untr", "mean.expHM.log", "mean.lnHM.untr", "mean.lnHM.log")
    names.DD <- c("diffVar", "disp.MDSeq.zi", "disp.MDSeq.nozi", 
                  "disp.expHM.untr", "disp.expHM.log", "disp.lnHM.untr", "disp.lnHM.log")
    names.DEDD <- c("diffVar_DEDD", "expHMM", "lnHMM", "voom_dV", "voom_MD", "voom_HM", 
                    "HM_dV", "HM_MD", "HM_HM")
  }
  
  for (i in c("DE", "DD", "DEDD")) {
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
                names.diffVar, names.MDSeq, names.expHM, names.lnHM, names.combined, 
                names.diffVar_DEDD)) {
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
        if (!grepl("_", j)) {
          results.DE$fpr[[j]] <- get(paste0("fpr.", j))
          results.DE$fdr[[j]] <- get(paste0("fdr.", j))
          results.DE$tpr[[j]] <- get(paste0("tpr.", j))
        }
      }
    }
    for (i in names.DD) {
      results.DD$auc[[i]] <- get(paste0("auc.", i))
      results.DD$pauc[[i]] <- get(paste0("pauc.", i))
      results.DD$mean.fdr[[i]] <- get(paste0("mean.fdr.", i))
      results.DD$mean.discoveries[[i]] <- get(paste0("mean.discoveries.", i))
      for (j in grep(i, names.fpr.fdr.tpr[-grep("DEDD", names.fpr.fdr.tpr)], value=T)) {
        results.DD$fpr[[j]] <- get(paste0("fpr.", j))
        results.DD$fdr[[j]] <- get(paste0("fdr.", j))
        results.DD$tpr[[j]] <- get(paste0("tpr.", j))
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
  
  folder <- "Results/GTEx blood updated results Sept 2020"
  saveRDS(results.DE, file=here(folder, paste0("DE.results.blood_", samples_per_group, "_DEDD.rds")))
  saveRDS(results.DD, file=here(folder, paste0("DD.results.blood_", samples_per_group, "_DEDD.rds")))
  saveRDS(results.DEDD, file=here(folder, paste0("DEDD.results.blood_", samples_per_group, "_DEDD.rds")))
}
