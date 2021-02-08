library(here)
library(compcodeR)
library(ROCR)
library(caret)

# First create variables to store results
{
  names.lnHMM <- c("lnHMM")
  names.lnHM <- c("mean.untr", "mean.log", "disp.untr", "disp.log")
  names.combined_DE <- c("mean.log")
  names.combined_DD <- c("disp.log")
  names.combined_long <- c("mean.log_disp.log")
  names.combined <- c("HM_HM")
  names.combined_all <- c("voom_HM", "HM_dV", "HM_MD", "HM_HM")
  results.list <- c(names.lnHMM, names.lnHM, names.combined_all)
}


for (samples_per_group in c(2, 5, 10, 20, 50)) {
  {
    for (i in c(names.lnHM, names.combined_all)) {
      for (j in c("fpr.", "fdr.", "tpr.")) {
        assign(paste0(j, i), rep(NA, 10))
      }
    }
    for (i in c(names.lnHMM)) {
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
  folder <- "Results/GTEx muscle long chain results Aug 2020"
  for (i in 1:10) {
    filename <- paste0("results.lnHM.muscle_", samples_per_group, "_set", i, "_DEDD.rds")
    res <- readRDS(here(folder, filename))
    DE <- factor(res$DE, levels=c("1", "0"))
    DD <- factor(res$DD, levels=c("1", "0"))
    DEDD <- factor(as.numeric(res$DE==1 | res$DD==1), levels=c("1", "0"))
    res_old <- readRDS(here("Results/GTEx muscle simulated DE, DD, DEDD results March 2020", 
                            paste0("results.muscle_", samples_per_group, "_set", i, "_DEDD.rds")))
    p.voom <- res_old$p.voom.tmm
    p.dV <- res_old$p.diffVar.tmm
    p.MD <- res_old$p.disp.MDSeq.zi.tmm
    rm(res_old)
    p.voom_HM <- pmin(p.voom, res$p.disp.log)
    p.HM_dV <- pmin(res$p.mean.log, p.dV)
    p.HM_MD <- pmin(res$p.mean.log, p.MD)
    rm(p.voom, p.dV, p.MD)
    for (j in c("voom_HM", "HM_dV", "HM_MD")) {
      assign(paste0("q.", j), p.adjust(get(paste0("p.", j)), method="BH"))
    }
    for (k in names.combined_DE) {
      for (l in names.combined_DD) {
        assign(paste0("p.", k, "_", l), pmin(get(paste0("p.", k), res), get(paste0("p.", l), res)))
        assign(paste0("q.", k, "_", l), p.adjust(get(paste0("p.", k, "_", l)), method="BH"))
      }
    }
    for (j in seq_len(length(names.combined))) {
      assign(paste0("p.", names.combined[j]), get(paste0("p.", names.combined_long[j])))
      assign(paste0("q.", names.combined[j]), get(paste0("q.", names.combined_long[j])))
    }
    for (j in c(names.lnHM)) {
      if (mean(is.na(get(paste0("p.",j), res))) < 1) { # added because no data for DSS for some runs
        assign(paste0("call.", j), factor(as.numeric(get(paste0("q.", j), res) < 0.05), 
                                          levels=c("1", "0")))
        if (grepl("disp", j)) {
          assign(paste0("pred.", j), prediction(1 - get(paste0("p.", j), res), DD))
          assign(paste0("fpr.", j), `[<-`(get(paste0("fpr.", j)), i,
                                          value=1 - specificity(get(paste0("call.", j)), DD)))
          assign(paste0("fdr.", j), `[<-`(get(paste0("fdr.", j)), i,
                                          value=1 - precision(get(paste0("call.", j)), DD)))
          assign(paste0("tpr.", j), `[<-`(get(paste0("tpr.", j)), i,
                                          value=sensitivity(get(paste0("call.", j)), DD)))
        }
        else {
          assign(paste0("pred.", j), prediction(1 - get(paste0("p.",j), res), DE))
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
    for (j in c(names.lnHMM)) {
      assign(paste0("call.point5.", j), factor(as.numeric(get("prob", res) > 0.5), 
                                               levels=c("1", "0")))
      assign(paste0("call.thr.", j), 
             factor(as.numeric(get("prob", res) > get("thr", res)), levels=c("1", "0")))
      assign(paste0("call.bfdr.", j), factor(as.numeric(get("bfdr", res) < 0.05), 
                                             levels=c("1", "0")))
      assign(paste0("pred.", j), prediction(get("prob", res), DEDD))
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
    for (j in names.combined_all) {
      assign(paste0("call.", j), factor(as.numeric(get(paste0("q.", j)) < 0.05), 
                                        levels=c("1", "0")))
      assign(paste0("pred.", j), 
             prediction(1 - get(paste0("p.", j))[which(!is.na(get(paste0("p.", j))))], 
                        DEDD[which(!is.na(get(paste0("p.", j))))]))
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
    names.DE <- c("mean.untr", "mean.log")
    names.DD <- c("disp.untr", "disp.log")
    names.DEDD <- c("lnHMM", "voom_HM", "HM_dV", "HM_MD", "HM_HM")
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
    for (i in c(names.lnHM, names.combined_all)) {
      names.fpr.fdr.tpr <- c(names.fpr.fdr.tpr, i)
    }
    for (i in c(names.lnHMM)) {
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
      results.DE$fpr[[i]] <- get(paste0("fpr.", i))
      results.DE$fdr[[i]] <- get(paste0("fdr.", i))
      results.DE$tpr[[i]] <- get(paste0("tpr.", i))
    }
    for (i in names.DD) {
      results.DD$auc[[i]] <- get(paste0("auc.", i))
      results.DD$pauc[[i]] <- get(paste0("pauc.", i))
      results.DD$mean.fdr[[i]] <- get(paste0("mean.fdr.", i))
      results.DD$mean.discoveries[[i]] <- get(paste0("mean.discoveries.", i))
      results.DD$fpr[[i]] <- get(paste0("fpr.", i))
      results.DD$fdr[[i]] <- get(paste0("fdr.", i))
      results.DD$tpr[[i]] <- get(paste0("tpr.", i))
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
  
  saveRDS(results.DE, file=here(folder, paste0("DE.results.lnHM.muscle_", samples_per_group, "_DEDD.rds")))
  saveRDS(results.DD, file=here(folder, paste0("DD.results.lnHM.muscle_", samples_per_group, "_DEDD.rds")))
  saveRDS(results.DEDD, file=here(folder, paste0("DEDD.results.lnHM.muscle_", samples_per_group, "_DEDD.rds")))
}
