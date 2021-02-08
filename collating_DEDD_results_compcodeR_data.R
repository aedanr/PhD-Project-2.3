library(here)
library(compcodeR)
library(ROCR)
library(caret)

# First create variables to store results
{
  names.edgeR <- c("edgeR.ql.tmm", "lfc1.edgeR.ql.tmm", "lfc2.edgeR.ql.tmm", 
                   "edgeR.lr.tmm", "lfc1.edgeR.lr.tmm", "lfc2.edgeR.lr.tmm", 
                   "edgeR.et.tmm", 
                   "edgeR.ql.rle", "lfc1.edgeR.ql.rle", "lfc2.edgeR.ql.rle", 
                   "edgeR.lr.rle", "lfc1.edgeR.lr.rle", "lfc2.edgeR.lr.rle", 
                   "edgeR.et.rle")
  names.DESeq2 <- c("DESeq2.noif.tmm", "lfc1.DESeq2.noif.tmm", "lfc2.DESeq2.noif.tmm", 
                    "DESeq2.if.tmm", "lfc1.DESeq2.if.tmm", "lfc2.DESeq2.if.tmm", 
                    "DESeq2.noif.rle", "lfc1.DESeq2.noif.rle", "lfc2.DESeq2.noif.rle", 
                    "DESeq2.if.rle", "lfc1.DESeq2.if.rle", "lfc2.DESeq2.if.rle")
  names.voom <- c("voom.tmm", "lfc1.voom.tmm", "lfc2.voom.tmm", 
                  "voom.rle", "lfc1.voom.rle", "lfc2.voom.rle")
  names.DSS <- c("DSS.tmm", "DSS.rle")
  names.baySeq <- c("baySeq.tmm", "baySeq.rle")
  names.MDSeq <- c("mean.MDSeq.zi.tmm", "mean.lfc1.MDSeq.zi.tmm", "mean.lfc2.MDSeq.zi.tmm", 
                   "mean.MDSeq.nozi.tmm", "mean.lfc1.MDSeq.nozi.tmm", "mean.lfc2.MDSeq.nozi.tmm", 
                   "disp.MDSeq.zi.tmm", "disp.lfc1.MDSeq.zi.tmm", "disp.lfc2.MDSeq.zi.tmm", 
                   "disp.MDSeq.nozi.tmm", "disp.lfc1.MDSeq.nozi.tmm", "disp.lfc2.MDSeq.nozi.tmm", 
                   "mean.MDSeq.zi.rle", "mean.lfc1.MDSeq.zi.rle", "mean.lfc2.MDSeq.zi.rle", 
                   "mean.MDSeq.nozi.rle", "mean.lfc1.MDSeq.nozi.rle", "mean.lfc2.MDSeq.nozi.rle", 
                   "disp.MDSeq.zi.rle", "disp.lfc1.MDSeq.zi.rle", "disp.lfc2.MDSeq.zi.rle", 
                   "disp.MDSeq.nozi.rle", "disp.lfc1.MDSeq.nozi.rle", "disp.lfc2.MDSeq.nozi.rle")
  names.expHMM <- c("expHMM.tmm", "expHMM.rle")
  names.expHM <- c("mean.expHM.untr.tmm", "mean.expHM.log.tmm", 
                   "mean.lfc1.expHM.tmm", "mean.lfc2.expHM.tmm", 
                   "disp.expHM.untr.tmm", "disp.expHM.log.tmm", 
                   "disp.lfc1.expHM.tmm", "disp.lfc2.expHM.tmm", 
                   "mean.expHM.untr.rle", "mean.expHM.log.rle", 
                   "mean.lfc1.expHM.rle", "mean.lfc2.expHM.rle", 
                   "disp.expHM.untr.rle", "disp.expHM.log.rle", 
                   "disp.lfc1.expHM.rle", "disp.lfc2.expHM.rle")
  names.lnHMM <- c("lnHMM.tmm", "lnHMM.rle")
  names.lnHM <- c("mean.lnHM.untr.tmm", "mean.lnHM.log.tmm", 
                  "mean.lfc1.lnHM.tmm", "mean.lfc2.lnHM.tmm", 
                  "disp.lnHM.untr.tmm", "disp.lnHM.log.tmm", 
                  "disp.lfc1.lnHM.tmm", "disp.lfc2.lnHM.tmm", 
                  "mean.lnHM.untr.rle", "mean.lnHM.log.rle", 
                  "mean.lfc1.lnHM.rle", "mean.lfc2.lnHM.rle", 
                  "disp.lnHM.untr.rle", "disp.lnHM.log.rle", 
                  "disp.lfc1.lnHM.rle", "disp.lfc2.lnHM.rle")
  names.combined_DE <- c("edgeR.ql", "mean.lnHM.untr")
  names.combined_DD <- c("disp.MDSeq.zi", "disp.lnHM.untr")
  names.combined_long <- c("edgeR.ql_disp.MDSeq.zi.tmm", "edgeR.ql_disp.lnHM.untr.tmm", 
                           "mean.lnHM.untr_disp.MDSeq.zi.tmm", "mean.lnHM.untr_disp.lnHM.untr.tmm", 
                           "edgeR.ql_disp.MDSeq.zi.rle", "edgeR.ql_disp.lnHM.untr.rle", 
                           "mean.lnHM.untr_disp.MDSeq.zi.rle", "mean.lnHM.untr_disp.lnHM.untr.rle")
  names.combined <- c("edgeR_MD.tmm", "edgeR_HM.tmm", "HM_MD.tmm", "HM.HM.tmm", 
                      "edgeR_MD.rle", "edgeR_HM.rle", "HM_MD.rle", "HM.HM.rle")
  # results.list <- c(names.edgeR, names.DESeq2, names.voom, names.DSS, names.baySeq, 
  #                   names.MDSeq, names.expHMM, names.expHM, names.lnHMM, names.lnHM, 
  #                   names.combined)
  results.list <- c(names.edgeR, names.DESeq2, names.voom, names.baySeq, 
                    names.MDSeq, names.expHMM, names.expHM, names.lnHMM, names.lnHM, 
                    names.combined)
}


for (samples_per_group in c(50)) {
  {
    for (i in c(names.edgeR, names.DESeq2, names.voom, 
                # names.DSS, 
                names.baySeq, 
                names.MDSeq, names.expHM, names.lnHM, names.combined)) {
      for (j in c("fpr.", "fdr.", "tpr.")) {
        assign(paste0(j, i), numeric(50))
      }
    }
    # for (i in names.DSS) {
    #   for (j in c("fpr.lfdr.", "fdr.lfdr.", "tpr.lfdr.")) {
    #     assign(paste0(j, i), numeric(50))
    #   }
    # }
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
    filename <- paste0("results.DEDD", samples_per_group, ".", i, ".rds")
    res <- readRDS(here(folder, filename))
    DE <- factor(res$DE, levels=c("1", "0"))
    DD <- factor(res$DD, levels=c("1", "0"))
    DEDD <- factor(as.numeric(res$DE==1 | res$DD==1), levels=c("1", "0"))
    lfcm1 <- factor(as.numeric(res$lfcm1), levels=c("1", "0"))
    lfcm2 <- factor(as.numeric(res$lfcm2), levels=c("1", "0"))
    lfcd1 <- factor(as.numeric(res$lfcd1), levels=c("1", "0"))
    lfcd2 <- factor(as.numeric(res$lfcd2), levels=c("1", "0"))
    for (j in c(".tmm", ".rle")) {
      for (k in names.combined_DE) {
        for (l in names.combined_DD) {
          assign(paste0("p.", k, "_", l, j), pmin(get(paste0("p.", k, j), res), get(paste0("p.", l, j), res)))
          assign(paste0("q.", k, "_", l, j), p.adjust(get(paste0("p.", k, "_", l, j)), method="BH"))
        }
      }
    }
    for (j in seq_len(length(names.combined))) {
      assign(paste0("p.", names.combined[j]), get(paste0("p.", names.combined_long[j])))
      assign(paste0("q.", names.combined[j]), get(paste0("q.", names.combined_long[j])))
    }
    for (j in c(names.edgeR, names.DESeq2, names.voom, 
                # names.DSS, 
                names.MDSeq, 
                names.expHM, names.lnHM)) {
      assign(paste0("call.", j), factor(as.numeric(get(paste0("q.", j), res) < 0.05), 
                                        levels=c("1", "0")))
      if (grepl("disp", j)) {
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
      }
      else
        if (grepl("lfc1", j)) {
          assign(paste0("pred.", j), prediction(1 - get(paste0("p.", j), res)[which(!is.na(get(paste0("p.", j), res)))], 
                                                lfcm1[which(!is.na(get(paste0("p.", j), res)))]))
          assign(paste0("fpr.", j), `[<-`(get(paste0("fpr.",j)), i,
                                          value=1 - specificity(get(paste0("call.", j)), lfcm1)))
          assign(paste0("fdr.", j), `[<-`(get(paste0("fdr.",j)), i,
                                          value=1 - precision(get(paste0("call.", j)), lfcm1)))
          assign(paste0("tpr.", j), `[<-`(get(paste0("tpr.",j)), i,
                                          value=sensitivity(get(paste0("call.", j)), lfcm1)))
        }
      else if (grepl("lfc2", j)) {
        assign(paste0("pred.", j), prediction(1 - get(paste0("p.", j), res)[which(!is.na(get(paste0("p.", j), res)))], 
                                              lfcm2[which(!is.na(get(paste0("p.", j), res)))]))
        assign(paste0("fpr.", j), `[<-`(get(paste0("fpr.", j)), i,
                                        value=1 - specificity(get(paste0("call.", j)), lfcm2)))
        assign(paste0("fdr.", j), `[<-`(get(paste0("fdr.", j)), i,
                                        value=1 - precision(get(paste0("call.", j)), lfcm2)))
        assign(paste0("tpr.", j), `[<-`(get(paste0("tpr.", j)), i,
                                        value=sensitivity(get(paste0("call.", j)), lfcm2)))
      }
      else {
        assign(paste0("pred.", j), prediction(1 - get(paste0("p.",j), res)[which(!is.na(get(paste0("p.", j), res)))], 
                                              DE[which(!is.na(get(paste0("p.", j), res)))]))
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
    # for (j in names.DSS) {
    #   assign(paste0("call.lfdr.", j), factor(as.numeric(get(paste0("lfdr.", j), res) < 0.05), 
    #                                          levels=c("1", "0")))
    #   assign(paste0("fpr.lfdr.", j), `[<-`(get(paste0("fpr.lfdr.",j)), i,
    #                                        value=1 - specificity(get(paste0("call.lfdr.", j)), DE)))
    #   assign(paste0("fdr.lfdr.", j), `[<-`(get(paste0("fdr.lfdr.",j)), i,
    #                                        value=1 - precision(get(paste0("call.lfdr.", j)), DE)))
    #   assign(paste0("tpr.lfdr.", j), `[<-`(get(paste0("tpr.lfdr.",j)), i,
    #                                        value=sensitivity(get(paste0("call.lfdr.", j)), DE)))
    # }
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
                                             value=1 - specificity(get(paste0("call.point5.",j)), DEDD)))
      assign(paste0("fpr.thr.", j), `[<-`(get(paste0("fpr.thr.", j)), i, 
                                          value=1 - specificity(get(paste0("call.thr.",j)), DEDD)))
      assign(paste0("fpr.bfdr.", j), `[<-`(get(paste0("fpr.bfdr.", j)), i, 
                                           value=1 - specificity(get(paste0("call.bfdr.",j)), DEDD)))
      assign(paste0("fdr.point5.", j), `[<-`(get(paste0("fdr.point5.", j)), i, 
                                             value=1 - precision(get(paste0("call.point5.",j)), DEDD)))
      assign(paste0("fdr.thr.", j), `[<-`(get(paste0("fdr.thr.", j)), i, 
                                          value=1 - precision(get(paste0("call.thr.",j)), DEDD)))
      assign(paste0("fdr.bfdr.", j), `[<-`(get(paste0("fdr.bfdr.", j)), i, 
                                           value=1 - precision(get(paste0("call.bfdr.",j)), DEDD)))
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
      assign(paste0("pred.", j), prediction(1 - get(paste0("p.", j))[which(!is.na(get(paste0("p.", j))))], 
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
    rm(list=c("i", "j", "filename", "res", "DE", "DD", "DEDD", "lfcm1", "lfcm2", "lfcd1", "lfcd2"))
  }
  
  # Re-format results to save by experiment type ####
  
  # Set up variables
  {
    names.DE <- c("edgeR.ql.tmm", "edgeR.lr.tmm", "edgeR.et.tmm", "DESeq2.noif.tmm", "DESeq2.if.tmm", 
                  "voom.tmm", 
                  # "DSS.tmm", 
                  "baySeq.tmm", "mean.MDSeq.zi.tmm", "mean.MDSeq.nozi.tmm", 
                  "mean.expHM.untr.tmm", "mean.expHM.log.tmm", "mean.lnHM.untr.tmm", "mean.lnHM.log.tmm", 
                  "edgeR.ql.rle", "edgeR.lr.rle", "edgeR.et.rle", "DESeq2.noif.rle", "DESeq2.if.rle", 
                  "voom.rle", 
                  # "DSS.rle", 
                  "baySeq.rle", "mean.MDSeq.zi.rle", "mean.MDSeq.nozi.rle", 
                  "mean.expHM.untr.rle", "mean.expHM.log.rle", "mean.lnHM.untr.rle", "mean.lnHM.log.rle")
    names.DE.lfc1 <- c("lfc1.edgeR.ql.tmm", "lfc1.edgeR.lr.tmm", "lfc1.DESeq2.noif.tmm", "lfc1.DESeq2.if.tmm", 
                       "lfc1.voom.tmm", "mean.lfc1.MDSeq.zi.tmm", "mean.lfc1.MDSeq.nozi.tmm", 
                       "mean.lfc1.expHM.tmm", "mean.lfc1.lnHM.tmm", 
                       "lfc1.edgeR.ql.rle", "lfc1.edgeR.lr.rle", "lfc1.DESeq2.noif.rle", "lfc1.DESeq2.if.rle", 
                       "lfc1.voom.rle", "mean.lfc1.MDSeq.zi.rle", "mean.lfc1.MDSeq.nozi.rle", 
                       "mean.lfc1.expHM.rle", "mean.lfc1.lnHM.rle")
    names.DE.lfc2 <- c("lfc2.edgeR.ql.tmm", "lfc2.edgeR.lr.tmm", "lfc2.DESeq2.noif.tmm", "lfc2.DESeq2.if.tmm", 
                       "lfc2.voom.tmm", "mean.lfc2.MDSeq.zi.tmm", "mean.lfc2.MDSeq.nozi.tmm", 
                       "mean.lfc2.expHM.tmm", "mean.lfc2.lnHM.tmm", 
                       "lfc2.edgeR.ql.rle", "lfc2.edgeR.lr.rle", "lfc2.DESeq2.noif.rle", "lfc2.DESeq2.if.rle", 
                       "lfc2.voom.rle", "mean.lfc2.MDSeq.zi.rle", "mean.lfc2.MDSeq.nozi.rle", 
                       "mean.lfc2.expHM.rle", "mean.lfc2.lnHM.rle")
    names.DD <- c("disp.MDSeq.zi.tmm", "disp.MDSeq.nozi.tmm", 
                  "disp.expHM.untr.tmm", "disp.expHM.log.tmm", "disp.lnHM.untr.tmm", "disp.lnHM.log.tmm", 
                  "disp.MDSeq.zi.rle", "disp.MDSeq.nozi.rle", 
                  "disp.expHM.untr.rle", "disp.expHM.log.rle", "disp.lnHM.untr.rle", "disp.lnHM.log.rle")
    names.DD.lfc1 <- c("disp.lfc1.MDSeq.zi.tmm", "disp.lfc1.MDSeq.nozi.tmm", 
                       "disp.lfc1.expHM.tmm", "disp.lfc1.lnHM.tmm", 
                       "disp.lfc1.MDSeq.zi.rle", "disp.lfc1.MDSeq.nozi.rle", 
                       "disp.lfc1.expHM.rle", "disp.lfc1.lnHM.rle")
    names.DD.lfc2 <- c("disp.lfc2.MDSeq.zi.tmm", "disp.lfc2.MDSeq.nozi.tmm", 
                       "disp.lfc2.expHM.tmm", "disp.lfc2.lnHM.tmm", 
                       "disp.lfc2.MDSeq.zi.rle", "disp.lfc2.MDSeq.nozi.rle", 
                       "disp.lfc2.expHM.rle", "disp.lfc2.lnHM.rle")
    names.DEDD <- c(names.expHMM, names.lnHMM, names.combined)
  }
  
  for (i in c("DE", "DE.lfc1", "DE.lfc2", "DD", "DD.lfc1", "DD.lfc2", "DEDD")) {
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
    for (i in c(names.edgeR, names.DESeq2, names.voom, 
                # names.DSS, 
                names.baySeq, 
                names.MDSeq, names.expHM, names.lnHM, names.combined)) {
      names.fpr.fdr.tpr <- c(names.fpr.fdr.tpr, i)
    }
    # for (i in names.DSS) {
    #   names.fpr.fdr.tpr <- c(names.fpr.fdr.tpr, paste0("lfdr.", i))
    # }
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
    for (i in names.DE.lfc1) {
      results.DE.lfc1$auc[[i]] <- get(paste0("auc.", i))
      results.DE.lfc1$pauc[[i]] <- get(paste0("pauc.", i))
      results.DE.lfc1$mean.fdr[[i]] <- get(paste0("mean.fdr.", i))
      results.DE.lfc1$mean.discoveries[[i]] <- get(paste0("mean.discoveries.", i))
      for (j in grep(i, names.fpr.fdr.tpr, value=T)) {
        results.DE.lfc1$fpr[[j]] <- get(paste0("fpr.", j))
        results.DE.lfc1$fdr[[j]] <- get(paste0("fdr.", j))
        results.DE.lfc1$tpr[[j]] <- get(paste0("tpr.", j))
      }
    }
    for (i in names.DE.lfc2) {
      results.DE.lfc2$auc[[i]] <- get(paste0("auc.", i))
      results.DE.lfc2$pauc[[i]] <- get(paste0("pauc.", i))
      results.DE.lfc2$mean.fdr[[i]] <- get(paste0("mean.fdr.", i))
      results.DE.lfc2$mean.discoveries[[i]] <- get(paste0("mean.discoveries.", i))
      for (j in grep(i, names.fpr.fdr.tpr, value=T)) {
        results.DE.lfc2$fpr[[j]] <- get(paste0("fpr.", j))
        results.DE.lfc2$fdr[[j]] <- get(paste0("fdr.", j))
        results.DE.lfc2$tpr[[j]] <- get(paste0("tpr.", j))
      }
    }
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
  
  # Fix results to ensure no extraneous results in each file ####
  # Extra results in some files due to imprecise naming meaning that matching to place 
  # results in the right files sometimes caught more data than it should have for fpr, 
  # fdr and tpr, e.g. "voom.lfc1" contains "voom", so lfc results for voom were included 
  # in DE results. There was also some doubling up, e.g. "noif.DESeq" 
  # contains "if.DESeq" as well as "noif.DESeq", but because results are filled by name, 
  # these duplications aren't in the results files as they're overwritten by the next 
  # match in the loop.
  {
    results.DE$fpr <- results.DE$fpr[, -grep("lfc", names(results.DE$fpr))]
    results.DE$fdr <- results.DE$fdr[, -grep("lfc", names(results.DE$fdr))]
    results.DE$tpr <- results.DE$tpr[, -grep("lfc", names(results.DE$tpr))]
  }
  
  saveRDS(results.DE, file=here(folder, paste0("DE.results.DEDD", samples_per_group, ".rds")))
  saveRDS(results.DE.lfc1, file=here(folder, paste0("DE.lfc1.results.DEDD", samples_per_group, ".rds")))
  saveRDS(results.DE.lfc2, file=here(folder, paste0("DE.lfc2.results.DEDD", samples_per_group, ".rds")))
  saveRDS(results.DD, file=here(folder, paste0("DD.results.DEDD", samples_per_group, ".rds")))
  saveRDS(results.DD.lfc1, file=here(folder, paste0("DD.lfc1.results.DEDD", samples_per_group, ".rds")))
  saveRDS(results.DD.lfc2, file=here(folder, paste0("DD.lfc2.results.DEDD", samples_per_group, ".rds")))
  saveRDS(results.DEDD, file=here(folder, paste0("DEDD.results.DEDD", samples_per_group, ".rds")))
}
