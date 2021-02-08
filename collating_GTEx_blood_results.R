library(here)
library(ROCR)
library(caret)
source(here('scripts','2019-05-03_bfdr_function.R'))
DD4 <- readRDS(here("recount data/GTEx/quarantine", "results.blood_4_DD.rds"))
DD10 <- readRDS(here("recount data/GTEx/quarantine", "results.blood_10_DD.rds"))
DD20 <- readRDS(here("recount data/GTEx/quarantine", "results.blood_20_DD.rds"))
DD100 <- readRDS(here("recount data/GTEx/quarantine", "results.blood_100_DD.rds"))
DE4 <- readRDS(here("recount data/GTEx/quarantine", "results.blood_4_DE.rds"))
DE10 <- readRDS(here("recount data/GTEx/quarantine", "results.blood_10_DE.rds"))
DE20 <- readRDS(here("recount data/GTEx/quarantine", "results.blood_20_DE.rds"))
DE100 <- readRDS(here("recount data/GTEx/quarantine", "results.blood_100_DE.rds"))
DEDD4 <- readRDS(here("recount data/GTEx/quarantine", "results.blood_4_DEDD.rds"))
DEDD10 <- readRDS(here("recount data/GTEx/quarantine", "results.blood_10_DEDD.rds"))
DEDD20 <- readRDS(here("recount data/GTEx/quarantine", "results.blood_20_DEDD.rds"))
DEDD100 <- readRDS(here("recount data/GTEx/quarantine", "results.blood_100_DEDD.rds"))


# First create variables to store results
{
  names.edgeR <- c('ql.edgeR', 'lr.edgeR', 'et.edgeR')
  names.DESeq2 <- c('noif.DESeq2', 'if.DESeq2')
  names.voom <- c('voom')
  # names.DSS <- c('notrend.DSS', 'trend.DSS')
  names.DSS <- 'notrend.DSS'
  names.baySeq <- 'baySeq'
  # names.MDSeq <- c('disp.zi.MDSeq', 'disp.nozi.MDSeq')
  names.MDSeq <- c('mean.zi.MDSeq', 'mean.nozi.MDSeq')
  # names.MDSeq <- c('mean.zi.MDSeq', 'mean.nozi.MDSeq', 'disp.zi.MDSeq', 'disp.nozi.MDSeq')
  # names.expHM <- c('expHMM', 'disp.expHM', 'ldisp.expHM')
  names.expHM <- c('expHMM', 'mean.expHM', 'lmean.expHM')
  # names.expHM <- c('expHMM', 'mean.expHM', 'lmean.expHM', 'disp.expHM', 'ldisp.expHM')
  # names.lnHM <- c('lnHMM', 'disp.lnHM', 'ldisp.lnHM')
  names.lnHM <- c('lnHMM', 'mean.lnHM', 'lmean.lnHM')
  # names.lnHM <- c('lnHMM', 'mean.lnHM', 'lmean.lnHM', 'disp.lnHM', 'ldisp.lnHM')
  # names.combined <- c('expHM_expHM', 'lnHM_lnHM', 'MDSeq_MDSeq', 'voom_expHM', 'voom_lnHM', 'voom_MDSeq')
  # results.list <- c(
  #   names.MDSeq,
  #   names.expHM, names.lnHM)
  results.list <- c(names.edgeR, names.DESeq2, names.voom, names.DSS, names.baySeq,
                    names.MDSeq, names.expHM, names.lnHM)
}

{
  for (i in c(
    names.edgeR, names.DESeq2, names.voom, names.DSS, names.baySeq,
    # names.combined
    # ,
    names.MDSeq, names.expHM[-1], names.lnHM[-1]
    )) {
    for (j in c('fpr.','fdr.','tpr.')) {
      assign(paste0(j,i), numeric(4))
    }
  }
  for (i in names.DSS) {
   for (j in c('fpr.lfdr.','fdr.lfdr.','tpr.lfdr.')) {
     assign(paste0(j,i), numeric(4))
   }
  }
  for (i in c(names.expHM[1],names.lnHM[1])) {
    for (j in c('fpr.5.','fdr.5.','tpr.5.','fpr.thr.','fdr.thr.','tpr.thr.', 
                'fpr.bfdr.','fdr.bfdr.','tpr.bfdr.')) {
      assign(paste0(j,i), numeric(4))
    }
  }
  for (i in c(results.list
  # , names.combined
  )) {
    assign(paste0('auc.',i), numeric(4))
    assign(paste0('discoveries.',i), matrix(nrow=4, ncol=20000))
    assign(paste0('false.discoveries.',i), matrix(nrow=4, ncol=20000))
  }
}


# Then put results in variables
for (i in c(4,10,20,100)) {
  if (i == 4) {k <- 1} else if (i == 10) {k <- 2} else if (i == 20) {k <- 3} else if (i == 100) {k <- 4}
  # res <- get(paste0('DD', i))
  res <- get(paste0('DE', i))
  # res <- get(paste0('DEDD', i))
  # DD <- factor(res$DD, levels=c('1','0'))
  DE <- factor(res$DE, levels=c('1','0'))
  # DEDD <- DD
  DEDD <- DE
  # DEDD <- factor(res$DEDD, levels=c('1','0'))
  thr.expHM <- sort(res$prob.expHM, decreasing=T)[round(nrow(res$counts) * res$prop.expHM)]
  thr.lnHM <- sort(res$prob.lnHM, decreasing=T)[round(nrow(res$counts) * res$prop.lnHM)]
  for (j in c('expHM', 'lnHM')) {
    assign(paste0('bfdr.', j), bfdr(get(paste0('prob.', j), res)))
  }
  # expHM_expHM <- pmin(res$p.mean.expHM, res$p.disp.expHM)
  # lnHM_lnHM <- pmin(res$p.mean.lnHM, res$p.disp.lnHM)
  # MDSeq_MDSeq <- pmin(res$p.mean.zi.MDSeq, res$p.disp.zi.MDSeq)
  # voom_expHM <- pmin(res$p.voom, res$p.disp.expHM)
  # voom_lnHM <- pmin(res$p.voom, res$p.disp.lnHM)
  # voom_MDSeq <- pmin(res$p.voom, res$p.disp.zi.MDSeq)
  # call.MDSeq_MDSeq <- factor(as.numeric(p.adjust(MDSeq_MDSeq, method="BH") < 0.05), levels=c("1", "0"))
  # call.lnHM_lnHM <- factor(as.numeric(p.adjust(lnHM_lnHM, method="BH") < 0.05), levels=c("1", "0"))
  # call.expHM_expHM <- factor(as.numeric(p.adjust(expHM_expHM, method="BH") < 0.05), levels=c("1", "0"))
  # call.voom_MDSeq <- factor(as.numeric(p.adjust(voom_MDSeq, method="BH") < 0.05), levels=c("1", "0"))
  # call.voom_lnHM <- factor(as.numeric(p.adjust(voom_lnHM, method="BH") < 0.05), levels=c("1", "0"))
  # call.voom_expHM <- factor(as.numeric(p.adjust(voom_expHM, method="BH") < 0.05), levels=c("1", "0"))
  for (j in c(
    names.edgeR, names.DESeq2, names.voom, names.DSS,
    names.MDSeq, names.expHM[-1], names.lnHM[-1])) {
    assign(paste0('call.',j), factor(as.numeric(get(paste0('q.',j), res) < 0.05), levels=c('1','0')))
    # if (grepl('disp',j)) {
      # assign(paste0('pred.',j), prediction(1-get(paste0('p.',j), res), DD))
      # assign(paste0('fpr.',j), `[<-`(get(paste0('fpr.',j)), k,
      #                                    value=1-specificity(get(paste0('call.',j)), DD)))
      # assign(paste0('fdr.',j), `[<-`(get(paste0('fdr.',j)), k,
      #                                    value=1-precision(get(paste0('call.',j)), DD)))
      # assign(paste0('tpr.',j), `[<-`(get(paste0('tpr.',j)), k,
      #                                    value=sensitivity(get(paste0('call.',j)), DD)))
    #  }
    # else {
      assign(paste0('pred.',j), prediction(1-get(paste0('p.',j), res), DE))
      assign(paste0('fpr.',j), `[<-`(get(paste0('fpr.',j)), k,
                                         value=1-specificity(get(paste0('call.',j)), DE)))
      assign(paste0('fdr.',j), `[<-`(get(paste0('fdr.',j)), k,
                                         value=1-precision(get(paste0('call.',j)), DE)))
      assign(paste0('tpr.',j), `[<-`(get(paste0('tpr.',j)), k,
                                         value=sensitivity(get(paste0('call.',j)), DE)))
    # }
    assign(paste0('auc.',j), `[<-`(get(paste0('auc.',j)), k,
                                   value=performance(get(paste0('pred.',j)), measure='auc')@y.values[[1]]))
    assign(paste0('false.discoveries.',j), `[<-`(get(paste0('false.discoveries.',j)), k,,
                                                 value=c(get(paste0('pred.',j))@fp[[1]],
                                                         rep(NA, 20000-length(get(paste0('pred.',j))@fp[[1]])))))
    assign(paste0('discoveries.',j), `[<-`(get(paste0('discoveries.',j)), k,,
                                           value=c(get(paste0('pred.',j))@n.pos.pred[[1]],
                                                   rep(NA, 20000-length(get(paste0('pred.',j))@n.pos.pred[[1]])))))
  }
  for (j in names.DSS) {
    assign(paste0('call.lfdr.',j), factor(as.numeric(get(paste0('lfdr.',j), res) < 0.05), levels=c('1','0')))
    assign(paste0('fpr.lfdr.',j), `[<-`(get(paste0('fpr.lfdr.',j)), k,
                                        value=1-specificity(get(paste0('call.lfdr.',j)), DE)))
    assign(paste0('fdr.lfdr.',j), `[<-`(get(paste0('fdr.lfdr.',j)), k,
                                        value=1-precision(get(paste0('call.lfdr.',j)), DE)))
    assign(paste0('tpr.lfdr.',j), `[<-`(get(paste0('tpr.lfdr.',j)), k,
                                        value=sensitivity(get(paste0('call.lfdr.',j)), DE)))
  }
  for (j in names.baySeq) {
    assign(paste0('call.',j), factor(as.numeric(get(paste0('q.',j), res) < 0.05), levels=c('1','0')))
    assign(paste0('pred.',j), prediction(get(paste0('prob.',j), res), DE))
    assign(paste0('fpr.',j), `[<-`(get(paste0('fpr.',j)), k,
                                       value=1-specificity(get(paste0('call.',j)), DE)))
    assign(paste0('fdr.',j), `[<-`(get(paste0('fdr.',j)), k,
                                       value=1-precision(get(paste0('call.',j)), DE)))
    assign(paste0('tpr.',j), `[<-`(get(paste0('tpr.',j)), k,
                                       value=sensitivity(get(paste0('call.',j)), DE)))
    assign(paste0('auc.',j), `[<-`(get(paste0('auc.',j)), k,
                                   value=performance(get(paste0('pred.',j)), measure='auc')@y.values[[1]]))
    assign(paste0('false.discoveries.',j), `[<-`(get(paste0('false.discoveries.',j)), k,,
                                                 value=c(get(paste0('pred.',j))@fp[[1]],
                                                         rep(NA, 20000-length(get(paste0('pred.',j))@fp[[1]])))))
    assign(paste0('discoveries.',j), `[<-`(get(paste0('discoveries.',j)), k,,
                                           value=c(get(paste0('pred.',j))@n.pos.pred[[1]],
                                                   rep(NA, 20000-length(get(paste0('pred.',j))@n.pos.pred[[1]])))))
  }
  for (j in c('expHM', 'lnHM')) {
    assign(paste0('call.5.', j), factor(as.numeric(get(paste0('prob.',j), res) > 0.5), levels=c('1','0')))
    assign(paste0('call.thr.', j), factor(as.numeric(get(paste0('prob.',j), res) > get(paste0('thr.',j))), 
                                         levels=c('1','0')))
    assign(paste0('call.bfdr.', j), factor(as.numeric(get(paste0('bfdr.', j)) < 0.05), levels=c('1','0')))
    assign(paste0('pred.', paste0(j, 'M')), prediction(get(paste0('prob.',j), res), DEDD))
    assign(paste0('fpr.5.', paste0(j, 'M')), `[<-`(get(paste0('fpr.5.', paste0(j, 'M'))), k, 
                                     value=1-specificity(get(paste0('call.5.',j)), DEDD)))
    assign(paste0('fpr.thr.', paste0(j, 'M')), `[<-`(get(paste0('fpr.thr.', paste0(j, 'M'))), k, 
                                       value=1-specificity(get(paste0('call.thr.',j)), DEDD)))
    assign(paste0('fpr.bfdr.', paste0(j, 'M')), `[<-`(get(paste0('fpr.bfdr.', paste0(j, 'M'))), k, 
                                        value=1-specificity(get(paste0('call.bfdr.',j)), DEDD)))
    assign(paste0('fdr.5.', paste0(j, 'M')), `[<-`(get(paste0('fdr.5.', paste0(j, 'M'))), k, 
                                     value=1-precision(get(paste0('call.5.',j)), DEDD)))
    assign(paste0('fdr.thr.', paste0(j, 'M')), `[<-`(get(paste0('fdr.thr.', paste0(j, 'M'))), k, 
                                       value=1-precision(get(paste0('call.thr.',j)), DEDD)))
    assign(paste0('fdr.bfdr.', paste0(j, 'M')), `[<-`(get(paste0('fdr.bfdr.', paste0(j, 'M'))), k, 
                                        value=1-precision(get(paste0('call.bfdr.',j)), DEDD)))
    assign(paste0('tpr.5.', paste0(j, 'M')), `[<-`(get(paste0('tpr.5.', paste0(j, 'M'))), k, 
                                     value=sensitivity(get(paste0('call.5.',j)), DEDD)))
    assign(paste0('tpr.thr.', paste0(j, 'M')), `[<-`(get(paste0('tpr.thr.', paste0(j, 'M'))), k, 
                                       value=sensitivity(get(paste0('call.thr.',j)), DEDD)))
    assign(paste0('tpr.bfdr.', paste0(j, 'M')), `[<-`(get(paste0('tpr.bfdr.', paste0(j, 'M'))), k, 
                                        value=sensitivity(get(paste0('call.bfdr.',j)), DEDD)))
    assign(paste0('auc.', paste0(j, 'M')), `[<-`(get(paste0('auc.', paste0(j, 'M'))), k, 
                                   value=performance(get(paste0('pred.', paste0(j, 'M'))), measure='auc')@y.values[[1]]))
    assign(paste0('false.discoveries.', paste0(j, 'M')), `[<-`(get(paste0('false.discoveries.', paste0(j, 'M'))), k,, 
                                                 value=c(get(paste0('pred.', paste0(j, 'M')))@fp[[1]], 
                                                         rep(NA, 20000-length(get(paste0('pred.', paste0(j, 'M')))@fp[[1]])))))
    assign(paste0('discoveries.', paste0(j, 'M')), `[<-`(get(paste0('discoveries.', paste0(j, 'M'))), k,, 
                                           value=c(get(paste0('pred.', paste0(j, 'M')))@n.pos.pred[[1]], 
                                                   rep(NA, 20000-length(get(paste0('pred.', paste0(j, 'M')))@n.pos.pred[[1]])))))
  }
  # for (j in names.combined) {
  #   assign(paste0('pred.',j), prediction(1-get(j), DEDD))
  #   assign(paste0('fpr.', j), `[<-`(get(paste0('fpr.', j)), k, 
  #                                     value=1-specificity(get(paste0('call.',j)), DEDD)))
  #   assign(paste0('fdr.', j), `[<-`(get(paste0('fdr.', j)), k, 
  #                                     value=1-precision(get(paste0('call.',j)), DEDD)))
  #   assign(paste0('tpr.', j), `[<-`(get(paste0('tpr.', j)), k, 
  #                                     value=sensitivity(get(paste0('call.',j)), DEDD)))
  #   assign(paste0('auc.', j), `[<-`(get(paste0('auc.', j)), k, 
  #                                   value=performance(get(paste0('pred.', j)), measure='auc')@y.values[[1]]))
  #   assign(paste0('false.discoveries.', j), `[<-`(get(paste0('false.discoveries.', j)), k,, 
  #                                                 value=c(get(paste0('pred.', j))@fp[[1]], 
  #                                                         rep(NA, 20000-length(get(paste0('pred.', j))@fp[[1]])))))
  #   assign(paste0('discoveries.', j), `[<-`(get(paste0('discoveries.', j)), k,, 
  #                                           value=c(get(paste0('pred.', j))@n.pos.pred[[1]], 
  #                                                   rep(NA, 20000-length(get(paste0('pred.', j))@n.pos.pred[[1]])))))
  # }
}


# Re-format results to save by experiment type ####

# Set up variables
{
  names.DE <- c('ql.edgeR', 'lr.edgeR', 'et.edgeR', 'noif.DESeq2', 'if.DESeq2', 'voom', 'notrend.DSS',
                 'baySeq', 'mean.zi.MDSeq', 'mean.nozi.MDSeq', 'mean.expHM', 'lmean.expHM', 'mean.lnHM', 'lmean.lnHM',
                'expHMM', 'lnHMM')
  # names.DD <- c('disp.zi.MDSeq', 'disp.nozi.MDSeq', 'disp.expHM', 'ldisp.expHM', 'disp.lnHM', 'ldisp.lnHM',
  # 'expHMM', 'lnHMM')
  # names.DEDD <- c('expHMM', 'lnHMM', names.combined)
}

for (i in c(
  'DE'
  # ,
  # 'DD'
  # ,
  # 'DEDD'
  )) {
  assign(paste0('results.',i), list())
  for (j in c('auc', 'fpr', 'fdr', 'tpr')) {
    assign(paste0('results.',i), `[[<-`(get(paste0('results.',i)), j, data.frame(matrix(nrow=4, ncol=0))))
  }
}

{names.fpr.fdr.tpr <- vector()
  for (i in c(
    names.edgeR, names.DESeq2, names.voom, names.DSS, names.baySeq,
    names.MDSeq, names.expHM[-1], names.lnHM[-1]
    # ,
    # names.combined
    )) {
      names.fpr.fdr.tpr <- c(names.fpr.fdr.tpr, i)
  }
  for (i in names.DSS) {
     names.fpr.fdr.tpr <- c(names.fpr.fdr.tpr, paste0('lfdr.', i))
  }
  for (i in c(names.expHM[1],names.lnHM[1])) {
    for (j in c('5.','thr.','bfdr.')) {
      names.fpr.fdr.tpr <- c(names.fpr.fdr.tpr,paste0(j,i))
    }
  }
}

# Put results into variables
{
  for (i in names.DE) {
   results.DE$auc[[i]] <- get(paste0('auc.',i))
   for (j in grep(i, names.fpr.fdr.tpr, value=T)) {
     results.DE$fpr[[j]] <- get(paste0('fpr.',j))
     results.DE$fdr[[j]] <- get(paste0('fdr.',j))
     results.DE$tpr[[j]] <- get(paste0('tpr.',j))
   }
  }
  # for (i in names.DD) {
  #   results.DD$auc[[i]] <- get(paste0('auc.',i))
  #   for (j in grep(i, names.fpr.fdr.tpr, value=T)) {
  #     results.DD$fpr[[j]] <- get(paste0('fpr.',j))
  #     results.DD$fdr[[j]] <- get(paste0('fdr.',j))
  #     results.DD$tpr[[j]] <- get(paste0('tpr.',j))
  #   }
  # }
  # for (i in names.DEDD) {
  #   results.DEDD$auc[[i]] <- get(paste0('auc.',i))
  #   for (j in grep(i, names.fpr.fdr.tpr, value=T)) {
  #     results.DEDD$fpr[[j]] <- get(paste0('fpr.',j))
  #     results.DEDD$fdr[[j]] <- get(paste0('fdr.',j))
  #     results.DEDD$tpr[[j]] <- get(paste0('tpr.',j))
  #   }
  # }
}


saveRDS(results.DD, file=here("recount data/GTEx/quarantine", "blood_DD_results.rds"))
saveRDS(results.DE, file=here("recount data/GTEx/quarantine", "blood_DE_results.rds"))
saveRDS(results.DEDD, file=here("recount data/GTEx/quarantine", "blood_DEDD_results.rds"))


results.DD <- readRDS(here("Results/GTEx blood artificial DD, DE, DEDD results Dec 2019", 
                           "blood_DD_results.rds"))
results.DE <- readRDS(here("Results/GTEx blood artificial DD, DE, DEDD results Dec 2019", 
                           "blood_DE_results.rds"))
results.DEDD <- readRDS(here("Results/GTEx blood artificial DD, DE, DEDD results Dec 2019", 
                             "blood_DEDD_results.rds"))





