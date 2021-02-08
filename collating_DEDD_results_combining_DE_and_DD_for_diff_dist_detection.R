library(here)
library(compcodeR)
library(ROCR)
library(caret)

# First create variables to store results
names <- c('expHMM', 'lnHMM', 'lnHM_lnHM', 'expHM_expHM', 'MDSeq_MDSeq', 'HM_edgeR')
for (i in names) {
  assign(paste0('fdr.',i), numeric(50))
  assign(paste0('tpr.',i), numeric(50))
  assign(paste0('auc.',i), numeric(50))
  assign(paste0('pauc.',i), numeric(50))
  assign(paste0('discoveries.',i), matrix(nrow=50, ncol=20000))
  assign(paste0('false.discoveries.',i), matrix(nrow=50, ncol=20000))
}

# Then import each run in turn and save results
for (i in 1:50) {
  # import <- paste0('results.DEDD50.',i,'.all.rds')
  import <- paste0('results.DEDD50.',i,'.rds')
  TMM <- readRDS(here('Results/DEDD compcodeR data results July-Aug 2019',import))
  import <- paste0('results.DEDD50.',i,'.DESeqnorm.rds')
  DESeq <- readRDS(here('Results/DEDD compcodeR data results DESeq norm Sept 2019',import))
  DEDD <- factor(as.numeric(TMM$DE==1 | TMM$DD==1), levels=c('1','0'))
  expHMM <- 1 - DESeq$prob.expHM
  lnHMM <- 1 - DESeq$prob.lnHM
  expHM_expHM <- pmin(DESeq$p.mean.expHM, DESeq$p.disp.expHM)
  lnHM_lnHM <- pmin(DESeq$p.mean.lnHM, DESeq$p.disp.lnHM)
  MDSeq_MDSeq <- pmin(DESeq$p.mean.zi.MDSeq, DESeq$p.disp.zi.MDSeq)
  HM_edgeR <- pmin(TMM$p.ql.edgeR, DESeq$p.disp.lnHM)
  thr.expHMM <- sort(DESeq$prob.expHM, decreasing=T)[round(nrow(DESeq$data@count.matrix)*DESeq$prop.expHM)]
  thr.lnHMM <- sort(DESeq$prob.lnHM, decreasing=T)[round(nrow(DESeq$data@count.matrix)*DESeq$prop.lnHM)]
  call.lnHM_lnHM <- factor(as.numeric(p.adjust(lnHM_lnHM, method="BH") < 0.05), levels=c("1", "0"))
  call.expHM_expHM <- factor(as.numeric(p.adjust(expHM_expHM, method="BH") < 0.05), levels=c("1", "0"))
  call.MDSeq_MDSeq <- factor(as.numeric(p.adjust(MDSeq_MDSeq, method="BH") < 0.05), levels=c("1", "0"))
  call.HM_edgeR <- factor(as.numeric(p.adjust(HM_edgeR, method="BH") < 0.05), levels=c("1", "0"))
  call.lnHMM <- factor(as.numeric(lnHMM < thr.lnHMM), levels=c("1", "0"))
  call.expHMM <- factor(as.numeric(expHMM < thr.expHMM), levels=c("1", "0"))
  for (j in names) {
    assign(paste0('fdr.', j), `[<-`(get(paste0('fdr.', j)), i, 
                                    value=1-precision(get(paste0('call.', j)), DEDD)))
    assign(paste0('tpr.', j), `[<-`(get(paste0('tpr.', j)), i, 
                                    value=sensitivity(get(paste0('call.', j)), DEDD)))
    assign(paste0('pred.', j), prediction(1 - get(j), DEDD))
    assign(paste0('auc.', j), `[<-`(get(paste0('auc.', j)), i, 
                                    value=performance(get(paste0('pred.', j)), 
                                                      measure='auc')@y.values[[1]]))
    assign(paste0('pauc.', j), 
           `[<-`(get(paste0('pauc.', j)), i, 
                 value=performance(get(paste0('pred.', j)), 
                                   measure='auc', 
                                   fpr.stop=0.05)@y.values[[1]]))
    assign(paste0('false.discoveries.', j), 
           `[<-`(get(paste0('false.discoveries.', j)), i,,
                 value=c(get(paste0('pred.', j))@fp[[1]], 
                         rep(NA, 20000 - length(get(paste0('pred.', j))@fp[[1]])))))
    assign(paste0('discoveries.', j), 
           `[<-`(get(paste0('discoveries.', j)), i,, 
                 value=c(get(paste0('pred.', j))@n.pos.pred[[1]], 
                         rep(NA, 20000 - length(get(paste0('pred.', j))@n.pos.pred[[1]])))))
  }
}

for (i in names) {
  assign(paste0('mean.discoveries.',i), colMeans(get(paste0('discoveries.',i))))
  assign(paste0('mean.fdr.',i), colMeans(get(paste0('false.discoveries.',i)) / 
                                           get(paste0('discoveries.',i))))
}

# Clean up temporary variables no longer needed
{
  rm(list=ls()[grep('^pred',ls())])
  rm(list=ls()[grep('^discoveries',ls())])
  rm(list=ls()[grep('^false.discoveries',ls())])
  rm(list=c('i','j','import','TMM','DESeq','DEDD', names))
}

# Re-format results to save by experiment type ####

# Set up variables
fdr <- data.frame(matrix(nrow=50, ncol=0))
tpr <- data.frame(matrix(nrow=50, ncol=0))
auc <- data.frame(matrix(nrow=50, ncol=0))
pauc <- data.frame(matrix(nrow=50, ncol=0))
mean.fdr <- data.frame(matrix(nrow=20000, ncol=0))
mean.discoveries <- data.frame(matrix(nrow=20000, ncol=0))

# Put results into variables
for (i in names) {
  fdr[[i]] <- get(paste0('fdr.', i))
  tpr[[i]] <- get(paste0('tpr.', i))
  auc[[i]] <- get(paste0('auc.', i))
  pauc[[i]] <- get(paste0('pauc.', i))
  mean.fdr[[i]] <- get(paste0('mean.fdr.', i))
  mean.discoveries[[i]] <- get(paste0('mean.discoveries.', i))
}

results <- list('fdr' = fdr, 
                'tpr' = tpr, 
                'auc' = auc, 
                'pauc' = pauc, 
                'mean.fdr' = mean.fdr, 
                'mean.discoveries' = mean.discoveries)

saveRDS(results, file=here('Results/Diff dist combining DE, DD predictions Nov 2019', 
                           'diff_dist_results_DEDD50.rds'))

