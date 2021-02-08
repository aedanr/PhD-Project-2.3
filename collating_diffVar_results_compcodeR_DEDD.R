# diffVar results on compcodeR data to add for consistency since included it for GTEx data
# Only TMM, no RLE.

library(here)
library(compcodeR)
library(ROCR)
library(caret)

folder <- "Results/diffVar compcodeR results Dec 2020"

for (samples.per.cond in c(2, 5, 10, 20, 50)) {
  
  auc <- numeric(50)
  pauc <- numeric(50)
  fpr <- numeric(50)
  fdr <- numeric(50)
  tpr <- numeric(50)
  discoveries <- matrix(nrow=50, ncol=20000)
  false.discoveries <- matrix(nrow=50, ncol=20000)
  
  auc_DEDD <- numeric(50)
  pauc_DEDD <- numeric(50)
  fpr_DEDD <- numeric(50)
  fdr_DEDD <- numeric(50)
  tpr_DEDD <- numeric(50)
  discoveries_DEDD <- matrix(nrow=50, ncol=20000)
  false.discoveries_DEDD <- matrix(nrow=50, ncol=20000)
  
  for (i in 1:50) {
    filename <- paste0("results.DEDD", samples.per.cond, ".", i, ".rds")
    res <- readRDS(here(folder, filename))
    DD <- factor(res$DD, levels=c("1", "0"))
    DEDD <- factor(res$DEDD, levels=c("1", "0"))
    
    call <- factor(as.numeric(res$q.diffVar.tmm < 0.05), levels=c("1", "0"))
    fpr[i] <- 1 - specificity(call, DD)
    fdr[i] <- 1 - precision(call, DD)
    tpr[i] <- sensitivity(call, DD)
    pred <- prediction(1 - res$p.diffVar.tmm[which(!is.na(res$p.diffVar.tmm))], DD[which(!is.na(res$p.diffVar.tmm))])
    auc[i] <- performance(pred, measure="auc")@y.values[[1]]
    pauc[i] <-performance(pred, measure="auc", fpr.stop=0.05)@y.values[[1]]
    discoveries[i,] <- c(pred@n.pos.pred[[1]], rep(NA, 20000 - length(pred@n.pos.pred[[1]])))
    false.discoveries[i,] <- c(pred@fp[[1]], rep(NA, 20000 - length(pred@fp[[1]])))
    
    fpr_DEDD[i] <- 1 - specificity(call, DEDD)
    fdr_DEDD[i] <- 1 - precision(call, DEDD)
    tpr_DEDD[i] <- sensitivity(call, DEDD)
    pred_DEDD <- prediction(1 - res$p.diffVar.tmm[which(!is.na(res$p.diffVar.tmm))], DEDD[which(!is.na(res$p.diffVar.tmm))])
    auc_DEDD[i] <- performance(pred_DEDD, measure="auc")@y.values[[1]]
    pauc_DEDD[i] <-performance(pred_DEDD, measure="auc", fpr.stop=0.05)@y.values[[1]]
    discoveries_DEDD[i,] <- c(pred_DEDD@n.pos.pred[[1]], rep(NA, 20000 - length(pred_DEDD@n.pos.pred[[1]])))
    false.discoveries_DEDD[i,] <- c(pred_DEDD@fp[[1]], rep(NA, 20000 - length(pred_DEDD@fp[[1]])))
  }
  
  mean.discoveries <- colMeans(discoveries)
  mean.fdr <- colMeans(false.discoveries / discoveries)

  mean.discoveries_DEDD <- colMeans(discoveries_DEDD)
  mean.fdr_DEDD <- colMeans(false.discoveries_DEDD / discoveries_DEDD)
  
  results.DD <- list(
    auc = auc, 
    pauc = pauc, 
    fpr = fpr, 
    fdr = fdr, 
    tpr = tpr, 
    mean.fdr = mean.fdr, 
    mean.discoveries = mean.discoveries
  )
  
  results.DEDD <- list(
    auc = auc_DEDD, 
    pauc = pauc_DEDD, 
    fpr = fpr_DEDD, 
    fdr = fdr_DEDD, 
    tpr = tpr_DEDD, 
    mean.fdr = mean.fdr_DEDD, 
    mean.discoveries = mean.discoveries_DEDD
  )
  
  saveRDS(results.DD, file=here(folder, paste0("DD.results.", "DEDD", samples.per.cond, ".rds")))
  saveRDS(results.DEDD, file=here(folder, paste0("DEDD.results.", "DEDD", samples.per.cond, ".rds")))
}
