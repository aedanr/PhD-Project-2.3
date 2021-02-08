library(here)
library(ROCR)
library(caret)

folder <- "Results/ExpVarQuant compcodeR results Dec 2020"

for (samples.per.cond in c(2,5,10,20,50)) {

  auc <- numeric(50)
  pauc <- numeric(50)
  fpr <- numeric(50)
  fdr <- numeric(50)
  tpr <- numeric(50)
  discoveries <- matrix(nrow=50, ncol=20000)
  false.discoveries <- matrix(nrow=50, ncol=20000)
  
  for (i in 1:50) {
    
    filename <- paste0("results.ExpVarQuantDD", samples.per.cond, ".", i, ".rds")
    res <- readRDS(here(folder, filename))
    DD <- factor(res$DD_no_NA, levels=c("1", "0"))
    
    pred <- prediction(1 - res$res_no_NA$p.cv, DD)
    auc[i] <- performance(pred, measure="auc")@y.values[[1]]
    pauc[i] <- performance(pred, measure="auc", fpr.stop=0.05)@y.values[[1]]
    call <- factor(as.numeric(res$res_no_NA$padj.cv < 0.05), levels=c("1", "0"))
    fpr[i] <- 1 - specificity(call, DD)
    fdr[i] <- 1 - precision(call, DD)
    tpr[i] <- sensitivity(call, DD)
    discoveries[i, ] <- c(pred@n.pos.pred[[1]], rep(NA, 20000 - length(pred@n.pos.pred[[1]])))
    false.discoveries[i, ] <- c(pred@fp[[1]], rep(NA, 20000 - length(pred@fp[[1]])))

  }
  
  mean.discoveries <- colMeans(discoveries)
  mean.fdr <- colMeans(false.discoveries / discoveries)
  
  results <- list(
    auc = auc, 
    pauc = pauc, 
    fpr = fpr, 
    fdr = fdr, 
    tpr = tpr, 
    mean.fdr = mean.fdr, 
    mean.discoveries = mean.discoveries
  )
  
  saveRDS(results, file=here(folder, paste0("DD.results.", "DD", samples.per.cond, ".rds")))
  
}

