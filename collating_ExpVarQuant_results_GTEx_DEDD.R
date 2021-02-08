library(here)
library(ROCR)
library(caret)

folder <- "Results/ExpVarQuant GTEx results Dec 2020"

for (samples.per.cond in c(2,5,10,20,50)) {

  auc <- numeric(10)
  pauc <- numeric(10)
  fpr <- numeric(10)
  fdr <- numeric(10)
  tpr <- numeric(10)
  discoveries <- matrix(nrow=10, ncol=20000)
  false.discoveries <- matrix(nrow=10, ncol=20000)
  
  for (i in 1:10) {
    
    filename <- paste0("results.ExpVarQuantblood_", samples.per.cond, "_set", i, "_DEDD.rds")
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
  
  saveRDS(results, file=here(folder, paste0("DD.results.blood_", samples.per.cond, "_DEDD.rds")))
  
}


## Quick comparison of results against HMs shows nearly identical performance.
## GAMLSS very slightly better on most measures but not FDR, and HMs have better 
## FDR curves for smaller samples.
folder <- "Results/ExpVarQuant GTEx results Dec 2020"
for (i in c(2,5,10,20,50)) {
  file <- paste0("DD.results.blood_", i, "_DEDD.rds")
  assign(paste0("res.", i), readRDS(here(folder, file)))
}
boxplot(cbind(res.2$auc, res.5$auc, res.10$auc, res.20$auc, res.50$auc))
colMeans(cbind(res.2$auc, res.5$auc, res.10$auc, res.20$auc, res.50$auc))
boxplot(cbind(res.2$pauc, res.5$pauc, res.10$pauc, res.20$pauc, res.50$pauc))
colMeans(cbind(res.2$pauc, res.5$pauc, res.10$pauc, res.20$pauc, res.50$pauc))
boxplot(cbind(res.2$fdr, res.5$fdr, res.10$fdr, res.20$fdr, res.50$fdr))
colMeans(cbind(res.2$fdr, res.5$fdr, res.10$fdr, res.20$fdr, res.50$fdr))
boxplot(cbind(res.2$tpr, res.5$tpr, res.10$tpr, res.20$tpr, res.50$tpr))
colMeans(cbind(res.2$tpr, res.5$tpr, res.10$tpr, res.20$tpr, res.50$tpr))

par(mfrow=c(3,2))
plot(res.2$mean.discoveries, res.2$mean.fdr, type='l', xlim=c(0,2000), ylim=c(0,0.9))
lines(mean.discoveries_blood_2$blood_lnHM.log, mean.fdr_blood_2$blood_lnHM.log, 
      lwd=2, col=col_vector[3], cex.axis=1.5)
abline(h=0.05, col='lightgrey')
plot(res.5$mean.discoveries, res.5$mean.fdr, type='l', xlim=c(0,2000), ylim=c(0,0.9))
lines(mean.discoveries_blood_5$blood_lnHM.log, mean.fdr_blood_5$blood_lnHM.log, 
      lwd=2, col=col_vector[3], cex.axis=1.5)
abline(h=0.05, col='lightgrey')
plot(res.10$mean.discoveries, res.10$mean.fdr, type='l', xlim=c(0,2000), ylim=c(0,0.9))
lines(mean.discoveries_blood_10$blood_lnHM.log, mean.fdr_blood_10$blood_lnHM.log, 
      lwd=2, col=col_vector[3], cex.axis=1.5)
abline(h=0.05, col='lightgrey')
plot(res.20$mean.discoveries, res.20$mean.fdr, type='l', xlim=c(0,2000), ylim=c(0,0.9))
lines(mean.discoveries_blood_20$blood_lnHM.log, mean.fdr_blood_20$blood_lnHM.log, 
      lwd=2, col=col_vector[3], cex.axis=1.5)
abline(h=0.05, col='lightgrey')
plot(res.50$mean.discoveries, res.50$mean.fdr, type='l', xlim=c(0,2000), ylim=c(0,0.9))
lines(mean.discoveries_blood_50$blood_lnHM.log, mean.fdr_blood_50$blood_lnHM.log, 
      lwd=2, col=col_vector[3], cex.axis=1.5)
abline(h=0.05, col='lightgrey')


## GAMLSS fails to converge for some genes
folder <- "Results/ExpVarQuant GTEx results Dec 2020"
for (i in c(2,5,10,20,50)) {
  for (j in 1:10) {
    for (k in c("blood", "muscle")) {
      file <- paste0("results.ExpVarQuant", k, "_", i, "_set", j, "_DEDD.rds")
      assign(paste0(k, "_", i, "_set", j), readRDS(here(folder, file)))
    }
  }
}

for (j in c(2,5,10,20,50)) {
  for (k in c("blood", "muscle")) {
    missing <- numeric(1)
    genes <- numeric(1)
    for (i in 1:10) {
      dat <- get(paste0(k, "_", j, "_set", i))
      missing <- missing + length(dat$DD) - length(dat$DD_no_NA)
      genes <- genes + length(dat$DD)
    }
    assign(paste0("missing_", k, "_", j), missing / 10)
    assign(paste0("genes_", k, "_", j), genes / 10)
  }
}

data.frame(rbind(
  c(missing_blood_2, genes_blood_2, missing_muscle_2, genes_muscle_2), 
  c(missing_blood_5, genes_blood_5, missing_muscle_5, genes_muscle_5), 
  c(missing_blood_10, genes_blood_10, missing_muscle_10, genes_muscle_10), 
  c(missing_blood_20, genes_blood_20, missing_muscle_20, genes_muscle_20), 
  c(missing_blood_50, genes_blood_50, missing_muscle_50, genes_muscle_50)
), 
row.names=c(2, 5, 10, 20, 50))

# average around 3000/18000 missing for 2 samples per group, but very few 
# for larger samples - even for 5, average below 40 missing.

