library(here)
library(RColorBrewer)

# Results initially saved by dataset, but want to assess by experiment type, 
# i.e. detecting differential expression, dispersion or distribution.

## Load data, create colour vector ####
folder <- "Results/compcodeR DE, DD, DEDD results Feb 2020"
# folder <- "Results/temp_quarantine"
for (i in c("DE", "DD", "DEDD")) {
  for (j in c("DEDD2", "DEDD5", "DEDD10", "DEDD20", "DEDD50")) {
    assign(paste0(i, ".results.", j), 
           readRDS(here(folder, paste0(i, ".results.", j, ".rds"))))
  }
}
for (i in c("DE", "DEDD")) {
  for (j in c("DE2", "DE5", "DE10", "DE20", "DE50")) {
    assign(paste0(i, ".results.", j), 
           readRDS(here(folder, paste0(i, ".results.", j, ".rds"))))
  }
}
for (i in c("DD", "DEDD")) {
  for (j in c("DD2", "DD5", "DD10", "DD20", "DD50")) {
    assign(paste0(i, ".results.", j), 
           readRDS(here(folder, paste0(i, ".results.", j, ".rds"))))
  }
}
rm(i,j)

n <- 23
qual_col_pals = brewer.pal.info[brewer.pal.info$category == "qual" & 
                                  brewer.pal.info$colorblind == T,]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# pie(rep(1,n), col=col_vector)
col_vector <- col_vector[c(3,4,9,19)]
legend=c("expHM", "expHM log", "lnHM", "lnHM log")

{
{
  DD.results.DD2$auc <- DD.results.DD2$auc[, grep("HM", names(DD.results.DD2$auc))]
  DD.results.DD2$pauc <- DD.results.DD2$pauc[, grep("HM", names(DD.results.DD2$pauc))]
  DD.results.DD2$fdr <- DD.results.DD2$fdr[, grep("HM", names(DD.results.DD2$fdr))]
  DD.results.DD2$tpr <- DD.results.DD2$tpr[, grep("HM", names(DD.results.DD2$tpr))]
  DD.results.DD2$mean.fdr <- DD.results.DD2$mean.fdr[, grep("HM", names(DD.results.DD2$mean.fdr))]
  DD.results.DD2$mean.discoveries <- DD.results.DD2$mean.discoveries[, grep("HM", names(DD.results.DD2$mean.discoveries))]
  DD.results.DD5$auc <- DD.results.DD5$auc[, grep("HM", names(DD.results.DD5$auc))]
  DD.results.DD5$pauc <- DD.results.DD5$pauc[, grep("HM", names(DD.results.DD5$pauc))]
  DD.results.DD5$fdr <- DD.results.DD5$fdr[, grep("HM", names(DD.results.DD5$fdr))]
  DD.results.DD5$tpr <- DD.results.DD5$tpr[, grep("HM", names(DD.results.DD5$tpr))]
  DD.results.DD5$mean.fdr <- DD.results.DD5$mean.fdr[, grep("HM", names(DD.results.DD5$mean.fdr))]
  DD.results.DD5$mean.discoveries <- DD.results.DD5$mean.discoveries[, grep("HM", names(DD.results.DD5$mean.discoveries))]
  DD.results.DD10$auc <- DD.results.DD10$auc[, grep("HM", names(DD.results.DD10$auc))]
  DD.results.DD10$pauc <- DD.results.DD10$pauc[, grep("HM", names(DD.results.DD10$pauc))]
  DD.results.DD10$fdr <- DD.results.DD10$fdr[, grep("HM", names(DD.results.DD10$fdr))]
  DD.results.DD10$tpr <- DD.results.DD10$tpr[, grep("HM", names(DD.results.DD10$tpr))]
  DD.results.DD10$mean.fdr <- DD.results.DD10$mean.fdr[, grep("HM", names(DD.results.DD10$mean.fdr))]
  DD.results.DD10$mean.discoveries <- DD.results.DD10$mean.discoveries[, grep("HM", names(DD.results.DD10$mean.discoveries))]
  DD.results.DD20$auc <- DD.results.DD20$auc[, grep("HM", names(DD.results.DD20$auc))]
  DD.results.DD20$pauc <- DD.results.DD20$pauc[, grep("HM", names(DD.results.DD20$pauc))]
  DD.results.DD20$fdr <- DD.results.DD20$fdr[, grep("HM", names(DD.results.DD20$fdr))]
  DD.results.DD20$tpr <- DD.results.DD20$tpr[, grep("HM", names(DD.results.DD20$tpr))]
  DD.results.DD20$mean.fdr <- DD.results.DD20$mean.fdr[, grep("HM", names(DD.results.DD20$mean.fdr))]
  DD.results.DD20$mean.discoveries <- DD.results.DD20$mean.discoveries[, grep("HM", names(DD.results.DD20$mean.discoveries))]
  DD.results.DD50$auc <- DD.results.DD50$auc[, grep("HM", names(DD.results.DD50$auc))]
  DD.results.DD50$pauc <- DD.results.DD50$pauc[, grep("HM", names(DD.results.DD50$pauc))]
  DD.results.DD50$fdr <- DD.results.DD50$fdr[, grep("HM", names(DD.results.DD50$fdr))]
  DD.results.DD50$tpr <- DD.results.DD50$tpr[, grep("HM", names(DD.results.DD50$tpr))]
  DD.results.DD50$mean.fdr <- DD.results.DD50$mean.fdr[, grep("HM", names(DD.results.DD50$mean.fdr))]
  DD.results.DD50$mean.discoveries <- DD.results.DD50$mean.discoveries[, grep("HM", names(DD.results.DD50$mean.discoveries))]
  DEDD.results.DD2$auc <- DEDD.results.DD2$auc[, grep("HM", names(DEDD.results.DD2$auc))]
  DEDD.results.DD2$pauc <- DEDD.results.DD2$pauc[, grep("HM", names(DEDD.results.DD2$pauc))]
  DEDD.results.DD2$fdr <- DEDD.results.DD2$fdr[, grep("HM", names(DEDD.results.DD2$fdr))]
  DEDD.results.DD2$tpr <- DEDD.results.DD2$tpr[, grep("HM", names(DEDD.results.DD2$tpr))]
  DEDD.results.DD2$mean.fdr <- DEDD.results.DD2$mean.fdr[, grep("HM", names(DEDD.results.DD2$mean.fdr))]
  DEDD.results.DD2$mean.discoveries <- DEDD.results.DD2$mean.discoveries[, grep("HM", names(DEDD.results.DD2$mean.discoveries))]
  DEDD.results.DD5$auc <- DEDD.results.DD5$auc[, grep("HM", names(DEDD.results.DD5$auc))]
  DEDD.results.DD5$pauc <- DEDD.results.DD5$pauc[, grep("HM", names(DEDD.results.DD5$pauc))]
  DEDD.results.DD5$fdr <- DEDD.results.DD5$fdr[, grep("HM", names(DEDD.results.DD5$fdr))]
  DEDD.results.DD5$tpr <- DEDD.results.DD5$tpr[, grep("HM", names(DEDD.results.DD5$tpr))]
  DEDD.results.DD5$mean.fdr <- DEDD.results.DD5$mean.fdr[, grep("HM", names(DEDD.results.DD5$mean.fdr))]
  DEDD.results.DD5$mean.discoveries <- DEDD.results.DD5$mean.discoveries[, grep("HM", names(DEDD.results.DD5$mean.discoveries))]
  DEDD.results.DD10$auc <- DEDD.results.DD10$auc[, grep("HM", names(DEDD.results.DD10$auc))]
  DEDD.results.DD10$pauc <- DEDD.results.DD10$pauc[, grep("HM", names(DEDD.results.DD10$pauc))]
  DEDD.results.DD10$fdr <- DEDD.results.DD10$fdr[, grep("HM", names(DEDD.results.DD10$fdr))]
  DEDD.results.DD10$tpr <- DEDD.results.DD10$tpr[, grep("HM", names(DEDD.results.DD10$tpr))]
  DEDD.results.DD10$mean.fdr <- DEDD.results.DD10$mean.fdr[, grep("HM", names(DEDD.results.DD10$mean.fdr))]
  DEDD.results.DD10$mean.discoveries <- DEDD.results.DD10$mean.discoveries[, grep("HM", names(DEDD.results.DD10$mean.discoveries))]
  DEDD.results.DD20$auc <- DEDD.results.DD20$auc[, grep("HM", names(DEDD.results.DD20$auc))]
  DEDD.results.DD20$pauc <- DEDD.results.DD20$pauc[, grep("HM", names(DEDD.results.DD20$pauc))]
  DEDD.results.DD20$fdr <- DEDD.results.DD20$fdr[, grep("HM", names(DEDD.results.DD20$fdr))]
  DEDD.results.DD20$tpr <- DEDD.results.DD20$tpr[, grep("HM", names(DEDD.results.DD20$tpr))]
  DEDD.results.DD20$mean.fdr <- DEDD.results.DD20$mean.fdr[, grep("HM", names(DEDD.results.DD20$mean.fdr))]
  DEDD.results.DD20$mean.discoveries <- DEDD.results.DD20$mean.discoveries[, grep("HM", names(DEDD.results.DD20$mean.discoveries))]
  DEDD.results.DD50$auc <- DEDD.results.DD50$auc[, grep("HM", names(DEDD.results.DD50$auc))]
  DEDD.results.DD50$pauc <- DEDD.results.DD50$pauc[, grep("HM", names(DEDD.results.DD50$pauc))]
  DEDD.results.DD50$fdr <- DEDD.results.DD50$fdr[, grep("HM", names(DEDD.results.DD50$fdr))]
  DEDD.results.DD50$tpr <- DEDD.results.DD50$tpr[, grep("HM", names(DEDD.results.DD50$tpr))]
  DEDD.results.DD50$mean.fdr <- DEDD.results.DD50$mean.fdr[, grep("HM", names(DEDD.results.DD50$mean.fdr))]
  DEDD.results.DD50$mean.discoveries <- DEDD.results.DD50$mean.discoveries[, grep("HM", names(DEDD.results.DD50$mean.discoveries))]
}

{
  DE.results.DE2$auc <- DE.results.DE2$auc[, grep("HM", names(DE.results.DE2$auc))]
  DE.results.DE2$pauc <- DE.results.DE2$pauc[, grep("HM", names(DE.results.DE2$pauc))]
  DE.results.DE2$fdr <- DE.results.DE2$fdr[, grep("HM", names(DE.results.DE2$fdr))]
  DE.results.DE2$tpr <- DE.results.DE2$tpr[, grep("HM", names(DE.results.DE2$tpr))]
  DE.results.DE2$mean.fdr <- DE.results.DE2$mean.fdr[, grep("HM", names(DE.results.DE2$mean.fdr))]
  DE.results.DE2$mean.discoveries <- DE.results.DE2$mean.discoveries[, grep("HM", names(DE.results.DE2$mean.discoveries))]
  DE.results.DE5$auc <- DE.results.DE5$auc[, grep("HM", names(DE.results.DE5$auc))]
  DE.results.DE5$pauc <- DE.results.DE5$pauc[, grep("HM", names(DE.results.DE5$pauc))]
  DE.results.DE5$fdr <- DE.results.DE5$fdr[, grep("HM", names(DE.results.DE5$fdr))]
  DE.results.DE5$tpr <- DE.results.DE5$tpr[, grep("HM", names(DE.results.DE5$tpr))]
  DE.results.DE5$mean.fdr <- DE.results.DE5$mean.fdr[, grep("HM", names(DE.results.DE5$mean.fdr))]
  DE.results.DE5$mean.discoveries <- DE.results.DE5$mean.discoveries[, grep("HM", names(DE.results.DE5$mean.discoveries))]
  DE.results.DE10$auc <- DE.results.DE10$auc[, grep("HM", names(DE.results.DE10$auc))]
  DE.results.DE10$pauc <- DE.results.DE10$pauc[, grep("HM", names(DE.results.DE10$pauc))]
  DE.results.DE10$fdr <- DE.results.DE10$fdr[, grep("HM", names(DE.results.DE10$fdr))]
  DE.results.DE10$tpr <- DE.results.DE10$tpr[, grep("HM", names(DE.results.DE10$tpr))]
  DE.results.DE10$mean.fdr <- DE.results.DE10$mean.fdr[, grep("HM", names(DE.results.DE10$mean.fdr))]
  DE.results.DE10$mean.discoveries <- DE.results.DE10$mean.discoveries[, grep("HM", names(DE.results.DE10$mean.discoveries))]
  DE.results.DE20$auc <- DE.results.DE20$auc[, grep("HM", names(DE.results.DE20$auc))]
  DE.results.DE20$pauc <- DE.results.DE20$pauc[, grep("HM", names(DE.results.DE20$pauc))]
  DE.results.DE20$fdr <- DE.results.DE20$fdr[, grep("HM", names(DE.results.DE20$fdr))]
  DE.results.DE20$tpr <- DE.results.DE20$tpr[, grep("HM", names(DE.results.DE20$tpr))]
  DE.results.DE20$mean.fdr <- DE.results.DE20$mean.fdr[, grep("HM", names(DE.results.DE20$mean.fdr))]
  DE.results.DE20$mean.discoveries <- DE.results.DE20$mean.discoveries[, grep("HM", names(DE.results.DE20$mean.discoveries))]
  DE.results.DE50$auc <- DE.results.DE50$auc[, grep("HM", names(DE.results.DE50$auc))]
  DE.results.DE50$pauc <- DE.results.DE50$pauc[, grep("HM", names(DE.results.DE50$pauc))]
  DE.results.DE50$fdr <- DE.results.DE50$fdr[, grep("HM", names(DE.results.DE50$fdr))]
  DE.results.DE50$tpr <- DE.results.DE50$tpr[, grep("HM", names(DE.results.DE50$tpr))]
  DE.results.DE50$mean.fdr <- DE.results.DE50$mean.fdr[, grep("HM", names(DE.results.DE50$mean.fdr))]
  DE.results.DE50$mean.discoveries <- DE.results.DE50$mean.discoveries[, grep("HM", names(DE.results.DE50$mean.discoveries))]
  DEDD.results.DE2$auc <- DEDD.results.DE2$auc[, grep("HM", names(DEDD.results.DE2$auc))]
  DEDD.results.DE2$pauc <- DEDD.results.DE2$pauc[, grep("HM", names(DEDD.results.DE2$pauc))]
  DEDD.results.DE2$fdr <- DEDD.results.DE2$fdr[, grep("HM", names(DEDD.results.DE2$fdr))]
  DEDD.results.DE2$tpr <- DEDD.results.DE2$tpr[, grep("HM", names(DEDD.results.DE2$tpr))]
  DEDD.results.DE2$mean.fdr <- DEDD.results.DE2$mean.fdr[, grep("HM", names(DEDD.results.DE2$mean.fdr))]
  DEDD.results.DE2$mean.discoveries <- DEDD.results.DE2$mean.discoveries[, grep("HM", names(DEDD.results.DE2$mean.discoveries))]
  DEDD.results.DE5$auc <- DEDD.results.DE5$auc[, grep("HM", names(DEDD.results.DE5$auc))]
  DEDD.results.DE5$pauc <- DEDD.results.DE5$pauc[, grep("HM", names(DEDD.results.DE5$pauc))]
  DEDD.results.DE5$fdr <- DEDD.results.DE5$fdr[, grep("HM", names(DEDD.results.DE5$fdr))]
  DEDD.results.DE5$tpr <- DEDD.results.DE5$tpr[, grep("HM", names(DEDD.results.DE5$tpr))]
  DEDD.results.DE5$mean.fdr <- DEDD.results.DE5$mean.fdr[, grep("HM", names(DEDD.results.DE5$mean.fdr))]
  DEDD.results.DE5$mean.discoveries <- DEDD.results.DE5$mean.discoveries[, grep("HM", names(DEDD.results.DE5$mean.discoveries))]
  DEDD.results.DE10$auc <- DEDD.results.DE10$auc[, grep("HM", names(DEDD.results.DE10$auc))]
  DEDD.results.DE10$pauc <- DEDD.results.DE10$pauc[, grep("HM", names(DEDD.results.DE10$pauc))]
  DEDD.results.DE10$fdr <- DEDD.results.DE10$fdr[, grep("HM", names(DEDD.results.DE10$fdr))]
  DEDD.results.DE10$tpr <- DEDD.results.DE10$tpr[, grep("HM", names(DEDD.results.DE10$tpr))]
  DEDD.results.DE10$mean.fdr <- DEDD.results.DE10$mean.fdr[, grep("HM", names(DEDD.results.DE10$mean.fdr))]
  DEDD.results.DE10$mean.discoveries <- DEDD.results.DE10$mean.discoveries[, grep("HM", names(DEDD.results.DE10$mean.discoveries))]
  DEDD.results.DE20$auc <- DEDD.results.DE20$auc[, grep("HM", names(DEDD.results.DE20$auc))]
  DEDD.results.DE20$pauc <- DEDD.results.DE20$pauc[, grep("HM", names(DEDD.results.DE20$pauc))]
  DEDD.results.DE20$fdr <- DEDD.results.DE20$fdr[, grep("HM", names(DEDD.results.DE20$fdr))]
  DEDD.results.DE20$tpr <- DEDD.results.DE20$tpr[, grep("HM", names(DEDD.results.DE20$tpr))]
  DEDD.results.DE20$mean.fdr <- DEDD.results.DE20$mean.fdr[, grep("HM", names(DEDD.results.DE20$mean.fdr))]
  DEDD.results.DE20$mean.discoveries <- DEDD.results.DE20$mean.discoveries[, grep("HM", names(DEDD.results.DE20$mean.discoveries))]
  DEDD.results.DE50$auc <- DEDD.results.DE50$auc[, grep("HM", names(DEDD.results.DE50$auc))]
  DEDD.results.DE50$pauc <- DEDD.results.DE50$pauc[, grep("HM", names(DEDD.results.DE50$pauc))]
  DEDD.results.DE50$fdr <- DEDD.results.DE50$fdr[, grep("HM", names(DEDD.results.DE50$fdr))]
  DEDD.results.DE50$tpr <- DEDD.results.DE50$tpr[, grep("HM", names(DEDD.results.DE50$tpr))]
  DEDD.results.DE50$mean.fdr <- DEDD.results.DE50$mean.fdr[, grep("HM", names(DEDD.results.DE50$mean.fdr))]
  DEDD.results.DE50$mean.discoveries <- DEDD.results.DE50$mean.discoveries[, grep("HM", names(DEDD.results.DE50$mean.discoveries))]
}

{
  DD.results.DEDD2$auc <- DD.results.DEDD2$auc[, grep("HM", names(DD.results.DEDD2$auc))]
  DD.results.DEDD2$pauc <- DD.results.DEDD2$pauc[, grep("HM", names(DD.results.DEDD2$pauc))]
  DD.results.DEDD2$fdr <- DD.results.DEDD2$fdr[, grep("HM", names(DD.results.DEDD2$fdr))]
  DD.results.DEDD2$tpr <- DD.results.DEDD2$tpr[, grep("HM", names(DD.results.DEDD2$tpr))]
  DD.results.DEDD2$mean.fdr <- DD.results.DEDD2$mean.fdr[, grep("HM", names(DD.results.DEDD2$mean.fdr))]
  DD.results.DEDD2$mean.discoveries <- DD.results.DEDD2$mean.discoveries[, grep("HM", names(DD.results.DEDD2$mean.discoveries))]
  DD.results.DEDD5$auc <- DD.results.DEDD5$auc[, grep("HM", names(DD.results.DEDD5$auc))]
  DD.results.DEDD5$pauc <- DD.results.DEDD5$pauc[, grep("HM", names(DD.results.DEDD5$pauc))]
  DD.results.DEDD5$fdr <- DD.results.DEDD5$fdr[, grep("HM", names(DD.results.DEDD5$fdr))]
  DD.results.DEDD5$tpr <- DD.results.DEDD5$tpr[, grep("HM", names(DD.results.DEDD5$tpr))]
  DD.results.DEDD5$mean.fdr <- DD.results.DEDD5$mean.fdr[, grep("HM", names(DD.results.DEDD5$mean.fdr))]
  DD.results.DEDD5$mean.discoveries <- DD.results.DEDD5$mean.discoveries[, grep("HM", names(DD.results.DEDD5$mean.discoveries))]
  DD.results.DEDD10$auc <- DD.results.DEDD10$auc[, grep("HM", names(DD.results.DEDD10$auc))]
  DD.results.DEDD10$pauc <- DD.results.DEDD10$pauc[, grep("HM", names(DD.results.DEDD10$pauc))]
  DD.results.DEDD10$fdr <- DD.results.DEDD10$fdr[, grep("HM", names(DD.results.DEDD10$fdr))]
  DD.results.DEDD10$tpr <- DD.results.DEDD10$tpr[, grep("HM", names(DD.results.DEDD10$tpr))]
  DD.results.DEDD10$mean.fdr <- DD.results.DEDD10$mean.fdr[, grep("HM", names(DD.results.DEDD10$mean.fdr))]
  DD.results.DEDD10$mean.discoveries <- DD.results.DEDD10$mean.discoveries[, grep("HM", names(DD.results.DEDD10$mean.discoveries))]
  DD.results.DEDD20$auc <- DD.results.DEDD20$auc[, grep("HM", names(DD.results.DEDD20$auc))]
  DD.results.DEDD20$pauc <- DD.results.DEDD20$pauc[, grep("HM", names(DD.results.DEDD20$pauc))]
  DD.results.DEDD20$fdr <- DD.results.DEDD20$fdr[, grep("HM", names(DD.results.DEDD20$fdr))]
  DD.results.DEDD20$tpr <- DD.results.DEDD20$tpr[, grep("HM", names(DD.results.DEDD20$tpr))]
  DD.results.DEDD20$mean.fdr <- DD.results.DEDD20$mean.fdr[, grep("HM", names(DD.results.DEDD20$mean.fdr))]
  DD.results.DEDD20$mean.discoveries <- DD.results.DEDD20$mean.discoveries[, grep("HM", names(DD.results.DEDD20$mean.discoveries))]
  DD.results.DEDD50$auc <- DD.results.DEDD50$auc[, grep("HM", names(DD.results.DEDD50$auc))]
  DD.results.DEDD50$pauc <- DD.results.DEDD50$pauc[, grep("HM", names(DD.results.DEDD50$pauc))]
  DD.results.DEDD50$fdr <- DD.results.DEDD50$fdr[, grep("HM", names(DD.results.DEDD50$fdr))]
  DD.results.DEDD50$tpr <- DD.results.DEDD50$tpr[, grep("HM", names(DD.results.DEDD50$tpr))]
  DD.results.DEDD50$mean.fdr <- DD.results.DEDD50$mean.fdr[, grep("HM", names(DD.results.DEDD50$mean.fdr))]
  DD.results.DEDD50$mean.discoveries <- DD.results.DEDD50$mean.discoveries[, grep("HM", names(DD.results.DEDD50$mean.discoveries))]
  DE.results.DEDD2$auc <- DE.results.DEDD2$auc[, grep("HM", names(DE.results.DEDD2$auc))]
  DE.results.DEDD2$pauc <- DE.results.DEDD2$pauc[, grep("HM", names(DE.results.DEDD2$pauc))]
  DE.results.DEDD2$fdr <- DE.results.DEDD2$fdr[, grep("HM", names(DE.results.DEDD2$fdr))]
  DE.results.DEDD2$tpr <- DE.results.DEDD2$tpr[, grep("HM", names(DE.results.DEDD2$tpr))]
  DE.results.DEDD2$mean.fdr <- DE.results.DEDD2$mean.fdr[, grep("HM", names(DE.results.DEDD2$mean.fdr))]
  DE.results.DEDD2$mean.discoveries <- DE.results.DEDD2$mean.discoveries[, grep("HM", names(DE.results.DEDD2$mean.discoveries))]
  DE.results.DEDD5$auc <- DE.results.DEDD5$auc[, grep("HM", names(DE.results.DEDD5$auc))]
  DE.results.DEDD5$pauc <- DE.results.DEDD5$pauc[, grep("HM", names(DE.results.DEDD5$pauc))]
  DE.results.DEDD5$fdr <- DE.results.DEDD5$fdr[, grep("HM", names(DE.results.DEDD5$fdr))]
  DE.results.DEDD5$tpr <- DE.results.DEDD5$tpr[, grep("HM", names(DE.results.DEDD5$tpr))]
  DE.results.DEDD5$mean.fdr <- DE.results.DEDD5$mean.fdr[, grep("HM", names(DE.results.DEDD5$mean.fdr))]
  DE.results.DEDD5$mean.discoveries <- DE.results.DEDD5$mean.discoveries[, grep("HM", names(DE.results.DEDD5$mean.discoveries))]
  DE.results.DEDD10$auc <- DE.results.DEDD10$auc[, grep("HM", names(DE.results.DEDD10$auc))]
  DE.results.DEDD10$pauc <- DE.results.DEDD10$pauc[, grep("HM", names(DE.results.DEDD10$pauc))]
  DE.results.DEDD10$fdr <- DE.results.DEDD10$fdr[, grep("HM", names(DE.results.DEDD10$fdr))]
  DE.results.DEDD10$tpr <- DE.results.DEDD10$tpr[, grep("HM", names(DE.results.DEDD10$tpr))]
  DE.results.DEDD10$mean.fdr <- DE.results.DEDD10$mean.fdr[, grep("HM", names(DE.results.DEDD10$mean.fdr))]
  DE.results.DEDD10$mean.discoveries <- DE.results.DEDD10$mean.discoveries[, grep("HM", names(DE.results.DEDD10$mean.discoveries))]
  DE.results.DEDD20$auc <- DE.results.DEDD20$auc[, grep("HM", names(DE.results.DEDD20$auc))]
  DE.results.DEDD20$pauc <- DE.results.DEDD20$pauc[, grep("HM", names(DE.results.DEDD20$pauc))]
  DE.results.DEDD20$fdr <- DE.results.DEDD20$fdr[, grep("HM", names(DE.results.DEDD20$fdr))]
  DE.results.DEDD20$tpr <- DE.results.DEDD20$tpr[, grep("HM", names(DE.results.DEDD20$tpr))]
  DE.results.DEDD20$mean.fdr <- DE.results.DEDD20$mean.fdr[, grep("HM", names(DE.results.DEDD20$mean.fdr))]
  DE.results.DEDD20$mean.discoveries <- DE.results.DEDD20$mean.discoveries[, grep("HM", names(DE.results.DEDD20$mean.discoveries))]
  DE.results.DEDD50$auc <- DE.results.DEDD50$auc[, grep("HM", names(DE.results.DEDD50$auc))]
  DE.results.DEDD50$pauc <- DE.results.DEDD50$pauc[, grep("HM", names(DE.results.DEDD50$pauc))]
  DE.results.DEDD50$fdr <- DE.results.DEDD50$fdr[, grep("HM", names(DE.results.DEDD50$fdr))]
  DE.results.DEDD50$tpr <- DE.results.DEDD50$tpr[, grep("HM", names(DE.results.DEDD50$tpr))]
  DE.results.DEDD50$mean.fdr <- DE.results.DEDD50$mean.fdr[, grep("HM", names(DE.results.DEDD50$mean.fdr))]
  DE.results.DEDD50$mean.discoveries <- DE.results.DEDD50$mean.discoveries[, grep("HM", names(DE.results.DEDD50$mean.discoveries))]
  DEDD.results.DEDD2$auc <- DEDD.results.DEDD2$auc[, grep("HMM", names(DEDD.results.DEDD2$auc))]
  DEDD.results.DEDD2$pauc <- DEDD.results.DEDD2$pauc[, grep("HMM", names(DEDD.results.DEDD2$pauc))]
  DEDD.results.DEDD2$fdr <- DEDD.results.DEDD2$fdr[, grep("HMM", names(DEDD.results.DEDD2$fdr))]
  DEDD.results.DEDD2$tpr <- DEDD.results.DEDD2$tpr[, grep("HMM", names(DEDD.results.DEDD2$tpr))]
  DEDD.results.DEDD2$mean.fdr <- DEDD.results.DEDD2$mean.fdr[, grep("HMM", names(DEDD.results.DEDD2$mean.fdr))]
  DEDD.results.DEDD2$mean.discoveries <- DEDD.results.DEDD2$mean.discoveries[
    , grep("HMM", names(DEDD.results.DEDD2$mean.discoveries))
    ]
  DEDD.results.DEDD5$auc <- DEDD.results.DEDD5$auc[, grep("HMM", names(DEDD.results.DEDD5$auc))]
  DEDD.results.DEDD5$pauc <- DEDD.results.DEDD5$pauc[, grep("HMM", names(DEDD.results.DEDD5$pauc))]
  DEDD.results.DEDD5$fdr <- DEDD.results.DEDD5$fdr[, grep("HMM", names(DEDD.results.DEDD5$fdr))]
  DEDD.results.DEDD5$tpr <- DEDD.results.DEDD5$tpr[, grep("HMM", names(DEDD.results.DEDD5$tpr))]
  DEDD.results.DEDD5$mean.fdr <- DEDD.results.DEDD5$mean.fdr[, grep("HMM", names(DEDD.results.DEDD5$mean.fdr))]
  DEDD.results.DEDD5$mean.discoveries <- DEDD.results.DEDD5$mean.discoveries[
    , grep("HMM", names(DEDD.results.DEDD5$mean.discoveries))
    ]
  DEDD.results.DEDD10$auc <- DEDD.results.DEDD10$auc[, grep("HMM", names(DEDD.results.DEDD10$auc))]
  DEDD.results.DEDD10$pauc <- DEDD.results.DEDD10$pauc[, grep("HMM", names(DEDD.results.DEDD10$pauc))]
  DEDD.results.DEDD10$fdr <- DEDD.results.DEDD10$fdr[, grep("HMM", names(DEDD.results.DEDD10$fdr))]
  DEDD.results.DEDD10$tpr <- DEDD.results.DEDD10$tpr[, grep("HMM", names(DEDD.results.DEDD10$tpr))]
  DEDD.results.DEDD10$mean.fdr <- DEDD.results.DEDD10$mean.fdr[, grep("HMM", names(DEDD.results.DEDD10$mean.fdr))]
  DEDD.results.DEDD10$mean.discoveries <- DEDD.results.DEDD10$mean.discoveries[
    , grep("HMM", names(DEDD.results.DEDD10$mean.discoveries))
    ]
  DEDD.results.DEDD20$auc <- DEDD.results.DEDD20$auc[, grep("HMM", names(DEDD.results.DEDD20$auc))]
  DEDD.results.DEDD20$pauc <- DEDD.results.DEDD20$pauc[, grep("HMM", names(DEDD.results.DEDD20$pauc))]
  DEDD.results.DEDD20$fdr <- DEDD.results.DEDD20$fdr[, grep("HMM", names(DEDD.results.DEDD20$fdr))]
  DEDD.results.DEDD20$tpr <- DEDD.results.DEDD20$tpr[, grep("HMM", names(DEDD.results.DEDD20$tpr))]
  DEDD.results.DEDD20$mean.fdr <- DEDD.results.DEDD20$mean.fdr[, grep("HMM", names(DEDD.results.DEDD20$mean.fdr))]
  DEDD.results.DEDD20$mean.discoveries <- DEDD.results.DEDD20$mean.discoveries[
    , grep("HMM", names(DEDD.results.DEDD20$mean.discoveries))
    ]
  DEDD.results.DEDD50$auc <- DEDD.results.DEDD50$auc[, grep("HMM", names(DEDD.results.DEDD50$auc))]
  DEDD.results.DEDD50$pauc <- DEDD.results.DEDD50$pauc[, grep("HMM", names(DEDD.results.DEDD50$pauc))]
  DEDD.results.DEDD50$fdr <- DEDD.results.DEDD50$fdr[, grep("HMM", names(DEDD.results.DEDD50$fdr))]
  DEDD.results.DEDD50$tpr <- DEDD.results.DEDD50$tpr[, grep("HMM", names(DEDD.results.DEDD50$tpr))]
  DEDD.results.DEDD50$mean.fdr <- DEDD.results.DEDD50$mean.fdr[, grep("HMM", names(DEDD.results.DEDD50$mean.fdr))]
  DEDD.results.DEDD50$mean.discoveries <- DEDD.results.DEDD50$mean.discoveries[
    , grep("HMM", names(DEDD.results.DEDD50$mean.discoveries))
    ]
}

{
  DD.results.DD2$auc <- DD.results.DD2$auc[, grep("tmm", names(DD.results.DD2$auc))]
  DD.results.DD2$pauc <- DD.results.DD2$pauc[, grep("tmm", names(DD.results.DD2$pauc))]
  DD.results.DD2$fdr <- DD.results.DD2$fdr[, grep("tmm", names(DD.results.DD2$fdr))]
  DD.results.DD2$tpr <- DD.results.DD2$tpr[, grep("tmm", names(DD.results.DD2$tpr))]
  DD.results.DD2$mean.fdr <- DD.results.DD2$mean.fdr[, grep("tmm", names(DD.results.DD2$mean.fdr))]
  DD.results.DD2$mean.discoveries <- DD.results.DD2$mean.discoveries[, grep("tmm", names(DD.results.DD2$mean.discoveries))]
  DD.results.DD5$auc <- DD.results.DD5$auc[, grep("tmm", names(DD.results.DD5$auc))]
  DD.results.DD5$pauc <- DD.results.DD5$pauc[, grep("tmm", names(DD.results.DD5$pauc))]
  DD.results.DD5$fdr <- DD.results.DD5$fdr[, grep("tmm", names(DD.results.DD5$fdr))]
  DD.results.DD5$tpr <- DD.results.DD5$tpr[, grep("tmm", names(DD.results.DD5$tpr))]
  DD.results.DD5$mean.fdr <- DD.results.DD5$mean.fdr[, grep("tmm", names(DD.results.DD5$mean.fdr))]
  DD.results.DD5$mean.discoveries <- DD.results.DD5$mean.discoveries[, grep("tmm", names(DD.results.DD5$mean.discoveries))]
  DD.results.DD10$auc <- DD.results.DD10$auc[, grep("tmm", names(DD.results.DD10$auc))]
  DD.results.DD10$pauc <- DD.results.DD10$pauc[, grep("tmm", names(DD.results.DD10$pauc))]
  DD.results.DD10$fdr <- DD.results.DD10$fdr[, grep("tmm", names(DD.results.DD10$fdr))]
  DD.results.DD10$tpr <- DD.results.DD10$tpr[, grep("tmm", names(DD.results.DD10$tpr))]
  DD.results.DD10$mean.fdr <- DD.results.DD10$mean.fdr[, grep("tmm", names(DD.results.DD10$mean.fdr))]
  DD.results.DD10$mean.discoveries <- DD.results.DD10$mean.discoveries[, grep("tmm", names(DD.results.DD10$mean.discoveries))]
  DD.results.DD20$auc <- DD.results.DD20$auc[, grep("tmm", names(DD.results.DD20$auc))]
  DD.results.DD20$pauc <- DD.results.DD20$pauc[, grep("tmm", names(DD.results.DD20$pauc))]
  DD.results.DD20$fdr <- DD.results.DD20$fdr[, grep("tmm", names(DD.results.DD20$fdr))]
  DD.results.DD20$tpr <- DD.results.DD20$tpr[, grep("tmm", names(DD.results.DD20$tpr))]
  DD.results.DD20$mean.fdr <- DD.results.DD20$mean.fdr[, grep("tmm", names(DD.results.DD20$mean.fdr))]
  DD.results.DD20$mean.discoveries <- DD.results.DD20$mean.discoveries[, grep("tmm", names(DD.results.DD20$mean.discoveries))]
  DD.results.DD50$auc <- DD.results.DD50$auc[, grep("tmm", names(DD.results.DD50$auc))]
  DD.results.DD50$pauc <- DD.results.DD50$pauc[, grep("tmm", names(DD.results.DD50$pauc))]
  DD.results.DD50$fdr <- DD.results.DD50$fdr[, grep("tmm", names(DD.results.DD50$fdr))]
  DD.results.DD50$tpr <- DD.results.DD50$tpr[, grep("tmm", names(DD.results.DD50$tpr))]
  DD.results.DD50$mean.fdr <- DD.results.DD50$mean.fdr[, grep("tmm", names(DD.results.DD50$mean.fdr))]
  DD.results.DD50$mean.discoveries <- DD.results.DD50$mean.discoveries[, grep("tmm", names(DD.results.DD50$mean.discoveries))]
  DEDD.results.DD2$auc <- DEDD.results.DD2$auc[, grep("tmm", names(DEDD.results.DD2$auc))]
  DEDD.results.DD2$pauc <- DEDD.results.DD2$pauc[, grep("tmm", names(DEDD.results.DD2$pauc))]
  DEDD.results.DD2$fdr <- DEDD.results.DD2$fdr[, grep("tmm", names(DEDD.results.DD2$fdr))]
  DEDD.results.DD2$tpr <- DEDD.results.DD2$tpr[, grep("tmm", names(DEDD.results.DD2$tpr))]
  DEDD.results.DD2$mean.fdr <- DEDD.results.DD2$mean.fdr[, grep("tmm", names(DEDD.results.DD2$mean.fdr))]
  DEDD.results.DD2$mean.discoveries <- DEDD.results.DD2$mean.discoveries[, grep("tmm", names(DEDD.results.DD2$mean.discoveries))]
  DEDD.results.DD5$auc <- DEDD.results.DD5$auc[, grep("tmm", names(DEDD.results.DD5$auc))]
  DEDD.results.DD5$pauc <- DEDD.results.DD5$pauc[, grep("tmm", names(DEDD.results.DD5$pauc))]
  DEDD.results.DD5$fdr <- DEDD.results.DD5$fdr[, grep("tmm", names(DEDD.results.DD5$fdr))]
  DEDD.results.DD5$tpr <- DEDD.results.DD5$tpr[, grep("tmm", names(DEDD.results.DD5$tpr))]
  DEDD.results.DD5$mean.fdr <- DEDD.results.DD5$mean.fdr[, grep("tmm", names(DEDD.results.DD5$mean.fdr))]
  DEDD.results.DD5$mean.discoveries <- DEDD.results.DD5$mean.discoveries[, grep("tmm", names(DEDD.results.DD5$mean.discoveries))]
  DEDD.results.DD10$auc <- DEDD.results.DD10$auc[, grep("tmm", names(DEDD.results.DD10$auc))]
  DEDD.results.DD10$pauc <- DEDD.results.DD10$pauc[, grep("tmm", names(DEDD.results.DD10$pauc))]
  DEDD.results.DD10$fdr <- DEDD.results.DD10$fdr[, grep("tmm", names(DEDD.results.DD10$fdr))]
  DEDD.results.DD10$tpr <- DEDD.results.DD10$tpr[, grep("tmm", names(DEDD.results.DD10$tpr))]
  DEDD.results.DD10$mean.fdr <- DEDD.results.DD10$mean.fdr[, grep("tmm", names(DEDD.results.DD10$mean.fdr))]
  DEDD.results.DD10$mean.discoveries <- DEDD.results.DD10$mean.discoveries[, grep("tmm", names(DEDD.results.DD10$mean.discoveries))]
  DEDD.results.DD20$auc <- DEDD.results.DD20$auc[, grep("tmm", names(DEDD.results.DD20$auc))]
  DEDD.results.DD20$pauc <- DEDD.results.DD20$pauc[, grep("tmm", names(DEDD.results.DD20$pauc))]
  DEDD.results.DD20$fdr <- DEDD.results.DD20$fdr[, grep("tmm", names(DEDD.results.DD20$fdr))]
  DEDD.results.DD20$tpr <- DEDD.results.DD20$tpr[, grep("tmm", names(DEDD.results.DD20$tpr))]
  DEDD.results.DD20$mean.fdr <- DEDD.results.DD20$mean.fdr[, grep("tmm", names(DEDD.results.DD20$mean.fdr))]
  DEDD.results.DD20$mean.discoveries <- DEDD.results.DD20$mean.discoveries[, grep("tmm", names(DEDD.results.DD20$mean.discoveries))]
  DEDD.results.DD50$auc <- DEDD.results.DD50$auc[, grep("tmm", names(DEDD.results.DD50$auc))]
  DEDD.results.DD50$pauc <- DEDD.results.DD50$pauc[, grep("tmm", names(DEDD.results.DD50$pauc))]
  DEDD.results.DD50$fdr <- DEDD.results.DD50$fdr[, grep("tmm", names(DEDD.results.DD50$fdr))]
  DEDD.results.DD50$tpr <- DEDD.results.DD50$tpr[, grep("tmm", names(DEDD.results.DD50$tpr))]
  DEDD.results.DD50$mean.fdr <- DEDD.results.DD50$mean.fdr[, grep("tmm", names(DEDD.results.DD50$mean.fdr))]
  DEDD.results.DD50$mean.discoveries <- DEDD.results.DD50$mean.discoveries[, grep("tmm", names(DEDD.results.DD50$mean.discoveries))]
}

{
  DE.results.DE2$auc <- DE.results.DE2$auc[, grep("tmm", names(DE.results.DE2$auc))]
  DE.results.DE2$pauc <- DE.results.DE2$pauc[, grep("tmm", names(DE.results.DE2$pauc))]
  DE.results.DE2$fdr <- DE.results.DE2$fdr[, grep("tmm", names(DE.results.DE2$fdr))]
  DE.results.DE2$tpr <- DE.results.DE2$tpr[, grep("tmm", names(DE.results.DE2$tpr))]
  DE.results.DE2$mean.fdr <- DE.results.DE2$mean.fdr[, grep("tmm", names(DE.results.DE2$mean.fdr))]
  DE.results.DE2$mean.discoveries <- DE.results.DE2$mean.discoveries[, grep("tmm", names(DE.results.DE2$mean.discoveries))]
  DE.results.DE5$auc <- DE.results.DE5$auc[, grep("tmm", names(DE.results.DE5$auc))]
  DE.results.DE5$pauc <- DE.results.DE5$pauc[, grep("tmm", names(DE.results.DE5$pauc))]
  DE.results.DE5$fdr <- DE.results.DE5$fdr[, grep("tmm", names(DE.results.DE5$fdr))]
  DE.results.DE5$tpr <- DE.results.DE5$tpr[, grep("tmm", names(DE.results.DE5$tpr))]
  DE.results.DE5$mean.fdr <- DE.results.DE5$mean.fdr[, grep("tmm", names(DE.results.DE5$mean.fdr))]
  DE.results.DE5$mean.discoveries <- DE.results.DE5$mean.discoveries[, grep("tmm", names(DE.results.DE5$mean.discoveries))]
  DE.results.DE10$auc <- DE.results.DE10$auc[, grep("tmm", names(DE.results.DE10$auc))]
  DE.results.DE10$pauc <- DE.results.DE10$pauc[, grep("tmm", names(DE.results.DE10$pauc))]
  DE.results.DE10$fdr <- DE.results.DE10$fdr[, grep("tmm", names(DE.results.DE10$fdr))]
  DE.results.DE10$tpr <- DE.results.DE10$tpr[, grep("tmm", names(DE.results.DE10$tpr))]
  DE.results.DE10$mean.fdr <- DE.results.DE10$mean.fdr[, grep("tmm", names(DE.results.DE10$mean.fdr))]
  DE.results.DE10$mean.discoveries <- DE.results.DE10$mean.discoveries[, grep("tmm", names(DE.results.DE10$mean.discoveries))]
  DE.results.DE20$auc <- DE.results.DE20$auc[, grep("tmm", names(DE.results.DE20$auc))]
  DE.results.DE20$pauc <- DE.results.DE20$pauc[, grep("tmm", names(DE.results.DE20$pauc))]
  DE.results.DE20$fdr <- DE.results.DE20$fdr[, grep("tmm", names(DE.results.DE20$fdr))]
  DE.results.DE20$tpr <- DE.results.DE20$tpr[, grep("tmm", names(DE.results.DE20$tpr))]
  DE.results.DE20$mean.fdr <- DE.results.DE20$mean.fdr[, grep("tmm", names(DE.results.DE20$mean.fdr))]
  DE.results.DE20$mean.discoveries <- DE.results.DE20$mean.discoveries[, grep("tmm", names(DE.results.DE20$mean.discoveries))]
  DE.results.DE50$auc <- DE.results.DE50$auc[, grep("tmm", names(DE.results.DE50$auc))]
  DE.results.DE50$pauc <- DE.results.DE50$pauc[, grep("tmm", names(DE.results.DE50$pauc))]
  DE.results.DE50$fdr <- DE.results.DE50$fdr[, grep("tmm", names(DE.results.DE50$fdr))]
  DE.results.DE50$tpr <- DE.results.DE50$tpr[, grep("tmm", names(DE.results.DE50$tpr))]
  DE.results.DE50$mean.fdr <- DE.results.DE50$mean.fdr[, grep("tmm", names(DE.results.DE50$mean.fdr))]
  DE.results.DE50$mean.discoveries <- DE.results.DE50$mean.discoveries[, grep("tmm", names(DE.results.DE50$mean.discoveries))]
  DEDD.results.DE2$auc <- DEDD.results.DE2$auc[, grep("tmm", names(DEDD.results.DE2$auc))]
  DEDD.results.DE2$pauc <- DEDD.results.DE2$pauc[, grep("tmm", names(DEDD.results.DE2$pauc))]
  DEDD.results.DE2$fdr <- DEDD.results.DE2$fdr[, grep("tmm", names(DEDD.results.DE2$fdr))]
  DEDD.results.DE2$tpr <- DEDD.results.DE2$tpr[, grep("tmm", names(DEDD.results.DE2$tpr))]
  DEDD.results.DE2$mean.fdr <- DEDD.results.DE2$mean.fdr[, grep("tmm", names(DEDD.results.DE2$mean.fdr))]
  DEDD.results.DE2$mean.discoveries <- DEDD.results.DE2$mean.discoveries[, grep("tmm", names(DEDD.results.DE2$mean.discoveries))]
  DEDD.results.DE5$auc <- DEDD.results.DE5$auc[, grep("tmm", names(DEDD.results.DE5$auc))]
  DEDD.results.DE5$pauc <- DEDD.results.DE5$pauc[, grep("tmm", names(DEDD.results.DE5$pauc))]
  DEDD.results.DE5$fdr <- DEDD.results.DE5$fdr[, grep("tmm", names(DEDD.results.DE5$fdr))]
  DEDD.results.DE5$tpr <- DEDD.results.DE5$tpr[, grep("tmm", names(DEDD.results.DE5$tpr))]
  DEDD.results.DE5$mean.fdr <- DEDD.results.DE5$mean.fdr[, grep("tmm", names(DEDD.results.DE5$mean.fdr))]
  DEDD.results.DE5$mean.discoveries <- DEDD.results.DE5$mean.discoveries[, grep("tmm", names(DEDD.results.DE5$mean.discoveries))]
  DEDD.results.DE10$auc <- DEDD.results.DE10$auc[, grep("tmm", names(DEDD.results.DE10$auc))]
  DEDD.results.DE10$pauc <- DEDD.results.DE10$pauc[, grep("tmm", names(DEDD.results.DE10$pauc))]
  DEDD.results.DE10$fdr <- DEDD.results.DE10$fdr[, grep("tmm", names(DEDD.results.DE10$fdr))]
  DEDD.results.DE10$tpr <- DEDD.results.DE10$tpr[, grep("tmm", names(DEDD.results.DE10$tpr))]
  DEDD.results.DE10$mean.fdr <- DEDD.results.DE10$mean.fdr[, grep("tmm", names(DEDD.results.DE10$mean.fdr))]
  DEDD.results.DE10$mean.discoveries <- DEDD.results.DE10$mean.discoveries[, grep("tmm", names(DEDD.results.DE10$mean.discoveries))]
  DEDD.results.DE20$auc <- DEDD.results.DE20$auc[, grep("tmm", names(DEDD.results.DE20$auc))]
  DEDD.results.DE20$pauc <- DEDD.results.DE20$pauc[, grep("tmm", names(DEDD.results.DE20$pauc))]
  DEDD.results.DE20$fdr <- DEDD.results.DE20$fdr[, grep("tmm", names(DEDD.results.DE20$fdr))]
  DEDD.results.DE20$tpr <- DEDD.results.DE20$tpr[, grep("tmm", names(DEDD.results.DE20$tpr))]
  DEDD.results.DE20$mean.fdr <- DEDD.results.DE20$mean.fdr[, grep("tmm", names(DEDD.results.DE20$mean.fdr))]
  DEDD.results.DE20$mean.discoveries <- DEDD.results.DE20$mean.discoveries[, grep("tmm", names(DEDD.results.DE20$mean.discoveries))]
  DEDD.results.DE50$auc <- DEDD.results.DE50$auc[, grep("tmm", names(DEDD.results.DE50$auc))]
  DEDD.results.DE50$pauc <- DEDD.results.DE50$pauc[, grep("tmm", names(DEDD.results.DE50$pauc))]
  DEDD.results.DE50$fdr <- DEDD.results.DE50$fdr[, grep("tmm", names(DEDD.results.DE50$fdr))]
  DEDD.results.DE50$tpr <- DEDD.results.DE50$tpr[, grep("tmm", names(DEDD.results.DE50$tpr))]
  DEDD.results.DE50$mean.fdr <- DEDD.results.DE50$mean.fdr[, grep("tmm", names(DEDD.results.DE50$mean.fdr))]
  DEDD.results.DE50$mean.discoveries <- DEDD.results.DE50$mean.discoveries[, grep("tmm", names(DEDD.results.DE50$mean.discoveries))]
}

{
  DD.results.DEDD2$auc <- DD.results.DEDD2$auc[, grep("tmm", names(DD.results.DEDD2$auc))]
  DD.results.DEDD2$pauc <- DD.results.DEDD2$pauc[, grep("tmm", names(DD.results.DEDD2$pauc))]
  DD.results.DEDD2$fdr <- DD.results.DEDD2$fdr[, grep("tmm", names(DD.results.DEDD2$fdr))]
  DD.results.DEDD2$tpr <- DD.results.DEDD2$tpr[, grep("tmm", names(DD.results.DEDD2$tpr))]
  DD.results.DEDD2$mean.fdr <- DD.results.DEDD2$mean.fdr[, grep("tmm", names(DD.results.DEDD2$mean.fdr))]
  DD.results.DEDD2$mean.discoveries <- DD.results.DEDD2$mean.discoveries[, grep("tmm", names(DD.results.DEDD2$mean.discoveries))]
  DD.results.DEDD5$auc <- DD.results.DEDD5$auc[, grep("tmm", names(DD.results.DEDD5$auc))]
  DD.results.DEDD5$pauc <- DD.results.DEDD5$pauc[, grep("tmm", names(DD.results.DEDD5$pauc))]
  DD.results.DEDD5$fdr <- DD.results.DEDD5$fdr[, grep("tmm", names(DD.results.DEDD5$fdr))]
  DD.results.DEDD5$tpr <- DD.results.DEDD5$tpr[, grep("tmm", names(DD.results.DEDD5$tpr))]
  DD.results.DEDD5$mean.fdr <- DD.results.DEDD5$mean.fdr[, grep("tmm", names(DD.results.DEDD5$mean.fdr))]
  DD.results.DEDD5$mean.discoveries <- DD.results.DEDD5$mean.discoveries[, grep("tmm", names(DD.results.DEDD5$mean.discoveries))]
  DD.results.DEDD10$auc <- DD.results.DEDD10$auc[, grep("tmm", names(DD.results.DEDD10$auc))]
  DD.results.DEDD10$pauc <- DD.results.DEDD10$pauc[, grep("tmm", names(DD.results.DEDD10$pauc))]
  DD.results.DEDD10$fdr <- DD.results.DEDD10$fdr[, grep("tmm", names(DD.results.DEDD10$fdr))]
  DD.results.DEDD10$tpr <- DD.results.DEDD10$tpr[, grep("tmm", names(DD.results.DEDD10$tpr))]
  DD.results.DEDD10$mean.fdr <- DD.results.DEDD10$mean.fdr[, grep("tmm", names(DD.results.DEDD10$mean.fdr))]
  DD.results.DEDD10$mean.discoveries <- DD.results.DEDD10$mean.discoveries[, grep("tmm", names(DD.results.DEDD10$mean.discoveries))]
  DD.results.DEDD20$auc <- DD.results.DEDD20$auc[, grep("tmm", names(DD.results.DEDD20$auc))]
  DD.results.DEDD20$pauc <- DD.results.DEDD20$pauc[, grep("tmm", names(DD.results.DEDD20$pauc))]
  DD.results.DEDD20$fdr <- DD.results.DEDD20$fdr[, grep("tmm", names(DD.results.DEDD20$fdr))]
  DD.results.DEDD20$tpr <- DD.results.DEDD20$tpr[, grep("tmm", names(DD.results.DEDD20$tpr))]
  DD.results.DEDD20$mean.fdr <- DD.results.DEDD20$mean.fdr[, grep("tmm", names(DD.results.DEDD20$mean.fdr))]
  DD.results.DEDD20$mean.discoveries <- DD.results.DEDD20$mean.discoveries[, grep("tmm", names(DD.results.DEDD20$mean.discoveries))]
  DD.results.DEDD50$auc <- DD.results.DEDD50$auc[, grep("tmm", names(DD.results.DEDD50$auc))]
  DD.results.DEDD50$pauc <- DD.results.DEDD50$pauc[, grep("tmm", names(DD.results.DEDD50$pauc))]
  DD.results.DEDD50$fdr <- DD.results.DEDD50$fdr[, grep("tmm", names(DD.results.DEDD50$fdr))]
  DD.results.DEDD50$tpr <- DD.results.DEDD50$tpr[, grep("tmm", names(DD.results.DEDD50$tpr))]
  DD.results.DEDD50$mean.fdr <- DD.results.DEDD50$mean.fdr[, grep("tmm", names(DD.results.DEDD50$mean.fdr))]
  DD.results.DEDD50$mean.discoveries <- DD.results.DEDD50$mean.discoveries[, grep("tmm", names(DD.results.DEDD50$mean.discoveries))]
  DE.results.DEDD2$auc <- DE.results.DEDD2$auc[, grep("tmm", names(DE.results.DEDD2$auc))]
  DE.results.DEDD2$pauc <- DE.results.DEDD2$pauc[, grep("tmm", names(DE.results.DEDD2$pauc))]
  DE.results.DEDD2$fdr <- DE.results.DEDD2$fdr[, grep("tmm", names(DE.results.DEDD2$fdr))]
  DE.results.DEDD2$tpr <- DE.results.DEDD2$tpr[, grep("tmm", names(DE.results.DEDD2$tpr))]
  DE.results.DEDD2$mean.fdr <- DE.results.DEDD2$mean.fdr[, grep("tmm", names(DE.results.DEDD2$mean.fdr))]
  DE.results.DEDD2$mean.discoveries <- DE.results.DEDD2$mean.discoveries[, grep("tmm", names(DE.results.DEDD2$mean.discoveries))]
  DE.results.DEDD5$auc <- DE.results.DEDD5$auc[, grep("tmm", names(DE.results.DEDD5$auc))]
  DE.results.DEDD5$pauc <- DE.results.DEDD5$pauc[, grep("tmm", names(DE.results.DEDD5$pauc))]
  DE.results.DEDD5$fdr <- DE.results.DEDD5$fdr[, grep("tmm", names(DE.results.DEDD5$fdr))]
  DE.results.DEDD5$tpr <- DE.results.DEDD5$tpr[, grep("tmm", names(DE.results.DEDD5$tpr))]
  DE.results.DEDD5$mean.fdr <- DE.results.DEDD5$mean.fdr[, grep("tmm", names(DE.results.DEDD5$mean.fdr))]
  DE.results.DEDD5$mean.discoveries <- DE.results.DEDD5$mean.discoveries[, grep("tmm", names(DE.results.DEDD5$mean.discoveries))]
  DE.results.DEDD10$auc <- DE.results.DEDD10$auc[, grep("tmm", names(DE.results.DEDD10$auc))]
  DE.results.DEDD10$pauc <- DE.results.DEDD10$pauc[, grep("tmm", names(DE.results.DEDD10$pauc))]
  DE.results.DEDD10$fdr <- DE.results.DEDD10$fdr[, grep("tmm", names(DE.results.DEDD10$fdr))]
  DE.results.DEDD10$tpr <- DE.results.DEDD10$tpr[, grep("tmm", names(DE.results.DEDD10$tpr))]
  DE.results.DEDD10$mean.fdr <- DE.results.DEDD10$mean.fdr[, grep("tmm", names(DE.results.DEDD10$mean.fdr))]
  DE.results.DEDD10$mean.discoveries <- DE.results.DEDD10$mean.discoveries[, grep("tmm", names(DE.results.DEDD10$mean.discoveries))]
  DE.results.DEDD20$auc <- DE.results.DEDD20$auc[, grep("tmm", names(DE.results.DEDD20$auc))]
  DE.results.DEDD20$pauc <- DE.results.DEDD20$pauc[, grep("tmm", names(DE.results.DEDD20$pauc))]
  DE.results.DEDD20$fdr <- DE.results.DEDD20$fdr[, grep("tmm", names(DE.results.DEDD20$fdr))]
  DE.results.DEDD20$tpr <- DE.results.DEDD20$tpr[, grep("tmm", names(DE.results.DEDD20$tpr))]
  DE.results.DEDD20$mean.fdr <- DE.results.DEDD20$mean.fdr[, grep("tmm", names(DE.results.DEDD20$mean.fdr))]
  DE.results.DEDD20$mean.discoveries <- DE.results.DEDD20$mean.discoveries[, grep("tmm", names(DE.results.DEDD20$mean.discoveries))]
  DE.results.DEDD50$auc <- DE.results.DEDD50$auc[, grep("tmm", names(DE.results.DEDD50$auc))]
  DE.results.DEDD50$pauc <- DE.results.DEDD50$pauc[, grep("tmm", names(DE.results.DEDD50$pauc))]
  DE.results.DEDD50$fdr <- DE.results.DEDD50$fdr[, grep("tmm", names(DE.results.DEDD50$fdr))]
  DE.results.DEDD50$tpr <- DE.results.DEDD50$tpr[, grep("tmm", names(DE.results.DEDD50$tpr))]
  DE.results.DEDD50$mean.fdr <- DE.results.DEDD50$mean.fdr[, grep("tmm", names(DE.results.DEDD50$mean.fdr))]
  DE.results.DEDD50$mean.discoveries <- DE.results.DEDD50$mean.discoveries[, grep("tmm", names(DE.results.DEDD50$mean.discoveries))]
  DEDD.results.DEDD2$auc <- DEDD.results.DEDD2$auc[, grep("tmm", names(DEDD.results.DEDD2$auc))]
  DEDD.results.DEDD2$pauc <- DEDD.results.DEDD2$pauc[, grep("tmm", names(DEDD.results.DEDD2$pauc))]
  DEDD.results.DEDD2$fdr <- DEDD.results.DEDD2$fdr[, grep("tmm", names(DEDD.results.DEDD2$fdr))]
  DEDD.results.DEDD2$tpr <- DEDD.results.DEDD2$tpr[, grep("tmm", names(DEDD.results.DEDD2$tpr))]
  DEDD.results.DEDD2$mean.fdr <- DEDD.results.DEDD2$mean.fdr[, grep("tmm", names(DEDD.results.DEDD2$mean.fdr))]
  DEDD.results.DEDD2$mean.discoveries <- DEDD.results.DEDD2$mean.discoveries[
    , grep("tmm", names(DEDD.results.DEDD2$mean.discoveries))
    ]
  DEDD.results.DEDD5$auc <- DEDD.results.DEDD5$auc[, grep("tmm", names(DEDD.results.DEDD5$auc))]
  DEDD.results.DEDD5$pauc <- DEDD.results.DEDD5$pauc[, grep("tmm", names(DEDD.results.DEDD5$pauc))]
  DEDD.results.DEDD5$fdr <- DEDD.results.DEDD5$fdr[, grep("tmm", names(DEDD.results.DEDD5$fdr))]
  DEDD.results.DEDD5$tpr <- DEDD.results.DEDD5$tpr[, grep("tmm", names(DEDD.results.DEDD5$tpr))]
  DEDD.results.DEDD5$mean.fdr <- DEDD.results.DEDD5$mean.fdr[, grep("tmm", names(DEDD.results.DEDD5$mean.fdr))]
  DEDD.results.DEDD5$mean.discoveries <- DEDD.results.DEDD5$mean.discoveries[
    , grep("tmm", names(DEDD.results.DEDD5$mean.discoveries))
    ]
  DEDD.results.DEDD10$auc <- DEDD.results.DEDD10$auc[, grep("tmm", names(DEDD.results.DEDD10$auc))]
  DEDD.results.DEDD10$pauc <- DEDD.results.DEDD10$pauc[, grep("tmm", names(DEDD.results.DEDD10$pauc))]
  DEDD.results.DEDD10$fdr <- DEDD.results.DEDD10$fdr[, grep("tmm", names(DEDD.results.DEDD10$fdr))]
  DEDD.results.DEDD10$tpr <- DEDD.results.DEDD10$tpr[, grep("tmm", names(DEDD.results.DEDD10$tpr))]
  DEDD.results.DEDD10$mean.fdr <- DEDD.results.DEDD10$mean.fdr[, grep("tmm", names(DEDD.results.DEDD10$mean.fdr))]
  DEDD.results.DEDD10$mean.discoveries <- DEDD.results.DEDD10$mean.discoveries[
    , grep("tmm", names(DEDD.results.DEDD10$mean.discoveries))
    ]
  DEDD.results.DEDD20$auc <- DEDD.results.DEDD20$auc[, grep("tmm", names(DEDD.results.DEDD20$auc))]
  DEDD.results.DEDD20$pauc <- DEDD.results.DEDD20$pauc[, grep("tmm", names(DEDD.results.DEDD20$pauc))]
  DEDD.results.DEDD20$fdr <- DEDD.results.DEDD20$fdr[, grep("tmm", names(DEDD.results.DEDD20$fdr))]
  DEDD.results.DEDD20$tpr <- DEDD.results.DEDD20$tpr[, grep("tmm", names(DEDD.results.DEDD20$tpr))]
  DEDD.results.DEDD20$mean.fdr <- DEDD.results.DEDD20$mean.fdr[, grep("tmm", names(DEDD.results.DEDD20$mean.fdr))]
  DEDD.results.DEDD20$mean.discoveries <- DEDD.results.DEDD20$mean.discoveries[
    , grep("tmm", names(DEDD.results.DEDD20$mean.discoveries))
    ]
  DEDD.results.DEDD50$auc <- DEDD.results.DEDD50$auc[, grep("tmm", names(DEDD.results.DEDD50$auc))]
  DEDD.results.DEDD50$pauc <- DEDD.results.DEDD50$pauc[, grep("tmm", names(DEDD.results.DEDD50$pauc))]
  DEDD.results.DEDD50$fdr <- DEDD.results.DEDD50$fdr[, grep("tmm", names(DEDD.results.DEDD50$fdr))]
  DEDD.results.DEDD50$tpr <- DEDD.results.DEDD50$tpr[, grep("tmm", names(DEDD.results.DEDD50$tpr))]
  DEDD.results.DEDD50$mean.fdr <- DEDD.results.DEDD50$mean.fdr[, grep("tmm", names(DEDD.results.DEDD50$mean.fdr))]
  DEDD.results.DEDD50$mean.discoveries <- DEDD.results.DEDD50$mean.discoveries[
    , grep("tmm", names(DEDD.results.DEDD50$mean.discoveries))
    ]
}
}

#################################
#### Differential dispersion ####
#################################

###########
### AUC ###
par(mfrow=c(2,5), mar=c(1,2,2,1), mgp=c(3,0.7,0))
boxplot(DD.results.DD2$auc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.5,0.9), main=paste0("AUC diff disp, DD2"), cex.main=1.5)
legend("topleft", fill=col_vector[1:4], bty='n', cex=1.5, ncol=2, legend=legend)
boxplot(DD.results.DD5$auc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.5,0.9), main=paste0("AUC diff disp, DD5"), cex.main=1.5)
boxplot(DD.results.DD10$auc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.5,0.9), main=paste0("AUC diff disp, DD10"), cex.main=1.5)
boxplot(DD.results.DD20$auc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.5,0.9), main=paste0("AUC diff disp, DD20"), cex.main=1.5)
boxplot(DD.results.DD50$auc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.5,0.9), main=paste0("AUC diff disp, DD50"), cex.main=1.5)
boxplot(DD.results.DEDD2$auc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.5,0.9), main=paste0("AUC diff disp, DEDD2"), cex.main=1.5)
boxplot(DD.results.DEDD5$auc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.5,0.9), main=paste0("AUC diff disp, DEDD5"), cex.main=1.5)
boxplot(DD.results.DEDD10$auc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.5,0.9), main=paste0("AUC diff disp, DEDD10"), cex.main=1.5)
boxplot(DD.results.DEDD20$auc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.5,0.9), main=paste0("AUC diff disp, DEDD50"), cex.main=1.5)
boxplot(DD.results.DEDD50$auc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.5,0.9), main=paste0("AUC diff disp, DEDD50"), cex.main=1.5)
rbind(colMeans(DD.results.DD2$auc), colMeans(DD.results.DD5$auc), colMeans(DD.results.DD10$auc), 
      colMeans(DD.results.DD20$auc), colMeans(DD.results.DD50$auc))
rbind(colMeans(DD.results.DEDD2$auc), colMeans(DD.results.DEDD5$auc), colMeans(DD.results.DEDD10$auc), 
      colMeans(DD.results.DEDD20$auc), colMeans(DD.results.DEDD50$auc))
# Essentially no differences in AUC, except for small samples, higher with log and for expHM.

###########
### pAUC ###
par(mfrow=c(2,5), mar=c(1,2,2,1), mgp=c(3,0.7,0))
boxplot(DD.results.DD2$pauc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.03), main=paste0("Partial AUC diff disp, DD2"), cex.main=1.5)
legend("topleft", fill=col_vector[1:4], bty='n', cex=1.5, ncol=2, legend=legend)
boxplot(DD.results.DD5$pauc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.03), main=paste0("Partial AUC diff disp, DD5"), cex.main=1.5)
boxplot(DD.results.DD10$pauc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.03), main=paste0("Partial AUC diff disp, DD10"), cex.main=1.5)
boxplot(DD.results.DD20$pauc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.03), main=paste0("Partial AUC diff disp, DD20"), cex.main=1.5)
boxplot(DD.results.DD50$pauc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.03), main=paste0("Partial AUC diff disp, DD50"), cex.main=1.5)
boxplot(DD.results.DEDD2$pauc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.03), main=paste0("Partial AUC diff disp, DEDD2"), cex.main=1.5)
boxplot(DD.results.DEDD5$pauc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.03), main=paste0("Partial AUC diff disp, DEDD5"), cex.main=1.5)
boxplot(DD.results.DEDD10$pauc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.03), main=paste0("Partial AUC diff disp, DEDD10"), cex.main=1.5)
boxplot(DD.results.DEDD20$pauc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.03), main=paste0("Partial AUC diff disp, DEDD20"), cex.main=1.5)
boxplot(DD.results.DEDD50$pauc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.03), main=paste0("Partial AUC diff disp, DEDD50"), cex.main=1.5)
rbind(colMeans(DD.results.DD2$pauc), colMeans(DD.results.DD5$pauc), colMeans(DD.results.DD10$pauc), 
      colMeans(DD.results.DD20$pauc), colMeans(DD.results.DD50$pauc))
rbind(colMeans(DD.results.DEDD2$pauc), colMeans(DD.results.DEDD5$pauc), colMeans(DD.results.DEDD10$pauc), 
      colMeans(DD.results.DEDD20$pauc), colMeans(DD.results.DEDD50$pauc))
# Essentially no differences in pAUC, except for small samples, higher with log and for expHM.

###########
### FDR ###
par(mfrow=c(2,5), mar=c(1,2,2,1), mgp=c(3,0.7,0))
boxplot(DD.results.DD2$fdr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.25), main=paste0("FDR diff disp, DD2"), cex.main=1.5)
abline(h=0.05, col='lightgrey')
legend("topleft", fill=col_vector[1:4], bty='n', cex=1.5, ncol=2, legend=legend)
boxplot(DD.results.DD5$fdr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.25), main=paste0("FDR diff disp, DD5"), cex.main=1.5)
abline(h=0.05, col='lightgrey')
boxplot(DD.results.DD10$fdr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.25), main=paste0("FDR diff disp, DD10"), cex.main=1.5)
abline(h=0.05, col='lightgrey')
boxplot(DD.results.DD20$fdr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.25), main=paste0("FDR diff disp, DD20"), cex.main=1.5)
abline(h=0.05, col='lightgrey')
boxplot(DD.results.DD50$fdr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.25), main=paste0("FDR diff disp, DD50"), cex.main=1.5)
abline(h=0.05, col='lightgrey')
boxplot(DD.results.DEDD2$fdr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.25), main=paste0("FDR diff disp, DEDD2"), cex.main=1.5)
abline(h=0.05, col='lightgrey')
boxplot(DD.results.DEDD5$fdr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.25), main=paste0("FDR diff disp, DEDD5"), cex.main=1.5)
abline(h=0.05, col='lightgrey')
boxplot(DD.results.DEDD10$fdr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.25), main=paste0("FDR diff disp, DEDD10"), cex.main=1.5)
abline(h=0.05, col='lightgrey')
boxplot(DD.results.DEDD20$fdr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.25), main=paste0("FDR diff disp, DEDD20"), cex.main=1.5)
abline(h=0.05, col='lightgrey')
boxplot(DD.results.DEDD50$fdr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.25), main=paste0("FDR diff disp, DEDD50"), cex.main=1.5)
abline(h=0.05, col='lightgrey')
# Mean FDRs:
rbind(colMeans(DD.results.DD2$fdr, na.rm=T), colMeans(DD.results.DD5$fdr, na.rm=T), colMeans(DD.results.DD10$fdr, na.rm=T), 
      colMeans(DD.results.DD20$fdr, na.rm=T), colMeans(DD.results.DD50$fdr, na.rm=T))
rbind(colMeans(DD.results.DEDD2$fdr, na.rm=T), colMeans(DD.results.DEDD5$fdr, na.rm=T), colMeans(DD.results.DEDD10$fdr, na.rm=T), 
      colMeans(DD.results.DEDD20$fdr, na.rm=T), colMeans(DD.results.DEDD50$fdr, na.rm=T))
# Mean squared distances from 0.05:
rbind(colMeans((DD.results.DD2$fdr - 0.05)^2, na.rm=T), colMeans((DD.results.DD5$fdr - 0.05)^2, na.rm=T), 
      colMeans((DD.results.DD10$fdr - 0.05)^2, na.rm=T), colMeans((DD.results.DD20$fdr - 0.05)^2, na.rm=T), 
      colMeans((DD.results.DD50$fdr - 0.05)^2, na.rm=T))
rbind(colMeans((DD.results.DEDD2$fdr - 0.05)^2, na.rm=T), colMeans((DD.results.DEDD5$fdr - 0.05)^2, na.rm=T), 
      colMeans((DD.results.DEDD10$fdr - 0.05)^2, na.rm=T), colMeans((DD.results.DEDD20$fdr - 0.05)^2, na.rm=T), 
      colMeans((DD.results.DEDD50$fdr - 0.05)^2, na.rm=T))
# Can only judge FDR for 20 and 50 samples per group. For 20, all are pretty bad, especially for DD only, lnHM better than expHM, 
# and no real difference between raw and log. For 50, lnHM generally closer to 0.05 than expHM, and less frequently over 0.05, and 
# untranslated generally closer to 0.05 than log, and less frequently over 0.05.

###########
### TPR ###
par(mfrow=c(2,5), mar=c(1,2,2,1), mgp=c(3,0.7,0))
boxplot(DD.results.DD2$tpr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.42), main=paste0("TPR diff disp, DD2"), cex.main=1.5)
legend("topleft", fill=col_vector[1:4], bty='n', cex=1.5, ncol=2, legend=legend)
boxplot(DD.results.DD5$tpr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.42), main=paste0("TPR diff disp, DD5"), cex.main=1.5)
boxplot(DD.results.DD10$tpr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.42), main=paste0("TPR diff disp, DD10"), cex.main=1.5)
boxplot(DD.results.DD20$tpr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.42), main=paste0("TPR diff disp, DD20"), cex.main=1.5)
boxplot(DD.results.DD50$tpr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.42), main=paste0("TPR diff disp, DD50"), cex.main=1.5)
boxplot(DD.results.DEDD2$tpr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.42), main=paste0("TPR diff disp, DEDD2"), cex.main=1.5)
boxplot(DD.results.DEDD5$tpr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.42), main=paste0("TPR diff disp, DEDD5"), cex.main=1.5)
boxplot(DD.results.DEDD10$tpr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.42), main=paste0("TPR diff disp, DEDD10"), cex.main=1.5)
boxplot(DD.results.DEDD20$tpr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.42), main=paste0("TPR diff disp, DEDD20"), cex.main=1.5)
boxplot(DD.results.DEDD50$tpr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.42), main=paste0("TPR diff disp, DEDD50"), cex.main=1.5)
rbind(colMeans(DD.results.DD2$tpr), colMeans(DD.results.DD5$tpr), colMeans(DD.results.DD10$tpr), 
      colMeans(DD.results.DD20$tpr), colMeans(DD.results.DD50$tpr))
rbind(colMeans(DD.results.DEDD2$tpr), colMeans(DD.results.DEDD5$tpr), colMeans(DD.results.DEDD10$tpr), 
      colMeans(DD.results.DEDD20$tpr), colMeans(DD.results.DEDD50$tpr))
# Can only judge TPR for 20 and 50 samples per group. Slightly higher for expHM than lnHM, and for log than untranslated, i.e. 
# higher TPR for situations where FDR control is poorer.

#############################
### False discovery plots ###
par(mfrow=c(2,5), mar=c(2,2,2,1), mgp=c(3,0.7,0))

plot(DD.results.DD2$mean.discoveries$disp.expHM.untr.tmm, DD.results.DD2$mean.fdr$disp.expHM.untr.tmm, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], main=paste0("DD2"), cex.main=1.5)
lines(DD.results.DD2$mean.discoveries$disp.expHM.log.tmm, DD.results.DD2$mean.fdr$disp.expHM.log.tmm, 
      col=col_vector[2])
lines(DD.results.DD2$mean.discoveries$disp.lnHM.untr.tmm, DD.results.DD2$mean.fdr$disp.lnHM.untr.tmm, 
      col=col_vector[3])
lines(DD.results.DD2$mean.discoveries$disp.lnHM.log.tmm, DD.results.DD2$mean.fdr$disp.lnHM.log.tmm, 
      col=col_vector[4])
abline(h=0.05, col='lightgrey')
plot(DD.results.DD5$mean.discoveries$disp.expHM.untr.tmm, DD.results.DD5$mean.fdr$disp.expHM.untr.tmm, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], main=paste0("DD5"), cex.main=1.5)
lines(DD.results.DD5$mean.discoveries$disp.expHM.log.tmm, DD.results.DD5$mean.fdr$disp.expHM.log.tmm, 
      col=col_vector[2])
lines(DD.results.DD5$mean.discoveries$disp.lnHM.untr.tmm, DD.results.DD5$mean.fdr$disp.lnHM.untr.tmm, 
      col=col_vector[3])
lines(DD.results.DD5$mean.discoveries$disp.lnHM.log.tmm, DD.results.DD5$mean.fdr$disp.lnHM.log.tmm, 
      col=col_vector[4])
abline(h=0.05, col='lightgrey')
plot(DD.results.DD10$mean.discoveries$disp.expHM.untr.tmm, DD.results.DD10$mean.fdr$disp.expHM.untr.tmm, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], main=paste0("DD10"), cex.main=1.5)
lines(DD.results.DD10$mean.discoveries$disp.expHM.log.tmm, DD.results.DD10$mean.fdr$disp.expHM.log.tmm, 
      col=col_vector[2])
lines(DD.results.DD10$mean.discoveries$disp.lnHM.untr.tmm, DD.results.DD10$mean.fdr$disp.lnHM.untr.tmm, 
      col=col_vector[3])
lines(DD.results.DD10$mean.discoveries$disp.lnHM.log.tmm, DD.results.DD10$mean.fdr$disp.lnHM.log.tmm, 
      col=col_vector[4])
abline(h=0.05, col='lightgrey')
plot(DD.results.DD20$mean.discoveries$disp.expHM.untr.tmm, DD.results.DD20$mean.fdr$disp.expHM.untr.tmm, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], main=paste0("DD20"), cex.main=1.5)
lines(DD.results.DD20$mean.discoveries$disp.expHM.log.tmm, DD.results.DD20$mean.fdr$disp.expHM.log.tmm, 
      col=col_vector[2])
lines(DD.results.DD20$mean.discoveries$disp.lnHM.untr.tmm, DD.results.DD20$mean.fdr$disp.lnHM.untr.tmm, 
      col=col_vector[3])
lines(DD.results.DD20$mean.fdr$disp.lnHM.untr.tmm, DD.results.DD20$mean.discoveries$disp.lnHM.log.tmm,  
      col=col_vector[4])
abline(h=0.05, col='lightgrey')
plot(DD.results.DD50$mean.discoveries$disp.expHM.untr.tmm, DD.results.DD50$mean.fdr$disp.expHM.untr.tmm, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], main=paste0("DD50"), cex.main=1.5)
lines(DD.results.DD50$mean.discoveries$disp.expHM.log.tmm, DD.results.DD50$mean.fdr$disp.expHM.log.tmm, 
      col=col_vector[2])
lines(DD.results.DD50$mean.discoveries$disp.lnHM.untr.tmm, DD.results.DD50$mean.fdr$disp.lnHM.untr.tmm, 
      col=col_vector[3])
lines(DD.results.DD50$mean.discoveries$disp.lnHM.log.tmm, DD.results.DD50$mean.fdr$disp.lnHM.log.tmm, 
       col=col_vector[4])
abline(h=0.05, col='lightgrey')
legend("topleft", bty='n', col=col_vector[1:4], lty=1, ncol=1, cex=1.5, legend=legend)

plot(DD.results.DEDD2$mean.discoveries$disp.expHM.untr.tmm, DD.results.DEDD2$mean.fdr$disp.expHM.untr.tmm, 
     type='l', ylim=c(0,0.85), xlim=c(0,2000), col=col_vector[1], main=paste0("DD2"), cex.main=1.5)
lines(DD.results.DEDD2$mean.discoveries$disp.expHM.log.tmm, DD.results.DEDD2$mean.fdr$disp.expHM.log.tmm, 
      col=col_vector[2])
lines(DD.results.DEDD2$mean.discoveries$disp.lnHM.untr.tmm, DD.results.DEDD2$mean.fdr$disp.lnHM.untr.tmm, 
      col=col_vector[3])
lines(DD.results.DEDD2$mean.discoveries$disp.lnHM.log.tmm, DD.results.DEDD2$mean.fdr$disp.lnHM.log.tmm, 
      col=col_vector[4])
abline(h=0.05, col='lightgrey')
plot(DD.results.DEDD5$mean.discoveries$disp.expHM.untr.tmm, DD.results.DEDD5$mean.fdr$disp.expHM.untr.tmm, 
     type='l', ylim=c(0,0.85), xlim=c(0,2000), col=col_vector[1], main=paste0("DD5"), cex.main=1.5)
lines(DD.results.DEDD5$mean.discoveries$disp.expHM.log.tmm, DD.results.DEDD5$mean.fdr$disp.expHM.log.tmm, 
      col=col_vector[2])
lines(DD.results.DEDD5$mean.discoveries$disp.lnHM.untr.tmm, DD.results.DEDD5$mean.fdr$disp.lnHM.untr.tmm, 
      col=col_vector[3])
lines(DD.results.DEDD5$mean.discoveries$disp.lnHM.log.tmm, DD.results.DEDD5$mean.fdr$disp.lnHM.log.tmm, 
      col=col_vector[4])
abline(h=0.05, col='lightgrey')
plot(DD.results.DEDD10$mean.discoveries$disp.expHM.untr.tmm, DD.results.DEDD10$mean.fdr$disp.expHM.untr.tmm, 
     type='l', ylim=c(0,0.85), xlim=c(0,2000), col=col_vector[1], main=paste0("DD10"), cex.main=1.5)
lines(DD.results.DEDD10$mean.discoveries$disp.expHM.log.tmm, DD.results.DEDD10$mean.fdr$disp.expHM.log.tmm, 
      col=col_vector[2])
lines(DD.results.DEDD10$mean.discoveries$disp.lnHM.untr.tmm, DD.results.DEDD10$mean.fdr$disp.lnHM.untr.tmm, 
      col=col_vector[3])
lines(DD.results.DEDD10$mean.discoveries$disp.lnHM.log.tmm, DD.results.DEDD10$mean.fdr$disp.lnHM.log.tmm, 
      col=col_vector[4])
abline(h=0.05, col='lightgrey')
plot(DD.results.DEDD20$mean.discoveries$disp.expHM.untr.tmm, DD.results.DEDD20$mean.fdr$disp.expHM.untr.tmm, 
     type='l', ylim=c(0,0.85), xlim=c(0,2000), col=col_vector[1], main=paste0("DD20"), cex.main=1.5)
lines(DD.results.DEDD20$mean.discoveries$disp.expHM.log.tmm, DD.results.DEDD20$mean.fdr$disp.expHM.log.tmm, 
      col=col_vector[2])
lines(DD.results.DEDD20$mean.discoveries$disp.lnHM.untr.tmm, DD.results.DEDD20$mean.fdr$disp.lnHM.untr.tmm, 
      col=col_vector[3])
lines(DD.results.DEDD20$mean.fdr$disp.lnHM.untr.tmm, DD.results.DEDD20$mean.discoveries$disp.lnHM.log.tmm,  
      col=col_vector[4])
abline(h=0.05, col='lightgrey')
plot(DD.results.DEDD50$mean.discoveries$disp.expHM.untr.tmm, DD.results.DEDD50$mean.fdr$disp.expHM.untr.tmm, 
     type='l', ylim=c(0,0.85), xlim=c(0,2000), col=col_vector[1], main=paste0("DD50"), cex.main=1.5)
lines(DD.results.DEDD50$mean.discoveries$disp.expHM.log.tmm, DD.results.DEDD50$mean.fdr$disp.expHM.log.tmm, 
      col=col_vector[2])
lines(DD.results.DEDD50$mean.discoveries$disp.lnHM.untr.tmm, DD.results.DEDD50$mean.fdr$disp.lnHM.untr.tmm, 
      col=col_vector[3])
lines(DD.results.DEDD50$mean.discoveries$disp.lnHM.log.tmm, DD.results.DEDD50$mean.fdr$disp.lnHM.log.tmm, 
      col=col_vector[4])
abline(h=0.05, col='lightgrey')
legend("topleft", bty='n', col=col_vector[1:4], lty=1, ncol=1, cex=1.5, legend=legend)
# FDR curves indistinguishable.


#################################
#### Differential expression ####
#################################

###########
### AUC ###
par(mfrow=c(2,5), mar=c(1,2,2,1), mgp=c(3,0.7,0))
boxplot(DE.results.DE2$auc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.78,1), main=paste0("AUC diff exp, DE2"), cex.main=1.5)
legend("topleft", fill=col_vector[1:4], bty='n', cex=1.5, ncol=2, legend=legend)
boxplot(DE.results.DE5$auc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.78,1), main=paste0("AUC diff exp, DE5"), cex.main=1.5)
boxplot(DE.results.DE10$auc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.78,1), main=paste0("AUC diff exp, DE10"), cex.main=1.5)
boxplot(DE.results.DE20$auc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.78,1), main=paste0("AUC diff exp, DE20"), cex.main=1.5)
boxplot(DE.results.DE50$auc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.78,1), main=paste0("AUC diff exp, DE50"), cex.main=1.5)
boxplot(DE.results.DEDD2$auc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.78,1), main=paste0("AUC diff exp, DEDD2"), cex.main=1.5)
boxplot(DE.results.DEDD5$auc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.78,1), main=paste0("AUC diff exp, DEDD5"), cex.main=1.5)
boxplot(DE.results.DEDD10$auc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.78,1), main=paste0("AUC diff exp, DEDD10"), cex.main=1.5)
boxplot(DE.results.DEDD20$auc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.78,1), main=paste0("AUC diff exp, DEDD50"), cex.main=1.5)
boxplot(DE.results.DEDD50$auc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.78,1), main=paste0("AUC diff exp, DEDD50"), cex.main=1.5)
rbind(colMeans(DE.results.DE2$auc), colMeans(DE.results.DE5$auc), colMeans(DE.results.DE10$auc), 
      colMeans(DE.results.DE20$auc), colMeans(DE.results.DE50$auc))
rbind(colMeans(DE.results.DEDD2$auc), colMeans(DE.results.DEDD5$auc), colMeans(DE.results.DEDD10$auc), 
      colMeans(DE.results.DEDD20$auc), colMeans(DE.results.DEDD50$auc))
# Essentially no differences in AUC except for small samples, for which lnHM is better than expHM and log is better than 
# untranslated.

###########
### pAUC ###
par(mfrow=c(2,5), mar=c(1,2,2,1), mgp=c(3,0.7,0))
boxplot(DE.results.DE2$pauc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.01,0.05), main=paste0("Partial AUC diff exp, DE2"), cex.main=1.5)
legend("topleft", fill=col_vector[1:4], bty='n', cex=1.5, ncol=2, legend=legend)
boxplot(DE.results.DE5$pauc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.01,0.05), main=paste0("Partial AUC diff exp, DE5"), cex.main=1.5)
boxplot(DE.results.DE10$pauc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.01,0.05), main=paste0("Partial AUC diff exp, DE10"), cex.main=1.5)
boxplot(DE.results.DE20$pauc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.01,0.05), main=paste0("Partial AUC diff exp, DE20"), cex.main=1.5)
boxplot(DE.results.DE50$pauc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.01,0.05), main=paste0("Partial AUC diff exp, DE50"), cex.main=1.5)
boxplot(DE.results.DEDD2$pauc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.01,0.05), main=paste0("Partial AUC diff exp, DEDD2"), cex.main=1.5)
boxplot(DE.results.DEDD5$pauc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.01,0.05), main=paste0("Partial AUC diff exp, DEDD5"), cex.main=1.5)
boxplot(DE.results.DEDD10$pauc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.01,0.05), main=paste0("Partial AUC diff exp, DEDD10"), cex.main=1.5)
boxplot(DE.results.DEDD20$pauc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.01,0.05), main=paste0("Partial AUC diff exp, DEDD20"), cex.main=1.5)
boxplot(DE.results.DEDD50$pauc, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.01,0.05), main=paste0("Partial AUC diff exp, DEDD50"), cex.main=1.5)
rbind(colMeans(DE.results.DE2$pauc), colMeans(DE.results.DE5$pauc), colMeans(DE.results.DE10$pauc), 
      colMeans(DE.results.DE20$pauc), colMeans(DE.results.DE50$pauc))
rbind(colMeans(DE.results.DEDD2$pauc), colMeans(DE.results.DEDD5$pauc), colMeans(DE.results.DEDD10$pauc), 
      colMeans(DE.results.DEDD20$pauc), colMeans(DE.results.DEDD50$pauc))
# Essentially no difference in pAUC except for small samples, for which lnHM is better than expHM and untranslated is better than 
# log.

###########
### FDR ###
par(mfrow=c(2,5), mar=c(1,2,2,1), mgp=c(3,0.7,0))
boxplot(DE.results.DE2$fdr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.15), main=paste0("FDR diff exp, DE2"), cex.main=1.5)
abline(h=0.05, col='lightgrey')
legend("topleft", fill=col_vector[1:4], bty='n', cex=1.5, ncol=2, legend=legend)
boxplot(DE.results.DE5$fdr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.15), main=paste0("FDR diff exp, DE5"), cex.main=1.5)
abline(h=0.05, col='lightgrey')
boxplot(DE.results.DE10$fdr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.15), main=paste0("FDR diff exp, DE10"), cex.main=1.5)
abline(h=0.05, col='lightgrey')
boxplot(DE.results.DE20$fdr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.15), main=paste0("FDR diff exp, DE20"), cex.main=1.5)
abline(h=0.05, col='lightgrey')
boxplot(DE.results.DE50$fdr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.15), main=paste0("FDR diff exp, DE50"), cex.main=1.5)
abline(h=0.05, col='lightgrey')
boxplot(DE.results.DEDD2$fdr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.15), main=paste0("FDR diff exp, DEDD2"), cex.main=1.5)
abline(h=0.05, col='lightgrey')
boxplot(DE.results.DEDD5$fdr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.15), main=paste0("FDR diff exp, DEDD5"), cex.main=1.5)
abline(h=0.05, col='lightgrey')
boxplot(DE.results.DEDD10$fdr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.15), main=paste0("FDR diff exp, DEDD10"), cex.main=1.5)
abline(h=0.05, col='lightgrey')
boxplot(DE.results.DEDD20$fdr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.15), main=paste0("FDR diff exp, DEDD20"), cex.main=1.5)
abline(h=0.05, col='lightgrey')
boxplot(DE.results.DEDD50$fdr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.15), main=paste0("FDR diff exp, DEDD50"), cex.main=1.5)
abline(h=0.05, col='lightgrey')
# Mean FDRs:
rbind(colMeans(DE.results.DE2$fdr, na.rm=T), colMeans(DE.results.DE5$fdr, na.rm=T), colMeans(DE.results.DE10$fdr, na.rm=T), 
      colMeans(DE.results.DE20$fdr, na.rm=T), colMeans(DE.results.DE50$fdr, na.rm=T))
rbind(colMeans(DE.results.DEDD2$fdr, na.rm=T), colMeans(DE.results.DEDD5$fdr, na.rm=T), colMeans(DE.results.DEDD10$fdr, na.rm=T), 
      colMeans(DE.results.DEDD20$fdr, na.rm=T), colMeans(DE.results.DEDD50$fdr, na.rm=T))
# Mean squared distances from 0.05:
rbind(colMeans((DE.results.DE2$fdr - 0.05)^2, na.rm=T), colMeans((DE.results.DE5$fdr - 0.05)^2, na.rm=T), 
      colMeans((DE.results.DE10$fdr - 0.05)^2, na.rm=T), colMeans((DE.results.DE20$fdr - 0.05)^2, na.rm=T), 
      colMeans((DE.results.DE50$fdr - 0.05)^2, na.rm=T))
rbind(colMeans((DE.results.DEDD2$fdr - 0.05)^2, na.rm=T), colMeans((DE.results.DEDD5$fdr - 0.05)^2, na.rm=T), 
      colMeans((DE.results.DEDD10$fdr - 0.05)^2, na.rm=T), colMeans((DE.results.DEDD20$fdr - 0.05)^2, na.rm=T), 
      colMeans((DE.results.DEDD50$fdr - 0.05)^2, na.rm=T))
# With differences in mean only, lnHM and untranslated give lower FDRs and closer on average to 0.05 than lnHM and log. With 
# differences in mean and dispersion, there's no real difference between expHM and lnHM, and log is generally closer to 0.05 than 
# untranslated but more often exceeds 0.05, although still with generally acceptable FDRs except for a small number of cases with 5 
# samples per group.

###########
### TPR ###
par(mfrow=c(2,5), mar=c(1,2,2,1), mgp=c(3,0.7,0))
boxplot(DE.results.DE2$tpr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.13,0.96), main=paste0("TPR diff exp, DE2"), cex.main=1.5)
legend("topleft", fill=col_vector[1:4], bty='n', cex=1.5, ncol=2, legend=legend)
boxplot(DE.results.DE5$tpr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.13,0.96), main=paste0("TPR diff exp, DE5"), cex.main=1.5)
boxplot(DE.results.DE10$tpr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.13,0.96), main=paste0("TPR diff exp, DE10"), cex.main=1.5)
boxplot(DE.results.DE20$tpr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.13,0.96), main=paste0("TPR diff exp, DE20"), cex.main=1.5)
boxplot(DE.results.DE50$tpr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.13,0.96), main=paste0("TPR diff exp, DE50"), cex.main=1.5)
boxplot(DE.results.DEDD2$tpr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.13,0.96), main=paste0("TPR diff exp, DEDD2"), cex.main=1.5)
boxplot(DE.results.DEDD5$tpr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.13,0.96), main=paste0("TPR diff exp, DEDD5"), cex.main=1.5)
boxplot(DE.results.DEDD10$tpr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.13,0.96), main=paste0("TPR diff exp, DEDD10"), cex.main=1.5)
boxplot(DE.results.DEDD20$tpr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.13,0.96), main=paste0("TPR diff exp, DEDD20"), cex.main=1.5)
boxplot(DE.results.DEDD50$tpr, names=NA, col=col_vector[1:4], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.13,0.96), main=paste0("TPR diff exp, DEDD50"), cex.main=1.5)
rbind(colMeans(DE.results.DE2$tpr), colMeans(DE.results.DE5$tpr), colMeans(DE.results.DE10$tpr), 
      colMeans(DE.results.DE20$tpr), colMeans(DE.results.DE50$tpr))
rbind(colMeans(DE.results.DEDD2$tpr), colMeans(DE.results.DEDD5$tpr), colMeans(DE.results.DEDD10$tpr), 
      colMeans(DE.results.DEDD20$tpr), colMeans(DE.results.DEDD50$tpr))
# Can only judge TPR for at least 5 samples per group. Always higher for log than untranslated, and higher for lnHM than expHM 
# except for 20 and 50 samples per group, where expHM is very slightly higher.

#############################
### False discovery plots ###
par(mfrow=c(2,5), mar=c(2,2,2,1), mgp=c(3,0.7,0))

plot(DE.results.DE2$mean.discoveries$mean.expHM.untr.tmm, DE.results.DE2$mean.fdr$mean.expHM.untr.tmm, 
     type='l', ylim=c(0,0.8), xlim=c(0,2000), col=col_vector[1], main=paste0("DD2"), cex.main=1.5)
lines(DE.results.DE2$mean.discoveries$mean.expHM.log.tmm, DE.results.DE2$mean.fdr$mean.expHM.log.tmm, 
      type='l', col=col_vector[2])
lines(DE.results.DE2$mean.discoveries$mean.lnHM.untr.tmm, DE.results.DE2$mean.fdr$mean.lnHM.untr.tmm, 
      type='l', col=col_vector[3])
lines(DE.results.DE2$mean.discoveries$mean.lnHM.log.tmm, DE.results.DE2$mean.fdr$mean.lnHM.log.tmm, 
      type='l', col=col_vector[4])
abline(h=0.05, col='lightgrey')
plot(DE.results.DE5$mean.discoveries$mean.expHM.untr.tmm, DE.results.DE5$mean.fdr$mean.expHM.untr.tmm, 
     type='l', ylim=c(0,0.8), xlim=c(0,2000), col=col_vector[1], main=paste0("DD5"), cex.main=1.5)
lines(DE.results.DE5$mean.discoveries$mean.expHM.log.tmm, DE.results.DE5$mean.fdr$mean.expHM.log.tmm, 
      type='l', col=col_vector[2])
lines(DE.results.DE5$mean.discoveries$mean.lnHM.untr.tmm, DE.results.DE5$mean.fdr$mean.lnHM.untr.tmm, 
      type='l', col=col_vector[3])
lines(DE.results.DE5$mean.discoveries$mean.lnHM.log.tmm, DE.results.DE5$mean.fdr$mean.lnHM.log.tmm, 
      type='l', col=col_vector[4])
abline(h=0.05, col='lightgrey')
plot(DE.results.DE10$mean.discoveries$mean.expHM.untr.tmm, DE.results.DE10$mean.fdr$mean.expHM.untr.tmm, 
     type='l', ylim=c(0,0.8), xlim=c(0,2000), col=col_vector[1], main=paste0("DD10"), cex.main=1.5)
lines(DE.results.DE10$mean.discoveries$mean.expHM.log.tmm, DE.results.DE10$mean.fdr$mean.expHM.log.tmm, 
      type='l', col=col_vector[2])
lines(DE.results.DE10$mean.discoveries$mean.lnHM.untr.tmm, DE.results.DE10$mean.fdr$mean.lnHM.untr.tmm, 
      type='l', col=col_vector[3])
lines(DE.results.DE10$mean.discoveries$mean.lnHM.log.tmm, DE.results.DE10$mean.fdr$mean.lnHM.log.tmm, 
      type='l', col=col_vector[4])
abline(h=0.05, col='lightgrey')
plot(DE.results.DE20$mean.discoveries$mean.expHM.untr.tmm, DE.results.DE20$mean.fdr$mean.expHM.untr.tmm, 
     type='l', ylim=c(0,0.8), xlim=c(0,2000), col=col_vector[1], main=paste0("DD20"), cex.main=1.5)
lines(DE.results.DE20$mean.discoveries$mean.expHM.log.tmm, DE.results.DE20$mean.fdr$mean.expHM.log.tmm, 
      type='l', col=col_vector[2])
lines(DE.results.DE20$mean.discoveries$mean.lnHM.untr.tmm, DE.results.DE20$mean.fdr$mean.lnHM.untr.tmm, 
      type='l', col=col_vector[3])
lines(DE.results.DE50$mean.fdr$mean.lnHM.untr.tmm, DE.results.DE50$mean.discoveries$mean.lnHM.log.tmm,  
      type='l', col=col_vector[4])
abline(h=0.05, col='lightgrey')
plot(DE.results.DE50$mean.discoveries$mean.expHM.untr.tmm, DE.results.DE50$mean.fdr$mean.expHM.untr.tmm, 
     type='l', ylim=c(0,0.8), xlim=c(0,2000), col=col_vector[1], main=paste0("DD50"), cex.main=1.5)
lines(DE.results.DE50$mean.discoveries$mean.expHM.log.tmm, DE.results.DE50$mean.fdr$mean.expHM.log.tmm, 
      type='l', col=col_vector[2])
lines(DE.results.DE50$mean.discoveries$mean.lnHM.untr.tmm, DE.results.DE50$mean.fdr$mean.lnHM.untr.tmm, 
      type='l', col=col_vector[3])
lines(DE.results.DE50$mean.discoveries$mean.lnHM.log.tmm, DE.results.DE50$mean.fdr$mean.lnHM.log.tmm, 
      type='l', col=col_vector[4])
abline(h=0.05, col='lightgrey')
legend("topleft", bty='n', col=col_vector[1:4], lty=1, ncol=1, cex=1.5, legend=legend)

plot(DE.results.DEDD2$mean.discoveries$mean.expHM.untr.tmm, DE.results.DEDD2$mean.fdr$mean.expHM.untr.tmm, 
     type='l', ylim=c(0,0.6), xlim=c(0,2000), col=col_vector[1], main=paste0("DEDD2"), cex.main=1.5)
lines(DE.results.DEDD2$mean.discoveries$mean.expHM.log.tmm, DE.results.DEDD2$mean.fdr$mean.expHM.log.tmm, 
      type='l', col=col_vector[2])
lines(DE.results.DEDD2$mean.discoveries$mean.lnHM.untr.tmm, DE.results.DEDD2$mean.fdr$mean.lnHM.untr.tmm, 
      type='l', col=col_vector[3])
lines(DE.results.DEDD2$mean.discoveries$mean.lnHM.log.tmm, DE.results.DEDD2$mean.fdr$mean.lnHM.log.tmm, 
      type='l', col=col_vector[4])
abline(h=0.05, col='lightgrey')
plot(DE.results.DEDD5$mean.discoveries$mean.expHM.untr.tmm, DE.results.DEDD5$mean.fdr$mean.expHM.untr.tmm, 
     type='l', ylim=c(0,0.6), xlim=c(0,2000), col=col_vector[1], main=paste0("DEDD5"), cex.main=1.5)
lines(DE.results.DEDD5$mean.discoveries$mean.expHM.log.tmm, DE.results.DEDD5$mean.fdr$mean.expHM.log.tmm, 
      type='l', col=col_vector[2])
lines(DE.results.DEDD5$mean.discoveries$mean.lnHM.untr.tmm, DE.results.DEDD5$mean.fdr$mean.lnHM.untr.tmm, 
      type='l', col=col_vector[3])
lines(DE.results.DEDD5$mean.discoveries$mean.lnHM.log.tmm, DE.results.DEDD5$mean.fdr$mean.lnHM.log.tmm, 
      type='l', col=col_vector[4])
abline(h=0.05, col='lightgrey')
plot(DE.results.DEDD10$mean.discoveries$mean.expHM.untr.tmm, DE.results.DEDD10$mean.fdr$mean.expHM.untr.tmm, 
     type='l', ylim=c(0,0.6), xlim=c(0,2000), col=col_vector[1], main=paste0("DEDD10"), cex.main=1.5)
lines(DE.results.DEDD10$mean.discoveries$mean.expHM.log.tmm, DE.results.DEDD10$mean.fdr$mean.expHM.log.tmm, 
      type='l', col=col_vector[2])
lines(DE.results.DEDD10$mean.discoveries$mean.lnHM.untr.tmm, DE.results.DEDD10$mean.fdr$mean.lnHM.untr.tmm, 
      type='l', col=col_vector[3])
lines(DE.results.DEDD10$mean.discoveries$mean.lnHM.log.tmm, DE.results.DEDD10$mean.fdr$mean.lnHM.log.tmm, 
      type='l', col=col_vector[4])
abline(h=0.05, col='lightgrey')
plot(DE.results.DEDD20$mean.discoveries$mean.expHM.untr.tmm, DE.results.DEDD20$mean.fdr$mean.expHM.untr.tmm, 
     type='l', ylim=c(0,0.6), xlim=c(0,2000), col=col_vector[1], main=paste0("DEDD20"), cex.main=1.5)
lines(DE.results.DEDD20$mean.discoveries$mean.expHM.log.tmm, DE.results.DEDD20$mean.fdr$mean.expHM.log.tmm, 
      type='l', col=col_vector[2])
lines(DE.results.DEDD20$mean.discoveries$mean.lnHM.untr.tmm, DE.results.DEDD20$mean.fdr$mean.lnHM.untr.tmm, 
      type='l', col=col_vector[3])
lines(DE.results.DEDD50$mean.fdr$mean.lnHM.untr.tmm, DE.results.DEDD50$mean.discoveries$mean.lnHM.log.tmm,  
      type='l', col=col_vector[4])
abline(h=0.05, col='lightgrey')
plot(DE.results.DEDD50$mean.discoveries$mean.expHM.untr.tmm, DE.results.DEDD50$mean.fdr$mean.expHM.untr.tmm, 
     type='l', ylim=c(0,0.6), xlim=c(0,2000), col=col_vector[1], main=paste0("DEDD50"), cex.main=1.5)
lines(DE.results.DEDD50$mean.discoveries$mean.expHM.log.tmm, DE.results.DEDD50$mean.fdr$mean.expHM.log.tmm, 
      type='l', col=col_vector[2])
lines(DE.results.DEDD50$mean.discoveries$mean.lnHM.untr.tmm, DE.results.DEDD50$mean.fdr$mean.lnHM.untr.tmm, 
      type='l', col=col_vector[3])
lines(DE.results.DEDD50$mean.discoveries$mean.lnHM.log.tmm, DE.results.DEDD50$mean.fdr$mean.lnHM.log.tmm, 
      type='l', col=col_vector[4])
abline(h=0.05, col='lightgrey')
legend("topleft", bty='n', col=col_vector[1:4], lty=1, ncol=1, cex=1.5, legend=legend)
# Virtually indistinguishable for 10 samples per group upwards. For 2 and 5, lnHM slightly better than expHM, and untranslated very 
# slightly better than log.


###################################
#### Differential distribution ####
###################################
legend <- c("expHMM", "lnHMM")

###########
### AUC ###
par(mfrow=c(3,5), mar=c(1,2,2,1), mgp=c(3,0.7,0))
boxplot(DEDD.results.DD2$auc, names=NA, col=col_vector[1:2], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.47,0.81), main=paste0("AUC diff exp, DD2"), cex.main=1.5)
legend("topleft", fill=col_vector[1:2], bty='n', cex=1.5, ncol=2, legend=legend)
boxplot(DEDD.results.DD5$auc, names=NA, col=col_vector[1:2], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.47,0.81), main=paste0("AUC diff exp, DD5"), cex.main=1.5)
boxplot(DEDD.results.DD10$auc, names=NA, col=col_vector[1:2], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.47,0.81), main=paste0("AUC diff exp, DD10"), cex.main=1.5)
boxplot(DEDD.results.DD20$auc, names=NA, col=col_vector[1:2], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.47,0.81), main=paste0("AUC diff exp, DD20"), cex.main=1.5)
boxplot(DEDD.results.DD50$auc, names=NA, col=col_vector[1:2], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.47,0.81), main=paste0("AUC diff exp, DD50"), cex.main=1.5)
boxplot(DEDD.results.DE2$auc, names=NA, col=col_vector[1:2], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.63,1), main=paste0("AUC diff exp, DE2"), cex.main=1.5)
boxplot(DEDD.results.DE5$auc, names=NA, col=col_vector[1:2], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.63,1), main=paste0("AUC diff exp, DE5"), cex.main=1.5)
boxplot(DEDD.results.DE10$auc, names=NA, col=col_vector[1:2], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.63,1), main=paste0("AUC diff exp, DE10"), cex.main=1.5)
boxplot(DEDD.results.DE20$auc, names=NA, col=col_vector[1:2], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.63,1), main=paste0("AUC diff exp, DE20"), cex.main=1.5)
boxplot(DEDD.results.DE50$auc, names=NA, col=col_vector[1:2], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.63,1), main=paste0("AUC diff exp, DE50"), cex.main=1.5)
boxplot(DEDD.results.DEDD2$auc, names=NA, col=col_vector[1:2], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.6,0.94), main=paste0("AUC diff exp, DEDD2"), cex.main=1.5)
boxplot(DEDD.results.DEDD5$auc, names=NA, col=col_vector[1:2], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.6,0.94), main=paste0("AUC diff exp, DEDD5"), cex.main=1.5)
boxplot(DEDD.results.DEDD10$auc, names=NA, col=col_vector[1:2], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.6,0.94), main=paste0("AUC diff exp, DEDD10"), cex.main=1.5)
boxplot(DEDD.results.DEDD20$auc, names=NA, col=col_vector[1:2], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.6,0.94), main=paste0("AUC diff exp, DEDD50"), cex.main=1.5)
boxplot(DEDD.results.DEDD50$auc, names=NA, col=col_vector[1:2], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.6,0.94), main=paste0("AUC diff exp, DEDD50"), cex.main=1.5)
rbind(colMeans(DEDD.results.DD2$auc), colMeans(DEDD.results.DD5$auc), colMeans(DEDD.results.DD10$auc), 
      colMeans(DEDD.results.DD20$auc), colMeans(DEDD.results.DD50$auc))
rbind(colMeans(DEDD.results.DE2$auc), colMeans(DEDD.results.DE5$auc), colMeans(DEDD.results.DE10$auc), 
      colMeans(DEDD.results.DE20$auc), colMeans(DEDD.results.DE50$auc))
rbind(colMeans(DEDD.results.DEDD2$auc), colMeans(DEDD.results.DEDD5$auc), colMeans(DEDD.results.DEDD10$auc), 
      colMeans(DEDD.results.DEDD20$auc), colMeans(DEDD.results.DEDD50$auc))
# AUC virtually always higher for lnHMM, but very small differences except for 2 to 10 samples per group and with differences in 
# dispersion only for 20 and 50 samples per group.

###########
### pAUC ###
par(mfrow=c(3,5), mar=c(1,2,2,1), mgp=c(3,0.7,0))
boxplot(DEDD.results.DD2$pauc, names=NA, col=col_vector[1:2], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.023), main=paste0("Partial AUC diff exp, DE2"), cex.main=1.5)
legend("topleft", fill=col_vector[1:2], bty='n', cex=1.5, ncol=2, legend=legend)
boxplot(DEDD.results.DD5$pauc, names=NA, col=col_vector[1:2], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.023), main=paste0("Partial AUC diff exp, DE5"), cex.main=1.5)
boxplot(DEDD.results.DD10$pauc, names=NA, col=col_vector[1:2], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.023), main=paste0("Partial AUC diff exp, DE10"), cex.main=1.5)
boxplot(DEDD.results.DD20$pauc, names=NA, col=col_vector[1:2], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.023), main=paste0("Partial AUC diff exp, DE20"), cex.main=1.5)
boxplot(DEDD.results.DD50$pauc, names=NA, col=col_vector[1:2], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.023), main=paste0("Partial AUC diff exp, DE50"), cex.main=1.5)
boxplot(DEDD.results.DE2$pauc, names=NA, col=col_vector[1:2], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.007,0.048), main=paste0("Partial AUC diff exp, DE2"), cex.main=1.5)
boxplot(DEDD.results.DE5$pauc, names=NA, col=col_vector[1:2], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.007,0.048), main=paste0("Partial AUC diff exp, DE5"), cex.main=1.5)
boxplot(DEDD.results.DE10$pauc, names=NA, col=col_vector[1:2], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.007,0.048), main=paste0("Partial AUC diff exp, DE10"), cex.main=1.5)
boxplot(DEDD.results.DE20$pauc, names=NA, col=col_vector[1:2], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.007,0.048), main=paste0("Partial AUC diff exp, DE20"), cex.main=1.5)
boxplot(DEDD.results.DE50$pauc, names=NA, col=col_vector[1:2], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.007,0.048), main=paste0("Partial AUC diff exp, DE50"), cex.main=1.5)
boxplot(DEDD.results.DEDD2$pauc, names=NA, col=col_vector[1:2], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.005,0.04), main=paste0("Partial AUC diff exp, DEDD2"), cex.main=1.5)
boxplot(DEDD.results.DEDD5$pauc, names=NA, col=col_vector[1:2], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.005,0.04), main=paste0("Partial AUC diff exp, DEDD5"), cex.main=1.5)
boxplot(DEDD.results.DEDD10$pauc, names=NA, col=col_vector[1:2], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.005,0.04), main=paste0("Partial AUC diff exp, DEDD10"), cex.main=1.5)
boxplot(DEDD.results.DEDD20$pauc, names=NA, col=col_vector[1:2], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.005,0.04), main=paste0("Partial AUC diff exp, DEDD20"), cex.main=1.5)
boxplot(DEDD.results.DEDD50$pauc, names=NA, col=col_vector[1:2], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0.005,0.04), main=paste0("Partial AUC diff exp, DEDD50"), cex.main=1.5)
rbind(colMeans(DEDD.results.DD2$pauc), colMeans(DEDD.results.DD5$pauc), colMeans(DEDD.results.DD10$pauc), 
      colMeans(DEDD.results.DD20$pauc), colMeans(DEDD.results.DD50$pauc))
rbind(colMeans(DEDD.results.DE2$pauc), colMeans(DEDD.results.DE5$pauc), colMeans(DEDD.results.DE10$pauc), 
      colMeans(DEDD.results.DE20$pauc), colMeans(DEDD.results.DE50$pauc))
rbind(colMeans(DEDD.results.DEDD2$pauc), colMeans(DEDD.results.DEDD5$pauc), colMeans(DEDD.results.DEDD10$pauc), 
      colMeans(DEDD.results.DEDD20$pauc), colMeans(DEDD.results.DEDD50$pauc))
# Essentially no differences in pAUC, except lnHMM lower than expHMM for 50 samples per group and slightly higher for 5 samples per 
# group.

###########
### FDR ###
par(mfrow=c(3,5), mar=c(1,2,2,1), mgp=c(3,0.7,0))
boxplot(DEDD.results.DD2$fdr, names=NA, col=col_vector[1:3], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.1), main=paste0("FDR diff exp, DD2"), cex.main=1.5)
abline(h=0.05, col='lightgrey')
legend("topleft", fill=col_vector[1:3], bty='n', cex=1.5, ncol=1, 
       legend=c("p = 0.5 threshold", "Posterior threshold", "bfdr"))
boxplot(DEDD.results.DD5$fdr, names=NA, col=col_vector[1:3], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.1), main=paste0("FDR diff exp, DD5"), cex.main=1.5)
abline(h=0.05, col='lightgrey')
boxplot(DEDD.results.DD10$fdr, names=NA, col=col_vector[1:3], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.1), main=paste0("FDR diff exp, DD10"), cex.main=1.5)
abline(h=0.05, col='lightgrey')
boxplot(DEDD.results.DD20$fdr, names=NA, col=col_vector[1:3], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.1), main=paste0("FDR diff exp, DD20"), cex.main=1.5)
abline(h=0.05, col='lightgrey')
boxplot(DEDD.results.DD50$fdr, names=NA, col=col_vector[1:3], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.1), main=paste0("FDR diff exp, DD50"), cex.main=1.5)
abline(h=0.05, col='lightgrey')
boxplot(DEDD.results.DE2$fdr, names=NA, col=col_vector[1:3], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.1), main=paste0("FDR diff exp, DE2"), cex.main=1.5)
abline(h=0.05, col='lightgrey')
boxplot(DEDD.results.DE5$fdr, names=NA, col=col_vector[1:3], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.1), main=paste0("FDR diff exp, DE5"), cex.main=1.5)
abline(h=0.05, col='lightgrey')
boxplot(DEDD.results.DE10$fdr, names=NA, col=col_vector[1:3], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.1), main=paste0("FDR diff exp, DE10"), cex.main=1.5)
abline(h=0.05, col='lightgrey')
boxplot(DEDD.results.DE20$fdr, names=NA, col=col_vector[1:3], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.1), main=paste0("FDR diff exp, DE20"), cex.main=1.5)
abline(h=0.05, col='lightgrey')
boxplot(DEDD.results.DE50$fdr, names=NA, col=col_vector[1:3], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.1), main=paste0("FDR diff exp, DE50"), cex.main=1.5)
abline(h=0.05, col='lightgrey')
boxplot(DEDD.results.DEDD2$fdr, names=NA, col=col_vector[1:3], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.1), main=paste0("FDR diff exp, DEDD2"), cex.main=1.5)
abline(h=0.05, col='lightgrey')
boxplot(DEDD.results.DEDD5$fdr, names=NA, col=col_vector[1:3], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.1), main=paste0("FDR diff exp, DEDD5"), cex.main=1.5)
abline(h=0.05, col='lightgrey')
boxplot(DEDD.results.DEDD10$fdr, names=NA, col=col_vector[1:3], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.1), main=paste0("FDR diff exp, DEDD10"), cex.main=1.5)
abline(h=0.05, col='lightgrey')
boxplot(DEDD.results.DEDD20$fdr, names=NA, col=col_vector[1:3], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.1), main=paste0("FDR diff exp, DEDD20"), cex.main=1.5)
abline(h=0.05, col='lightgrey')
boxplot(DEDD.results.DEDD50$fdr, names=NA, col=col_vector[1:3], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.1), main=paste0("FDR diff exp, DEDD50"), cex.main=1.5)
abline(h=0.05, col='lightgrey')
rbind(colMeans(DEDD.results.DEDD2$fdr, na.rm=T), colMeans(DEDD.results.DEDD5$fdr, na.rm=T), 
      colMeans(DEDD.results.DEDD10$fdr, na.rm=T), colMeans(DEDD.results.DEDD20$fdr, na.rm=T), 
      colMeans(DEDD.results.DEDD50$fdr, na.rm=T))
rbind(colMeans(DEDD.results.DE2$fdr, na.rm=T), colMeans(DEDD.results.DE5$fdr, na.rm=T), 
      colMeans(DEDD.results.DE10$fdr, na.rm=T), colMeans(DEDD.results.DE20$fdr, na.rm=T), 
      colMeans(DEDD.results.DE50$fdr, na.rm=T))
rbind(colMeans(DEDD.results.DEDD2$fdr, na.rm=T), colMeans(DEDD.results.DEDD5$fdr, na.rm=T), 
      colMeans(DEDD.results.DEDD10$fdr, na.rm=T), colMeans(DEDD.results.DEDD20$fdr, na.rm=T), 
      colMeans(DEDD.results.DEDD50$fdr, na.rm=T))
# Can't really judge FDRs at all as they're all extremely low except for 2 samples per group and for 5 and 10 samples per group with 
# differences in dispersion only, where they're either undefined or close to 1. lnHMM looks to generally be better than expHMM in 
# that with 10 or more samples per group the FDRs are slightly higher, but still far below 0.05. Posterior threshold also generally 
# gives FDRs slightly closer to 0.05 than BFDR and 0.5 posterior probability threshold, but posterior threshold isn't aiming at any 
# particular FDR, so can really just say that it isn't giving problematically high FDRs. BFDR is trying to control FDR at 0.05, and 
# is very conservative.

###########
### TPR ###
par(mfrow=c(3,5), mar=c(1,2,2,1), mgp=c(3,0.7,0))
boxplot(DEDD.results.DD2$tpr, names=NA, col=col_vector[1:3], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.1), main=paste0("TPR diff exp, DE2"), cex.main=1.5)
legend("topleft", fill=col_vector[1:3], bty='n', cex=1.5, ncol=1, 
       legend=c("p = 0.5 threshold", "Posterior threshold", "bfdr"))
boxplot(DEDD.results.DD5$tpr, names=NA, col=col_vector[1:3], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.1), main=paste0("TPR diff exp, DE5"), cex.main=1.5)
boxplot(DEDD.results.DD10$tpr, names=NA, col=col_vector[1:3], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.1), main=paste0("TPR diff exp, DE10"), cex.main=1.5)
boxplot(DEDD.results.DD20$tpr, names=NA, col=col_vector[1:3], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.1), main=paste0("TPR diff exp, DE20"), cex.main=1.5)
boxplot(DEDD.results.DD50$tpr, names=NA, col=col_vector[1:3], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.1), main=paste0("TPR diff exp, DE50"), cex.main=1.5)
boxplot(DEDD.results.DE2$tpr, names=NA, col=col_vector[1:3], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.92), main=paste0("TPR diff exp, DE2"), cex.main=1.5)
boxplot(DEDD.results.DE5$tpr, names=NA, col=col_vector[1:3], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.92), main=paste0("TPR diff exp, DE5"), cex.main=1.5)
boxplot(DEDD.results.DE10$tpr, names=NA, col=col_vector[1:3], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.92), main=paste0("TPR diff exp, DE10"), cex.main=1.5)
boxplot(DEDD.results.DE20$tpr, names=NA, col=col_vector[1:3], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.92), main=paste0("TPR diff exp, DE20"), cex.main=1.5)
boxplot(DEDD.results.DE50$tpr, names=NA, col=col_vector[1:3], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.92), main=paste0("TPR diff exp, DE50"), cex.main=1.5)
boxplot(DEDD.results.DEDD2$tpr, names=NA, col=col_vector[1:3], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.68), main=paste0("TPR diff exp, DEDD2"), cex.main=1.5)
boxplot(DEDD.results.DEDD5$tpr, names=NA, col=col_vector[1:3], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.68), main=paste0("TPR diff exp, DEDD5"), cex.main=1.5)
boxplot(DEDD.results.DEDD10$tpr, names=NA, col=col_vector[1:3], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.68), main=paste0("TPR diff exp, DEDD10"), cex.main=1.5)
boxplot(DEDD.results.DEDD20$tpr, names=NA, col=col_vector[1:3], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.68), main=paste0("TPR diff exp, DEDD20"), cex.main=1.5)
boxplot(DEDD.results.DEDD50$tpr, names=NA, col=col_vector[1:3], pch=20, cex.axis=1.2, 
        xaxt='n', ylim=c(0,0.68), main=paste0("TPR diff exp, DEDD50"), cex.main=1.5)
rbind(colMeans(DEDD.results.DD2$tpr), colMeans(DEDD.results.DD5$tpr), colMeans(DEDD.results.DD10$tpr), 
      colMeans(DEDD.results.DD20$tpr), colMeans(DEDD.results.DD50$tpr))
rbind(colMeans(DEDD.results.DE2$tpr), colMeans(DEDD.results.DE5$tpr), colMeans(DEDD.results.DE10$tpr), 
      colMeans(DEDD.results.DE20$tpr), colMeans(DEDD.results.DE50$tpr))
rbind(colMeans(DEDD.results.DEDD2$tpr), colMeans(DEDD.results.DEDD5$tpr), colMeans(DEDD.results.DEDD10$tpr), 
      colMeans(DEDD.results.DEDD20$tpr), colMeans(DEDD.results.DEDD50$tpr))
# lnHMM always gives higher TPRs than expHMM. Posterior threshold gives higher power than other methods in all situations except for 
# 50 samples per group with expHMM with differences in mean only or mean and dispersion, where BFDR is slightly higher, but lnHMM 
# with posterior threshold is always higher than expHMM with any method. BFDR is far worse than the others for 2, 50 and 10 samples 
# per group, slightly worse for 20 samples per group, and better than 0.5 posterior probability threshold for 50 samples per group.

#############################
### False discovery plots ###
par(mfrow=c(3,5), mar=c(2,2,2,1), mgp=c(3,0.7,0))

plot(DEDD.results.DD2$mean.discoveries$lnHMM.tmm, DEDD.results.DD2$mean.fdr$lnHMM.tmm, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], main="DD2", cex.main=1.5)
lines(DEDD.results.DD2$mean.discoveries$expHMM.tmm, DEDD.results.DD2$mean.fdr$expHMM.tmm, 
      type='l', col=col_vector[3])
abline(h=0.05, col='lightgrey')
plot(DEDD.results.DD5$mean.discoveries$lnHMM.tmm, DEDD.results.DD5$mean.fdr$lnHMM.tmm, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], main="DD5", cex.main=1.5)
lines(DEDD.results.DD5$mean.discoveries$expHMM.tmm, DEDD.results.DD5$mean.fdr$expHMM.tmm, 
      type='l', col=col_vector[3])
abline(h=0.05, col='lightgrey')
plot(DEDD.results.DD10$mean.discoveries$lnHMM.tmm, DEDD.results.DD10$mean.fdr$lnHMM.tmm, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], main="DD10", cex.main=1.5)
lines(DEDD.results.DD10$mean.discoveries$expHMM.tmm, DEDD.results.DD10$mean.fdr$expHMM.tmm, 
      type='l', col=col_vector[3])
abline(h=0.05, col='lightgrey')
plot(DEDD.results.DD20$mean.discoveries$lnHMM.tmm, DEDD.results.DD20$mean.fdr$lnHMM.tmm, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], main="DD20", cex.main=1.5)
lines(DEDD.results.DD20$mean.discoveries$expHMM.tmm, DEDD.results.DD20$mean.fdr$expHMM.tmm, 
      type='l', col=col_vector[3])
abline(h=0.05, col='lightgrey')
plot(DEDD.results.DD50$mean.discoveries$lnHMM.tmm, DEDD.results.DD50$mean.fdr$lnHMM.tmm, 
     type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1], main="DD50", cex.main=1.5)
lines(DEDD.results.DD50$mean.discoveries$expHMM.tmm, DEDD.results.DD50$mean.fdr$expHMM.tmm, 
      type='l', col=col_vector[3])
abline(h=0.05, col='lightgrey')
legend("topright", bty='n', col=col_vector[c(1,3)], lty=1, legend=legend, cex=1.5)

plot(DEDD.results.DE2$mean.discoveries$lnHMM.tmm, DEDD.results.DE2$mean.fdr$lnHMM.tmm, 
     type='l', ylim=c(0,0.9), xlim=c(0,2000), col=col_vector[1], main="DE2", cex.main=1.5)
lines(DEDD.results.DE2$mean.discoveries$expHMM.tmm, DEDD.results.DE2$mean.fdr$expHMM.tmm, 
      type='l', col=col_vector[3])
abline(h=0.05, col='lightgrey')
plot(DEDD.results.DE5$mean.discoveries$lnHMM.tmm, DEDD.results.DE5$mean.fdr$lnHMM.tmm, 
     type='l', ylim=c(0,0.9), xlim=c(0,2000), col=col_vector[1], main="DE5", cex.main=1.5)
lines(DEDD.results.DE5$mean.discoveries$expHMM.tmm, DEDD.results.DE5$mean.fdr$expHMM.tmm, 
      type='l', col=col_vector[3])
abline(h=0.05, col='lightgrey')
plot(DEDD.results.DE10$mean.discoveries$lnHMM.tmm, DEDD.results.DE10$mean.fdr$lnHMM.tmm, 
     type='l', ylim=c(0,0.9), xlim=c(0,2000), col=col_vector[1], main="DE10", cex.main=1.5)
lines(DEDD.results.DE10$mean.discoveries$expHMM.tmm, DEDD.results.DE10$mean.fdr$expHMM.tmm, 
      type='l', col=col_vector[3])
abline(h=0.05, col='lightgrey')
plot(DEDD.results.DE20$mean.discoveries$lnHMM.tmm, DEDD.results.DE20$mean.fdr$lnHMM.tmm, 
     type='l', ylim=c(0,0.9), xlim=c(0,2000), col=col_vector[1], main="DE20", cex.main=1.5)
lines(DEDD.results.DE20$mean.discoveries$expHMM.tmm, DEDD.results.DE20$mean.fdr$expHMM.tmm, 
      type='l', col=col_vector[3])
abline(h=0.05, col='lightgrey')
plot(DEDD.results.DE50$mean.discoveries$lnHMM.tmm, DEDD.results.DE50$mean.fdr$lnHMM.tmm, 
     type='l', ylim=c(0,0.9), xlim=c(0,2000), col=col_vector[1], main="DE50", cex.main=1.5)
lines(DEDD.results.DE50$mean.discoveries$expHMM.tmm, DEDD.results.DE50$mean.fdr$expHMM.tmm, 
      type='l', col=col_vector[3])
abline(h=0.05, col='lightgrey')
legend("topright", bty='n', col=col_vector[c(1,3)], lty=1, legend=legend, cex=1.5)

plot(DEDD.results.DEDD2$mean.discoveries$expHMM.tmm, DEDD.results.DEDD2$mean.fdr$expHMM.tmm, 
     type='l', ylim=c(0,0.65), xlim=c(0,2000), col=col_vector[1], main="DEDD2", cex.main=1.5)
lines(DEDD.results.DEDD2$mean.discoveries$lnHMM.tmm, DEDD.results.DEDD2$mean.fdr$lnHMM.tmm, 
      type='l', col=col_vector[3])
abline(h=0.05, col='lightgrey')
plot(DEDD.results.DEDD5$mean.discoveries$expHMM.tmm, DEDD.results.DEDD5$mean.fdr$expHMM.tmm, 
     type='l', ylim=c(0,0.65), xlim=c(0,2000), col=col_vector[1], main="DEDD5", cex.main=1.5)
lines(DEDD.results.DEDD5$mean.discoveries$lnHMM.tmm, DEDD.results.DEDD5$mean.fdr$lnHMM.tmm, 
      type='l', col=col_vector[3])
abline(h=0.05, col='lightgrey')
plot(DEDD.results.DEDD10$mean.discoveries$expHMM.tmm, DEDD.results.DEDD10$mean.fdr$expHMM.tmm, 
     type='l', ylim=c(0,0.65), xlim=c(0,2000), col=col_vector[1], main="DEDD10", cex.main=1.5)
lines(DEDD.results.DEDD10$mean.discoveries$lnHMM.tmm, DEDD.results.DEDD10$mean.fdr$lnHMM.tmm, 
      type='l', col=col_vector[3])
abline(h=0.05, col='lightgrey')
plot(DEDD.results.DEDD20$mean.discoveries$expHMM.tmm, DEDD.results.DEDD20$mean.fdr$expHMM.tmm, 
     type='l', ylim=c(0,0.65), xlim=c(0,2000), col=col_vector[1], main="DEDD20", cex.main=1.5)
lines(DEDD.results.DEDD20$mean.discoveries$lnHMM.tmm, DEDD.results.DEDD20$mean.fdr$lnHMM.tmm, 
      type='l', col=col_vector[3])
abline(h=0.05, col='lightgrey')
plot(DEDD.results.DEDD50$mean.discoveries$expHMM.tmm, DEDD.results.DEDD50$mean.fdr$expHMM.tmm, 
     type='l', ylim=c(0,0.65), xlim=c(0,2000), col=col_vector[1], main="DEDD50", cex.main=1.5)
lines(DEDD.results.DEDD50$mean.discoveries$lnHMM.tmm, DEDD.results.DEDD50$mean.fdr$lnHMM.tmm, 
      type='l', col=col_vector[3])
abline(h=0.05, col='lightgrey')
legend("topright", bty='n', col=col_vector[c(1,3)], lty=1, legend=legend, cex=1.5)
# FDRs almost always lower for lnHMM than expHMM. Very little difference at low end in all cases except with 2 samples per group. At 
# higher end, lnHMM better than expHMM with differences in dispersion only or mean only, but expHMM better than lnHMM with 
# differences in mean and dispersion.
