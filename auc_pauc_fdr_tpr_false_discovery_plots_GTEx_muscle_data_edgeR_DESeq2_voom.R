library(here)
library(RColorBrewer)

## Load data, create colour vector ####
folder <- "Results/GTEx muscle simulated DE, DD, DEDD results Feb 2020/edgeR, DESeq2, voom only"
for (j in c("2_DE", "5_DE", "10_DE", "20_DE", "50_DE")) {
  assign(paste0("DE.results.muscle_", j), 
         readRDS(here(folder, paste0("DE.results.muscle_", j, ".rds"))))
}
rm(j)

n <- 23
qual_col_pals = brewer.pal.info[brewer.pal.info$category == "qual" & 
                                  brewer.pal.info$colorblind == T,]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, 
                           rownames(qual_col_pals)))
# pie(rep(1,n), col=col_vector)
col_vector <- col_vector[-c(7,10,11,12,20)]


#################################
#### Differential expression ####
#################################

##############
#### AUC #####
par(mfrow=c(2,3), mar=c(0.5,2,5,1), mgp=c(3,0.7,0))
boxplot(DE.results.muscle_2_DE$auc, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0.65,1.1), 
        main=paste0("AUC diff exp, 2 samples per group\n", 
                    "Differences in mean only"), 
        cex.main=1.7)
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", 
                "DESeq2 no IF", "DESeq2 IF", "voom"))
boxplot(DE.results.muscle_5_DE$auc, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0.65,1.1), 
        main=paste0("AUC diff exp, 5 samples per group\n", 
                    "Differences in mean only"), 
        cex.main=1.7)
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", 
                "DESeq2 no IF", "DESeq2 IF", "voom"))
boxplot(DE.results.muscle_10_DE$auc, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0.65,1.1), 
        main=paste0("AUC diff exp, 10 samples per group\n", 
                    "Differences in mean only"), 
        cex.main=1.7)
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", 
                "DESeq2 no IF", "DESeq2 IF", "voom"))
boxplot(DE.results.muscle_20_DE$auc, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0.65,1.1), 
        main=paste0("AUC diff exp, 20 samples per group\n", 
                    "Differences in mean only"), 
        cex.main=1.7)
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", 
                "DESeq2 no IF", "DESeq2 IF", "voom"))
boxplot(DE.results.muscle_50_DE$auc, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0.65,1.1), 
        main=paste0("AUC diff exp, 50 samples per group\n", 
                    "Differences in mean only"), 
        cex.main=1.7)
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.2, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", 
                "DESeq2 no IF", "DESeq2 IF", "voom"))
rbind(colMeans(DE.results.muscle_2_DE$auc), 
      colMeans(DE.results.muscle_5_DE$auc), 
      colMeans(DE.results.muscle_10_DE$auc), 
      colMeans(DE.results.muscle_20_DE$auc), 
      colMeans(DE.results.muscle_50_DE$auc))


##############
#### pAUC ####
par(mfrow=c(2,3), mar=c(0.5,2,5,1), mgp=c(3,0.7,0))
boxplot(DE.results.muscle_2_DE$pauc, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,0.065), 
        main=paste0("Partial AUC diff exp, 2 samples per group\n", 
                    "Differences in mean only"), 
        cex.main=1.7)
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.4, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", 
                "DESeq2 no IF", "DESeq2 IF", "voom"))
boxplot(DE.results.muscle_5_DE$pauc, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,0.065), 
        main=paste0("Partial AUC diff exp, 5 samples per group\n", 
                    "Differences in mean only"), 
        cex.main=1.7)
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.4, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", 
                "DESeq2 no IF", "DESeq2 IF", "voom"))
boxplot(DE.results.muscle_10_DE$pauc, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,0.065), 
        main=paste0("Partial AUC diff exp, 10 samples per group\n", 
                    "Differences in mean only"), 
        cex.main=1.7)
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.4, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", 
                "DESeq2 no IF", "DESeq2 IF", "voom"))
boxplot(DE.results.muscle_20_DE$pauc, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,0.065), 
        main=paste0("Partial AUC diff exp, 20 samples per group\n", 
                    "Differences in mean only"), 
        cex.main=1.7)
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.4, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", 
                "DESeq2 no IF", "DESeq2 IF", "voom"))
boxplot(DE.results.muscle_50_DE$pauc, names=NA, col=col_vector[1:6], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,0.065), 
        main=paste0("Partial AUC diff exp, 50 samples per group\n", 
                    "Differences in mean only"), 
        cex.main=1.7)
legend("topleft", fill=col_vector[1:6], bty='n', cex=1.4, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", 
                "DESeq2 no IF", "DESeq2 IF", "voom"))
rbind(colMeans(DE.results.muscle_2_DE$pauc), 
      colMeans(DE.results.muscle_5_DE$pauc), 
      colMeans(DE.results.muscle_10_DE$pauc), 
      colMeans(DE.results.muscle_20_DE$pauc), 
      colMeans(DE.results.muscle_50_DE$pauc))


#############
#### FDR ####
par(mfrow=c(2,3), mar=c(0.5,2,5,1), mgp=c(3,0.7,0))
boxplot(DE.results.muscle_2_DE$fdr, names=NA, col=col_vector[1:15], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.3), 
        main=paste0("FDR diff exp, 2 samples per group\n", 
                    "Differences in mean only"), 
        cex.main=1.7)
abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[1:15], bty='n', cex=1.4, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", 
                "DESeq2 no IF", "DESeq2 IF", "voom"))

boxplot(DE.results.muscle_5_DE$fdr, names=NA, col=col_vector[1:15], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.3), 
        main=paste0("FDR diff exp, 5 samples per group\n", 
                    "Differences in mean only"), 
        cex.main=1.7)
abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[1:15], bty='n', cex=1.4, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", 
                "DESeq2 no IF", "DESeq2 IF", "voom"))

boxplot(DE.results.muscle_10_DE$fdr, names=NA, col=col_vector[1:15], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.3), 
        main=paste0("FDR diff exp, 10 samples per group\n", 
                    "Differences in mean only"), 
        cex.main=1.7)
abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[1:15], bty='n', cex=1.4, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", 
                "DESeq2 no IF", "DESeq2 IF", "voom"))

boxplot(DE.results.muscle_20_DE$fdr, names=NA, col=col_vector[1:15], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.3), 
        main=paste0("FDR diff exp, 20 samples per group\n", 
                    "Differences in mean only"), 
        cex.main=1.7)
abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[1:15], bty='n', cex=1.4, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", 
                "DESeq2 no IF", "DESeq2 IF", "voom"))

boxplot(DE.results.muscle_50_DE$fdr, names=NA, col=col_vector[1:15], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.3), 
        main=paste0("FDR diff exp, 50 samples per group\n", 
                    "Differences in mean only"), 
        cex.main=1.7)
abline(h=0.05, col="lightgrey")
legend("topleft", fill=col_vector[1:15], bty='n', cex=1.4, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", 
                "DESeq2 no IF", "DESeq2 IF", "voom"))
rbind(colMeans(DE.results.muscle_2_DE$fdr), 
      colMeans(DE.results.muscle_5_DE$fdr), 
      colMeans(DE.results.muscle_10_DE$fdr), 
      colMeans(DE.results.muscle_20_DE$fdr), 
      colMeans(DE.results.muscle_50_DE$fdr))


#############
#### TPR ####
par(mfrow=c(2,3), mar=c(0.5,2,5,1), mgp=c(3,0.7,0))
boxplot(DE.results.muscle_2_DE$tpr, names=NA, col=col_vector[1:15], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.25), 
        main=paste0("TPR diff exp, 2 samples per group\n", 
                    "Differences in mean only"), 
        cex.main=1.7)
legend("topleft", fill=col_vector[1:15], bty='n', cex=1.4, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", 
                "DESeq2 no IF", "DESeq2 IF", "voom"))
boxplot(DE.results.muscle_5_DE$tpr, names=NA, col=col_vector[1:15], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.25), 
        main=paste0("TPR diff exp, 5 samples per group\n", 
                    "Differences in mean only"), 
        cex.main=1.7)
legend("topleft", fill=col_vector[1:15], bty='n', cex=1.4, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", 
                "DESeq2 no IF", "DESeq2 IF", "voom"))
boxplot(DE.results.muscle_10_DE$tpr, names=NA, col=col_vector[1:15], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.25), 
        main=paste0("TPR diff exp, 10 samples per group\n", 
                    "Differences in mean only"), 
        cex.main=1.7)
legend("topleft", fill=col_vector[1:15], bty='n', cex=1.4, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", 
                "DESeq2 no IF", "DESeq2 IF", "voom"))
boxplot(DE.results.muscle_20_DE$tpr, names=NA, col=col_vector[1:15], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.25), 
        main=paste0("TPR diff exp, 20 samples per group\n", 
                    "Differences in mean only"), 
        cex.main=1.7)
legend("topleft", fill=col_vector[1:15], bty='n', cex=1.4, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", 
                "DESeq2 no IF", "DESeq2 IF", "voom"))
boxplot(DE.results.muscle_50_DE$tpr, names=NA, col=col_vector[1:15], 
        pch=20, cex.axis=1.7, xaxt='n', ylim=c(0,1.25), 
        main=paste0("TPR diff exp, 50 samples per group\n", 
                    "Differences in mean only"), 
        cex.main=1.7)
legend("topleft", fill=col_vector[1:15], bty='n', cex=1.4, ncol=3, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", 
                "DESeq2 no IF", "DESeq2 IF", "voom"))
rbind(colMeans(DE.results.muscle_2_DE$tpr), colMeans(DE.results.muscle_5_DE$tpr), 
      colMeans(DE.results.muscle_10_DE$tpr), colMeans(DE.results.muscle_20_DE$tpr), 
      colMeans(DE.results.muscle_50_DE$tpr))


###############################
#### False discovery plots ####
par(mfrow=c(2,3), mar=c(2,1.5,2.5,0.5), mgp=c(3,0.5,0))
plot(DE.results.muscle_2_DE$mean.discoveries$edgeR.ql.tmm, 
     DE.results.muscle_2_DE$mean.fdr$edgeR.ql.tmm, 
     type='l', ylim=c(0,0.7), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("2 samples per group\n", 
                 "Differences in mean only"))
lines(DE.results.muscle_2_DE$mean.discoveries$edgeR.lr.tmm, 
      DE.results.muscle_2_DE$mean.fdr$edgeR.lr.tmm, 
      type='l', col=col_vector[2])
lines(DE.results.muscle_2_DE$mean.discoveries$edgeR.et.tmm, 
      DE.results.muscle_2_DE$mean.fdr$edgeR.et.tmm, 
      type='l', col=col_vector[3])
lines(DE.results.muscle_2_DE$mean.discoveries$DESeq2.noif.tmm, 
      DE.results.muscle_2_DE$mean.fdr$DESeq2.noif.tmm, 
      type='l', col=col_vector[4])
lines(DE.results.muscle_2_DE$mean.discoveries$DESeq2.if.tmm, 
      DE.results.muscle_2_DE$mean.fdr$DESeq2.if.tmm, 
      type='l', col=col_vector[5])
lines(DE.results.muscle_2_DE$mean.discoveries$voom.tmm, 
      DE.results.muscle_2_DE$mean.fdr$voom.tmm, 
      type='l', col=col_vector[6])
abline(h=0.05, col='lightgrey')
legend("topright", bty='n', col=col_vector[1:6], lty=1, ncol=2, lwd=2, cex=1.2, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", 
                "DESeq2 no IF", "DESeq2 IF", "voom"))
plot(DE.results.muscle_5_DE$mean.discoveries$edgeR.ql.tmm, 
     DE.results.muscle_5_DE$mean.fdr$edgeR.ql.tmm, 
     type='l', ylim=c(0,0.7), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("5 samples per group\n", 
                 "Differences in mean only"))
lines(DE.results.muscle_5_DE$mean.discoveries$edgeR.lr.tmm, 
      DE.results.muscle_5_DE$mean.fdr$edgeR.lr.tmm, 
      type='l', col=col_vector[2])
lines(DE.results.muscle_5_DE$mean.discoveries$edgeR.et.tmm, 
      DE.results.muscle_5_DE$mean.fdr$edgeR.et.tmm, 
      type='l', col=col_vector[3])
lines(DE.results.muscle_5_DE$mean.discoveries$DESeq2.noif.tmm, 
      DE.results.muscle_5_DE$mean.fdr$DESeq2.noif.tmm, 
      type='l', col=col_vector[4])
lines(DE.results.muscle_5_DE$mean.discoveries$DESeq2.if.tmm, 
      DE.results.muscle_5_DE$mean.fdr$DESeq2.if.tmm, 
      type='l', col=col_vector[5])
lines(DE.results.muscle_5_DE$mean.discoveries$voom.tmm, 
      DE.results.muscle_5_DE$mean.fdr$voom.tmm, 
      type='l', col=col_vector[6])
abline(h=0.05, col='lightgrey')
plot(DE.results.muscle_10_DE$mean.discoveries$edgeR.ql.tmm, 
     DE.results.muscle_10_DE$mean.fdr$edgeR.ql.tmm, 
     type='l', ylim=c(0,0.7), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("10 samples per group\n", 
                 "Differences in mean only"))
lines(DE.results.muscle_10_DE$mean.discoveries$edgeR.lr.tmm, 
      DE.results.muscle_10_DE$mean.fdr$edgeR.lr.tmm, 
      type='l', col=col_vector[2])
lines(DE.results.muscle_10_DE$mean.discoveries$edgeR.et.tmm, 
      DE.results.muscle_10_DE$mean.fdr$edgeR.et.tmm, 
      type='l', col=col_vector[3])
lines(DE.results.muscle_10_DE$mean.discoveries$DESeq2.noif.tmm, 
      DE.results.muscle_10_DE$mean.fdr$DESeq2.noif.tmm, 
      type='l', col=col_vector[4])
lines(DE.results.muscle_10_DE$mean.discoveries$DESeq2.if.tmm, 
      DE.results.muscle_10_DE$mean.fdr$DESeq2.if.tmm, 
      type='l', col=col_vector[5])
lines(DE.results.muscle_10_DE$mean.discoveries$voom.tmm, 
      DE.results.muscle_10_DE$mean.fdr$voom.tmm, 
      type='l', col=col_vector[6])
abline(h=0.05, col='lightgrey')
plot(DE.results.muscle_20_DE$mean.discoveries$edgeR.ql.tmm, 
     DE.results.muscle_20_DE$mean.fdr$edgeR.ql.tmm, 
     type='l', ylim=c(0,0.7), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("20 samples per group\n", 
                 "Differences in mean only"))
lines(DE.results.muscle_20_DE$mean.discoveries$edgeR.lr.tmm, 
      DE.results.muscle_20_DE$mean.fdr$edgeR.lr.tmm, 
      type='l', col=col_vector[2])
lines(DE.results.muscle_20_DE$mean.discoveries$edgeR.et.tmm, 
      DE.results.muscle_20_DE$mean.fdr$edgeR.et.tmm, 
      type='l', col=col_vector[3])
lines(DE.results.muscle_20_DE$mean.discoveries$DESeq2.noif.tmm, 
      DE.results.muscle_20_DE$mean.fdr$DESeq2.noif.tmm, 
      type='l', col=col_vector[4])
lines(DE.results.muscle_20_DE$mean.discoveries$DESeq2.if.tmm, 
      DE.results.muscle_20_DE$mean.fdr$DESeq2.if.tmm, 
      type='l', col=col_vector[5])
lines(DE.results.muscle_20_DE$mean.discoveries$voom.tmm, 
      DE.results.muscle_20_DE$mean.fdr$voom.tmm, 
      type='l', col=col_vector[6])
abline(h=0.05, col='lightgrey')
plot(DE.results.muscle_50_DE$mean.discoveries$edgeR.ql.tmm, 
     DE.results.muscle_50_DE$mean.fdr$edgeR.ql.tmm, 
     type='l', ylim=c(0,0.7), xlim=c(0,2000), col=col_vector[1], 
     main=paste0("50 samples per group\n", 
                 "Differences in mean only"))
lines(DE.results.muscle_50_DE$mean.discoveries$edgeR.lr.tmm, 
      DE.results.muscle_50_DE$mean.fdr$edgeR.lr.tmm, 
      type='l', col=col_vector[2])
lines(DE.results.muscle_50_DE$mean.discoveries$edgeR.et.tmm, 
      DE.results.muscle_50_DE$mean.fdr$edgeR.et.tmm, 
      type='l', col=col_vector[3])
lines(DE.results.muscle_50_DE$mean.discoveries$DESeq2.noif.tmm, 
      DE.results.muscle_50_DE$mean.fdr$DESeq2.noif.tmm, 
      type='l', col=col_vector[4])
lines(DE.results.muscle_50_DE$mean.discoveries$DESeq2.if.tmm, 
      DE.results.muscle_50_DE$mean.fdr$DESeq2.if.tmm, 
      type='l', col=col_vector[5])
lines(DE.results.muscle_50_DE$mean.discoveries$voom.tmm, 
      DE.results.muscle_50_DE$mean.fdr$voom.tmm, 
      type='l', col=col_vector[6])
abline(h=0.05, col='lightgrey')



