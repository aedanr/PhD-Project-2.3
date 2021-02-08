library(here)
library(RColorBrewer)

## Load data, create colour vector ####
folder <- "Results/GTEx combined results Sept 2020"
for (i in c("2", "5", "10", "20", "50")) {
  assign(paste0("DEDD.results.DEDD", i), 
         readRDS(here(folder, paste0("DEDD.results.DEDD", i, ".rds"))))
}
rm(i,folder)

qual_col_pals = brewer.pal.info[brewer.pal.info$category == "qual" & 
                                  brewer.pal.info$colorblind == T,]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector <- col_vector[c(1:3,9,11,26)]
rm(qual_col_pals)


## Keep only diffVar, voom/lnHM.log, lnHMM (posterior threshold and BFDR).
## Separate objects for each metric and for blood and muscle
for (metric in c("fdr", "tpr")) {
  for (size in c("2", "5", "10", "20", "50")) {
    assign(
      paste0(metric, "_blood_", size), 
      get(paste0("DEDD.results.DEDD", size))[[metric]][
        , c("blood_dV", "blood_voom_HM", "blood_lnHMM.thr", "blood_lnHMM.bfdr")
      ]
    )
    assign(
      paste0(metric, "_muscle_", size), 
      get(paste0("DEDD.results.DEDD", size))[[metric]][
        , c("muscle_dV", "muscle_voom_HM", "muscle_lnHMM.thr", "muscle_lnHMM.bfdr")
      ]
    )
  }
}


positions = rep(1:20)
offsets <- c(rep(0, 4), rep(2, 4), rep(4, 4), rep(6, 4), rep(8, 4))

### FDR and TPR blood ####
par(mfrow=c(2,1), mar=c(2,2,2,0.5), mgp=c(3,0.7,0))
plot(NA, xlim=c(1,28), ylim=c(0,1), xaxt="n", yaxt="n", 
     main="False discovery rate", cex.main=1.5)
lines(c(0.5,4.5), c(0.05,0.05), col="lightgrey")
lines(c(6.5,10.5), c(0.05,0.05), col="lightgrey")
lines(c(12.5,16.5), c(0.05,0.05), col="lightgrey")
lines(c(18.5,22.5), c(0.05,0.05), col="lightgrey")
lines(c(24.5,28.5), c(0.05,0.05), col="lightgrey")
boxplot(cbind(fdr_blood_2, fdr_blood_5, fdr_blood_10, fdr_blood_20, fdr_blood_50), 
        col=col_vector[1:4], pch=20, cex.axis=1.5, xaxt='n', 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(2.5,8.5,14.5,20.5,26.5), labels=c(2,5,10,20,50), tick=F, cex.axis=1.5)
legend("topright", fill=col_vector[1:6], bty='n', cex=1.5, ncol=1, 
       legend=c("diffVar", "Hybrid", "HMM, posterior threshold", "HMM, BFDR"))
plot(NA, xlim=c(1,28), ylim=c(0,1), xaxt="n", yaxt="n", 
     main="Sensitivity", cex.main=1.5)
boxplot(cbind(tpr_blood_2, tpr_blood_5, tpr_blood_10, tpr_blood_20, tpr_blood_50), 
        col=col_vector[1:4], pch=20, cex.axis=1.5, xaxt='n', 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(2.5,8.5,14.5,20.5,26.5), labels=c(2,5,10,20,50), tick=F, cex.axis=1.5)

rbind(colMeans(cbind(fdr_blood_2)), colMeans(cbind(fdr_blood_5)), colMeans(cbind(fdr_blood_10)), 
      colMeans(cbind(fdr_blood_20)), colMeans(cbind(fdr_blood_50)))
rbind(colMeans(cbind(tpr_blood_2)), colMeans(cbind(tpr_blood_5)), colMeans(cbind(tpr_blood_10)), 
      colMeans(cbind(tpr_blood_20)), colMeans(cbind(tpr_blood_50)))


### FDR and TPR muscle ####
par(mfrow=c(2,1), mar=c(2,2,2,0.5), mgp=c(3,0.7,0))
plot(NA, xlim=c(1,28), ylim=c(0,1), xaxt="n", yaxt="n", 
     main="False discovery rate", cex.main=1.5)
lines(c(0.5,4.5), c(0.05,0.05), col="lightgrey")
lines(c(6.5,10.5), c(0.05,0.05), col="lightgrey")
lines(c(12.5,16.5), c(0.05,0.05), col="lightgrey")
lines(c(18.5,22.5), c(0.05,0.05), col="lightgrey")
lines(c(24.5,28.5), c(0.05,0.05), col="lightgrey")
boxplot(cbind(fdr_muscle_2, fdr_muscle_5, fdr_muscle_10, fdr_muscle_20, fdr_muscle_50), 
        col=col_vector[1:4], pch=20, cex.axis=1.5, xaxt='n', 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(2.5,8.5,14.5,20.5,26.5), labels=c(2,5,10,20,50), tick=F, cex.axis=1.5)
legend("topright", fill=col_vector[1:6], bty='n', cex=1.5, ncol=1, 
       legend=c("diffVar", "Hybrid", "HMM, posterior threshold", "HMM, BFDR"))
plot(NA, xlim=c(1,28), ylim=c(0,1), xaxt="n", yaxt="n", 
     main="Sensitivity", cex.main=1.5)
boxplot(cbind(tpr_muscle_2, tpr_muscle_5, tpr_muscle_10, tpr_muscle_20, tpr_muscle_50), 
        col=col_vector[1:4], pch=20, cex.axis=1.5, xaxt='n', 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(2.5,8.5,14.5,20.5,26.5), labels=c(2,5,10,20,50), tick=F, cex.axis=1.5)

rbind(colMeans(cbind(fdr_muscle_2)), colMeans(cbind(fdr_muscle_5)), colMeans(cbind(fdr_muscle_10)), 
      colMeans(cbind(fdr_muscle_20)), colMeans(cbind(fdr_muscle_50)))
rbind(colMeans(cbind(tpr_muscle_2)), colMeans(cbind(tpr_muscle_5)), colMeans(cbind(tpr_muscle_10)), 
      colMeans(cbind(tpr_muscle_20)), colMeans(cbind(tpr_muscle_50)))


