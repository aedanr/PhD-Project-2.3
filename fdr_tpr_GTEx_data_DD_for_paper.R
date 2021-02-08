library(here)
library(RColorBrewer)

## Load data, create colour vector ####
folder <- "Results/GTEx combined results Dec 2020"
for (i in c("2", "5", "10", "20", "50")) {
  assign(paste0("DD.results.DEDD", i), 
         readRDS(here(folder, paste0("DD.results.DEDD", i, ".rds"))))
}
rm(i,folder)

qual_col_pals = brewer.pal.info[brewer.pal.info$category == "qual" & 
                                  brewer.pal.info$colorblind == T,]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# pie(rep(1,length(col_vector)), col=col_vector)
col_vector <- col_vector[c(1:3,9,11,26)]
rm(qual_col_pals)


## Keep only one version of each method: MDSeq ZI, lnHM log.
## Separate objects for each metric and for blood and muscle
for (metric in c("fdr", "tpr")) {
  for (size in c("2", "5", "10", "20", "50")) {
    assign(
      paste0(metric, "_blood_", size), 
      get(paste0("DD.results.DEDD", size))[[metric]][
        , c("blood_MD.zi", "blood_EVQ", "blood_lnHM.log")
      ]
    )
    assign(
      paste0(metric, "_muscle_", size), 
      get(paste0("DD.results.DEDD", size))[[metric]][
        , c("muscle_MD.zi", "muscle_EVQ", "muscle_lnHM.log")
      ]
    )
  }
}
rm(metric,size)

positions = rep(1:15)
offsets <- c(rep(0, 3), rep(1, 3), rep(2, 3), rep(3, 3), rep(4, 3))

### FDR and TPR blood ####
par(mfrow=c(2,1), mar=c(2,2,2,0.5), mgp=c(3,0.7,0))
plot(NA, xlim=c(0.75,19.25), ylim=c(0,0.93), xaxt="n", yaxt="n", 
     main="False discovery rate", cex.main=1.5)
lines(c(0.25,3.75), c(0.05,0.05), col="lightgrey")
lines(c(4.25,7.75), c(0.05,0.05), col="lightgrey")
lines(c(8.25,11.75), c(0.05,0.05), col="lightgrey")
lines(c(12.25,15.75), c(0.05,0.05), col="lightgrey")
lines(c(16.25,19.75), c(0.05,0.05), col="lightgrey")
boxplot(cbind(fdr_blood_2, fdr_blood_5, fdr_blood_10, fdr_blood_20, fdr_blood_50), 
        col=col_vector[1:3], pch=20, cex.axis=1.5, xaxt='n', 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(2,6,10,14,18), labels=c(2,5,10,20,50), tick=F, cex.axis=1.5)
plot(NA, xlim=c(0.75,19.25), ylim=c(0,0.65), xaxt="n", yaxt="n", 
     main="Sensitivity", cex.main=1.5)
boxplot(cbind(tpr_blood_2, tpr_blood_5, tpr_blood_10, tpr_blood_20, tpr_blood_50), 
        col=col_vector[1:3], pch=20, cex.axis=1.5, xaxt='n', 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(2,6,10,14,18), labels=c(2,5,10,20,50), tick=F, cex.axis=1.5)
legend("topleft", fill=col_vector[1:3], bty='n', cex=1.5, ncol=1, 
       legend=c("MDSeq", "GAMLSS", "HM"))

rbind(colMeans(cbind(fdr_blood_2)), colMeans(cbind(fdr_blood_5)), colMeans(cbind(fdr_blood_10)), 
      colMeans(cbind(fdr_blood_20)), colMeans(cbind(fdr_blood_50)))
rbind(colMeans(cbind(tpr_blood_2)), colMeans(cbind(tpr_blood_5)), colMeans(cbind(tpr_blood_10)), 
      colMeans(cbind(tpr_blood_20)), colMeans(cbind(tpr_blood_50)))


### FDR and TPR muscle ####
par(mfrow=c(2,1), mar=c(2,2,2,0.5), mgp=c(3,0.7,0))
plot(NA, xlim=c(0.75,19.25), ylim=c(0,0.93), xaxt="n", yaxt="n", 
     main="False discovery rate", cex.main=1.5)
lines(c(0.25,3.75), c(0.05,0.05), col="lightgrey")
lines(c(4.25,7.75), c(0.05,0.05), col="lightgrey")
lines(c(8.25,11.75), c(0.05,0.05), col="lightgrey")
lines(c(12.25,15.75), c(0.05,0.05), col="lightgrey")
lines(c(16.25,19.75), c(0.05,0.05), col="lightgrey")
boxplot(cbind(fdr_muscle_2, fdr_muscle_5, fdr_muscle_10, fdr_muscle_20, fdr_muscle_50), 
        col=col_vector[1:3], pch=20, cex.axis=1.5, xaxt='n', 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(2,6,10,14,18), labels=c(2,5,10,20,50), tick=F, cex.axis=1.5)
plot(NA, xlim=c(0.75,19.25), ylim=c(0,0.6), xaxt="n", yaxt="n", 
     main="Sensitivity", cex.main=1.5)
boxplot(cbind(tpr_muscle_2, tpr_muscle_5, tpr_muscle_10, tpr_muscle_20, tpr_muscle_50), 
        col=col_vector[1:3], pch=20, cex.axis=1.5, xaxt='n', 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(2,6,10,14,18), labels=c(2,5,10,20,50), tick=F, cex.axis=1.5)
legend("topleft", fill=col_vector[1:3], bty='n', cex=1.5, ncol=1, 
       legend=c("MDSeq", "GAMLSS", "HM"))

rbind(colMeans(cbind(fdr_muscle_2)), colMeans(cbind(fdr_muscle_5)), colMeans(cbind(fdr_muscle_10)), 
      colMeans(cbind(fdr_muscle_20)), colMeans(cbind(fdr_muscle_50)))
rbind(colMeans(cbind(tpr_muscle_2)), colMeans(cbind(tpr_muscle_5)), colMeans(cbind(tpr_muscle_10)), 
      colMeans(cbind(tpr_muscle_20)), colMeans(cbind(tpr_muscle_50)))

