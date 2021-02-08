library(here)
library(RColorBrewer)

## Load data, create colour vector ####
folder <- "Results/GTEx DD lfc1 results Jan 2021"
for (i in c("2", "5", "10", "20", "50")) {
  for (j in c("blood", "muscle")) {
    assign(paste0("DD.lfc1.results.", j, ".DEDD", i), 
           readRDS(here(folder, paste0("DD.lfc1.results.", j, "_", i, "_DEDD.rds"))))
  }
}
rm(i,j,folder)

qual_col_pals = brewer.pal.info[brewer.pal.info$category == "qual" & 
                                  brewer.pal.info$colorblind == T,]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector <- col_vector[c(1:3,9,11,26)]
rm(qual_col_pals)


positions = rep(1:10)
offsets <- c(rep(0, 2), rep(1, 2), rep(2, 2), rep(3, 2), rep(4, 2))

# Blood
par(mfrow=c(2,1), mar=c(2,2,2,0.5), mgp=c(3,0.7,0))
plot(NA, xlim=c(0.75,14.25), ylim=c(0,1), xaxt="n", yaxt="n", 
     main="False discovery rate", cex.main=1.5)
lines(c(0.25,2.75), c(0.05,0.05), col="lightgrey")
lines(c(3.25,5.75), c(0.05,0.05), col="lightgrey")
lines(c(6.25,8.75), c(0.05,0.05), col="lightgrey")
lines(c(9.25,11.75), c(0.05,0.05), col="lightgrey")
lines(c(12.25,14.75), c(0.05,0.05), col="lightgrey")
boxplot(cbind(DD.lfc1.results.blood.DEDD2$fdr, 
              DD.lfc1.results.blood.DEDD5$fdr, 
              DD.lfc1.results.blood.DEDD10$fdr, 
              DD.lfc1.results.blood.DEDD20$fdr, 
              DD.lfc1.results.blood.DEDD50$fdr), 
        col=col_vector[1:2], pch=20, cex.axis=1.5, xaxt='n', 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(1.5,4.5,7.5,10.5,13.5), labels=c(2,5,10,20,50), tick=F, cex.axis=1.5)
plot(NA, xlim=c(0.75,14.25), ylim=c(0,0.43), xaxt="n", yaxt="n", 
     main="Sensitivity", cex.main=1.5)
boxplot(cbind(DD.lfc1.results.blood.DEDD2$tpr, 
              DD.lfc1.results.blood.DEDD5$tpr, 
              DD.lfc1.results.blood.DEDD10$tpr, 
              DD.lfc1.results.blood.DEDD20$tpr, 
              DD.lfc1.results.blood.DEDD50$tpr), 
        col=col_vector[1:2], pch=20, cex.axis=1.5, xaxt='n', 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(1.5,4.5,7.5,10.5,13.5), labels=c(2,5,10,20,50), tick=F, cex.axis=1.5)
legend("topleft", fill=col_vector[1:3], bty='n', cex=1.5, ncol=1, 
       legend=c("MDSeq", "HM"))

# Muscle
par(mfrow=c(2,1), mar=c(2,2,2,0.5), mgp=c(3,0.7,0))
plot(NA, xlim=c(0.75,14.25), ylim=c(0,1), xaxt="n", yaxt="n", 
     main="False discovery rate", cex.main=1.5)
lines(c(0.25,2.75), c(0.05,0.05), col="lightgrey")
lines(c(3.25,5.75), c(0.05,0.05), col="lightgrey")
lines(c(6.25,8.75), c(0.05,0.05), col="lightgrey")
lines(c(9.25,11.75), c(0.05,0.05), col="lightgrey")
lines(c(12.25,14.75), c(0.05,0.05), col="lightgrey")
boxplot(cbind(DD.lfc1.results.muscle.DEDD2$fdr, 
              DD.lfc1.results.muscle.DEDD5$fdr, 
              DD.lfc1.results.muscle.DEDD10$fdr, 
              DD.lfc1.results.muscle.DEDD20$fdr, 
              DD.lfc1.results.muscle.DEDD50$fdr), 
        col=col_vector[1:2], pch=20, cex.axis=1.5, xaxt='n', 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(1.5,4.5,7.5,10.5,13.5), labels=c(2,5,10,20,50), tick=F, cex.axis=1.5)
plot(NA, xlim=c(0.75,14.25), ylim=c(0,0.43), xaxt="n", yaxt="n", 
     main="Sensitivity", cex.main=1.5)
boxplot(cbind(DD.lfc1.results.muscle.DEDD2$tpr, 
              DD.lfc1.results.muscle.DEDD5$tpr, 
              DD.lfc1.results.muscle.DEDD10$tpr, 
              DD.lfc1.results.muscle.DEDD20$tpr, 
              DD.lfc1.results.muscle.DEDD50$tpr), 
        col=col_vector[1:2], pch=20, cex.axis=1.5, xaxt='n', 
        at=positions + offsets, lwd=0.5, cex=0.5, 
        add=T)
axis(side=1, at=c(1.5,4.5,7.5,10.5,13.5), labels=c(2,5,10,20,50), tick=F, cex.axis=1.5)
legend("topleft", fill=col_vector[1:3], bty='n', cex=1.5, ncol=1, 
       legend=c("MDSeq", "HM"))

