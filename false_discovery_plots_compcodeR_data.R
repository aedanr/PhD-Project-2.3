library(here)
library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual' & 
                                  brewer.pal.info$colorblind == T,]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# n <- 15
# pie(rep(1,n), col=col_vector)
# rm(n)
col_vector <- col_vector[-c(1,5,7,10,12,15)]

for (i in c('DE', 'DD', 'DEDD')) {
  for (j in c('DEDD2', 'DEDD5', 'DEDD10', 'DEDD20')) {
    assign(paste0(i,'.results.TMM.',j), 
           readRDS(here('Results/DEDD compcodeR data results July-Aug 2019', 
                        paste0(i,'.results.',j,'.TMM.rds'))))
  }
}

for (i in c('DE', 'DEDD')) {
  for (j in c('DE2', 'DE5', 'DE10', 'DE20')) {
    assign(paste0(i,'.results.TMM.',j), 
           readRDS(here('Results/DE compcodeR data results July-Aug 2019', 
                        paste0(i,'.results.',j,'.TMM.rds'))))
  }
}

for (i in c('DD', 'DEDD')) {
  for (j in c('DD2', 'DD5', 'DD10', 'DD20')) {
    assign(paste0(i,'.results.TMM.',j), 
           readRDS(here('Results/DD compcodeR data results July-Aug 2019', 
                        paste0(i,'.results.',j,'.TMM.rds'))))
  }
}

for (i in c('DE', 'DD', 'DEDD')) {
  for (j in c('DEDD2', 'DEDD5', 'DEDD10', 'DEDD20')) {
    assign(paste0(i,'.results.DESeq.',j), 
           readRDS(here('Results/DEDD compcodeR data results DESeq norm Sept 2019', 
                        paste0(i,'.results.',j,'.DESeqnorm.rds'))))
  }
}

for (i in c('DE', 'DEDD')) {
  for (j in c('DE2', 'DE5', 'DE10', 'DE20')) {
    assign(paste0(i,'.results.DESeq.',j), 
           readRDS(here('Results/DE compcodeR data results DESeq norm Aug-Sept 2019', 
                        paste0(i,'.results.',j,'.DESeqnorm.rds'))))
  }
}

for (i in c('DD', 'DEDD')) {
  for (j in c('DD2', 'DD5', 'DD10', 'DD20')) {
    assign(paste0(i,'.results.DESeq.',j), 
           readRDS(here('Results/DD compcodeR data results DESeq norm Sept 2019', 
                        paste0(i,'.results.',j,'.DESeqnorm.rds'))))
  }
}
rm(i,j)


# DD2,5,10,20 TMM,DESeq DD: MDSeq(2), HMs(4)
# DD2,5,10,20 TMM,DESeq DEDD: HMMs(2)
# DE2,5,10,20 TMM,DESeq: edgeR(3), DESeq2(2), voom(1), DSS(1), baySeq(1), MDSeq(2), HMs(4)
# DE2,5,10,20 TMM,DESeq DEDD: HMMs(2)
# DEDD2,5,10,20 TMM,DESeq DD: MDSeq(2), HMs(4)
# DEDD2,5,10,20 TMM,DESeq DE: edgeR(3), DESeq2(2), voom(1), DSS(1), baySeq(1), MDSeq(2), HMs(4)
# DEDD2,5,10,20 TMM,DESeq DEDD: HMMs(2)


# Remove results not in all comparisons (ShrinkBayes, DSS trend), redundant versions (noif.DESeq) 
# and FDR versions that are clearly worse than others (BH for DSS, 0.5 threshold for baySeq).
dim(DD.results.TMM.DD20$mean.fdr)
# 20000 6
dim(DEDD.results.TMM.DD20$mean.fdr)
# 20000 2
DE.results.TMM.DE2$mean.fdr <- subset(DE.results.TMM.DE2$mean.fdr, 
                                      select = -c(trend.DSS, noif.DESeq))
DE.results.TMM.DE2$mean.discoveries <- subset(DE.results.TMM.DE2$mean.discoveries, 
                                              select = -c(trend.DSS, noif.DESeq))
DE.results.TMM.DE5$mean.fdr <- subset(DE.results.TMM.DE5$mean.fdr, 
                                      select = -c(trend.DSS, noif.DESeq))
DE.results.TMM.DE5$mean.discoveries <- subset(DE.results.TMM.DE5$mean.discoveries, 
                                              select = -c(trend.DSS, noif.DESeq))
DE.results.TMM.DE10$mean.fdr <- subset(DE.results.TMM.DE10$mean.fdr, 
                                       select = -noif.DESeq)
DE.results.TMM.DE10$mean.discoveries <- subset(DE.results.TMM.DE10$mean.discoveries, 
                                               select = -noif.DESeq)
DE.results.TMM.DE20$mean.fdr <- subset(DE.results.TMM.DE20$mean.fdr, 
                                       select = -noif.DESeq)
DE.results.TMM.DE20$mean.discoveries <- subset(DE.results.TMM.DE20$mean.discoveries, 
                                               select = -noif.DESeq)
DE.results.DESeq.DE2$mean.fdr <- subset(DE.results.DESeq.DE2$mean.fdr, 
                                        select = -noif.DESeq)
DE.results.DESeq.DE2$mean.discoveries <- subset(DE.results.DESeq.DE2$mean.discoveries, 
                                                select = -noif.DESeq)
DE.results.DESeq.DE5$mean.fdr <- subset(DE.results.DESeq.DE5$mean.fdr, 
                                        select = -noif.DESeq)
DE.results.DESeq.DE5$mean.discoveries <- subset(DE.results.DESeq.DE5$mean.discoveries, 
                                                select = -noif.DESeq)
DE.results.DESeq.DE10$mean.fdr <- subset(DE.results.DESeq.DE10$mean.fdr, 
                                         select = -noif.DESeq)
DE.results.DESeq.DE10$mean.discoveries <- subset(DE.results.DESeq.DE10$mean.discoveries, 
                                                 select = -noif.DESeq)
DE.results.DESeq.DE20$mean.fdr <- subset(DE.results.DESeq.DE20$mean.fdr, 
                                         select = -noif.DESeq)
DE.results.DESeq.DE20$mean.discoveries <- subset(DE.results.DESeq.DE20$mean.discoveries, 
                                                 select = -noif.DESeq)
dim(DE.results.TMM.DE20$mean.fdr)
# 20000 13
dim(DEDD.results.TMM.DE20$mean.fdr)
# 20000 2

dim(DD.results.TMM.DEDD20$mean.fdr)
# 20000 6
DE.results.TMM.DEDD2$mean.fdr <- subset(DE.results.TMM.DEDD2$mean.fdr, 
                                        select = -c(trend.DSS, mix.ShrinkBayes, 
                                                    np.ShrinkBayes, noif.DESeq))
DE.results.TMM.DEDD2$mean.discoveries <- subset(DE.results.TMM.DEDD2$mean.discoveries, 
                                                select = -c(trend.DSS, mix.ShrinkBayes, 
                                                            np.ShrinkBayes, noif.DESeq))
DE.results.TMM.DEDD5$mean.fdr <- subset(DE.results.TMM.DEDD5$mean.fdr, 
                                        select = -noif.DESeq)
DE.results.TMM.DEDD5$mean.discoveries <- subset(DE.results.TMM.DEDD5$mean.discoveries, 
                                                select = -noif.DESeq)
DE.results.TMM.DEDD10$mean.fdr <- subset(DE.results.TMM.DEDD10$mean.fdr, 
                                         select = -noif.DESeq)
DE.results.TMM.DEDD10$mean.discoveries <- subset(DE.results.TMM.DEDD10$mean.discoveries, 
                                                 select = -noif.DESeq)
DE.results.TMM.DEDD20$mean.fdr <- subset(DE.results.TMM.DEDD20$mean.fdr, 
                                         select = -noif.DESeq)
DE.results.TMM.DEDD20$mean.discoveries <- subset(DE.results.TMM.DEDD20$mean.discoveries, 
                                                 select = -noif.DESeq)
DE.results.DESeq.DEDD2$mean.fdr <- subset(DE.results.DESeq.DEDD2$mean.fdr, 
                                          select = -noif.DESeq)
DE.results.DESeq.DEDD2$mean.discoveries <- subset(DE.results.DESeq.DEDD2$mean.discoveries, 
                                                  select = -noif.DESeq)
DE.results.DESeq.DEDD5$mean.fdr <- subset(DE.results.DESeq.DEDD5$mean.fdr, 
                                          select = -noif.DESeq)
DE.results.DESeq.DEDD5$mean.discoveries <- subset(DE.results.DESeq.DEDD5$mean.discoveries, 
                                                  select = -noif.DESeq)
DE.results.DESeq.DEDD10$mean.fdr <- subset(DE.results.DESeq.DEDD10$mean.fdr, 
                                           select = -noif.DESeq)
DE.results.DESeq.DEDD10$mean.discoveries <- subset(DE.results.DESeq.DEDD10$mean.discoveries, 
                                                   select = -noif.DESeq)
DE.results.DESeq.DEDD20$mean.fdr <- subset(DE.results.DESeq.DEDD20$mean.fdr, 
                                           select = -noif.DESeq)
DE.results.DESeq.DEDD20$mean.discoveries <- subset(DE.results.DESeq.DEDD20$mean.discoveries, 
                                                   select = -noif.DESeq)
dim(DE.results.TMM.DEDD20$mean.fdr)
# 20000 13
dim(DEDD.results.TMM.DEDD20$mean.fdr)
# 20000 2


# DD2,5,10,20 DD: MDSeq(2), HMs(4)
# DD2,5,10,20 DEDD: HMMs(2)
# DE2,5,10,20 DE: edgeR(3), DESeq2(1), voom(1), DSS(1), baySeq(1), MDSeq(2), HMs(4)
# DE2,5,10,20 DEDD: HMMs(2)
# DEDD2,5,10,20 DD: MDSeq(2), HMs(4)
# DEDD2,5,10,20 DE: edgeR(3), DESeq2(1), voom(1), DSS(1), baySeq(1), MDSeq(2), HMs(4)
# DEDD2,5,10,20 DEDD: HMMs(2)

#################################
#### Differential dispersion ####
#################################

## DD2 ##
par(mfrow=c(2,6), mar=c(2,2,1,1), mgp=c(3,0.7,0))
plot(DD.results.TMM.DD2$mean.discoveries$disp.zi.MDSeq, 
     DD.results.TMM.DD2$mean.fdr$disp.zi.MDSeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='MDSeq +ZI, TMM', cex=1.2)
plot(DD.results.TMM.DD2$mean.discoveries$disp.nozi.MDSeq, 
     DD.results.TMM.DD2$mean.fdr$disp.nozi.MDSeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='MDSeq -ZI, TMM', cex=1.2)
plot(DD.results.TMM.DD2$mean.discoveries$disp.expHM, 
     DD.results.TMM.DD2$mean.fdr$disp.expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='expHM, TMM', cex=1.2)
plot(DD.results.TMM.DD2$mean.discoveries$ldisp.expHM, 
     DD.results.TMM.DD2$mean.fdr$ldisp.expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='expHM log, TMM', cex=1.2)
plot(DD.results.TMM.DD2$mean.discoveries$disp.lnHM, 
     DD.results.TMM.DD2$mean.fdr$disp.lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='lnHM, TMM', cex=1.2)
plot(DD.results.TMM.DD2$mean.discoveries$ldisp.lnHM, 
     DD.results.TMM.DD2$mean.fdr$ldisp.lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='lnHM log, TMM', cex=1.2)
plot(DD.results.DESeq.DD2$mean.discoveries$disp.zi.MDSeq, 
     DD.results.DESeq.DD2$mean.fdr$disp.zi.MDSeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='MDSeq +ZI, DESeq norm', cex=1.2)
plot(DD.results.DESeq.DD2$mean.discoveries$disp.nozi.MDSeq, 
     DD.results.DESeq.DD2$mean.fdr$disp.nozi.MDSeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='MDSeq -ZI, DESeq norm', cex=1.2)
plot(DD.results.DESeq.DD2$mean.discoveries$disp.expHM, 
     DD.results.DESeq.DD2$mean.fdr$disp.expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='expHM, DESeq norm', cex=1.2)
plot(DD.results.DESeq.DD2$mean.discoveries$ldisp.expHM, 
     DD.results.DESeq.DD2$mean.fdr$ldisp.expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='expHM log, DESeq norm', cex=1.2)
plot(DD.results.DESeq.DD2$mean.discoveries$disp.lnHM, 
     DD.results.DESeq.DD2$mean.fdr$disp.lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='lnHM, DESeq norm', cex=1.2)
plot(DD.results.DESeq.DD2$mean.discoveries$ldisp.lnHM, 
     DD.results.DESeq.DD2$mean.fdr$ldisp.lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='lnHM log, DESeq norm', cex=1.2)

## DD5 ##
par(mfrow=c(2,6), mar=c(2,2,1,1), mgp=c(3,0.7,0))
plot(DD.results.TMM.DD5$mean.discoveries$disp.zi.MDSeq, 
     DD.results.TMM.DD5$mean.fdr$disp.zi.MDSeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='MDSeq +ZI, TMM', cex=1.2)
plot(DD.results.TMM.DD5$mean.discoveries$disp.nozi.MDSeq, 
     DD.results.TMM.DD5$mean.fdr$disp.nozi.MDSeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='MDSeq -ZI, TMM', cex=1.2)
plot(DD.results.TMM.DD5$mean.discoveries$disp.expHM, 
     DD.results.TMM.DD5$mean.fdr$disp.expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='expHM, TMM', cex=1.2)
plot(DD.results.TMM.DD5$mean.discoveries$ldisp.expHM, 
     DD.results.TMM.DD5$mean.fdr$ldisp.expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='expHM log, TMM', cex=1.2)
plot(DD.results.TMM.DD5$mean.discoveries$disp.lnHM, 
     DD.results.TMM.DD5$mean.fdr$disp.lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='lnHM, TMM', cex=1.2)
plot(DD.results.TMM.DD5$mean.discoveries$ldisp.lnHM, 
     DD.results.TMM.DD5$mean.fdr$ldisp.lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='lnHM log, TMM', cex=1.2)
plot(DD.results.DESeq.DD5$mean.discoveries$disp.zi.MDSeq, 
     DD.results.DESeq.DD5$mean.fdr$disp.zi.MDSeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='MDSeq +ZI, DESeq norm', cex=1.2)
plot(DD.results.DESeq.DD5$mean.discoveries$disp.nozi.MDSeq, 
     DD.results.DESeq.DD5$mean.fdr$disp.nozi.MDSeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='MDSeq -ZI, DESeq norm', cex=1.2)
plot(DD.results.DESeq.DD5$mean.discoveries$disp.expHM, 
     DD.results.DESeq.DD5$mean.fdr$disp.expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='expHM, DESeq norm', cex=1.2)
plot(DD.results.DESeq.DD5$mean.discoveries$ldisp.expHM, 
     DD.results.DESeq.DD5$mean.fdr$ldisp.expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='expHM log, DESeq norm', cex=1.2)
plot(DD.results.DESeq.DD5$mean.discoveries$disp.lnHM, 
     DD.results.DESeq.DD5$mean.fdr$disp.lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='lnHM, DESeq norm', cex=1.2)
plot(DD.results.DESeq.DD5$mean.discoveries$ldisp.lnHM, 
     DD.results.DESeq.DD5$mean.fdr$ldisp.lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='lnHM log, DESeq norm', cex=1.2)

## DD10 ##
par(mfrow=c(2,6), mar=c(2,2,1,1), mgp=c(3,0.7,0))
plot(DD.results.TMM.DD10$mean.discoveries$disp.zi.MDSeq, 
     DD.results.TMM.DD10$mean.fdr$disp.zi.MDSeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='MDSeq +ZI, TMM', cex=1.2)
plot(DD.results.TMM.DD10$mean.discoveries$disp.nozi.MDSeq, 
     DD.results.TMM.DD10$mean.fdr$disp.nozi.MDSeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='MDSeq -ZI, TMM', cex=1.2)
plot(DD.results.TMM.DD10$mean.discoveries$disp.expHM, 
     DD.results.TMM.DD10$mean.fdr$disp.expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='expHM, TMM', cex=1.2)
plot(DD.results.TMM.DD10$mean.discoveries$ldisp.expHM, 
     DD.results.TMM.DD10$mean.fdr$ldisp.expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='expHM log, TMM', cex=1.2)
plot(DD.results.TMM.DD10$mean.discoveries$disp.lnHM, 
     DD.results.TMM.DD10$mean.fdr$disp.lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='lnHM, TMM', cex=1.2)
plot(DD.results.TMM.DD10$mean.discoveries$ldisp.lnHM, 
     DD.results.TMM.DD10$mean.fdr$ldisp.lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='lnHM log, TMM', cex=1.2)
plot(DD.results.DESeq.DD10$mean.discoveries$disp.zi.MDSeq, 
     DD.results.DESeq.DD10$mean.fdr$disp.zi.MDSeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='MDSeq +ZI, DESeq norm', cex=1.2)
plot(DD.results.DESeq.DD10$mean.discoveries$disp.nozi.MDSeq, 
     DD.results.DESeq.DD10$mean.fdr$disp.nozi.MDSeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='MDSeq -ZI, DESeq norm', cex=1.2)
plot(DD.results.DESeq.DD10$mean.discoveries$disp.expHM, 
     DD.results.DESeq.DD10$mean.fdr$disp.expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='expHM, DESeq norm', cex=1.2)
plot(DD.results.DESeq.DD10$mean.discoveries$ldisp.expHM, 
     DD.results.DESeq.DD10$mean.fdr$ldisp.expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='expHM log, DESeq norm', cex=1.2)
plot(DD.results.DESeq.DD10$mean.discoveries$disp.lnHM, 
     DD.results.DESeq.DD10$mean.fdr$disp.lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='lnHM, DESeq norm', cex=1.2)
plot(DD.results.DESeq.DD10$mean.discoveries$ldisp.lnHM, 
     DD.results.DESeq.DD10$mean.fdr$ldisp.lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='lnHM log, DESeq norm', cex=1.2)

## DD20 ##
par(mfrow=c(2,6), mar=c(2,2,1,1), mgp=c(3,0.7,0))
plot(DD.results.TMM.DD20$mean.discoveries$disp.zi.MDSeq, 
     DD.results.TMM.DD20$mean.fdr$disp.zi.MDSeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='MDSeq +ZI, TMM', cex=1.2)
plot(DD.results.TMM.DD20$mean.discoveries$disp.nozi.MDSeq, 
     DD.results.TMM.DD20$mean.fdr$disp.nozi.MDSeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='MDSeq -ZI, TMM', cex=1.2)
plot(DD.results.TMM.DD20$mean.discoveries$disp.expHM, 
     DD.results.TMM.DD20$mean.fdr$disp.expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='expHM, TMM', cex=1.2)
plot(DD.results.TMM.DD20$mean.discoveries$ldisp.expHM, 
     DD.results.TMM.DD20$mean.fdr$ldisp.expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='expHM log, TMM', cex=1.2)
plot(DD.results.TMM.DD20$mean.discoveries$disp.lnHM, 
     DD.results.TMM.DD20$mean.fdr$disp.lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='lnHM, TMM', cex=1.2)
plot(DD.results.TMM.DD20$mean.discoveries$ldisp.lnHM, 
     DD.results.TMM.DD20$mean.fdr$ldisp.lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='lnHM log, TMM', cex=1.2)
plot(DD.results.DESeq.DD20$mean.discoveries$disp.zi.MDSeq, 
     DD.results.DESeq.DD20$mean.fdr$disp.zi.MDSeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='MDSeq +ZI, DESeq norm', cex=1.2)
plot(DD.results.DESeq.DD20$mean.discoveries$disp.nozi.MDSeq, 
     DD.results.DESeq.DD20$mean.fdr$disp.nozi.MDSeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='MDSeq -ZI, DESeq norm', cex=1.2)
plot(DD.results.DESeq.DD20$mean.discoveries$disp.expHM, 
     DD.results.DESeq.DD20$mean.fdr$disp.expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='expHM, DESeq norm', cex=1.2)
plot(DD.results.DESeq.DD20$mean.discoveries$ldisp.expHM, 
     DD.results.DESeq.DD20$mean.fdr$ldisp.expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='expHM log, DESeq norm', cex=1.2)
plot(DD.results.DESeq.DD20$mean.discoveries$disp.lnHM, 
     DD.results.DESeq.DD20$mean.fdr$disp.lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='lnHM, DESeq norm', cex=1.2)
plot(DD.results.DESeq.DD20$mean.discoveries$ldisp.lnHM, 
     DD.results.DESeq.DD20$mean.fdr$ldisp.lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='lnHM log, DESeq norm', cex=1.2)

## DEDD2 ##
par(mfrow=c(2,6), mar=c(2,2,1,1), mgp=c(3,0.7,0))
plot(DD.results.TMM.DEDD2$mean.discoveries$disp.zi.MDSeq, 
     DD.results.TMM.DEDD2$mean.fdr$disp.zi.MDSeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='MDSeq +ZI, TMM', cex=1.2)
plot(DD.results.TMM.DEDD2$mean.discoveries$disp.nozi.MDSeq, 
     DD.results.TMM.DEDD2$mean.fdr$disp.nozi.MDSeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='MDSeq -ZI, TMM', cex=1.2)
plot(DD.results.TMM.DEDD2$mean.discoveries$disp.expHM, 
     DD.results.TMM.DEDD2$mean.fdr$disp.expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='expHM, TMM', cex=1.2)
plot(DD.results.TMM.DEDD2$mean.discoveries$ldisp.expHM, 
     DD.results.TMM.DEDD2$mean.fdr$ldisp.expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='expHM log, TMM', cex=1.2)
plot(DD.results.TMM.DEDD2$mean.discoveries$disp.lnHM, 
     DD.results.TMM.DEDD2$mean.fdr$disp.lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='lnHM, TMM', cex=1.2)
plot(DD.results.TMM.DEDD2$mean.discoveries$ldisp.lnHM, 
     DD.results.TMM.DEDD2$mean.fdr$ldisp.lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='lnHM log, TMM', cex=1.2)
plot(DD.results.DESeq.DEDD2$mean.discoveries$disp.zi.MDSeq, 
     DD.results.DESeq.DEDD2$mean.fdr$disp.zi.MDSeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='MDSeq +ZI, DESeq norm', cex=1.2)
plot(DD.results.DESeq.DEDD2$mean.discoveries$disp.nozi.MDSeq, 
     DD.results.DESeq.DEDD2$mean.fdr$disp.nozi.MDSeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='MDSeq -ZI, DESeq norm', cex=1.2)
plot(DD.results.DESeq.DEDD2$mean.discoveries$disp.expHM, 
     DD.results.DESeq.DEDD2$mean.fdr$disp.expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='expHM, DESeq norm', cex=1.2)
plot(DD.results.DESeq.DEDD2$mean.discoveries$ldisp.expHM, 
     DD.results.DESeq.DEDD2$mean.fdr$ldisp.expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='expHM log, DESeq norm', cex=1.2)
plot(DD.results.DESeq.DEDD2$mean.discoveries$disp.lnHM, 
     DD.results.DESeq.DEDD2$mean.fdr$disp.lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='lnHM, DESeq norm', cex=1.2)
plot(DD.results.DESeq.DEDD2$mean.discoveries$ldisp.lnHM, 
     DD.results.DESeq.DEDD2$mean.fdr$ldisp.lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='lnHM log, DESeq norm', cex=1.2)

## DEDD5 ##
par(mfrow=c(2,6), mar=c(2,2,1,1), mgp=c(3,0.7,0))
plot(DD.results.TMM.DEDD5$mean.discoveries$disp.zi.MDSeq, 
     DD.results.TMM.DEDD5$mean.fdr$disp.zi.MDSeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='MDSeq +ZI, TMM', cex=1.2)
plot(DD.results.TMM.DEDD5$mean.discoveries$disp.nozi.MDSeq, 
     DD.results.TMM.DEDD5$mean.fdr$disp.nozi.MDSeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='MDSeq -ZI, TMM', cex=1.2)
plot(DD.results.TMM.DEDD5$mean.discoveries$disp.expHM, 
     DD.results.TMM.DEDD5$mean.fdr$disp.expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='expHM, TMM', cex=1.2)
plot(DD.results.TMM.DEDD5$mean.discoveries$ldisp.expHM, 
     DD.results.TMM.DEDD5$mean.fdr$ldisp.expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='expHM log, TMM', cex=1.2)
plot(DD.results.TMM.DEDD5$mean.discoveries$disp.lnHM, 
     DD.results.TMM.DEDD5$mean.fdr$disp.lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='lnHM, TMM', cex=1.2)
plot(DD.results.TMM.DEDD5$mean.discoveries$ldisp.lnHM, 
     DD.results.TMM.DEDD5$mean.fdr$ldisp.lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='lnHM log, TMM', cex=1.2)
plot(DD.results.DESeq.DEDD5$mean.discoveries$disp.zi.MDSeq, 
     DD.results.DESeq.DEDD5$mean.fdr$disp.zi.MDSeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='MDSeq +ZI, DESeq norm', cex=1.2)
plot(DD.results.DESeq.DEDD5$mean.discoveries$disp.nozi.MDSeq, 
     DD.results.DESeq.DEDD5$mean.fdr$disp.nozi.MDSeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='MDSeq -ZI, DESeq norm', cex=1.2)
plot(DD.results.DESeq.DEDD5$mean.discoveries$disp.expHM, 
     DD.results.DESeq.DEDD5$mean.fdr$disp.expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='expHM, DESeq norm', cex=1.2)
plot(DD.results.DESeq.DEDD5$mean.discoveries$ldisp.expHM, 
     DD.results.DESeq.DEDD5$mean.fdr$ldisp.expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='expHM log, DESeq norm', cex=1.2)
plot(DD.results.DESeq.DEDD5$mean.discoveries$disp.lnHM, 
     DD.results.DESeq.DEDD5$mean.fdr$disp.lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='lnHM, DESeq norm', cex=1.2)
plot(DD.results.DESeq.DEDD5$mean.discoveries$ldisp.lnHM, 
     DD.results.DESeq.DEDD5$mean.fdr$ldisp.lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='lnHM log, DESeq norm', cex=1.2)

## DEDD10 ##
par(mfrow=c(2,6), mar=c(2,2,1,1), mgp=c(3,0.7,0))
plot(DD.results.TMM.DEDD10$mean.discoveries$disp.zi.MDSeq, 
     DD.results.TMM.DEDD10$mean.fdr$disp.zi.MDSeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='MDSeq +ZI, TMM', cex=1.2)
plot(DD.results.TMM.DEDD10$mean.discoveries$disp.nozi.MDSeq, 
     DD.results.TMM.DEDD10$mean.fdr$disp.nozi.MDSeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='MDSeq -ZI, TMM', cex=1.2)
plot(DD.results.TMM.DEDD10$mean.discoveries$disp.expHM, 
     DD.results.TMM.DEDD10$mean.fdr$disp.expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='expHM, TMM', cex=1.2)
plot(DD.results.TMM.DEDD10$mean.discoveries$ldisp.expHM, 
     DD.results.TMM.DEDD10$mean.fdr$ldisp.expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='expHM log, TMM', cex=1.2)
plot(DD.results.TMM.DEDD10$mean.discoveries$disp.lnHM, 
     DD.results.TMM.DEDD10$mean.fdr$disp.lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='lnHM, TMM', cex=1.2)
plot(DD.results.TMM.DEDD10$mean.discoveries$ldisp.lnHM, 
     DD.results.TMM.DEDD10$mean.fdr$ldisp.lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='lnHM log, TMM', cex=1.2)
plot(DD.results.DESeq.DEDD10$mean.discoveries$disp.zi.MDSeq, 
     DD.results.DESeq.DEDD10$mean.fdr$disp.zi.MDSeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='MDSeq +ZI, DESeq norm', cex=1.2)
plot(DD.results.DESeq.DEDD10$mean.discoveries$disp.nozi.MDSeq, 
     DD.results.DESeq.DEDD10$mean.fdr$disp.nozi.MDSeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='MDSeq -ZI, DESeq norm', cex=1.2)
plot(DD.results.DESeq.DEDD10$mean.discoveries$disp.expHM, 
     DD.results.DESeq.DEDD10$mean.fdr$disp.expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='expHM, DESeq norm', cex=1.2)
plot(DD.results.DESeq.DEDD10$mean.discoveries$ldisp.expHM, 
     DD.results.DESeq.DEDD10$mean.fdr$ldisp.expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='expHM log, DESeq norm', cex=1.2)
plot(DD.results.DESeq.DEDD10$mean.discoveries$disp.lnHM, 
     DD.results.DESeq.DEDD10$mean.fdr$disp.lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='lnHM, DESeq norm', cex=1.2)
plot(DD.results.DESeq.DEDD10$mean.discoveries$ldisp.lnHM, 
     DD.results.DESeq.DEDD10$mean.fdr$ldisp.lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='lnHM log, DESeq norm', cex=1.2)

## DEDD20 ##
par(mfrow=c(2,6), mar=c(2,2,1,1), mgp=c(3,0.7,0))
plot(DD.results.TMM.DEDD20$mean.discoveries$disp.zi.MDSeq, 
     DD.results.TMM.DEDD20$mean.fdr$disp.zi.MDSeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='MDSeq +ZI, TMM', cex=1.2)
plot(DD.results.TMM.DEDD20$mean.discoveries$disp.nozi.MDSeq, 
     DD.results.TMM.DEDD20$mean.fdr$disp.nozi.MDSeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='MDSeq -ZI, TMM', cex=1.2)
plot(DD.results.TMM.DEDD20$mean.discoveries$disp.expHM, 
     DD.results.TMM.DEDD20$mean.fdr$disp.expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='expHM, TMM', cex=1.2)
plot(DD.results.TMM.DEDD20$mean.discoveries$ldisp.expHM, 
     DD.results.TMM.DEDD20$mean.fdr$ldisp.expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='expHM log, TMM', cex=1.2)
plot(DD.results.TMM.DEDD20$mean.discoveries$disp.lnHM, 
     DD.results.TMM.DEDD20$mean.fdr$disp.lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='lnHM, TMM', cex=1.2)
plot(DD.results.TMM.DEDD20$mean.discoveries$ldisp.lnHM, 
     DD.results.TMM.DEDD20$mean.fdr$ldisp.lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='lnHM log, TMM', cex=1.2)
plot(DD.results.DESeq.DEDD20$mean.discoveries$disp.zi.MDSeq, 
     DD.results.DESeq.DEDD20$mean.fdr$disp.zi.MDSeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='MDSeq +ZI, DESeq norm', cex=1.2)
plot(DD.results.DESeq.DEDD20$mean.discoveries$disp.nozi.MDSeq, 
     DD.results.DESeq.DEDD20$mean.fdr$disp.nozi.MDSeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='MDSeq -ZI, DESeq norm', cex=1.2)
plot(DD.results.DESeq.DEDD20$mean.discoveries$disp.expHM, 
     DD.results.DESeq.DEDD20$mean.fdr$disp.expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='expHM, DESeq norm', cex=1.2)
plot(DD.results.DESeq.DEDD20$mean.discoveries$ldisp.expHM, 
     DD.results.DESeq.DEDD20$mean.fdr$ldisp.expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='expHM log, DESeq norm', cex=1.2)
plot(DD.results.DESeq.DEDD20$mean.discoveries$disp.lnHM, 
     DD.results.DESeq.DEDD20$mean.fdr$disp.lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='lnHM, DESeq norm', cex=1.2)
plot(DD.results.DESeq.DEDD20$mean.discoveries$ldisp.lnHM, 
     DD.results.DESeq.DEDD20$mean.fdr$ldisp.lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='lnHM log, DESeq norm', cex=1.2)


#################################
#### Differential expression ####
#################################

## DE2 ##
par(mfrow=c(4,7), mar=c(2,2,1,1), mgp=c(3,0.7,0))
plot(DE.results.TMM.DE2$mean.discoveries$ql.edgeR, 
     DE.results.TMM.DE2$mean.fdr$ql.edgeR, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='edgeR QL, TMM')
plot(DE.results.TMM.DE2$mean.discoveries$lr.edgeR, 
     DE.results.TMM.DE2$mean.fdr$lr.edgeR, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='edgeR LR, TMM')
plot(DE.results.TMM.DE2$mean.discoveries$et.edgeR, 
     DE.results.TMM.DE2$mean.fdr$et.edgeR, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='edgeR ET, TMM')
plot(DE.results.TMM.DE2$mean.discoveries$if.DESeq, 
     DE.results.TMM.DE2$mean.fdr$if.DESeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='DESeq2, TMM')
plot(DE.results.TMM.DE2$mean.discoveries$voom, 
     DE.results.TMM.DE2$mean.fdr$voom, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='limma, TMM')
plot(DE.results.TMM.DE2$mean.discoveries$notrend.DSS, 
     DE.results.TMM.DE2$mean.fdr$notrend.DSS, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='DSS, TMM')
plot(DE.results.TMM.DE2$mean.discoveries$baySeq, 
     DE.results.TMM.DE2$mean.fdr$baySeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='baySeq, TMM')
plot(DE.results.TMM.DE2$mean.discoveries$mean.zi.MDSeq, 
     DE.results.TMM.DE2$mean.fdr$mean.zi.MDSeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='MDSeq +ZI, TMM')
plot(DE.results.TMM.DE2$mean.discoveries$mean.nozi.MDSeq, 
     DE.results.TMM.DE2$mean.fdr$mean.nozi.MDSeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='MDSeq -ZI')
plot(DE.results.TMM.DE2$mean.discoveries$mean.expHM, 
     DE.results.TMM.DE2$mean.fdr$mean.expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='expHM, TMM')
plot(DE.results.TMM.DE2$mean.discoveries$lmean.expHM, 
     DE.results.TMM.DE2$mean.fdr$lmean.expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='expHM log, TMM')
plot(DE.results.TMM.DE2$mean.discoveries$mean.lnHM, 
     DE.results.TMM.DE2$mean.fdr$mean.lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='lnHM, TMM')
plot(DE.results.TMM.DE2$mean.discoveries$lmean.lnHM, 
     DE.results.TMM.DE2$mean.fdr$lmean.lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='lnHM log, TMM')
plot.new()
plot(DE.results.DESeq.DE2$mean.discoveries$ql.edgeR, 
     DE.results.DESeq.DE2$mean.fdr$ql.edgeR, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='edgeR QL, DESeq norm')
plot(DE.results.DESeq.DE2$mean.discoveries$lr.edgeR, 
     DE.results.DESeq.DE2$mean.fdr$lr.edgeR, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='edgeR LR, DESeq norm')
plot(DE.results.DESeq.DE2$mean.discoveries$et.edgeR, 
     DE.results.DESeq.DE2$mean.fdr$et.edgeR, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='edgeR ET, DESeq norm')
plot(DE.results.DESeq.DE2$mean.discoveries$if.DESeq, 
     DE.results.DESeq.DE2$mean.fdr$if.DESeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='DESeq2, DESeq norm')
plot(DE.results.DESeq.DE2$mean.discoveries$voom, 
     DE.results.DESeq.DE2$mean.fdr$voom, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='limma, DESeq norm')
plot(DE.results.DESeq.DE2$mean.discoveries$notrend.DSS, 
     DE.results.DESeq.DE2$mean.fdr$notrend.DSS, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='DSS, DESeq norm')
plot(DE.results.DESeq.DE2$mean.discoveries$baySeq, 
     DE.results.DESeq.DE2$mean.fdr$baySeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='baySeq, DESeq norm')
plot(DE.results.DESeq.DE2$mean.discoveries$mean.zi.MDSeq, 
     DE.results.DESeq.DE2$mean.fdr$mean.zi.MDSeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='MDSeq +ZI, DESeq norm')
plot(DE.results.DESeq.DE2$mean.discoveries$mean.nozi.MDSeq, 
     DE.results.DESeq.DE2$mean.fdr$mean.nozi.MDSeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='MDSeq -ZI')
plot(DE.results.DESeq.DE2$mean.discoveries$mean.expHM, 
     DE.results.DESeq.DE2$mean.fdr$mean.expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='expHM, DESeq norm')
plot(DE.results.DESeq.DE2$mean.discoveries$lmean.expHM, 
     DE.results.DESeq.DE2$mean.fdr$lmean.expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='expHM log, DESeq norm')
plot(DE.results.DESeq.DE2$mean.discoveries$mean.lnHM, 
     DE.results.DESeq.DE2$mean.fdr$mean.lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='lnHM, DESeq norm')
plot(DE.results.DESeq.DE2$mean.discoveries$lmean.lnHM, 
     DE.results.DESeq.DE2$mean.fdr$lmean.lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='lnHM log, DESeq norm')
plot.new()

## DE5 ##
par(mfrow=c(4,7), mar=c(2,2,1,1), mgp=c(3,0.7,0))
plot(DE.results.TMM.DE5$mean.discoveries$ql.edgeR, 
     DE.results.TMM.DE5$mean.fdr$ql.edgeR, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='edgeR QL, TMM')
plot(DE.results.TMM.DE5$mean.discoveries$lr.edgeR, 
     DE.results.TMM.DE5$mean.fdr$lr.edgeR, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='edgeR LR, TMM')
plot(DE.results.TMM.DE5$mean.discoveries$et.edgeR, 
     DE.results.TMM.DE5$mean.fdr$et.edgeR, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='edgeR ET, TMM')
plot(DE.results.TMM.DE5$mean.discoveries$if.DESeq, 
     DE.results.TMM.DE5$mean.fdr$if.DESeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='DESeq2, TMM')
plot(DE.results.TMM.DE5$mean.discoveries$voom, 
     DE.results.TMM.DE5$mean.fdr$voom, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='limma, TMM')
plot(DE.results.TMM.DE5$mean.discoveries$notrend.DSS, 
     DE.results.TMM.DE5$mean.fdr$notrend.DSS, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='DSS, TMM')
plot(DE.results.TMM.DE5$mean.discoveries$baySeq, 
     DE.results.TMM.DE5$mean.fdr$baySeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='baySeq, TMM')
plot(DE.results.TMM.DE5$mean.discoveries$mean.zi.MDSeq, 
     DE.results.TMM.DE5$mean.fdr$mean.zi.MDSeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='MDSeq +ZI, TMM')
plot(DE.results.TMM.DE5$mean.discoveries$mean.nozi.MDSeq, 
     DE.results.TMM.DE5$mean.fdr$mean.nozi.MDSeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='MDSeq -ZI')
plot(DE.results.TMM.DE5$mean.discoveries$mean.expHM, 
     DE.results.TMM.DE5$mean.fdr$mean.expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='expHM, TMM')
plot(DE.results.TMM.DE5$mean.discoveries$lmean.expHM, 
     DE.results.TMM.DE5$mean.fdr$lmean.expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='expHM log, TMM')
plot(DE.results.TMM.DE5$mean.discoveries$mean.lnHM, 
     DE.results.TMM.DE5$mean.fdr$mean.lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='lnHM, TMM')
plot(DE.results.TMM.DE5$mean.discoveries$lmean.lnHM, 
     DE.results.TMM.DE5$mean.fdr$lmean.lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='lnHM log, TMM')
plot.new()
plot(DE.results.DESeq.DE5$mean.discoveries$ql.edgeR, 
     DE.results.DESeq.DE5$mean.fdr$ql.edgeR, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='edgeR QL, DESeq norm')
plot(DE.results.DESeq.DE5$mean.discoveries$lr.edgeR, 
     DE.results.DESeq.DE5$mean.fdr$lr.edgeR, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='edgeR LR, DESeq norm')
plot(DE.results.DESeq.DE5$mean.discoveries$et.edgeR, 
     DE.results.DESeq.DE5$mean.fdr$et.edgeR, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='edgeR ET, DESeq norm')
plot(DE.results.DESeq.DE5$mean.discoveries$if.DESeq, 
     DE.results.DESeq.DE5$mean.fdr$if.DESeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='DESeq2, DESeq norm')
plot(DE.results.DESeq.DE5$mean.discoveries$voom, 
     DE.results.DESeq.DE5$mean.fdr$voom, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='limma, DESeq norm')
plot(DE.results.DESeq.DE5$mean.discoveries$notrend.DSS, 
     DE.results.DESeq.DE5$mean.fdr$notrend.DSS, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='DSS, DESeq norm')
plot(DE.results.DESeq.DE5$mean.discoveries$baySeq, 
     DE.results.DESeq.DE5$mean.fdr$baySeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='baySeq, DESeq norm')
plot(DE.results.DESeq.DE5$mean.discoveries$mean.zi.MDSeq, 
     DE.results.DESeq.DE5$mean.fdr$mean.zi.MDSeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='MDSeq +ZI, DESeq norm')
plot(DE.results.DESeq.DE5$mean.discoveries$mean.nozi.MDSeq, 
     DE.results.DESeq.DE5$mean.fdr$mean.nozi.MDSeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='MDSeq -ZI')
plot(DE.results.DESeq.DE5$mean.discoveries$mean.expHM, 
     DE.results.DESeq.DE5$mean.fdr$mean.expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='expHM, DESeq norm')
plot(DE.results.DESeq.DE5$mean.discoveries$lmean.expHM, 
     DE.results.DESeq.DE5$mean.fdr$lmean.expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='expHM log, DESeq norm')
plot(DE.results.DESeq.DE5$mean.discoveries$mean.lnHM, 
     DE.results.DESeq.DE5$mean.fdr$mean.lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='lnHM, DESeq norm')
plot(DE.results.DESeq.DE5$mean.discoveries$lmean.lnHM, 
     DE.results.DESeq.DE5$mean.fdr$lmean.lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='lnHM log, DESeq norm')
plot.new()

## DE10 ##
par(mfrow=c(4,7), mar=c(2,2,1,1), mgp=c(3,0.7,0))
plot(DE.results.TMM.DE10$mean.discoveries$ql.edgeR, 
     DE.results.TMM.DE10$mean.fdr$ql.edgeR, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='edgeR QL, TMM')
plot(DE.results.TMM.DE10$mean.discoveries$lr.edgeR, 
     DE.results.TMM.DE10$mean.fdr$lr.edgeR, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='edgeR LR, TMM')
plot(DE.results.TMM.DE10$mean.discoveries$et.edgeR, 
     DE.results.TMM.DE10$mean.fdr$et.edgeR, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='edgeR ET, TMM')
plot(DE.results.TMM.DE10$mean.discoveries$if.DESeq, 
     DE.results.TMM.DE10$mean.fdr$if.DESeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='DESeq2, TMM')
plot(DE.results.TMM.DE10$mean.discoveries$voom, 
     DE.results.TMM.DE10$mean.fdr$voom, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='limma, TMM')
plot(DE.results.TMM.DE10$mean.discoveries$notrend.DSS, 
     DE.results.TMM.DE10$mean.fdr$notrend.DSS, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='DSS, TMM')
plot(DE.results.TMM.DE10$mean.discoveries$baySeq, 
     DE.results.TMM.DE10$mean.fdr$baySeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='baySeq, TMM')
plot(DE.results.TMM.DE10$mean.discoveries$mean.zi.MDSeq, 
     DE.results.TMM.DE10$mean.fdr$mean.zi.MDSeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='MDSeq +ZI, TMM')
plot(DE.results.TMM.DE10$mean.discoveries$mean.nozi.MDSeq, 
     DE.results.TMM.DE10$mean.fdr$mean.nozi.MDSeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='MDSeq -ZI')
plot(DE.results.TMM.DE10$mean.discoveries$mean.expHM, 
     DE.results.TMM.DE10$mean.fdr$mean.expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='expHM, TMM')
plot(DE.results.TMM.DE10$mean.discoveries$lmean.expHM, 
     DE.results.TMM.DE10$mean.fdr$lmean.expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='expHM log, TMM')
plot(DE.results.TMM.DE10$mean.discoveries$mean.lnHM, 
     DE.results.TMM.DE10$mean.fdr$mean.lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='lnHM, TMM')
plot(DE.results.TMM.DE10$mean.discoveries$lmean.lnHM, 
     DE.results.TMM.DE10$mean.fdr$lmean.lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='lnHM log, TMM')
plot.new()
plot(DE.results.DESeq.DE10$mean.discoveries$ql.edgeR, 
     DE.results.DESeq.DE10$mean.fdr$ql.edgeR, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='edgeR QL, DESeq norm')
plot(DE.results.DESeq.DE10$mean.discoveries$lr.edgeR, 
     DE.results.DESeq.DE10$mean.fdr$lr.edgeR, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='edgeR LR, DESeq norm')
plot(DE.results.DESeq.DE10$mean.discoveries$et.edgeR, 
     DE.results.DESeq.DE10$mean.fdr$et.edgeR, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='edgeR ET, DESeq norm')
plot(DE.results.DESeq.DE10$mean.discoveries$if.DESeq, 
     DE.results.DESeq.DE10$mean.fdr$if.DESeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='DESeq2, DESeq norm')
plot(DE.results.DESeq.DE10$mean.discoveries$voom, 
     DE.results.DESeq.DE10$mean.fdr$voom, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='limma, DESeq norm')
plot(DE.results.DESeq.DE10$mean.discoveries$notrend.DSS, 
     DE.results.DESeq.DE10$mean.fdr$notrend.DSS, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='DSS, DESeq norm')
plot(DE.results.DESeq.DE10$mean.discoveries$baySeq, 
     DE.results.DESeq.DE10$mean.fdr$baySeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='baySeq, DESeq norm')
plot(DE.results.DESeq.DE10$mean.discoveries$mean.zi.MDSeq, 
     DE.results.DESeq.DE10$mean.fdr$mean.zi.MDSeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='MDSeq +ZI, DESeq norm')
plot(DE.results.DESeq.DE10$mean.discoveries$mean.nozi.MDSeq, 
     DE.results.DESeq.DE10$mean.fdr$mean.nozi.MDSeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='MDSeq -ZI')
plot(DE.results.DESeq.DE10$mean.discoveries$mean.expHM, 
     DE.results.DESeq.DE10$mean.fdr$mean.expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='expHM, DESeq norm')
plot(DE.results.DESeq.DE10$mean.discoveries$lmean.expHM, 
     DE.results.DESeq.DE10$mean.fdr$lmean.expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='expHM log, DESeq norm')
plot(DE.results.DESeq.DE10$mean.discoveries$mean.lnHM, 
     DE.results.DESeq.DE10$mean.fdr$mean.lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='lnHM, DESeq norm')
plot(DE.results.DESeq.DE10$mean.discoveries$lmean.lnHM, 
     DE.results.DESeq.DE10$mean.fdr$lmean.lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='lnHM log, DESeq norm')
plot.new()

## DE20 ##
par(mfrow=c(4,7), mar=c(2,2,1,1), mgp=c(3,0.7,0))
plot(DE.results.TMM.DE20$mean.discoveries$ql.edgeR, 
     DE.results.TMM.DE20$mean.fdr$ql.edgeR, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='edgeR QL, TMM')
plot(DE.results.TMM.DE20$mean.discoveries$lr.edgeR, 
     DE.results.TMM.DE20$mean.fdr$lr.edgeR, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='edgeR LR, TMM')
plot(DE.results.TMM.DE20$mean.discoveries$et.edgeR, 
     DE.results.TMM.DE20$mean.fdr$et.edgeR, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='edgeR ET, TMM')
plot(DE.results.TMM.DE20$mean.discoveries$if.DESeq, 
     DE.results.TMM.DE20$mean.fdr$if.DESeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='DESeq2, TMM')
plot(DE.results.TMM.DE20$mean.discoveries$voom, 
     DE.results.TMM.DE20$mean.fdr$voom, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='limma, TMM')
plot(DE.results.TMM.DE20$mean.discoveries$notrend.DSS, 
     DE.results.TMM.DE20$mean.fdr$notrend.DSS, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='DSS, TMM')
plot(DE.results.TMM.DE20$mean.discoveries$baySeq, 
     DE.results.TMM.DE20$mean.fdr$baySeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='baySeq, TMM')
plot(DE.results.TMM.DE20$mean.discoveries$mean.zi.MDSeq, 
     DE.results.TMM.DE20$mean.fdr$mean.zi.MDSeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='MDSeq +ZI, TMM')
plot(DE.results.TMM.DE20$mean.discoveries$mean.nozi.MDSeq, 
     DE.results.TMM.DE20$mean.fdr$mean.nozi.MDSeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='MDSeq -ZI')
plot(DE.results.TMM.DE20$mean.discoveries$mean.expHM, 
     DE.results.TMM.DE20$mean.fdr$mean.expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='expHM, TMM')
plot(DE.results.TMM.DE20$mean.discoveries$lmean.expHM, 
     DE.results.TMM.DE20$mean.fdr$lmean.expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='expHM log, TMM')
plot(DE.results.TMM.DE20$mean.discoveries$mean.lnHM, 
     DE.results.TMM.DE20$mean.fdr$mean.lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='lnHM, TMM')
plot(DE.results.TMM.DE20$mean.discoveries$lmean.lnHM, 
     DE.results.TMM.DE20$mean.fdr$lmean.lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='lnHM log, TMM')
plot.new()
plot(DE.results.DESeq.DE20$mean.discoveries$ql.edgeR, 
     DE.results.DESeq.DE20$mean.fdr$ql.edgeR, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='edgeR QL, DESeq norm')
plot(DE.results.DESeq.DE20$mean.discoveries$lr.edgeR, 
     DE.results.DESeq.DE20$mean.fdr$lr.edgeR, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='edgeR LR, DESeq norm')
plot(DE.results.DESeq.DE20$mean.discoveries$et.edgeR, 
     DE.results.DESeq.DE20$mean.fdr$et.edgeR, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='edgeR ET, DESeq norm')
plot(DE.results.DESeq.DE20$mean.discoveries$if.DESeq, 
     DE.results.DESeq.DE20$mean.fdr$if.DESeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='DESeq2, DESeq norm')
plot(DE.results.DESeq.DE20$mean.discoveries$voom, 
     DE.results.DESeq.DE20$mean.fdr$voom, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='limma, DESeq norm')
plot(DE.results.DESeq.DE20$mean.discoveries$notrend.DSS, 
     DE.results.DESeq.DE20$mean.fdr$notrend.DSS, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='DSS, DESeq norm')
plot(DE.results.DESeq.DE20$mean.discoveries$baySeq, 
     DE.results.DESeq.DE20$mean.fdr$baySeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='baySeq, DESeq norm')
plot(DE.results.DESeq.DE20$mean.discoveries$mean.zi.MDSeq, 
     DE.results.DESeq.DE20$mean.fdr$mean.zi.MDSeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='MDSeq +ZI, DESeq norm')
plot(DE.results.DESeq.DE20$mean.discoveries$mean.nozi.MDSeq, 
     DE.results.DESeq.DE20$mean.fdr$mean.nozi.MDSeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='MDSeq -ZI')
plot(DE.results.DESeq.DE20$mean.discoveries$mean.expHM, 
     DE.results.DESeq.DE20$mean.fdr$mean.expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='expHM, DESeq norm')
plot(DE.results.DESeq.DE20$mean.discoveries$lmean.expHM, 
     DE.results.DESeq.DE20$mean.fdr$lmean.expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='expHM log, DESeq norm')
plot(DE.results.DESeq.DE20$mean.discoveries$mean.lnHM, 
     DE.results.DESeq.DE20$mean.fdr$mean.lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='lnHM, DESeq norm')
plot(DE.results.DESeq.DE20$mean.discoveries$lmean.lnHM, 
     DE.results.DESeq.DE20$mean.fdr$lmean.lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='lnHM log, DESeq norm')
plot.new()

## DEDD2 ##
par(mfrow=c(4,7), mar=c(2,2,1,1), mgp=c(3,0.7,0))
plot(DE.results.TMM.DEDD2$mean.discoveries$ql.edgeR, 
     DE.results.TMM.DEDD2$mean.fdr$ql.edgeR, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='edgeR QL, TMM')
plot(DE.results.TMM.DEDD2$mean.discoveries$lr.edgeR, 
     DE.results.TMM.DEDD2$mean.fdr$lr.edgeR, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='edgeR LR, TMM')
plot(DE.results.TMM.DEDD2$mean.discoveries$et.edgeR, 
     DE.results.TMM.DEDD2$mean.fdr$et.edgeR, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='edgeR ET, TMM')
plot(DE.results.TMM.DEDD2$mean.discoveries$if.DESeq, 
     DE.results.TMM.DEDD2$mean.fdr$if.DESeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='DESeq2, TMM')
plot(DE.results.TMM.DEDD2$mean.discoveries$voom, 
     DE.results.TMM.DEDD2$mean.fdr$voom, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='limma, TMM')
plot(DE.results.TMM.DEDD2$mean.discoveries$notrend.DSS, 
     DE.results.TMM.DEDD2$mean.fdr$notrend.DSS, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='DSS, TMM')
plot(DE.results.TMM.DEDD2$mean.discoveries$baySeq, 
     DE.results.TMM.DEDD2$mean.fdr$baySeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='baySeq, TMM')
plot(DE.results.TMM.DEDD2$mean.discoveries$mean.zi.MDSeq, 
     DE.results.TMM.DEDD2$mean.fdr$mean.zi.MDSeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='MDSeq +ZI, TMM')
plot(DE.results.TMM.DEDD2$mean.discoveries$mean.nozi.MDSeq, 
     DE.results.TMM.DEDD2$mean.fdr$mean.nozi.MDSeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='MDSeq -ZI')
plot(DE.results.TMM.DEDD2$mean.discoveries$mean.expHM, 
     DE.results.TMM.DEDD2$mean.fdr$mean.expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='expHM, TMM')
plot(DE.results.TMM.DEDD2$mean.discoveries$lmean.expHM, 
     DE.results.TMM.DEDD2$mean.fdr$lmean.expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='expHM log, TMM')
plot(DE.results.TMM.DEDD2$mean.discoveries$mean.lnHM, 
     DE.results.TMM.DEDD2$mean.fdr$mean.lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='lnHM, TMM')
plot(DE.results.TMM.DEDD2$mean.discoveries$lmean.lnHM, 
     DE.results.TMM.DEDD2$mean.fdr$lmean.lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='lnHM log, TMM')
plot.new()
plot(DE.results.DESeq.DEDD2$mean.discoveries$ql.edgeR, 
     DE.results.DESeq.DEDD2$mean.fdr$ql.edgeR, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='edgeR QL, DESeq norm')
plot(DE.results.DESeq.DEDD2$mean.discoveries$lr.edgeR, 
     DE.results.DESeq.DEDD2$mean.fdr$lr.edgeR, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='edgeR LR, DESeq norm')
plot(DE.results.DESeq.DEDD2$mean.discoveries$et.edgeR, 
     DE.results.DESeq.DEDD2$mean.fdr$et.edgeR, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='edgeR ET, DESeq norm')
plot(DE.results.DESeq.DEDD2$mean.discoveries$if.DESeq, 
     DE.results.DESeq.DEDD2$mean.fdr$if.DESeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='DESeq2, DESeq norm')
plot(DE.results.DESeq.DEDD2$mean.discoveries$voom, 
     DE.results.DESeq.DEDD2$mean.fdr$voom, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='limma, DESeq norm')
plot(DE.results.DESeq.DEDD2$mean.discoveries$notrend.DSS, 
     DE.results.DESeq.DEDD2$mean.fdr$notrend.DSS, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='DSS, DESeq norm')
plot(DE.results.DESeq.DEDD2$mean.discoveries$baySeq, 
     DE.results.DESeq.DEDD2$mean.fdr$baySeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='baySeq, DESeq norm')
plot(DE.results.DESeq.DEDD2$mean.discoveries$mean.zi.MDSeq, 
     DE.results.DESeq.DEDD2$mean.fdr$mean.zi.MDSeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='MDSeq +ZI, DESeq norm')
plot(DE.results.DESeq.DEDD2$mean.discoveries$mean.nozi.MDSeq, 
     DE.results.DESeq.DEDD2$mean.fdr$mean.nozi.MDSeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='MDSeq -ZI')
plot(DE.results.DESeq.DEDD2$mean.discoveries$mean.expHM, 
     DE.results.DESeq.DEDD2$mean.fdr$mean.expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='expHM, DESeq norm')
plot(DE.results.DESeq.DEDD2$mean.discoveries$lmean.expHM, 
     DE.results.DESeq.DEDD2$mean.fdr$lmean.expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='expHM log, DESeq norm')
plot(DE.results.DESeq.DEDD2$mean.discoveries$mean.lnHM, 
     DE.results.DESeq.DEDD2$mean.fdr$mean.lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='lnHM, DESeq norm')
plot(DE.results.DESeq.DEDD2$mean.discoveries$lmean.lnHM, 
     DE.results.DESeq.DEDD2$mean.fdr$lmean.lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='lnHM log, DESeq norm')
plot.new()

## DEDD5 ##
par(mfrow=c(4,7), mar=c(2,2,1,1), mgp=c(3,0.7,0))
plot(DE.results.TMM.DEDD5$mean.discoveries$ql.edgeR, 
     DE.results.TMM.DEDD5$mean.fdr$ql.edgeR, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='edgeR QL, TMM')
plot(DE.results.TMM.DEDD5$mean.discoveries$lr.edgeR, 
     DE.results.TMM.DEDD5$mean.fdr$lr.edgeR, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='edgeR LR, TMM')
plot(DE.results.TMM.DEDD5$mean.discoveries$et.edgeR, 
     DE.results.TMM.DEDD5$mean.fdr$et.edgeR, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='edgeR ET, TMM')
plot(DE.results.TMM.DEDD5$mean.discoveries$if.DESeq, 
     DE.results.TMM.DEDD5$mean.fdr$if.DESeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='DESeq2, TMM')
plot(DE.results.TMM.DEDD5$mean.discoveries$voom, 
     DE.results.TMM.DEDD5$mean.fdr$voom, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='limma, TMM')
plot(DE.results.TMM.DEDD5$mean.discoveries$notrend.DSS, 
     DE.results.TMM.DEDD5$mean.fdr$notrend.DSS, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='DSS, TMM')
plot(DE.results.TMM.DEDD5$mean.discoveries$baySeq, 
     DE.results.TMM.DEDD5$mean.fdr$baySeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='baySeq, TMM')
plot(DE.results.TMM.DEDD5$mean.discoveries$mean.zi.MDSeq, 
     DE.results.TMM.DEDD5$mean.fdr$mean.zi.MDSeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='MDSeq +ZI, TMM')
plot(DE.results.TMM.DEDD5$mean.discoveries$mean.nozi.MDSeq, 
     DE.results.TMM.DEDD5$mean.fdr$mean.nozi.MDSeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='MDSeq -ZI')
plot(DE.results.TMM.DEDD5$mean.discoveries$mean.expHM, 
     DE.results.TMM.DEDD5$mean.fdr$mean.expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='expHM, TMM')
plot(DE.results.TMM.DEDD5$mean.discoveries$lmean.expHM, 
     DE.results.TMM.DEDD5$mean.fdr$lmean.expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='expHM log, TMM')
plot(DE.results.TMM.DEDD5$mean.discoveries$mean.lnHM, 
     DE.results.TMM.DEDD5$mean.fdr$mean.lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='lnHM, TMM')
plot(DE.results.TMM.DEDD5$mean.discoveries$lmean.lnHM, 
     DE.results.TMM.DEDD5$mean.fdr$lmean.lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='lnHM log, TMM')
plot.new()
plot(DE.results.DESeq.DEDD5$mean.discoveries$ql.edgeR, 
     DE.results.DESeq.DEDD5$mean.fdr$ql.edgeR, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='edgeR QL, DESeq norm')
plot(DE.results.DESeq.DEDD5$mean.discoveries$lr.edgeR, 
     DE.results.DESeq.DEDD5$mean.fdr$lr.edgeR, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='edgeR LR, DESeq norm')
plot(DE.results.DESeq.DEDD5$mean.discoveries$et.edgeR, 
     DE.results.DESeq.DEDD5$mean.fdr$et.edgeR, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='edgeR ET, DESeq norm')
plot(DE.results.DESeq.DEDD5$mean.discoveries$if.DESeq, 
     DE.results.DESeq.DEDD5$mean.fdr$if.DESeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='DESeq2, DESeq norm')
plot(DE.results.DESeq.DEDD5$mean.discoveries$voom, 
     DE.results.DESeq.DEDD5$mean.fdr$voom, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='limma, DESeq norm')
plot(DE.results.DESeq.DEDD5$mean.discoveries$notrend.DSS, 
     DE.results.DESeq.DEDD5$mean.fdr$notrend.DSS, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='DSS, DESeq norm')
plot(DE.results.DESeq.DEDD5$mean.discoveries$baySeq, 
     DE.results.DESeq.DEDD5$mean.fdr$baySeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='baySeq, DESeq norm')
plot(DE.results.DESeq.DEDD5$mean.discoveries$mean.zi.MDSeq, 
     DE.results.DESeq.DEDD5$mean.fdr$mean.zi.MDSeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='MDSeq +ZI, DESeq norm')
plot(DE.results.DESeq.DEDD5$mean.discoveries$mean.nozi.MDSeq, 
     DE.results.DESeq.DEDD5$mean.fdr$mean.nozi.MDSeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='MDSeq -ZI')
plot(DE.results.DESeq.DEDD5$mean.discoveries$mean.expHM, 
     DE.results.DESeq.DEDD5$mean.fdr$mean.expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='expHM, DESeq norm')
plot(DE.results.DESeq.DEDD5$mean.discoveries$lmean.expHM, 
     DE.results.DESeq.DEDD5$mean.fdr$lmean.expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='expHM log, DESeq norm')
plot(DE.results.DESeq.DEDD5$mean.discoveries$mean.lnHM, 
     DE.results.DESeq.DEDD5$mean.fdr$mean.lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='lnHM, DESeq norm')
plot(DE.results.DESeq.DEDD5$mean.discoveries$lmean.lnHM, 
     DE.results.DESeq.DEDD5$mean.fdr$lmean.lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='lnHM log, DESeq norm')
plot.new()

## DEDD10 ##
par(mfrow=c(4,7), mar=c(2,2,1,1), mgp=c(3,0.7,0))
plot(DE.results.TMM.DEDD10$mean.discoveries$ql.edgeR, 
     DE.results.TMM.DEDD10$mean.fdr$ql.edgeR, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='edgeR QL, TMM')
plot(DE.results.TMM.DEDD10$mean.discoveries$lr.edgeR, 
     DE.results.TMM.DEDD10$mean.fdr$lr.edgeR, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='edgeR LR, TMM')
plot(DE.results.TMM.DEDD10$mean.discoveries$et.edgeR, 
     DE.results.TMM.DEDD10$mean.fdr$et.edgeR, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='edgeR ET, TMM')
plot(DE.results.TMM.DEDD10$mean.discoveries$if.DESeq, 
     DE.results.TMM.DEDD10$mean.fdr$if.DESeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='DESeq2, TMM')
plot(DE.results.TMM.DEDD10$mean.discoveries$voom, 
     DE.results.TMM.DEDD10$mean.fdr$voom, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='limma, TMM')
plot(DE.results.TMM.DEDD10$mean.discoveries$notrend.DSS, 
     DE.results.TMM.DEDD10$mean.fdr$notrend.DSS, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='DSS, TMM')
plot(DE.results.TMM.DEDD10$mean.discoveries$baySeq, 
     DE.results.TMM.DEDD10$mean.fdr$baySeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='baySeq, TMM')
plot(DE.results.TMM.DEDD10$mean.discoveries$mean.zi.MDSeq, 
     DE.results.TMM.DEDD10$mean.fdr$mean.zi.MDSeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='MDSeq +ZI, TMM')
plot(DE.results.TMM.DEDD10$mean.discoveries$mean.nozi.MDSeq, 
     DE.results.TMM.DEDD10$mean.fdr$mean.nozi.MDSeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='MDSeq -ZI')
plot(DE.results.TMM.DEDD10$mean.discoveries$mean.expHM, 
     DE.results.TMM.DEDD10$mean.fdr$mean.expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='expHM, TMM')
plot(DE.results.TMM.DEDD10$mean.discoveries$lmean.expHM, 
     DE.results.TMM.DEDD10$mean.fdr$lmean.expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='expHM log, TMM')
plot(DE.results.TMM.DEDD10$mean.discoveries$mean.lnHM, 
     DE.results.TMM.DEDD10$mean.fdr$mean.lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='lnHM, TMM')
plot(DE.results.TMM.DEDD10$mean.discoveries$lmean.lnHM, 
     DE.results.TMM.DEDD10$mean.fdr$lmean.lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='lnHM log, TMM')
plot.new()
plot(DE.results.DESeq.DEDD10$mean.discoveries$ql.edgeR, 
     DE.results.DESeq.DEDD10$mean.fdr$ql.edgeR, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='edgeR QL, DESeq norm')
plot(DE.results.DESeq.DEDD10$mean.discoveries$lr.edgeR, 
     DE.results.DESeq.DEDD10$mean.fdr$lr.edgeR, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='edgeR LR, DESeq norm')
plot(DE.results.DESeq.DEDD10$mean.discoveries$et.edgeR, 
     DE.results.DESeq.DEDD10$mean.fdr$et.edgeR, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='edgeR ET, DESeq norm')
plot(DE.results.DESeq.DEDD10$mean.discoveries$if.DESeq, 
     DE.results.DESeq.DEDD10$mean.fdr$if.DESeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='DESeq2, DESeq norm')
plot(DE.results.DESeq.DEDD10$mean.discoveries$voom, 
     DE.results.DESeq.DEDD10$mean.fdr$voom, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='limma, DESeq norm')
plot(DE.results.DESeq.DEDD10$mean.discoveries$notrend.DSS, 
     DE.results.DESeq.DEDD10$mean.fdr$notrend.DSS, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='DSS, DESeq norm')
plot(DE.results.DESeq.DEDD10$mean.discoveries$baySeq, 
     DE.results.DESeq.DEDD10$mean.fdr$baySeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='baySeq, DESeq norm')
plot(DE.results.DESeq.DEDD10$mean.discoveries$mean.zi.MDSeq, 
     DE.results.DESeq.DEDD10$mean.fdr$mean.zi.MDSeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='MDSeq +ZI, DESeq norm')
plot(DE.results.DESeq.DEDD10$mean.discoveries$mean.nozi.MDSeq, 
     DE.results.DESeq.DEDD10$mean.fdr$mean.nozi.MDSeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='MDSeq -ZI')
plot(DE.results.DESeq.DEDD10$mean.discoveries$mean.expHM, 
     DE.results.DESeq.DEDD10$mean.fdr$mean.expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='expHM, DESeq norm')
plot(DE.results.DESeq.DEDD10$mean.discoveries$lmean.expHM, 
     DE.results.DESeq.DEDD10$mean.fdr$lmean.expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='expHM log, DESeq norm')
plot(DE.results.DESeq.DEDD10$mean.discoveries$mean.lnHM, 
     DE.results.DESeq.DEDD10$mean.fdr$mean.lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='lnHM, DESeq norm')
plot(DE.results.DESeq.DEDD10$mean.discoveries$lmean.lnHM, 
     DE.results.DESeq.DEDD10$mean.fdr$lmean.lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='lnHM log, DESeq norm')
plot.new()

## DEDD20 ##
par(mfrow=c(4,7), mar=c(2,2,1,1), mgp=c(3,0.7,0))
plot(DE.results.TMM.DEDD20$mean.discoveries$ql.edgeR, 
     DE.results.TMM.DEDD20$mean.fdr$ql.edgeR, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='edgeR QL, TMM')
plot(DE.results.TMM.DEDD20$mean.discoveries$lr.edgeR, 
     DE.results.TMM.DEDD20$mean.fdr$lr.edgeR, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='edgeR LR, TMM')
plot(DE.results.TMM.DEDD20$mean.discoveries$et.edgeR, 
     DE.results.TMM.DEDD20$mean.fdr$et.edgeR, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='edgeR ET, TMM')
plot(DE.results.TMM.DEDD20$mean.discoveries$if.DESeq, 
     DE.results.TMM.DEDD20$mean.fdr$if.DESeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='DESeq2, TMM')
plot(DE.results.TMM.DEDD20$mean.discoveries$voom, 
     DE.results.TMM.DEDD20$mean.fdr$voom, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='limma, TMM')
plot(DE.results.TMM.DEDD20$mean.discoveries$notrend.DSS, 
     DE.results.TMM.DEDD20$mean.fdr$notrend.DSS, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='DSS, TMM')
plot(DE.results.TMM.DEDD20$mean.discoveries$baySeq, 
     DE.results.TMM.DEDD20$mean.fdr$baySeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='baySeq, TMM')
plot(DE.results.TMM.DEDD20$mean.discoveries$mean.zi.MDSeq, 
     DE.results.TMM.DEDD20$mean.fdr$mean.zi.MDSeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='MDSeq +ZI, TMM')
plot(DE.results.TMM.DEDD20$mean.discoveries$mean.nozi.MDSeq, 
     DE.results.TMM.DEDD20$mean.fdr$mean.nozi.MDSeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='MDSeq -ZI')
plot(DE.results.TMM.DEDD20$mean.discoveries$mean.expHM, 
     DE.results.TMM.DEDD20$mean.fdr$mean.expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='expHM, TMM')
plot(DE.results.TMM.DEDD20$mean.discoveries$lmean.expHM, 
     DE.results.TMM.DEDD20$mean.fdr$lmean.expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='expHM log, TMM')
plot(DE.results.TMM.DEDD20$mean.discoveries$mean.lnHM, 
     DE.results.TMM.DEDD20$mean.fdr$mean.lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='lnHM, TMM')
plot(DE.results.TMM.DEDD20$mean.discoveries$lmean.lnHM, 
     DE.results.TMM.DEDD20$mean.fdr$lmean.lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='lnHM log, TMM')
plot.new()
plot(DE.results.DESeq.DEDD20$mean.discoveries$ql.edgeR, 
     DE.results.DESeq.DEDD20$mean.fdr$ql.edgeR, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='edgeR QL, DESeq norm')
plot(DE.results.DESeq.DEDD20$mean.discoveries$lr.edgeR, 
     DE.results.DESeq.DEDD20$mean.fdr$lr.edgeR, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='edgeR LR, DESeq norm')
plot(DE.results.DESeq.DEDD20$mean.discoveries$et.edgeR, 
     DE.results.DESeq.DEDD20$mean.fdr$et.edgeR, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='edgeR ET, DESeq norm')
plot(DE.results.DESeq.DEDD20$mean.discoveries$if.DESeq, 
     DE.results.DESeq.DEDD20$mean.fdr$if.DESeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='DESeq2, DESeq norm')
plot(DE.results.DESeq.DEDD20$mean.discoveries$voom, 
     DE.results.DESeq.DEDD20$mean.fdr$voom, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='limma, DESeq norm')
plot(DE.results.DESeq.DEDD20$mean.discoveries$notrend.DSS, 
     DE.results.DESeq.DEDD20$mean.fdr$notrend.DSS, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='DSS, DESeq norm')
plot(DE.results.DESeq.DEDD20$mean.discoveries$baySeq, 
     DE.results.DESeq.DEDD20$mean.fdr$baySeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='baySeq, DESeq norm')
plot(DE.results.DESeq.DEDD20$mean.discoveries$mean.zi.MDSeq, 
     DE.results.DESeq.DEDD20$mean.fdr$mean.zi.MDSeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='MDSeq +ZI, DESeq norm')
plot(DE.results.DESeq.DEDD20$mean.discoveries$mean.nozi.MDSeq, 
     DE.results.DESeq.DEDD20$mean.fdr$mean.nozi.MDSeq, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='MDSeq -ZI')
plot(DE.results.DESeq.DEDD20$mean.discoveries$mean.expHM, 
     DE.results.DESeq.DEDD20$mean.fdr$mean.expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='expHM, DESeq norm')
plot(DE.results.DESeq.DEDD20$mean.discoveries$lmean.expHM, 
     DE.results.DESeq.DEDD20$mean.fdr$lmean.expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='expHM log, DESeq norm')
plot(DE.results.DESeq.DEDD20$mean.discoveries$mean.lnHM, 
     DE.results.DESeq.DEDD20$mean.fdr$mean.lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='lnHM, DESeq norm')
plot(DE.results.DESeq.DEDD20$mean.discoveries$lmean.lnHM, 
     DE.results.DESeq.DEDD20$mean.fdr$lmean.lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomright", bty='n', legend='lnHM log, DESeq norm')
plot.new()

###################################
#### Differential distribution ####
###################################

## DD2, DE2, DEDD2 ##
par(mfcol=c(2,6), mar=c(2,2,1,1), mgp=c(3,0.7,0))
plot(DEDD.results.TMM.DD2$mean.discoveries$expHM, 
     DEDD.results.TMM.DD2$mean.fdr$expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomleft", bty='n', legend=c('Disp', 'expHM, TMM'))
plot(DEDD.results.DESeq.DD2$mean.discoveries$expHM, 
     DEDD.results.DESeq.DD2$mean.fdr$expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomleft", bty='n', legend=c('Disp', 'expHM, DESeq norm'))
plot(DEDD.results.TMM.DD2$mean.discoveries$lnHM, 
     DEDD.results.TMM.DD2$mean.fdr$lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomleft", bty='n', legend=c('Disp', 'lnHM, TMM'))
plot(DEDD.results.DESeq.DD2$mean.discoveries$lnHM, 
     DEDD.results.DESeq.DD2$mean.fdr$lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomleft", bty='n', legend=c('Disp', 'lnHM, DESeq norm'))
plot(DEDD.results.TMM.DE2$mean.discoveries$expHM, 
     DEDD.results.TMM.DE2$mean.fdr$expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomleft", bty='n', legend=c('Mean', 'expHM, TMM'))
plot(DEDD.results.DESeq.DE2$mean.discoveries$expHM, 
     DEDD.results.DESeq.DE2$mean.fdr$expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomleft", bty='n', legend=c('Mean', 'expHM, DESeq norm'))
plot(DEDD.results.TMM.DE2$mean.discoveries$lnHM, 
     DEDD.results.TMM.DE2$mean.fdr$lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomleft", bty='n', legend=c('Mean', 'lnHM, TMM'))
plot(DEDD.results.DESeq.DE2$mean.discoveries$lnHM, 
     DEDD.results.DESeq.DE2$mean.fdr$lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomleft", bty='n', legend=c('Mean', 'lnHM, DESeq norm'))
plot(DEDD.results.TMM.DEDD2$mean.discoveries$expHM, 
     DEDD.results.TMM.DEDD2$mean.fdr$expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomleft", bty='n', legend=c('Both', 'expHM, TMM'))
plot(DEDD.results.DESeq.DEDD2$mean.discoveries$expHM, 
     DEDD.results.DESeq.DEDD2$mean.fdr$expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomleft", bty='n', legend=c('Both', 'expHM, DESeq norm'))
plot(DEDD.results.TMM.DEDD2$mean.discoveries$lnHM, 
     DEDD.results.TMM.DEDD2$mean.fdr$lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomleft", bty='n', legend=c('Both', 'lnHM, TMM'))
plot(DEDD.results.DESeq.DEDD2$mean.discoveries$lnHM, 
     DEDD.results.DESeq.DEDD2$mean.fdr$lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomleft", bty='n', legend=c('Both', 'lnHM, DESeq norm'))

## DD5, DE5, DEDD5 ##
par(mfcol=c(2,6), mar=c(2,2,1,1), mgp=c(3,0.7,0))
plot(DEDD.results.TMM.DD5$mean.discoveries$expHM, 
     DEDD.results.TMM.DD5$mean.fdr$expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomleft", bty='n', legend=c('Disp', 'expHM, TMM'))
plot(DEDD.results.DESeq.DD5$mean.discoveries$expHM, 
     DEDD.results.DESeq.DD5$mean.fdr$expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomleft", bty='n', legend=c('Disp', 'expHM, DESeq norm'))
plot(DEDD.results.TMM.DD5$mean.discoveries$lnHM, 
     DEDD.results.TMM.DD5$mean.fdr$lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomleft", bty='n', legend=c('Disp', 'lnHM, TMM'))
plot(DEDD.results.DESeq.DD5$mean.discoveries$lnHM, 
     DEDD.results.DESeq.DD5$mean.fdr$lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomleft", bty='n', legend=c('Disp', 'lnHM, DESeq norm'))
plot(DEDD.results.TMM.DE5$mean.discoveries$expHM, 
     DEDD.results.TMM.DE5$mean.fdr$expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("topleft", bty='n', legend=c('Mean', 'expHM, TMM'))
plot(DEDD.results.DESeq.DE5$mean.discoveries$expHM, 
     DEDD.results.DESeq.DE5$mean.fdr$expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("topleft", bty='n', legend=c('Mean', 'expHM, DESeq norm'))
plot(DEDD.results.TMM.DE5$mean.discoveries$lnHM, 
     DEDD.results.TMM.DE5$mean.fdr$lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("topleft", bty='n', legend=c('Mean', 'lnHM, TMM'))
plot(DEDD.results.DESeq.DE5$mean.discoveries$lnHM, 
     DEDD.results.DESeq.DE5$mean.fdr$lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("topleft", bty='n', legend=c('Mean', 'lnHM, DESeq norm'))
plot(DEDD.results.TMM.DEDD5$mean.discoveries$expHM, 
     DEDD.results.TMM.DEDD5$mean.fdr$expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("topleft", bty='n', legend=c('Both', 'expHM, TMM'))
plot(DEDD.results.DESeq.DEDD5$mean.discoveries$expHM, 
     DEDD.results.DESeq.DEDD5$mean.fdr$expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("topleft", bty='n', legend=c('Both', 'expHM, DESeq norm'))
plot(DEDD.results.TMM.DEDD5$mean.discoveries$lnHM, 
     DEDD.results.TMM.DEDD5$mean.fdr$lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("topleft", bty='n', legend=c('Both', 'lnHM, TMM'))
plot(DEDD.results.DESeq.DEDD5$mean.discoveries$lnHM, 
     DEDD.results.DESeq.DEDD5$mean.fdr$lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("topleft", bty='n', legend=c('Both', 'lnHM, DESeq norm'))

## DD10, DE10, DEDD10 ##
par(mfcol=c(2,6), mar=c(2,2,1,1), mgp=c(3,0.7,0))
plot(DEDD.results.TMM.DD10$mean.discoveries$expHM, 
     DEDD.results.TMM.DD10$mean.fdr$expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomleft", bty='n', legend=c('Disp', 'expHM, TMM'))
plot(DEDD.results.DESeq.DD10$mean.discoveries$expHM, 
     DEDD.results.DESeq.DD10$mean.fdr$expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomleft", bty='n', legend=c('Disp', 'expHM, DESeq norm'))
plot(DEDD.results.TMM.DD10$mean.discoveries$lnHM, 
     DEDD.results.TMM.DD10$mean.fdr$lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomleft", bty='n', legend=c('Disp', 'lnHM, TMM'))
plot(DEDD.results.DESeq.DD10$mean.discoveries$lnHM, 
     DEDD.results.DESeq.DD10$mean.fdr$lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("bottomleft", bty='n', legend=c('Disp', 'lnHM, DESeq norm'))
plot(DEDD.results.TMM.DE10$mean.discoveries$expHM, 
     DEDD.results.TMM.DE10$mean.fdr$expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("topleft", bty='n', legend=c('Mean', 'expHM, TMM'))
plot(DEDD.results.DESeq.DE10$mean.discoveries$expHM, 
     DEDD.results.DESeq.DE10$mean.fdr$expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("topleft", bty='n', legend=c('Mean', 'expHM, DESeq norm'))
plot(DEDD.results.TMM.DE10$mean.discoveries$lnHM, 
     DEDD.results.TMM.DE10$mean.fdr$lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("topleft", bty='n', legend=c('Mean', 'lnHM, TMM'))
plot(DEDD.results.DESeq.DE10$mean.discoveries$lnHM, 
     DEDD.results.DESeq.DE10$mean.fdr$lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("topleft", bty='n', legend=c('Mean', 'lnHM, DESeq norm'))
plot(DEDD.results.TMM.DEDD10$mean.discoveries$expHM, 
     DEDD.results.TMM.DEDD10$mean.fdr$expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("topleft", bty='n', legend=c('Both', 'expHM, TMM'))
plot(DEDD.results.DESeq.DEDD10$mean.discoveries$expHM, 
     DEDD.results.DESeq.DEDD10$mean.fdr$expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("topleft", bty='n', legend=c('Both', 'expHM, DESeq norm'))
plot(DEDD.results.TMM.DEDD10$mean.discoveries$lnHM, 
     DEDD.results.TMM.DEDD10$mean.fdr$lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("topleft", bty='n', legend=c('Both', 'lnHM, TMM'))
plot(DEDD.results.DESeq.DEDD10$mean.discoveries$lnHM, 
     DEDD.results.DESeq.DEDD10$mean.fdr$lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("topleft", bty='n', legend=c('Both', 'lnHM, DESeq norm'))

## DD20, DE20, DEDD20 ##
par(mfcol=c(2,6), mar=c(2,2,1,1), mgp=c(3,0.7,0))
plot(DEDD.results.TMM.DD20$mean.discoveries$expHM, 
     DEDD.results.TMM.DD20$mean.fdr$expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("topleft", bty='n', legend=c('Disp', 'expHM, TMM'))
plot(DEDD.results.DESeq.DD20$mean.discoveries$expHM, 
     DEDD.results.DESeq.DD20$mean.fdr$expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("topleft", bty='n', legend=c('Disp', 'expHM, DESeq norm'))
plot(DEDD.results.TMM.DD20$mean.discoveries$lnHM, 
     DEDD.results.TMM.DD20$mean.fdr$lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("topleft", bty='n', legend=c('Disp', 'lnHM, TMM'))
plot(DEDD.results.DESeq.DD20$mean.discoveries$lnHM, 
     DEDD.results.DESeq.DD20$mean.fdr$lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("topleft", bty='n', legend=c('Disp', 'lnHM, DESeq norm'))
plot(DEDD.results.TMM.DE20$mean.discoveries$expHM, 
     DEDD.results.TMM.DE20$mean.fdr$expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("topleft", bty='n', legend=c('Mean', 'expHM, TMM'))
plot(DEDD.results.DESeq.DE20$mean.discoveries$expHM, 
     DEDD.results.DESeq.DE20$mean.fdr$expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("topleft", bty='n', legend=c('Mean', 'expHM, DESeq norm'))
plot(DEDD.results.TMM.DE20$mean.discoveries$lnHM, 
     DEDD.results.TMM.DE20$mean.fdr$lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("topleft", bty='n', legend=c('Mean', 'lnHM, TMM'))
plot(DEDD.results.DESeq.DE20$mean.discoveries$lnHM, 
     DEDD.results.DESeq.DE20$mean.fdr$lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("topleft", bty='n', legend=c('Mean', 'lnHM, DESeq norm'))
plot(DEDD.results.TMM.DEDD20$mean.discoveries$expHM, 
     DEDD.results.TMM.DEDD20$mean.fdr$expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("topleft", bty='n', legend=c('Both', 'expHM, TMM'))
plot(DEDD.results.DESeq.DEDD20$mean.discoveries$expHM, 
     DEDD.results.DESeq.DEDD20$mean.fdr$expHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("topleft", bty='n', legend=c('Both', 'expHM, DESeq norm'))
plot(DEDD.results.TMM.DEDD20$mean.discoveries$lnHM, 
     DEDD.results.TMM.DEDD20$mean.fdr$lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("topleft", bty='n', legend=c('Both', 'lnHM, TMM'))
plot(DEDD.results.DESeq.DEDD20$mean.discoveries$lnHM, 
     DEDD.results.DESeq.DEDD20$mean.fdr$lnHM, type='l', ylim=c(0,1), 
     xlab='Number of discoveries', ylab='FDR')
legend("topleft", bty='n', legend=c('Both', 'lnHM, DESeq norm'))


## Compare variants DE within methods ####
par(mfrow=c(4,4), mar=c(2,2,1,1), mgp=c(3,0.7,0))
for (i in c('DE2', 'DE5', 'DE10', 'DE20')) {
  plot(get(paste0('DE.results.TMM.',i))$mean.discoveries$ql.edgeR, 
       get(paste0('DE.results.TMM.',i))$mean.fdr$ql.edgeR, 
       type='l', ylim=c(0,0.8), xlim=c(0,2000), col=col_vector[1])
  lines(get(paste0('DE.results.TMM.',i))$mean.discoveries$lr.edgeR, 
        get(paste0('DE.results.TMM.',i))$mean.fdr$lr.edgeR, 
        type='l', col=col_vector[2])
  lines(get(paste0('DE.results.TMM.',i))$mean.discoveries$et.edgeR, 
        get(paste0('DE.results.TMM.',i))$mean.fdr$et.edgeR, 
        type='l', col=col_vector[3])
  legend("bottomright", bty='n', legend=c('QL','LR','ET'), col=col_vector[1:3], lty=1)
}
for (i in c('DEDD2', 'DEDD5', 'DEDD10', 'DEDD20')) {
  plot(get(paste0('DE.results.TMM.',i))$mean.discoveries$ql.edgeR, 
       get(paste0('DE.results.TMM.',i))$mean.fdr$ql.edgeR, 
       type='l', ylim=c(0,0.8), xlim=c(0,2000), col=col_vector[1])
  lines(get(paste0('DE.results.TMM.',i))$mean.discoveries$lr.edgeR, 
        get(paste0('DE.results.TMM.',i))$mean.fdr$lr.edgeR, 
        type='l', col=col_vector[2])
  legend("topleft", bty='n', legend=c('QL','LR','ET'), col=col_vector[1:3], lty=1)
}
# for (i in c('DE2', 'DE5', 'DE10', 'DE20')) {
#   plot(get(paste0('DE.results.DESeq.',i))$mean.discoveries$mean.zi.MDSeq, 
#        get(paste0('DE.results.DESeq.',i))$mean.fdr$mean.zi.MDSeq, 
#        type='l', ylim=c(0,0.8), xlim=c(0,2000), col=col_vector[1])
#   lines(get(paste0('DE.results.DESeq.',i))$mean.discoveries$mean.nozi.MDSeq, 
#         get(paste0('DE.results.DESeq.',i))$mean.fdr$mean.nozi.MDSeq, 
#         type='l', col=col_vector[2])
#   legend("bottomright", bty='n', legend=c('+ZI','-ZI'), col=col_vector[1:2], lty=1)
# }
# for (i in c('DEDD2', 'DEDD5', 'DEDD10', 'DEDD20')) {
#   plot(get(paste0('DE.results.DESeq.',i))$mean.discoveries$mean.zi.MDSeq, 
#        get(paste0('DE.results.DESeq.',i))$mean.fdr$mean.zi.MDSeq, 
#        type='l', ylim=c(0,0.8), xlim=c(0,2000), col=col_vector[1])
#   lines(get(paste0('DE.results.DESeq.',i))$mean.discoveries$mean.nozi.MDSeq, 
#         get(paste0('DE.results.DESeq.',i))$mean.fdr$mean.nozi.MDSeq, 
#         type='l', col=col_vector[2])
#   legend("topleft", bty='n', legend=c('+ZI','-ZI'), col=col_vector[1:2], lty=1)
# }
for (i in c('DE2', 'DE5', 'DE10', 'DE20')) {
  plot(get(paste0('DE.results.DESeq.',i))$mean.discoveries$mean.expHM, 
       get(paste0('DE.results.DESeq.',i))$mean.fdr$mean.expHM, 
       type='l', ylim=c(0,0.8), xlim=c(0,2000), col=col_vector[1])
  lines(get(paste0('DE.results.DESeq.',i))$mean.discoveries$lmean.expHM, 
        get(paste0('DE.results.DESeq.',i))$mean.fdr$lmean.expHM, 
        type='l', col=col_vector[2])
  lines(get(paste0('DE.results.DESeq.',i))$mean.discoveries$mean.lnHM, 
        get(paste0('DE.results.DESeq.',i))$mean.fdr$mean.lnHM, 
        type='l', col=col_vector[3])
  lines(get(paste0('DE.results.DESeq.',i))$mean.discoveries$lmean.lnHM, 
        get(paste0('DE.results.DESeq.',i))$mean.fdr$lmean.lnHM, 
        type='l', col=col_vector[4])
  legend("bottomright", bty='n', legend=c('expHM','expHM log','lnHM','lnHM log'), 
         col=col_vector[1:4], lty=1)
}
for (i in c('DEDD2', 'DEDD5', 'DEDD10', 'DEDD20')) {
  plot(get(paste0('DE.results.DESeq.',i))$mean.discoveries$mean.expHM, 
       get(paste0('DE.results.DESeq.',i))$mean.fdr$mean.expHM, 
       type='l', ylim=c(0,0.8), xlim=c(0,2000), col=col_vector[1])
  lines(get(paste0('DE.results.DESeq.',i))$mean.discoveries$lmean.expHM, 
        get(paste0('DE.results.DESeq.',i))$mean.fdr$lmean.expHM, 
        type='l', col=col_vector[2])
  lines(get(paste0('DE.results.DESeq.',i))$mean.discoveries$mean.lnHM, 
        get(paste0('DE.results.DESeq.',i))$mean.fdr$mean.lnHM, 
        type='l', col=col_vector[3])
  lines(get(paste0('DE.results.DESeq.',i))$mean.discoveries$lmean.lnHM, 
        get(paste0('DE.results.DESeq.',i))$mean.fdr$lmean.lnHM, 
        type='l', col=col_vector[4])
  legend("topleft", bty='n', legend=c('expHM','expHM log','lnHM','lnHM log'), 
         col=col_vector[1:4], lty=1)
}


## Compare all DE methods on same plots ####
par(mfrow=c(2,4), mar=c(2,2,1,1), mgp=c(3,0.7,0))
for (i in c('DE2', 'DE5', 'DE10', 'DE20')) {
  plot(get(paste0('DE.results.TMM.',i))$mean.discoveries$ql.edgeR, 
       get(paste0('DE.results.TMM.',i))$mean.fdr$ql.edgeR, 
       type='l', ylim=c(0,0.8), xlim=c(0,2000), col=col_vector[1])
  lines(get(paste0('DE.results.DESeq.',i))$mean.discoveries$if.DESeq, 
        get(paste0('DE.results.DESeq.',i))$mean.fdr$if.DESeq, 
        type='l', col=col_vector[2])
  lines(get(paste0('DE.results.TMM.',i))$mean.discoveries$voom, 
        get(paste0('DE.results.TMM.',i))$mean.fdr$voom, 
        type='l', col=col_vector[3])
  lines(get(paste0('DE.results.DESeq.',i))$mean.discoveries$notrend.DSS, 
        get(paste0('DE.results.DESeq.',i))$mean.fdr$notrend.DSS, 
        type='l', col=col_vector[4])
  lines(get(paste0('DE.results.DESeq.',i))$mean.discoveries$baySeq, 
        get(paste0('DE.results.DESeq.',i))$mean.fdr$baySeq, 
        type='l', col=col_vector[5])
  lines(get(paste0('DE.results.DESeq.',i))$mean.discoveries$mean.zi.MDSeq, 
        get(paste0('DE.results.DESeq.',i))$mean.fdr$mean.zi.MDSeq, 
        type='l', col=col_vector[6])
  lines(get(paste0('DE.results.DESeq.',i))$mean.discoveries$mean.lnHM, 
        get(paste0('DE.results.DESeq.',i))$mean.fdr$mean.lnHM, 
        type='l', col=col_vector[7])
  legend("bottomright", bty='n', col=col_vector[1:7], lty=1, ncol=2, 
         legend=c('edgeR', 'DESeq', 'limma', 'DSS', 'baySeq', 'MDSeq', 'HM'))
}
for (i in c('DEDD2', 'DEDD5', 'DEDD10', 'DEDD20')) {
  plot(get(paste0('DE.results.TMM.',i))$mean.discoveries$ql.edgeR, 
       get(paste0('DE.results.TMM.',i))$mean.fdr$ql.edgeR, 
       type='l', ylim=c(0,0.8), xlim=c(0,2000), col=col_vector[1])
  lines(get(paste0('DE.results.DESeq.',i))$mean.discoveries$if.DESeq, 
        get(paste0('DE.results.DESeq.',i))$mean.fdr$if.DESeq, 
        type='l', col=col_vector[2])
  lines(get(paste0('DE.results.TMM.',i))$mean.discoveries$voom, 
        get(paste0('DE.results.TMM.',i))$mean.fdr$voom, 
        type='l', col=col_vector[3])
  lines(get(paste0('DE.results.DESeq.',i))$mean.discoveries$notrend.DSS, 
        get(paste0('DE.results.DESeq.',i))$mean.fdr$notrend.DSS, 
        type='l', col=col_vector[4])
  lines(get(paste0('DE.results.DESeq.',i))$mean.discoveries$baySeq, 
        get(paste0('DE.results.DESeq.',i))$mean.fdr$baySeq, 
        type='l', col=col_vector[5])
  lines(get(paste0('DE.results.DESeq.',i))$mean.discoveries$mean.zi.MDSeq, 
        get(paste0('DE.results.DESeq.',i))$mean.fdr$mean.zi.MDSeq, 
        type='l', col=col_vector[6])
  lines(get(paste0('DE.results.DESeq.',i))$mean.discoveries$mean.lnHM, 
        get(paste0('DE.results.DESeq.',i))$mean.fdr$mean.lnHM, 
        type='l', col=col_vector[7])
  legend("topleft", bty='n', col=col_vector[1:7], lty=1, ncol=2, 
         legend=c('edgeR', 'DESeq', 'limma', 'DSS', 'baySeq', 'MDSeq', 'HM'))
}

# Compare all DE methods with different scales for each scenario to best show differences ####
par(mfrow=c(2,4), mar=c(2,2,1,1), mgp=c(3,0.7,0))
# DE2
plot(DE.results.TMM.DE2$mean.discoveries$ql.edgeR, 
     DE.results.TMM.DE2$mean.fdr$ql.edgeR, 
     type='l', ylim=c(0,0.6), xlim=c(0,250), col=col_vector[1])
lines(DE.results.DESeq.DE2$mean.discoveries$if.DESeq, 
      DE.results.DESeq.DE2$mean.fdr$if.DESeq, 
      type='l', col=col_vector[2])
lines(DE.results.TMM.DE2$mean.discoveries$voom, 
      DE.results.TMM.DE2$mean.fdr$voom, 
      type='l', col=col_vector[3])
lines(DE.results.DESeq.DE2$mean.discoveries$notrend.DSS, 
      DE.results.DESeq.DE2$mean.fdr$notrend.DSS, 
      type='l', col=col_vector[4])
lines(DE.results.DESeq.DE2$mean.discoveries$baySeq, 
      DE.results.DESeq.DE2$mean.fdr$baySeq, 
      type='l', col=col_vector[5])
lines(DE.results.DESeq.DE2$mean.discoveries$mean.zi.MDSeq, 
      DE.results.DESeq.DE2$mean.fdr$mean.zi.MDSeq, 
      type='l', col=col_vector[6])
lines(DE.results.DESeq.DE2$mean.discoveries$mean.lnHM, 
      DE.results.DESeq.DE2$mean.fdr$mean.lnHM, 
      type='l', col=col_vector[7])
legend("bottomright", bty='n', col=col_vector[1:7], lty=1, 
       legend=c('edgeR', 'DESeq', 'limma', 'DSS', 'baySeq', 'MDSeq', 'HM'))
# DE5
plot(DE.results.TMM.DE5$mean.discoveries$ql.edgeR,
     DE.results.TMM.DE5$mean.fdr$ql.edgeR, 
     type='l', ylim=c(0,0.06), xlim=c(0,300), col=col_vector[1])
lines(DE.results.DESeq.DE5$mean.discoveries$if.DESeq, 
      DE.results.DESeq.DE5$mean.fdr$if.DESeq, 
      type='l', col=col_vector[2])
lines(DE.results.TMM.DE5$mean.discoveries$voom, 
      DE.results.TMM.DE5$mean.fdr$voom, 
      type='l', col=col_vector[3])
lines(DE.results.DESeq.DE5$mean.discoveries$notrend.DSS, 
      DE.results.DESeq.DE5$mean.fdr$notrend.DSS, 
      type='l', col=col_vector[4])
lines(DE.results.DESeq.DE5$mean.discoveries$baySeq, 
      DE.results.DESeq.DE5$mean.fdr$baySeq, 
      type='l', col=col_vector[5])
lines(DE.results.DESeq.DE5$mean.discoveries$mean.zi.MDSeq, 
      DE.results.DESeq.DE5$mean.fdr$mean.zi.MDSeq, 
      type='l', col=col_vector[6])
lines(DE.results.DESeq.DE5$mean.discoveries$mean.lnHM, 
      DE.results.DESeq.DE5$mean.fdr$mean.lnHM, 
      type='l', col=col_vector[7])
abline(h=0.05, col='lightgrey')
legend("bottomright", bty='n', col=col_vector[1:7], lty=1, 
       legend=c('edgeR', 'DESeq', 'limma', 'DSS', 'baySeq', 'MDSeq', 'HM'))
# DE10
plot(DE.results.TMM.DE10$mean.discoveries$ql.edgeR,
     DE.results.TMM.DE10$mean.fdr$ql.edgeR, 
     type='l', ylim=c(0,0.06), xlim=c(300,550), col=col_vector[1])
lines(DE.results.DESeq.DE10$mean.discoveries$if.DESeq, 
      DE.results.DESeq.DE10$mean.fdr$if.DESeq, 
      type='l', col=col_vector[2])
lines(DE.results.TMM.DE10$mean.discoveries$voom, 
      DE.results.TMM.DE10$mean.fdr$voom, 
      type='l', col=col_vector[3])
lines(DE.results.DESeq.DE10$mean.discoveries$notrend.DSS, 
      DE.results.DESeq.DE10$mean.fdr$notrend.DSS, 
      type='l', col=col_vector[4])
lines(DE.results.DESeq.DE10$mean.discoveries$baySeq, 
      DE.results.DESeq.DE10$mean.fdr$baySeq, 
      type='l', col=col_vector[5])
lines(DE.results.DESeq.DE10$mean.discoveries$mean.zi.MDSeq, 
      DE.results.DESeq.DE10$mean.fdr$mean.zi.MDSeq, 
      type='l', col=col_vector[6])
lines(DE.results.DESeq.DE10$mean.discoveries$mean.lnHM, 
      DE.results.DESeq.DE10$mean.fdr$mean.lnHM, 
      type='l', col=col_vector[7])
abline(h=0.05, col='lightgrey')
legend("bottomright", bty='n', col=col_vector[1:7], lty=1, 
       legend=c('edgeR', 'DESeq', 'limma', 'DSS', 'baySeq', 'MDSeq', 'HM'))
# DE20
plot(DE.results.TMM.DE20$mean.discoveries$ql.edgeR, 
     DE.results.TMM.DE20$mean.fdr$ql.edgeR, 
     type='l', ylim=c(0,0.06), xlim=c(500,700), col=col_vector[1])
lines(DE.results.DESeq.DE20$mean.discoveries$if.DESeq, 
      DE.results.DESeq.DE20$mean.fdr$if.DESeq, 
      type='l', col=col_vector[2])
lines(DE.results.TMM.DE20$mean.discoveries$voom, 
      DE.results.TMM.DE20$mean.fdr$voom, 
      type='l', col=col_vector[3])
lines(DE.results.DESeq.DE20$mean.discoveries$notrend.DSS, 
      DE.results.DESeq.DE20$mean.fdr$notrend.DSS, 
      type='l', col=col_vector[4])
lines(DE.results.DESeq.DE20$mean.discoveries$baySeq, 
      DE.results.DESeq.DE20$mean.fdr$baySeq, 
      type='l', col=col_vector[5])
lines(DE.results.DESeq.DE20$mean.discoveries$mean.zi.MDSeq, 
      DE.results.DESeq.DE20$mean.fdr$mean.zi.MDSeq, 
      type='l', col=col_vector[6])
lines(DE.results.DESeq.DE20$mean.discoveries$mean.lnHM, 
      DE.results.DESeq.DE20$mean.fdr$mean.lnHM, 
      type='l', col=col_vector[7])
abline(h=0.05, col='lightgrey')
legend("bottomright", bty='n', col=col_vector[1:7], lty=1, 
       legend=c('edgeR', 'DESeq', 'limma', 'DSS', 'baySeq', 'MDSeq', 'HM'))
# DEDD2
plot(DE.results.TMM.DEDD2$mean.discoveries$ql.edgeR, 
     DE.results.TMM.DEDD2$mean.fdr$ql.edgeR, 
     type='l', ylim=c(0,0.4), xlim=c(0,300), col=col_vector[1])
lines(DE.results.DESeq.DEDD2$mean.discoveries$if.DESeq, 
      DE.results.DESeq.DEDD2$mean.fdr$if.DESeq, 
      type='l', col=col_vector[2])
lines(DE.results.TMM.DEDD2$mean.discoveries$voom, 
      DE.results.TMM.DEDD2$mean.fdr$voom, 
      type='l', col=col_vector[3])
lines(DE.results.DESeq.DEDD2$mean.discoveries$notrend.DSS, 
      DE.results.DESeq.DEDD2$mean.fdr$notrend.DSS, 
      type='l', col=col_vector[4])
lines(DE.results.DESeq.DEDD2$mean.discoveries$baySeq, 
      DE.results.DESeq.DEDD2$mean.fdr$baySeq, 
      type='l', col=col_vector[5])
lines(DE.results.DESeq.DEDD2$mean.discoveries$mean.zi.MDSeq, 
      DE.results.DESeq.DEDD2$mean.fdr$mean.zi.MDSeq, 
      type='l', col=col_vector[6])
lines(DE.results.DESeq.DEDD2$mean.discoveries$mean.lnHM, 
      DE.results.DESeq.DEDD2$mean.fdr$mean.lnHM, 
      type='l', col=col_vector[7])
legend("bottomright", bty='n', col=col_vector[1:7], lty=1, 
       legend=c('edgeR', 'DESeq', 'limma', 'DSS', 'baySeq', 'MDSeq', 'HM'))
# DEDD5
plot(DE.results.TMM.DEDD5$mean.discoveries$ql.edgeR,
     DE.results.TMM.DEDD5$mean.fdr$ql.edgeR, 
     type='l', ylim=c(0,0.06), xlim=c(0,650), col=col_vector[1])
lines(DE.results.DESeq.DEDD5$mean.discoveries$if.DESeq, 
      DE.results.DESeq.DEDD5$mean.fdr$if.DESeq, 
      type='l', col=col_vector[2])
lines(DE.results.TMM.DEDD5$mean.discoveries$voom, 
      DE.results.TMM.DEDD5$mean.fdr$voom, 
      type='l', col=col_vector[3])
lines(DE.results.DESeq.DEDD5$mean.discoveries$notrend.DSS, 
      DE.results.DESeq.DEDD5$mean.fdr$notrend.DSS, 
      type='l', col=col_vector[4])
lines(DE.results.DESeq.DEDD5$mean.discoveries$baySeq, 
      DE.results.DESeq.DEDD5$mean.fdr$baySeq, 
      type='l', col=col_vector[5])
lines(DE.results.DESeq.DEDD5$mean.discoveries$mean.zi.MDSeq, 
      DE.results.DESeq.DEDD5$mean.fdr$mean.zi.MDSeq, 
      type='l', col=col_vector[6])
lines(DE.results.DESeq.DEDD5$mean.discoveries$mean.lnHM, 
      DE.results.DESeq.DEDD5$mean.fdr$mean.lnHM, 
      type='l', col=col_vector[7])
abline(h=0.05, col='lightgrey')
legend("bottomright", bty='n', col=col_vector[1:7], lty=1, 
       legend=c('edgeR', 'DESeq', 'limma', 'DSS', 'baySeq', 'MDSeq', 'HM'))
# DEDD10
plot(DE.results.TMM.DEDD10$mean.discoveries$ql.edgeR,
     DE.results.TMM.DEDD10$mean.fdr$ql.edgeR, 
     type='l', ylim=c(0,0.06), xlim=c(600,1200), col=col_vector[1])
lines(DE.results.DESeq.DEDD10$mean.discoveries$if.DESeq, 
      DE.results.DESeq.DEDD10$mean.fdr$if.DESeq, 
      type='l', col=col_vector[2])
lines(DE.results.TMM.DEDD10$mean.discoveries$voom, 
      DE.results.TMM.DEDD10$mean.fdr$voom, 
      type='l', col=col_vector[3])
lines(DE.results.DESeq.DEDD10$mean.discoveries$notrend.DSS, 
      DE.results.DESeq.DEDD10$mean.fdr$notrend.DSS, 
      type='l', col=col_vector[4])
lines(DE.results.DESeq.DEDD10$mean.discoveries$baySeq, 
      DE.results.DESeq.DEDD10$mean.fdr$baySeq, 
      type='l', col=col_vector[5])
lines(DE.results.DESeq.DEDD10$mean.discoveries$mean.zi.MDSeq, 
      DE.results.DESeq.DEDD10$mean.fdr$mean.zi.MDSeq, 
      type='l', col=col_vector[6])
lines(DE.results.DESeq.DEDD10$mean.discoveries$mean.lnHM, 
      DE.results.DESeq.DEDD10$mean.fdr$mean.lnHM, 
      type='l', col=col_vector[7])
abline(h=0.05, col='lightgrey')
legend("bottomright", bty='n', col=col_vector[1:7], lty=1, 
       legend=c('edgeR', 'DESeq', 'limma', 'DSS', 'baySeq', 'MDSeq', 'HM'))
# DEDD20
plot(DE.results.TMM.DEDD20$mean.discoveries$ql.edgeR, 
     DE.results.TMM.DEDD20$mean.fdr$ql.edgeR, 
     type='l', ylim=c(0,0.06), xlim=c(800,1500), col=col_vector[1])
lines(DE.results.DESeq.DEDD20$mean.discoveries$if.DESeq, 
      DE.results.DESeq.DEDD20$mean.fdr$if.DESeq, 
      type='l', col=col_vector[2])
lines(DE.results.TMM.DEDD20$mean.discoveries$voom, 
      DE.results.TMM.DEDD20$mean.fdr$voom, 
      type='l', col=col_vector[3])
lines(DE.results.DESeq.DEDD20$mean.discoveries$notrend.DSS, 
      DE.results.DESeq.DEDD20$mean.fdr$notrend.DSS, 
      type='l', col=col_vector[4])
lines(DE.results.DESeq.DEDD20$mean.discoveries$baySeq, 
      DE.results.DESeq.DEDD20$mean.fdr$baySeq, 
      type='l', col=col_vector[5])
lines(DE.results.DESeq.DEDD20$mean.discoveries$mean.zi.MDSeq, 
      DE.results.DESeq.DEDD20$mean.fdr$mean.zi.MDSeq, 
      type='l', col=col_vector[6])
lines(DE.results.DESeq.DEDD20$mean.discoveries$mean.lnHM, 
      DE.results.DESeq.DEDD20$mean.fdr$mean.lnHM, 
      type='l', col=col_vector[7])
abline(h=0.05, col='lightgrey')
legend("bottomright", bty='n', col=col_vector[1:7], lty=1, 
       legend=c('edgeR', 'DESeq', 'limma', 'DSS', 'baySeq', 'MDSeq', 'HM'))


## Compare all differential dispersion methods on same plots ####
par(mfrow=c(2,4), mar=c(2,2,1,1), mgp=c(3,0.7,0))
for (i in c('DD2', 'DD5', 'DD10', 'DD20')) {
  plot(get(paste0('DD.results.DESeq.',i))$mean.discoveries$disp.zi.MDSeq, 
       get(paste0('DD.results.DESeq.',i))$mean.fdr$disp.zi.MDSeq, 
       type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1])
  lines(get(paste0('DD.results.DESeq.',i))$mean.discoveries$disp.nozi.MDSeq, 
        get(paste0('DD.results.DESeq.',i))$mean.fdr$disp.nozi.MDSeq, 
        type='l', col=col_vector[2])
  lines(get(paste0('DD.results.DESeq.',i))$mean.discoveries$disp.expHM,
        get(paste0('DD.results.DESeq.',i))$mean.fdr$disp.expHM,
        type='l', col=col_vector[3])
  lines(get(paste0('DD.results.DESeq.',i))$mean.discoveries$ldisp.expHM,
        get(paste0('DD.results.DESeq.',i))$mean.fdr$ldisp.expHM,
        type='l', col=col_vector[4])
  lines(get(paste0('DD.results.DESeq.',i))$mean.discoveries$disp.lnHM,
        get(paste0('DD.results.DESeq.',i))$mean.fdr$disp.lnHM,
        type='l', col=col_vector[5])
  lines(get(paste0('DD.results.DESeq.',i))$mean.discoveries$ldisp.lnHM,
        get(paste0('DD.results.DESeq.',i))$mean.fdr$ldisp.lnHM,
        type='l', col=col_vector[6])
  legend("bottomright", bty='n', col=col_vector[1:6], lty=1, ncol=3, 
         legend=c('MDSeq +ZI', 'MDSeq -ZI', 'expHM', 'expHM log', 'lnHM', 'lnHM log'))
}
for (i in c('DEDD2', 'DEDD5', 'DEDD10', 'DEDD20')) {
  plot(get(paste0('DD.results.DESeq.',i))$mean.discoveries$disp.zi.MDSeq, 
       get(paste0('DD.results.DESeq.',i))$mean.fdr$disp.zi.MDSeq, 
       type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1])
  lines(get(paste0('DD.results.DESeq.',i))$mean.discoveries$disp.nozi.MDSeq, 
        get(paste0('DD.results.DESeq.',i))$mean.fdr$disp.nozi.MDSeq, 
        type='l', col=col_vector[2])
  lines(get(paste0('DD.results.DESeq.',i))$mean.discoveries$disp.expHM,
        get(paste0('DD.results.DESeq.',i))$mean.fdr$disp.expHM,
        type='l', col=col_vector[3])
  lines(get(paste0('DD.results.DESeq.',i))$mean.discoveries$ldisp.expHM,
        get(paste0('DD.results.DESeq.',i))$mean.fdr$ldisp.expHM,
        type='l', col=col_vector[4])
  lines(get(paste0('DD.results.DESeq.',i))$mean.discoveries$disp.lnHM,
        get(paste0('DD.results.DESeq.',i))$mean.fdr$disp.lnHM,
        type='l', col=col_vector[5])
  lines(get(paste0('DD.results.DESeq.',i))$mean.discoveries$ldisp.lnHM,
        get(paste0('DD.results.DESeq.',i))$mean.fdr$ldisp.lnHM,
        type='l', col=col_vector[6])
  legend("topleft", bty='n', col=col_vector[1:6], lty=1, ncol=3, 
         legend=c('MDSeq +ZI', 'MDSeq -ZI', 'expHM', 'expHM log', 'lnHM', 'lnHM log'))
}


# Compare all differential dispersion methods with different scales for each scenario ####
par(mfrow=c(2,4), mar=c(2,2,1,1), mgp=c(3,0.7,0))
# DD2
plot(DD.results.DESeq.DD2$mean.discoveries$disp.zi.MDSeq, 
     DD.results.DESeq.DD2$mean.fdr$disp.zi.MDSeq, 
     type='l', ylim=c(0.8,0.95), xlim=c(0,500), col=col_vector[1])
lines(DD.results.DESeq.DD2$mean.discoveries$disp.nozi.MDSeq, 
      DD.results.DESeq.DD2$mean.fdr$disp.nozi.MDSeq, 
      type='l', col=col_vector[2])
lines(DD.results.DESeq.DD2$mean.discoveries$disp.expHM,
      DD.results.DESeq.DD2$mean.fdr$disp.expHM,
      type='l', col=col_vector[3])
lines(DD.results.DESeq.DD2$mean.discoveries$ldisp.expHM,
      DD.results.DESeq.DD2$mean.fdr$ldisp.expHM,
      type='l', col=col_vector[4])
lines(DD.results.DESeq.DD2$mean.discoveries$disp.lnHM,
      DD.results.DESeq.DD2$mean.fdr$disp.lnHM,
      type='l', col=col_vector[5])
lines(DD.results.DESeq.DD2$mean.discoveries$ldisp.lnHM,
      DD.results.DESeq.DD2$mean.fdr$ldisp.lnHM,
      type='l', col=col_vector[6])
legend("bottomright", bty='n', col=col_vector[1:6], lty=1, ncol=3, 
       legend=c('MDSeq +ZI', 'MDSeq -ZI', 'expHM', 'expHM log', 'lnHM', 'lnHM log'))
# DD5
plot(DD.results.DESeq.DD5$mean.discoveries$disp.zi.MDSeq, 
     DD.results.DESeq.DD5$mean.fdr$disp.zi.MDSeq, 
     type='l', ylim=c(0.6,0.9), xlim=c(0,500), col=col_vector[1])
lines(DD.results.DESeq.DD5$mean.discoveries$disp.nozi.MDSeq, 
      DD.results.DESeq.DD5$mean.fdr$disp.nozi.MDSeq, 
      type='l', col=col_vector[2])
lines(DD.results.DESeq.DD5$mean.discoveries$disp.expHM,
      DD.results.DESeq.DD5$mean.fdr$disp.expHM,
      type='l', col=col_vector[3])
lines(DD.results.DESeq.DD5$mean.discoveries$ldisp.expHM,
      DD.results.DESeq.DD5$mean.fdr$ldisp.expHM,
      type='l', col=col_vector[4])
lines(DD.results.DESeq.DD5$mean.discoveries$disp.lnHM,
      DD.results.DESeq.DD5$mean.fdr$disp.lnHM,
      type='l', col=col_vector[5])
lines(DD.results.DESeq.DD5$mean.discoveries$ldisp.lnHM,
      DD.results.DESeq.DD5$mean.fdr$ldisp.lnHM,
      type='l', col=col_vector[6])
legend("bottomright", bty='n', col=col_vector[1:6], lty=1, ncol=3, 
       legend=c('MDSeq +ZI', 'MDSeq -ZI', 'expHM', 'expHM log', 'lnHM', 'lnHM log'))
# DD10
plot(DD.results.DESeq.DD10$mean.discoveries$disp.zi.MDSeq, 
     DD.results.DESeq.DD10$mean.fdr$disp.zi.MDSeq, 
     type='l', ylim=c(0.3,0.8), xlim=c(0,500), col=col_vector[1])
lines(DD.results.DESeq.DD10$mean.discoveries$disp.nozi.MDSeq, 
      DD.results.DESeq.DD10$mean.fdr$disp.nozi.MDSeq, 
      type='l', col=col_vector[2])
lines(DD.results.DESeq.DD10$mean.discoveries$disp.expHM,
      DD.results.DESeq.DD10$mean.fdr$disp.expHM,
      type='l', col=col_vector[3])
lines(DD.results.DESeq.DD10$mean.discoveries$ldisp.expHM,
      DD.results.DESeq.DD10$mean.fdr$ldisp.expHM,
      type='l', col=col_vector[4])
lines(DD.results.DESeq.DD10$mean.discoveries$disp.lnHM,
      DD.results.DESeq.DD10$mean.fdr$disp.lnHM,
      type='l', col=col_vector[5])
lines(DD.results.DESeq.DD10$mean.discoveries$ldisp.lnHM,
      DD.results.DESeq.DD10$mean.fdr$ldisp.lnHM,
      type='l', col=col_vector[6])
legend("bottomright", bty='n', col=col_vector[1:6], lty=1, ncol=3, 
       legend=c('MDSeq +ZI', 'MDSeq -ZI', 'expHM', 'expHM log', 'lnHM', 'lnHM log'))
# DD20
plot(DD.results.DESeq.DD20$mean.discoveries$disp.zi.MDSeq, 
     DD.results.DESeq.DD20$mean.fdr$disp.zi.MDSeq, 
     type='l', ylim=c(0.1,0.6), xlim=c(0,500), col=col_vector[1])
lines(DD.results.DESeq.DD20$mean.discoveries$disp.nozi.MDSeq, 
      DD.results.DESeq.DD20$mean.fdr$disp.nozi.MDSeq, 
      type='l', col=col_vector[2])
lines(DD.results.DESeq.DD20$mean.discoveries$disp.expHM,
      DD.results.DESeq.DD20$mean.fdr$disp.expHM,
      type='l', col=col_vector[3])
lines(DD.results.DESeq.DD20$mean.discoveries$ldisp.expHM,
      DD.results.DESeq.DD20$mean.fdr$ldisp.expHM,
      type='l', col=col_vector[4])
lines(DD.results.DESeq.DD20$mean.discoveries$disp.lnHM,
      DD.results.DESeq.DD20$mean.fdr$disp.lnHM,
      type='l', col=col_vector[5])
lines(DD.results.DESeq.DD20$mean.discoveries$ldisp.lnHM,
      DD.results.DESeq.DD20$mean.fdr$ldisp.lnHM,
      type='l', col=col_vector[6])
legend("bottomright", bty='n', col=col_vector[1:6], lty=1, ncol=3, 
       legend=c('MDSeq +ZI', 'MDSeq -ZI', 'expHM', 'expHM log', 'lnHM', 'lnHM log'))
# DEDD2
plot(DD.results.DESeq.DEDD2$mean.discoveries$disp.zi.MDSeq, 
     DD.results.DESeq.DEDD2$mean.fdr$disp.zi.MDSeq, 
     type='l', ylim=c(0.5,0.9), xlim=c(0,500), col=col_vector[1])
lines(DD.results.DESeq.DEDD2$mean.discoveries$disp.nozi.MDSeq, 
      DD.results.DESeq.DEDD2$mean.fdr$disp.nozi.MDSeq, 
      type='l', col=col_vector[2])
lines(DD.results.DESeq.DEDD2$mean.discoveries$disp.expHM,
      DD.results.DESeq.DEDD2$mean.fdr$disp.expHM,
      type='l', col=col_vector[3])
lines(DD.results.DESeq.DEDD2$mean.discoveries$ldisp.expHM,
      DD.results.DESeq.DEDD2$mean.fdr$ldisp.expHM,
      type='l', col=col_vector[4])
lines(DD.results.DESeq.DEDD2$mean.discoveries$disp.lnHM,
      DD.results.DESeq.DEDD2$mean.fdr$disp.lnHM,
      type='l', col=col_vector[5])
lines(DD.results.DESeq.DEDD2$mean.discoveries$ldisp.lnHM,
      DD.results.DESeq.DEDD2$mean.fdr$ldisp.lnHM,
      type='l', col=col_vector[6])
legend("bottomright", bty='n', col=col_vector[1:6], lty=1, ncol=3, 
       legend=c('MDSeq +ZI', 'MDSeq -ZI', 'expHM', 'expHM log', 'lnHM', 'lnHM log'))
# DEDD5
plot(DD.results.DESeq.DEDD5$mean.discoveries$disp.zi.MDSeq, 
     DD.results.DESeq.DEDD5$mean.fdr$disp.zi.MDSeq, 
     type='l', ylim=c(0.35,0.75), xlim=c(0,500), col=col_vector[1])
lines(DD.results.DESeq.DEDD5$mean.discoveries$disp.nozi.MDSeq, 
      DD.results.DESeq.DEDD5$mean.fdr$disp.nozi.MDSeq, 
      type='l', col=col_vector[2])
lines(DD.results.DESeq.DEDD5$mean.discoveries$disp.expHM,
      DD.results.DESeq.DEDD5$mean.fdr$disp.expHM,
      type='l', col=col_vector[3])
lines(DD.results.DESeq.DEDD5$mean.discoveries$ldisp.expHM,
      DD.results.DESeq.DEDD5$mean.fdr$ldisp.expHM,
      type='l', col=col_vector[4])
lines(DD.results.DESeq.DEDD5$mean.discoveries$disp.lnHM,
      DD.results.DESeq.DEDD5$mean.fdr$disp.lnHM,
      type='l', col=col_vector[5])
lines(DD.results.DESeq.DEDD5$mean.discoveries$ldisp.lnHM,
      DD.results.DESeq.DEDD5$mean.fdr$ldisp.lnHM,
      type='l', col=col_vector[6])
legend("bottomright", bty='n', col=col_vector[1:6], lty=1, ncol=3, 
       legend=c('MDSeq +ZI', 'MDSeq -ZI', 'expHM', 'expHM log', 'lnHM', 'lnHM log'))
# DEDD10
plot(DD.results.DESeq.DEDD10$mean.discoveries$disp.zi.MDSeq, 
     DD.results.DESeq.DEDD10$mean.fdr$disp.zi.MDSeq, 
     type='l', ylim=c(0.2,0.6), xlim=c(0,500), col=col_vector[1])
lines(DD.results.DESeq.DEDD10$mean.discoveries$disp.nozi.MDSeq, 
      DD.results.DESeq.DEDD10$mean.fdr$disp.nozi.MDSeq, 
      type='l', col=col_vector[2])
lines(DD.results.DESeq.DEDD10$mean.discoveries$disp.expHM,
      DD.results.DESeq.DEDD10$mean.fdr$disp.expHM,
      type='l', col=col_vector[3])
lines(DD.results.DESeq.DEDD10$mean.discoveries$ldisp.expHM,
      DD.results.DESeq.DEDD10$mean.fdr$ldisp.expHM,
      type='l', col=col_vector[4])
lines(DD.results.DESeq.DEDD10$mean.discoveries$disp.lnHM,
      DD.results.DESeq.DEDD10$mean.fdr$disp.lnHM,
      type='l', col=col_vector[5])
lines(DD.results.DESeq.DEDD10$mean.discoveries$ldisp.lnHM,
      DD.results.DESeq.DEDD10$mean.fdr$ldisp.lnHM,
      type='l', col=col_vector[6])
legend("bottomright", bty='n', col=col_vector[1:6], lty=1, ncol=3, 
       legend=c('MDSeq +ZI', 'MDSeq -ZI', 'expHM', 'expHM log', 'lnHM', 'lnHM log'))
# DEDD20
plot(DD.results.DESeq.DEDD20$mean.discoveries$disp.zi.MDSeq, 
     DD.results.DESeq.DEDD20$mean.fdr$disp.zi.MDSeq, 
     type='l', ylim=c(0.05,0.35), xlim=c(0,500), col=col_vector[1])
lines(DD.results.DESeq.DEDD20$mean.discoveries$disp.nozi.MDSeq, 
      DD.results.DESeq.DEDD20$mean.fdr$disp.nozi.MDSeq, 
      type='l', col=col_vector[2])
lines(DD.results.DESeq.DEDD20$mean.discoveries$disp.expHM,
      DD.results.DESeq.DEDD20$mean.fdr$disp.expHM,
      type='l', col=col_vector[3])
lines(DD.results.DESeq.DEDD20$mean.discoveries$ldisp.expHM,
      DD.results.DESeq.DEDD20$mean.fdr$ldisp.expHM,
      type='l', col=col_vector[4])
lines(DD.results.DESeq.DEDD20$mean.discoveries$disp.lnHM,
      DD.results.DESeq.DEDD20$mean.fdr$disp.lnHM,
      type='l', col=col_vector[5])
lines(DD.results.DESeq.DEDD20$mean.discoveries$ldisp.lnHM,
      DD.results.DESeq.DEDD20$mean.fdr$ldisp.lnHM,
      type='l', col=col_vector[6])
legend("bottomright", bty='n', col=col_vector[1:6], lty=1, ncol=3, 
       legend=c('MDSeq +ZI', 'MDSeq -ZI', 'expHM', 'expHM log', 'lnHM', 'lnHM log'))


## Compare all differential distribution methods on same plots ####
par(mfrow=c(3,4), mar=c(2,2,1,1), mgp=c(3,0.7,0))
for (i in c('DD2', 'DD5', 'DD10', 'DD20')) {
  plot(get(paste0('DEDD.results.DESeq.',i))$mean.discoveries$expHM, 
       get(paste0('DEDD.results.DESeq.',i))$mean.fdr$expHM, 
       type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1])
  lines(get(paste0('DEDD.results.DESeq.',i))$mean.discoveries$lnHM, 
        get(paste0('DEDD.results.DESeq.',i))$mean.fdr$lnHM, 
        type='l', col=col_vector[2])
  legend("bottomright", bty='n', col=col_vector[1:2], lty=1, ncol=1, 
         legend=c('expHM', 'lnHM'))
}
for (i in c('DE2', 'DE5', 'DE10', 'DE20')) {
  plot(get(paste0('DEDD.results.DESeq.',i))$mean.discoveries$expHM, 
       get(paste0('DEDD.results.DESeq.',i))$mean.fdr$expHM, 
       type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1])
  lines(get(paste0('DEDD.results.DESeq.',i))$mean.discoveries$lnHM, 
        get(paste0('DEDD.results.DESeq.',i))$mean.fdr$lnHM, 
        type='l', col=col_vector[2])
  legend("bottomright", bty='n', col=col_vector[1:2], lty=1, ncol=1, 
         legend=c('expHM', 'lnHM'))
}
for (i in c('DEDD2', 'DEDD5', 'DEDD10', 'DEDD20')) {
  plot(get(paste0('DEDD.results.DESeq.',i))$mean.discoveries$expHM, 
       get(paste0('DEDD.results.DESeq.',i))$mean.fdr$expHM, 
       type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1])
  lines(get(paste0('DEDD.results.DESeq.',i))$mean.discoveries$lnHM, 
        get(paste0('DEDD.results.DESeq.',i))$mean.fdr$lnHM, 
        type='l', col=col_vector[2])
  legend("topleft", bty='n', col=col_vector[1:2], lty=1, ncol=1, 
         legend=c('expHM', 'lnHM'))
}


# Compare all differential distribution methods with different scales for each scenario ####
par(mfrow=c(3,4), mar=c(2,2,1,1), mgp=c(3,0.7,0))
# DD2
plot(DEDD.results.DESeq.DD2$mean.discoveries$expHM, 
     DEDD.results.DESeq.DD2$mean.fdr$expHM, 
     type='l', ylim=c(0.82,0.94), xlim=c(0,1000), col=col_vector[1])
lines(DEDD.results.DESeq.DD2$mean.discoveries$lnHM, 
      DEDD.results.DESeq.DD2$mean.fdr$lnHM, 
      type='l', col=col_vector[2])
legend("bottomright", bty='n', col=col_vector[1:2], lty=1, ncol=1, 
       legend=c('expHM', 'lnHM'))
# DD5
plot(DEDD.results.DESeq.DD5$mean.discoveries$expHM, 
     DEDD.results.DESeq.DD5$mean.fdr$expHM, 
     type='l', ylim=c(0.84,0.93), xlim=c(0,800), col=col_vector[1])
lines(DEDD.results.DESeq.DD5$mean.discoveries$lnHM, 
      DEDD.results.DESeq.DD5$mean.fdr$lnHM, 
      type='l', col=col_vector[2])
legend("bottomright", bty='n', col=col_vector[1:2], lty=1, ncol=1, 
       legend=c('expHM', 'lnHM'))
# DD10
plot(DEDD.results.DESeq.DD10$mean.discoveries$expHM, 
     DEDD.results.DESeq.DD10$mean.fdr$expHM, 
     type='l', ylim=c(0.63,0.85), xlim=c(0,800), col=col_vector[1])
lines(DEDD.results.DESeq.DD10$mean.discoveries$lnHM, 
      DEDD.results.DESeq.DD10$mean.fdr$lnHM, 
      type='l', col=col_vector[2])
legend("bottomright", bty='n', col=col_vector[1:2], lty=1, ncol=1, 
       legend=c('expHM', 'lnHM'))
# DD20
plot(DEDD.results.DESeq.DD20$mean.discoveries$expHM, 
     DEDD.results.DESeq.DD20$mean.fdr$expHM, 
     type='l', ylim=c(0,0.4), xlim=c(0,400), col=col_vector[1])
lines(DEDD.results.DESeq.DD20$mean.discoveries$lnHM, 
      DEDD.results.DESeq.DD20$mean.fdr$lnHM, 
      type='l', col=col_vector[2])
legend("bottomright", bty='n', col=col_vector[1:2], lty=1, ncol=1, 
       legend=c('expHM', 'lnHM'))
# DE2
plot(DEDD.results.DESeq.DE2$mean.discoveries$expHM, 
     DEDD.results.DESeq.DE2$mean.fdr$expHM, 
     type='l', ylim=c(0.7,0.9), xlim=c(0,1200), col=col_vector[1])
lines(DEDD.results.DESeq.DE2$mean.discoveries$lnHM, 
      DEDD.results.DESeq.DE2$mean.fdr$lnHM, 
      type='l', col=col_vector[2])
legend("bottomright", bty='n', col=col_vector[1:2], lty=1, ncol=1, 
       legend=c('expHM', 'lnHM'))
# DE5
plot(DEDD.results.DESeq.DE5$mean.discoveries$expHM, 
     DEDD.results.DESeq.DE5$mean.fdr$expHM, 
     type='l', ylim=c(0,0.06), xlim=c(0,160), col=col_vector[1])
lines(DEDD.results.DESeq.DE5$mean.discoveries$lnHM, 
      DEDD.results.DESeq.DE5$mean.fdr$lnHM, 
      type='l', col=col_vector[2])
abline(h=0.05, col='lightgrey')
legend("topleft", bty='n', col=col_vector[1:2], lty=1, ncol=1, 
       legend=c('expHM', 'lnHM'))
# DE10
plot(DEDD.results.DESeq.DE10$mean.discoveries$expHM, 
     DEDD.results.DESeq.DE10$mean.fdr$expHM, 
     type='l', ylim=c(0,0.06), xlim=c(150,500), col=col_vector[1])
lines(DEDD.results.DESeq.DE10$mean.discoveries$lnHM, 
      DEDD.results.DESeq.DE10$mean.fdr$lnHM, 
      type='l', col=col_vector[2])
abline(h=0.05, col='lightgrey')
legend("topleft", bty='n', col=col_vector[1:2], lty=1, ncol=1, 
       legend=c('expHM', 'lnHM'))
# DE20
plot(DEDD.results.DESeq.DE20$mean.discoveries$expHM, 
     DEDD.results.DESeq.DE20$mean.fdr$expHM, 
     type='l', ylim=c(0,0.06), xlim=c(400,650), col=col_vector[1])
lines(DEDD.results.DESeq.DE20$mean.discoveries$lnHM, 
      DEDD.results.DESeq.DE20$mean.fdr$lnHM, 
      type='l', col=col_vector[2])
abline(h=0.05, col='lightgrey')
legend("topleft", bty='n', col=col_vector[1:2], lty=1, ncol=1, 
       legend=c('expHM', 'lnHM'))
# DEDD2
plot(DEDD.results.DESeq.DEDD2$mean.discoveries$expHM, 
     DEDD.results.DESeq.DEDD2$mean.fdr$expHM, 
     type='l', ylim=c(0.47,0.64), xlim=c(0,1000), col=col_vector[1])
lines(DEDD.results.DESeq.DEDD2$mean.discoveries$lnHM, 
      DEDD.results.DESeq.DEDD2$mean.fdr$lnHM, 
      type='l', col=col_vector[2])
legend("bottomright", bty='n', col=col_vector[1:2], lty=1, ncol=1, 
       legend=c('expHM', 'lnHM'))
# DEDD5
plot(DEDD.results.DESeq.DEDD5$mean.discoveries$expHM, 
     DEDD.results.DESeq.DEDD5$mean.fdr$expHM, 
     type='l', ylim=c(0,0.06), xlim=c(0,500), col=col_vector[1])
lines(DEDD.results.DESeq.DEDD5$mean.discoveries$lnHM, 
      DEDD.results.DESeq.DEDD5$mean.fdr$lnHM, 
      type='l', col=col_vector[2])
abline(h=0.05, col='lightgrey')
legend("topleft", bty='n', col=col_vector[1:2], lty=1, ncol=1, 
       legend=c('expHM', 'lnHM'))
# DEDD10
plot(DEDD.results.DESeq.DEDD10$mean.discoveries$expHM, 
     DEDD.results.DESeq.DEDD10$mean.fdr$expHM, 
     type='l', ylim=c(0,0.06), xlim=c(400,1100), col=col_vector[1])
lines(DEDD.results.DESeq.DEDD10$mean.discoveries$lnHM, 
      DEDD.results.DESeq.DEDD10$mean.fdr$lnHM, 
      type='l', col=col_vector[2])
abline(h=0.05, col='lightgrey')
legend("topleft", bty='n', col=col_vector[1:2], lty=1, ncol=1, 
       legend=c('expHM', 'lnHM'))
# DEDD20
plot(DEDD.results.DESeq.DEDD20$mean.discoveries$expHM, 
     DEDD.results.DESeq.DEDD20$mean.fdr$expHM, 
     type='l', ylim=c(0,0.06), xlim=c(900,1500), col=col_vector[1])
lines(DEDD.results.DESeq.DEDD20$mean.discoveries$lnHM, 
      DEDD.results.DESeq.DEDD20$mean.fdr$lnHM, 
      type='l', col=col_vector[2])
abline(h=0.05, col='lightgrey')
legend("topleft", bty='n', col=col_vector[1:2], lty=1, ncol=1, 
       legend=c('expHM', 'lnHM'))


## Compare HMMs with DD methods on DD datasets ####
par(mfrow=c(2,4), mar=c(2,2,1,1), mgp=c(3,0.7,0))
# All on same scale
for (i in c('DD2', 'DD5', 'DD10', 'DD20')) {
  plot(get(paste0('DD.results.DESeq.',i))$mean.discoveries$disp.zi.MDSeq, 
       get(paste0('DD.results.DESeq.',i))$mean.fdr$disp.zi.MDSeq, 
       type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1])
  lines(get(paste0('DD.results.DESeq.',i))$mean.discoveries$disp.expHM, 
        get(paste0('DD.results.DESeq.',i))$mean.fdr$disp.expHM, 
        type='l', col=col_vector[2])
  lines(get(paste0('DD.results.DESeq.',i))$mean.discoveries$ldisp.expHM, 
        get(paste0('DD.results.DESeq.',i))$mean.fdr$ldisp.expHM, 
        type='l', col=col_vector[3])
  lines(get(paste0('DD.results.DESeq.',i))$mean.discoveries$disp.lnHM, 
        get(paste0('DD.results.DESeq.',i))$mean.fdr$disp.lnHM, 
        type='l', col=col_vector[4])
  lines(get(paste0('DD.results.DESeq.',i))$mean.discoveries$ldisp.lnHM, 
        get(paste0('DD.results.DESeq.',i))$mean.fdr$ldisp.lnHM, 
        type='l', col=col_vector[5])
  lines(get(paste0('DEDD.results.DESeq.',i))$mean.discoveries$expHM, 
        get(paste0('DEDD.results.DESeq.',i))$mean.fdr$expHM, 
        type='l', col=col_vector[6])
  lines(get(paste0('DEDD.results.DESeq.',i))$mean.discoveries$lnHM, 
        get(paste0('DEDD.results.DESeq.',i))$mean.fdr$lnHM, 
        type='l', col=col_vector[7])
  legend("bottomright", bty='n', col=col_vector[1:7], lty=1, ncol=2, 
         legend=c('MDSeq +ZI', 'expHM', 'expHM log', 'lnHM', 'lnHM log', 
                  'expHMM', 'lnHMM'))
}
# DD2
plot(DD.results.DESeq.DD2$mean.discoveries$disp.zi.MDSeq, 
     DD.results.DESeq.DD2$mean.fdr$disp.zi.MDSeq, 
     type='l', ylim=c(0.82,0.94), xlim=c(0,1000), col=col_vector[1])
lines(DD.results.DESeq.DD2$mean.discoveries$disp.expHM, 
      DD.results.DESeq.DD2$mean.fdr$disp.expHM, 
      type='l', col=col_vector[2])
lines(DD.results.DESeq.DD2$mean.discoveries$ldisp.expHM, 
      DD.results.DESeq.DD2$mean.fdr$ldisp.expHM, 
      type='l', col=col_vector[3])
lines(DD.results.DESeq.DD2$mean.discoveries$disp.lnHM, 
      DD.results.DESeq.DD2$mean.fdr$disp.lnHM, 
      type='l', col=col_vector[4])
lines(DD.results.DESeq.DD2$mean.discoveries$ldisp.lnHM, 
      DD.results.DESeq.DD2$mean.fdr$ldisp.lnHM, 
      type='l', col=col_vector[5])
lines(DEDD.results.DESeq.DD2$mean.discoveries$expHM, 
      DEDD.results.DESeq.DD2$mean.fdr$expHM, 
      type='l', col=col_vector[6])
lines(DEDD.results.DESeq.DD2$mean.discoveries$lnHM, 
      DEDD.results.DESeq.DD2$mean.fdr$lnHM, 
      type='l', col=col_vector[7])
legend("bottomright", bty='n', col=col_vector[1:7], lty=1, ncol=2, 
       legend=c('MDSeq +ZI', 'expHM', 'expHM log', 'lnHM', 'lnHM log', 
                'expHMM', 'lnHMM'))
# DD5
plot(DD.results.DESeq.DD5$mean.discoveries$disp.zi.MDSeq, 
     DD.results.DESeq.DD5$mean.fdr$disp.zi.MDSeq, 
     type='l', ylim=c(0.6,0.93), xlim=c(0,800), col=col_vector[1])
lines(DD.results.DESeq.DD5$mean.discoveries$disp.expHM, 
      DD.results.DESeq.DD5$mean.fdr$disp.expHM, 
      type='l', col=col_vector[2])
lines(DD.results.DESeq.DD5$mean.discoveries$ldisp.expHM, 
      DD.results.DESeq.DD5$mean.fdr$ldisp.expHM, 
      type='l', col=col_vector[3])
lines(DD.results.DESeq.DD5$mean.discoveries$disp.lnHM, 
      DD.results.DESeq.DD5$mean.fdr$disp.lnHM, 
      type='l', col=col_vector[4])
lines(DD.results.DESeq.DD5$mean.discoveries$ldisp.lnHM, 
      DD.results.DESeq.DD5$mean.fdr$ldisp.lnHM, 
      type='l', col=col_vector[5])
lines(DEDD.results.DESeq.DD5$mean.discoveries$expHM, 
      DEDD.results.DESeq.DD5$mean.fdr$expHM, 
      type='l', col=col_vector[6])
lines(DEDD.results.DESeq.DD5$mean.discoveries$lnHM, 
      DEDD.results.DESeq.DD5$mean.fdr$lnHM, 
      type='l', col=col_vector[7])
legend("bottomright", bty='n', col=col_vector[1:7], lty=1, ncol=2, 
       legend=c('MDSeq +ZI', 'expHM', 'expHM log', 'lnHM', 'lnHM log', 
                'expHMM', 'lnHMM'))
# DD10
plot(DD.results.DESeq.DD10$mean.discoveries$disp.zi.MDSeq, 
     DD.results.DESeq.DD10$mean.fdr$disp.zi.MDSeq, 
     type='l', ylim=c(0.4,0.83), xlim=c(0,800), col=col_vector[1])
lines(DD.results.DESeq.DD10$mean.discoveries$disp.expHM, 
      DD.results.DESeq.DD10$mean.fdr$disp.expHM, 
      type='l', col=col_vector[2])
lines(DD.results.DESeq.DD10$mean.discoveries$ldisp.expHM, 
      DD.results.DESeq.DD10$mean.fdr$ldisp.expHM, 
      type='l', col=col_vector[3])
lines(DD.results.DESeq.DD10$mean.discoveries$disp.lnHM, 
      DD.results.DESeq.DD10$mean.fdr$disp.lnHM, 
      type='l', col=col_vector[4])
lines(DD.results.DESeq.DD10$mean.discoveries$ldisp.lnHM, 
      DD.results.DESeq.DD10$mean.fdr$ldisp.lnHM, 
      type='l', col=col_vector[5])
lines(DEDD.results.DESeq.DD10$mean.discoveries$expHM, 
      DEDD.results.DESeq.DD10$mean.fdr$expHM, 
      type='l', col=col_vector[6])
lines(DEDD.results.DESeq.DD10$mean.discoveries$lnHM, 
      DEDD.results.DESeq.DD10$mean.fdr$lnHM, 
      type='l', col=col_vector[7])
legend("bottomright", bty='n', col=col_vector[1:7], lty=1, ncol=2, 
       legend=c('MDSeq +ZI', 'expHM', 'expHM log', 'lnHM', 'lnHM log', 
                'expHMM', 'lnHMM'))
# DD20
plot(DD.results.DESeq.DD20$mean.discoveries$disp.zi.MDSeq, 
     DD.results.DESeq.DD20$mean.fdr$disp.zi.MDSeq, 
     type='l', ylim=c(0.1,0.5), xlim=c(0,400), col=col_vector[1])
lines(DD.results.DESeq.DD20$mean.discoveries$disp.expHM, 
      DD.results.DESeq.DD20$mean.fdr$disp.expHM, 
      type='l', col=col_vector[2])
lines(DD.results.DESeq.DD20$mean.discoveries$ldisp.expHM, 
      DD.results.DESeq.DD20$mean.fdr$ldisp.expHM, 
      type='l', col=col_vector[3])
lines(DD.results.DESeq.DD20$mean.discoveries$disp.lnHM, 
      DD.results.DESeq.DD20$mean.fdr$disp.lnHM, 
      type='l', col=col_vector[4])
lines(DD.results.DESeq.DD20$mean.discoveries$ldisp.lnHM, 
      DD.results.DESeq.DD20$mean.fdr$ldisp.lnHM, 
      type='l', col=col_vector[5])
lines(DEDD.results.DESeq.DD20$mean.discoveries$expHM, 
      DEDD.results.DESeq.DD20$mean.fdr$expHM, 
      type='l', col=col_vector[6])
lines(DEDD.results.DESeq.DD20$mean.discoveries$lnHM, 
      DEDD.results.DESeq.DD20$mean.fdr$lnHM, 
      type='l', col=col_vector[7])
legend("bottomright", bty='n', col=col_vector[1:7], lty=1, ncol=2, 
       legend=c('MDSeq +ZI', 'expHM', 'expHM log', 'lnHM', 'lnHM log', 
                'expHMM', 'lnHMM'))

## Compare HMMs with DE methods on DE datasets ####
par(mfrow=c(2,4), mar=c(2,2,1,1), mgp=c(3,0.7,0))
# All on same scale
for (i in c('DE2', 'DE5', 'DE10', 'DE20')) {
  plot(get(paste0('DE.results.TMM.',i))$mean.discoveries$ql.edgeR, 
       get(paste0('DE.results.TMM.',i))$mean.fdr$ql.edgeR, 
       type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1])
  lines(get(paste0('DE.results.DESeq.',i))$mean.discoveries$if.DESeq, 
        get(paste0('DE.results.DESeq.',i))$mean.fdr$if.DESeq, 
        type='l', col=col_vector[2])
  lines(get(paste0('DE.results.DESeq.',i))$mean.discoveries$notrend.DSS, 
        get(paste0('DE.results.DESeq.',i))$mean.fdr$notrend.DSS, 
        type='l', col=col_vector[3])
  lines(get(paste0('DE.results.DESeq.',i))$mean.discoveries$mean.zi.MDSeq, 
        get(paste0('DE.results.DESeq.',i))$mean.fdr$mean.zi.MDSeq, 
        type='l', col=col_vector[4])
  lines(get(paste0('DE.results.DESeq.',i))$mean.discoveries$mean.lnHM, 
        get(paste0('DE.results.DESeq.',i))$mean.fdr$mean.lnHM, 
        type='l', col=col_vector[5])
  lines(get(paste0('DEDD.results.DESeq.',i))$mean.discoveries$expHM, 
        get(paste0('DEDD.results.DESeq.',i))$mean.fdr$expHM, 
        type='l', col=col_vector[6])
  lines(get(paste0('DEDD.results.DESeq.',i))$mean.discoveries$lnHM, 
        get(paste0('DEDD.results.DESeq.',i))$mean.fdr$lnHM, 
        type='l', col=col_vector[7])
  abline(h=0.05, col='lightgrey')
  legend("topright", bty='n', col=col_vector[1:7], lty=1, ncol=2, 
         legend=c('edgeR', 'DESeq', 'DSS', 'MDSeq', 'lnHM', 'expHMM', 
                  'lnHMM'))
}
# DE2
plot(DE.results.TMM.DE2$mean.discoveries$ql.edgeR, 
     DE.results.TMM.DE2$mean.fdr$ql.edgeR, 
     type='l', ylim=c(0.1,0.8), xlim=c(0,1200), col=col_vector[1])
lines(DE.results.DESeq.DE2$mean.discoveries$if.DESeq, 
      DE.results.DESeq.DE2$mean.fdr$if.DESeq, 
      type='l', col=col_vector[2])
lines(DE.results.DESeq.DE2$mean.discoveries$notrend.DSS, 
      DE.results.DESeq.DE2$mean.fdr$notrend.DSS, 
      type='l', col=col_vector[3])
lines(DE.results.DESeq.DE2$mean.discoveries$mean.zi.MDSeq, 
      DE.results.DESeq.DE2$mean.fdr$mean.zi.MDSeq, 
      type='l', col=col_vector[4])
lines(DE.results.DESeq.DE2$mean.discoveries$mean.lnHM, 
      DE.results.DESeq.DE2$mean.fdr$mean.lnHM, 
      type='l', col=col_vector[5])
lines(DEDD.results.DESeq.DE2$mean.discoveries$expHM, 
      DEDD.results.DESeq.DE2$mean.fdr$expHM, 
      type='l', col=col_vector[6])
lines(DEDD.results.DESeq.DE2$mean.discoveries$lnHM, 
      DEDD.results.DESeq.DE2$mean.fdr$lnHM, 
      type='l', col=col_vector[7])
legend("bottomright", bty='n', col=col_vector[1:7], lty=1, ncol=2, 
       legend=c('edgeR', 'DESeq', 'DSS', 'MDSeq', 'lnHM', 'expHMM', 
                'lnHMM'))
# DE5
plot(DE.results.TMM.DE5$mean.discoveries$ql.edgeR, 
     DE.results.TMM.DE5$mean.fdr$ql.edgeR, 
     type='l', ylim=c(0,0.06), xlim=c(0,300), col=col_vector[1])
lines(DE.results.DESeq.DE5$mean.discoveries$if.DESeq, 
      DE.results.DESeq.DE5$mean.fdr$if.DESeq, 
      type='l', col=col_vector[2])
lines(DE.results.DESeq.DE5$mean.discoveries$notrend.DSS, 
      DE.results.DESeq.DE5$mean.fdr$notrend.DSS, 
      type='l', col=col_vector[3])
lines(DE.results.DESeq.DE5$mean.discoveries$mean.zi.MDSeq, 
      DE.results.DESeq.DE5$mean.fdr$mean.zi.MDSeq, 
      type='l', col=col_vector[4])
lines(DE.results.DESeq.DE5$mean.discoveries$mean.lnHM, 
      DE.results.DESeq.DE5$mean.fdr$mean.lnHM, 
      type='l', col=col_vector[5])
lines(DEDD.results.DESeq.DE5$mean.discoveries$expHM, 
      DEDD.results.DESeq.DE5$mean.fdr$expHM, 
      type='l', col=col_vector[6])
lines(DEDD.results.DESeq.DE5$mean.discoveries$lnHM, 
      DEDD.results.DESeq.DE5$mean.fdr$lnHM, 
      type='l', col=col_vector[7])
abline(h=0.05, col='lightgrey')
legend("bottomright", bty='n', col=col_vector[1:7], lty=1, ncol=2, 
       legend=c('edgeR', 'DESeq', 'DSS', 'MDSeq', 'lnHM', 'expHMM', 
                'lnHMM'))
# DE10
plot(DE.results.TMM.DE10$mean.discoveries$ql.edgeR, 
     DE.results.TMM.DE10$mean.fdr$ql.edgeR, 
     type='l', ylim=c(0,0.06), xlim=c(150,550), col=col_vector[1])
lines(DE.results.DESeq.DE10$mean.discoveries$if.DESeq, 
      DE.results.DESeq.DE10$mean.fdr$if.DESeq, 
      type='l', col=col_vector[2])
lines(DE.results.DESeq.DE10$mean.discoveries$notrend.DSS, 
      DE.results.DESeq.DE10$mean.fdr$notrend.DSS, 
      type='l', col=col_vector[3])
lines(DE.results.DESeq.DE10$mean.discoveries$mean.zi.MDSeq, 
      DE.results.DESeq.DE10$mean.fdr$mean.zi.MDSeq, 
      type='l', col=col_vector[4])
lines(DE.results.DESeq.DE10$mean.discoveries$mean.lnHM, 
      DE.results.DESeq.DE10$mean.fdr$mean.lnHM, 
      type='l', col=col_vector[5])
lines(DEDD.results.DESeq.DE10$mean.discoveries$expHM, 
      DEDD.results.DESeq.DE10$mean.fdr$expHM, 
      type='l', col=col_vector[6])
lines(DEDD.results.DESeq.DE10$mean.discoveries$lnHM, 
      DEDD.results.DESeq.DE10$mean.fdr$lnHM, 
      type='l', col=col_vector[7])
abline(h=0.05, col='lightgrey')
legend("topleft", bty='n', col=col_vector[1:7], lty=1, ncol=2, 
       legend=c('edgeR', 'DESeq', 'DSS', 'MDSeq', 'lnHM', 'expHMM', 
                'lnHMM'))
# DE20
plot(DE.results.TMM.DE20$mean.discoveries$ql.edgeR, 
     DE.results.TMM.DE20$mean.fdr$ql.edgeR, 
     type='l', ylim=c(0,0.06), xlim=c(450,750), col=col_vector[1])
lines(DE.results.DESeq.DE20$mean.discoveries$if.DESeq, 
      DE.results.DESeq.DE20$mean.fdr$if.DESeq, 
      type='l', col=col_vector[2])
lines(DE.results.DESeq.DE20$mean.discoveries$notrend.DSS, 
      DE.results.DESeq.DE20$mean.fdr$notrend.DSS, 
      type='l', col=col_vector[3])
lines(DE.results.DESeq.DE20$mean.discoveries$mean.zi.MDSeq, 
      DE.results.DESeq.DE20$mean.fdr$mean.zi.MDSeq, 
      type='l', col=col_vector[4])
lines(DE.results.DESeq.DE20$mean.discoveries$mean.lnHM, 
      DE.results.DESeq.DE20$mean.fdr$mean.lnHM, 
      type='l', col=col_vector[5])
lines(DEDD.results.DESeq.DE20$mean.discoveries$expHM, 
      DEDD.results.DESeq.DE20$mean.fdr$expHM, 
      type='l', col=col_vector[6])
lines(DEDD.results.DESeq.DE20$mean.discoveries$lnHM, 
      DEDD.results.DESeq.DE20$mean.fdr$lnHM, 
      type='l', col=col_vector[7])
abline(h=0.05, col='lightgrey')
legend("topleft", bty='n', col=col_vector[1:7], lty=1, ncol=2, 
       legend=c('edgeR', 'DESeq', 'DSS', 'MDSeq', 'lnHM', 'expHMM', 
                'lnHMM'))


## Compare diff dist detection from HMMs and combined DE/DD methods ####
# expHMM, lnHMM, (lnHM+lnHM), (expHM+expHM), (MDSeq+MDSeq), (edgeR QL+lnHM)
combined.DEDD.results.DEDD2 <- readRDS(here('Results/Diff dist combining DE, DD predictions Nov 2019', 
                                            'diff_dist_results_DEDD2.rds'))
combined.DEDD.results.DEDD5 <- readRDS(here('Results/Diff dist combining DE, DD predictions Nov 2019', 
                                            'diff_dist_results_DEDD5.rds'))
combined.DEDD.results.DEDD10 <- readRDS(here('Results/Diff dist combining DE, DD predictions Nov 2019', 
                                             'diff_dist_results_DEDD10.rds'))
combined.DEDD.results.DEDD20 <- readRDS(here('Results/Diff dist combining DE, DD predictions Nov 2019', 
                                             'diff_dist_results_DEDD20.rds'))

par(mfrow=c(2,4), mar=c(2,2,1,1), mgp=c(3,0.7,0))
# All on same scale
for (i in c('DEDD2', 'DEDD5', 'DEDD10', 'DEDD20')) {
  plot(get(paste0('combined.DEDD.results.',i))$mean.discoveries$lnHMM, 
       get(paste0('combined.DEDD.results.',i))$mean.fdr$lnHMM, 
       type='l', ylim=c(0,1), xlim=c(0,2000), col=col_vector[1])
  lines(get(paste0('combined.DEDD.results.',i))$mean.discoveries$expHMM, 
        get(paste0('combined.DEDD.results.',i))$mean.fdr$expHMM, 
        type='l', col=col_vector[2])
  lines(get(paste0('combined.DEDD.results.',i))$mean.discoveries$lnHM_lnHM, 
        get(paste0('combined.DEDD.results.',i))$mean.fdr$lnHM_lnHM, 
        type='l', col=col_vector[3])
  lines(get(paste0('combined.DEDD.results.',i))$mean.discoveries$expHM_expHM, 
        get(paste0('combined.DEDD.results.',i))$mean.fdr$expHM_expHM, 
        type='l', col=col_vector[4])
  lines(get(paste0('combined.DEDD.results.',i))$mean.discoveries$MDSeq_MDSeq, 
        get(paste0('combined.DEDD.results.',i))$mean.fdr$MDSeq_MDSeq, 
        type='l', col=col_vector[5])
  lines(get(paste0('combined.DEDD.results.',i))$mean.discoveries$HM_edgeR, 
        get(paste0('combined.DEDD.results.',i))$mean.fdr$HM_edgeR, 
        type='l', col=col_vector[6])
  abline(h=0.05, col='lightgrey')
  legend("topright", bty='n', col=col_vector[1:7], lty=1, ncol=2, 
         legend=c('lnHMM', 'expHMM', 'lnHM, lnHM', 'expHM, expHM', 'MDSeq, MDSeq', 
                  'edgeR, HM'))
}
# DEDD2
plot(combined.DEDD.results.DEDD2$mean.discoveries$lnHMM, 
     combined.DEDD.results.DEDD2$mean.fdr$lnHMM, 
     type='l', ylim=c(0.2,0.65), xlim=c(0,1300), col=col_vector[1])
lines(combined.DEDD.results.DEDD2$mean.discoveries$expHMM, 
      combined.DEDD.results.DEDD2$mean.fdr$expHMM, 
      type='l', col=col_vector[2])
lines(combined.DEDD.results.DEDD2$mean.discoveries$lnHM_lnHM, 
      combined.DEDD.results.DEDD2$mean.fdr$lnHM_lnHM, 
      type='l', col=col_vector[3])
lines(combined.DEDD.results.DEDD2$mean.discoveries$expHM_expHM, 
      combined.DEDD.results.DEDD2$mean.fdr$expHM_expHM, 
      type='l', col=col_vector[4])
lines(combined.DEDD.results.DEDD2$mean.discoveries$MDSeq_MDSeq, 
      combined.DEDD.results.DEDD2$mean.fdr$MDSeq_MDSeq, 
      type='l', col=col_vector[5])
lines(combined.DEDD.results.DEDD2$mean.discoveries$HM_edgeR, 
      combined.DEDD.results.DEDD2$mean.fdr$HM_edgeR, 
      type='l', col=col_vector[6])
legend("bottomright", bty='n', col=col_vector[1:7], lty=1, ncol=2, 
       legend=c('lnHMM', 'expHMM', 'lnHM, lnHM', 'expHM, expHM', 'MDSeq, MDSeq', 
                'edgeR, HM'))
# DEDD5
plot(combined.DEDD.results.DEDD5$mean.discoveries$lnHMM, 
     combined.DEDD.results.DEDD5$mean.fdr$lnHMM, 
     type='l', ylim=c(0,0.06), xlim=c(0,550), col=col_vector[1])
lines(combined.DEDD.results.DEDD5$mean.discoveries$expHMM, 
      combined.DEDD.results.DEDD5$mean.fdr$expHMM, 
      type='l', col=col_vector[2])
lines(combined.DEDD.results.DEDD5$mean.discoveries$lnHM_lnHM, 
      combined.DEDD.results.DEDD5$mean.fdr$lnHM_lnHM, 
      type='l', col=col_vector[3])
lines(combined.DEDD.results.DEDD5$mean.discoveries$expHM_expHM, 
      combined.DEDD.results.DEDD5$mean.fdr$expHM_expHM, 
      type='l', col=col_vector[4])
lines(combined.DEDD.results.DEDD5$mean.discoveries$MDSeq_MDSeq, 
      combined.DEDD.results.DEDD5$mean.fdr$MDSeq_MDSeq, 
      type='l', col=col_vector[5])
lines(combined.DEDD.results.DEDD5$mean.discoveries$HM_edgeR, 
      combined.DEDD.results.DEDD5$mean.fdr$HM_edgeR, 
      type='l', col=col_vector[6])
abline(h=0.05, col='lightgrey')
legend("topleft", bty='n', col=col_vector[1:7], lty=1, ncol=2, 
       legend=c('lnHMM', 'expHMM', 'lnHM, lnHM', 'expHM, expHM', 'MDSeq, MDSeq', 
                'edgeR, HM'))
# DEDD10
plot(combined.DEDD.results.DEDD10$mean.discoveries$lnHMM, 
     combined.DEDD.results.DEDD10$mean.fdr$lnHMM, 
     type='l', ylim=c(0,0.06), xlim=c(500,1200), col=col_vector[1])
lines(combined.DEDD.results.DEDD10$mean.discoveries$expHMM, 
      combined.DEDD.results.DEDD10$mean.fdr$expHMM, 
      type='l', col=col_vector[2])
lines(combined.DEDD.results.DEDD10$mean.discoveries$lnHM_lnHM, 
      combined.DEDD.results.DEDD10$mean.fdr$lnHM_lnHM, 
      type='l', col=col_vector[3])
lines(combined.DEDD.results.DEDD10$mean.discoveries$expHM_expHM, 
      combined.DEDD.results.DEDD10$mean.fdr$expHM_expHM, 
      type='l', col=col_vector[4])
lines(combined.DEDD.results.DEDD10$mean.discoveries$MDSeq_MDSeq, 
      combined.DEDD.results.DEDD10$mean.fdr$MDSeq_MDSeq, 
      type='l', col=col_vector[5])
lines(combined.DEDD.results.DEDD10$mean.discoveries$HM_edgeR, 
      combined.DEDD.results.DEDD10$mean.fdr$HM_edgeR, 
      type='l', col=col_vector[6])
abline(h=0.05, col='lightgrey')
legend("topleft", bty='n', col=col_vector[1:7], lty=1, ncol=2, 
       legend=c('lnHMM', 'expHMM', 'lnHM, lnHM', 'expHM, expHM', 'MDSeq, MDSeq', 
                'edgeR, HM'))
# DEDD20
plot(combined.DEDD.results.DEDD20$mean.discoveries$lnHMM, 
     combined.DEDD.results.DEDD20$mean.fdr$lnHMM, 
     type='l', ylim=c(0,0.06), xlim=c(1000,1600), col=col_vector[1])
lines(combined.DEDD.results.DEDD20$mean.discoveries$expHMM, 
      combined.DEDD.results.DEDD20$mean.fdr$expHMM, 
      type='l', col=col_vector[2])
lines(combined.DEDD.results.DEDD20$mean.discoveries$lnHM_lnHM, 
      combined.DEDD.results.DEDD20$mean.fdr$lnHM_lnHM, 
      type='l', col=col_vector[3])
lines(combined.DEDD.results.DEDD20$mean.discoveries$expHM_expHM, 
      combined.DEDD.results.DEDD20$mean.fdr$expHM_expHM, 
      type='l', col=col_vector[4])
lines(combined.DEDD.results.DEDD20$mean.discoveries$MDSeq_MDSeq, 
      combined.DEDD.results.DEDD20$mean.fdr$MDSeq_MDSeq, 
      type='l', col=col_vector[5])
lines(combined.DEDD.results.DEDD20$mean.discoveries$HM_edgeR, 
      combined.DEDD.results.DEDD20$mean.fdr$HM_edgeR, 
      type='l', col=col_vector[6])
abline(h=0.05, col='lightgrey')
legend("topleft", bty='n', col=col_vector[1:7], lty=1, ncol=2, 
       legend=c('lnHMM', 'expHMM', 'lnHM, lnHM', 'expHM, expHM', 'MDSeq, MDSeq', 
                'edgeR, HM'))





