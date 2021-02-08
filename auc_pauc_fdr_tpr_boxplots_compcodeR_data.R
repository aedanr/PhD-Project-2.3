library(here)
library(RColorBrewer)

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

for (j in c('DEDD2', 'DEDD5', 'DEDD10', 'DEDD20')) {
  assign(paste0('combined_DEDD.results.',j), 
         readRDS(here('Results/Diff dist combining DE, DD predictions Nov 2019', 
                      paste0('diff_dist_results_',j,'.rds'))))
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
# To use "5.baySeq" in subset(), need to change name to not start with numeric.
names(DE.results.TMM.DE2$fdr)[which(names(DE.results.TMM.DE2$fdr) == 
                                      '5.baySeq')] <- 'thr5.baySeq'
names(DE.results.TMM.DE5$fdr)[which(names(DE.results.TMM.DE5$fdr) == 
                                      '5.baySeq')] <- 'thr5.baySeq'
names(DE.results.TMM.DE10$fdr)[which(names(DE.results.TMM.DE10$fdr) == 
                                       '5.baySeq')] <- 'thr5.baySeq'
names(DE.results.TMM.DE20$fdr)[which(names(DE.results.TMM.DE20$fdr) == 
                                       '5.baySeq')] <- 'thr5.baySeq'
names(DE.results.DESeq.DE2$fdr)[which(names(DE.results.DESeq.DE2$fdr) == 
                                        '5.baySeq')] <- 'thr5.baySeq'
names(DE.results.DESeq.DE5$fdr)[which(names(DE.results.DESeq.DE5$fdr) == 
                                        '5.baySeq')] <- 'thr5.baySeq'
names(DE.results.DESeq.DE10$fdr)[which(names(DE.results.DESeq.DE10$fdr) == 
                                         '5.baySeq')] <- 'thr5.baySeq'
names(DE.results.DESeq.DE20$fdr)[which(names(DE.results.DESeq.DE20$fdr) == 
                                         '5.baySeq')] <- 'thr5.baySeq'
names(DE.results.TMM.DEDD2$fdr)[which(names(DE.results.TMM.DEDD2$fdr) == 
                                        '5.baySeq')] <- 'thr5.baySeq'
names(DE.results.TMM.DEDD5$fdr)[which(names(DE.results.TMM.DEDD5$fdr) == 
                                        '5.baySeq')] <- 'thr5.baySeq'
names(DE.results.TMM.DEDD10$fdr)[which(names(DE.results.TMM.DEDD10$fdr) == 
                                         '5.baySeq')] <- 'thr5.baySeq'
names(DE.results.TMM.DEDD20$fdr)[which(names(DE.results.TMM.DEDD20$fdr) == 
                                         '5.baySeq')] <- 'thr5.baySeq'
names(DE.results.DESeq.DEDD2$fdr)[which(names(DE.results.DESeq.DEDD2$fdr) == 
                                          '5.baySeq')] <- 'thr5.baySeq'
names(DE.results.DESeq.DEDD5$fdr)[which(names(DE.results.DESeq.DEDD5$fdr) == 
                                          '5.baySeq')] <- 'thr5.baySeq'
names(DE.results.DESeq.DEDD10$fdr)[which(names(DE.results.DESeq.DEDD10$fdr) == 
                                           '5.baySeq')] <- 'thr5.baySeq'
names(DE.results.DESeq.DEDD20$fdr)[which(names(DE.results.DESeq.DEDD20$fdr) == 
                                           '5.baySeq')] <- 'thr5.baySeq'
names(DE.results.TMM.DE2$tpr)[which(names(DE.results.TMM.DE2$tpr) == 
                                      '5.baySeq')] <- 'thr5.baySeq'
names(DE.results.TMM.DE5$tpr)[which(names(DE.results.TMM.DE5$tpr) == 
                                      '5.baySeq')] <- 'thr5.baySeq'
names(DE.results.TMM.DE10$tpr)[which(names(DE.results.TMM.DE10$tpr) == 
                                       '5.baySeq')] <- 'thr5.baySeq'
names(DE.results.TMM.DE20$tpr)[which(names(DE.results.TMM.DE20$tpr) == 
                                       '5.baySeq')] <- 'thr5.baySeq'
names(DE.results.DESeq.DE2$tpr)[which(names(DE.results.DESeq.DE2$tpr) == 
                                        '5.baySeq')] <- 'thr5.baySeq'
names(DE.results.DESeq.DE5$tpr)[which(names(DE.results.DESeq.DE5$tpr) == 
                                        '5.baySeq')] <- 'thr5.baySeq'
names(DE.results.DESeq.DE10$tpr)[which(names(DE.results.DESeq.DE10$tpr) == 
                                         '5.baySeq')] <- 'thr5.baySeq'
names(DE.results.DESeq.DE20$tpr)[which(names(DE.results.DESeq.DE20$tpr) == 
                                         '5.baySeq')] <- 'thr5.baySeq'
names(DE.results.TMM.DEDD2$tpr)[which(names(DE.results.TMM.DEDD2$tpr) == 
                                        '5.baySeq')] <- 'thr5.baySeq'
names(DE.results.TMM.DEDD5$tpr)[which(names(DE.results.TMM.DEDD5$tpr) == 
                                        '5.baySeq')] <- 'thr5.baySeq'
names(DE.results.TMM.DEDD10$tpr)[which(names(DE.results.TMM.DEDD10$tpr) == 
                                         '5.baySeq')] <- 'thr5.baySeq'
names(DE.results.TMM.DEDD20$tpr)[which(names(DE.results.TMM.DEDD20$tpr) == 
                                         '5.baySeq')] <- 'thr5.baySeq'
names(DE.results.DESeq.DEDD2$tpr)[which(names(DE.results.DESeq.DEDD2$tpr) == 
                                          '5.baySeq')] <- 'thr5.baySeq'
names(DE.results.DESeq.DEDD5$tpr)[which(names(DE.results.DESeq.DEDD5$tpr) == 
                                          '5.baySeq')] <- 'thr5.baySeq'
names(DE.results.DESeq.DEDD10$tpr)[which(names(DE.results.DESeq.DEDD10$tpr) == 
                                           '5.baySeq')] <- 'thr5.baySeq'
names(DE.results.DESeq.DEDD20$tpr)[which(names(DE.results.DESeq.DEDD20$tpr) == 
                                           '5.baySeq')] <- 'thr5.baySeq'

dim(DD.results.TMM.DD20$auc)
# 50 6
dim(DEDD.results.TMM.DD20$auc)
# 50 2

DE.results.TMM.DE2$auc <- subset(DE.results.TMM.DE2$auc, 
                                 select = -c(trend.DSS, noif.DESeq))
DE.results.TMM.DE2$pauc <- subset(DE.results.TMM.DE2$pauc, 
                                 select = -c(trend.DSS, noif.DESeq))
DE.results.TMM.DE2$fdr <- subset(DE.results.TMM.DE2$fdr, 
                                 select = -c(fdr.trend.DSS, lfdr.trend.DSS, bh.trend.DSS, 
                                             fdr.noif.DESeq, bh.notrend.DSS, thr5.baySeq))
DE.results.TMM.DE2$tpr <- subset(DE.results.TMM.DE2$tpr, 
                                 select = -c(fdr.trend.DSS, lfdr.trend.DSS, bh.trend.DSS, 
                                             fdr.noif.DESeq, bh.notrend.DSS, thr5.baySeq))
DE.results.TMM.DE5$auc <- subset(DE.results.TMM.DE5$auc, 
                                 select = -c(trend.DSS, noif.DESeq))
DE.results.TMM.DE5$pauc <- subset(DE.results.TMM.DE5$pauc, 
                                 select = -c(trend.DSS, noif.DESeq))
DE.results.TMM.DE5$fdr <- subset(DE.results.TMM.DE5$fdr, 
                                 select = -c(fdr.trend.DSS, lfdr.trend.DSS, bh.trend.DSS, 
                                             fdr.noif.DESeq, bh.notrend.DSS, thr5.baySeq))
DE.results.TMM.DE5$tpr <- subset(DE.results.TMM.DE5$tpr, 
                                 select = -c(fdr.trend.DSS, lfdr.trend.DSS, bh.trend.DSS, 
                                             fdr.noif.DESeq, bh.notrend.DSS, thr5.baySeq))
DE.results.TMM.DE10$auc <- subset(DE.results.TMM.DE10$auc, 
                                  select = -noif.DESeq)
DE.results.TMM.DE10$pauc <- subset(DE.results.TMM.DE10$pauc, 
                                  select = -noif.DESeq)
DE.results.TMM.DE10$fdr <- subset(DE.results.TMM.DE10$fdr, 
                                  select = -c(fdr.noif.DESeq, bh.notrend.DSS, thr5.baySeq))
DE.results.TMM.DE10$tpr <- subset(DE.results.TMM.DE10$tpr, 
                                  select = -c(fdr.noif.DESeq, bh.notrend.DSS, thr5.baySeq))
DE.results.TMM.DE20$auc <- subset(DE.results.TMM.DE20$auc, 
                                  select = -noif.DESeq)
DE.results.TMM.DE20$pauc <- subset(DE.results.TMM.DE20$pauc, 
                                  select = -noif.DESeq)
DE.results.TMM.DE20$fdr <- subset(DE.results.TMM.DE20$fdr, 
                                  select = -c(fdr.noif.DESeq, bh.notrend.DSS, thr5.baySeq))
DE.results.TMM.DE20$tpr <- subset(DE.results.TMM.DE20$tpr, 
                                  select = -c(fdr.noif.DESeq, bh.notrend.DSS, thr5.baySeq))
DE.results.DESeq.DE2$auc <- subset(DE.results.DESeq.DE2$auc, 
                                   select = -noif.DESeq)
DE.results.DESeq.DE2$pauc <- subset(DE.results.DESeq.DE2$pauc, 
                                   select = -noif.DESeq)
DE.results.DESeq.DE2$fdr <- subset(DE.results.DESeq.DE2$fdr, 
                                   select = -c(fdr.noif.DESeq, bh.notrend.DSS, thr5.baySeq))
DE.results.DESeq.DE2$tpr <- subset(DE.results.DESeq.DE2$tpr, 
                                   select = -c(fdr.noif.DESeq, bh.notrend.DSS, thr5.baySeq))
DE.results.DESeq.DE5$auc <- subset(DE.results.DESeq.DE5$auc, 
                                   select = -noif.DESeq)
DE.results.DESeq.DE5$pauc <- subset(DE.results.DESeq.DE5$pauc, 
                                   select = -noif.DESeq)
DE.results.DESeq.DE5$fdr <- subset(DE.results.DESeq.DE5$fdr, 
                                   select = -c(fdr.noif.DESeq, bh.notrend.DSS, thr5.baySeq))
DE.results.DESeq.DE5$tpr <- subset(DE.results.DESeq.DE5$tpr, 
                                   select = -c(fdr.noif.DESeq, bh.notrend.DSS, thr5.baySeq))
DE.results.DESeq.DE10$auc <- subset(DE.results.DESeq.DE10$auc, 
                                    select = -noif.DESeq)
DE.results.DESeq.DE10$pauc <- subset(DE.results.DESeq.DE10$pauc, 
                                    select = -noif.DESeq)
DE.results.DESeq.DE10$fdr <- subset(DE.results.DESeq.DE10$fdr, 
                                    select = -c(fdr.noif.DESeq, bh.notrend.DSS, thr5.baySeq))
DE.results.DESeq.DE10$tpr <- subset(DE.results.DESeq.DE10$tpr, 
                                    select = -c(fdr.noif.DESeq, bh.notrend.DSS, thr5.baySeq))
DE.results.DESeq.DE20$auc <- subset(DE.results.DESeq.DE20$auc, 
                                    select = -noif.DESeq)
DE.results.DESeq.DE20$pauc <- subset(DE.results.DESeq.DE20$pauc, 
                                    select = -noif.DESeq)
DE.results.DESeq.DE20$fdr <- subset(DE.results.DESeq.DE20$fdr, 
                                    select = -c(fdr.noif.DESeq, bh.notrend.DSS, thr5.baySeq))
DE.results.DESeq.DE20$tpr <- subset(DE.results.DESeq.DE20$tpr, 
                                    select = -c(fdr.noif.DESeq, bh.notrend.DSS, thr5.baySeq))
dim(DE.results.TMM.DE20$auc)
# 50 13
dim(DEDD.results.TMM.DE20$auc)
# 50 2

dim(DD.results.TMM.DEDD20$auc)
# 50 6
DE.results.TMM.DEDD2$auc <- subset(DE.results.TMM.DEDD2$auc, 
                                   select = -c(trend.DSS, mix.ShrinkBayes, np.ShrinkBayes, 
                                               noif.DESeq))
DE.results.TMM.DEDD2$pauc <- subset(DE.results.TMM.DEDD2$pauc, 
                                   select = -c(trend.DSS, mix.ShrinkBayes, np.ShrinkBayes, 
                                               noif.DESeq))
DE.results.TMM.DEDD2$fdr <- subset(DE.results.TMM.DEDD2$fdr, 
                                   select = -c(fdr.trend.DSS, lfdr.trend.DSS, bh.trend.DSS, 
                                               lfdr.mix.ShrinkBayes, bfdr.mix.ShrinkBayes, 
                                               lfdr.np.ShrinkBayes, bfdr.np.ShrinkBayes, 
                                               fdr.noif.DESeq, bh.notrend.DSS, thr5.baySeq))
DE.results.TMM.DEDD2$tpr <- subset(DE.results.TMM.DEDD2$tpr, 
                                   select = -c(fdr.trend.DSS, lfdr.trend.DSS, bh.trend.DSS, 
                                               lfdr.mix.ShrinkBayes, bfdr.mix.ShrinkBayes, 
                                               lfdr.np.ShrinkBayes, bfdr.np.ShrinkBayes, 
                                               fdr.noif.DESeq, bh.notrend.DSS, thr5.baySeq))
DE.results.TMM.DEDD5$auc <- subset(DE.results.TMM.DEDD5$auc, 
                                   select = -noif.DESeq)
DE.results.TMM.DEDD5$pauc <- subset(DE.results.TMM.DEDD5$pauc, 
                                   select = -noif.DESeq)
DE.results.TMM.DEDD5$fdr <- subset(DE.results.TMM.DEDD5$fdr, 
                                   select = -c(fdr.noif.DESeq, bh.notrend.DSS, thr5.baySeq))
DE.results.TMM.DEDD5$tpr <- subset(DE.results.TMM.DEDD5$tpr, 
                                   select = -c(fdr.noif.DESeq, bh.notrend.DSS, thr5.baySeq))
DE.results.TMM.DEDD10$auc <- subset(DE.results.TMM.DEDD10$auc, 
                                    select = -noif.DESeq)
DE.results.TMM.DEDD10$pauc <- subset(DE.results.TMM.DEDD10$pauc, 
                                    select = -noif.DESeq)
DE.results.TMM.DEDD10$fdr <- subset(DE.results.TMM.DEDD10$fdr, 
                                    select = -c(fdr.noif.DESeq, bh.notrend.DSS, thr5.baySeq))
DE.results.TMM.DEDD10$tpr <- subset(DE.results.TMM.DEDD10$tpr, 
                                    select = -c(fdr.noif.DESeq, bh.notrend.DSS, thr5.baySeq))
DE.results.TMM.DEDD20$auc <- subset(DE.results.TMM.DEDD20$auc, 
                                    select = -noif.DESeq)
DE.results.TMM.DEDD20$pauc <- subset(DE.results.TMM.DEDD20$pauc, 
                                    select = -noif.DESeq)
DE.results.TMM.DEDD20$fdr <- subset(DE.results.TMM.DEDD20$fdr, 
                                    select = -c(fdr.noif.DESeq, bh.notrend.DSS, thr5.baySeq))
DE.results.TMM.DEDD20$tpr <- subset(DE.results.TMM.DEDD20$tpr, 
                                    select = -c(fdr.noif.DESeq, bh.notrend.DSS, thr5.baySeq))
DE.results.DESeq.DEDD2$auc <- subset(DE.results.DESeq.DEDD2$auc, 
                                     select = -noif.DESeq)
DE.results.DESeq.DEDD2$pauc <- subset(DE.results.DESeq.DEDD2$pauc, 
                                     select = -noif.DESeq)
DE.results.DESeq.DEDD2$fdr <- subset(DE.results.DESeq.DEDD2$fdr, 
                                     select = -c(fdr.noif.DESeq, bh.notrend.DSS, thr5.baySeq))
DE.results.DESeq.DEDD2$tpr <- subset(DE.results.DESeq.DEDD2$tpr, 
                                     select = -c(fdr.noif.DESeq, bh.notrend.DSS, thr5.baySeq))
DE.results.DESeq.DEDD5$auc <- subset(DE.results.DESeq.DEDD5$auc, 
                                     select = -noif.DESeq)
DE.results.DESeq.DEDD5$pauc <- subset(DE.results.DESeq.DEDD5$pauc, 
                                     select = -noif.DESeq)
DE.results.DESeq.DEDD5$fdr <- subset(DE.results.DESeq.DEDD5$fdr, 
                                     select = -c(fdr.noif.DESeq, bh.notrend.DSS, thr5.baySeq))
DE.results.DESeq.DEDD5$tpr <- subset(DE.results.DESeq.DEDD5$tpr, 
                                     select = -c(fdr.noif.DESeq, bh.notrend.DSS, thr5.baySeq))
DE.results.DESeq.DEDD10$auc <- subset(DE.results.DESeq.DEDD10$auc, 
                                      select = -noif.DESeq)
DE.results.DESeq.DEDD10$pauc <- subset(DE.results.DESeq.DEDD10$pauc, 
                                      select = -noif.DESeq)
DE.results.DESeq.DEDD10$fdr <- subset(DE.results.DESeq.DEDD10$fdr, 
                                      select = -c(fdr.noif.DESeq, bh.notrend.DSS, thr5.baySeq))
DE.results.DESeq.DEDD10$tpr <- subset(DE.results.DESeq.DEDD10$tpr, 
                                      select = -c(fdr.noif.DESeq, bh.notrend.DSS, thr5.baySeq))
DE.results.DESeq.DEDD20$auc <- subset(DE.results.DESeq.DEDD20$auc, 
                                      select = -noif.DESeq)
DE.results.DESeq.DEDD20$pauc <- subset(DE.results.DESeq.DEDD20$pauc, 
                                      select = -noif.DESeq)
DE.results.DESeq.DEDD20$fdr <- subset(DE.results.DESeq.DEDD20$fdr, 
                                      select = -c(fdr.noif.DESeq, bh.notrend.DSS, thr5.baySeq))
DE.results.DESeq.DEDD20$tpr <- subset(DE.results.DESeq.DEDD20$tpr, 
                                      select = -c(fdr.noif.DESeq, bh.notrend.DSS, thr5.baySeq))
dim(DE.results.TMM.DEDD20$auc)
# 50 13
dim(DEDD.results.TMM.DEDD20$auc)
# 50 2


# remove raw scores from FDRs and TPRs
DD.results.TMM.DD2$fdr <- 
  DD.results.TMM.DD2$fdr[,-grep('raw', names(DD.results.TMM.DD2$fdr))]
DE.results.TMM.DE2$fdr <- 
  DE.results.TMM.DE2$fdr[,-grep('raw', names(DE.results.TMM.DE2$fdr))]
DE.results.TMM.DEDD2$fdr <- 
  DE.results.TMM.DEDD2$fdr[,-grep('raw', names(DE.results.TMM.DEDD2$fdr))]
DD.results.TMM.DEDD2$fdr <- 
  DD.results.TMM.DEDD2$fdr[,-grep('raw', names(DD.results.TMM.DEDD2$fdr))]
DD.results.TMM.DD5$fdr <- 
  DD.results.TMM.DD5$fdr[,-grep('raw', names(DD.results.TMM.DD5$fdr))]
DE.results.TMM.DE5$fdr <- 
  DE.results.TMM.DE5$fdr[,-grep('raw', names(DE.results.TMM.DE5$fdr))]
DE.results.TMM.DEDD5$fdr <- 
  DE.results.TMM.DEDD5$fdr[,-grep('raw', names(DE.results.TMM.DEDD5$fdr))]
DD.results.TMM.DEDD5$fdr <- 
  DD.results.TMM.DEDD5$fdr[,-grep('raw', names(DD.results.TMM.DEDD5$fdr))]
DD.results.TMM.DD10$fdr <- 
  DD.results.TMM.DD10$fdr[,-grep('raw', names(DD.results.TMM.DD10$fdr))]
DE.results.TMM.DE10$fdr <- 
  DE.results.TMM.DE10$fdr[,-grep('raw', names(DE.results.TMM.DE10$fdr))]
DE.results.TMM.DEDD10$fdr <- 
  DE.results.TMM.DEDD10$fdr[,-grep('raw', names(DE.results.TMM.DEDD10$fdr))]
DD.results.TMM.DEDD10$fdr <- 
  DD.results.TMM.DEDD10$fdr[,-grep('raw', names(DD.results.TMM.DEDD10$fdr))]
DD.results.TMM.DD20$fdr <- 
  DD.results.TMM.DD20$fdr[,-grep('raw', names(DD.results.TMM.DD20$fdr))]
DE.results.TMM.DE20$fdr <- 
  DE.results.TMM.DE20$fdr[,-grep('raw', names(DE.results.TMM.DE20$fdr))]
DE.results.TMM.DEDD20$fdr <- 
  DE.results.TMM.DEDD20$fdr[,-grep('raw', names(DE.results.TMM.DEDD20$fdr))]
DD.results.TMM.DEDD20$fdr <- 
  DD.results.TMM.DEDD20$fdr[,-grep('raw', names(DD.results.TMM.DEDD20$fdr))]
DD.results.DESeq.DD2$fdr <- 
  DD.results.DESeq.DD2$fdr[,-grep('raw', names(DD.results.DESeq.DD2$fdr))]
DE.results.DESeq.DE2$fdr <- 
  DE.results.DESeq.DE2$fdr[,-grep('raw', names(DE.results.DESeq.DE2$fdr))]
DE.results.DESeq.DEDD2$fdr <- 
  DE.results.DESeq.DEDD2$fdr[,-grep('raw', names(DE.results.DESeq.DEDD2$fdr))]
DD.results.DESeq.DEDD2$fdr <- 
  DD.results.DESeq.DEDD2$fdr[,-grep('raw', names(DD.results.DESeq.DEDD2$fdr))]
DD.results.DESeq.DD5$fdr <- 
  DD.results.DESeq.DD5$fdr[,-grep('raw', names(DD.results.DESeq.DD5$fdr))]
DE.results.DESeq.DE5$fdr <- 
  DE.results.DESeq.DE5$fdr[,-grep('raw', names(DE.results.DESeq.DE5$fdr))]
DE.results.DESeq.DEDD5$fdr <- 
  DE.results.DESeq.DEDD5$fdr[,-grep('raw', names(DE.results.DESeq.DEDD5$fdr))]
DD.results.DESeq.DEDD5$fdr <- 
  DD.results.DESeq.DEDD5$fdr[,-grep('raw', names(DD.results.DESeq.DEDD5$fdr))]
DD.results.DESeq.DD10$fdr <- 
  DD.results.DESeq.DD10$fdr[,-grep('raw', names(DD.results.DESeq.DD10$fdr))]
DE.results.DESeq.DE10$fdr <- 
  DE.results.DESeq.DE10$fdr[,-grep('raw', names(DE.results.DESeq.DE10$fdr))]
DE.results.DESeq.DEDD10$fdr <- 
  DE.results.DESeq.DEDD10$fdr[,-grep('raw', names(DE.results.DESeq.DEDD10$fdr))]
DD.results.DESeq.DEDD10$fdr <- 
  DD.results.DESeq.DEDD10$fdr[,-grep('raw', names(DD.results.DESeq.DEDD10$fdr))]
DD.results.DESeq.DD20$fdr <- 
  DD.results.DESeq.DD20$fdr[,-grep('raw', names(DD.results.DESeq.DD20$fdr))]
DE.results.DESeq.DE20$fdr <- 
  DE.results.DESeq.DE20$fdr[,-grep('raw', names(DE.results.DESeq.DE20$fdr))]
DE.results.DESeq.DEDD20$fdr <- 
  DE.results.DESeq.DEDD20$fdr[,-grep('raw', names(DE.results.DESeq.DEDD20$fdr))]
DD.results.DESeq.DEDD20$fdr <- 
  DD.results.DESeq.DEDD20$fdr[,-grep('raw', names(DD.results.DESeq.DEDD20$fdr))]
DD.results.TMM.DD2$tpr <- 
  DD.results.TMM.DD2$tpr[,-grep('raw', names(DD.results.TMM.DD2$tpr))]
DE.results.TMM.DE2$tpr <- 
  DE.results.TMM.DE2$tpr[,-grep('raw', names(DE.results.TMM.DE2$tpr))]
DE.results.TMM.DEDD2$tpr <- 
  DE.results.TMM.DEDD2$tpr[,-grep('raw', names(DE.results.TMM.DEDD2$tpr))]
DD.results.TMM.DEDD2$tpr <- 
  DD.results.TMM.DEDD2$tpr[,-grep('raw', names(DD.results.TMM.DEDD2$tpr))]
DD.results.TMM.DD5$tpr <- 
  DD.results.TMM.DD5$tpr[,-grep('raw', names(DD.results.TMM.DD5$tpr))]
DE.results.TMM.DE5$tpr <- 
  DE.results.TMM.DE5$tpr[,-grep('raw', names(DE.results.TMM.DE5$tpr))]
DE.results.TMM.DEDD5$tpr <- 
  DE.results.TMM.DEDD5$tpr[,-grep('raw', names(DE.results.TMM.DEDD5$tpr))]
DD.results.TMM.DEDD5$tpr <- 
  DD.results.TMM.DEDD5$tpr[,-grep('raw', names(DD.results.TMM.DEDD5$tpr))]
DD.results.TMM.DD10$tpr <- 
  DD.results.TMM.DD10$tpr[,-grep('raw', names(DD.results.TMM.DD10$tpr))]
DE.results.TMM.DE10$tpr <- 
  DE.results.TMM.DE10$tpr[,-grep('raw', names(DE.results.TMM.DE10$tpr))]
DE.results.TMM.DEDD10$tpr <- 
  DE.results.TMM.DEDD10$tpr[,-grep('raw', names(DE.results.TMM.DEDD10$tpr))]
DD.results.TMM.DEDD10$tpr <- 
  DD.results.TMM.DEDD10$tpr[,-grep('raw', names(DD.results.TMM.DEDD10$tpr))]
DD.results.TMM.DD20$tpr <- 
  DD.results.TMM.DD20$tpr[,-grep('raw', names(DD.results.TMM.DD20$tpr))]
DE.results.TMM.DE20$tpr <- 
  DE.results.TMM.DE20$tpr[,-grep('raw', names(DE.results.TMM.DE20$tpr))]
DE.results.TMM.DEDD20$tpr <- 
  DE.results.TMM.DEDD20$tpr[,-grep('raw', names(DE.results.TMM.DEDD20$tpr))]
DD.results.TMM.DEDD20$tpr <- 
  DD.results.TMM.DEDD20$tpr[,-grep('raw', names(DD.results.TMM.DEDD20$tpr))]
DD.results.DESeq.DD2$tpr <- 
  DD.results.DESeq.DD2$tpr[,-grep('raw', names(DD.results.DESeq.DD2$tpr))]
DE.results.DESeq.DE2$tpr <- 
  DE.results.DESeq.DE2$tpr[,-grep('raw', names(DE.results.DESeq.DE2$tpr))]
DE.results.DESeq.DEDD2$tpr <- 
  DE.results.DESeq.DEDD2$tpr[,-grep('raw', names(DE.results.DESeq.DEDD2$tpr))]
DD.results.DESeq.DEDD2$tpr <- 
  DD.results.DESeq.DEDD2$tpr[,-grep('raw', names(DD.results.DESeq.DEDD2$tpr))]
DD.results.DESeq.DD5$tpr <- 
  DD.results.DESeq.DD5$tpr[,-grep('raw', names(DD.results.DESeq.DD5$tpr))]
DE.results.DESeq.DE5$tpr <- 
  DE.results.DESeq.DE5$tpr[,-grep('raw', names(DE.results.DESeq.DE5$tpr))]
DE.results.DESeq.DEDD5$tpr <- 
  DE.results.DESeq.DEDD5$tpr[,-grep('raw', names(DE.results.DESeq.DEDD5$tpr))]
DD.results.DESeq.DEDD5$tpr <- 
  DD.results.DESeq.DEDD5$tpr[,-grep('raw', names(DD.results.DESeq.DEDD5$tpr))]
DD.results.DESeq.DD10$tpr <- 
  DD.results.DESeq.DD10$tpr[,-grep('raw', names(DD.results.DESeq.DD10$tpr))]
DE.results.DESeq.DE10$tpr <- 
  DE.results.DESeq.DE10$tpr[,-grep('raw', names(DE.results.DESeq.DE10$tpr))]
DE.results.DESeq.DEDD10$tpr <- 
  DE.results.DESeq.DEDD10$tpr[,-grep('raw', names(DE.results.DESeq.DEDD10$tpr))]
DD.results.DESeq.DEDD10$tpr <- 
  DD.results.DESeq.DEDD10$tpr[,-grep('raw', names(DD.results.DESeq.DEDD10$tpr))]
DD.results.DESeq.DD20$tpr <- 
  DD.results.DESeq.DD20$tpr[,-grep('raw', names(DD.results.DESeq.DD20$tpr))]
DE.results.DESeq.DE20$tpr <- 
  DE.results.DESeq.DE20$tpr[,-grep('raw', names(DE.results.DESeq.DE20$tpr))]
DE.results.DESeq.DEDD20$tpr <- 
  DE.results.DESeq.DEDD20$tpr[,-grep('raw', names(DE.results.DESeq.DEDD20$tpr))]
DD.results.DESeq.DEDD20$tpr <- 
  DD.results.DESeq.DEDD20$tpr[,-grep('raw', names(DD.results.DESeq.DEDD20$tpr))]
dim(DD.results.TMM.DD20$fdr)
# 50 14
dim(DEDD.results.TMM.DD20$fdr)
# 50 6
dim(DE.results.TMM.DE20$fdr)
# 50 22
dim(DEDD.results.TMM.DE20$fdr)
# 50 6
dim(DD.results.TMM.DEDD20$fdr)
# 50 14
dim(DE.results.TMM.DEDD20$fdr)
# 50 22
dim(DEDD.results.TMM.DEDD20$fdr)
# 50 6


# DD2,5,10,20 DD: MDSeq(2), HMs(4)
# DD2,5,10,20 DEDD: HMMs(2)
# DE2,5,10,20 DE: edgeR(3), DESeq2(1), voom(1), DSS(1), baySeq(1), MDSeq(2), HMs(4)
# DE2,5,10,20 DEDD: HMMs(2)
# DEDD2,5,10,20 DD: MDSeq(2), HMs(4)
# DEDD2,5,10,20 DE: edgeR(3), DESeq2(1), voom(1), DSS(1), baySeq(1), MDSeq(2), HMs(4)
# DEDD2,5,10,20 DEDD: HMMs(2)

n <- 23
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual' & 
                                  brewer.pal.info$colorblind == T,]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# pie(rep(1,n), col=col_vector)
col_vector <- col_vector[-c(7,10,11,12,20)]

#################################
#### Differential dispersion ####
#################################

###########
### AUC ###
par(mfrow=c(2,2), mar=c(1,2,3,1), mgp=c(3,0.7,0))
boxplot(cbind(DD.results.TMM.DD2$auc, DD.results.DESeq.DD2$auc), 
        names=NA, col=col_vector[1:6], xaxt='n', ylim=c(0.45,0.9), 
        pch=20, main=paste0("AUC for differential dispersion, 2 samples per group, TMM (left)\n" , 
                            "or DESeq (right) normalisation, differences in dispersion only"))
abline(h=seq(0.48,0.56,0.02), col='lightgrey'); lines(c(6.5,6.5), c(0.48,0.56), col='darkgrey')
legend("top", fill=col_vector[1:6], ncol=3, bty='n', cex=1.2, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "exp HM", "exp HM log", "ln HM", "ln HM log"))
boxplot(cbind(DD.results.TMM.DD5$auc, DD.results.DESeq.DD5$auc), 
        names=NA, col=col_vector[1:6], xaxt='n', ylim=c(0.45,0.9), 
        pch=20, main=paste0("AUC for differential dispersion, 5 samples per group, TMM (left)\n" , 
                            "or DESeq (right) normalisation, differences in dispersion only"))
abline(h=seq(0.5,0.62,0.02), col='lightgrey'); lines(c(6.5,6.5), c(0.5,0.62), col='darkgrey')
legend("top", fill=col_vector[1:6], ncol=3, bty='n', cex=1.2, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "exp HM", "exp HM log", "ln HM", "ln HM log"))
boxplot(cbind(DD.results.TMM.DD10$auc, DD.results.DESeq.DD10$auc), 
        names=NA, col=col_vector[1:6], xaxt='n', ylim=c(0.45,0.9), 
        pch=20, main=paste0("AUC for differential dispersion, 10 samples per group, TMM (left)\n" , 
                            "or DESeq (right) normalisation, differences in dispersion only"))
abline(h=seq(0.56,0.7,0.02), col='lightgrey'); lines(c(6.5,6.5), c(0.56,0.7), col='darkgrey')
legend("top", fill=col_vector[1:6], ncol=3, bty='n', cex=1.2, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "exp HM", "exp HM log", "ln HM", "ln HM log"))
boxplot(cbind(DD.results.TMM.DD20$auc, DD.results.DESeq.DD20$auc), 
        names=NA, col=col_vector[1:6], xaxt='n', ylim=c(0.45,0.9), 
        pch=20, main=paste0("AUC for differential dispersion, 20 samples per group, TMM (left)\n" , 
                            "or DESeq (right) normalisation, differences in dispersion only"))
abline(h=seq(0.62,0.8,0.02), col='lightgrey'); lines(c(6.5,6.5), c(0.62,0.8), col='darkgrey')
legend("top", fill=col_vector[1:6], ncol=3, bty='n', cex=1.2, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "exp HM", "exp HM log", "ln HM", "ln HM log"))
colMeans(cbind(DD.results.TMM.DD2$auc, DD.results.DESeq.DD2$auc))
colMeans(cbind(DD.results.TMM.DD5$auc, DD.results.DESeq.DD5$auc))
colMeans(cbind(DD.results.TMM.DD10$auc, DD.results.DESeq.DD10$auc))
colMeans(cbind(DD.results.TMM.DD20$auc, DD.results.DESeq.DD20$auc))
# HMs better than MDSeq, improvement increases with sample size.
# MDSeq and HMs better with DESeq norm, improvement increases with sample size.

par(mfrow=c(2,2), mar=c(1,2,3,1), mgp=c(3,0.7,0))
boxplot(cbind(DD.results.TMM.DEDD2$auc, DD.results.DESeq.DEDD2$auc), 
        names=NA, col=col_vector[1:6], xaxt='n', ylim=c(0.45,0.9), 
        pch=20, main=paste0("AUC for differential dispersion, 2 samples per group, TMM (left)\n" , 
                            "or DESeq (right) normalisation, differences in mean and dispersion"))
abline(h=seq(0.48,0.56,0.02), col='lightgrey'); lines(c(6.5,6.5), c(0.48,0.56), col='darkgrey')
legend("top", fill=col_vector[1:6], ncol=3, bty='n', cex=1.2, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "exp HM", "exp HM log", "ln HM", "ln HM log"))
boxplot(cbind(DD.results.TMM.DEDD5$auc, DD.results.DESeq.DEDD5$auc), 
        names=NA, col=col_vector[1:6], xaxt='n', ylim=c(0.45,0.9), 
        pch=20, main=paste0("AUC for differential dispersion, 5 samples per group, TMM (left)\n" , 
                            "or DESeq (right) normalisation, differences in mean and dispersion"))
abline(h=seq(0.52,0.62,0.02), col='lightgrey'); lines(c(6.5,6.5), c(0.52,0.62), col='darkgrey')
legend("top", fill=col_vector[1:6], ncol=3, bty='n', cex=1.2, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "exp HM", "exp HM log", "ln HM", "ln HM log"))
boxplot(cbind(DD.results.TMM.DEDD10$auc, DD.results.DESeq.DEDD10$auc), 
        names=NA, col=col_vector[1:6], xaxt='n', ylim=c(0.45,0.9), 
        pch=20, main=paste0("AUC for differential dispersion, 10 samples per group, TMM (left)\n" , 
                            "or DESeq (right) normalisation, differences in mean and dispersion"))
abline(h=seq(0.56,0.7,0.02), col='lightgrey'); lines(c(6.5,6.5), c(0.56,0.7), col='darkgrey')
legend("top", fill=col_vector[1:6], ncol=3, bty='n', cex=1.2, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "exp HM", "exp HM log", "ln HM", "ln HM log"))
boxplot(cbind(DD.results.TMM.DEDD20$auc, DD.results.DESeq.DEDD20$auc), 
        names=NA, col=col_vector[1:6], xaxt='n', ylim=c(0.45,0.9), xlim=c(0.5,12.5), 
        pch=20, main=paste0("AUC for differential dispersion, 10 samples per group, TMM (left)\n" , 
                            "or DESeq (right) normalisation, differences in mean and dispersion"))
abline(h=seq(0.66,0.78,0.02), col='lightgrey'); lines(c(6.5,6.5), c(0.66,0.78), col='darkgrey')
legend("top", fill=col_vector[1:6], ncol=3, bty='n', cex=1.2, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "exp HM", "exp HM log", "ln HM", "ln HM log"))
colMeans(cbind(DD.results.TMM.DEDD2$auc, DD.results.DESeq.DEDD2$auc))
colMeans(cbind(DD.results.TMM.DEDD5$auc, DD.results.DESeq.DEDD5$auc))
colMeans(cbind(DD.results.TMM.DEDD10$auc, DD.results.DESeq.DEDD10$auc))
colMeans(cbind(DD.results.TMM.DEDD20$auc, DD.results.DESeq.DEDD20$auc))
# For TMM, MDSeq slightly better for some samples sizes, HMs more variable.
# For DESeq norm, very similar for smaller samples sizes, HMs clearly better for 20.
# Both better with DESeq norm, improvement increasing with sample size.


###########
### pAUC ###
par(mfrow=c(2,2), mar=c(1,2,3,1), mgp=c(3,0.7,0))
boxplot(cbind(DD.results.TMM.DD2$pauc, DD.results.DESeq.DD2$pauc), 
        names=NA, col=col_vector[1:6], xaxt='n', ylim=c(0,0.02), 
        pch=20, main=paste0("pAUC for differential dispersion, 2 samples per group, TMM (left)\n" , 
                            "or DESeq (right) normalisation, differences in dispersion only"))
abline(h=seq(0.0005,0.003,0.0005), col='lightgrey'); lines(c(6.5,6.5), c(0,0.017), col='darkgrey')
legend("top", fill=col_vector[1:6], ncol=3, bty='n', cex=1.2, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "exp HM", "exp HM log", "ln HM", "ln HM log"))
boxplot(cbind(DD.results.TMM.DD5$pauc, DD.results.DESeq.DD5$pauc), 
        names=NA, col=col_vector[1:6], xaxt='n', ylim=c(0,0.02), 
        pch=20, main=paste0("pAUC for differential dispersion, 5 samples per group, TMM (left)\n" , 
                            "or DESeq (right) normalisation, differences in dispersion only"))
abline(h=seq(0.001,0.005,0.001), col='lightgrey'); lines(c(6.5,6.5), c(0,0.017), col='darkgrey')
legend("top", fill=col_vector[1:6], ncol=3, bty='n', cex=1.2, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "exp HM", "exp HM log", "ln HM", "ln HM log"))
boxplot(cbind(DD.results.TMM.DD10$pauc, DD.results.DESeq.DD10$pauc), 
        names=NA, col=col_vector[1:6], xaxt='n', ylim=c(0,0.02), 
        pch=20, main=paste0("pAUC for differential dispersion, 10 samples per group, TMM (left)\n" , 
                            "or DESeq (right) normalisation, differences in dispersion only"))
abline(h=seq(0.002,0.009,0.001), col='lightgrey'); lines(c(6.5,6.5), c(0,0.017), col='darkgrey')
legend("top", fill=col_vector[1:6], ncol=3, bty='n', cex=1.2, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "exp HM", "exp HM log", "ln HM", "ln HM log"))
boxplot(cbind(DD.results.TMM.DD20$pauc, DD.results.DESeq.DD20$pauc), 
        names=NA, col=col_vector[1:6], xaxt='n', ylim=c(0,0.02), 
        pch=20, main=paste0("pAUC for differential dispersion, 20 samples per group, TMM (left)\n" , 
                            "or DESeq (right) normalisation, differences in dispersion only"))
abline(h=seq(0.005,0.017,0.001), col='lightgrey'); lines(c(6.5,6.5), c(0,0.017), col='darkgrey')
legend("top", fill=col_vector[1:6], ncol=3, bty='n', cex=1.2, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "exp HM", "exp HM log", "ln HM", "ln HM log"))
colMeans(cbind(DD.results.TMM.DD2$pauc, DD.results.DESeq.DD2$pauc))
colMeans(cbind(DD.results.TMM.DD5$pauc, DD.results.DESeq.DD5$pauc))
colMeans(cbind(DD.results.TMM.DD10$pauc, DD.results.DESeq.DD10$pauc))
colMeans(cbind(DD.results.TMM.DD20$pauc, DD.results.DESeq.DD20$pauc))

par(mfrow=c(2,2), mar=c(1,2,3,1), mgp=c(3,0.7,0))
boxplot(cbind(DD.results.TMM.DEDD2$pauc, DD.results.DESeq.DEDD2$pauc), 
        names=NA, col=col_vector[1:6], xaxt='n', ylim=c(0,0.02), 
        pch=20, main=paste0("pAUC for differential dispersion, 2 samples per group, TMM (left)\n" , 
                            "or DESeq (right) normalisation, differences in mean and dispersion"))
abline(h=seq(0.0005,0.003,0.0005), col='lightgrey'); lines(c(6.5,6.5), c(0,0.017), col='darkgrey')
legend("top", fill=col_vector[1:6], ncol=3, bty='n', cex=1.2, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "exp HM", "exp HM log", "ln HM", "ln HM log"))
boxplot(cbind(DD.results.TMM.DEDD5$pauc, DD.results.DESeq.DEDD5$pauc), 
        names=NA, col=col_vector[1:6], xaxt='n', ylim=c(0,0.02), 
        pch=20, main=paste0("pAUC for differential dispersion, 5 samples per group, TMM (left)\n" , 
                            "or DESeq (right) normalisation, differences in mean and dispersion"))
abline(h=seq(0.001,0.005,0.001), col='lightgrey'); lines(c(6.5,6.5), c(0,0.017), col='darkgrey')
legend("top", fill=col_vector[1:6], ncol=3, bty='n', cex=1.2, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "exp HM", "exp HM log", "ln HM", "ln HM log"))
boxplot(cbind(DD.results.TMM.DEDD10$pauc, DD.results.DESeq.DEDD10$pauc), 
        names=NA, col=col_vector[1:6], xaxt='n', ylim=c(0,0.02), 
        pch=20, main=paste0("pAUC for differential dispersion, 10 samples per group, TMM (left)\n" , 
                            "or DESeq (right) normalisation, differences in mean and dispersion"))
abline(h=seq(0.002,0.009,0.001), col='lightgrey'); lines(c(6.5,6.5), c(0,0.017), col='darkgrey')
legend("top", fill=col_vector[1:6], ncol=3, bty='n', cex=1.2, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "exp HM", "exp HM log", "ln HM", "ln HM log"))
boxplot(cbind(DD.results.TMM.DEDD20$pauc, DD.results.DESeq.DEDD20$pauc), 
        names=NA, col=col_vector[1:6], xaxt='n', ylim=c(0,0.02), xlim=c(0.5,12.5), 
        pch=20, main=paste0("pAUC for differential dispersion, 10 samples per group, TMM (left)\n" , 
                            "or DESeq (right) normalisation, differences in mean and dispersion"))
abline(h=seq(0.005,0.017,0.001), col='lightgrey'); lines(c(6.5,6.5), c(0,0.017), col='darkgrey')
legend("top", fill=col_vector[1:6], ncol=3, bty='n', cex=1.2, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "exp HM", "exp HM log", "ln HM", "ln HM log"))
colMeans(cbind(DD.results.TMM.DEDD2$pauc, DD.results.DESeq.DEDD2$pauc))
colMeans(cbind(DD.results.TMM.DEDD5$pauc, DD.results.DESeq.DEDD5$pauc))
colMeans(cbind(DD.results.TMM.DEDD10$pauc, DD.results.DESeq.DEDD10$pauc))
colMeans(cbind(DD.results.TMM.DEDD20$pauc, DD.results.DESeq.DEDD20$pauc))


###########
### FDR ###
par(mfrow=c(2,2), mar=c(1,2,3,1), mgp=c(3,0.7,0))
boxplot(cbind(DD.results.TMM.DD2$fdr, DD.results.DESeq.DD2$fdr), 
        ylim=c(0,1.2), names=NA, col=col_vector[1:14], xaxt='n', 
        pch=20, main=paste0("FDR for differential dispersion, 2 samples per group, TMM (left)\n" , 
                            "or DESeq (right) normalisation, differences in dispersion only"))
abline(h=0.05, col='lightgrey'); lines(c(14.5,14.5), c(0,1), col='darkgrey')
legend("top", fill=col_vector[1:14], ncol=4, bty='n', cex=0.9, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "exp HM (BH)", "exp HM (BY)", "exp HM (q)", 
                "exp HM log (BH)", "exp HM log (BY)", "exp HM log (q)", "ln HM (BH)", "ln HM (BY)", 
                "ln HM (q)", "ln HM log (BH)", "ln HM log (BY)", "ln HM log (q)"))
boxplot(cbind(DD.results.TMM.DD5$fdr, DD.results.DESeq.DD5$fdr), 
        ylim=c(0,1.2), names=NA, col=col_vector[1:14], xaxt='n', 
        pch=20, main=paste0("FDR for differential dispersion, 5 samples per group, TMM (left)\n" , 
                            "or DESeq (right) normalisation, differences in dispersion only"))
abline(h=0.05, col='lightgrey'); lines(c(14.5,14.5), c(0,1), col='darkgrey')
legend("top", fill=col_vector[1:14], ncol=4, bty='n', cex=0.9, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "exp HM (BH)", "exp HM (BY)", "exp HM (q)", 
                "exp HM log (BH)", "exp HM log (BY)", "exp HM log (q)", "ln HM (BH)", "ln HM (BY)", 
                "ln HM (q)", "ln HM log (BH)", "ln HM log (BY)", "ln HM log (q)"))
boxplot(cbind(DD.results.TMM.DD10$fdr, DD.results.DESeq.DD10$fdr), 
        ylim=c(0,1.2), names=NA, col=col_vector[1:14], xaxt='n', 
        pch=20, main=paste0("FDR for differential dispersion, 10 samples per group, TMM (left)\n" , 
                            "or DESeq (right) normalisation, differences in dispersion only"))
abline(h=0.05, col='lightgrey'); lines(c(14.5,14.5), c(0,1), col='darkgrey')
legend("top", fill=col_vector[1:14], ncol=4, bty='n', cex=0.9, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "exp HM (BH)", "exp HM (BY)", "exp HM (q)", 
                "exp HM log (BH)", "exp HM log (BY)", "exp HM log (q)", "ln HM (BH)", "ln HM (BY)", 
                "ln HM (q)", "ln HM log (BH)", "ln HM log (BY)", "ln HM log (q)"))
boxplot(cbind(DD.results.TMM.DD20$fdr, DD.results.DESeq.DD20$fdr), 
        ylim=c(0,1.2), names=NA, col=col_vector[1:14], xaxt='n', 
        pch=20, main=paste0("FDR for differential dispersion, 20 samples per group, TMM (left)\n" , 
                            "or DESeq (right) normalisation, differences in dispersion only"))
abline(h=0.05, col='lightgrey'); lines(c(14.5,14.5), c(0,1), col='darkgrey')
legend("top", fill=col_vector[1:14], ncol=4, bty='n', cex=0.9, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "exp HM (BH)", "exp HM (BY)", "exp HM (q)", 
                "exp HM log (BH)", "exp HM log (BY)", "exp HM log (q)", "ln HM (BH)", "ln HM (BY)", 
                "ln HM (q)", "ln HM log (BH)", "ln HM log (BY)", "ln HM log (q)"))
colMeans(cbind(DD.results.TMM.DD2$fdr, DD.results.DESeq.DD2$fdr), na.rm=T)
colMeans(cbind(DD.results.TMM.DD5$fdr, DD.results.DESeq.DD5$fdr), na.rm=T)
colMeans(cbind(DD.results.TMM.DD10$fdr, DD.results.DESeq.DD10$fdr), na.rm=T)
colMeans(cbind(DD.results.TMM.DD20$fdr, DD.results.DESeq.DD20$fdr), na.rm=T)
# HMs no calls for 2, 5, or 10 samples.
# High FDRs for 20 samples, but better and less variable than MDSeq.
# DESeq norm better for MDSeq and HMs.
# No calls for HMs with BY. BH far worse than q with TMM for expHM, but similar for lnHM and with 
# DESeq norm.

par(mfrow=c(2,2), mar=c(1,2,3,1), mgp=c(3,0.7,0))
boxplot(cbind(DD.results.TMM.DEDD2$fdr, DD.results.DESeq.DEDD2$fdr), 
        ylim=c(0,1.2), names=NA, col=col_vector[1:14], xaxt='n', 
        pch=20, main=paste0("FDR for differential dispersion, 2 samples per group, TMM (left)\n" , 
                            "or DESeq (right) normalisation, differences in mean and dispersion"))
abline(h=0.05, col='lightgrey'); lines(c(14.5,14.5), c(0,1), col='darkgrey')
legend("top", fill=col_vector[1:14], ncol=4, bty='n', cex=0.9, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "exp HM (BH)", "exp HM (BY)", "exp HM (q)", 
                "exp HM log (BH)", "exp HM log (BY)", "exp HM log (q)", "ln HM (BH)", "ln HM (BY)", 
                "ln HM (q)", "ln HM log (BH)", "ln HM log (BY)", "ln HM log (q)"))
boxplot(cbind(DD.results.TMM.DEDD5$fdr, DD.results.DESeq.DEDD5$fdr), 
        ylim=c(0,1.2), names=NA, col=col_vector[1:14], xaxt='n', 
        pch=20, main=paste0("FDR for differential dispersion, 5 samples per group, TMM (left)\n" , 
                            "or DESeq (right) normalisation, differences in mean and dispersion"))
abline(h=0.05, col='lightgrey'); lines(c(14.5,14.5), c(0,1), col='darkgrey')
legend("top", fill=col_vector[1:14], ncol=4, bty='n', cex=0.9, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "exp HM (BH)", "exp HM (BY)", "exp HM (q)", 
                "exp HM log (BH)", "exp HM log (BY)", "exp HM log (q)", "ln HM (BH)", "ln HM (BY)", 
                "ln HM (q)", "ln HM log (BH)", "ln HM log (BY)", "ln HM log (q)"))
boxplot(cbind(DD.results.TMM.DEDD10$fdr, DD.results.DESeq.DEDD10$fdr), 
        ylim=c(0,1.2), names=NA, col=col_vector[1:14], xaxt='n', 
        pch=20, main=paste0("FDR for differential dispersion, 10 samples per group, TMM (left)\n" , 
                            "or DESeq (right) normalisation, differences in mean and dispersion"))
abline(h=0.05, col='lightgrey'); lines(c(14.5,14.5), c(0,1), col='darkgrey')
legend("top", fill=col_vector[1:14], ncol=4, bty='n', cex=0.9, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "exp HM (BH)", "exp HM (BY)", "exp HM (q)", 
                "exp HM log (BH)", "exp HM log (BY)", "exp HM log (q)", "ln HM (BH)", "ln HM (BY)", 
                "ln HM (q)", "ln HM log (BH)", "ln HM log (BY)", "ln HM log (q)"))
boxplot(cbind(DD.results.TMM.DEDD20$fdr, DD.results.DESeq.DEDD20$fdr), 
        ylim=c(0,1.2), xlim=c(0.5,28.5), names=NA, col=col_vector[1:14], xaxt='n', 
        pch=20, main=paste0("FDR for differential dispersion, 20 samples per group, TMM (left)\n" , 
                            "or DESeq (right) normalisation, differences in mean and dispersion"))
abline(h=0.05, col='lightgrey'); lines(c(14.5,14.5), c(0,1), col='darkgrey')
legend("top", fill=col_vector[1:14], ncol=4, bty='n', cex=0.9, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "exp HM (BH)", "exp HM (BY)", "exp HM (q)", 
                "exp HM log (BH)", "exp HM log (BY)", "exp HM log (q)", "ln HM (BH)", "ln HM (BY)", 
                "ln HM (q)", "ln HM log (BH)", "ln HM log (BY)", "ln HM log (q)"))
colMeans(cbind(DD.results.TMM.DEDD2$fdr, DD.results.DESeq.DEDD2$fdr), na.rm=T)
colMeans(cbind(DD.results.TMM.DEDD5$fdr, DD.results.DESeq.DEDD5$fdr), na.rm=T)
colMeans(cbind(DD.results.TMM.DEDD10$fdr, DD.results.DESeq.DEDD10$fdr), na.rm=T)
colMeans(cbind(DD.results.TMM.DEDD20$fdr, DD.results.DESeq.DEDD20$fdr), na.rm=T)
# HMs nearly no calls for 2, 5, or 10 samples. For 20 samples, better than MDSeq, and better than 
# DD only data.
# DESeq norm better for MDSeq and HMs; means not much betterfor lnHM, but variance lower.
# No calls for HMs with BY. BH slightly better than q. Close to but still above 0.05 for lnHM.


###########
### TPR ###
par(mfrow=c(2,2), mar=c(1,2,3,1), mgp=c(3,0.7,0))
boxplot(cbind(DD.results.TMM.DD2$tpr, DD.results.DESeq.DD2$tpr), ylim=c(0,1.2), 
        names=NA, col=col_vector[1:14], xaxt='n', 
        pch=20, main=paste0("TPR for differential dispersion, 2 samples per group, TMM (left)\n" , 
                            "or DESeq (right) normalisation, differences in dispersion only"))
legend("top", fill=col_vector[1:14], ncol=4, bty='n', cex=0.9, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "exp HM (BH)", "exp HM (BY)", "exp HM (q)", 
                "exp HM log (BH)", "exp HM log (BY)", "exp HM log (q)", "ln HM (BH)", "ln HM (BY)", 
                "ln HM (q)", "ln HM log (BH)", "ln HM log (BY)", "ln HM log (q)"))
lines(c(14.5,14.5), c(0,1), col='darkgrey')
boxplot(cbind(DD.results.TMM.DD5$tpr, DD.results.DESeq.DD5$tpr), ylim=c(0,1.2), 
        names=NA, col=col_vector[1:14], xaxt='n', 
        pch=20, main=paste0("TPR for differential dispersion, 5 samples per group, TMM (left)\n" , 
                            "or DESeq (right) normalisation, differences in dispersion only"))
legend("top", fill=col_vector[1:14], ncol=4, bty='n', cex=0.9, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "exp HM (BH)", "exp HM (BY)", "exp HM (q)", 
                "exp HM log (BH)", "exp HM log (BY)", "exp HM log (q)", "ln HM (BH)", "ln HM (BY)", 
                "ln HM (q)", "ln HM log (BH)", "ln HM log (BY)", "ln HM log (q)"))
lines(c(14.5,14.5), c(0,1), col='darkgrey')
boxplot(cbind(DD.results.TMM.DD10$tpr, DD.results.DESeq.DD10$tpr), ylim=c(0,1.2), 
        names=NA, col=col_vector[1:14], xaxt='n', 
        pch=20, main=paste0("TPR for differential dispersion, 10 samples per group, TMM (left)\n" , 
                            "or DESeq (right) normalisation, differences in dispersion only"))
legend("top", fill=col_vector[1:14], ncol=4, bty='n', cex=0.9, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "exp HM (BH)", "exp HM (BY)", "exp HM (q)", 
                "exp HM log (BH)", "exp HM log (BY)", "exp HM log (q)", "ln HM (BH)", "ln HM (BY)", 
                "ln HM (q)", "ln HM log (BH)", "ln HM log (BY)", "ln HM log (q)"))
lines(c(14.5,14.5), c(0,1), col='darkgrey')
boxplot(cbind(DD.results.TMM.DD20$tpr, DD.results.DESeq.DD20$tpr), ylim=c(0,1.2), 
        names=NA, col=col_vector[1:14], xaxt='n', 
        pch=20, main=paste0("TPR for differential dispersion, 20 samples per group, TMM (left)\n" , 
                            "or DESeq (right) normalisation, differences in dispersion only"))
legend("top", fill=col_vector[1:14], ncol=4, bty='n', cex=0.9, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "exp HM (BH)", "exp HM (BY)", "exp HM (q)", 
                "exp HM log (BH)", "exp HM log (BY)", "exp HM log (q)", "ln HM (BH)", "ln HM (BY)", 
                "ln HM (q)", "ln HM log (BH)", "ln HM log (BY)", "ln HM log (q)"))
lines(c(14.5,14.5), c(0,1), col='darkgrey')
colMeans(cbind(DD.results.TMM.DD2$tpr, DD.results.DESeq.DD2$tpr), na.rm=T)
colMeans(cbind(DD.results.TMM.DD5$tpr, DD.results.DESeq.DD5$tpr), na.rm=T)
colMeans(cbind(DD.results.TMM.DD10$tpr, DD.results.DESeq.DD10$tpr), na.rm=T)
colMeans(cbind(DD.results.TMM.DD20$tpr, DD.results.DESeq.DD20$tpr), na.rm=T)
# HMs no calls for 2, 5, or 10 samples. Very low power for 20 samples, but better than MDSeq.
# For MDSeq with TMM, power decreases with increasing sample size; inconsistent for DESeq norm.
# DESeq norm always better than TMM for HMs.

par(mfrow=c(2,2), mar=c(1,2,3,1), mgp=c(3,0.7,0))
boxplot(cbind(DD.results.TMM.DEDD2$tpr, DD.results.DESeq.DEDD2$tpr), 
        ylim=c(0,1.2), names=NA, col=col_vector[1:14], xaxt='n', 
        pch=20, main=paste0("TPR for differential dispersion, 2 samples per group, TMM (left)\n" , 
                            "or DESeq (right) normalisation, differences in mean and dispersion"))
legend("top", fill=col_vector[1:14], ncol=4, bty='n', cex=0.9, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "exp HM (BH)", "exp HM (BY)", "exp HM (q)", 
                "exp HM log (BH)", "exp HM log (BY)", "exp HM log (q)", "ln HM (BH)", "ln HM (BY)", 
                "ln HM (q)", "ln HM log (BH)", "ln HM log (BY)", "ln HM log (q)"))
abline(h=seq(0,0.08,0.02), col='lightgrey'); lines(c(14.5,14.5), c(0,1), col='darkgrey')
boxplot(cbind(DD.results.TMM.DEDD5$tpr, DD.results.DESeq.DEDD5$tpr), 
        ylim=c(0,1.2), names=NA, col=col_vector[1:14], xaxt='n', 
        pch=20, main=paste0("TPR for differential dispersion, 5 samples per group, TMM (left)\n" , 
                            "or DESeq (right) normalisation, differences in mean and dispersion"))
legend("top", fill=col_vector[1:14], ncol=4, bty='n', cex=0.9, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "exp HM (BH)", "exp HM (BY)", "exp HM (q)", 
                "exp HM log (BH)", "exp HM log (BY)", "exp HM log (q)", "ln HM (BH)", "ln HM (BY)", 
                "ln HM (q)", "ln HM log (BH)", "ln HM log (BY)", "ln HM log (q)"))
abline(h=seq(0,0.04,0.02), col='lightgrey'); lines(c(14.5,14.5), c(0,1), col='darkgrey')
boxplot(cbind(DD.results.TMM.DEDD10$tpr, DD.results.DESeq.DEDD10$tpr), 
        ylim=c(0,1.2), names=NA, col=col_vector[1:14], xaxt='n', 
        pch=20, main=paste0("TPR for differential dispersion, 10 samples per group, TMM (left)\n" , 
                            "or DESeq (right) normalisation, differences in mean and dispersion"))
legend("top", fill=col_vector[1:14], ncol=4, bty='n', cex=0.9, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "exp HM (BH)", "exp HM (BY)", "exp HM (q)", 
                "exp HM log (BH)", "exp HM log (BY)", "exp HM log (q)", "ln HM (BH)", "ln HM (BY)", 
                "ln HM (q)", "ln HM log (BH)", "ln HM log (BY)", "ln HM log (q)"))
abline(h=seq(0,0.04,0.02), col='lightgrey'); lines(c(14.5,14.5), c(0,1), col='darkgrey')
boxplot(cbind(DD.results.TMM.DEDD20$tpr, DD.results.DESeq.DEDD20$tpr), 
        ylim=c(0,1.2), xlim=c(0.5,28.5), names=NA, col=col_vector[1:14], xaxt='n', 
        pch=20, main=paste0("TPR for differential dispersion, 20 samples per group, TMM (left)\n" , 
                            "or DESeq (right) normalisation, differences in mean and dispersion"))
legend("top", fill=col_vector[1:14], ncol=4, bty='n', cex=0.9, 
       legend=c("MDSeq ZI", "MDSeq no ZI", "exp HM (BH)", "exp HM (BY)", "exp HM (q)", 
                "exp HM log (BH)", "exp HM log (BY)", "exp HM log (q)", "ln HM (BH)", "ln HM (BY)", 
                "ln HM (q)", "ln HM log (BH)", "ln HM log (BY)", "ln HM log (q)"))
abline(h=seq(0,0.12,0.02), col='lightgrey'); lines(c(14.5,14.5), c(0,1), col='darkgrey')
colMeans(cbind(DD.results.TMM.DEDD2$tpr, DD.results.DESeq.DEDD2$tpr), na.rm=T)
colMeans(cbind(DD.results.TMM.DEDD5$tpr, DD.results.DESeq.DEDD5$tpr), na.rm=T)
colMeans(cbind(DD.results.TMM.DEDD10$tpr, DD.results.DESeq.DEDD10$tpr), na.rm=T)
colMeans(cbind(DD.results.TMM.DEDD20$tpr, DD.results.DESeq.DEDD20$tpr), na.rm=T)
# HMs nearly no calls for 2, 5, or 10 samples. Very low power for 20 samples, and lower than MDSeq.
# Inconsistent relationship between power and sample size for MDSeq.
# DESeq norm always better than TMM for HMs.


#################################
#### Differential expression ####
#################################

###########
### AUC ###
par(mfrow=c(2,2), mar=c(1,2,3,1), mgp=c(3,0.7,0))
boxplot(cbind(DE.results.TMM.DE2$auc, DE.results.DESeq.DE2$auc), 
        names=NA, col=col_vector[1:13], xaxt='n', ylim=c(0.6,1.07), 
        pch=20, main=paste0("AUC for differential expression, 2 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in mean only"))
abline(h=seq(0.6,0.85,0.05), col='lightgrey'); lines(c(13.5,13.5), c(0.6,0.85), col='darkgrey')
legend("top", fill=col_vector[1:13], ncol=5, bty='n', cex=0.9, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq", "voom", "DSS", "baySeq", 
                "MDSeq ZI", "MDSeq no ZI", "exp HM", "exp HM log", "ln HM", "ln HM log"))
boxplot(cbind(DE.results.TMM.DE5$auc, DE.results.DESeq.DE5$auc), 
        names=NA, col=col_vector[1:13], xaxt='n', ylim=c(0.6,1.07), 
        pch=20, main=paste0("AUC for differential expression, 5 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in mean only"))
abline(h=seq(0.75,0.95,0.05), col='lightgrey'); lines(c(13.5,13.5), c(0.75,0.95), col='darkgrey')
legend("top", fill=col_vector[1:13], ncol=5, bty='n', cex=0.9, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq", "voom", "DSS", "baySeq", 
                "MDSeq ZI", "MDSeq no ZI", "exp HM", "exp HM log", "ln HM", "ln HM log"))
boxplot(cbind(DE.results.TMM.DE10$auc, DE.results.DESeq.DE10$auc), 
        names=NA, col=col_vector[1:13], xaxt='n', ylim=c(0.6,1.07), 
        pch=20, main=paste0("AUC for differential expression, 5 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in mean only"))
abline(h=seq(0.82,0.98,0.02), col='lightgrey'); lines(c(13.5,13.5), c(0.82,0.98), col='darkgrey')
legend("top", fill=col_vector[1:13], ncol=5, bty='n', cex=0.9, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq", "voom", "DSS", "baySeq", 
                "MDSeq ZI", "MDSeq no ZI", "exp HM", "exp HM log", "ln HM", "ln HM log"))
boxplot(cbind(DE.results.TMM.DE20$auc, DE.results.DESeq.DE20$auc), 
        names=NA, col=col_vector[1:13], xaxt='n', ylim=c(0.6,1.07), xlim=c(0.5,26.5), 
        pch=20, main=paste0("AUC for differential expression, 20 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in mean only"))
abline(h=seq(0.88,1,0.02), col='lightgrey'); lines(c(13.5,13.5), c(0.88,1), col='darkgrey')
legend("top", fill=col_vector[1:13], ncol=5, bty='n', cex=0.9, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq", "voom", "DSS", "baySeq", 
                "MDSeq ZI", "MDSeq no ZI", "exp HM", "exp HM log", "ln HM", "ln HM log"))
colMeans(cbind(DE.results.TMM.DE2$auc, DE.results.DESeq.DE2$auc))
colMeans(cbind(DE.results.TMM.DE5$auc, DE.results.DESeq.DE5$auc))
colMeans(cbind(DE.results.TMM.DE10$auc, DE.results.DESeq.DE10$auc))
colMeans(cbind(DE.results.TMM.DE20$auc, DE.results.DESeq.DE20$auc))
# lnHM with DESeq norm equal to edgeR with TMM, similar to DSS and DESeq with DESeq norm, better 
# than all others.

par(mfrow=c(2,2), mar=c(1,2,3,1), mgp=c(3,0.7,0))
boxplot(cbind(DE.results.TMM.DEDD2$auc, DE.results.DESeq.DEDD2$auc), 
        names=NA, col=col_vector[1:13], xaxt='n', ylim=c(0.6,1.07), 
        pch=20, main=paste0("AUC for differential expression, 2 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in mean and dispersion"))
abline(h=seq(0.6,0.85,0.05), col='lightgrey'); lines(c(13.5,13.5), c(0.6,0.85), col='darkgrey')
legend("top", fill=col_vector[1:13], ncol=5, bty='n', cex=0.9, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq", "voom", "DSS", "baySeq", 
                "MDSeq ZI", "MDSeq no ZI", "exp HM", "exp HM log", "ln HM", "ln HM log"))
boxplot(cbind(DE.results.TMM.DEDD5$auc, DE.results.DESeq.DEDD5$auc), 
        names=NA, col=col_vector[1:13], xaxt='n', ylim=c(0.6,1.07), 
        pch=20, main=paste0("AUC for differential expression, 5 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in mean and dispersion"))
abline(h=seq(0.7,0.95,0.05), col='lightgrey'); lines(c(13.5,13.5), c(0.7,0.95), col='darkgrey')
legend("top", fill=col_vector[1:13], ncol=5, bty='n', cex=0.9, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq", "voom", "DSS", "baySeq", 
                "MDSeq ZI", "MDSeq no ZI", "exp HM", "exp HM log", "ln HM", "ln HM log"))
boxplot(cbind(DE.results.TMM.DEDD10$auc, DE.results.DESeq.DEDD10$auc), 
        names=NA, col=col_vector[1:13], xaxt='n', ylim=c(0.6,1.07), 
        pch=20, main=paste0("AUC for differential expression, 10 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in mean and dispersion"))
abline(h=seq(0.84,0.96,0.02), col='lightgrey'); lines(c(13.5,13.5), c(0.84,0.96), col='darkgrey')
legend("top", fill=col_vector[1:13], ncol=5, bty='n', cex=0.9, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq", "voom", "DSS", "baySeq", 
                "MDSeq ZI", "MDSeq no ZI", "exp HM", "exp HM log", "ln HM", "ln HM log"))
boxplot(cbind(DE.results.TMM.DEDD20$auc, DE.results.DESeq.DEDD20$auc), 
        names=NA, col=col_vector[1:13], xaxt='n', ylim=c(0.6,1.07), xlim=c(0.5,26.5), 
        pch=20, main=paste0("AUC for differential expression, 20 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in mean and dispersion"))
abline(h=seq(0.88,1,0.02), col='lightgrey'); lines(c(13.5,13.5), c(0.88,1), col='darkgrey')
legend("top", fill=col_vector[1:13], ncol=5, bty='n', cex=0.9, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq", "voom", "DSS", "baySeq", 
                "MDSeq ZI", "MDSeq no ZI", "exp HM", "exp HM log", "ln HM", "ln HM log"))
colMeans(cbind(DE.results.TMM.DEDD2$auc, DE.results.DESeq.DEDD2$auc))
colMeans(cbind(DE.results.TMM.DEDD5$auc, DE.results.DESeq.DEDD5$auc))
colMeans(cbind(DE.results.TMM.DEDD10$auc, DE.results.DESeq.DEDD10$auc))
colMeans(cbind(DE.results.TMM.DEDD20$auc, DE.results.DESeq.DEDD20$auc))
# HMs with DESeq norm equal to edgeR with TMM, similar to DSS and DESeq with DESeq norm, better 
# than all others.


###########
### pAUC ###
par(mfrow=c(2,2), mar=c(1,2,3,1), mgp=c(3,0.7,0))
boxplot(cbind(DE.results.TMM.DE2$pauc, DE.results.DESeq.DE2$pauc), 
        names=NA, col=col_vector[1:13], xaxt='n', ylim=c(0,0.055), 
        pch=20, main=paste0("pAUC for differential expression, 2 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in mean only"))
abline(h=seq(0.004,0.018,0.002), col='lightgrey'); lines(c(13.5,13.5), c(0,0.046), col='darkgrey')
legend("top", fill=col_vector[1:13], ncol=5, bty='n', cex=0.9, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq", "voom", "DSS", "baySeq", 
                "MDSeq ZI", "MDSeq no ZI", "exp HM", "exp HM log", "ln HM", "ln HM log"))
boxplot(cbind(DE.results.TMM.DE5$pauc, DE.results.DESeq.DE5$pauc), 
        names=NA, col=col_vector[1:13], xaxt='n', ylim=c(0,0.055), 
        pch=20, main=paste0("pAUC for differential expression, 5 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in mean only"))
abline(h=seq(0.016,0.032,0.002), col='lightgrey'); lines(c(13.5,13.5), c(0,0.046), col='darkgrey')
legend("top", fill=col_vector[1:13], ncol=5, bty='n', cex=0.9, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq", "voom", "DSS", "baySeq", 
                "MDSeq ZI", "MDSeq no ZI", "exp HM", "exp HM log", "ln HM", "ln HM log"))
boxplot(cbind(DE.results.TMM.DE10$pauc, DE.results.DESeq.DE10$pauc), 
        names=NA, col=col_vector[1:13], xaxt='n', ylim=c(0,0.055), 
        pch=20, main=paste0("pAUC for differential expression, 5 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in mean only"))
abline(h=seq(0.026,0.042,0.002), col='lightgrey'); lines(c(13.5,13.5), c(0,0.046), col='darkgrey')
legend("top", fill=col_vector[1:13], ncol=5, bty='n', cex=0.9, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq", "voom", "DSS", "baySeq", 
                "MDSeq ZI", "MDSeq no ZI", "exp HM", "exp HM log", "ln HM", "ln HM log"))
boxplot(cbind(DE.results.TMM.DE20$pauc, DE.results.DESeq.DE20$pauc), 
        names=NA, col=col_vector[1:13], xaxt='n', ylim=c(0,0.055), 
        pch=20, main=paste0("pAUC for differential expression, 20 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in mean only"))
abline(h=seq(0.032,0.046,0.002), col='lightgrey'); lines(c(13.5,13.5), c(0,0.046), col='darkgrey')
legend("top", fill=col_vector[1:13], ncol=5, bty='n', cex=0.9, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq", "voom", "DSS", "baySeq", 
                "MDSeq ZI", "MDSeq no ZI", "exp HM", "exp HM log", "ln HM", "ln HM log"))
colMeans(cbind(DE.results.TMM.DE2$pauc, DE.results.DESeq.DE2$pauc))
colMeans(cbind(DE.results.TMM.DE5$pauc, DE.results.DESeq.DE5$pauc))
colMeans(cbind(DE.results.TMM.DE10$pauc, DE.results.DESeq.DE10$pauc))
colMeans(cbind(DE.results.TMM.DE20$pauc, DE.results.DESeq.DE20$pauc))

par(mfrow=c(2,2), mar=c(1,2,3,1), mgp=c(3,0.7,0))
boxplot(cbind(DE.results.TMM.DEDD2$pauc, DE.results.DESeq.DEDD2$pauc), 
        names=NA, col=col_vector[1:13], xaxt='n', ylim=c(0,0.055), 
        pch=20, main=paste0("pAUC for differential expression, 2 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in mean and dispersion"))
abline(h=seq(0.004,0.018,0.002), col='lightgrey'); lines(c(13.5,13.5), c(0,0.046), col='darkgrey')
legend("top", fill=col_vector[1:13], ncol=5, bty='n', cex=0.9, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq", "voom", "DSS", "baySeq", 
                "MDSeq ZI", "MDSeq no ZI", "exp HM", "exp HM log", "ln HM", "ln HM log"))
boxplot(cbind(DE.results.TMM.DEDD5$pauc, DE.results.DESeq.DEDD5$pauc), 
        names=NA, col=col_vector[1:13], xaxt='n', ylim=c(0,0.055), 
        pch=20, main=paste0("pAUC for differential expression, 5 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in mean and dispersion"))
abline(h=seq(0.014,0.032,0.002), col='lightgrey'); lines(c(13.5,13.5), c(0,0.046), col='darkgrey')
legend("top", fill=col_vector[1:13], ncol=5, bty='n', cex=0.9, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq", "voom", "DSS", "baySeq", 
                "MDSeq ZI", "MDSeq no ZI", "exp HM", "exp HM log", "ln HM", "ln HM log"))
boxplot(cbind(DE.results.TMM.DEDD10$pauc, DE.results.DESeq.DEDD10$pauc), 
        names=NA, col=col_vector[1:13], xaxt='n', ylim=c(0,0.055), 
        pch=20, main=paste0("pAUC for differential expression, 10 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in mean and dispersion"))
abline(h=seq(0.026,0.042,0.002), col='lightgrey'); lines(c(13.5,13.5), c(0,0.046), col='darkgrey')
legend("top", fill=col_vector[1:13], ncol=5, bty='n', cex=0.9, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq", "voom", "DSS", "baySeq", 
                "MDSeq ZI", "MDSeq no ZI", "exp HM", "exp HM log", "ln HM", "ln HM log"))
boxplot(cbind(DE.results.TMM.DEDD20$pauc, DE.results.DESeq.DEDD20$pauc), 
        names=NA, col=col_vector[1:13], xaxt='n', ylim=c(0,0.055), xlim=c(0.5,26.5), 
        pch=20, main=paste0("pAUC for differential expression, 20 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in mean and dispersion"))
abline(h=seq(0.03,0.046,0.002), col='lightgrey'); lines(c(13.5,13.5), c(0,0.046), col='darkgrey')
legend("top", fill=col_vector[1:13], ncol=5, bty='n', cex=0.9, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq", "voom", "DSS", "baySeq", 
                "MDSeq ZI", "MDSeq no ZI", "exp HM", "exp HM log", "ln HM", "ln HM log"))
colMeans(cbind(DE.results.TMM.DEDD2$pauc, DE.results.DESeq.DEDD2$pauc))
colMeans(cbind(DE.results.TMM.DEDD5$pauc, DE.results.DESeq.DEDD5$pauc))
colMeans(cbind(DE.results.TMM.DEDD10$pauc, DE.results.DESeq.DEDD10$pauc))
colMeans(cbind(DE.results.TMM.DEDD20$pauc, DE.results.DESeq.DEDD20$pauc))


###########
### FDR ###
par(mfrow=c(2,2), mar=c(1,2,3,1), mgp=c(3,0.7,0))
boxplot(cbind(DE.results.TMM.DE2$fdr, DE.results.DESeq.DE2$fdr), 
        names=NA, col=col_vector[1:22], xaxt='n', ylim=c(0,1.3), 
        pch=20, main=paste0("FDR for differential expression, 2 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in mean only"))
abline(h=0.05, col='lightgrey'); lines(c(22.5,22.5), c(0,1), col='darkgrey')
legend("top", fill=col_vector[1:22], ncol=5, bty='n', cex=0.8, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq", "voom", "DSS (FDR)", "DSS (LFDR)", 
                "baySeq", "MDSeq ZI", "MDSeq no ZI", "exp HM (BH)", "exp HM (BY)", "exp HM (q)", 
                "exp HM log (BH)", "exp HM log (BY)", "exp HM log (q)", "ln HM (BH)", "ln HM (BY)", 
                "ln HM (q)", "ln HM log (BH)", "ln HM log (BY)", "ln HM log (q)"))
boxplot(cbind(DE.results.TMM.DE5$fdr, DE.results.DESeq.DE5$fdr), 
        names=NA, col=col_vector[1:22], xaxt='n', ylim=c(0,1.3), 
        pch=20, main=paste0("FDR for differential expression, 5 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in mean only"))
abline(h=0.05, col='lightgrey'); lines(c(22.5,22.5), c(0,1), col='darkgrey')
legend("top", fill=col_vector[1:22], ncol=5, bty='n', cex=0.8, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq", "voom", "DSS (FDR)", "DSS (LFDR)", 
                "baySeq", "MDSeq ZI", "MDSeq no ZI", "exp HM (BH)", "exp HM (BY)", "exp HM (q)", 
                "exp HM log (BH)", "exp HM log (BY)", "exp HM log (q)", "ln HM (BH)", "ln HM (BY)", 
                "ln HM (q)", "ln HM log (BH)", "ln HM log (BY)", "ln HM log (q)"))
boxplot(cbind(DE.results.TMM.DE10$fdr, DE.results.DESeq.DE10$fdr), 
        names=NA, col=col_vector[1:22], xaxt='n', ylim=c(0,1.3), 
        pch=20, main=paste0("FDR for differential expression, 10 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in mean only"))
abline(h=0.05, col='lightgrey'); lines(c(22.5,22.5), c(0,1), col='darkgrey')
legend("top", fill=col_vector[1:22], ncol=5, bty='n', cex=0.8, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq", "voom", "DSS (FDR)", "DSS (LFDR)", 
                "baySeq", "MDSeq ZI", "MDSeq no ZI", "exp HM (BH)", "exp HM (BY)", "exp HM (q)", 
                "exp HM log (BH)", "exp HM log (BY)", "exp HM log (q)", "ln HM (BH)", "ln HM (BY)", 
                "ln HM (q)", "ln HM log (BH)", "ln HM log (BY)", "ln HM log (q)"))
boxplot(cbind(DE.results.TMM.DE20$fdr, DE.results.DESeq.DE20$fdr), 
        names=NA, col=col_vector[1:22], xaxt='n', ylim=c(0,1.3), 
        pch=20, main=paste0("FDR for differential expression, 20 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in mean only"))
abline(h=0.05, col='lightgrey'); lines(c(22.5,22.5), c(0,1), col='darkgrey')
legend("top", fill=col_vector[1:22], ncol=5, bty='n', cex=0.8, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq", "voom", "DSS (FDR)", "DSS (LFDR)", 
                "baySeq", "MDSeq ZI", "MDSeq no ZI", "exp HM (BH)", "exp HM (BY)", "exp HM (q)", 
                "exp HM log (BH)", "exp HM log (BY)", "exp HM log (q)", "ln HM (BH)", "ln HM (BY)", 
                "ln HM (q)", "ln HM log (BH)", "ln HM log (BY)", "ln HM log (q)"))
colMeans(cbind(DE.results.TMM.DE2$fdr, DE.results.DESeq.DE2$fdr), na.rm=T)
colMeans(cbind(DE.results.TMM.DE5$fdr, DE.results.DESeq.DE5$fdr), na.rm=T)
colMeans(cbind(DE.results.TMM.DE10$fdr, DE.results.DESeq.DE10$fdr), na.rm=T)
colMeans(cbind(DE.results.TMM.DE20$fdr, DE.results.DESeq.DE20$fdr), na.rm=T)
# edgeR QL, voom with TMM and baySeq for small samples best for rarely exceeding 0.05 
# without being extremely conservative.
# HMs with BH for DESeq norm and edgeR LR and ET for TMM closest to 0.05 on average.

par(mfrow=c(2,2), mar=c(1,2,3,1), mgp=c(3,0.7,0))
boxplot(cbind(DE.results.TMM.DEDD2$fdr, DE.results.DESeq.DEDD2$fdr), 
        names=NA, col=col_vector[1:22], xaxt='n', ylim=c(0,1.3), 
        pch=20, main=paste0("FDR for differential expression, 2 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in mean and disperison"))
abline(h=0.05, col='lightgrey'); lines(c(22.5,22.5), c(0,1), col='darkgrey')
legend("top", fill=col_vector[1:22], ncol=5, bty='n', cex=0.8, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq", "voom", "DSS (FDR)", "DSS (LFDR)", 
                "baySeq", "MDSeq ZI", "MDSeq no ZI", "exp HM (BH)", "exp HM (BY)", "exp HM (q)", 
                "exp HM log (BH)", "exp HM log (BY)", "exp HM log (q)", "ln HM (BH)", "ln HM (BY)", 
                "ln HM (q)", "ln HM log (BH)", "ln HM log (BY)", "ln HM log (q)"))
boxplot(cbind(DE.results.TMM.DEDD5$fdr, DE.results.DESeq.DEDD5$fdr), 
        names=NA, col=col_vector[1:22], xaxt='n', ylim=c(0,1.3), 
        pch=20, main=paste0("FDR for differential expression, 5 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in mean and disperison"))
abline(h=0.05, col='lightgrey'); lines(c(22.5,22.5), c(0,1), col='darkgrey')
legend("top", fill=col_vector[1:22], ncol=5, bty='n', cex=0.8, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq", "voom", "DSS (FDR)", "DSS (LFDR)", 
                "baySeq", "MDSeq ZI", "MDSeq no ZI", "exp HM (BH)", "exp HM (BY)", "exp HM (q)", 
                "exp HM log (BH)", "exp HM log (BY)", "exp HM log (q)", "ln HM (BH)", "ln HM (BY)", 
                "ln HM (q)", "ln HM log (BH)", "ln HM log (BY)", "ln HM log (q)"))
boxplot(cbind(DE.results.TMM.DEDD10$fdr, DE.results.DESeq.DEDD10$fdr), 
        names=NA, col=col_vector[1:22], xaxt='n', ylim=c(0,1.3), 
        pch=20, main=paste0("FDR for differential expression, 10 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in mean and disperison"))
abline(h=0.05, col='lightgrey'); lines(c(22.5,22.5), c(0,1), col='darkgrey')
legend("top", fill=col_vector[1:22], ncol=5, bty='n', cex=0.8, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq", "voom", "DSS (FDR)", "DSS (LFDR)", 
                "baySeq", "MDSeq ZI", "MDSeq no ZI", "exp HM (BH)", "exp HM (BY)", "exp HM (q)", 
                "exp HM log (BH)", "exp HM log (BY)", "exp HM log (q)", "ln HM (BH)", "ln HM (BY)", 
                "ln HM (q)", "ln HM log (BH)", "ln HM log (BY)", "ln HM log (q)"))
boxplot(cbind(DE.results.TMM.DEDD20$fdr, DE.results.DESeq.DEDD20$fdr), 
        names=NA, col=col_vector[1:22], xaxt='n', ylim=c(0,1.3), 
        pch=20, main=paste0("FDR for differential expression, 20 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in mean and disperison"))
abline(h=0.05, col='lightgrey'); lines(c(22.5,22.5), c(0,1), col='darkgrey')
legend("top", fill=col_vector[1:22], ncol=5, bty='n', cex=0.8, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq", "voom", "DSS (FDR)", "DSS (LFDR)", 
                "baySeq", "MDSeq ZI", "MDSeq no ZI", "exp HM (BH)", "exp HM (BY)", "exp HM (q)", 
                "exp HM log (BH)", "exp HM log (BY)", "exp HM log (q)", "ln HM (BH)", "ln HM (BY)", 
                "ln HM (q)", "ln HM log (BH)", "ln HM log (BY)", "ln HM log (q)"))
colMeans(cbind(DE.results.TMM.DEDD2$fdr, DE.results.DESeq.DEDD2$fdr), na.rm=T)
colMeans(cbind(DE.results.TMM.DEDD5$fdr, DE.results.DESeq.DEDD5$fdr), na.rm=T)
colMeans(cbind(DE.results.TMM.DEDD10$fdr, DE.results.DESeq.DEDD10$fdr), na.rm=T)
colMeans(cbind(DE.results.TMM.DEDD20$fdr, DE.results.DESeq.DEDD20$fdr), na.rm=T)
# HMs with BH, DSS FDR for DESeq norm, and edgeR QL for either normalisation all similar.
# HMs with BH for DESeq norm closest to 0.05 on average.


###########
### TPR ###
par(mfrow=c(2,2), mar=c(1,2,3,1), mgp=c(3,0.7,0))
boxplot(cbind(DE.results.TMM.DE2$tpr, DE.results.DESeq.DE2$tpr), 
        names=NA, col=col_vector[1:22], xaxt='n', ylim=c(0,1.15), 
        pch=20, main=paste0("TPR for differential expression, 2 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in mean only"))
abline(h=seq(0,0.5,0.05), col='lightgrey'); lines(c(22.5,22.5), c(0,0.95), col='darkgrey')
legend("top", fill=col_vector[1:22], ncol=5, bty='n', cex=0.8, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq", "voom", "DSS (FDR)", "DSS (LFDR)", 
                "baySeq", "MDSeq ZI", "MDSeq no ZI", "exp HM (BH)", "exp HM (BY)", "exp HM (q)", 
                "exp HM log (BH)", "exp HM log (BY)", "exp HM log (q)", "ln HM (BH)", "ln HM (BY)", 
                "ln HM (q)", "ln HM log (BH)", "ln HM log (BY)", "ln HM log (q)"))
boxplot(cbind(DE.results.TMM.DE5$tpr, DE.results.DESeq.DE5$tpr), 
        names=NA, col=col_vector[1:22], xaxt='n', ylim=c(0,1.15), 
        pch=20, main=paste0("TPR for differential expression, 5 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in mean only"))
abline(h=seq(0,0.55,0.05), col='lightgrey'); lines(c(22.5,22.5), c(0,0.95), col='darkgrey')
legend("top", fill=col_vector[1:22], ncol=5, bty='n', cex=0.8, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq", "voom", "DSS (FDR)", "DSS (LFDR)", 
                "baySeq", "MDSeq ZI", "MDSeq no ZI", "exp HM (BH)", "exp HM (BY)", "exp HM (q)", 
                "exp HM log (BH)", "exp HM log (BY)", "exp HM log (q)", "ln HM (BH)", "ln HM (BY)", 
                "ln HM (q)", "ln HM log (BH)", "ln HM log (BY)", "ln HM log (q)"))
boxplot(cbind(DE.results.TMM.DE10$tpr, DE.results.DESeq.DE10$tpr), 
        names=NA, col=col_vector[1:22], xaxt='n', ylim=c(0,1.15), 
        pch=20, main=paste0("TPR for differential expression, 10 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in mean only"))
abline(h=seq(0.35,0.75,0.05), col='lightgrey'); lines(c(22.5,22.5), c(0,0.95), col='darkgrey')
legend("top", fill=col_vector[1:22], ncol=5, bty='n', cex=0.8, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq", "voom", "DSS (FDR)", "DSS (LFDR)", 
                "baySeq", "MDSeq ZI", "MDSeq no ZI", "exp HM (BH)", "exp HM (BY)", "exp HM (q)", 
                "exp HM log (BH)", "exp HM log (BY)", "exp HM log (q)", "ln HM (BH)", "ln HM (BY)", 
                "ln HM (q)", "ln HM log (BH)", "ln HM log (BY)", "ln HM log (q)"))
boxplot(cbind(DE.results.TMM.DE20$tpr, DE.results.DESeq.DE20$tpr), 
        names=NA, col=col_vector[1:22], xaxt='n', ylim=c(0,1.15), 
        pch=20, main=paste0("TPR for differential expression, 20 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in mean only"))
abline(h=seq(0.6,0.95,0.05), col='lightgrey'); lines(c(22.5,22.5), c(0,0.95), col='darkgrey')
legend("top", fill=col_vector[1:22], ncol=5, bty='n', cex=0.8, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq", "voom", "DSS (FDR)", "DSS (LFDR)", 
                "baySeq", "MDSeq ZI", "MDSeq no ZI", "exp HM (BH)", "exp HM (BY)", "exp HM (q)", 
                "exp HM log (BH)", "exp HM log (BY)", "exp HM log (q)", "ln HM (BH)", "ln HM (BY)", 
                "ln HM (q)", "ln HM log (BH)", "ln HM log (BY)", "ln HM log (q)"))
colMeans(cbind(DE.results.TMM.DE2$tpr, DE.results.DESeq.DE2$tpr), na.rm=T)
colMeans(cbind(DE.results.TMM.DE5$tpr, DE.results.DESeq.DE5$tpr), na.rm=T)
colMeans(cbind(DE.results.TMM.DE10$tpr, DE.results.DESeq.DE10$tpr), na.rm=T)
colMeans(cbind(DE.results.TMM.DE20$tpr, DE.results.DESeq.DE20$tpr), na.rm=T)
# MDSeq, DESeq with DESeq norm best for 5 samples, then edgeR LR and ET with TMM.
# DESeq with DESeq norm best for 10, then edgeR LR and ET with TMM, then MDSeq and HMs with DESeq 
# norm, edgeR QL with TMM, DSS FDR with DESeq norm.
# Most similar for 20 samples. HMs and DESeq with DESeq norm best, then edgeR and DSS FDR with TMM.

par(mfrow=c(2,2), mar=c(1,2,3,1), mgp=c(3,0.7,0))
boxplot(cbind(DE.results.TMM.DEDD2$tpr, DE.results.DESeq.DEDD2$tpr), 
        names=NA, col=col_vector[1:22], xaxt='n', ylim=c(0,1.15), 
        pch=20, main=paste0("TPR for differential expression, 2 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in mean and dispersion"))
abline(h=seq(0,0.5,0.05), col='lightgrey'); lines(c(22.5,22.5), c(0,0.95), col='darkgrey')
legend("top", fill=col_vector[1:22], ncol=5, bty='n', cex=0.8, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq", "voom", "DSS (FDR)", "DSS (LFDR)", 
                "baySeq", "MDSeq ZI", "MDSeq no ZI", "exp HM (BH)", "exp HM (BY)", "exp HM (q)", 
                "exp HM log (BH)", "exp HM log (BY)", "exp HM log (q)", "ln HM (BH)", "ln HM (BY)", 
                "ln HM (q)", "ln HM log (BH)", "ln HM log (BY)", "ln HM log (q)"))
boxplot(cbind(DE.results.TMM.DEDD5$tpr, DE.results.DESeq.DEDD5$tpr), 
        names=NA, col=col_vector[1:22], xaxt='n', ylim=c(0,1.15), 
        pch=20, main=paste0("TPR for differential expression, 2 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in mean and dispersion"))
abline(h=seq(0,0.55,0.05), col='lightgrey'); lines(c(22.5,22.5), c(0,0.95), col='darkgrey')
legend("top", fill=col_vector[1:22], ncol=5, bty='n', cex=0.8, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq", "voom", "DSS (FDR)", "DSS (LFDR)", 
                "baySeq", "MDSeq ZI", "MDSeq no ZI", "exp HM (BH)", "exp HM (BY)", "exp HM (q)", 
                "exp HM log (BH)", "exp HM log (BY)", "exp HM log (q)", "ln HM (BH)", "ln HM (BY)", 
                "ln HM (q)", "ln HM log (BH)", "ln HM log (BY)", "ln HM log (q)"))
boxplot(cbind(DE.results.TMM.DEDD10$tpr, DE.results.DESeq.DEDD10$tpr), 
        names=NA, col=col_vector[1:22], xaxt='n', ylim=c(0,1.15), 
        pch=20, main=paste0("TPR for differential expression, 2 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in mean and dispersion"))
abline(h=seq(0.3,0.75,0.05), col='lightgrey'); lines(c(22.5,22.5), c(0,0.95), col='darkgrey')
legend("top", fill=col_vector[1:22], ncol=5, bty='n', cex=0.8, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq", "voom", "DSS (FDR)", "DSS (LFDR)", 
                "baySeq", "MDSeq ZI", "MDSeq no ZI", "exp HM (BH)", "exp HM (BY)", "exp HM (q)", 
                "exp HM log (BH)", "exp HM log (BY)", "exp HM log (q)", "ln HM (BH)", "ln HM (BY)", 
                "ln HM (q)", "ln HM log (BH)", "ln HM log (BY)", "ln HM log (q)"))
boxplot(cbind(DE.results.TMM.DEDD20$tpr, DE.results.DESeq.DEDD20$tpr), 
        names=NA, col=col_vector[1:22], xaxt='n', ylim=c(0,1.15), 
        pch=20, main=paste0("TPR for differential expression, 2 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in mean and dispersion"))
abline(h=seq(0.6,0.95,0.05), col='lightgrey'); lines(c(22.5,22.5), c(0,0.95), col='darkgrey')
legend("top", fill=col_vector[1:22], ncol=5, bty='n', cex=0.8, 
       legend=c("edgeR QL", "edgeR LR", "edgeR ET", "DESeq", "voom", "DSS (FDR)", "DSS (LFDR)", 
                "baySeq", "MDSeq ZI", "MDSeq no ZI", "exp HM (BH)", "exp HM (BY)", "exp HM (q)", 
                "exp HM log (BH)", "exp HM log (BY)", "exp HM log (q)", "ln HM (BH)", "ln HM (BY)", 
                "ln HM (q)", "ln HM log (BH)", "ln HM log (BY)", "ln HM log (q)"))
colMeans(cbind(DE.results.TMM.DEDD2$tpr, DE.results.DESeq.DEDD2$tpr), na.rm=T)
colMeans(cbind(DE.results.TMM.DEDD5$tpr, DE.results.DESeq.DEDD5$tpr), na.rm=T)
colMeans(cbind(DE.results.TMM.DEDD10$tpr, DE.results.DESeq.DEDD10$tpr), na.rm=T)
colMeans(cbind(DE.results.TMM.DEDD20$tpr, DE.results.DESeq.DEDD20$tpr), na.rm=T)
# MDSeq, DESeq with DESeq norm best for 5 samples, then edgeR LR, ET with TMM.
# DESeq with DESeq norm and edgeR LR, ET with TMM best for 10, then HMs with DESeq norm, edgeR QL 
# with TMM, MDSeq with DESeq norm, DSS FDR with TMM.
# HMs with either normalisation, edgeR with TMM, DESeq with DESeq norm, DSS with either 
# normalisation all similar for 20 samples.


###################################
#### Differential distribution ####
###################################

###########
### AUC ###
par(mfrow=c(2,2), mar=c(1,2,3,1), mgp=c(3,0.7,0))
boxplot(cbind(DEDD.results.TMM.DD2$auc, DEDD.results.DESeq.DD2$auc), 
        names=NA, col=col_vector[1:2], xaxt='n', ylim=c(0.45,1.05), 
        pch=20, main=paste0("AUC for differential distribution, 2 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in dispersion only"))
abline(h=seq(0.46, 0.54, 0.02), col='lightgrey'); lines(c(2.5,2.5), c(0,1), col='darkgrey')
legend("top", fill=col_vector[1:2], ncol=3, bty='n', cex=1.2, 
       legend=c("exp HMM", "ln HMM"))
boxplot(cbind(DEDD.results.TMM.DD5$auc, DEDD.results.DESeq.DD5$auc), 
        names=NA, col=col_vector[1:2], xaxt='n', ylim=c(0.45,1.05), 
        pch=20, main=paste0("AUC for differential distribution, 5 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in dispersion only"))
abline(h=seq(0.48, 0.54, 0.02), col='lightgrey'); lines(c(2.5,2.5), c(0,1), col='darkgrey')
legend("top", fill=col_vector[1:2], ncol=3, bty='n', cex=1.2, 
       legend=c("exp HMM", "ln HMM"))
boxplot(cbind(DEDD.results.TMM.DD10$auc, DEDD.results.DESeq.DD10$auc), 
        names=NA, col=col_vector[1:2], xaxt='n', ylim=c(0.45,1.05), 
        pch=20, main=paste0("AUC for differential distribution, 10 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in dispersion only"))
abline(h=seq(0.5, 0.56, 0.02), col='lightgrey'); lines(c(2.5,2.5), c(0,1), col='darkgrey')
legend("top", fill=col_vector[1:2], ncol=3, bty='n', cex=1.2, 
       legend=c("exp HMM", "ln HMM"))
boxplot(cbind(DEDD.results.TMM.DD20$auc, DEDD.results.DESeq.DD20$auc), 
        names=NA, col=col_vector[1:2], xaxt='n', ylim=c(0.45,1.05), 
        pch=20, main=paste0("AUC for differential distribution, 20 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in dispersion only"))
abline(h=seq(0.52, 0.62, 0.02), col='lightgrey'); lines(c(2.5,2.5), c(0,1), col='darkgrey')
legend("top", fill=col_vector[1:2], ncol=3, bty='n', cex=1.2, 
       legend=c("exp HMM", "ln HMM"))
colMeans(cbind(DEDD.results.TMM.DD2$auc, DEDD.results.DESeq.DD2$auc))
colMeans(cbind(DEDD.results.TMM.DD5$auc, DEDD.results.DESeq.DD5$auc))
colMeans(cbind(DEDD.results.TMM.DD10$auc, DEDD.results.DESeq.DD10$auc))
colMeans(cbind(DEDD.results.TMM.DD20$auc, DEDD.results.DESeq.DD20$auc))
# Better with DESeq norm, improvement increases with sample size.
# Very little difference between expHM and lnHM.

par(mfrow=c(2,2), mar=c(1,2,3,1), mgp=c(3,0.7,0))
boxplot(cbind(DEDD.results.TMM.DE2$auc, DEDD.results.DESeq.DE2$auc), 
        names=NA, col=col_vector[1:2], xaxt='n', ylim=c(0.45,1.05), 
        pch=20, main=paste0("AUC for differential distribution, 2 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in mean only"))
abline(h=seq(0.55, 0.75, 0.05), col='lightgrey'); lines(c(2.5,2.5), c(0,1), col='darkgrey')
legend("top", fill=col_vector[1:2], ncol=3, bty='n', cex=1.2, 
       legend=c("exp HMM", "ln HMM"))
boxplot(cbind(DEDD.results.TMM.DE5$auc, DEDD.results.DESeq.DE5$auc), 
        names=NA, col=col_vector[1:2], xaxt='n', ylim=c(0.45,1.05), 
        pch=20, main=paste0("AUC for differential distribution, 5 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in mean only"))
abline(h=seq(0.7, 0.9, 0.05), col='lightgrey'); lines(c(2.5,2.5), c(0,1), col='darkgrey')
legend("top", fill=col_vector[1:2], ncol=3, bty='n', cex=1.2, 
       legend=c("exp HMM", "ln HMM"))
boxplot(cbind(DEDD.results.TMM.DE10$auc, DEDD.results.DESeq.DE10$auc), 
        names=NA, col=col_vector[1:2], xaxt='n', ylim=c(0.45,1.05), 
        pch=20, main=paste0("AUC for differential distribution, 10 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in mean only"))
abline(h=seq(0.78, 0.92, 0.02), col='lightgrey'); lines(c(2.5,2.5), c(0,1), col='darkgrey')
legend("top", fill=col_vector[1:2], ncol=3, bty='n', cex=1.2, 
       legend=c("exp HMM", "ln HMM"))
boxplot(cbind(DEDD.results.TMM.DE20$auc, DEDD.results.DESeq.DE20$auc), 
        names=NA, col=col_vector[1:2], xaxt='n', ylim=c(0.45,1.05), 
        pch=20, main=paste0("AUC for differential distribution, 20 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in mean only"))
abline(h=seq(0.88, 0.96, 0.02), col='lightgrey'); lines(c(2.5,2.5), c(0,1), col='darkgrey')
legend("top", fill=col_vector[1:2], ncol=3, bty='n', cex=1.2, 
       legend=c("exp HMM", "ln HMM"))
colMeans(cbind(DEDD.results.TMM.DE2$auc, DEDD.results.DESeq.DE2$auc))
colMeans(cbind(DEDD.results.TMM.DE5$auc, DEDD.results.DESeq.DE5$auc))
colMeans(cbind(DEDD.results.TMM.DE10$auc, DEDD.results.DESeq.DE10$auc))
colMeans(cbind(DEDD.results.TMM.DE20$auc, DEDD.results.DESeq.DE20$auc))
# Better with DESeq norm but improvement decreases with sample size, presumably because for bigger 
# samples DE is easy to detect.
# lnHM clearly better than expHM for small samples, similar for large, presumably for same reason.

par(mfrow=c(2,2), mar=c(1,2,3,1), mgp=c(3,0.7,0))
boxplot(cbind(DEDD.results.TMM.DEDD2$auc, DEDD.results.DESeq.DEDD2$auc), 
        names=NA, col=col_vector[1:2], xaxt='n', ylim=c(0.45,1.05), 
        pch=20, main=paste0("AUC for differential distribution, 2 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in mean and dispersion"))
abline(h=seq(0.56, 0.68, 0.02), col='lightgrey'); lines(c(2.5,2.5), c(0,1), col='darkgrey')
legend("top", fill=col_vector[1:2], ncol=3, bty='n', cex=1.2, 
       legend=c("exp HMM", "ln HMM"))
boxplot(cbind(DEDD.results.TMM.DEDD5$auc, DEDD.results.DESeq.DEDD5$auc), 
        names=NA, col=col_vector[1:2], xaxt='n', ylim=c(0.45,1.05), 
        pch=20, main=paste0("AUC for differential distribution, 5 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in mean and dispersion"))
abline(h=seq(0.64, 0.76, 0.02), col='lightgrey'); lines(c(2.5,2.5), c(0,1), col='darkgrey')
legend("top", fill=col_vector[1:2], ncol=3, bty='n', cex=1.2, 
       legend=c("exp HMM", "ln HMM"))
boxplot(cbind(DEDD.results.TMM.DEDD10$auc, DEDD.results.DESeq.DEDD10$auc), 
        names=NA, col=col_vector[1:2], xaxt='n', ylim=c(0.45,1.05), 
        pch=20, main=paste0("AUC for differential distribution, 10 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in mean and dispersion"))
abline(h=seq(0.74, 0.82, 0.02), col='lightgrey'); lines(c(2.5,2.5), c(0,1), col='darkgrey')
legend("top", fill=col_vector[1:2], ncol=3, bty='n', cex=1.2, 
       legend=c("exp HMM", "ln HMM"))
boxplot(cbind(DEDD.results.TMM.DEDD20$auc, DEDD.results.DESeq.DEDD20$auc), 
        names=NA, col=col_vector[1:2], xaxt='n', ylim=c(0.45,1.05), 
        pch=20, main=paste0("AUC for differential distribution, 20 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in mean and dispersion"))
abline(h=seq(0.8, 0.88, 0.02), col='lightgrey'); lines(c(2.5,2.5), c(0,1), col='darkgrey')
legend("top", fill=col_vector[1:2], ncol=3, bty='n', cex=1.2, 
       legend=c("exp HMM", "ln HMM"))
colMeans(cbind(DEDD.results.TMM.DEDD2$auc, DEDD.results.DESeq.DEDD2$auc))
colMeans(cbind(DEDD.results.TMM.DEDD5$auc, DEDD.results.DESeq.DEDD5$auc))
colMeans(cbind(DEDD.results.TMM.DEDD10$auc, DEDD.results.DESeq.DEDD10$auc))
colMeans(cbind(DEDD.results.TMM.DEDD20$auc, DEDD.results.DESeq.DEDD20$auc))
# Better with DESeq norm, no clear trend in improvement with increasing sample size.
# Higher variance for small samples with TMM, but similar between normalisations for bigger samples.
# lhHM better than expHM for small samples, similar or worse for bigger samples.


############
### pAUC ###
par(mfrow=c(2,2), mar=c(1,2,3,1), mgp=c(3,0.7,0))
boxplot(cbind(DEDD.results.TMM.DD2$pauc, DEDD.results.DESeq.DD2$pauc), 
        names=NA, col=col_vector[1:2], xaxt='n', ylim=c(0,0.05), 
        pch=20, main=paste0("pAUC for differential distribution, 2 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in dispersion only"))
abline(h=seq(0,0.003,0.001), col='lightgrey'); lines(c(2.5,2.5), c(0,0.045), col='darkgrey')
legend("top", fill=col_vector[1:2], ncol=3, bty='n', cex=1.2, 
       legend=c("exp HMM", "ln HMM"))
boxplot(cbind(DEDD.results.TMM.DD5$pauc, DEDD.results.DESeq.DD5$pauc), 
        names=NA, col=col_vector[1:2], xaxt='n', ylim=c(0,0.05), 
        pch=20, main=paste0("pAUC for differential distribution, 5 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in dispersion only"))
abline(h=seq(0,0.003,0.001), col='lightgrey'); lines(c(2.5,2.5), c(0,0.045), col='darkgrey')
legend("top", fill=col_vector[1:2], ncol=3, bty='n', cex=1.2, 
       legend=c("exp HMM", "ln HMM"))
boxplot(cbind(DEDD.results.TMM.DD10$pauc, DEDD.results.DESeq.DD10$pauc), 
        names=NA, col=col_vector[1:2], xaxt='n', ylim=c(0,0.05), 
        pch=20, main=paste0("pAUC for differential distribution, 10 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in dispersion only"))
abline(h=seq(0.001,0.005,0.001), col='lightgrey'); lines(c(2.5,2.5), c(0,0.045), col='darkgrey')
legend("top", fill=col_vector[1:2], ncol=3, bty='n', cex=1.2, 
       legend=c("exp HMM", "ln HMM"))
boxplot(cbind(DEDD.results.TMM.DD20$pauc, DEDD.results.DESeq.DD20$pauc), 
        names=NA, col=col_vector[1:2], xaxt='n', ylim=c(0,0.05), 
        pch=20, main=paste0("pAUC for differential distribution, 20 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in dispersion only"))
abline(h=seq(0.002,0.01,0.002), col='lightgrey'); lines(c(2.5,2.5), c(0,0.045), col='darkgrey')
legend("top", fill=col_vector[1:2], ncol=3, bty='n', cex=1.2, 
       legend=c("exp HMM", "ln HMM"))
colMeans(cbind(DEDD.results.TMM.DD2$pauc, DEDD.results.DESeq.DD2$pauc))
colMeans(cbind(DEDD.results.TMM.DD5$pauc, DEDD.results.DESeq.DD5$pauc))
colMeans(cbind(DEDD.results.TMM.DD10$pauc, DEDD.results.DESeq.DD10$pauc))
colMeans(cbind(DEDD.results.TMM.DD20$pauc, DEDD.results.DESeq.DD20$pauc))

par(mfrow=c(2,2), mar=c(1,2,3,1), mgp=c(3,0.7,0))
boxplot(cbind(DEDD.results.TMM.DE2$pauc, DEDD.results.DESeq.DE2$pauc), 
        names=NA, col=col_vector[1:2], xaxt='n', ylim=c(0,0.05), 
        pch=20, main=paste0("pAUC for differential distribution, 2 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in mean only"))
abline(h=seq(0.002,0.01,0.002), col='lightgrey'); lines(c(2.5,2.5), c(0,0.045), col='darkgrey')
legend("top", fill=col_vector[1:2], ncol=3, bty='n', cex=1.2, 
       legend=c("exp HMM", "ln HMM"))
boxplot(cbind(DEDD.results.TMM.DE5$pauc, DEDD.results.DESeq.DE5$pauc), 
        names=NA, col=col_vector[1:2], xaxt='n', ylim=c(0,0.05), 
        pch=20, main=paste0("pAUC for differential distribution, 5 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in mean only"))
abline(h=seq(0.014,0.028,0.002), col='lightgrey'); lines(c(2.5,2.5), c(0,0.045), col='darkgrey')
legend("top", fill=col_vector[1:2], ncol=3, bty='n', cex=1.2, 
       legend=c("exp HMM", "ln HMM"))
boxplot(cbind(DEDD.results.TMM.DE10$pauc, DEDD.results.DESeq.DE10$pauc), 
        names=NA, col=col_vector[1:2], xaxt='n', ylim=c(0,0.05), 
        pch=20, main=paste0("pAUC for differential distribution, 10 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in mean only"))
abline(h=seq(0.024,0.038,0.002), col='lightgrey'); lines(c(2.5,2.5), c(0,0.045), col='darkgrey')
legend("top", fill=col_vector[1:2], ncol=3, bty='n', cex=1.2, 
       legend=c("exp HMM", "ln HMM"))
boxplot(cbind(DEDD.results.TMM.DE20$pauc, DEDD.results.DESeq.DE20$pauc), 
        names=NA, col=col_vector[1:2], xaxt='n', ylim=c(0,0.05), 
        pch=20, main=paste0("pAUC for differential distribution, 20 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in mean only"))
abline(h=seq(0.034,0.044,0.002), col='lightgrey'); lines(c(2.5,2.5), c(0,0.045), col='darkgrey')
legend("top", fill=col_vector[1:2], ncol=3, bty='n', cex=1.2, 
       legend=c("exp HMM", "ln HMM"))
colMeans(cbind(DEDD.results.TMM.DE2$pauc, DEDD.results.DESeq.DE2$pauc))
colMeans(cbind(DEDD.results.TMM.DE5$pauc, DEDD.results.DESeq.DE5$pauc))
colMeans(cbind(DEDD.results.TMM.DE10$pauc, DEDD.results.DESeq.DE10$pauc))
colMeans(cbind(DEDD.results.TMM.DE20$pauc, DEDD.results.DESeq.DE20$pauc))

par(mfrow=c(2,2), mar=c(1,2,3,1), mgp=c(3,0.7,0))
boxplot(cbind(DEDD.results.TMM.DEDD2$pauc, DEDD.results.DESeq.DEDD2$pauc), 
        names=NA, col=col_vector[1:2], xaxt='n', ylim=c(0,0.05), 
        pch=20, main=paste0("pAUC for differential distribution, 2 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in mean and dispersion"))
abline(h=seq(0.003,0.007,0.001), col='lightgrey'); lines(c(2.5,2.5), c(0,0.045), col='darkgrey')
legend("top", fill=col_vector[1:2], ncol=3, bty='n', cex=1.2, 
       legend=c("exp HMM", "ln HMM"))
boxplot(cbind(DEDD.results.TMM.DEDD5$pauc, DEDD.results.DESeq.DEDD5$pauc), 
        names=NA, col=col_vector[1:2], xaxt='n', ylim=c(0,0.05), 
        pch=20, main=paste0("pAUC for differential distribution, 5 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in mean and dispersion"))
abline(h=seq(0.01,0.02,0.002), col='lightgrey'); lines(c(2.5,2.5), c(0,0.045), col='darkgrey')
legend("top", fill=col_vector[1:2], ncol=3, bty='n', cex=1.2, 
       legend=c("exp HMM", "ln HMM"))
boxplot(cbind(DEDD.results.TMM.DEDD10$pauc, DEDD.results.DESeq.DEDD10$pauc), 
        names=NA, col=col_vector[1:2], xaxt='n', ylim=c(0,0.05), 
        pch=20, main=paste0("pAUC for differential distribution, 10 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in mean and dispersion"))
abline(h=seq(0.018,0.026,0.002), col='lightgrey'); lines(c(2.5,2.5), c(0,0.045), col='darkgrey')
legend("top", fill=col_vector[1:2], ncol=3, bty='n', cex=1.2, 
       legend=c("exp HMM", "ln HMM"))
boxplot(cbind(DEDD.results.TMM.DEDD20$pauc, DEDD.results.DESeq.DEDD20$pauc), 
        names=NA, col=col_vector[1:2], xaxt='n', ylim=c(0,0.05), 
        pch=20, main=paste0("pAUC for differential distribution, 20 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in mean and dispersion"))
abline(h=seq(0.024,0.034,0.002), col='lightgrey'); lines(c(2.5,2.5), c(0,0.045), col='darkgrey')
legend("top", fill=col_vector[1:2], ncol=3, bty='n', cex=1.2, 
       legend=c("exp HMM", "ln HMM"))
colMeans(cbind(DEDD.results.TMM.DEDD2$pauc, DEDD.results.DESeq.DEDD2$pauc))
colMeans(cbind(DEDD.results.TMM.DEDD5$pauc, DEDD.results.DESeq.DEDD5$pauc))
colMeans(cbind(DEDD.results.TMM.DEDD10$pauc, DEDD.results.DESeq.DEDD10$pauc))
colMeans(cbind(DEDD.results.TMM.DEDD20$pauc, DEDD.results.DESeq.DEDD20$pauc))


###########
### FDR ###
par(mfrow=c(2,2), mar=c(1,2,3,1), mgp=c(3,0.7,0))
boxplot(cbind(DEDD.results.TMM.DD2$fdr, DEDD.results.DESeq.DD2$fdr), 
        names=NA, col=col_vector[1:6], xaxt='n', ylim=c(0,1.2), 
        pch=20, main=paste0("FDR for differential distribution, 2 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in dispersion only"))
abline(h=0.05, col='lightgrey'); lines(c(6.5,6.5), c(0,1), col='darkgrey')
legend("top", fill=col_vector[1:6], ncol=4, bty='n', cex=0.9, 
       legend=c("exp HMM (0.5)", "exp HMM (thr)", "exp HMM (BFDR)", "ln HMM (0.5)", "ln HMM (thr)", 
                "ln HMM (BFDR)"))
boxplot(cbind(DEDD.results.TMM.DD5$fdr, DEDD.results.DESeq.DD5$fdr), 
        names=NA, col=col_vector[1:6], xaxt='n', ylim=c(0,1.2), 
        pch=20, main=paste0("FDR for differential distribution, 5 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in dispersion only"))
abline(h=0.05, col='lightgrey'); lines(c(6.5,6.5), c(0,1), col='darkgrey')
legend("top", fill=col_vector[1:6], ncol=4, bty='n', cex=0.9, 
       legend=c("exp HMM (0.5)", "exp HMM (thr)", "exp HMM (BFDR)", "ln HMM (0.5)", "ln HMM (thr)", 
                "ln HMM (BFDR)"))
boxplot(cbind(DEDD.results.TMM.DD10$fdr, DEDD.results.DESeq.DD10$fdr), 
        names=NA, col=col_vector[1:6], xaxt='n', ylim=c(0,1.2), 
        pch=20, main=paste0("FDR for differential distribution, 10 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in dispersion only"))
abline(h=0.05, col='lightgrey'); lines(c(6.5,6.5), c(0,1), col='darkgrey')
legend("top", fill=col_vector[1:6], ncol=4, bty='n', cex=0.9, 
       legend=c("exp HMM (0.5)", "exp HMM (thr)", "exp HMM (BFDR)", "ln HMM (0.5)", "ln HMM (thr)", 
                "ln HMM (BFDR)"))
boxplot(cbind(DEDD.results.TMM.DD20$fdr, DEDD.results.DESeq.DD20$fdr), 
        names=NA, col=col_vector[1:6], xaxt='n', ylim=c(0,1.2), 
        pch=20, main=paste0("FDR for differential distribution, 20 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in dispersion only"))
abline(h=0.05, col='lightgrey'); lines(c(6.5,6.5), c(0,1), col='darkgrey')
legend("top", fill=col_vector[1:6], ncol=4, bty='n', cex=0.9, 
       legend=c("exp HMM (0.5)", "exp HMM (thr)", "exp HMM (BFDR)", "ln HMM (0.5)", "ln HMM (thr)", 
                "ln HMM (BFDR)"))
colMeans(cbind(DEDD.results.TMM.DD2$fdr, DEDD.results.DESeq.DD2$fdr), na.rm=T)
colMeans(cbind(DEDD.results.TMM.DD5$fdr, DEDD.results.DESeq.DD5$fdr), na.rm=T)
colMeans(cbind(DEDD.results.TMM.DD10$fdr, DEDD.results.DESeq.DD10$fdr), na.rm=T)
colMeans(cbind(DEDD.results.TMM.DD20$fdr, DEDD.results.DESeq.DD20$fdr), na.rm=T)
# Hard to take much from this. FDRs really unreliable for these sample sizes for differences in 
# dispersion only. Looks like posterior threshold is probably best, and all look better with 
# DESeq norm than TMM.

par(mfrow=c(2,2), mar=c(1,2,3,1), mgp=c(3,0.7,0))
boxplot(cbind(DEDD.results.TMM.DE2$fdr, DEDD.results.DESeq.DE2$fdr), 
        names=NA, col=col_vector[1:6], xaxt='n', ylim=c(0,1.2), 
        pch=20, main=paste0("FDR for differential distribution, 2 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in mean only"))
abline(h=0.05, col='lightgrey'); lines(c(6.5,6.5), c(0,1), col='darkgrey')
legend("top", fill=col_vector[1:6], ncol=4, bty='n', cex=0.9, 
       legend=c("exp HMM (0.5)", "exp HMM (thr)", "exp HMM (BFDR)", "ln HMM (0.5)", "ln HMM (thr)", 
                "ln HMM (BFDR)"))
boxplot(cbind(DEDD.results.TMM.DE5$fdr, DEDD.results.DESeq.DE5$fdr), 
        names=NA, col=col_vector[1:6], xaxt='n', ylim=c(0,1.2), 
        pch=20, main=paste0("FDR for differential distribution, 5 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in mean only"))
abline(h=0.05, col='lightgrey'); lines(c(6.5,6.5), c(0,1), col='darkgrey')
legend("top", fill=col_vector[1:6], ncol=4, bty='n', cex=0.9, 
       legend=c("exp HMM (0.5)", "exp HMM (thr)", "exp HMM (BFDR)", "ln HMM (0.5)", "ln HMM (thr)", 
                "ln HMM (BFDR)"))
boxplot(cbind(DEDD.results.TMM.DE10$fdr, DEDD.results.DESeq.DE10$fdr), 
        names=NA, col=col_vector[1:6], xaxt='n', ylim=c(0,1.2), 
        pch=20, main=paste0("FDR for differential distribution, 10 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in mean only"))
abline(h=0.05, col='lightgrey'); lines(c(6.5,6.5), c(0,1), col='darkgrey')
legend("top", fill=col_vector[1:6], ncol=4, bty='n', cex=0.9, 
       legend=c("exp HMM (0.5)", "exp HMM (thr)", "exp HMM (BFDR)", "ln HMM (0.5)", "ln HMM (thr)", 
                "ln HMM (BFDR)"))
boxplot(cbind(DEDD.results.TMM.DE20$fdr, DEDD.results.DESeq.DE20$fdr), 
        names=NA, col=col_vector[1:6], xaxt='n', ylim=c(0,1.2), 
        pch=20, main=paste0("FDR for differential distribution, 20 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in mean only"))
abline(h=0.05, col='lightgrey'); lines(c(6.5,6.5), c(0,1), col='darkgrey')
legend("top", fill=col_vector[1:6], ncol=4, bty='n', cex=0.9, 
       legend=c("exp HMM (0.5)", "exp HMM (thr)", "exp HMM (BFDR)", "ln HMM (0.5)", "ln HMM (thr)", 
                "ln HMM (BFDR)"))
colMeans(cbind(DEDD.results.TMM.DE2$fdr, DEDD.results.DESeq.DE2$fdr), na.rm=T)
colMeans(cbind(DEDD.results.TMM.DE5$fdr, DEDD.results.DESeq.DE5$fdr), na.rm=T)
colMeans(cbind(DEDD.results.TMM.DE10$fdr, DEDD.results.DESeq.DE10$fdr), na.rm=T)
colMeans(cbind(DEDD.results.TMM.DE20$fdr, DEDD.results.DESeq.DE20$fdr), na.rm=T)
# Not much difference between normalisation methods for bigger samples. No method gives good 
# FDRs - posterior and 0.5 threshold get extremely conservative as sample size increases.
# BFDR always zero.

par(mfrow=c(2,2), mar=c(1,2,3,1), mgp=c(3,0.7,0))
boxplot(cbind(DEDD.results.TMM.DEDD2$fdr, DEDD.results.DESeq.DEDD2$fdr), 
        names=NA, col=col_vector[1:6], xaxt='n', ylim=c(0,1.2), 
        pch=20, main=paste0("FDR for differential distribution, 2 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in mean and dispersion"))
abline(h=0.05, col='lightgrey'); lines(c(6.5,6.5), c(0,1), col='darkgrey')
legend("top", fill=col_vector[1:6], ncol=4, bty='n', cex=0.9, 
       legend=c("exp HMM (0.5)", "exp HMM (thr)", "exp HMM (BFDR)", "ln HMM (0.5)", "ln HMM (thr)", 
                "ln HMM (BFDR)"))
boxplot(cbind(DEDD.results.TMM.DEDD5$fdr, DEDD.results.DESeq.DEDD5$fdr), 
        names=NA, col=col_vector[1:6], xaxt='n', ylim=c(0,1.2), 
        pch=20, main=paste0("FDR for differential distribution, 5 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in mean and dispersion"))
abline(h=0.05, col='lightgrey'); lines(c(6.5,6.5), c(0,1), col='darkgrey')
legend("top", fill=col_vector[1:6], ncol=4, bty='n', cex=0.9, 
       legend=c("exp HMM (0.5)", "exp HMM (thr)", "exp HMM (BFDR)", "ln HMM (0.5)", "ln HMM (thr)", 
                "ln HMM (BFDR)"))
boxplot(cbind(DEDD.results.TMM.DEDD10$fdr, DEDD.results.DESeq.DEDD10$fdr), 
        names=NA, col=col_vector[1:6], xaxt='n', ylim=c(0,1.2), 
        pch=20, main=paste0("FDR for differential distribution, 10 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in mean and dispersion"))
abline(h=0.05, col='lightgrey'); lines(c(6.5,6.5), c(0,1), col='darkgrey')
legend("top", fill=col_vector[1:6], ncol=4, bty='n', cex=0.9, 
       legend=c("exp HMM (0.5)", "exp HMM (thr)", "exp HMM (BFDR)", "ln HMM (0.5)", "ln HMM (thr)", 
                "ln HMM (BFDR)"))
boxplot(cbind(DEDD.results.TMM.DEDD20$fdr, DEDD.results.DESeq.DEDD20$fdr), 
        names=NA, col=col_vector[1:6], xaxt='n', ylim=c(0,1.2), 
        pch=20, main=paste0("FDR for differential distribution, 20 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in mean and dispersion"))
abline(h=0.05, col='lightgrey'); lines(c(6.5,6.5), c(0,1), col='darkgrey')
legend("top", fill=col_vector[1:6], ncol=4, bty='n', cex=0.9, 
       legend=c("exp HMM (0.5)", "exp HMM (thr)", "exp HMM (BFDR)", "ln HMM (0.5)", "ln HMM (thr)", 
                "ln HMM (BFDR)"))
colMeans(cbind(DEDD.results.TMM.DEDD2$fdr, DEDD.results.DESeq.DEDD2$fdr), na.rm=T)
colMeans(cbind(DEDD.results.TMM.DEDD5$fdr, DEDD.results.DESeq.DEDD5$fdr), na.rm=T)
colMeans(cbind(DEDD.results.TMM.DEDD10$fdr, DEDD.results.DESeq.DEDD10$fdr), na.rm=T)
colMeans(cbind(DEDD.results.TMM.DEDD20$fdr, DEDD.results.DESeq.DEDD20$fdr), na.rm=T)
# Similar to DE only: not much difference by normalisation, BFDR all zero, threshold methods 
# increasingly conservative as sample size increases.


###########
### TPR ###
par(mfrow=c(2,2), mar=c(1,2,3,1), mgp=c(3,0.7,0))
boxplot(cbind(DEDD.results.TMM.DD2$tpr, DEDD.results.DESeq.DD2$tpr), 
        names=NA, col=col_vector[1:6], xaxt='n', ylim=c(0,1.2), 
        pch=20, main=paste0("TPR for differential distribution, 2 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in dispersion only"))
lines(c(6.5,6.5), c(0,1), col='darkgrey')
legend("top", fill=col_vector[1:6], ncol=4, bty='n', cex=0.9, 
       legend=c("exp HMM (0.5)", "exp HMM (thr)", "exp HMM (BFDR)", "ln HMM (0.5)", "ln HMM (thr)", 
                "ln HMM (BFDR)"))
boxplot(cbind(DEDD.results.TMM.DD5$tpr, DEDD.results.DESeq.DD5$tpr), 
        names=NA, col=col_vector[1:6], xaxt='n', ylim=c(0,1.2), 
        pch=20, main=paste0("TPR for differential distribution, 5 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in dispersion only"))
lines(c(6.5,6.5), c(0,1), col='darkgrey')
legend("top", fill=col_vector[1:6], ncol=4, bty='n', cex=0.9, 
       legend=c("exp HMM (0.5)", "exp HMM (thr)", "exp HMM (BFDR)", "ln HMM (0.5)", "ln HMM (thr)", 
                "ln HMM (BFDR)"))
boxplot(cbind(DEDD.results.TMM.DD10$tpr, DEDD.results.DESeq.DD10$tpr), 
        names=NA, col=col_vector[1:6], xaxt='n', ylim=c(0,1.2), 
        pch=20, main=paste0("TPR for differential distribution, 10 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in dispersion only"))
lines(c(6.5,6.5), c(0,1), col='darkgrey')
legend("top", fill=col_vector[1:6], ncol=4, bty='n', cex=0.9, 
       legend=c("exp HMM (0.5)", "exp HMM (thr)", "exp HMM (BFDR)", "ln HMM (0.5)", "ln HMM (thr)", 
                "ln HMM (BFDR)"))
boxplot(cbind(DEDD.results.TMM.DD20$tpr, DEDD.results.DESeq.DD20$tpr), 
        names=NA, col=col_vector[1:6], xaxt='n', ylim=c(0,1.2), 
        pch=20, main=paste0("TPR for differential distribution, 20 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in dispersion only"))
lines(c(6.5,6.5), c(0,1), col='darkgrey')
legend("top", fill=col_vector[1:6], ncol=4, bty='n', cex=0.9, 
       legend=c("exp HMM (0.5)", "exp HMM (thr)", "exp HMM (BFDR)", "ln HMM (0.5)", "ln HMM (thr)", 
                "ln HMM (BFDR)"))
colMeans(cbind(DEDD.results.TMM.DD2$tpr, DEDD.results.DESeq.DD2$tpr), na.rm=T)
colMeans(cbind(DEDD.results.TMM.DD5$tpr, DEDD.results.DESeq.DD5$tpr), na.rm=T)
colMeans(cbind(DEDD.results.TMM.DD10$tpr, DEDD.results.DESeq.DD10$tpr), na.rm=T)
colMeans(cbind(DEDD.results.TMM.DD20$tpr, DEDD.results.DESeq.DD20$tpr), na.rm=T)
# All virtually indistinguishable from zero but DESeq norm always higher than TMM; BFDR all zero.

par(mfrow=c(2,2), mar=c(1,2,3,1), mgp=c(3,0.7,0))
boxplot(cbind(DEDD.results.TMM.DE2$tpr, DEDD.results.DESeq.DE2$tpr), 
        names=NA, col=col_vector[1:6], xaxt='n', ylim=c(0,1.2), 
        pch=20, main=paste0("TPR for differential distribution, 2 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in mean only"))
abline(h=seq(0,0.05,0.05), col='lightgrey'); lines(c(6.5,6.5), c(0,1), col='darkgrey')
legend("top", fill=col_vector[1:6], ncol=4, bty='n', cex=0.9, 
       legend=c("exp HMM (0.5)", "exp HMM (thr)", "exp HMM (BFDR)", "ln HMM (0.5)", "ln HMM (thr)", 
                "ln HMM (BFDR)"))
boxplot(cbind(DEDD.results.TMM.DE5$tpr, DEDD.results.DESeq.DE5$tpr), 
        names=NA, col=col_vector[1:6], xaxt='n', ylim=c(0,1.2), 
        pch=20, main=paste0("TPR for differential distribution, 5 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in mean only"))
abline(h=seq(0,0.2,0.05), col='lightgrey'); lines(c(6.5,6.5), c(0,1), col='darkgrey')
legend("top", fill=col_vector[1:6], ncol=4, bty='n', cex=0.9, 
       legend=c("exp HMM (0.5)", "exp HMM (thr)", "exp HMM (BFDR)", "ln HMM (0.5)", "ln HMM (thr)", 
                "ln HMM (BFDR)"))
boxplot(cbind(DEDD.results.TMM.DE10$tpr, DEDD.results.DESeq.DE10$tpr), 
        names=NA, col=col_vector[1:6], xaxt='n', ylim=c(0,1.2), 
        pch=20, main=paste0("TPR for differential distribution, 10 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in mean only"))
abline(h=seq(0.15,0.5,0.05), col='lightgrey'); lines(c(6.5,6.5), c(0,1), col='darkgrey')
legend("top", fill=col_vector[1:6], ncol=4, bty='n', cex=0.9, 
       legend=c("exp HMM (0.5)", "exp HMM (thr)", "exp HMM (BFDR)", "ln HMM (0.5)", "ln HMM (thr)", 
                "ln HMM (BFDR)"))
boxplot(cbind(DEDD.results.TMM.DE20$tpr, DEDD.results.DESeq.DE20$tpr), 
        names=NA, col=col_vector[1:6], xaxt='n', ylim=c(0,1.2), 
        pch=20, main=paste0("TPR for differential distribution, 20 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in mean only"))
abline(h=seq(0.45,0.75,0.05), col='lightgrey'); lines(c(6.5,6.5), c(0,1), col='darkgrey')
legend("top", fill=col_vector[1:6], ncol=4, bty='n', cex=0.9, 
       legend=c("exp HMM (0.5)", "exp HMM (thr)", "exp HMM (BFDR)", "ln HMM (0.5)", "ln HMM (thr)", 
                "ln HMM (BFDR)"))
colMeans(cbind(DEDD.results.TMM.DE2$tpr, DEDD.results.DESeq.DE2$tpr), na.rm=T)
colMeans(cbind(DEDD.results.TMM.DE5$tpr, DEDD.results.DESeq.DE5$tpr), na.rm=T)
colMeans(cbind(DEDD.results.TMM.DE10$tpr, DEDD.results.DESeq.DE10$tpr), na.rm=T)
colMeans(cbind(DEDD.results.TMM.DE20$tpr, DEDD.results.DESeq.DE20$tpr), na.rm=T)
# DESeq norm consistently better than TMM.
# Posterior threshold consistently better than 0.5. BFDR all zero.
# lnHM consistently better than expHM.
# Very low power for 2 or 5 samples, reasonable for 10, good for 20.

par(mfrow=c(2,2), mar=c(1,2,3,1), mgp=c(3,0.7,0))
boxplot(cbind(DEDD.results.TMM.DEDD2$tpr, DEDD.results.DESeq.DEDD2$tpr), 
        names=NA, col=col_vector[1:6], xaxt='n', ylim=c(0,1.2), 
        pch=20, main=paste0("TPR for differential distribution, 2 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in mean and dispersion"))
abline(h=seq(0,0.02,0.02), col='lightgrey'); lines(c(6.5,6.5), c(0,1), col='darkgrey')
legend("top", fill=col_vector[1:6], ncol=4, bty='n', cex=0.9, 
       legend=c("exp HMM (0.5)", "exp HMM (thr)", "exp HMM (BFDR)", "ln HMM (0.5)", "ln HMM (thr)", 
                "ln HMM (BFDR)"))
boxplot(cbind(DEDD.results.TMM.DEDD5$tpr, DEDD.results.DESeq.DEDD5$tpr), 
        names=NA, col=col_vector[1:6], xaxt='n', ylim=c(0,1.2), 
        pch=20, main=paste0("TPR for differential distribution, 5 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in mean and dispersion"))
abline(h=seq(0,0.15,0.05), col='lightgrey'); lines(c(6.5,6.5), c(0,1), col='darkgrey')
legend("top", fill=col_vector[1:6], ncol=4, bty='n', cex=0.9, 
       legend=c("exp HMM (0.5)", "exp HMM (thr)", "exp HMM (BFDR)", "ln HMM (0.5)", "ln HMM (thr)", 
                "ln HMM (BFDR)"))
boxplot(cbind(DEDD.results.TMM.DEDD10$tpr, DEDD.results.DESeq.DEDD10$tpr), 
        names=NA, col=col_vector[1:6], xaxt='n', ylim=c(0,1.2), 
        pch=20, main=paste0("TPR for differential distribution, 10 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in mean and dispersion"))
abline(h=seq(0.1,0.3,0.05), col='lightgrey'); lines(c(6.5,6.5), c(0,1), col='darkgrey')
legend("top", fill=col_vector[1:6], ncol=4, bty='n', cex=0.9, 
       legend=c("exp HMM (0.5)", "exp HMM (thr)", "exp HMM (BFDR)", "ln HMM (0.5)", "ln HMM (thr)", 
                "ln HMM (BFDR)"))
boxplot(cbind(DEDD.results.TMM.DEDD20$tpr, DEDD.results.DESeq.DEDD20$tpr), 
        names=NA, col=col_vector[1:6], xaxt='n', ylim=c(0,1.2), 
        pch=20, main=paste0("TPR for differential distribution, 20 samples per group, TMM (left)\n", 
                            "or DESeq (right) normalisation, differences in mean and dispersion"))
abline(h=seq(0.35,0.55,0.05), col='lightgrey'); lines(c(6.5,6.5), c(0,1), col='darkgrey')
legend("top", fill=col_vector[1:6], ncol=4, bty='n', cex=0.9, 
       legend=c("exp HMM (0.5)", "exp HMM (thr)", "exp HMM (BFDR)", "ln HMM (0.5)", "ln HMM (thr)", 
                "ln HMM (BFDR)"))
colMeans(cbind(DEDD.results.TMM.DEDD2$tpr, DEDD.results.DESeq.DEDD2$tpr), na.rm=T)
colMeans(cbind(DEDD.results.TMM.DEDD5$tpr, DEDD.results.DESeq.DEDD5$tpr), na.rm=T)
colMeans(cbind(DEDD.results.TMM.DEDD10$tpr, DEDD.results.DESeq.DEDD10$tpr), na.rm=T)
colMeans(cbind(DEDD.results.TMM.DEDD20$tpr, DEDD.results.DESeq.DEDD20$tpr), na.rm=T)
# DESeq norm consistently better than TMM.
# Posterior threshold consistently better than 0.5. BFDR all zero.
# lnHM consistently better than expHM.
# Very low power for 2 or 5 samples, reasonable for 10 and 20.


###############################################################
#### Differential distribution combining DE and DD results ####
###############################################################

###########
### AUC ###
par(mfrow=c(2,2), mar=c(1,2,3,1), mgp=c(3,0.7,0))
boxplot(combined_DEDD.results.DEDD2$auc, 
        names=NA, col=col_vector[1:2], xaxt='n', ylim=c(0.45,1.05), 
        pch=20, main=paste0("AUC for differential distribution, 2 samples per group,\n", 
                            "mixture model v combined DE & differential dispersion"))
abline(h=seq(0.60, 0.74, 0.02), col='lightgrey')
legend("top", fill=col_vector[1:6], ncol=3, bty='n', cex=1.2, 
       legend=c("exp HMM", "ln HMM", "lnHM_lnHM", "expHM_expHM", "MDSeq_MDSeq", "HM_edgeR"))
boxplot(combined_DEDD.results.DEDD5$auc, 
        names=NA, col=col_vector[1:2], xaxt='n', ylim=c(0.45,1.05), 
        pch=20, main=paste0("AUC for differential distribution, 5 samples per group,\n", 
                            "mixture model v combined DE & differential dispersion"))
abline(h=seq(0.70, 0.82, 0.02), col='lightgrey')
legend("top", fill=col_vector[1:6], ncol=3, bty='n', cex=1.2, 
       legend=c("exp HMM", "ln HMM", "lnHM_lnHM", "expHM_expHM", "MDSeq_MDSeq", "HM_edgeR"))
boxplot(combined_DEDD.results.DEDD10$auc, 
        names=NA, col=col_vector[1:2], xaxt='n', ylim=c(0.45,1.05), 
        pch=20, main=paste0("AUC for differential distribution, 10 samples per group,\n", 
                            "mixture model v combined DE & differential dispersion"))
abline(h=seq(0.78, 0.86, 0.02), col='lightgrey')
legend("top", fill=col_vector[1:6], ncol=3, bty='n', cex=1.2, 
       legend=c("exp HMM", "ln HMM", "lnHM_lnHM", "expHM_expHM", "MDSeq_MDSeq", "HM_edgeR"))
boxplot(combined_DEDD.results.DEDD20$auc, 
        names=NA, col=col_vector[1:2], xaxt='n', ylim=c(0.45,1.05), 
        pch=20, main=paste0("AUC for differential distribution, 20 samples per group,\n", 
                            "mixture model v combined DE & differential dispersion"))
abline(h=seq(0.84, 0.92, 0.02), col='lightgrey')
legend("top", fill=col_vector[1:6], ncol=3, bty='n', cex=1.2, 
       legend=c("exp HMM", "ln HMM", "lnHM_lnHM", "expHM_expHM", "MDSeq_MDSeq", "HM_edgeR"))
colMeans(combined_DEDD.results.DEDD2$auc)
colMeans(combined_DEDD.results.DEDD5$auc)
colMeans(combined_DEDD.results.DEDD10$auc)
colMeans(combined_DEDD.results.DEDD20$auc)
# HMMs always worse than combined methods. MDSeq always worse than other combined methods, 
# but close for 20 samples per group. HMs very similar to HM/edgeR, but interestingly, 
# lnHM for DE and diff disp usually better on average than lnHM for diff disp with edgeR Ql 
# for DE. expHM very slightly worse than lnHM and edgeR combinations for 2 and 5 samples per 
# groups, but very slightly better for 10 and 20.


############
### pAUC ###
par(mfrow=c(2,2), mar=c(1,2,3,1), mgp=c(3,0.7,0))
boxplot(combined_DEDD.results.DEDD2$pauc, 
        names=NA, col=col_vector[1:2], xaxt='n', ylim=c(0,0.05), 
        pch=20, main=paste0("pAUC for differential distribution, 2 samples per group,\n", 
                            "mixture model v combined DE & differential dispersion"))
abline(h=seq(0.005,0.011,0.001), col='lightgrey')
legend("top", fill=col_vector[1:6], ncol=3, bty='n', cex=1.2, 
       legend=c("exp HMM", "ln HMM", "lnHM_lnHM", "expHM_expHM", "MDSeq_MDSeq", "HM_edgeR"))
boxplot(combined_DEDD.results.DEDD5$pauc, 
        names=NA, col=col_vector[1:2], xaxt='n', ylim=c(0,0.05), 
        pch=20, main=paste0("pAUC for differential distribution, 5 samples per group,\n", 
                            "mixture model v combined DE & differential dispersion"))
abline(h=seq(0.014,0.022,0.002), col='lightgrey')
legend("top", fill=col_vector[1:6], ncol=3, bty='n', cex=1.2, 
       legend=c("exp HMM", "ln HMM", "lnHM_lnHM", "expHM_expHM", "MDSeq_MDSeq", "HM_edgeR"))
boxplot(combined_DEDD.results.DEDD10$pauc, 
        names=NA, col=col_vector[1:2], xaxt='n', ylim=c(0,0.05), 
        pch=20, main=paste0("pAUC for differential distribution, 10 samples per group,\n", 
                            "mixture model v combined DE & differential dispersion"))
abline(h=seq(0.022,0.028,0.002), col='lightgrey')
legend("top", fill=col_vector[1:6], ncol=3, bty='n', cex=1.2, 
       legend=c("exp HMM", "ln HMM", "lnHM_lnHM", "expHM_expHM", "MDSeq_MDSeq", "HM_edgeR"))
boxplot(combined_DEDD.results.DEDD20$pauc, 
        names=NA, col=col_vector[1:2], xaxt='n', ylim=c(0,0.05), 
        pch=20, main=paste0("pAUC for differential distribution, 20 samples per group,\n", 
                            "mixture model v combined DE & differential dispersion"))
abline(h=seq(0.030,0.034,0.002), col='lightgrey')
legend("top", fill=col_vector[1:6], ncol=3, bty='n', cex=1.2, 
       legend=c("exp HMM", "ln HMM", "lnHM_lnHM", "expHM_expHM", "MDSeq_MDSeq", "HM_edgeR"))
colMeans(combined_DEDD.results.DEDD2$pauc)
colMeans(combined_DEDD.results.DEDD5$pauc)
colMeans(combined_DEDD.results.DEDD10$pauc)
colMeans(combined_DEDD.results.DEDD20$pauc)
# HMMs always lowest, followed by MDSeq, as for AUC, but unlike AUC, lnHM is consistently 
# higher than expHM, and lnHM combined with edgeR is consistently the highest.


###########
### FDR ###
par(mfrow=c(2,2), mar=c(1,2,3,1), mgp=c(3,0.7,0))
boxplot(combined_DEDD.results.DEDD2$fdr, 
        names=NA, col=col_vector[1:6], xaxt='n', ylim=c(0,1.2), 
        pch=20, main=paste0("FDR for differential distribution, 2 samples per group,\n", 
                            "mixture model v combined DE & differential dispersion"))
abline(h=0.05, col='lightgrey')
legend("top", fill=col_vector[1:6], ncol=3, bty='n', cex=1.2, 
       legend=c("exp HMM", "ln HMM", "lnHM_lnHM", "expHM_expHM", "MDSeq_MDSeq", "HM_edgeR"))
boxplot(combined_DEDD.results.DEDD5$fdr, 
        names=NA, col=col_vector[1:6], xaxt='n', ylim=c(0,1.2), 
        pch=20, main=paste0("FDR for differential distribution, 5 samples per group,\n", 
                            "mixture model v combined DE & differential dispersion"))
abline(h=0.05, col='lightgrey')
legend("top", fill=col_vector[1:6], ncol=3, bty='n', cex=1.2, 
       legend=c("exp HMM", "ln HMM", "lnHM_lnHM", "expHM_expHM", "MDSeq_MDSeq", "HM_edgeR"))
boxplot(combined_DEDD.results.DEDD10$fdr, 
        names=NA, col=col_vector[1:6], xaxt='n', ylim=c(0,1.2), 
        pch=20, main=paste0("FDR for differential distribution, 10 samples per group,\n", 
                            "mixture model v combined DE & differential dispersion"))
abline(h=0.05, col='lightgrey')
legend("top", fill=col_vector[1:6], ncol=3, bty='n', cex=1.2, 
       legend=c("exp HMM", "ln HMM", "lnHM_lnHM", "expHM_expHM", "MDSeq_MDSeq", "HM_edgeR"))
boxplot(combined_DEDD.results.DEDD20$fdr, 
        names=NA, col=col_vector[1:6], xaxt='n', ylim=c(0,1.2), 
        pch=20, main=paste0("FDR for differential distribution, 20 samples per group,\n", 
                            "mixture model v combined DE & differential dispersion"))
abline(h=0.05, col='lightgrey')
legend("top", fill=col_vector[1:6], ncol=3, bty='n', cex=1.2, 
       legend=c("exp HMM", "ln HMM", "lnHM_lnHM", "expHM_expHM", "MDSeq_MDSeq", "HM_edgeR"))
colMeans(combined_DEDD.results.DEDD2$fdr, na.rm=T)
colMeans(combined_DEDD.results.DEDD5$fdr, na.rm=T)
colMeans(combined_DEDD.results.DEDD10$fdr, na.rm=T)
colMeans(combined_DEDD.results.DEDD20$fdr, na.rm=T)
# No discoviers for any methods except MDSeq for 2 samples per group. For other sample 
# sizes, extremely low FDRs for HMMs and getting lower with bigger samples. Reasonable FDRs 
# for combined methods other than MDSeq but seeming to increase with sample size - below 0.05 
# for 5 and 10 samples per group, above for 20. lnHM combined with edgeR always closer to 0.05
# than lnHM alone and expHM - i.e. bigger when below 0.05 and smaller when above 0.05.


###########
### TPR ###
par(mfrow=c(2,2), mar=c(1,2,3,1), mgp=c(3,0.7,0))
boxplot(combined_DEDD.results.DEDD2$tpr, 
        names=NA, col=col_vector[1:6], xaxt='n', ylim=c(0,1.2), 
        pch=20, main=paste0("TPR for differential distribution, 2 samples per group,\n", 
                            "mixture model v combined DE & differential dispersion"))
abline(h=0, col='lightgrey')
legend("top", fill=col_vector[1:6], ncol=3, bty='n', cex=1.2, 
       legend=c("lnHM_lnHM", "expHM_expHM", "MDSeq_MDSeq", "HM_edgeR", "exp HMM", "ln HMM"))
boxplot(combined_DEDD.results.DEDD5$tpr, 
        names=NA, col=col_vector[1:6], xaxt='n', ylim=c(0,1.2), 
        pch=20, main=paste0("TPR for differential distribution, 5 samples per group,\n", 
                            "mixture model v combined DE & differential dispersion"))
abline(h=seq(0,0.50,0.05), col='lightgrey')
legend("top", fill=col_vector[1:6], ncol=3, bty='n', cex=1.2, 
       legend=c("lnHM_lnHM", "expHM_expHM", "MDSeq_MDSeq", "HM_edgeR", "exp HMM", "ln HMM"))
boxplot(combined_DEDD.results.DEDD10$tpr, 
        names=NA, col=col_vector[1:6], xaxt='n', ylim=c(0,1.2), 
        pch=20, main=paste0("TPR for differential distribution, 10 samples per group,\n", 
                            "mixture model v combined DE & differential dispersion"))
abline(h=seq(0.20,0.60,0.05), col='lightgrey')
legend("top", fill=col_vector[1:6], ncol=3, bty='n', cex=1.2, 
       legend=c("lnHM_lnHM", "expHM_expHM", "MDSeq_MDSeq", "HM_edgeR", "exp HMM", "ln HMM"))
boxplot(combined_DEDD.results.DEDD20$tpr, 
        names=NA, col=col_vector[1:6], xaxt='n', ylim=c(0,1.2), 
        pch=20, main=paste0("TPR for differential distribution, 20 samples per group,\n", 
                            "mixture model v combined DE & differential dispersion"))
abline(h=seq(0.40,0.70,0.05), col='lightgrey')
legend("top", fill=col_vector[1:6], ncol=3, bty='n', cex=1.2, 
       legend=c("lnHM_lnHM", "expHM_expHM", "MDSeq_MDSeq", "HM_edgeR", "exp HMM", "ln HMM"))
colMeans(combined_DEDD.results.DEDD2$tpr, na.rm=T)
colMeans(combined_DEDD.results.DEDD5$tpr, na.rm=T)
colMeans(combined_DEDD.results.DEDD10$tpr, na.rm=T)
colMeans(combined_DEDD.results.DEDD20$tpr, na.rm=T)
# Highest power for MDSeq, but since it always has unacceptable FDRs, shouldn't be considered. 
# Among other methods, HMMs far lower than combined methods, and lnHM combined with edgeR 
# generally highest power as well as having FDRs closest to 0.05; slightly lower power than 
# HMs for 20 samples per group, but only by a very small amount, and with better FDRs. lnHM 
# has higher power than expHM for 5 and 10 samples per group despite having lower FDRs, and 
# only very slightly lower power for 20 samples per group with better FDRs. Overall, between 
# FDRs and TPRs, lnHM-edgeR clearly best, followed by lnHM then expHM.

