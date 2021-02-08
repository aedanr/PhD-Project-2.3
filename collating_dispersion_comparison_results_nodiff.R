library(here)
library(compcodeR)

for (i in c('.nodiff2','.nodiff5','.nodiff10', '.nodiff20')) {
  for (j in c('raw.','norm.')) {
    for (k in c('edgeR.trend', 'edgeR.tag', 'DESeq2', 'DSS.notrend', 'DSS.trend')) {
      assign(paste0('mse.', j, k, i), numeric(50))
    }
    for (l in c('lnHM', 'expHM')) {
      assign(paste0('mse.', j, l, i, '.overall_est'), numeric(50))
      assign(paste0('mse.', j, l, i, '.group1_est'), numeric(50))
      assign(paste0('mse.', j, l, i, '.group2_est'), numeric(50))
    }
  }
}
rm(i,j,k,l)

for (n in 1:50) {
  for (i in c('.nodiff2','.nodiff5','.nodiff10', '.nodiff20')) {
    res <- readRDS(here('Results/Dispersion estimation Dec 2020', 
                        paste0('disp.results', i, '.', n, '.rds')))
    true.disps <- res$data@variable.annotations$truedispersions.S1
    edgeR.trend <- res$disps.edgeR.trended
    edgeR.tag <- res$disps.edgeR.tagwise
    DESeq2 <- res$disps.DESeq
    DSS.notrend <- res$disps.DSS.notrend
    DSS.trend <- res$disps.DSS.trend
    expHM.group1_est <- res$disps.expHM.1
    expHM.group2_est <- res$disps.expHM.2
    expHM.overall_est <- res$disps.expHM
    lnHM.group1_est <- res$disps.lnHM.1
    lnHM.group2_est <- res$disps.lnHM.2
    lnHM.overall_est <- res$disps.lnHM
    for (k in c('edgeR.trend', 'edgeR.tag', 'DESeq2', 'DSS.notrend', 'DSS.trend')) {
        assign(paste0('mse.raw.', k, i), `[<-`(get(paste0('mse.raw.', k, i)), n, 
                                               value=mean((get(k) - true.disps)^2)))
        assign(paste0('mse.norm.', k, i), `[<-`(get(paste0('mse.norm.', k, i)), n, 
                                                value=mean(((get(k) - true.disps) / true.disps)^2)))
      }
      for (l in c('lnHM', 'expHM')) {
        assign(paste0('mse.raw.', l, i, '.overall_est'), `[<-`(get(paste0('mse.raw.', l, i, '.overall_est')), n, 
                     value=mean((get(paste0(l, '.overall_est')) - true.disps)^2)))
        assign(paste0('mse.raw.', l, i, '.group1_est'), `[<-`(get(paste0('mse.raw.', l, i, '.group1_est')), n, 
                     value=mean((get(paste0(l, '.group1_est')) - true.disps)^2)))
        assign(paste0('mse.raw.', l, i, '.group2_est'), `[<-`(get(paste0('mse.raw.', l, i, '.group2_est')), n, 
                     value=mean((get(paste0(l, '.group2_est')) - true.disps)^2)))
        assign(paste0('mse.norm.', l, i, '.overall_est'), `[<-`(get(paste0('mse.norm.', l, i, '.overall_est')), n, 
                     value=mean(((get(paste0(l, '.overall_est')) - true.disps) / 
                                   true.disps)^2)))
        assign(paste0('mse.norm.', l, i, '.group1_est'), `[<-`(get(paste0('mse.norm.', l, i, '.group1_est')), n, 
                     value=mean(((get(paste0(l, '.group1_est')) - true.disps) / 
                                   true.disps)^2)))
        assign(paste0('mse.norm.', l, i, '.group2_est'), `[<-`(get(paste0('mse.norm.', l, i, '.group2_est')), n, 
                     value=mean(((get(paste0(l, '.group2_est')) - true.disps) / 
                                   true.disps)^2)))
      }
  }
}
rm(res,i,k,l,n,true.disps)

nodiff2 <- list(
  raw.MSEs = data.frame(
    edgeR.tag = mse.raw.edgeR.tag.nodiff2, 
    edgeR.trend = mse.raw.edgeR.trend.nodiff2, 
    DESeq2 = mse.raw.DESeq2.nodiff2, 
    DSS.trend = mse.raw.DSS.trend.nodiff2,
    DSS.notrend = mse.raw.DSS.notrend.nodiff2, 
    expHM.1 =  mse.raw.expHM.nodiff2.group1_est, 
    expHM.2 =  mse.raw.expHM.nodiff2.group2_est, 
    expHM.overall = mse.raw.expHM.nodiff2.overall_est, 
    lnHM.1 =  mse.raw.lnHM.nodiff2.group1_est, 
    lnHM.2 =  mse.raw.lnHM.nodiff2.group2_est, 
    lnHM.overall = mse.raw.lnHM.nodiff2.overall_est
  ), 
  normalised.MSEs = data.frame(
    edgeR.tag = mse.norm.edgeR.tag.nodiff2, 
    edgeR.trend = mse.norm.edgeR.trend.nodiff2, 
    DESeq2 = mse.norm.DESeq2.nodiff2, 
    DSS.trend = mse.norm.DSS.trend.nodiff2,
    DSS.notrend = mse.norm.DSS.notrend.nodiff2, 
    expHM.1 =  mse.norm.expHM.nodiff2.group1_est, 
    expHM.2 =  mse.norm.expHM.nodiff2.group2_est, 
    expHM.overall = mse.norm.expHM.nodiff2.overall_est, 
    lnHM.1 =  mse.norm.lnHM.nodiff2.group1_est, 
    lnHM.2 =  mse.norm.lnHM.nodiff2.group2_est, 
    lnHM.overall = mse.norm.lnHM.nodiff2.overall_est
  )
)

nodiff5 <- list(
  raw.MSEs = data.frame(
    edgeR.tag = mse.raw.edgeR.tag.nodiff5, 
    edgeR.trend = mse.raw.edgeR.trend.nodiff5, 
    DESeq2 = mse.raw.DESeq2.nodiff5, 
    DSS.trend = mse.raw.DSS.trend.nodiff5,
    DSS.notrend = mse.raw.DSS.notrend.nodiff5, 
    expHM.1 =  mse.raw.expHM.nodiff5.group1_est, 
    expHM.2 =  mse.raw.expHM.nodiff5.group2_est, 
    expHM.overall = mse.raw.expHM.nodiff5.overall_est, 
    lnHM.1 =  mse.raw.lnHM.nodiff5.group1_est, 
    lnHM.2 =  mse.raw.lnHM.nodiff5.group2_est, 
    lnHM.overall = mse.raw.lnHM.nodiff5.overall_est
  ), 
  normalised.MSEs = data.frame(
    edgeR.tag = mse.norm.edgeR.tag.nodiff5, 
    edgeR.trend = mse.norm.edgeR.trend.nodiff5, 
    DESeq2 = mse.norm.DESeq2.nodiff5, 
    DSS.trend = mse.norm.DSS.trend.nodiff5,
    DSS.notrend = mse.norm.DSS.notrend.nodiff5, 
    expHM.1 =  mse.norm.expHM.nodiff5.group1_est, 
    expHM.2 =  mse.norm.expHM.nodiff5.group2_est, 
    expHM.overall = mse.norm.expHM.nodiff5.overall_est, 
    lnHM.1 =  mse.norm.lnHM.nodiff5.group1_est, 
    lnHM.2 =  mse.norm.lnHM.nodiff5.group2_est, 
    lnHM.overall = mse.norm.lnHM.nodiff5.overall_est
  )
)

nodiff10 <- list(
  raw.MSEs = data.frame(
    edgeR.tag = mse.raw.edgeR.tag.nodiff10, 
    edgeR.trend = mse.raw.edgeR.trend.nodiff10, 
    DESeq2 = mse.raw.DESeq2.nodiff10, 
    DSS.trend = mse.raw.DSS.trend.nodiff10,
    DSS.notrend = mse.raw.DSS.notrend.nodiff10, 
    expHM.1 =  mse.raw.expHM.nodiff10.group1_est, 
    expHM.2 =  mse.raw.expHM.nodiff10.group2_est, 
    expHM.overall = mse.raw.expHM.nodiff10.overall_est, 
    lnHM.1 =  mse.raw.lnHM.nodiff10.group1_est, 
    lnHM.2 =  mse.raw.lnHM.nodiff10.group2_est, 
    lnHM.overall = mse.raw.lnHM.nodiff10.overall_est
  ), 
  normalised.MSEs = data.frame(
    edgeR.tag = mse.norm.edgeR.tag.nodiff10, 
    edgeR.trend = mse.norm.edgeR.trend.nodiff10, 
    DESeq2 = mse.norm.DESeq2.nodiff10, 
    DSS.trend = mse.norm.DSS.trend.nodiff10,
    DSS.notrend = mse.norm.DSS.notrend.nodiff10, 
    expHM.1 =  mse.norm.expHM.nodiff10.group1_est, 
    expHM.2 =  mse.norm.expHM.nodiff10.group2_est, 
    expHM.overall = mse.norm.expHM.nodiff10.overall_est, 
    lnHM.1 =  mse.norm.lnHM.nodiff10.group1_est, 
    lnHM.2 =  mse.norm.lnHM.nodiff10.group2_est, 
    lnHM.overall = mse.norm.lnHM.nodiff10.overall_est
  )
)

nodiff20 <- list(
  raw.MSEs = data.frame(
    edgeR.tag = mse.raw.edgeR.tag.nodiff20, 
    edgeR.trend = mse.raw.edgeR.trend.nodiff20, 
    DESeq2 = mse.raw.DESeq2.nodiff20, 
    DSS.trend = mse.raw.DSS.trend.nodiff20,
    DSS.notrend = mse.raw.DSS.notrend.nodiff20, 
    expHM.1 =  mse.raw.expHM.nodiff20.group1_est, 
    expHM.2 =  mse.raw.expHM.nodiff20.group2_est, 
    expHM.overall = mse.raw.expHM.nodiff20.overall_est, 
    lnHM.1 =  mse.raw.lnHM.nodiff20.group1_est, 
    lnHM.2 =  mse.raw.lnHM.nodiff20.group2_est, 
    lnHM.overall = mse.raw.lnHM.nodiff20.overall_est
  ), 
  normalised.MSEs = data.frame(
    edgeR.tag = mse.norm.edgeR.tag.nodiff20, 
    edgeR.trend = mse.norm.edgeR.trend.nodiff20, 
    DESeq2 = mse.norm.DESeq2.nodiff20, 
    DSS.trend = mse.norm.DSS.trend.nodiff20,
    DSS.notrend = mse.norm.DSS.notrend.nodiff20, 
    expHM.1 =  mse.norm.expHM.nodiff20.group1_est, 
    expHM.2 =  mse.norm.expHM.nodiff20.group2_est, 
    expHM.overall = mse.norm.expHM.nodiff20.overall_est, 
    lnHM.1 =  mse.norm.lnHM.nodiff20.group1_est, 
    lnHM.2 =  mse.norm.lnHM.nodiff20.group2_est, 
    lnHM.overall = mse.norm.lnHM.nodiff20.overall_est
  )
)

saveRDS(nodiff2, file=here('Results/Dispersion estimation Dec 2020','mse.disp.nodiff2.rds'))
saveRDS(nodiff5, file=here('Results/Dispersion estimation Dec 2020','mse.disp.nodiff5.rds'))
saveRDS(nodiff10, file=here('Results/Dispersion estimation Dec 2020','mse.disp.nodiff10.rds'))
saveRDS(nodiff20, file=here('Results/Dispersion estimation Dec 2020','mse.disp.nodiff20.rds'))





### below notes are from DEDD dispersion estimation with TMM (from original analysis Oct 2019) ###
### seem to apply similarly to DE with TMM ###
between <- function(x, min, max) {
  return(x > min & x <= max)
}
c(mean(between(true.disps, exp(-15), exp(-3))), 
  mean(between(true.disps, exp(-3), exp(-2.5))), 
  mean(between(true.disps, exp(-2.5), exp(-2))), 
  mean(between(true.disps, exp(-2), exp(-1.5))), 
  mean(between(true.disps, exp(-1.5), exp(-1))), 
  mean(between(true.disps, exp(-1), exp(0))), 
  mean(between(true.disps, exp(0), exp(2))))
g1.1 <- which(true.disps < exp(-3))
g1.2 <- which(between(true.disps, exp(-3), exp(-2.5)))
g1.3 <- which(between(true.disps, exp(-2.5), exp(-2)))
g1.4 <- which(between(true.disps, exp(-2), exp(-1.5)))
g1.5 <- which(between(true.disps, exp(-1.5), exp(-1)))
g1.6 <- which(between(true.disps, exp(-1), exp(0)))
g1.7 <- which(true.disps >= exp(0))

rbind(c(mean((lnHM.overall_est[g1.1] - true.disps[g1.1])^2), 
        mean((lnHM.overall_est[g1.2] - true.disps[g1.2])^2), 
        mean((lnHM.overall_est[g1.3] - true.disps[g1.3])^2), 
        mean((lnHM.overall_est[g1.4] - true.disps[g1.4])^2), 
        mean((lnHM.overall_est[g1.5] - true.disps[g1.5])^2), 
        mean((lnHM.overall_est[g1.6] - true.disps[g1.6])^2), 
        mean((lnHM.overall_est[g1.7] - true.disps[g1.7])^2)), 
      c(mean((edgeR.tag[g1.1] - true.disps[g1.1])^2), 
        mean((edgeR.tag[g1.2] - true.disps[g1.2])^2), 
        mean((edgeR.tag[g1.3] - true.disps[g1.3])^2), 
        mean((edgeR.tag[g1.4] - true.disps[g1.4])^2), 
        mean((edgeR.tag[g1.5] - true.disps[g1.5])^2), 
        mean((edgeR.tag[g1.6] - true.disps[g1.6])^2), 
        mean((edgeR.tag[g1.7] - true.disps[g1.7])^2)))
# DESeq norm: lnHM better than edgeR on most ranges of dispersions, and much better for smaller dispersions
# TMM: edgeR better than lnHM except for extreme highest dispersions

par(mfrow=c(2,2), mar=c(2,2,1,1), mgp=c(2,1,0))
plot(log(true.disps), 
     abs(lnHM.overall_est - true.disps)^2, pch=20, col='red', 
     xlim=c(-5,1), ylim=c(0,0.8))
plot(log(true.disps), 
     abs(edgeR.tag - true.disps)^2, pch=20, col='blue', 
     xlim=c(-5,1), ylim=c(0,0.8))
plot(log(true.disps), 
     abs(DESeq2 - true.disps)^2, pch=20, col='yellow', 
     xlim=c(-5,1), ylim=c(0,0.8))
plot(log(true.disps), 
     abs(lnHM.overall_est - true.disps) - 
       abs(edgeR.tag - true.disps), pch=20)
plot(log(true.disps), 
     abs(lnHM.overall_est - true.disps)^2 - 
       abs(edgeR.tag - true.disps)^2, pch=20)
mean(abs(lnHM.overall_est - true.disps) - 
       abs(edgeR.tag - true.disps))
mean(abs(lnHM.overall_est - true.disps)^2 - 
       abs(edgeR.tag - true.disps)^2)
mean(abs(lnHM.overall_est - true.disps) > 
       abs(edgeR.tag - true.disps))
mean(abs(lnHM.overall_est - true.disps) > 
       abs(DESeq2 - true.disps))
plot(edgeR.tag, lnHM.overall_est, pch=20)
mean(lnHM.overall_est > edgeR.tag)
plot(log(edgeR.tag), log(lnHM.overall_est), pch=20)
mean(abs(lnHM.overall_est - true.disps) > 
         abs(edgeR.tag - true.disps))
mean(abs(lnHM.overall_est - true.disps) > 
         abs(DESeq2 - true.disps))

plot(true.disps, abs(lnHM.overall_est - true.disps) - abs(edgeR.tag - true.disps), pch=20)

rbind(c('<0.1','0.1-0.2','0.2-0.3','0.3-0.4','0.4-0.5','0.5-1','1-2','>2'), 
      round(c(mean(abs(lnHM.overall_est - true.disps)[which(true.disps < 0.1)] > 
                     abs(edgeR.tag - true.disps)[which(true.disps < 0.1)]),
              mean(abs(lnHM.overall_est - true.disps)[which(between(true.disps, 0.1, 0.2))] > 
                     abs(edgeR.tag - true.disps)[which(between(true.disps, 0.1, 0.2))]),
              mean(abs(lnHM.overall_est - true.disps)[which(between(true.disps, 0.2, 0.3))] > 
                     abs(edgeR.tag - true.disps)[which(between(true.disps, 0.2, 0.3))]),
              mean(abs(lnHM.overall_est - true.disps)[which(between(true.disps, 0.3, 0.4))] > 
                     abs(edgeR.tag - true.disps)[which(between(true.disps, 0.3, 0.4))]),
              mean(abs(lnHM.overall_est - true.disps)[which(between(true.disps, 0.4, 0.5))] > 
                     abs(edgeR.tag - true.disps)[which(between(true.disps, 0.4, 0.5))]),
              mean(abs(lnHM.overall_est - true.disps)[which(between(true.disps, 0.5, 1))] > 
                     abs(edgeR.tag - true.disps)[which(between(true.disps, 0.5, 1))]),
              mean(abs(lnHM.overall_est - true.disps)[which(between(true.disps, 1, 2))] > 
                     abs(edgeR.tag - true.disps)[which(between(true.disps, 1, 2))]),
              mean(abs(lnHM.overall_est - true.disps)[which(true.disps > 2)] > 
                     abs(edgeR.tag - true.disps)[which(true.disps > 2)])), 2))
# DESeqnorm: lnHM error > edgeR.tag error around half the time for dispersions > 0.2, almost never for < 0.2
# TMM: lnHM error > edgeR.tag error around half the time for dispersions > 0.2, almost all the time for < 0.2
mean(log(edgeR.tag))
mean(log(lnHM.overall_est))
mean(log(true.disps))      
hist(log(edgeR.tag))
hist(log(lnHM.overall_est))      
plot(density(log(edgeR.tag)), col='blue', ylim=c(0,0.6))
lines(density(log(lnHM.overall_est)), col='red')
abline(v=mean(log(true.disps)))
median(edgeR.tag)
median(lnHM.overall_est)
median(true.disps)
lines(density(log(true.disps)), col='orange')
# DESeqnorm: lnHM density matches true density almost exactly, edgeR.tag ok for high dispersions but squeezed to right, 
# and mean and median for lnHM close to true, but edgeR overestimates
# TMM: edgeR.tag density closely matches true density, lnHM good for high dispersions but slightly squeezed to right, 
# and mean and median for edgeR.tag close to true, but lnHM overestimates


