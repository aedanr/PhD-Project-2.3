library(here)
library(compcodeR)

for (i in c('.DE2','.DE5','.DE10', '.DE20')) {
  for (j in c('raw.','norm.')) {
    for (k in c('edgeR.trend', 'edgeR.tag', 'DESeq2', 'DSS.notrend', 'DSS.trend')) {
      assign(paste0('mse.', j, k, i, '.noDE'), numeric(50))
      assign(paste0('mse.', j, k, i, '.DE'), numeric(50))
      assign(paste0('mse.', j, k, i, '.all_genes'), numeric(50))
    }
    for (l in c('lnHM', 'expHM')) {
      assign(paste0('mse.', j, l, i, '.overall_est.noDE'), numeric(50))
      assign(paste0('mse.', j, l, i, '.group1_est.noDE'), numeric(50))
      assign(paste0('mse.', j, l, i, '.group2_est.noDE'), numeric(50))
      assign(paste0('mse.', j, l, i, '.overall_est.DE'), numeric(50))
      assign(paste0('mse.', j, l, i, '.group1_est.DE'), numeric(50))
      assign(paste0('mse.', j, l, i, '.group2_est.DE'), numeric(50))
      assign(paste0('mse.', j, l, i, '.overall_est.all_genes'), numeric(50))
      assign(paste0('mse.', j, l, i, '.group1_est.all_genes'), numeric(50))
      assign(paste0('mse.', j, l, i, '.group2_est.all_genes'), numeric(50))
    }
  }
}
rm(i,j,k,l)

## MSEs ####
for (n in 1:50) {
  for (i in c('.DE2', '.DE5', '.DE10', '.DE20')) {
    res <- readRDS(here('Results/Dispersion estimation Dec 2020', 
                        paste0('disp.results', i, '.', n, '.rds')))
    DE <- res$data@variable.annotations$differential.expression
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
        assign(paste0('mse.raw.', k, i, '.noDE'), 
               `[<-`(get(paste0('mse.raw.', k, i, '.noDE')), n, 
                     value=mean((get(k) - true.disps)[which(DE == 0)]^2)))
        assign(paste0('mse.raw.', k, i, '.DE'), 
               `[<-`(get(paste0('mse.raw.', k, i, '.DE')), n, 
                     value=mean((get(k) - true.disps)[which(DE == 1)]^2)))
        assign(paste0('mse.raw.', k, i, '.all_genes'), 
               `[<-`(get(paste0('mse.raw.', k, i, '.all_genes')), n, 
                     value=mean((get(k) - true.disps)^2)))
        assign(paste0('mse.norm.', k, i, '.noDE'), 
               `[<-`(get(paste0('mse.norm.', k, i, '.noDE')), n, 
                     value=mean(((get(k) - true.disps) / true.disps)[which(DE == 0)]^2)))
        assign(paste0('mse.norm.', k, i, '.DE'), 
               `[<-`(get(paste0('mse.norm.', k, i, '.DE')), n, 
                     value=mean(((get(k) - true.disps) / true.disps)[which(DE == 1)]^2)))
        assign(paste0('mse.norm.', k, i, '.all_genes'), 
               `[<-`(get(paste0('mse.norm.', k, i, '.all_genes')), n, 
                     value=mean(((get(k) - true.disps) / true.disps)^2)))
      }
      for (l in c('lnHM', 'expHM')) {
        assign(paste0('mse.raw.', l, i, '.overall_est.noDE'), 
               `[<-`(get(paste0('mse.raw.', l, i, '.overall_est.noDE')), n, 
                     value=mean((get(paste0(l, '.overall_est')) - true.disps)[which(DE == 0)]^2)))
        assign(paste0('mse.raw.', l, i, '.group1_est.noDE'), 
               `[<-`(get(paste0('mse.raw.', l, i, '.group1_est.noDE')), n, 
                     value=mean((get(paste0(l, '.group1_est')) - true.disps)[which(DE == 0)]^2)))
        assign(paste0('mse.raw.', l, i, '.group2_est.noDE'), 
               `[<-`(get(paste0('mse.raw.', l, i, '.group2_est.noDE')), n, 
                     value=mean((get(paste0(l, '.group2_est')) - true.disps)[which(DE == 0)]^2)))
        assign(paste0('mse.raw.', l, i, '.overall_est.DE'), 
               `[<-`(get(paste0('mse.raw.', l, i, '.overall_est.DE')), n, 
                     value=mean((get(paste0(l, '.overall_est')) - true.disps)[which(DE == 1)]^2)))
        assign(paste0('mse.raw.', l, i, '.group1_est.DE'), 
               `[<-`(get(paste0('mse.raw.', l, i, '.group1_est.DE')), n, 
                     value=mean((get(paste0(l, '.group1_est')) - true.disps)[which(DE == 1)]^2)))
        assign(paste0('mse.raw.', l, i, '.group2_est.DE'), 
               `[<-`(get(paste0('mse.raw.', l, i, '.group2_est.DE')), n, 
                     value=mean((get(paste0(l, '.group2_est')) - true.disps)[which(DE == 1)]^2)))
        assign(paste0('mse.raw.', l, i, '.overall_est.all_genes'), 
               `[<-`(get(paste0('mse.raw.', l, i, '.overall_est.all_genes')), n, 
                     value=mean((get(paste0(l, '.overall_est')) - true.disps)^2)))
        assign(paste0('mse.raw.', l, i, '.group1_est.all_genes'), 
               `[<-`(get(paste0('mse.raw.', l, i, '.group1_est.all_genes')), n, 
                     value=mean((get(paste0(l, '.group1_est')) - true.disps)^2)))
        assign(paste0('mse.raw.', l, i, '.group2_est.all_genes'), 
               `[<-`(get(paste0('mse.raw.', l, i, '.group2_est.all_genes')), n, 
                     value=mean((get(paste0(l, '.group2_est')) - true.disps)^2)))
        assign(paste0('mse.norm.', l, i, '.overall_est.noDE'), 
               `[<-`(get(paste0('mse.norm.', l, i, '.overall_est.noDE')), n, 
                     value=mean(((get(paste0(l, '.overall_est')) - true.disps) / 
                                   true.disps)[which(DE == 0)]^2)))
        assign(paste0('mse.norm.', l, i, '.group1_est.noDE'), 
               `[<-`(get(paste0('mse.norm.', l, i, '.group1_est.noDE')), n, 
                     value=mean(((get(paste0(l, '.group1_est')) - true.disps) / 
                                   true.disps)[which(DE == 0)]^2)))
        assign(paste0('mse.norm.', l, i, '.group2_est.noDE'), 
               `[<-`(get(paste0('mse.norm.', l, i, '.group2_est.noDE')), n, 
                     value=mean(((get(paste0(l, '.group2_est')) - true.disps) / 
                                   true.disps)[which(DE == 0)]^2)))
        assign(paste0('mse.norm.', l, i, '.overall_est.DE'), 
               `[<-`(get(paste0('mse.norm.', l, i, '.overall_est.DE')), n, 
                     value=mean(((get(paste0(l, '.overall_est')) - true.disps) / 
                                   true.disps)[which(DE == 1)]^2)))
        assign(paste0('mse.norm.', l, i, '.group1_est.DE'), 
               `[<-`(get(paste0('mse.norm.', l, i, '.group1_est.DE')), n, 
                     value=mean(((get(paste0(l, '.group1_est')) - true.disps) / 
                                   true.disps)[which(DE == 1)]^2)))
        assign(paste0('mse.norm.', l, i, '.group2_est.DE'), 
               `[<-`(get(paste0('mse.norm.', l, i, '.group2_est.DE')), n, 
                     value=mean(((get(paste0(l, '.group2_est')) - true.disps) / 
                                   true.disps)[which(DE == 1)]^2)))
        assign(paste0('mse.norm.', l, i, '.overall_est.all_genes'), 
               `[<-`(get(paste0('mse.norm.', l, i, '.overall_est.all_genes')), n, 
                     value=mean(((get(paste0(l, '.overall_est')) - true.disps) / 
                                   true.disps)^2)))
        assign(paste0('mse.norm.', l, i, '.group1_est.all_genes'), 
               `[<-`(get(paste0('mse.norm.', l, i, '.group1_est.all_genes')), n, 
                     value=mean(((get(paste0(l, '.group1_est')) - true.disps) / 
                                   true.disps)^2)))
        assign(paste0('mse.norm.', l, i, '.group2_est.all_genes'), 
               `[<-`(get(paste0('mse.norm.', l, i, '.group2_est.all_genes')), n, 
                     value=mean(((get(paste0(l, '.group2_est')) - true.disps) / 
                                   true.disps)^2)))
      }
  }
}
rm(i,k,l,n,true.disps)

DE2 <- list(
  raw.MSEs = list(
    noDE = data.frame(
      edgeR.tag = mse.raw.edgeR.tag.DE2.noDE, 
      edgeR.trend = mse.raw.edgeR.trend.DE2.noDE, 
      DESeq2 = mse.raw.DESeq2.DE2.noDE, 
      DSS.trend = mse.raw.DSS.trend.DE2.noDE,
      DSS.notrend = mse.raw.DSS.notrend.DE2.noDE, 
      expHM.1 =  mse.raw.expHM.DE2.group1_est.noDE, 
      expHM.2 =  mse.raw.expHM.DE2.group2_est.noDE, 
      expHM.overall = mse.raw.expHM.DE2.overall_est.noDE, 
      lnHM.1 =  mse.raw.lnHM.DE2.group1_est.noDE, 
      lnHM.2 =  mse.raw.lnHM.DE2.group2_est.noDE, 
      lnHM.overall = mse.raw.lnHM.DE2.overall_est.noDE
    ), 
    DE = data.frame(
      edgeR.tag = mse.raw.edgeR.tag.DE2.DE, 
      edgeR.trend = mse.raw.edgeR.trend.DE2.DE, 
      DESeq2 = mse.raw.DESeq2.DE2.DE, 
      DSS.trend = mse.raw.DSS.trend.DE2.DE, 
      DSS.notrend = mse.raw.DSS.notrend.DE2.DE, 
      expHM.1 =  mse.raw.expHM.DE2.group1_est.DE, 
      expHM.2 =  mse.raw.expHM.DE2.group2_est.DE, 
      expHM.overall = mse.raw.expHM.DE2.overall_est.DE, 
      lnHM.1 =  mse.raw.lnHM.DE2.group1_est.DE, 
      lnHM.2 =  mse.raw.lnHM.DE2.group2_est.DE, 
      lnHM.overall = mse.raw.lnHM.DE2.overall_est.DE
    ), 
    all_genes = data.frame(
      edgeR.tag = mse.raw.edgeR.tag.DE2.all_genes, 
      edgeR.trend = mse.raw.edgeR.trend.DE2.all_genes, 
      DESeq2 = mse.raw.DESeq2.DE2.all_genes, 
      DSS.trend = mse.raw.DSS.trend.DE2.all_genes, 
      DSS.notrend = mse.raw.DSS.notrend.DE2.all_genes, 
      expHM.1 =  mse.raw.expHM.DE2.group1_est.all_genes, 
      expHM.2 =  mse.raw.expHM.DE2.group2_est.all_genes, 
      expHM.overall = mse.raw.expHM.DE2.overall_est.all_genes, 
      lnHM.1 =  mse.raw.lnHM.DE2.group1_est.all_genes, 
      lnHM.2 =  mse.raw.lnHM.DE2.group2_est.all_genes, 
      lnHM.overall = mse.raw.lnHM.DE2.overall_est.all_genes
    )
  ), 
  normalised.MSEs = list(
    noDE = data.frame(
      edgeR.tag = mse.norm.edgeR.tag.DE2.noDE, 
      edgeR.trend = mse.norm.edgeR.trend.DE2.noDE, 
      DESeq2 = mse.norm.DESeq2.DE2.noDE, 
      DSS.trend = mse.norm.DSS.trend.DE2.noDE,
      DSS.notrend = mse.norm.DSS.notrend.DE2.noDE, 
      expHM.1 =  mse.norm.expHM.DE2.group1_est.noDE, 
      expHM.2 =  mse.norm.expHM.DE2.group2_est.noDE, 
      expHM.overall = mse.norm.expHM.DE2.overall_est.noDE, 
      lnHM.1 =  mse.norm.lnHM.DE2.group1_est.noDE, 
      lnHM.2 =  mse.norm.lnHM.DE2.group2_est.noDE, 
      lnHM.overall = mse.norm.lnHM.DE2.overall_est.noDE
    ), 
    DE = data.frame(
      edgeR.tag = mse.norm.edgeR.tag.DE2.DE, 
      edgeR.trend = mse.norm.edgeR.trend.DE2.DE, 
      DESeq2 = mse.norm.DESeq2.DE2.DE, 
      DSS.trend = mse.norm.DSS.trend.DE2.DE, 
      DSS.notrend = mse.norm.DSS.notrend.DE2.DE, 
      expHM.1 =  mse.norm.expHM.DE2.group1_est.DE, 
      expHM.2 =  mse.norm.expHM.DE2.group2_est.DE, 
      expHM.overall = mse.norm.expHM.DE2.overall_est.DE, 
      lnHM.1 =  mse.norm.lnHM.DE2.group1_est.DE, 
      lnHM.2 =  mse.norm.lnHM.DE2.group2_est.DE, 
      lnHM.overall = mse.norm.lnHM.DE2.overall_est.DE
    ), 
    all_genes = data.frame(
      edgeR.tag = mse.norm.edgeR.tag.DE2.all_genes, 
      edgeR.trend = mse.norm.edgeR.trend.DE2.all_genes, 
      DESeq2 = mse.norm.DESeq2.DE2.all_genes, 
      DSS.trend = mse.norm.DSS.trend.DE2.all_genes, 
      DSS.notrend = mse.norm.DSS.notrend.DE2.all_genes, 
      expHM.1 =  mse.norm.expHM.DE2.group1_est.all_genes, 
      expHM.2 =  mse.norm.expHM.DE2.group2_est.all_genes, 
      expHM.overall = mse.norm.expHM.DE2.overall_est.all_genes, 
      lnHM.1 =  mse.norm.lnHM.DE2.group1_est.all_genes, 
      lnHM.2 =  mse.norm.lnHM.DE2.group2_est.all_genes, 
      lnHM.overall = mse.norm.lnHM.DE2.overall_est.all_genes
    )
  )
)

DE5 <- list(
  raw.MSEs = list(
    noDE = data.frame(
      edgeR.tag = mse.raw.edgeR.tag.DE5.noDE, 
      edgeR.trend = mse.raw.edgeR.trend.DE5.noDE, 
      DESeq2 = mse.raw.DESeq2.DE5.noDE, 
      DSS.trend = mse.raw.DSS.trend.DE5.noDE, 
      DSS.notrend = mse.raw.DSS.notrend.DE5.noDE, 
      expHM.1 =  mse.raw.expHM.DE5.group1_est.noDE, 
      expHM.2 =  mse.raw.expHM.DE5.group2_est.noDE, 
      expHM.overall = mse.raw.expHM.DE5.overall_est.noDE, 
      lnHM.1 =  mse.raw.lnHM.DE5.group1_est.noDE, 
      lnHM.2 =  mse.raw.lnHM.DE5.group2_est.noDE, 
      lnHM.overall = mse.raw.lnHM.DE5.overall_est.noDE
    ), 
    DE = data.frame(
      edgeR.tag = mse.raw.edgeR.tag.DE5.DE, 
      edgeR.trend = mse.raw.edgeR.trend.DE5.DE, 
      DESeq2 = mse.raw.DESeq2.DE5.DE, 
      DSS.trend = mse.raw.DSS.trend.DE5.DE, 
      DSS.notrend = mse.raw.DSS.notrend.DE5.DE, 
      expHM.1 =  mse.raw.expHM.DE5.group1_est.DE, 
      expHM.2 =  mse.raw.expHM.DE5.group2_est.DE, 
      expHM.overall = mse.raw.expHM.DE5.overall_est.DE, 
      lnHM.1 =  mse.raw.lnHM.DE5.group1_est.DE, 
      lnHM.2 =  mse.raw.lnHM.DE5.group2_est.DE, 
      lnHM.overall = mse.raw.lnHM.DE5.overall_est.DE
    ), 
    all_genes = data.frame(
      edgeR.tag = mse.raw.edgeR.tag.DE5.all_genes, 
      edgeR.trend = mse.raw.edgeR.trend.DE5.all_genes, 
      DESeq2 = mse.raw.DESeq2.DE5.all_genes, 
      DSS.trend = mse.raw.DSS.trend.DE5.all_genes, 
      DSS.notrend = mse.raw.DSS.notrend.DE5.all_genes, 
      expHM.1 =  mse.raw.expHM.DE5.group1_est.all_genes, 
      expHM.2 =  mse.raw.expHM.DE5.group2_est.all_genes, 
      expHM.overall = mse.raw.expHM.DE5.overall_est.all_genes, 
      lnHM.1 =  mse.raw.lnHM.DE5.group1_est.all_genes, 
      lnHM.2 =  mse.raw.lnHM.DE5.group2_est.all_genes, 
      lnHM.overall = mse.raw.lnHM.DE5.overall_est.all_genes
    )
  ), 
  normalised.MSEs = list(
    noDE = data.frame(
      edgeR.tag = mse.norm.edgeR.tag.DE5.noDE, 
      edgeR.trend = mse.norm.edgeR.trend.DE5.noDE, 
      DESeq2 = mse.norm.DESeq2.DE5.noDE, 
      DSS.trend = mse.norm.DSS.trend.DE5.noDE, 
      DSS.notrend = mse.norm.DSS.notrend.DE5.noDE, 
      expHM.1 =  mse.norm.expHM.DE5.group1_est.noDE, 
      expHM.2 =  mse.norm.expHM.DE5.group2_est.noDE, 
      expHM.overall = mse.norm.expHM.DE5.overall_est.noDE, 
      lnHM.1 =  mse.norm.lnHM.DE5.group1_est.noDE, 
      lnHM.2 =  mse.norm.lnHM.DE5.group2_est.noDE, 
      lnHM.overall = mse.norm.lnHM.DE5.overall_est.noDE
    ), 
    DE = data.frame(
      edgeR.tag = mse.norm.edgeR.tag.DE5.DE, 
      edgeR.trend = mse.norm.edgeR.trend.DE5.DE, 
      DESeq2 = mse.norm.DESeq2.DE5.DE, 
      DSS.trend = mse.norm.DSS.trend.DE5.DE, 
      DSS.notrend = mse.norm.DSS.notrend.DE5.DE, 
      expHM.1 =  mse.norm.expHM.DE5.group1_est.DE, 
      expHM.2 =  mse.norm.expHM.DE5.group2_est.DE, 
      expHM.overall = mse.norm.expHM.DE5.overall_est.DE, 
      lnHM.1 =  mse.norm.lnHM.DE5.group1_est.DE, 
      lnHM.2 =  mse.norm.lnHM.DE5.group2_est.DE, 
      lnHM.overall = mse.norm.lnHM.DE5.overall_est.DE
    ), 
    all_genes = data.frame(
      edgeR.tag = mse.norm.edgeR.tag.DE5.all_genes, 
      edgeR.trend = mse.norm.edgeR.trend.DE5.all_genes, 
      DESeq2 = mse.norm.DESeq2.DE5.all_genes, 
      DSS.trend = mse.norm.DSS.trend.DE5.all_genes, 
      DSS.notrend = mse.norm.DSS.notrend.DE5.all_genes, 
      expHM.1 =  mse.norm.expHM.DE5.group1_est.all_genes, 
      expHM.2 =  mse.norm.expHM.DE5.group2_est.all_genes, 
      expHM.overall = mse.norm.expHM.DE5.overall_est.all_genes, 
      lnHM.1 =  mse.norm.lnHM.DE5.group1_est.all_genes, 
      lnHM.2 =  mse.norm.lnHM.DE5.group2_est.all_genes, 
      lnHM.overall = mse.norm.lnHM.DE5.overall_est.all_genes
    )
  )
)

DE10 <- list(
  raw.MSEs = list(
    noDE = data.frame(
      edgeR.tag = mse.raw.edgeR.tag.DE10.noDE, 
      edgeR.trend = mse.raw.edgeR.trend.DE10.noDE, 
      DESeq2 = mse.raw.DESeq2.DE10.noDE, 
      DSS.trend = mse.raw.DSS.trend.DE10.noDE, 
      DSS.notrend = mse.raw.DSS.notrend.DE10.noDE, 
      expHM.1 =  mse.raw.expHM.DE10.group1_est.noDE, 
      expHM.2 =  mse.raw.expHM.DE10.group2_est.noDE, 
      expHM.overall = mse.raw.expHM.DE10.overall_est.noDE, 
      lnHM.1 =  mse.raw.lnHM.DE10.group1_est.noDE, 
      lnHM.2 =  mse.raw.lnHM.DE10.group2_est.noDE, 
      lnHM.overall = mse.raw.lnHM.DE10.overall_est.noDE
    ), 
    DE = data.frame(
      edgeR.tag = mse.raw.edgeR.tag.DE10.DE, 
      edgeR.trend = mse.raw.edgeR.trend.DE10.DE, 
      DESeq2 = mse.raw.DESeq2.DE10.DE, 
      DSS.trend = mse.raw.DSS.trend.DE10.DE, 
      DSS.notrend = mse.raw.DSS.notrend.DE10.DE, 
      expHM.1 =  mse.raw.expHM.DE10.group1_est.DE, 
      expHM.2 =  mse.raw.expHM.DE10.group2_est.DE, 
      expHM.overall = mse.raw.expHM.DE10.overall_est.DE, 
      lnHM.1 =  mse.raw.lnHM.DE10.group1_est.DE, 
      lnHM.2 =  mse.raw.lnHM.DE10.group2_est.DE, 
      lnHM.overall = mse.raw.lnHM.DE10.overall_est.DE
    ), 
    all_genes = data.frame(
      edgeR.tag = mse.raw.edgeR.tag.DE10.all_genes, 
      edgeR.trend = mse.raw.edgeR.trend.DE10.all_genes, 
      DESeq2 = mse.raw.DESeq2.DE10.all_genes, 
      DSS.trend = mse.raw.DSS.trend.DE10.all_genes, 
      DSS.notrend = mse.raw.DSS.notrend.DE10.all_genes, 
      expHM.1 =  mse.raw.expHM.DE10.group1_est.all_genes, 
      expHM.2 =  mse.raw.expHM.DE10.group2_est.all_genes, 
      expHM.overall = mse.raw.expHM.DE10.overall_est.all_genes, 
      lnHM.1 =  mse.raw.lnHM.DE10.group1_est.all_genes, 
      lnHM.2 =  mse.raw.lnHM.DE10.group2_est.all_genes, 
      lnHM.overall = mse.raw.lnHM.DE10.overall_est.all_genes
    )
  ), 
  normalised.MSEs = list(
    noDE = data.frame(
      edgeR.tag = mse.norm.edgeR.tag.DE10.noDE, 
      edgeR.trend = mse.norm.edgeR.trend.DE10.noDE, 
      DESeq2 = mse.norm.DESeq2.DE10.noDE, 
      DSS.trend = mse.norm.DSS.trend.DE10.noDE, 
      DSS.notrend = mse.norm.DSS.notrend.DE10.noDE, 
      expHM.1 =  mse.norm.expHM.DE10.group1_est.noDE, 
      expHM.2 =  mse.norm.expHM.DE10.group2_est.noDE, 
      expHM.overall = mse.norm.expHM.DE10.overall_est.noDE, 
      lnHM.1 =  mse.norm.lnHM.DE10.group1_est.noDE, 
      lnHM.2 =  mse.norm.lnHM.DE10.group2_est.noDE, 
      lnHM.overall = mse.norm.lnHM.DE10.overall_est.noDE
    ), 
    DE = data.frame(
      edgeR.tag = mse.norm.edgeR.tag.DE10.DE, 
      edgeR.trend = mse.norm.edgeR.trend.DE10.DE, 
      DESeq2 = mse.norm.DESeq2.DE10.DE, 
      DSS.trend = mse.norm.DSS.trend.DE10.DE, 
      DSS.notrend = mse.norm.DSS.notrend.DE10.DE, 
      expHM.1 =  mse.norm.expHM.DE10.group1_est.DE, 
      expHM.2 =  mse.norm.expHM.DE10.group2_est.DE, 
      expHM.overall = mse.norm.expHM.DE10.overall_est.DE, 
      lnHM.1 =  mse.norm.lnHM.DE10.group1_est.DE, 
      lnHM.2 =  mse.norm.lnHM.DE10.group2_est.DE, 
      lnHM.overall = mse.norm.lnHM.DE10.overall_est.DE
    ), 
    all_genes = data.frame(
      edgeR.tag = mse.norm.edgeR.tag.DE10.all_genes, 
      edgeR.trend = mse.norm.edgeR.trend.DE10.all_genes, 
      DESeq2 = mse.norm.DESeq2.DE10.all_genes, 
      DSS.trend = mse.norm.DSS.trend.DE10.all_genes, 
      DSS.notrend = mse.norm.DSS.notrend.DE10.all_genes, 
      expHM.1 =  mse.norm.expHM.DE10.group1_est.all_genes, 
      expHM.2 =  mse.norm.expHM.DE10.group2_est.all_genes, 
      expHM.overall = mse.norm.expHM.DE10.overall_est.all_genes, 
      lnHM.1 =  mse.norm.lnHM.DE10.group1_est.all_genes, 
      lnHM.2 =  mse.norm.lnHM.DE10.group2_est.all_genes, 
      lnHM.overall = mse.norm.lnHM.DE10.overall_est.all_genes
    )
  )
)

DE20 <- list(
  raw.MSEs = list(
    noDE = data.frame(
      edgeR.tag = mse.raw.edgeR.tag.DE20.noDE, 
      edgeR.trend = mse.raw.edgeR.trend.DE20.noDE, 
      DESeq2 = mse.raw.DESeq2.DE20.noDE, 
      DSS.trend = mse.raw.DSS.trend.DE20.noDE, 
      DSS.notrend = mse.raw.DSS.notrend.DE20.noDE, 
      expHM.1 =  mse.raw.expHM.DE20.group1_est.noDE, 
      expHM.2 =  mse.raw.expHM.DE20.group2_est.noDE, 
      expHM.overall = mse.raw.expHM.DE20.overall_est.noDE, 
      lnHM.1 =  mse.raw.lnHM.DE20.group1_est.noDE, 
      lnHM.2 =  mse.raw.lnHM.DE20.group2_est.noDE, 
      lnHM.overall = mse.raw.lnHM.DE20.overall_est.noDE
    ), 
    DE = data.frame(
      edgeR.tag = mse.raw.edgeR.tag.DE20.DE, 
      edgeR.trend = mse.raw.edgeR.trend.DE20.DE, 
      DESeq2 = mse.raw.DESeq2.DE20.DE, 
      DSS.trend = mse.raw.DSS.trend.DE20.DE, 
      DSS.notrend = mse.raw.DSS.notrend.DE20.DE, 
      expHM.1 =  mse.raw.expHM.DE20.group1_est.DE, 
      expHM.2 =  mse.raw.expHM.DE20.group2_est.DE, 
      expHM.overall = mse.raw.expHM.DE20.overall_est.DE, 
      lnHM.1 =  mse.raw.lnHM.DE20.group1_est.DE, 
      lnHM.2 =  mse.raw.lnHM.DE20.group2_est.DE, 
      lnHM.overall = mse.raw.lnHM.DE20.overall_est.DE
    ), 
    all_genes = data.frame(
      edgeR.tag = mse.raw.edgeR.tag.DE20.all_genes, 
      edgeR.trend = mse.raw.edgeR.trend.DE20.all_genes, 
      DESeq2 = mse.raw.DESeq2.DE20.all_genes, 
      DSS.trend = mse.raw.DSS.trend.DE20.all_genes, 
      DSS.notrend = mse.raw.DSS.notrend.DE20.all_genes, 
      expHM.1 =  mse.raw.expHM.DE20.group1_est.all_genes, 
      expHM.2 =  mse.raw.expHM.DE20.group2_est.all_genes, 
      expHM.overall = mse.raw.expHM.DE20.overall_est.all_genes, 
      lnHM.1 =  mse.raw.lnHM.DE20.group1_est.all_genes, 
      lnHM.2 =  mse.raw.lnHM.DE20.group2_est.all_genes, 
      lnHM.overall = mse.raw.lnHM.DE20.overall_est.all_genes
    )
  ), 
  normalised.MSEs = list(
    noDE = data.frame(
      edgeR.tag = mse.norm.edgeR.tag.DE20.noDE, 
      edgeR.trend = mse.norm.edgeR.trend.DE20.noDE, 
      DESeq2 = mse.norm.DESeq2.DE20.noDE, 
      DSS.trend = mse.norm.DSS.trend.DE20.noDE, 
      DSS.notrend = mse.norm.DSS.notrend.DE20.noDE, 
      expHM.1 =  mse.norm.expHM.DE20.group1_est.noDE, 
      expHM.2 =  mse.norm.expHM.DE20.group2_est.noDE, 
      expHM.overall = mse.norm.expHM.DE20.overall_est.noDE, 
      lnHM.1 =  mse.norm.lnHM.DE20.group1_est.noDE, 
      lnHM.2 =  mse.norm.lnHM.DE20.group2_est.noDE, 
      lnHM.overall = mse.norm.lnHM.DE20.overall_est.noDE
    ), 
    DE = data.frame(
      edgeR.tag = mse.norm.edgeR.tag.DE20.DE, 
      edgeR.trend = mse.norm.edgeR.trend.DE20.DE, 
      DESeq2 = mse.norm.DESeq2.DE20.DE, 
      DSS.trend = mse.norm.DSS.trend.DE20.DE, 
      DSS.notrend = mse.norm.DSS.notrend.DE20.DE, 
      expHM.1 =  mse.norm.expHM.DE20.group1_est.DE, 
      expHM.2 =  mse.norm.expHM.DE20.group2_est.DE, 
      expHM.overall = mse.norm.expHM.DE20.overall_est.DE, 
      lnHM.1 =  mse.norm.lnHM.DE20.group1_est.DE, 
      lnHM.2 =  mse.norm.lnHM.DE20.group2_est.DE, 
      lnHM.overall = mse.norm.lnHM.DE20.overall_est.DE
    ), 
    all_genes = data.frame(
      edgeR.tag = mse.norm.edgeR.tag.DE20.all_genes, 
      edgeR.trend = mse.norm.edgeR.trend.DE20.all_genes, 
      DESeq2 = mse.norm.DESeq2.DE20.all_genes, 
      DSS.trend = mse.norm.DSS.trend.DE20.all_genes, 
      DSS.notrend = mse.norm.DSS.notrend.DE20.all_genes, 
      expHM.1 =  mse.norm.expHM.DE20.group1_est.all_genes, 
      expHM.2 =  mse.norm.expHM.DE20.group2_est.all_genes, 
      expHM.overall = mse.norm.expHM.DE20.overall_est.all_genes, 
      lnHM.1 =  mse.norm.lnHM.DE20.group1_est.all_genes, 
      lnHM.2 =  mse.norm.lnHM.DE20.group2_est.all_genes, 
      lnHM.overall = mse.norm.lnHM.DE20.overall_est.all_genes
    )
  )
)

saveRDS(DE2, file=here('Results/Dispersion estimation Dec 2020','mse.disp.DE2.rds'))
saveRDS(DE5, file=here('Results/Dispersion estimation Dec 2020','mse.disp.DE5.rds'))
saveRDS(DE10, file=here('Results/Dispersion estimation Dec 2020','mse.disp.DE10.rds'))
saveRDS(DE20, file=here('Results/Dispersion estimation Dec 2020','mse.disp.DE20.rds'))
# one huge dispersion outlier in DE2.1 makes estimates extremely high for that gene, so much so 
# that MSEs over all genes for DE2.1 are around 7 instead of 0.2-0.4 like the rest
# not sure what to do with this. clearly not representative but surely legitimately sampled by 
# compcodeR like all other data.
# save files for now anyway but will need to deal with this if I use these results (which I might 
# not - nodiff is the main result of interest, although really do want to also look at dispersion 
# estimation in presence of DE)
# (Above notes from original analysis Oct 2019, but still look relevant)


## Exploration of differences between methods and densities of predictions ####
## (From original analysis Oct 2019)
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
# lnHM better than edgeR on most ranges of dispersions

rbind(c(mean(edgeR.tag - mean(true.disps)), 
  mean(edgeR.trend - mean(true.disps)), 
  mean(DESeq2 - mean(true.disps)), 
  mean(DSS.notrend - mean(true.disps)), 
  mean(expHM.overall_est - mean(true.disps)), 
  mean(lnHM.overall_est) - mean(true.disps)), 
c(mean((edgeR.tag - true.disps)^2), 
  mean((edgeR.trend - true.disps)^2), 
  mean((DESeq2 - true.disps)^2), 
  mean((DSS.notrend - true.disps)^2), 
  mean((expHM.overall_est - true.disps)^2), 
  mean((lnHM.overall_est - true.disps)^2)))
mean(true.disps)
# approximate relationship between mean of estimates for each method and MSE

par(mfrow=c(3,2), mar=c(2,2,1,1), mgp=c(2,1,0))
plot(log(true.disps), 
     abs(lnHM.overall_est - true.disps)^2, pch=20, col='red', 
     xlim=c(-14,1.7), ylim=c(0,0.8))
# squared error v log disp for lnHM
plot(log(true.disps), 
     abs(edgeR.tag - true.disps)^2, pch=20, col='blue', 
     xlim=c(-14,1.7), ylim=c(0,0.8))
# squared error v disp for edgeR tag
plot(log(true.disps), 
     abs(DESeq2 - true.disps)^2, pch=20, col='yellow', 
     xlim=c(-14,1.7), ylim=c(0,0.8))
# squared error v disp for DESeq2
plot(log(true.disps), 
     abs(lnHM.overall_est - true.disps) - 
       abs(edgeR.tag - true.disps), pch=20)
# difference in absolute errors between lnHM and edgeR tag v disp
plot(log(true.disps), 
     abs(lnHM.overall_est - true.disps)^2 - 
       abs(edgeR.tag - true.disps)^2, pch=20)
# difference in squared errors between lnHM and edgeR tag v disp
mean(abs(lnHM.overall_est - true.disps) - 
       abs(edgeR.tag - true.disps))
mean(abs(lnHM.overall_est - true.disps)^2 - 
       abs(edgeR.tag - true.disps)^2)
mean(abs(lnHM.overall_est - true.disps) > 
       abs(edgeR.tag - true.disps))
mean(abs(lnHM.overall_est - true.disps) > 
       abs(DESeq2 - true.disps))
par(mfrow=c(2,2), mar=c(2,2,1,1), mgp=c(2,1,0))
plot(edgeR.tag, lnHM.overall_est, pch=20)
# lnHM estimates v edgeR tag estimates
mean(lnHM.overall_est > edgeR.tag)
# edgeR estimates bigger for most genes
plot(log(edgeR.tag), log(lnHM.overall_est), pch=20)
# lnHM estimates v edgeR tag estimates on log scale
# lnHM estimates cover wider range
plot(log(edgeR.tag[which(DE == 0)]), log(lnHM.overall_est[which(DE == 0)]), pch=20, xlim=c(-3,2), ylim=c(-5,2))
# lnHM estimates v edgeR tag estimates on log scale for non-DE genes
plot(log(edgeR.tag[which(DE == 1)]), log(lnHM.overall_est[which(DE == 1)]), pch=20, xlim=c(-3,2), ylim=c(-5,2))
# lnHM estimates v edgeR tag estimates on log scale for DE genes
# both have narrower range when there is DE, but could be because small number of DE genes (5%)
c(mean(abs(lnHM.overall_est - true.disps) > 
         abs(edgeR.tag - true.disps)), 
  mean(abs(lnHM.overall_est - true.disps)[which(DE==0)] > 
         abs(edgeR.tag - true.disps)[which(DE==0)]), 
  mean(abs(lnHM.overall_est - true.disps)[which(DE==1)] > 
         abs(edgeR.tag - true.disps)[which(DE==1)]))
# lnHM errors smaller than edgeR tag for ~75% of non-DE genes, bigger for 70% of DE genes
c(mean(abs(lnHM.overall_est - true.disps) > 
         abs(DESeq2 - true.disps)), 
  mean(abs(lnHM.overall_est - true.disps)[which(DE==0)] > 
         abs(DESeq2 - true.disps)[which(DE==0)]), 
  mean(abs(lnHM.overall_est - true.disps)[which(DE==1)] > 
         abs(DESeq2 - true.disps)[which(DE==1)]))
# lnHM errors bigger than DESeq2 errors for just over 50% of non-DE genes, 80% of DE genes

plot(true.disps, abs(lnHM.overall_est - true.disps) - abs(edgeR.tag - true.disps), 
     pch=20, xlim=c(0,8), ylim=c(-1,1))
# difference in errors between lnHM and edgeR tag v disp
plot(log(true.disps), abs(lnHM.overall_est - true.disps) - abs(edgeR.tag - true.disps), pch=20, 
     xlim=c(-16,3), ylim=c(-1,1))
# difference in errors between lnHM and edgeR tag v log disp
plot(true.disps[which(DE == 0)], xlim=c(0,8), ylim=c(-1,1), 
     abs(lnHM.overall_est - true.disps)[which(DE == 0)] - abs(edgeR.tag - true.disps)[which(DE == 0)], pch=20)
# difference in errors between lnHM and edgeR tag v disp for non-DE genes
plot(log(true.disps[which(DE == 0)]), xlim=c(-16,3), ylim=c(-1,1), 
     abs(lnHM.overall_est - true.disps)[which(DE == 0)] - abs(edgeR.tag - true.disps)[which(DE == 0)], pch=20)
# difference in errors between lnHM and edgeR tag v log disp for non-DE genes
plot(true.disps[which(DE == 1)], xlim=c(0,8), ylim=c(-1,1), 
     abs(lnHM.overall_est - true.disps)[which(DE == 1)] - abs(edgeR.tag - true.disps)[which(DE == 1)], pch=20)
# difference in errors between lnHM and edgeR tag v disp for DE genes
plot(log(true.disps[which(DE == 1)]), xlim=c(-16,3), ylim=c(-1,1), 
     abs(lnHM.overall_est - true.disps)[which(DE == 1)] - abs(edgeR.tag - true.disps)[which(DE == 1)], pch=20)
# difference in errors between lnHM and edgeR tag v log disp for DE genes

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
# Overall, lnHM almost always better than edgeR.tag for dispersion < 0.2, and both similar for > 0.2.
rbind(c('<0.1','0.1-0.2','0.2-0.3','0.3-0.4','0.4-0.5','0.5-1','1-2','>2'), 
      round(c(mean(abs(lnHM.overall_est - true.disps)[which(true.disps < 0.1 & DE == 1)] > 
          abs(edgeR.tag - true.disps)[which(true.disps < 0.1 & DE == 1)]),
        mean(abs(lnHM.overall_est - true.disps)[which(between(true.disps, 0.1, 0.2) & DE == 1)] > 
            abs(edgeR.tag - true.disps)[which(between(true.disps, 0.1, 0.2) & DE == 1)]),
        mean(abs(lnHM.overall_est - true.disps)[which(between(true.disps, 0.2, 0.3) & DE == 1)] > 
            abs(edgeR.tag - true.disps)[which(between(true.disps, 0.2, 0.3) & DE == 1)]),
        mean(abs(lnHM.overall_est - true.disps)[which(between(true.disps, 0.3, 0.4) & DE == 1)] > 
            abs(edgeR.tag - true.disps)[which(between(true.disps, 0.3, 0.4) & DE == 1)]),
        mean(abs(lnHM.overall_est - true.disps)[which(between(true.disps, 0.4, 0.5) & DE == 1)] > 
            abs(edgeR.tag - true.disps)[which(between(true.disps, 0.4, 0.5) & DE == 1)]),
        mean(abs(lnHM.overall_est - true.disps)[which(between(true.disps, 0.5, 1) & DE == 1)] > 
            abs(edgeR.tag - true.disps)[which(between(true.disps, 0.5, 1) & DE == 1)]),
        mean(abs(lnHM.overall_est - true.disps)[which(between(true.disps, 1, 2) & DE == 1)] > 
            abs(edgeR.tag - true.disps)[which(between(true.disps, 1, 2) & DE == 1)]),
        mean(abs(lnHM.overall_est - true.disps)[which(true.disps > 2 & DE == 1)] > 
            abs(edgeR.tag - true.disps)[which(true.disps > 2 & DE == 1)])), 2))
# For DE genes, lnHM almost always worse than edgeR.tag for dispersion < 0.5, both simliar for > 0.5
mean(log(edgeR.tag)) 
mean(log(lnHM.overall_est))
mean(log(true.disps))      
par(mfrow=c(2,2), mar=c(2,2,1,1), mgp=c(2,1,0))
hist(log(edgeR.tag))
hist(log(lnHM.overall_est))      
plot(density(log(edgeR.tag[which(DE == 0)])), col='blue', ylim=c(0,0.7), xlim=c(-5,2))
lines(density(log(lnHM.overall_est[which(DE == 0)])), col='red')
lines(density(log(true.disps[which(DE == 0)])), col='orange')
abline(v=mean(log(true.disps)))
median(edgeR.tag)
median(lnHM.overall_est)
median(true.disps)
plot(density(log(edgeR.tag[which(DE == 1)])), col='blue', ylim=c(0,0.8), xlim=c(-5,2))
lines(density(log(true.disps[which(DE == 1)])), col='orange')
lines(density(log(lnHM.overall_est[which(DE == 1)])), col='red')
abline(v=mean(log(true.disps[which(DE == 1)])))

# lnHM very close to true distribution for non-DE genes, edgeR.tag similar at high end but 
# overestimates at low end. For DE genes, edgeR.tag pattern is the same, but lnHM is sort of 
# like a really extreme version of the edgeR.tag pattern, or just generally massively 
# overestimating except at very high end.

rbind(colMeans(DE2$raw.MSEs$DE), colMeans(DE5$raw.MSEs$DE), 
      colMeans(DE10$raw.MSEs$DE), colMeans(DE20$raw.MSEs$DE))
names(res)
slotNames(res$data)
dim(res$data@count.matrix)
names(res$data@variable.annotations)
lfc <- res$data@variable.annotations$truelog2foldchanges
plot(lfc)
par(mfrow=c(2,2), mar=c(2,2,1,1), mgp=c(2,1,0))
plot(abs(lfc), abs(lnHM.overall_est - true.disps))
plot(abs(lfc), abs(edgeR.tag - true.disps))
plot(abs(lfc), log(lnHM.overall_est), ylim=c(-5,2))
plot(abs(lfc), log(edgeR.tag), ylim=c(-5,2))
mval <- res$data@variable.annotations$M.value
plot(abs(mval), abs(lnHM.overall_est - true.disps))
plot(abs(mval), abs(edgeR.tag - true.disps))
plot(abs(mval), log(lnHM.overall_est), ylim=c(-5,2))
plot(abs(mval), log(edgeR.tag), ylim=c(-5,2))


## Create and save density plots for one simulation per sample size ####
