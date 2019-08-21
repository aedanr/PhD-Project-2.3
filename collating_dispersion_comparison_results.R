library(here)
library(compcodeR)

for (i in c('.DEDD2','.DEDD5','.DEDD10')) {
  for (j in c('raw.','norm.')) {
    for (k in c('edgeR.trend', 'edgeR.tag', 'DESeq2', 'DSS.notrend', 'DSS.trend')) {
      assign(paste0('mse.', j, k, i, '.noDD.noDE'), numeric(50))
      assign(paste0('mse.', j, k, i, '.noDD.DE'), numeric(50))
      assign(paste0('mse.', j, k, i, '.noDD.all_genes'), numeric(50))
      assign(paste0('mse.', j, k, i, '.DD.group1.noDE'), numeric(50))
      assign(paste0('mse.', j, k, i, '.DD.group1.DE'), numeric(50))
      assign(paste0('mse.', j, k, i, '.DD.group1.all_genes'), numeric(50))
      assign(paste0('mse.', j, k, i, '.DD.group2.noDE'), numeric(50))
      assign(paste0('mse.', j, k, i, '.DD.group2.DE'), numeric(50))
      assign(paste0('mse.', j, k, i, '.DD.group2.all_genes'), numeric(50))
    }
    for (l in c('lnHM', 'expHM')) {
      assign(paste0('mse.', j, l, i, '.noDD.overall_est.noDE'), numeric(50))
      assign(paste0('mse.', j, l, i, '.noDD.group1_est.noDE'), numeric(50))
      assign(paste0('mse.', j, l, i, '.noDD.group2_est.noDE'), numeric(50))
      assign(paste0('mse.', j, l, i, '.noDD.overall_est.DE'), numeric(50))
      assign(paste0('mse.', j, l, i, '.noDD.group1_est.DE'), numeric(50))
      assign(paste0('mse.', j, l, i, '.noDD.group2_est.DE'), numeric(50))
      assign(paste0('mse.', j, l, i, '.noDD.overall_est.all_genes'), numeric(50))
      assign(paste0('mse.', j, l, i, '.noDD.group1_est.all_genes'), numeric(50))
      assign(paste0('mse.', j, l, i, '.noDD.group2_est.all_genes'), numeric(50))
      assign(paste0('mse.', j, l, i, '.DD.overall_est.group1.noDE'), numeric(50))
      assign(paste0('mse.', j, l, i, '.DD.overall_est.group1.DE'), numeric(50))
      assign(paste0('mse.', j, l, i, '.DD.overall_est.group1.all_genes'), numeric(50))
      assign(paste0('mse.', j, l, i, '.DD.overall_est.group2.noDE'), numeric(50))
      assign(paste0('mse.', j, l, i, '.DD.overall_est.group2.DE'), numeric(50))
      assign(paste0('mse.', j, l, i, '.DD.overall_est.group2.all_genes'), numeric(50))
      assign(paste0('mse.', j, l, i, '.DD.group1_est.group1.noDE'), numeric(50))
      assign(paste0('mse.', j, l, i, '.DD.group1_est.group1.DE'), numeric(50))
      assign(paste0('mse.', j, l, i, '.DD.group1_est.group1.all_genes'), numeric(50))
      assign(paste0('mse.', j, l, i, '.DD.group2_est.group2.noDE'), numeric(50))
      assign(paste0('mse.', j, l, i, '.DD.group2_est.group2.DE'), numeric(50))
      assign(paste0('mse.', j, l, i, '.DD.group2_est.group2.all_genes'), numeric(50))
    }
  }
}

# for (i in c('.DE2','.DE5','.DE10')) {
  for (j in c('raw.','norm.')) {
    for (k in c('edgeR.trend', 'edgeR.tag', 'DESeq2', 'DSS.notrend', 'DSS.trend')) {
      assign(paste0('mse.', j, k, i, '.noDE'), numeric(50))
      assign(paste0('mse.', j, k, i, '.DE'), numeric(50))
      assign(paste0('mse.', j, k, i, '.overall'), numeric(50))
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

# for (i in c('.DD2','.DD5','.DD10')) {
  for (j in c('raw.','norm.')) {
    for (k in c('edgeR.trend', 'edgeR.tag', 'DESeq2', 'DSS.notrend', 'DSS.trend')) {
      assign(paste0('mse.', j, k, i, '.noDD'), numeric(50))
      assign(paste0('mse.', j, k, i, '.DD.group1'), numeric(50))
      assign(paste0('mse.', j, k, i, '.DD.group2'), numeric(50))
    }
    for (l in c('lnHM', 'expHM')) {
      assign(paste0('mse.', j, l, i, '.noDD.overall_est'), numeric(50))
      assign(paste0('mse.', j, l, i, '.noDD.group1_est'), numeric(50))
      assign(paste0('mse.', j, l, i, '.noDD.group2_est'), numeric(50))
      assign(paste0('mse.', j, l, i, '.DD.overall_est.group1'), numeric(50))
      assign(paste0('mse.', j, l, i, '.DD.overall_est.group2'), numeric(50))
      assign(paste0('mse.', j, l, i, '.DD.group1_est.group1'), numeric(50))
      assign(paste0('mse.', j, l, i, '.DD.group2_est.group2'), numeric(50))
    }
  }
}

for (n in 1:50) {
  for (i in c('.DEDD2', '.DEDD5', '.DEDD10')) {
    res <- readRDS(here('Results/Dispersion estimation results Aug 2019', 
                        paste0('disp.results', i, '.', n, '.rds')))
    DD <- res$data@variable.annotations$differential.dispersion
    DE <- res$data@variable.annotations$differential.expression
    true.disps.group1 <- res$data@variable.annotations$truedispersions.S1
    true.disps.group2 <- res$data@variable.annotations$truedispersions.S2
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
        assign(paste0('mse.raw.', k, i, '.noDD.noDE'), 
               `[<-`(get(paste0('mse.raw.', k, i, '.noDD.noDE')), n, 
                     value=mean((get(k) - true.disps.group1)[which(DD == 0 & DE == 0)]^2)))
        assign(paste0('mse.raw.', k, i, '.noDD.DE'), 
               `[<-`(get(paste0('mse.raw.', k, i, '.noDD.DE')), n, 
                     value=mean((get(k) - true.disps.group1)[which(DD == 0 & DE == 1)]^2)))
        assign(paste0('mse.raw.', k, i, '.noDD.all_genes'), 
               `[<-`(get(paste0('mse.raw.', k, i, '.noDD.all_genes')), n, 
                     value=mean((get(k) - true.disps.group1)[which(DD == 0)]^2)))
        assign(paste0('mse.raw.', k, i, '.DD.group1.noDE'), 
               `[<-`(get(paste0('mse.raw.', k, i, '.DD.group1.noDE')), n, 
                     value=mean((get(k) - true.disps.group1)[which(DD == 1 & DE == 0)]^2)))
        assign(paste0('mse.raw.', k, i, '.DD.group1.DE'), 
               `[<-`(get(paste0('mse.raw.', k, i, '.DD.group1.DE')), n, 
                     value=mean((get(k) - true.disps.group1)[which(DD == 1 & DE == 1)]^2)))
        assign(paste0('mse.raw.', k, i, '.DD.group1.all_genes'), 
               `[<-`(get(paste0('mse.raw.', k, i, '.DD.group1.all_genes')), n, 
                     value=mean((get(k) - true.disps.group1)[which(DD == 1)]^2)))
        assign(paste0('mse.raw.', k, i, '.DD.group2.noDE'), 
               `[<-`(get(paste0('mse.raw.', k, i, '.DD.group2.noDE')), n, 
                     value=mean((get(k) - true.disps.group2)[which(DD == 1 & DE == 0)]^2)))
        assign(paste0('mse.raw.', k, i, '.DD.group2.DE'), 
               `[<-`(get(paste0('mse.raw.', k, i, '.DD.group2.DE')), n, 
                     value=mean((get(k) - true.disps.group2)[which(DD == 1 & DE == 1)]^2)))
        assign(paste0('mse.raw.', k, i, '.DD.group2.all_genes'), 
               `[<-`(get(paste0('mse.raw.', k, i, '.DD.group2.all_genes')), n, 
                     value=mean((get(k) - true.disps.group2)[which(DD == 1)]^2)))
        assign(paste0('mse.norm.', k, i, '.noDD.noDE'), 
               `[<-`(get(paste0('mse.norm.', k, i, '.noDD.noDE')), n, 
                     value=mean(((get(k) - true.disps.group1) / true.disps.group1)[which(DD == 0 & DE == 0)]^2)))
        assign(paste0('mse.norm.', k, i, '.noDD.DE'), 
               `[<-`(get(paste0('mse.norm.', k, i, '.noDD.DE')), n, 
                     value=mean(((get(k) - true.disps.group1) / true.disps.group1)[which(DD == 0 & DE == 1)]^2)))
        assign(paste0('mse.norm.', k, i, '.noDD.all_genes'), 
               `[<-`(get(paste0('mse.norm.', k, i, '.noDD.all_genes')), n, 
                     value=mean(((get(k) - true.disps.group1) / true.disps.group1)[which(DD == 0)]^2)))
        assign(paste0('mse.norm.', k, i, '.DD.group1.noDE'), 
               `[<-`(get(paste0('mse.norm.', k, i, '.DD.group1.noDE')), n, 
                     value=mean(((get(k) - true.disps.group1) / true.disps.group1)[which(DD == 1 & DE == 0)]^2)))
        assign(paste0('mse.norm.', k, i, '.DD.group1.DE'), 
               `[<-`(get(paste0('mse.norm.', k, i, '.DD.group1.DE')), n, 
                     value=mean(((get(k) - true.disps.group1) / true.disps.group1)[which(DD == 1 & DE == 1)]^2)))
        assign(paste0('mse.norm.', k, i, '.DD.group1.all_genes'), 
               `[<-`(get(paste0('mse.norm.', k, i, '.DD.group1.all_genes')), n, 
                     value=mean(((get(k) - true.disps.group1) / true.disps.group1)[which(DD == 1)]^2)))
        assign(paste0('mse.norm.', k, i, '.DD.group2.noDE'), 
               `[<-`(get(paste0('mse.norm.', k, i, '.DD.group2.noDE')), n, 
                     value=mean(((get(k) - true.disps.group2) / true.disps.group2)[which(DD == 1 & DE == 0)]^2)))
        assign(paste0('mse.norm.', k, i, '.DD.group2.DE'), 
               `[<-`(get(paste0('mse.norm.', k, i, '.DD.group2.DE')), n, 
                     value=mean(((get(k) - true.disps.group2) / true.disps.group2)[which(DD == 1 & DE == 1)]^2)))
        assign(paste0('mse.norm.', k, i, '.DD.group2.all_genes'), 
               `[<-`(get(paste0('mse.norm.', k, i, '.DD.group2.all_genes')), n, 
                     value=mean(((get(k) - true.disps.group2) / true.disps.group2)[which(DD == 1)]^2)))
      }
      for (l in c('lnHM', 'expHM')) {
        assign(paste0('mse.raw.', l, i, '.noDD.overall_est.noDE'), 
               `[<-`(get(paste0('mse.raw.', l, i, '.noDD.overall_est.noDE')), n, 
                     value=mean((get(paste0(l, '.overall_est')) - true.disps.group1)[which(DD == 0 & DE == 0)]^2)))
        assign(paste0('mse.raw.', l, i, '.noDD.group1_est.noDE'), 
               `[<-`(get(paste0('mse.raw.', l, i, '.noDD.group1_est.noDE')), n, 
                     value=mean((get(paste0(l, '.group1_est')) - true.disps.group1)[which(DD == 0 & DE == 0)]^2)))
        assign(paste0('mse.raw.', l, i, '.noDD.group2_est.noDE'), 
               `[<-`(get(paste0('mse.raw.', l, i, '.noDD.group2_est.noDE')), n, 
                     value=mean((get(paste0(l, '.group2_est')) - true.disps.group1)[which(DD == 0 & DE == 0)]^2)))
        assign(paste0('mse.raw.', l, i, '.noDD.overall_est.DE'), 
               `[<-`(get(paste0('mse.raw.', l, i, '.noDD.overall_est.DE')), n, 
                     value=mean((get(paste0(l, '.overall_est')) - true.disps.group1)[which(DD == 0 & DE == 1)]^2)))
        assign(paste0('mse.raw.', l, i, '.noDD.group1_est.DE'), 
               `[<-`(get(paste0('mse.raw.', l, i, '.noDD.group1_est.DE')), n, 
                     value=mean((get(paste0(l, '.group1_est')) - true.disps.group1)[which(DD == 0 & DE == 1)]^2)))
        assign(paste0('mse.raw.', l, i, '.noDD.group2_est.DE'), 
               `[<-`(get(paste0('mse.raw.', l, i, '.noDD.group2_est.DE')), n, 
                     value=mean((get(paste0(l, '.group2_est')) - true.disps.group1)[which(DD == 0 & DE == 1)]^2)))
        assign(paste0('mse.raw.', l, i, '.noDD.overall_est.all_genes'), 
               `[<-`(get(paste0('mse.raw.', l, i, '.noDD.overall_est.all_genes')), n, 
                     value=mean((get(paste0(l, '.overall_est')) - true.disps.group1)[which(DD == 0)]^2)))
        assign(paste0('mse.raw.', l, i, '.noDD.group1_est.all_genes'), 
               `[<-`(get(paste0('mse.raw.', l, i, '.noDD.group1_est.all_genes')), n, 
                     value=mean((get(paste0(l, '.group1_est')) - true.disps.group1)[which(DD == 0)]^2)))
        assign(paste0('mse.raw.', l, i, '.noDD.group2_est.all_genes'), 
               `[<-`(get(paste0('mse.raw.', l, i, '.noDD.group2_est.all_genes')), n, 
                     value=mean((get(paste0(l, '.group2_est')) - true.disps.group1)[which(DD == 0)]^2)))
        assign(paste0('mse.raw.', l, i, '.DD.overall_est.group1.noDE'), 
               `[<-`(get(paste0('mse.raw.', l, i, '.DD.overall_est.group1.noDE')), n, 
                     value=mean((get(paste0(l, '.overall_est')) - true.disps.group1)[which(DD == 1 & DE == 0)]^2)))
        assign(paste0('mse.raw.', l, i, '.DD.overall_est.group1.DE'), 
               `[<-`(get(paste0('mse.raw.', l, i, '.DD.overall_est.group1.DE')), n, 
                     value=mean((get(paste0(l, '.overall_est')) - true.disps.group1)[which(DD == 1 & DE == 1)]^2)))
        assign(paste0('mse.raw.', l, i, '.DD.overall_est.group1.all_genes'), 
               `[<-`(get(paste0('mse.raw.', l, i, '.DD.overall_est.group1.all_genes')), n, 
                     value=mean((get(paste0(l, '.overall_est')) - true.disps.group1)[which(DD == 1)]^2)))
        assign(paste0('mse.raw.', l, i, '.DD.overall_est.group2.noDE'), 
               `[<-`(get(paste0('mse.raw.', l, i, '.DD.overall_est.group2.noDE')), n, 
                     value=mean((get(paste0(l, '.overall_est')) - true.disps.group2)[which(DD == 1 & DE == 0)]^2)))
        assign(paste0('mse.raw.', l, i, '.DD.overall_est.group2.DE'), 
               `[<-`(get(paste0('mse.raw.', l, i, '.DD.overall_est.group2.DE')), n, 
                     value=mean((get(paste0(l, '.overall_est')) - true.disps.group2)[which(DD == 1 & DE == 1)]^2)))
        assign(paste0('mse.raw.', l, i, '.DD.overall_est.group2.all_genes'), 
               `[<-`(get(paste0('mse.raw.', l, i, '.DD.overall_est.group2.all_genes')), n, 
                     value=mean((get(paste0(l, '.overall_est')) - true.disps.group2)[which(DD == 1)]^2)))
        assign(paste0('mse.raw.', l, i, '.DD.group1_est.group1.noDE'), 
               `[<-`(get(paste0('mse.raw.', l, i, '.DD.group1_est.group1.noDE')), n, 
                     value=mean((get(paste0(l, '.group1_est')) - true.disps.group1)[which(DD == 1 & DE == 0)]^2)))
        assign(paste0('mse.raw.', l, i, '.DD.group1_est.group1.DE'), 
               `[<-`(get(paste0('mse.raw.', l, i, '.DD.group1_est.group1.DE')), n, 
                     value=mean((get(paste0(l, '.group1_est')) - true.disps.group1)[which(DD == 1 & DE == 1)]^2)))
        assign(paste0('mse.raw.', l, i, '.DD.group1_est.group1.all_genes'), 
               `[<-`(get(paste0('mse.raw.', l, i, '.DD.group1_est.group1.all_genes')), n, 
                     value=mean((get(paste0(l, '.group1_est')) - true.disps.group1)[which(DD == 1)]^2)))
        assign(paste0('mse.raw.', l, i, '.DD.group2_est.group2.noDE'), 
               `[<-`(get(paste0('mse.raw.', l, i, '.DD.group2_est.group2.noDE')), n, 
                     value=mean((get(paste0(l, '.group2_est')) - true.disps.group2)[which(DD == 1 & DE == 0)]^2)))
        assign(paste0('mse.raw.', l, i, '.DD.group2_est.group2.DE'), 
               `[<-`(get(paste0('mse.raw.', l, i, '.DD.group2_est.group2.DE')), n, 
                     value=mean((get(paste0(l, '.group2_est')) - true.disps.group2)[which(DD == 1 & DE == 1)]^2)))
        assign(paste0('mse.raw.', l, i, '.DD.group2_est.group2.all_genes'), 
               `[<-`(get(paste0('mse.raw.', l, i, '.DD.group2_est.group2.all_genes')), n, 
                     value=mean((get(paste0(l, '.group2_est')) - true.disps.group2)[which(DD == 1)]^2)))
        assign(paste0('mse.norm.', l, i, '.noDD.overall_est.noDE'), 
               `[<-`(get(paste0('mse.norm.', l, i, '.noDD.overall_est.noDE')), n, 
                     value=mean(((get(paste0(l, '.overall_est')) - true.disps.group1) / 
                                   true.disps.group1)[which(DD == 0 & DE == 0)]^2)))
        assign(paste0('mse.norm.', l, i, '.noDD.group1_est.noDE'), 
               `[<-`(get(paste0('mse.norm.', l, i, '.noDD.group1_est.noDE')), n, 
                     value=mean(((get(paste0(l, '.group1_est')) - true.disps.group1) / 
                                   true.disps.group1)[which(DD == 0 & DE == 0)]^2)))
        assign(paste0('mse.norm.', l, i, '.noDD.group2_est.noDE'), 
               `[<-`(get(paste0('mse.norm.', l, i, '.noDD.group2_est.noDE')), n, 
                     value=mean(((get(paste0(l, '.group2_est')) - true.disps.group1) / 
                                   true.disps.group1)[which(DD == 0 & DE == 0)]^2)))
        assign(paste0('mse.norm.', l, i, '.noDD.overall_est.DE'), 
               `[<-`(get(paste0('mse.norm.', l, i, '.noDD.overall_est.DE')), n, 
                     value=mean(((get(paste0(l, '.overall_est')) - true.disps.group1) / 
                                   true.disps.group1)[which(DD == 0 & DE == 1)]^2)))
        assign(paste0('mse.norm.', l, i, '.noDD.group1_est.DE'), 
               `[<-`(get(paste0('mse.norm.', l, i, '.noDD.group1_est.DE')), n, 
                     value=mean(((get(paste0(l, '.group1_est')) - true.disps.group1) / 
                                   true.disps.group1)[which(DD == 0 & DE == 1)]^2)))
        assign(paste0('mse.norm.', l, i, '.noDD.group2_est.DE'), 
               `[<-`(get(paste0('mse.norm.', l, i, '.noDD.group2_est.DE')), n, 
                     value=mean(((get(paste0(l, '.group2_est')) - true.disps.group1) / 
                                   true.disps.group1)[which(DD == 0 & DE == 1)]^2)))
        assign(paste0('mse.norm.', l, i, '.noDD.overall_est.all_genes'), 
               `[<-`(get(paste0('mse.norm.', l, i, '.noDD.overall_est.all_genes')), n, 
                     value=mean(((get(paste0(l, '.overall_est')) - true.disps.group1) / 
                                   true.disps.group1)[which(DD == 0)]^2)))
        assign(paste0('mse.norm.', l, i, '.noDD.group1_est.all_genes'), 
               `[<-`(get(paste0('mse.norm.', l, i, '.noDD.group1_est.all_genes')), n, 
                     value=mean(((get(paste0(l, '.group1_est')) - true.disps.group1) / 
                                   true.disps.group1)[which(DD == 0)]^2)))
        assign(paste0('mse.norm.', l, i, '.noDD.group2_est.all_genes'), 
               `[<-`(get(paste0('mse.norm.', l, i, '.noDD.group2_est.all_genes')), n, 
                     value=mean(((get(paste0(l, '.group2_est')) - true.disps.group1) / 
                                   true.disps.group1)[which(DD == 0)]^2)))
        assign(paste0('mse.norm.', l, i, '.DD.overall_est.group1.noDE'), 
               `[<-`(get(paste0('mse.norm.', l, i, '.DD.overall_est.group1.noDE')), n, 
                     value=mean(((get(paste0(l, '.overall_est')) - true.disps.group1) / 
                                   true.disps.group1)[which(DD == 1 & DE == 0)]^2)))
        assign(paste0('mse.norm.', l, i, '.DD.overall_est.group1.DE'), 
               `[<-`(get(paste0('mse.norm.', l, i, '.DD.overall_est.group1.DE')), n, 
                     value=mean(((get(paste0(l, '.overall_est')) - true.disps.group1) / 
                                   true.disps.group1)[which(DD == 1 & DE == 1)]^2)))
        assign(paste0('mse.norm.', l, i, '.DD.overall_est.group1.all_genes'), 
               `[<-`(get(paste0('mse.norm.', l, i, '.DD.overall_est.group1.all_genes')), n, 
                     value=mean(((get(paste0(l, '.overall_est')) - true.disps.group1) / 
                                   true.disps.group1)[which(DD == 1)]^2)))
        assign(paste0('mse.norm.', l, i, '.DD.overall_est.group2.noDE'), 
               `[<-`(get(paste0('mse.norm.', l, i, '.DD.overall_est.group2.noDE')), n, 
                     value=mean(((get(paste0(l, '.overall_est')) - true.disps.group2) / 
                                   true.disps.group2)[which(DD == 1 & DE == 0)]^2)))
        assign(paste0('mse.norm.', l, i, '.DD.overall_est.group2.DE'), 
               `[<-`(get(paste0('mse.norm.', l, i, '.DD.overall_est.group2.DE')), n, 
                     value=mean(((get(paste0(l, '.overall_est')) - true.disps.group2) / 
                                   true.disps.group2)[which(DD == 1 & DE == 1)]^2)))
        assign(paste0('mse.norm.', l, i, '.DD.overall_est.group2.all_genes'), 
               `[<-`(get(paste0('mse.norm.', l, i, '.DD.overall_est.group2.all_genes')), n, 
                     value=mean(((get(paste0(l, '.overall_est')) - true.disps.group2) / 
                                   true.disps.group2)[which(DD == 1)]^2)))
        assign(paste0('mse.norm.', l, i, '.DD.group1_est.group1.noDE'), 
               `[<-`(get(paste0('mse.norm.', l, i, '.DD.group1_est.group1.noDE')), n, 
                     value=mean(((get(paste0(l, '.group1_est')) - true.disps.group1) / 
                                   true.disps.group1)[which(DD == 1 & DE == 0)]^2)))
        assign(paste0('mse.norm.', l, i, '.DD.group1_est.group1.DE'), 
               `[<-`(get(paste0('mse.norm.', l, i, '.DD.group1_est.group1.DE')), n, 
                     value=mean(((get(paste0(l, '.group1_est')) - true.disps.group1) / 
                                   true.disps.group1)[which(DD == 1 & DE == 1)]^2)))
        assign(paste0('mse.norm.', l, i, '.DD.group1_est.group1.all_genes'), 
               `[<-`(get(paste0('mse.norm.', l, i, '.DD.group1_est.group1.all_genes')), n, 
                     value=mean(((get(paste0(l, '.group1_est')) - true.disps.group1) / 
                                   true.disps.group1)[which(DD == 1)]^2)))
        assign(paste0('mse.norm.', l, i, '.DD.group2_est.group2.noDE'), 
               `[<-`(get(paste0('mse.norm.', l, i, '.DD.group2_est.group2.noDE')), n, 
                     value=mean(((get(paste0(l, '.group2_est')) - true.disps.group2) / 
                                   true.disps.group2)[which(DD == 1 & DE == 0)]^2)))
        assign(paste0('mse.norm.', l, i, '.DD.group2_est.group2.DE'), 
               `[<-`(get(paste0('mse.norm.', l, i, '.DD.group2_est.group2.DE')), n, 
                     value=mean(((get(paste0(l, '.group2_est')) - true.disps.group2) / 
                                   true.disps.group2)[which(DD == 1 & DE == 1)]^2)))
        assign(paste0('mse.norm.', l, i, '.DD.group2_est.group2.all_genes'), 
               `[<-`(get(paste0('mse.norm.', l, i, '.DD.group2_est.group2.all_genes')), n, 
                     value=mean(((get(paste0(l, '.group2_est')) - true.disps.group2) / 
                                   true.disps.group2)[which(DD == 1)]^2)))
      }
  }
}


DEDD2 <- list(
  raw.MSEs = list(
    noDD.noDE = data.frame(
      edgeR.tag = mse.raw.edgeR.tag.DEDD2.noDD.noDE, 
      edgeR.trend = mse.raw.edgeR.trend.DEDD2.noDD.noDE, 
      DESeq2 = mse.raw.DESeq2.DEDD2.noDD.noDE, 
      DSS.trend = mse.raw.DSS.trend.DEDD2.noDD.noDE, 
      DSS.notrend = mse.raw.DSS.notrend.DEDD2.noDD.noDE, 
      expHM.group1 = mse.raw.expHM.DEDD2.noDD.group1_est.noDE, 
      expHM.group2 = mse.raw.expHM.DEDD2.noDD.group2_est.noDE, 
      expHM.overall = mse.raw.expHM.DEDD2.noDD.overall_est.noDE, 
      lnHM.group1 = mse.raw.lnHM.DEDD2.noDD.group1_est.noDE, 
      lnHM.group2 = mse.raw.lnHM.DEDD2.noDD.group2_est.noDE, 
      lnHM.overall = mse.raw.lnHM.DEDD2.noDD.overall_est.noDE
    ), 
    noDD.DE = data.frame(
      edgeR.tag = mse.raw.edgeR.tag.DEDD2.noDD.DE, 
      edgeR.trend = mse.raw.edgeR.trend.DEDD2.noDD.DE, 
      DESeq2 = mse.raw.DESeq2.DEDD2.noDD.DE, 
      DSS.trend = mse.raw.DSS.trend.DEDD2.noDD.DE, 
      DSS.notrend = mse.raw.DSS.notrend.DEDD2.noDD.DE, 
      expHM.group1 = mse.raw.expHM.DEDD2.noDD.group1_est.DE, 
      expHM.group2 = mse.raw.expHM.DEDD2.noDD.group2_est.DE, 
      expHM.overall = mse.raw.expHM.DEDD2.noDD.overall_est.DE, 
      lnHM.group1 = mse.raw.lnHM.DEDD2.noDD.group1_est.DE, 
      lnHM.group2 = mse.raw.lnHM.DEDD2.noDD.group2_est.DE, 
      lnHM.overall = mse.raw.lnHM.DEDD2.noDD.overall_est.DE
    ), 
    noDD.all_genes = data.frame(
      edgeR.tag = mse.raw.edgeR.tag.DEDD2.noDD.all_genes, 
      edgeR.trend = mse.raw.edgeR.trend.DEDD2.noDD.all_genes, 
      DESeq2 = mse.raw.DESeq2.DEDD2.noDD.all_genes, 
      DSS.trend = mse.raw.DSS.trend.DEDD2.noDD.all_genes, 
      DSS.notrend = mse.raw.DSS.notrend.DEDD2.noDD.all_genes, 
      expHM.group1 = mse.raw.expHM.DEDD2.noDD.group1_est.all_genes, 
      expHM.group2 = mse.raw.expHM.DEDD2.noDD.group2_est.all_genes, 
      expHM.overall = mse.raw.expHM.DEDD2.noDD.overall_est.all_genes, 
      lnHM.group1 = mse.raw.lnHM.DEDD2.noDD.group1_est.all_genes, 
      lnHM.group2 = mse.raw.lnHM.DEDD2.noDD.group2_est.all_genes, 
      lnHM.overall = mse.raw.lnHM.DEDD2.noDD.overall_est.all_genes
    ),
    DD.noDE = data.frame(
      edgeR.tag = mse.raw.edgeR.tag.DEDD2.DD.group1.noDE, 
      edgeR.tag = mse.raw.edgeR.tag.DEDD2.DD.group2.noDE, 
      edgeR.trend = mse.raw.edgeR.trend.DEDD2.DD.group1.noDE, 
      edgeR.trend = mse.raw.edgeR.trend.DEDD2.DD.group2.noDE, 
      DESeq2 = mse.raw.DESeq2.DEDD2.DD.group1.noDE, 
      DESeq2 = mse.raw.DESeq2.DEDD2.DD.group2.noDE, 
      DSS.trend = mse.raw.DSS.trend.DEDD2.DD.group1.noDE, 
      DSS.trend = mse.raw.DSS.trend.DEDD2.DD.group2.noDE, 
      DSS.notrend = mse.raw.DSS.notrend.DEDD2.DD.group1.noDE, 
      DSS.notrend = mse.raw.DSS.notrend.DEDD2.DD.group2.noDE, 
      expHM.group1 = mse.raw.expHM.DEDD2.DD.group1_est.group1.noDE, 
      expHM.group2 = mse.raw.expHM.DEDD2.DD.group2_est.group2.noDE, 
      expHM.overall = mse.raw.expHM.DEDD2.DD.overall_est.group1.noDE, 
      expHM.overall = mse.raw.expHM.DEDD2.DD.overall_est.group2.noDE, 
      lnHM.group1 = mse.raw.lnHM.DEDD2.DD.group1_est.group1.noDE, 
      lnHM.group2 = mse.raw.lnHM.DEDD2.DD.group2_est.group2.noDE, 
      lnHM.overall = mse.raw.lnHM.DEDD2.DD.overall_est.group1.noDE, 
      lnHM.overall = mse.raw.lnHM.DEDD2.DD.overall_est.group2.noDE
    ), 
    DD.DE = data.frame(
      edgeR.tag = mse.raw.edgeR.tag.DEDD2.DD.group1.DE, 
      edgeR.tag = mse.raw.edgeR.tag.DEDD2.DD.group2.DE, 
      edgeR.trend = mse.raw.edgeR.trend.DEDD2.DD.group1.DE, 
      edgeR.trend = mse.raw.edgeR.trend.DEDD2.DD.group2.DE, 
      DESeq2 = mse.raw.DESeq2.DEDD2.DD.group1.DE, 
      DESeq2 = mse.raw.DESeq2.DEDD2.DD.group2.DE, 
      DSS.trend = mse.raw.DSS.trend.DEDD2.DD.group1.DE, 
      DSS.trend = mse.raw.DSS.trend.DEDD2.DD.group2.DE, 
      DSS.notrend = mse.raw.DSS.notrend.DEDD2.DD.group1.DE, 
      DSS.notrend = mse.raw.DSS.notrend.DEDD2.DD.group2.DE, 
      expHM.group1 = mse.raw.expHM.DEDD2.DD.group1_est.group1.DE, 
      expHM.group2 = mse.raw.expHM.DEDD2.DD.group2_est.group2.DE, 
      expHM.overall = mse.raw.expHM.DEDD2.DD.overall_est.group1.DE, 
      expHM.overall = mse.raw.expHM.DEDD2.DD.overall_est.group2.DE, 
      lnHM.group1 = mse.raw.lnHM.DEDD2.DD.group1_est.group1.DE, 
      lnHM.group2 = mse.raw.lnHM.DEDD2.DD.group2_est.group2.DE, 
      lnHM.overall = mse.raw.lnHM.DEDD2.DD.overall_est.group1.DE, 
      lnHM.overall = mse.raw.lnHM.DEDD2.DD.overall_est.group2.DE
    ), 
    DD.all_genes = data.frame(
      edgeR.tag = mse.raw.edgeR.tag.DEDD2.DD.group1.all_genes, 
      edgeR.tag = mse.raw.edgeR.tag.DEDD2.DD.group2.all_genes, 
      edgeR.trend = mse.raw.edgeR.trend.DEDD2.DD.group1.all_genes, 
      edgeR.trend = mse.raw.edgeR.trend.DEDD2.DD.group2.all_genes, 
      DESeq2 = mse.raw.DESeq2.DEDD2.DD.group1.all_genes, 
      DESeq2 = mse.raw.DESeq2.DEDD2.DD.group2.all_genes, 
      DSS.trend = mse.raw.DSS.trend.DEDD2.DD.group1.all_genes, 
      DSS.trend = mse.raw.DSS.trend.DEDD2.DD.group2.all_genes, 
      DSS.notrend = mse.raw.DSS.notrend.DEDD2.DD.group1.all_genes, 
      DSS.notrend = mse.raw.DSS.notrend.DEDD2.DD.group2.all_genes, 
      expHM.group1 = mse.raw.expHM.DEDD2.DD.group1_est.group1.all_genes, 
      expHM.group2 = mse.raw.expHM.DEDD2.DD.group2_est.group2.all_genes, 
      expHM.overall = mse.raw.expHM.DEDD2.DD.overall_est.group1.all_genes, 
      expHM.overall = mse.raw.expHM.DEDD2.DD.overall_est.group2.all_genes, 
      lnHM.group1 = mse.raw.lnHM.DEDD2.DD.group1_est.group1.all_genes, 
      lnHM.group2 = mse.raw.lnHM.DEDD2.DD.group2_est.group2.all_genes, 
      lnHM.overall = mse.raw.lnHM.DEDD2.DD.overall_est.group1.all_genes, 
      lnHM.overall = mse.raw.lnHM.DEDD2.DD.overall_est.group2.all_genes
    )
  ), 
  normalised.MSEs = list(
    noDD.noDE = data.frame(
      edgeR.tag = mse.norm.edgeR.tag.DEDD2.noDD.noDE, 
      edgeR.trend = mse.norm.edgeR.trend.DEDD2.noDD.noDE, 
      DESeq2 = mse.norm.DESeq2.DEDD2.noDD.noDE, 
      DSS.trend = mse.norm.DSS.trend.DEDD2.noDD.noDE, 
      DSS.notrend = mse.norm.DSS.notrend.DEDD2.noDD.noDE, 
      expHM.group1 = mse.norm.expHM.DEDD2.noDD.group1_est.noDE, 
      expHM.group2 = mse.norm.expHM.DEDD2.noDD.group2_est.noDE, 
      expHM.overall = mse.norm.expHM.DEDD2.noDD.overall_est.noDE, 
      lnHM.group1 = mse.norm.lnHM.DEDD2.noDD.group1_est.noDE, 
      lnHM.group2 = mse.norm.lnHM.DEDD2.noDD.group2_est.noDE, 
      lnHM.overall = mse.norm.lnHM.DEDD2.noDD.overall_est.noDE
    ), 
    noDD.DE = data.frame(
      edgeR.tag = mse.norm.edgeR.tag.DEDD2.noDD.DE, 
      edgeR.trend = mse.norm.edgeR.trend.DEDD2.noDD.DE, 
      DESeq2 = mse.norm.DESeq2.DEDD2.noDD.DE, 
      DSS.trend = mse.norm.DSS.trend.DEDD2.noDD.DE, 
      DSS.notrend = mse.norm.DSS.notrend.DEDD2.noDD.DE, 
      expHM.group1 = mse.norm.expHM.DEDD2.noDD.group1_est.DE, 
      expHM.group2 = mse.norm.expHM.DEDD2.noDD.group2_est.DE, 
      expHM.overall = mse.norm.expHM.DEDD2.noDD.overall_est.DE, 
      lnHM.group1 = mse.norm.lnHM.DEDD2.noDD.group1_est.DE, 
      lnHM.group2 = mse.norm.lnHM.DEDD2.noDD.group2_est.DE, 
      lnHM.overall = mse.norm.lnHM.DEDD2.noDD.overall_est.DE
    ), 
    noDD.all_genes = data.frame(
      edgeR.tag = mse.norm.edgeR.tag.DEDD2.noDD.all_genes, 
      edgeR.trend = mse.norm.edgeR.trend.DEDD2.noDD.all_genes, 
      DESeq2 = mse.norm.DESeq2.DEDD2.noDD.all_genes, 
      DSS.trend = mse.norm.DSS.trend.DEDD2.noDD.all_genes, 
      DSS.notrend = mse.norm.DSS.notrend.DEDD2.noDD.all_genes, 
      expHM.group1 = mse.norm.expHM.DEDD2.noDD.group1_est.all_genes, 
      expHM.group2 = mse.norm.expHM.DEDD2.noDD.group2_est.all_genes, 
      expHM.overall = mse.norm.expHM.DEDD2.noDD.overall_est.all_genes, 
      lnHM.group1 = mse.norm.lnHM.DEDD2.noDD.group1_est.all_genes, 
      lnHM.group2 = mse.norm.lnHM.DEDD2.noDD.group2_est.all_genes, 
      lnHM.overall = mse.norm.lnHM.DEDD2.noDD.overall_est.all_genes
    ),
    DD.noDE = data.frame(
      edgeR.tag = mse.norm.edgeR.tag.DEDD2.DD.group1.noDE, 
      edgeR.tag = mse.norm.edgeR.tag.DEDD2.DD.group2.noDE, 
      edgeR.trend = mse.norm.edgeR.trend.DEDD2.DD.group1.noDE, 
      edgeR.trend = mse.norm.edgeR.trend.DEDD2.DD.group2.noDE, 
      DESeq2 = mse.norm.DESeq2.DEDD2.DD.group1.noDE, 
      DESeq2 = mse.norm.DESeq2.DEDD2.DD.group2.noDE, 
      DSS.trend = mse.norm.DSS.trend.DEDD2.DD.group1.noDE, 
      DSS.trend = mse.norm.DSS.trend.DEDD2.DD.group2.noDE, 
      DSS.notrend = mse.norm.DSS.notrend.DEDD2.DD.group1.noDE, 
      DSS.notrend = mse.norm.DSS.notrend.DEDD2.DD.group2.noDE, 
      expHM.group1 = mse.norm.expHM.DEDD2.DD.group1_est.group1.noDE, 
      expHM.group2 = mse.norm.expHM.DEDD2.DD.group2_est.group2.noDE, 
      expHM.overall = mse.norm.expHM.DEDD2.DD.overall_est.group1.noDE, 
      expHM.overall = mse.norm.expHM.DEDD2.DD.overall_est.group2.noDE, 
      lnHM.group1 = mse.norm.lnHM.DEDD2.DD.group1_est.group1.noDE, 
      lnHM.group2 = mse.norm.lnHM.DEDD2.DD.group2_est.group2.noDE, 
      lnHM.overall = mse.norm.lnHM.DEDD2.DD.overall_est.group1.noDE, 
      lnHM.overall = mse.norm.lnHM.DEDD2.DD.overall_est.group2.noDE
    ), 
    DD.DE = data.frame(
      edgeR.tag = mse.norm.edgeR.tag.DEDD2.DD.group1.DE, 
      edgeR.tag = mse.norm.edgeR.tag.DEDD2.DD.group2.DE, 
      edgeR.trend = mse.norm.edgeR.trend.DEDD2.DD.group1.DE, 
      edgeR.trend = mse.norm.edgeR.trend.DEDD2.DD.group2.DE, 
      DESeq2 = mse.norm.DESeq2.DEDD2.DD.group1.DE, 
      DESeq2 = mse.norm.DESeq2.DEDD2.DD.group2.DE, 
      DSS.trend = mse.norm.DSS.trend.DEDD2.DD.group1.DE, 
      DSS.trend = mse.norm.DSS.trend.DEDD2.DD.group2.DE, 
      DSS.notrend = mse.norm.DSS.notrend.DEDD2.DD.group1.DE, 
      DSS.notrend = mse.norm.DSS.notrend.DEDD2.DD.group2.DE, 
      expHM.group1 = mse.norm.expHM.DEDD2.DD.group1_est.group1.DE, 
      expHM.group2 = mse.norm.expHM.DEDD2.DD.group2_est.group2.DE, 
      expHM.overall = mse.norm.expHM.DEDD2.DD.overall_est.group1.DE, 
      expHM.overall = mse.norm.expHM.DEDD2.DD.overall_est.group2.DE, 
      lnHM.group1 = mse.norm.lnHM.DEDD2.DD.group1_est.group1.DE, 
      lnHM.group2 = mse.norm.lnHM.DEDD2.DD.group2_est.group2.DE, 
      lnHM.overall = mse.norm.lnHM.DEDD2.DD.overall_est.group1.DE, 
      lnHM.overall = mse.norm.lnHM.DEDD2.DD.overall_est.group2.DE
    ), 
    DD.all_genes = data.frame(
      edgeR.tag = mse.norm.edgeR.tag.DEDD2.DD.group1.all_genes, 
      edgeR.tag = mse.norm.edgeR.tag.DEDD2.DD.group2.all_genes, 
      edgeR.trend = mse.norm.edgeR.trend.DEDD2.DD.group1.all_genes, 
      edgeR.trend = mse.norm.edgeR.trend.DEDD2.DD.group2.all_genes, 
      DESeq2 = mse.norm.DESeq2.DEDD2.DD.group1.all_genes, 
      DESeq2 = mse.norm.DESeq2.DEDD2.DD.group2.all_genes, 
      DSS.trend = mse.norm.DSS.trend.DEDD2.DD.group1.all_genes, 
      DSS.trend = mse.norm.DSS.trend.DEDD2.DD.group2.all_genes, 
      DSS.notrend = mse.norm.DSS.notrend.DEDD2.DD.group1.all_genes, 
      DSS.notrend = mse.norm.DSS.notrend.DEDD2.DD.group2.all_genes, 
      expHM.group1 = mse.norm.expHM.DEDD2.DD.group1_est.group1.all_genes, 
      expHM.group2 = mse.norm.expHM.DEDD2.DD.group2_est.group2.all_genes, 
      expHM.overall = mse.norm.expHM.DEDD2.DD.overall_est.group1.all_genes, 
      expHM.overall = mse.norm.expHM.DEDD2.DD.overall_est.group2.all_genes, 
      lnHM.group1 = mse.norm.lnHM.DEDD2.DD.group1_est.group1.all_genes, 
      lnHM.group2 = mse.norm.lnHM.DEDD2.DD.group2_est.group2.all_genes, 
      lnHM.overall = mse.norm.lnHM.DEDD2.DD.overall_est.group1.all_genes, 
      lnHM.overall = mse.norm.lnHM.DEDD2.DD.overall_est.group2.all_genes
    )
  )
)





