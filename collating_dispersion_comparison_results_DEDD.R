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
    lnHM.overall_est <- res$disps.lnHM.1
    lnHM.overall_est <- res$disps.lnHM.2
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
      expHM.1 =  mse.raw.expHM.DEDD2.noDD.group1_est.noDE, 
      expHM.2 =  mse.raw.expHM.DEDD2.noDD.group2_est.noDE, 
      expHM.overall = mse.raw.expHM.DEDD2.noDD.overall_est.noDE, 
      lnHM.1 =  mse.raw.lnHM.DEDD2.noDD.group1_est.noDE, 
      lnHM.2 =  mse.raw.lnHM.DEDD2.noDD.group2_est.noDE, 
      lnHM.overall = mse.raw.lnHM.DEDD2.noDD.overall_est.noDE
    ), 
    noDD.DE = data.frame(
      edgeR.tag = mse.raw.edgeR.tag.DEDD2.noDD.DE, 
      edgeR.trend = mse.raw.edgeR.trend.DEDD2.noDD.DE, 
      DESeq2 = mse.raw.DESeq2.DEDD2.noDD.DE, 
      DSS.trend = mse.raw.DSS.trend.DEDD2.noDD.DE, 
      DSS.notrend = mse.raw.DSS.notrend.DEDD2.noDD.DE, 
      expHM.1 =  mse.raw.expHM.DEDD2.noDD.group1_est.DE, 
      expHM.2 =  mse.raw.expHM.DEDD2.noDD.group2_est.DE, 
      expHM.overall = mse.raw.expHM.DEDD2.noDD.overall_est.DE, 
      lnHM.1 =  mse.raw.lnHM.DEDD2.noDD.group1_est.DE, 
      lnHM.2 =  mse.raw.lnHM.DEDD2.noDD.group2_est.DE, 
      lnHM.overall = mse.raw.lnHM.DEDD2.noDD.overall_est.DE
    ), 
    noDD.all_genes = data.frame(
      edgeR.tag = mse.raw.edgeR.tag.DEDD2.noDD.all_genes, 
      edgeR.trend = mse.raw.edgeR.trend.DEDD2.noDD.all_genes, 
      DESeq2 = mse.raw.DESeq2.DEDD2.noDD.all_genes, 
      DSS.trend = mse.raw.DSS.trend.DEDD2.noDD.all_genes, 
      DSS.notrend = mse.raw.DSS.notrend.DEDD2.noDD.all_genes, 
      expHM.1 =  mse.raw.expHM.DEDD2.noDD.group1_est.all_genes, 
      expHM.2 =  mse.raw.expHM.DEDD2.noDD.group2_est.all_genes, 
      expHM.overall = mse.raw.expHM.DEDD2.noDD.overall_est.all_genes, 
      lnHM.1 =  mse.raw.lnHM.DEDD2.noDD.group1_est.all_genes, 
      lnHM.2 =  mse.raw.lnHM.DEDD2.noDD.group2_est.all_genes, 
      lnHM.overall = mse.raw.lnHM.DEDD2.noDD.overall_est.all_genes
    ),
    DD.noDE = data.frame(
      edgeR.tag.1 =  mse.raw.edgeR.tag.DEDD2.DD.group1.noDE, 
      edgeR.tag.2 =  mse.raw.edgeR.tag.DEDD2.DD.group2.noDE, 
      edgeR.trend.1 =  mse.raw.edgeR.trend.DEDD2.DD.group1.noDE, 
      edgeR.trend.2 =  mse.raw.edgeR.trend.DEDD2.DD.group2.noDE, 
      DESeq2.1 =  mse.raw.DESeq2.DEDD2.DD.group1.noDE, 
      DESeq2.2 =  mse.raw.DESeq2.DEDD2.DD.group2.noDE, 
      DSS.trend.1 =  mse.raw.DSS.trend.DEDD2.DD.group1.noDE, 
      DSS.trend.2 =  mse.raw.DSS.trend.DEDD2.DD.group2.noDE, 
      DSS.notrend.1 =  mse.raw.DSS.notrend.DEDD2.DD.group1.noDE, 
      DSS.notrend.2 =  mse.raw.DSS.notrend.DEDD2.DD.group2.noDE, 
      expHM.1 =  mse.raw.expHM.DEDD2.DD.group1_est.group1.noDE, 
      expHM.2 =  mse.raw.expHM.DEDD2.DD.group2_est.group2.noDE, 
      expHM.overall.1 =  mse.raw.expHM.DEDD2.DD.overall_est.group1.noDE, 
      expHM.overall.2 =  mse.raw.expHM.DEDD2.DD.overall_est.group2.noDE, 
      lnHM.1 =  mse.raw.lnHM.DEDD2.DD.group1_est.group1.noDE, 
      lnHM.2 =  mse.raw.lnHM.DEDD2.DD.group2_est.group2.noDE, 
      lnHM.overall.1 =  mse.raw.lnHM.DEDD2.DD.overall_est.group1.noDE, 
      lnHM.overall.2 =  mse.raw.lnHM.DEDD2.DD.overall_est.group2.noDE
    ), 
    DD.DE = data.frame(
      edgeR.tag.1 =  mse.raw.edgeR.tag.DEDD2.DD.group1.DE, 
      edgeR.tag.2 =  mse.raw.edgeR.tag.DEDD2.DD.group2.DE, 
      edgeR.trend.1 =  mse.raw.edgeR.trend.DEDD2.DD.group1.DE, 
      edgeR.trend.2 =  mse.raw.edgeR.trend.DEDD2.DD.group2.DE, 
      DESeq2.1 =  mse.raw.DESeq2.DEDD2.DD.group1.DE, 
      DESeq2.2 =  mse.raw.DESeq2.DEDD2.DD.group2.DE, 
      DSS.trend.1 =  mse.raw.DSS.trend.DEDD2.DD.group1.DE, 
      DSS.trend.2 =  mse.raw.DSS.trend.DEDD2.DD.group2.DE, 
      DSS.notrend.1 =  mse.raw.DSS.notrend.DEDD2.DD.group1.DE, 
      DSS.notrend.2 =  mse.raw.DSS.notrend.DEDD2.DD.group2.DE, 
      expHM.1 =  mse.raw.expHM.DEDD2.DD.group1_est.group1.DE, 
      expHM.2 =  mse.raw.expHM.DEDD2.DD.group2_est.group2.DE, 
      expHM.overall.1 =  mse.raw.expHM.DEDD2.DD.overall_est.group1.DE, 
      expHM.overall.2 =  mse.raw.expHM.DEDD2.DD.overall_est.group2.DE, 
      lnHM.1 =  mse.raw.lnHM.DEDD2.DD.group1_est.group1.DE, 
      lnHM.2 =  mse.raw.lnHM.DEDD2.DD.group2_est.group2.DE, 
      lnHM.overall.1 =  mse.raw.lnHM.DEDD2.DD.overall_est.group1.DE, 
      lnHM.overall.2 =  mse.raw.lnHM.DEDD2.DD.overall_est.group2.DE
    ), 
    DD.all_genes = data.frame(
      edgeR.tag.1 =  mse.raw.edgeR.tag.DEDD2.DD.group1.all_genes, 
      edgeR.tag.2 =  mse.raw.edgeR.tag.DEDD2.DD.group2.all_genes, 
      edgeR.trend.1 =  mse.raw.edgeR.trend.DEDD2.DD.group1.all_genes, 
      edgeR.trend.2 =  mse.raw.edgeR.trend.DEDD2.DD.group2.all_genes, 
      DESeq2.1 =  mse.raw.DESeq2.DEDD2.DD.group1.all_genes, 
      DESeq2.2 =  mse.raw.DESeq2.DEDD2.DD.group2.all_genes, 
      DSS.trend.1 =  mse.raw.DSS.trend.DEDD2.DD.group1.all_genes, 
      DSS.trend.2 =  mse.raw.DSS.trend.DEDD2.DD.group2.all_genes, 
      DSS.notrend.1 =  mse.raw.DSS.notrend.DEDD2.DD.group1.all_genes, 
      DSS.notrend.2 =  mse.raw.DSS.notrend.DEDD2.DD.group2.all_genes, 
      expHM.1 =  mse.raw.expHM.DEDD2.DD.group1_est.group1.all_genes, 
      expHM.2 =  mse.raw.expHM.DEDD2.DD.group2_est.group2.all_genes, 
      expHM.overall.1 =  mse.raw.expHM.DEDD2.DD.overall_est.group1.all_genes, 
      expHM.overall.2 =  mse.raw.expHM.DEDD2.DD.overall_est.group2.all_genes, 
      lnHM.1 =  mse.raw.lnHM.DEDD2.DD.group1_est.group1.all_genes, 
      lnHM.2 =  mse.raw.lnHM.DEDD2.DD.group2_est.group2.all_genes, 
      lnHM.overall.1 =  mse.raw.lnHM.DEDD2.DD.overall_est.group1.all_genes, 
      lnHM.overall.2 =  mse.raw.lnHM.DEDD2.DD.overall_est.group2.all_genes
    )
  ), 
  normalised.MSEs = list(
    noDD.noDE = data.frame(
      edgeR.tag = mse.norm.edgeR.tag.DEDD2.noDD.noDE, 
      edgeR.trend = mse.norm.edgeR.trend.DEDD2.noDD.noDE, 
      DESeq2 = mse.norm.DESeq2.DEDD2.noDD.noDE, 
      DSS.trend = mse.norm.DSS.trend.DEDD2.noDD.noDE, 
      DSS.notrend = mse.norm.DSS.notrend.DEDD2.noDD.noDE, 
      expHM.1 =  mse.norm.expHM.DEDD2.noDD.group1_est.noDE, 
      expHM.2 =  mse.norm.expHM.DEDD2.noDD.group2_est.noDE, 
      expHM.overall = mse.norm.expHM.DEDD2.noDD.overall_est.noDE, 
      lnHM.1 =  mse.norm.lnHM.DEDD2.noDD.group1_est.noDE, 
      lnHM.2 =  mse.norm.lnHM.DEDD2.noDD.group2_est.noDE, 
      lnHM.overall = mse.norm.lnHM.DEDD2.noDD.overall_est.noDE
    ), 
    noDD.DE = data.frame(
      edgeR.tag = mse.norm.edgeR.tag.DEDD2.noDD.DE, 
      edgeR.trend = mse.norm.edgeR.trend.DEDD2.noDD.DE, 
      DESeq2 = mse.norm.DESeq2.DEDD2.noDD.DE, 
      DSS.trend = mse.norm.DSS.trend.DEDD2.noDD.DE, 
      DSS.notrend = mse.norm.DSS.notrend.DEDD2.noDD.DE, 
      expHM.1 =  mse.norm.expHM.DEDD2.noDD.group1_est.DE, 
      expHM.2 =  mse.norm.expHM.DEDD2.noDD.group2_est.DE, 
      expHM.overall = mse.norm.expHM.DEDD2.noDD.overall_est.DE, 
      lnHM.1 =  mse.norm.lnHM.DEDD2.noDD.group1_est.DE, 
      lnHM.2 =  mse.norm.lnHM.DEDD2.noDD.group2_est.DE, 
      lnHM.overall = mse.norm.lnHM.DEDD2.noDD.overall_est.DE
    ), 
    noDD.all_genes = data.frame(
      edgeR.tag = mse.norm.edgeR.tag.DEDD2.noDD.all_genes, 
      edgeR.trend = mse.norm.edgeR.trend.DEDD2.noDD.all_genes, 
      DESeq2 = mse.norm.DESeq2.DEDD2.noDD.all_genes, 
      DSS.trend = mse.norm.DSS.trend.DEDD2.noDD.all_genes, 
      DSS.notrend = mse.norm.DSS.notrend.DEDD2.noDD.all_genes, 
      expHM.1 =  mse.norm.expHM.DEDD2.noDD.group1_est.all_genes, 
      expHM.2 =  mse.norm.expHM.DEDD2.noDD.group2_est.all_genes, 
      expHM.overall = mse.norm.expHM.DEDD2.noDD.overall_est.all_genes, 
      lnHM.1 =  mse.norm.lnHM.DEDD2.noDD.group1_est.all_genes, 
      lnHM.2 =  mse.norm.lnHM.DEDD2.noDD.group2_est.all_genes, 
      lnHM.overall = mse.norm.lnHM.DEDD2.noDD.overall_est.all_genes
    ),
    DD.noDE = data.frame(
      edgeR.tag.1 =  mse.norm.edgeR.tag.DEDD2.DD.group1.noDE, 
      edgeR.tag.2 =  mse.norm.edgeR.tag.DEDD2.DD.group2.noDE, 
      edgeR.trend.1 =  mse.norm.edgeR.trend.DEDD2.DD.group1.noDE, 
      edgeR.trend.2 =  mse.norm.edgeR.trend.DEDD2.DD.group2.noDE, 
      DESeq2.1 =  mse.norm.DESeq2.DEDD2.DD.group1.noDE, 
      DESeq2.2 =  mse.norm.DESeq2.DEDD2.DD.group2.noDE, 
      DSS.trend.1 =  mse.norm.DSS.trend.DEDD2.DD.group1.noDE, 
      DSS.trend.2 =  mse.norm.DSS.trend.DEDD2.DD.group2.noDE, 
      DSS.notrend.1 =  mse.norm.DSS.notrend.DEDD2.DD.group1.noDE, 
      DSS.notrend.2 =  mse.norm.DSS.notrend.DEDD2.DD.group2.noDE, 
      expHM.1 =  mse.norm.expHM.DEDD2.DD.group1_est.group1.noDE, 
      expHM.2 =  mse.norm.expHM.DEDD2.DD.group2_est.group2.noDE, 
      expHM.overall.1 =  mse.norm.expHM.DEDD2.DD.overall_est.group1.noDE, 
      expHM.overall.2 =  mse.norm.expHM.DEDD2.DD.overall_est.group2.noDE, 
      lnHM.1 =  mse.norm.lnHM.DEDD2.DD.group1_est.group1.noDE, 
      lnHM.2 =  mse.norm.lnHM.DEDD2.DD.group2_est.group2.noDE, 
      lnHM.overall.1 =  mse.norm.lnHM.DEDD2.DD.overall_est.group1.noDE, 
      lnHM.overall.2 =  mse.norm.lnHM.DEDD2.DD.overall_est.group2.noDE
    ), 
    DD.DE = data.frame(
      edgeR.tag.1 =  mse.norm.edgeR.tag.DEDD2.DD.group1.DE, 
      edgeR.tag.2 =  mse.norm.edgeR.tag.DEDD2.DD.group2.DE, 
      edgeR.trend.1 =  mse.norm.edgeR.trend.DEDD2.DD.group1.DE, 
      edgeR.trend.2 =  mse.norm.edgeR.trend.DEDD2.DD.group2.DE, 
      DESeq2.1 =  mse.norm.DESeq2.DEDD2.DD.group1.DE, 
      DESeq2.2 =  mse.norm.DESeq2.DEDD2.DD.group2.DE, 
      DSS.trend.1 =  mse.norm.DSS.trend.DEDD2.DD.group1.DE, 
      DSS.trend.2 =  mse.norm.DSS.trend.DEDD2.DD.group2.DE, 
      DSS.notrend.1 =  mse.norm.DSS.notrend.DEDD2.DD.group1.DE, 
      DSS.notrend.2 =  mse.norm.DSS.notrend.DEDD2.DD.group2.DE, 
      expHM.1 =  mse.norm.expHM.DEDD2.DD.group1_est.group1.DE, 
      expHM.2 =  mse.norm.expHM.DEDD2.DD.group2_est.group2.DE, 
      expHM.overall.1 =  mse.norm.expHM.DEDD2.DD.overall_est.group1.DE, 
      expHM.overall.2 =  mse.norm.expHM.DEDD2.DD.overall_est.group2.DE, 
      lnHM.1 =  mse.norm.lnHM.DEDD2.DD.group1_est.group1.DE, 
      lnHM.2 =  mse.norm.lnHM.DEDD2.DD.group2_est.group2.DE, 
      lnHM.overall.1 =  mse.norm.lnHM.DEDD2.DD.overall_est.group1.DE, 
      lnHM.overall.2 =  mse.norm.lnHM.DEDD2.DD.overall_est.group2.DE
    ), 
    DD.all_genes = data.frame(
      edgeR.tag.1 =  mse.norm.edgeR.tag.DEDD2.DD.group1.all_genes, 
      edgeR.tag.2 =  mse.norm.edgeR.tag.DEDD2.DD.group2.all_genes, 
      edgeR.trend.1 =  mse.norm.edgeR.trend.DEDD2.DD.group1.all_genes, 
      edgeR.trend.2 =  mse.norm.edgeR.trend.DEDD2.DD.group2.all_genes, 
      DESeq2.1 =  mse.norm.DESeq2.DEDD2.DD.group1.all_genes, 
      DESeq2.2 =  mse.norm.DESeq2.DEDD2.DD.group2.all_genes, 
      DSS.trend.1 =  mse.norm.DSS.trend.DEDD2.DD.group1.all_genes, 
      DSS.trend.2 =  mse.norm.DSS.trend.DEDD2.DD.group2.all_genes, 
      DSS.notrend.1 =  mse.norm.DSS.notrend.DEDD2.DD.group1.all_genes, 
      DSS.notrend.2 =  mse.norm.DSS.notrend.DEDD2.DD.group2.all_genes, 
      expHM.1 =  mse.norm.expHM.DEDD2.DD.group1_est.group1.all_genes, 
      expHM.2 =  mse.norm.expHM.DEDD2.DD.group2_est.group2.all_genes, 
      expHM.overall.1 =  mse.norm.expHM.DEDD2.DD.overall_est.group1.all_genes, 
      expHM.overall.2 =  mse.norm.expHM.DEDD2.DD.overall_est.group2.all_genes, 
      lnHM.1 =  mse.norm.lnHM.DEDD2.DD.group1_est.group1.all_genes, 
      lnHM.2 =  mse.norm.lnHM.DEDD2.DD.group2_est.group2.all_genes, 
      lnHM.overall.1 =  mse.norm.lnHM.DEDD2.DD.overall_est.group1.all_genes, 
      lnHM.overall.2 =  mse.norm.lnHM.DEDD2.DD.overall_est.group2.all_genes
    )
  )
)


DEDD5 <- list(
  raw.MSEs = list(
    noDD.noDE = data.frame(
      edgeR.tag = mse.raw.edgeR.tag.DEDD5.noDD.noDE, 
      edgeR.trend = mse.raw.edgeR.trend.DEDD5.noDD.noDE, 
      DESeq2 = mse.raw.DESeq2.DEDD5.noDD.noDE, 
      DSS.trend = mse.raw.DSS.trend.DEDD5.noDD.noDE, 
      DSS.notrend = mse.raw.DSS.notrend.DEDD5.noDD.noDE, 
      expHM.1 =  mse.raw.expHM.DEDD5.noDD.group1_est.noDE, 
      expHM.2 =  mse.raw.expHM.DEDD5.noDD.group2_est.noDE, 
      expHM.overall = mse.raw.expHM.DEDD5.noDD.overall_est.noDE, 
      lnHM.1 =  mse.raw.lnHM.DEDD5.noDD.group1_est.noDE, 
      lnHM.2 =  mse.raw.lnHM.DEDD5.noDD.group2_est.noDE, 
      lnHM.overall = mse.raw.lnHM.DEDD5.noDD.overall_est.noDE
    ), 
    noDD.DE = data.frame(
      edgeR.tag = mse.raw.edgeR.tag.DEDD5.noDD.DE, 
      edgeR.trend = mse.raw.edgeR.trend.DEDD5.noDD.DE, 
      DESeq2 = mse.raw.DESeq2.DEDD5.noDD.DE, 
      DSS.trend = mse.raw.DSS.trend.DEDD5.noDD.DE, 
      DSS.notrend = mse.raw.DSS.notrend.DEDD5.noDD.DE, 
      expHM.1 =  mse.raw.expHM.DEDD5.noDD.group1_est.DE, 
      expHM.2 =  mse.raw.expHM.DEDD5.noDD.group2_est.DE, 
      expHM.overall = mse.raw.expHM.DEDD5.noDD.overall_est.DE, 
      lnHM.1 =  mse.raw.lnHM.DEDD5.noDD.group1_est.DE, 
      lnHM.2 =  mse.raw.lnHM.DEDD5.noDD.group2_est.DE, 
      lnHM.overall = mse.raw.lnHM.DEDD5.noDD.overall_est.DE
    ), 
    noDD.all_genes = data.frame(
      edgeR.tag = mse.raw.edgeR.tag.DEDD5.noDD.all_genes, 
      edgeR.trend = mse.raw.edgeR.trend.DEDD5.noDD.all_genes, 
      DESeq2 = mse.raw.DESeq2.DEDD5.noDD.all_genes, 
      DSS.trend = mse.raw.DSS.trend.DEDD5.noDD.all_genes, 
      DSS.notrend = mse.raw.DSS.notrend.DEDD5.noDD.all_genes, 
      expHM.1 =  mse.raw.expHM.DEDD5.noDD.group1_est.all_genes, 
      expHM.2 =  mse.raw.expHM.DEDD5.noDD.group2_est.all_genes, 
      expHM.overall = mse.raw.expHM.DEDD5.noDD.overall_est.all_genes, 
      lnHM.1 =  mse.raw.lnHM.DEDD5.noDD.group1_est.all_genes, 
      lnHM.2 =  mse.raw.lnHM.DEDD5.noDD.group2_est.all_genes, 
      lnHM.overall = mse.raw.lnHM.DEDD5.noDD.overall_est.all_genes
    ),
    DD.noDE = data.frame(
      edgeR.tag.1 =  mse.raw.edgeR.tag.DEDD5.DD.group1.noDE, 
      edgeR.tag.2 =  mse.raw.edgeR.tag.DEDD5.DD.group2.noDE, 
      edgeR.trend.1 =  mse.raw.edgeR.trend.DEDD5.DD.group1.noDE, 
      edgeR.trend.2 =  mse.raw.edgeR.trend.DEDD5.DD.group2.noDE, 
      DESeq2.1 =  mse.raw.DESeq2.DEDD5.DD.group1.noDE, 
      DESeq2.2 =  mse.raw.DESeq2.DEDD5.DD.group2.noDE, 
      DSS.trend.1 =  mse.raw.DSS.trend.DEDD5.DD.group1.noDE, 
      DSS.trend.2 =  mse.raw.DSS.trend.DEDD5.DD.group2.noDE, 
      DSS.notrend.1 =  mse.raw.DSS.notrend.DEDD5.DD.group1.noDE, 
      DSS.notrend.2 =  mse.raw.DSS.notrend.DEDD5.DD.group2.noDE, 
      expHM.1 =  mse.raw.expHM.DEDD5.DD.group1_est.group1.noDE, 
      expHM.2 =  mse.raw.expHM.DEDD5.DD.group2_est.group2.noDE, 
      expHM.overall.1 =  mse.raw.expHM.DEDD5.DD.overall_est.group1.noDE, 
      expHM.overall.2 =  mse.raw.expHM.DEDD5.DD.overall_est.group2.noDE, 
      lnHM.1 =  mse.raw.lnHM.DEDD5.DD.group1_est.group1.noDE, 
      lnHM.2 =  mse.raw.lnHM.DEDD5.DD.group2_est.group2.noDE, 
      lnHM.overall.1 =  mse.raw.lnHM.DEDD5.DD.overall_est.group1.noDE, 
      lnHM.overall.2 =  mse.raw.lnHM.DEDD5.DD.overall_est.group2.noDE
    ), 
    DD.DE = data.frame(
      edgeR.tag.1 =  mse.raw.edgeR.tag.DEDD5.DD.group1.DE, 
      edgeR.tag.2 =  mse.raw.edgeR.tag.DEDD5.DD.group2.DE, 
      edgeR.trend.1 =  mse.raw.edgeR.trend.DEDD5.DD.group1.DE, 
      edgeR.trend.2 =  mse.raw.edgeR.trend.DEDD5.DD.group2.DE, 
      DESeq2.1 =  mse.raw.DESeq2.DEDD5.DD.group1.DE, 
      DESeq2.2 =  mse.raw.DESeq2.DEDD5.DD.group2.DE, 
      DSS.trend.1 =  mse.raw.DSS.trend.DEDD5.DD.group1.DE, 
      DSS.trend.2 =  mse.raw.DSS.trend.DEDD5.DD.group2.DE, 
      DSS.notrend.1 =  mse.raw.DSS.notrend.DEDD5.DD.group1.DE, 
      DSS.notrend.2 =  mse.raw.DSS.notrend.DEDD5.DD.group2.DE, 
      expHM.1 =  mse.raw.expHM.DEDD5.DD.group1_est.group1.DE, 
      expHM.2 =  mse.raw.expHM.DEDD5.DD.group2_est.group2.DE, 
      expHM.overall.1 =  mse.raw.expHM.DEDD5.DD.overall_est.group1.DE, 
      expHM.overall.2 =  mse.raw.expHM.DEDD5.DD.overall_est.group2.DE, 
      lnHM.1 =  mse.raw.lnHM.DEDD5.DD.group1_est.group1.DE, 
      lnHM.2 =  mse.raw.lnHM.DEDD5.DD.group2_est.group2.DE, 
      lnHM.overall.1 =  mse.raw.lnHM.DEDD5.DD.overall_est.group1.DE, 
      lnHM.overall.2 =  mse.raw.lnHM.DEDD5.DD.overall_est.group2.DE
    ), 
    DD.all_genes = data.frame(
      edgeR.tag.1 =  mse.raw.edgeR.tag.DEDD5.DD.group1.all_genes, 
      edgeR.tag.2 =  mse.raw.edgeR.tag.DEDD5.DD.group2.all_genes, 
      edgeR.trend.1 =  mse.raw.edgeR.trend.DEDD5.DD.group1.all_genes, 
      edgeR.trend.2 =  mse.raw.edgeR.trend.DEDD5.DD.group2.all_genes, 
      DESeq2.1 =  mse.raw.DESeq2.DEDD5.DD.group1.all_genes, 
      DESeq2.2 =  mse.raw.DESeq2.DEDD5.DD.group2.all_genes, 
      DSS.trend.1 =  mse.raw.DSS.trend.DEDD5.DD.group1.all_genes, 
      DSS.trend.2 =  mse.raw.DSS.trend.DEDD5.DD.group2.all_genes, 
      DSS.notrend.1 =  mse.raw.DSS.notrend.DEDD5.DD.group1.all_genes, 
      DSS.notrend.2 =  mse.raw.DSS.notrend.DEDD5.DD.group2.all_genes, 
      expHM.1 =  mse.raw.expHM.DEDD5.DD.group1_est.group1.all_genes, 
      expHM.2 =  mse.raw.expHM.DEDD5.DD.group2_est.group2.all_genes, 
      expHM.overall.1 =  mse.raw.expHM.DEDD5.DD.overall_est.group1.all_genes, 
      expHM.overall.2 =  mse.raw.expHM.DEDD5.DD.overall_est.group2.all_genes, 
      lnHM.1 =  mse.raw.lnHM.DEDD5.DD.group1_est.group1.all_genes, 
      lnHM.2 =  mse.raw.lnHM.DEDD5.DD.group2_est.group2.all_genes, 
      lnHM.overall.1 =  mse.raw.lnHM.DEDD5.DD.overall_est.group1.all_genes, 
      lnHM.overall.2 =  mse.raw.lnHM.DEDD5.DD.overall_est.group2.all_genes
    )
  ), 
  normalised.MSEs = list(
    noDD.noDE = data.frame(
      edgeR.tag = mse.norm.edgeR.tag.DEDD5.noDD.noDE, 
      edgeR.trend = mse.norm.edgeR.trend.DEDD5.noDD.noDE, 
      DESeq2 = mse.norm.DESeq2.DEDD5.noDD.noDE, 
      DSS.trend = mse.norm.DSS.trend.DEDD5.noDD.noDE, 
      DSS.notrend = mse.norm.DSS.notrend.DEDD5.noDD.noDE, 
      expHM.1 =  mse.norm.expHM.DEDD5.noDD.group1_est.noDE, 
      expHM.2 =  mse.norm.expHM.DEDD5.noDD.group2_est.noDE, 
      expHM.overall = mse.norm.expHM.DEDD5.noDD.overall_est.noDE, 
      lnHM.1 =  mse.norm.lnHM.DEDD5.noDD.group1_est.noDE, 
      lnHM.2 =  mse.norm.lnHM.DEDD5.noDD.group2_est.noDE, 
      lnHM.overall = mse.norm.lnHM.DEDD5.noDD.overall_est.noDE
    ), 
    noDD.DE = data.frame(
      edgeR.tag = mse.norm.edgeR.tag.DEDD5.noDD.DE, 
      edgeR.trend = mse.norm.edgeR.trend.DEDD5.noDD.DE, 
      DESeq2 = mse.norm.DESeq2.DEDD5.noDD.DE, 
      DSS.trend = mse.norm.DSS.trend.DEDD5.noDD.DE, 
      DSS.notrend = mse.norm.DSS.notrend.DEDD5.noDD.DE, 
      expHM.1 =  mse.norm.expHM.DEDD5.noDD.group1_est.DE, 
      expHM.2 =  mse.norm.expHM.DEDD5.noDD.group2_est.DE, 
      expHM.overall = mse.norm.expHM.DEDD5.noDD.overall_est.DE, 
      lnHM.1 =  mse.norm.lnHM.DEDD5.noDD.group1_est.DE, 
      lnHM.2 =  mse.norm.lnHM.DEDD5.noDD.group2_est.DE, 
      lnHM.overall = mse.norm.lnHM.DEDD5.noDD.overall_est.DE
    ), 
    noDD.all_genes = data.frame(
      edgeR.tag = mse.norm.edgeR.tag.DEDD5.noDD.all_genes, 
      edgeR.trend = mse.norm.edgeR.trend.DEDD5.noDD.all_genes, 
      DESeq2 = mse.norm.DESeq2.DEDD5.noDD.all_genes, 
      DSS.trend = mse.norm.DSS.trend.DEDD5.noDD.all_genes, 
      DSS.notrend = mse.norm.DSS.notrend.DEDD5.noDD.all_genes, 
      expHM.1 =  mse.norm.expHM.DEDD5.noDD.group1_est.all_genes, 
      expHM.2 =  mse.norm.expHM.DEDD5.noDD.group2_est.all_genes, 
      expHM.overall = mse.norm.expHM.DEDD5.noDD.overall_est.all_genes, 
      lnHM.1 =  mse.norm.lnHM.DEDD5.noDD.group1_est.all_genes, 
      lnHM.2 =  mse.norm.lnHM.DEDD5.noDD.group2_est.all_genes, 
      lnHM.overall = mse.norm.lnHM.DEDD5.noDD.overall_est.all_genes
    ),
    DD.noDE = data.frame(
      edgeR.tag.1 =  mse.norm.edgeR.tag.DEDD5.DD.group1.noDE, 
      edgeR.tag.2 =  mse.norm.edgeR.tag.DEDD5.DD.group2.noDE, 
      edgeR.trend.1 =  mse.norm.edgeR.trend.DEDD5.DD.group1.noDE, 
      edgeR.trend.2 =  mse.norm.edgeR.trend.DEDD5.DD.group2.noDE, 
      DESeq2.1 =  mse.norm.DESeq2.DEDD5.DD.group1.noDE, 
      DESeq2.2 =  mse.norm.DESeq2.DEDD5.DD.group2.noDE, 
      DSS.trend.1 =  mse.norm.DSS.trend.DEDD5.DD.group1.noDE, 
      DSS.trend.2 =  mse.norm.DSS.trend.DEDD5.DD.group2.noDE, 
      DSS.notrend.1 =  mse.norm.DSS.notrend.DEDD5.DD.group1.noDE, 
      DSS.notrend.2 =  mse.norm.DSS.notrend.DEDD5.DD.group2.noDE, 
      expHM.1 =  mse.norm.expHM.DEDD5.DD.group1_est.group1.noDE, 
      expHM.2 =  mse.norm.expHM.DEDD5.DD.group2_est.group2.noDE, 
      expHM.overall.1 =  mse.norm.expHM.DEDD5.DD.overall_est.group1.noDE, 
      expHM.overall.2 =  mse.norm.expHM.DEDD5.DD.overall_est.group2.noDE, 
      lnHM.1 =  mse.norm.lnHM.DEDD5.DD.group1_est.group1.noDE, 
      lnHM.2 =  mse.norm.lnHM.DEDD5.DD.group2_est.group2.noDE, 
      lnHM.overall.1 =  mse.norm.lnHM.DEDD5.DD.overall_est.group1.noDE, 
      lnHM.overall.2 =  mse.norm.lnHM.DEDD5.DD.overall_est.group2.noDE
    ), 
    DD.DE = data.frame(
      edgeR.tag.1 =  mse.norm.edgeR.tag.DEDD5.DD.group1.DE, 
      edgeR.tag.2 =  mse.norm.edgeR.tag.DEDD5.DD.group2.DE, 
      edgeR.trend.1 =  mse.norm.edgeR.trend.DEDD5.DD.group1.DE, 
      edgeR.trend.2 =  mse.norm.edgeR.trend.DEDD5.DD.group2.DE, 
      DESeq2.1 =  mse.norm.DESeq2.DEDD5.DD.group1.DE, 
      DESeq2.2 =  mse.norm.DESeq2.DEDD5.DD.group2.DE, 
      DSS.trend.1 =  mse.norm.DSS.trend.DEDD5.DD.group1.DE, 
      DSS.trend.2 =  mse.norm.DSS.trend.DEDD5.DD.group2.DE, 
      DSS.notrend.1 =  mse.norm.DSS.notrend.DEDD5.DD.group1.DE, 
      DSS.notrend.2 =  mse.norm.DSS.notrend.DEDD5.DD.group2.DE, 
      expHM.1 =  mse.norm.expHM.DEDD5.DD.group1_est.group1.DE, 
      expHM.2 =  mse.norm.expHM.DEDD5.DD.group2_est.group2.DE, 
      expHM.overall.1 =  mse.norm.expHM.DEDD5.DD.overall_est.group1.DE, 
      expHM.overall.2 =  mse.norm.expHM.DEDD5.DD.overall_est.group2.DE, 
      lnHM.1 =  mse.norm.lnHM.DEDD5.DD.group1_est.group1.DE, 
      lnHM.2 =  mse.norm.lnHM.DEDD5.DD.group2_est.group2.DE, 
      lnHM.overall.1 =  mse.norm.lnHM.DEDD5.DD.overall_est.group1.DE, 
      lnHM.overall.2 =  mse.norm.lnHM.DEDD5.DD.overall_est.group2.DE
    ), 
    DD.all_genes = data.frame(
      edgeR.tag.1 =  mse.norm.edgeR.tag.DEDD5.DD.group1.all_genes, 
      edgeR.tag.2 =  mse.norm.edgeR.tag.DEDD5.DD.group2.all_genes, 
      edgeR.trend.1 =  mse.norm.edgeR.trend.DEDD5.DD.group1.all_genes, 
      edgeR.trend.2 =  mse.norm.edgeR.trend.DEDD5.DD.group2.all_genes, 
      DESeq2.1 =  mse.norm.DESeq2.DEDD5.DD.group1.all_genes, 
      DESeq2.2 =  mse.norm.DESeq2.DEDD5.DD.group2.all_genes, 
      DSS.trend.1 =  mse.norm.DSS.trend.DEDD5.DD.group1.all_genes, 
      DSS.trend.2 =  mse.norm.DSS.trend.DEDD5.DD.group2.all_genes, 
      DSS.notrend.1 =  mse.norm.DSS.notrend.DEDD5.DD.group1.all_genes, 
      DSS.notrend.2 =  mse.norm.DSS.notrend.DEDD5.DD.group2.all_genes, 
      expHM.1 =  mse.norm.expHM.DEDD5.DD.group1_est.group1.all_genes, 
      expHM.2 =  mse.norm.expHM.DEDD5.DD.group2_est.group2.all_genes, 
      expHM.overall.1 =  mse.norm.expHM.DEDD5.DD.overall_est.group1.all_genes, 
      expHM.overall.2 =  mse.norm.expHM.DEDD5.DD.overall_est.group2.all_genes, 
      lnHM.1 =  mse.norm.lnHM.DEDD5.DD.group1_est.group1.all_genes, 
      lnHM.2 =  mse.norm.lnHM.DEDD5.DD.group2_est.group2.all_genes, 
      lnHM.overall.1 =  mse.norm.lnHM.DEDD5.DD.overall_est.group1.all_genes, 
      lnHM.overall.2 =  mse.norm.lnHM.DEDD5.DD.overall_est.group2.all_genes
    )
  )
)


DEDD10 <- list(
  raw.MSEs = list(
    noDD.noDE = data.frame(
      edgeR.tag = mse.raw.edgeR.tag.DEDD10.noDD.noDE, 
      edgeR.trend = mse.raw.edgeR.trend.DEDD10.noDD.noDE, 
      DESeq2 = mse.raw.DESeq2.DEDD10.noDD.noDE, 
      DSS.notrend = mse.raw.DSS.notrend.DEDD10.noDD.noDE, 
      expHM.1 =  mse.raw.expHM.DEDD10.noDD.group1_est.noDE, 
      expHM.2 =  mse.raw.expHM.DEDD10.noDD.group2_est.noDE, 
      expHM.overall = mse.raw.expHM.DEDD10.noDD.overall_est.noDE, 
      lnHM.1 =  mse.raw.lnHM.DEDD10.noDD.group1_est.noDE, 
      lnHM.2 =  mse.raw.lnHM.DEDD10.noDD.group2_est.noDE, 
      lnHM.overall = mse.raw.lnHM.DEDD10.noDD.overall_est.noDE
    ), 
    noDD.DE = data.frame(
      edgeR.tag = mse.raw.edgeR.tag.DEDD10.noDD.DE, 
      edgeR.trend = mse.raw.edgeR.trend.DEDD10.noDD.DE, 
      DESeq2 = mse.raw.DESeq2.DEDD10.noDD.DE, 
      DSS.notrend = mse.raw.DSS.notrend.DEDD10.noDD.DE, 
      expHM.1 =  mse.raw.expHM.DEDD10.noDD.group1_est.DE, 
      expHM.2 =  mse.raw.expHM.DEDD10.noDD.group2_est.DE, 
      expHM.overall = mse.raw.expHM.DEDD10.noDD.overall_est.DE, 
      lnHM.1 =  mse.raw.lnHM.DEDD10.noDD.group1_est.DE, 
      lnHM.2 =  mse.raw.lnHM.DEDD10.noDD.group2_est.DE, 
      lnHM.overall = mse.raw.lnHM.DEDD10.noDD.overall_est.DE
    ), 
    noDD.all_genes = data.frame(
      edgeR.tag = mse.raw.edgeR.tag.DEDD10.noDD.all_genes, 
      edgeR.trend = mse.raw.edgeR.trend.DEDD10.noDD.all_genes, 
      DESeq2 = mse.raw.DESeq2.DEDD10.noDD.all_genes, 
      DSS.notrend = mse.raw.DSS.notrend.DEDD10.noDD.all_genes, 
      expHM.1 =  mse.raw.expHM.DEDD10.noDD.group1_est.all_genes, 
      expHM.2 =  mse.raw.expHM.DEDD10.noDD.group2_est.all_genes, 
      expHM.overall = mse.raw.expHM.DEDD10.noDD.overall_est.all_genes, 
      lnHM.1 =  mse.raw.lnHM.DEDD10.noDD.group1_est.all_genes, 
      lnHM.2 =  mse.raw.lnHM.DEDD10.noDD.group2_est.all_genes, 
      lnHM.overall = mse.raw.lnHM.DEDD10.noDD.overall_est.all_genes
    ),
    DD.noDE = data.frame(
      edgeR.tag.1 =  mse.raw.edgeR.tag.DEDD10.DD.group1.noDE, 
      edgeR.tag.2 =  mse.raw.edgeR.tag.DEDD10.DD.group2.noDE, 
      edgeR.trend.1 =  mse.raw.edgeR.trend.DEDD10.DD.group1.noDE, 
      edgeR.trend.2 =  mse.raw.edgeR.trend.DEDD10.DD.group2.noDE, 
      DESeq2.1 =  mse.raw.DESeq2.DEDD10.DD.group1.noDE, 
      DESeq2.2 =  mse.raw.DESeq2.DEDD10.DD.group2.noDE, 
      DSS.notrend.1 =  mse.raw.DSS.notrend.DEDD10.DD.group1.noDE, 
      DSS.notrend.2 =  mse.raw.DSS.notrend.DEDD10.DD.group2.noDE, 
      expHM.1 =  mse.raw.expHM.DEDD10.DD.group1_est.group1.noDE, 
      expHM.2 =  mse.raw.expHM.DEDD10.DD.group2_est.group2.noDE, 
      expHM.overall.1 =  mse.raw.expHM.DEDD10.DD.overall_est.group1.noDE, 
      expHM.overall.2 =  mse.raw.expHM.DEDD10.DD.overall_est.group2.noDE, 
      lnHM.1 =  mse.raw.lnHM.DEDD10.DD.group1_est.group1.noDE, 
      lnHM.2 =  mse.raw.lnHM.DEDD10.DD.group2_est.group2.noDE, 
      lnHM.overall.1 =  mse.raw.lnHM.DEDD10.DD.overall_est.group1.noDE, 
      lnHM.overall.2 =  mse.raw.lnHM.DEDD10.DD.overall_est.group2.noDE
    ), 
    DD.DE = data.frame(
      edgeR.tag.1 =  mse.raw.edgeR.tag.DEDD10.DD.group1.DE, 
      edgeR.tag.2 =  mse.raw.edgeR.tag.DEDD10.DD.group2.DE, 
      edgeR.trend.1 =  mse.raw.edgeR.trend.DEDD10.DD.group1.DE, 
      edgeR.trend.2 =  mse.raw.edgeR.trend.DEDD10.DD.group2.DE, 
      DESeq2.1 =  mse.raw.DESeq2.DEDD10.DD.group1.DE, 
      DESeq2.2 =  mse.raw.DESeq2.DEDD10.DD.group2.DE, 
      DSS.notrend.1 =  mse.raw.DSS.notrend.DEDD10.DD.group1.DE, 
      DSS.notrend.2 =  mse.raw.DSS.notrend.DEDD10.DD.group2.DE, 
      expHM.1 =  mse.raw.expHM.DEDD10.DD.group1_est.group1.DE, 
      expHM.2 =  mse.raw.expHM.DEDD10.DD.group2_est.group2.DE, 
      expHM.overall.1 =  mse.raw.expHM.DEDD10.DD.overall_est.group1.DE, 
      expHM.overall.2 =  mse.raw.expHM.DEDD10.DD.overall_est.group2.DE, 
      lnHM.1 =  mse.raw.lnHM.DEDD10.DD.group1_est.group1.DE, 
      lnHM.2 =  mse.raw.lnHM.DEDD10.DD.group2_est.group2.DE, 
      lnHM.overall.1 =  mse.raw.lnHM.DEDD10.DD.overall_est.group1.DE, 
      lnHM.overall.2 =  mse.raw.lnHM.DEDD10.DD.overall_est.group2.DE
    ), 
    DD.all_genes = data.frame(
      edgeR.tag.1 =  mse.raw.edgeR.tag.DEDD10.DD.group1.all_genes, 
      edgeR.tag.2 =  mse.raw.edgeR.tag.DEDD10.DD.group2.all_genes, 
      edgeR.trend.1 =  mse.raw.edgeR.trend.DEDD10.DD.group1.all_genes, 
      edgeR.trend.2 =  mse.raw.edgeR.trend.DEDD10.DD.group2.all_genes, 
      DESeq2.1 =  mse.raw.DESeq2.DEDD10.DD.group1.all_genes, 
      DESeq2.2 =  mse.raw.DESeq2.DEDD10.DD.group2.all_genes, 
      DSS.notrend.1 =  mse.raw.DSS.notrend.DEDD10.DD.group1.all_genes, 
      DSS.notrend.2 =  mse.raw.DSS.notrend.DEDD10.DD.group2.all_genes, 
      expHM.1 =  mse.raw.expHM.DEDD10.DD.group1_est.group1.all_genes, 
      expHM.2 =  mse.raw.expHM.DEDD10.DD.group2_est.group2.all_genes, 
      expHM.overall.1 =  mse.raw.expHM.DEDD10.DD.overall_est.group1.all_genes, 
      expHM.overall.2 =  mse.raw.expHM.DEDD10.DD.overall_est.group2.all_genes, 
      lnHM.1 =  mse.raw.lnHM.DEDD10.DD.group1_est.group1.all_genes, 
      lnHM.2 =  mse.raw.lnHM.DEDD10.DD.group2_est.group2.all_genes, 
      lnHM.overall.1 =  mse.raw.lnHM.DEDD10.DD.overall_est.group1.all_genes, 
      lnHM.overall.2 =  mse.raw.lnHM.DEDD10.DD.overall_est.group2.all_genes
    )
  ), 
  normalised.MSEs = list(
    noDD.noDE = data.frame(
      edgeR.tag = mse.norm.edgeR.tag.DEDD10.noDD.noDE, 
      edgeR.trend = mse.norm.edgeR.trend.DEDD10.noDD.noDE, 
      DESeq2 = mse.norm.DESeq2.DEDD10.noDD.noDE, 
      DSS.notrend = mse.norm.DSS.notrend.DEDD10.noDD.noDE, 
      expHM.1 =  mse.norm.expHM.DEDD10.noDD.group1_est.noDE, 
      expHM.2 =  mse.norm.expHM.DEDD10.noDD.group2_est.noDE, 
      expHM.overall = mse.norm.expHM.DEDD10.noDD.overall_est.noDE, 
      lnHM.1 =  mse.norm.lnHM.DEDD10.noDD.group1_est.noDE, 
      lnHM.2 =  mse.norm.lnHM.DEDD10.noDD.group2_est.noDE, 
      lnHM.overall = mse.norm.lnHM.DEDD10.noDD.overall_est.noDE
    ), 
    noDD.DE = data.frame(
      edgeR.tag = mse.norm.edgeR.tag.DEDD10.noDD.DE, 
      edgeR.trend = mse.norm.edgeR.trend.DEDD10.noDD.DE, 
      DESeq2 = mse.norm.DESeq2.DEDD10.noDD.DE, 
      DSS.notrend = mse.norm.DSS.notrend.DEDD10.noDD.DE, 
      expHM.1 =  mse.norm.expHM.DEDD10.noDD.group1_est.DE, 
      expHM.2 =  mse.norm.expHM.DEDD10.noDD.group2_est.DE, 
      expHM.overall = mse.norm.expHM.DEDD10.noDD.overall_est.DE, 
      lnHM.1 =  mse.norm.lnHM.DEDD10.noDD.group1_est.DE, 
      lnHM.2 =  mse.norm.lnHM.DEDD10.noDD.group2_est.DE, 
      lnHM.overall = mse.norm.lnHM.DEDD10.noDD.overall_est.DE
    ), 
    noDD.all_genes = data.frame(
      edgeR.tag = mse.norm.edgeR.tag.DEDD10.noDD.all_genes, 
      edgeR.trend = mse.norm.edgeR.trend.DEDD10.noDD.all_genes, 
      DESeq2 = mse.norm.DESeq2.DEDD10.noDD.all_genes, 
      DSS.notrend = mse.norm.DSS.notrend.DEDD10.noDD.all_genes, 
      expHM.1 =  mse.norm.expHM.DEDD10.noDD.group1_est.all_genes, 
      expHM.2 =  mse.norm.expHM.DEDD10.noDD.group2_est.all_genes, 
      expHM.overall = mse.norm.expHM.DEDD10.noDD.overall_est.all_genes, 
      lnHM.1 =  mse.norm.lnHM.DEDD10.noDD.group1_est.all_genes, 
      lnHM.2 =  mse.norm.lnHM.DEDD10.noDD.group2_est.all_genes, 
      lnHM.overall = mse.norm.lnHM.DEDD10.noDD.overall_est.all_genes
    ),
    DD.noDE = data.frame(
      edgeR.tag.1 =  mse.norm.edgeR.tag.DEDD10.DD.group1.noDE, 
      edgeR.tag.2 =  mse.norm.edgeR.tag.DEDD10.DD.group2.noDE, 
      edgeR.trend.1 =  mse.norm.edgeR.trend.DEDD10.DD.group1.noDE, 
      edgeR.trend.2 =  mse.norm.edgeR.trend.DEDD10.DD.group2.noDE, 
      DESeq2.1 =  mse.norm.DESeq2.DEDD10.DD.group1.noDE, 
      DESeq2.2 =  mse.norm.DESeq2.DEDD10.DD.group2.noDE, 
      DSS.notrend.1 =  mse.norm.DSS.notrend.DEDD10.DD.group1.noDE, 
      DSS.notrend.2 =  mse.norm.DSS.notrend.DEDD10.DD.group2.noDE, 
      expHM.1 =  mse.norm.expHM.DEDD10.DD.group1_est.group1.noDE, 
      expHM.2 =  mse.norm.expHM.DEDD10.DD.group2_est.group2.noDE, 
      expHM.overall.1 =  mse.norm.expHM.DEDD10.DD.overall_est.group1.noDE, 
      expHM.overall.2 =  mse.norm.expHM.DEDD10.DD.overall_est.group2.noDE, 
      lnHM.1 =  mse.norm.lnHM.DEDD10.DD.group1_est.group1.noDE, 
      lnHM.2 =  mse.norm.lnHM.DEDD10.DD.group2_est.group2.noDE, 
      lnHM.overall.1 =  mse.norm.lnHM.DEDD10.DD.overall_est.group1.noDE, 
      lnHM.overall.2 =  mse.norm.lnHM.DEDD10.DD.overall_est.group2.noDE
    ), 
    DD.DE = data.frame(
      edgeR.tag.1 =  mse.norm.edgeR.tag.DEDD10.DD.group1.DE, 
      edgeR.tag.2 =  mse.norm.edgeR.tag.DEDD10.DD.group2.DE, 
      edgeR.trend.1 =  mse.norm.edgeR.trend.DEDD10.DD.group1.DE, 
      edgeR.trend.2 =  mse.norm.edgeR.trend.DEDD10.DD.group2.DE, 
      DESeq2.1 =  mse.norm.DESeq2.DEDD10.DD.group1.DE, 
      DESeq2.2 =  mse.norm.DESeq2.DEDD10.DD.group2.DE, 
      DSS.notrend.1 =  mse.norm.DSS.notrend.DEDD10.DD.group1.DE, 
      DSS.notrend.2 =  mse.norm.DSS.notrend.DEDD10.DD.group2.DE, 
      expHM.1 =  mse.norm.expHM.DEDD10.DD.group1_est.group1.DE, 
      expHM.2 =  mse.norm.expHM.DEDD10.DD.group2_est.group2.DE, 
      expHM.overall.1 =  mse.norm.expHM.DEDD10.DD.overall_est.group1.DE, 
      expHM.overall.2 =  mse.norm.expHM.DEDD10.DD.overall_est.group2.DE, 
      lnHM.1 =  mse.norm.lnHM.DEDD10.DD.group1_est.group1.DE, 
      lnHM.2 =  mse.norm.lnHM.DEDD10.DD.group2_est.group2.DE, 
      lnHM.overall.1 =  mse.norm.lnHM.DEDD10.DD.overall_est.group1.DE, 
      lnHM.overall.2 =  mse.norm.lnHM.DEDD10.DD.overall_est.group2.DE
    ), 
    DD.all_genes = data.frame(
      edgeR.tag.1 =  mse.norm.edgeR.tag.DEDD10.DD.group1.all_genes, 
      edgeR.tag.2 =  mse.norm.edgeR.tag.DEDD10.DD.group2.all_genes, 
      edgeR.trend.1 =  mse.norm.edgeR.trend.DEDD10.DD.group1.all_genes, 
      edgeR.trend.2 =  mse.norm.edgeR.trend.DEDD10.DD.group2.all_genes, 
      DESeq2.1 =  mse.norm.DESeq2.DEDD10.DD.group1.all_genes, 
      DESeq2.2 =  mse.norm.DESeq2.DEDD10.DD.group2.all_genes, 
      DSS.notrend.1 =  mse.norm.DSS.notrend.DEDD10.DD.group1.all_genes, 
      DSS.notrend.2 =  mse.norm.DSS.notrend.DEDD10.DD.group2.all_genes, 
      expHM.1 =  mse.norm.expHM.DEDD10.DD.group1_est.group1.all_genes, 
      expHM.2 =  mse.norm.expHM.DEDD10.DD.group2_est.group2.all_genes, 
      expHM.overall.1 =  mse.norm.expHM.DEDD10.DD.overall_est.group1.all_genes, 
      expHM.overall.2 =  mse.norm.expHM.DEDD10.DD.overall_est.group2.all_genes, 
      lnHM.1 =  mse.norm.lnHM.DEDD10.DD.group1_est.group1.all_genes, 
      lnHM.2 =  mse.norm.lnHM.DEDD10.DD.group2_est.group2.all_genes, 
      lnHM.overall.1 =  mse.norm.lnHM.DEDD10.DD.overall_est.group1.all_genes, 
      lnHM.overall.2 =  mse.norm.lnHM.DEDD10.DD.overall_est.group2.all_genes
    )
  )
)

saveRDS(DEDD2, file=here('Results/Dispersion estimation results Aug 2019','mse.disp.DEDD2.rds'))
saveRDS(DEDD5, file=here('Results/Dispersion estimation results Aug 2019','mse.disp.DEDD5.rds'))
saveRDS(DEDD10, file=here('Results/Dispersion estimation results Aug 2019','mse.disp.DEDD10.rds'))




# Looking at some results for last analysis (for convenience, since these values are 
# still in workspace) - DEDD10.50
c(mean((edgeR.tag - true.disps.group1)^2), mean((edgeR.tag - true.disps.group2)^2))
# [1] 0.02715370 0.03261726
# overall higher MSE for group 2, but not by much
c(mean((edgeR.tag - true.disps.group1)[which(DD == 0)]^2), 
  mean((edgeR.tag - true.disps.group2)[which(DD == 0)]^2))
# [1] 0.02280115 0.02280115
# identical for group 1 genes, trivially since estimates and true values are same for 
# each group for these genes
c(mean((edgeR.tag - true.disps.group1)[which(DD == 1)]^2), 
  mean((edgeR.tag - true.disps.group2)[which(DD == 1)]^2))
# [1] 0.06596018 0.12013575
# much higher MSE for group 2 for DD genes
c(mean((edgeR.tag - true.disps.group1)[which(DD == 1 & true.disps.group1 > 0.01 & 
                                               true.disps.group1 < 1)]^2), 
  mean((edgeR.tag - true.disps.group2)[which(DD == 1 & true.disps.group2 > 0.01 & 
                                               true.disps.group2 < 1)]^2))
# [1] 0.03619965 0.03094867
# comparable MSEs (in fact lower in group 2) for DE genes restricted to the same range 
# of true dispersions for each group.
# Demonstration that poorer results for dispersion estimation in group 2 are 
# caused by the distribution of dispersions compared to group 1 - more extreme 
# low and high values.


between <- function(x, min, max) {
  return(x > min & x <= max)
}
c(mean(between(true.disps.group1, exp(-15), exp(-3))), 
  mean(between(true.disps.group1, exp(-3), exp(-2.5))), 
  mean(between(true.disps.group1, exp(-2.5), exp(-2))), 
  mean(between(true.disps.group1, exp(-2), exp(-1.5))), 
  mean(between(true.disps.group1, exp(-1.5), exp(-1))), 
  mean(between(true.disps.group1, exp(-1), exp(0))), 
  mean(between(true.disps.group1, exp(0), exp(2))))
c(mean(between(true.disps.group2, exp(-15), exp(-3))), 
  mean(between(true.disps.group2, exp(-3), exp(-2.5))), 
  mean(between(true.disps.group2, exp(-2.5), exp(-2))), 
  mean(between(true.disps.group2, exp(-2), exp(-1.5))), 
  mean(between(true.disps.group2, exp(-1.5), exp(-1))), 
  mean(between(true.disps.group2, exp(-1), exp(0))), 
  mean(between(true.disps.group2, exp(0), exp(2))))
g1.1 <- which(true.disps.group1 < exp(-3) & DD = 0)
g1.2 <- which(between(true.disps.group1, exp(-3), exp(-2.5)) & DD == 0)
g1.3 <- which(between(true.disps.group1, exp(-2.5), exp(-2)) & DD == 0)
g1.4 <- which(between(true.disps.group1, exp(-2), exp(-1.5)) & DD == 0)
g1.5 <- which(between(true.disps.group1, exp(-1.5), exp(-1)) & DD == 0)
g1.6 <- which(between(true.disps.group1, exp(-1), exp(0)) & DD == 0)
g1.7 <- which(true.disps.group1 >= exp(0) & DD == 0)
g2.1 <- which(true.disps.group2 < exp(-3) & DD == 0)
g2.2 <- which(between(true.disps.group2, exp(-3), exp(-2.5)) & DD == 0)
g2.3 <- which(between(true.disps.group2, exp(-2.5), exp(-2)) & DD == 0)
g2.4 <- which(between(true.disps.group2, exp(-2), exp(-1.5)) & DD == 0)
g2.5 <- which(between(true.disps.group2, exp(-1.5), exp(-1)) & DD == 0)
g2.6 <- which(between(true.disps.group2, exp(-1), exp(0)) & DD == 0)
g2.7 <- which(true.disps.group2 >= exp(0) & DD == 0)

rbind(c(mean((lnHM.overall_est[g1.1] - true.disps.group1[g1.1])^2), 
        mean((lnHM.overall_est[g1.2] - true.disps.group1[g1.2])^2), 
        mean((lnHM.overall_est[g1.3] - true.disps.group1[g1.3])^2), 
        mean((lnHM.overall_est[g1.4] - true.disps.group1[g1.4])^2), 
        mean((lnHM.overall_est[g1.5] - true.disps.group1[g1.5])^2), 
        mean((lnHM.overall_est[g1.6] - true.disps.group1[g1.6])^2), 
        mean((lnHM.overall_est[g1.7] - true.disps.group1[g1.7])^2)), 
      c(mean((lnHM.overall_est[g2.1] - true.disps.group2[g2.1])^2), 
        mean((lnHM.overall_est[g2.2] - true.disps.group2[g2.2])^2), 
        mean((lnHM.overall_est[g2.3] - true.disps.group2[g2.3])^2), 
        mean((lnHM.overall_est[g2.4] - true.disps.group2[g2.4])^2), 
        mean((lnHM.overall_est[g2.5] - true.disps.group2[g2.5])^2), 
        mean((lnHM.overall_est[g2.6] - true.disps.group2[g2.6])^2), 
        mean((lnHM.overall_est[g2.7] - true.disps.group2[g2.7])^2)), 
      c(mean((edgeR.tag[g1.1] - true.disps.group1[g1.1])^2), 
        mean((edgeR.tag[g1.2] - true.disps.group1[g1.2])^2), 
        mean((edgeR.tag[g1.3] - true.disps.group1[g1.3])^2), 
        mean((edgeR.tag[g1.4] - true.disps.group1[g1.4])^2), 
        mean((edgeR.tag[g1.5] - true.disps.group1[g1.5])^2), 
        mean((edgeR.tag[g1.6] - true.disps.group1[g1.6])^2), 
        mean((edgeR.tag[g1.7] - true.disps.group1[g1.7])^2)), 
      c(mean((edgeR.tag[g2.1] - true.disps.group2[g2.1])^2), 
        mean((edgeR.tag[g2.2] - true.disps.group2[g2.2])^2), 
        mean((edgeR.tag[g2.3] - true.disps.group2[g2.3])^2), 
        mean((edgeR.tag[g2.4] - true.disps.group2[g2.4])^2), 
        mean((edgeR.tag[g2.5] - true.disps.group2[g2.5])^2), 
        mean((edgeR.tag[g2.6] - true.disps.group2[g2.6])^2), 
        mean((edgeR.tag[g2.7] - true.disps.group2[g2.7])^2)))
# lnHM better than edgeR only for largest few dispersions, and much worse
# for most of range

rbind(c(mean(edgeR.tag[which(DD==0)] - mean(true.disps.group1[which(DD == 0)])), 
  mean(edgeR.trend[which(DD==0)] - mean(true.disps.group1[which(DD == 0)])), 
  mean(DESeq2[which(DD==0)] - mean(true.disps.group1[which(DD == 0)])), 
  mean(DSS.notrend[which(DD==0)] - mean(true.disps.group1[which(DD == 0)])), 
  mean(expHM.overall_est[which(DD==0)] - mean(true.disps.group1[which(DD == 0)])), 
  mean(lnHM.overall_est[which(DD==0)]) - mean(true.disps.group1[which(DD == 0)])), 
c(mean((edgeR.tag - true.disps.group1)[which(DD == 0)]^2), 
  mean((edgeR.trend - true.disps.group1)[which(DD == 0)]^2), 
  mean((DESeq2 - true.disps.group1)[which(DD == 0)]^2), 
  mean((DSS.notrend - true.disps.group1)[which(DD == 0)]^2), 
  mean((expHM.overall_est - true.disps.group1)[which(DD == 0)]^2), 
  mean((lnHM.overall_est - true.disps.group1)[which(DD == 0)]^2)))
mean(true.disps.group1[which(DD == 0)])
# approximate relationship between mean of estimates for each method and MSE

plot(log(true.disps.group1)[which(DD==0)], 
     abs(lnHM.overall_est - true.disps.group1)[which(DD==0)]^2, pch=20, col='red', 
     xlim=c(-5,1), ylim=c(0,0.8))
plot(log(true.disps.group1)[which(DD==0)], 
     abs(edgeR.tag - true.disps.group1)[which(DD==0)]^2, pch=20, col='blue', 
     xlim=c(-5,1), ylim=c(0,0.8))
plot(log(true.disps.group1)[which(DD==0)], 
     abs(DESeq2 - true.disps.group1)[which(DD==0)]^2, pch=20, col='yellow', 
     xlim=c(-5,1), ylim=c(0,0.8))
plot(log(true.disps.group1)[which(DD == 0)], 
     abs(lnHM.overall_est - true.disps.group1)[which(DD==0)] - 
       abs(edgeR.tag - true.disps.group1)[which(DD==0)], pch=20)
plot(log(true.disps.group1)[which(DD == 0)], 
     abs(lnHM.overall_est - true.disps.group1)[which(DD==0)]^2 - 
       abs(edgeR.tag - true.disps.group1)[which(DD==0)]^2, pch=20)
mean(abs(lnHM.overall_est - true.disps.group1)[which(DD==0)] - 
       abs(edgeR.tag - true.disps.group1)[which(DD==0)])
mean(abs(lnHM.overall_est - true.disps.group1)[which(DD==0)]^2 - 
       abs(edgeR.tag - true.disps.group1)[which(DD==0)]^2)
mean(abs(lnHM.overall_est - true.disps.group1)[which(DD==0)] > 
       abs(edgeR.tag - true.disps.group1)[which(DD==0)])
mean(abs(lnHM.overall_est - true.disps.group1)[which(DD==0)] > 
       abs(DESeq2 - true.disps.group1)[which(DD==0)])
plot(edgeR.tag, lnHM.overall_est, pch=20)
plot(edgeR.tag[which(DD == 0)], lnHM.overall_est[which(DD == 0)], pch=20)
mean(lnHM.overall_est > edgeR.tag)
plot(log(edgeR.tag[which(DD == 0)]), log(lnHM.overall_est[which(DD == 0)]), pch=20)
plot(log(edgeR.tag[which(DD == 1)]), log(lnHM.overall_est[which(DD == 1)]), pch=20)
plot(log(edgeR.tag[which(DE == 0)]), log(lnHM.overall_est[which(DE == 0)]), pch=20)
plot(log(edgeR.tag[which(DE == 1)]), log(lnHM.overall_est[which(DE == 1)]), pch=20)
plot(log(edgeR.tag[which(DE == 1)]), log(lnHM.group1_est[which(DE == 1)]), pch=20)
c(mean(abs(lnHM.overall_est - true.disps.group1)[which(DD==0)] > 
         abs(edgeR.tag - true.disps.group1)[which(DD==0)]), 
  mean(abs(lnHM.overall_est - true.disps.group1)[which(DD==1)] > 
         abs(edgeR.tag - true.disps.group1)[which(DD==1)]), 
  mean(abs(lnHM.overall_est - true.disps.group1)[which(DE==0)] > 
         abs(edgeR.tag - true.disps.group1)[which(DE==0)]), 
  mean(abs(lnHM.overall_est - true.disps.group1)[which(DE==1)] > 
         abs(edgeR.tag - true.disps.group1)[which(DE==1)]))
c(mean(abs(lnHM.overall_est - true.disps.group1)[which(DD==0)] > 
         abs(DESeq2 - true.disps.group1)[which(DD==0)]), 
  mean(abs(lnHM.overall_est - true.disps.group1)[which(DD==1)] > 
         abs(DESeq2 - true.disps.group1)[which(DD==1)]), 
  mean(abs(lnHM.overall_est - true.disps.group1)[which(DE==0)] > 
         abs(DESeq2 - true.disps.group1)[which(DE==0)]), 
  mean(abs(lnHM.overall_est - true.disps.group1)[which(DE==1)] > 
         abs(DESeq2 - true.disps.group1)[which(DE==1)]))

plot(true.disps.group1[which(DD == 0)], 
     abs(lnHM.overall_est - true.disps.group1)[which(DD==0)] - 
       abs(edgeR.tag - true.disps.group1)[which(DD==0)], pch=20)

rbind(c('<0.1','0.1-0.2','0.2-0.3','0.3-0.4','0.4-0.5','0.5-1','1-2','>2'), 
      round(c(mean(abs(lnHM.overall_est - true.disps.group1)[
        which(DD==0 & true.disps.group1 < 0.1)] > 
          abs(edgeR.tag - true.disps.group1)[
            which(DD==0 & true.disps.group1 < 0.1)]),
        mean(abs(lnHM.overall_est - true.disps.group1)[
          which(DD==0 & 
                  between(true.disps.group1, 0.1, 0.2))] > 
            abs(edgeR.tag - true.disps.group1)[
              which(DD==0 & 
                      between(true.disps.group1, 0.1, 0.2))]),
        mean(abs(lnHM.overall_est - true.disps.group1)[
          which(DD==0 & 
                  between(true.disps.group1, 0.2, 0.3))] > 
            abs(edgeR.tag - true.disps.group1)[
              which(DD==0 & 
                      between(true.disps.group1, 0.2, 0.3))]),
        mean(abs(lnHM.overall_est - true.disps.group1)[
          which(DD==0 & 
                  between(true.disps.group1, 0.3, 0.4))] > 
            abs(edgeR.tag - true.disps.group1)[
              which(DD==0 & 
                      between(true.disps.group1, 0.3, 0.4))]),
        mean(abs(lnHM.overall_est - true.disps.group1)[
          which(DD==0 & 
                  between(true.disps.group1, 0.4, 0.5))] > 
            abs(edgeR.tag - true.disps.group1)[
              which(DD==0 & 
                      between(true.disps.group1, 0.4, 0.5))]),
        mean(abs(lnHM.overall_est - true.disps.group1)[
          which(DD==0 & 
                  between(true.disps.group1, 0.5, 1))] > 
            abs(edgeR.tag - true.disps.group1)[
              which(DD==0 & 
                      between(true.disps.group1, 0.5, 1))]),
        mean(abs(lnHM.overall_est - true.disps.group1)[
          which(DD==0 & 
                  between(true.disps.group1, 1, 2))] > 
            abs(edgeR.tag - true.disps.group1)[
              which(DD==0 & 
                      between(true.disps.group1, 1, 2))]),
        mean(abs(lnHM.overall_est - true.disps.group1)[
          which(DD==0 & true.disps.group1 > 2)] > 
            abs(edgeR.tag - true.disps.group1)[
              which(DD==0 & true.disps.group1 > 2)])), 2))
mean(log(edgeR.tag)[which(DD==0)]) 
mean(log(lnHM.overall_est)[which(DD==0)])
mean(log(true.disps.group1)[which(DD==0)])      
hist(log(edgeR.tag[which(DD==0)]))
hist(log(lnHM.overall_est[which(DD==0)]))      
plot(density(log(edgeR.tag[which(DD==0)])), col='blue', ylim=c(0,0.6))
lines(density(log(lnHM.overall_est[which(DD==0)])), col='red')
abline(v=mean(log(true.disps.group1[which(DD==0)])))
median(edgeR.tag[which(DD==0)])
median(lnHM.overall_est[which(DD==0)])
median(true.disps.group1[which(DD==0)])
lines(density(log(true.disps.group1[which(DD==0)])), col='orange')

# Most obvious pattern is that lnHM is better than edgeR for larger dispersions, 
# particularly between about 0.5 and 2, and for dispersions below about 0.2, 
# lnHM is much worse. Could be too much shrinkage, and density plots show that 
# distribution of estimates for lnHM is narrower than true distribution, while 
# edgeR is similar.
# Other main observation from density plots is that mode of lnHM estimates matches 
# mean of true values almost exactly. Could be coincidence but could possibly also 
# be related to use of mode in MCMC algorithm. Should look into this. From plots, 
# edgeR estimates don't look great as peak of density is away from true mean, but 
# this is misleading, and I may have been misled in the same way in using the mode
# in my algorithm.


