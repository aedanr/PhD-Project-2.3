library(gamlss)
library(edgeR)
library(here)
library(compcodeR)
library(ROCR)

for (samples.per.cond in c(2,5,10,20,50)) {
  
  group <- factor(c(rep(1, samples.per.cond), rep(2, samples.per.cond)))
  design <- model.matrix(~group)
  folder <- paste0("recount data/GTEx/blood_", samples.per.cond, "_samples_per_group")
  
  for (i in 1:10) {
    
    filename <- paste0("blood_", samples.per.cond, "_set", i, "_DEDD")
    data <- readRDS(here(folder, paste0(filename, ".rds")))
    counts <- data$counts
    DD <- data$DDindex
    nf <- calcNormFactors(counts, method="TMM")
    Counts_edgeR <- DGEList(counts=counts, norm.factors=nf, group=group)
    Counts_edgeR <- estimateDisp(Counts_edgeR, design)
    rm(data, nf)
    
    # code modified from https://github.com/Vityay/ExpVarQuant/blob/master/ExpVarQuant.R
    
    # Uncomment this code when working in Rstudio to prevent it from locking up
    # (prevents output being shown, otherwise huge number of iteration updates printed)
    sink(file = "diversionOfOutput.txt", append = TRUE)
    # use closeAllConnections() to undo not showing output
    
    Counts_edgeR$CPM <- cpm.DGEList(Counts_edgeR)
    ofs <- log(Counts_edgeR$samples$lib.size * Counts_edgeR$samples$norm.factors)
    Counts_edgeR$samples$offset <- ofs
    rm(ofs)
    
    # Estimate group effects on mean and overdispersion parameters with GAMLSS for each gene.
    gene_i <- seq_along(Counts_edgeR$counts[,1])
    gamlss_NB <- lapply(gene_i, function(i) {
      
      # For each gene (i) a table containing: x - a factor of interest (age); y - RNA-seq. counts and offset (ofs) is created.
      dat <- data.frame(
        x = Counts_edgeR$samples$group,
        y = Counts_edgeR$counts[i,],
        ofs = Counts_edgeR$samples$offset
      )
      # x is releveled to use group 1 as a reference.
      dat$x <- relevel(dat$x, ref = c("1"))
      
      # Fit negative binomial (family = NBI()) GAMLSS model, which accounts for age effects on mean and overdispersion (non-Poisson noise).
      # fo = y~0+x+offset(ofs) specifies model for mean and sigma.fo=~0+x for overdispersion, offset - offset(ofs) normalize counts to library size. sigma.start = 0.1 provides starting value for overdispersion estimation (default is 1). n.cyc – number of fitting algorithm cycles, see help(gamlss).
      # In some cases, fitting of NB model may fail and tryCatch(..., warning= function(w) NULL, error= function(e) NULL) will return NULL as a result.
      m0 <- tryCatch(
        gamlss(fo = y ~ 0+x+offset(ofs), sigma.fo = ~ 0+x, data=dat,
               family = NBI(), sigma.start = 0.1, n.cyc = 100),
        warning= function(w) NULL, error= function(e) NULL
      )
      
      # Fit reduced model by omitting age factor from the estimation of overdispersion: sigma.fo = ~ 1. In essence, this model corresponds to the GLM model implemented in edgeR.
      m1 <- tryCatch(
        gamlss(fo = y ~ 0+x+offset(ofs), sigma.fo = ~ 1, data=dat,
               family = NBI(), sigma.start = 0.1, n.cyc = 100),
        warning= function(w) NULL, error= function(e) NULL
      )
      
      # Fit reduced model by omitting age factor from the estimation of mean: fo = y ~ offset(ofs).
      m2 <- tryCatch(
        gamlss(fo = y ~ offset(ofs), sigma.fo = ~ 0+x, data=dat,
               family = NBI(), sigma.start = 0.1, n.cyc = 100),
        warning= function(w) NULL, error= function(e) NULL
      )
      
      # Fit null model.
      m3 <- tryCatch(
        gamlss(fo = y ~ offset(ofs), sigma.fo = ~ 1, data=dat,
               family = NBI(), sigma.start = 0.1, n.cyc = 100),
        warning= function(w) NULL, error= function(e) NULL
      )
      
      # Create data frame res to store the results.
      res <- data.frame(
        cpm.1 = NA,
        cpm.2 = NA,
        LR.cpm = NA,
        p_gamlss.cpm = NA,
        p_glm.cpm = NA,
        CV.1 = NA,
        CV.2 = NA,
        LR.cv = NA,
        p.cv = NA
      )
      
      # Because fitting of the NB model may fail for some genes, check whether all models were fitted successfully. 
      if(!any(sapply(list(m0,m1,m2,m3), is.null))) 
      {
        # Write GAMLSS estimations of gene’s mean (CPM) counts from the m0 model.
        res$cpm.1 = exp(m0$mu.coefficients+log(1e06))[[1]]
        res$cpm.2 = exp(m0$mu.coefficients+log(1e06))[[2]]
        
        # Calculate log2 ratio for changes in CPMs between old and young mice.
        res$LR.cpm = log2(exp(m0$mu.coefficients+log(1e06))[[2]]/exp(m0$mu.coefficients+log(1e06))[[1]])
        
        # GAMLSS log-likelihood ratio (LR) test for a significance of an age effect on gene’s mean (CPM) counts.
        # p_gamlss.cpm – p value of LR test statistic: D_μ=-2log⁡[L(μ_0,α_j  ┤|  X_ij)/L(μ_j,α_j  ┤|  X_ij), comparing m0 and m2 models.
        res$p_gamlss.cpm = pchisq(2*(logLik(m0)-logLik(m2)), df=m0$df.fit-m2$df.fit, lower=F)[[1]]
        
        # GLM log-likelihood ratio (LR) test for a significance of age effect on gene’s mean (CPM) counts.
        # p_glm.cpm – p value of LR test statistic: D_(μ_GLM )=-2log⁡[L(μ_0,α_0  ┤|  X_ij)/L(μ_j,α_0  ┤|  X_ij), comparing m1 and m3 models.
        res$p_glm.cpm = pchisq(2*(logLik(m1)-logLik(m3)), df=m1$df.fit-m3$df.fit, lower=F)[[1]]
        
        # Write GAMLSS estimations of gene’s non-Poisson noise from the m0 model: cv(μ)=√α.
        res$CV.1  = sqrt(exp(m0$sigma.coefficients)[[1]])
        res$CV.2 = sqrt(exp(m0$sigma.coefficients)[[2]])
        
        # Calculate log2 ratio for changes in cv(μ) between old and young mice.
        res$LR.cv = log2(sqrt(exp(m0$sigma.coefficients)[[2]])/sqrt(exp(m0$sigma.coefficients)[[1]]))
        
        # GAMLSS log-likelihood ratio (LR) test for a significance of an age effect on non-Poisson noise.
        # p.cv – p value of LR test statistic: D_α=-2log⁡[L(μ_j,α_0  ┤|  X_ij)/L(μ_j,α_j  ┤|  X_ij), comparing m0 and m1 models.
        res$p.cv = pchisq(2*(logLik(m0)-logLik(m1)), df=m0$df.fit-m1$df.fit, lower=F)[[1]]
      }
      res
    })
    
    closeAllConnections()
    
    # Transform list gamlss_NB containing GAMLSS estimations for each gene to data frame
    gamlss_NB <- do.call(rbind, gamlss_NB)
    rownames(gamlss_NB) <- rownames(Counts_edgeR$counts)[gene_i]
    
    # Because GAMLSS fitting might fail for some genes or yield inflated estimates of overdispersion, the results have to be cleaned.
    # First, remove genes, for which GAMLSS model has failed.
    res_clean <- na.omit(gamlss_NB)
    
    # Second, remove genes, for which estimates of cv(μ)=√α were either inflated > 3 or close to Poisson < 10-3.
    # idx <- gamlss_NB_clean$CV.1 > 3 | gamlss_NB_clean$CV.1 < 1e-03 | gamlss_NB_clean$CV.2 > 3 | gamlss_NB_clean$CV.2 < 1e-03
    # gamlss_NB_clean <- gamlss_NB_clean[!idx,]
    
    # Finally, calculate false discovery rates to account for multiple hypothesis testing.
    res_clean$padj_gamlss.cpm <- p.adjust(res_clean$p_gamlss.cpm, "fdr")
    res_clean$padj_glm.cpm <- p.adjust(res_clean$p_glm.cpm, "fdr")
    res_clean$padj.cv <- p.adjust(res_clean$p.cv, "fdr")
    
    results <- list(
      counts = counts, 
      DD = DD, 
      full_res = gamlss_NB, 
      res_no_NA = res_clean, 
      DD_no_NA = DD[which(rownames(counts) %in% rownames(res_clean))]
    )
    
    saveRDS(results, here("Results/ExpVarQuant GTEx results Dec 2020", paste0("results.ExpVarQuant", filename, ".rds")))
    
  }
  
  rm(group, design)
  
}

