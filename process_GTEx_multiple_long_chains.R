library(coda)
library(here)

method <- "lnHM"
size <- "50"
source <- paste0("Results/GTEx muscle long chain results Aug 2020/Raw data/", size, " samples per group/", method)
temp <- paste0("Results/GTEx muscle long chain results Aug 2020/Intermediate data/", size, " samples per group/", method)
data.source <- "Results/GTEx muscle simulated DE, DD, DEDD results March 2020"
comb <- paste0("Results/GTEx muscle long chain results Aug 2020/Combined chain results/", size, " samples per group/", method)
final <- "Results/GTEx muscle long chain results Aug 2020/"

## Import results and save components individually ####

# Run 1
for (set in 7:10) {
  dat <- readRDS(here(source, paste0("raw.", method, ".muscle_", size, "_set", set, "_DEDD.rds")))
  means1_run1 <- as.matrix(dat$means1)
  saveRDS(means1_run1, here(temp, paste0("means1_run1.", method, ".muscle_", size, "_set", set, "_DEDD.rds")))
  rm(means1_run1)
  means2_run1 <- as.matrix(dat$means2)
  saveRDS(means2_run1, here(temp, paste0("means2_run1.", method, ".muscle_", size, "_set", set, "_DEDD.rds")))
  rm(means2_run1)
  disps1_run1 <- as.matrix(dat$disps1)
  saveRDS(disps1_run1, here(temp, paste0("disps1_run1.", method, ".muscle_", size, "_set", set, "_DEDD.rds")))
  rm(disps1_run1)
  disps2_run1 <- as.matrix(dat$disps2)
  saveRDS(disps2_run1, here(temp, paste0("disps2_run1.", method, ".muscle_", size, "_set", set, "_DEDD.rds")))
  rm(disps2_run1)
  prob_run1 <- as.matrix(dat$indicators)
  saveRDS(prob_run1, here(temp, paste0("prob_run1.", method, ".muscle_", size, "_set", set, "_DEDD.rds")))
  rm(prob_run1)
  post.prop_run1 <- as.matrix(dat$proportion)
  saveRDS(post.prop_run1, here(temp, paste0("post.prop_run1.", method, ".muscle_", size, "_set", set, "_DEDD.rds")))
  rm(post.prop_run1)
  rm(dat)
  gc()
}

# Run 2
for (set in 7:10) {
  dat <- readRDS(here(source, paste0("raw.", method, ".muscle_", size, "_set", set, "_DEDD_2.rds")))
  means1_run2 <- as.matrix(dat$means1)
  saveRDS(means1_run2, here(temp, paste0("means1_run2.", method, ".muscle_", size, "_set", set, "_DEDD.rds")))
  rm(means1_run2)
  means2_run2 <- as.matrix(dat$means2)
  saveRDS(means2_run2, here(temp, paste0("means2_run2.", method, ".muscle_", size, "_set", set, "_DEDD.rds")))
  rm(means2_run2)
  disps1_run2 <- as.matrix(dat$disps1)
  saveRDS(disps1_run2, here(temp, paste0("disps1_run2.", method, ".muscle_", size, "_set", set, "_DEDD.rds")))
  rm(disps1_run2)
  disps2_run2 <- as.matrix(dat$disps2)
  saveRDS(disps2_run2, here(temp, paste0("disps2_run2.", method, ".muscle_", size, "_set", set, "_DEDD.rds")))
  rm(disps2_run2)
  prob_run2 <- as.matrix(dat$indicators)
  saveRDS(prob_run2, here(temp, paste0("prob_run2.", method, ".muscle_", size, "_set", set, "_DEDD.rds")))
  rm(prob_run2)
  post.prop_run2 <- as.matrix(dat$proportion)
  saveRDS(post.prop_run2, here(temp, paste0("post.prop_run2.", method, ".muscle_", size, "_set", set, "_DEDD.rds")))
  rm(post.prop_run2)
  rm(dat)
  gc()
}

# Run 3
for (set in 7:10) {
  dat <- readRDS(here(source, paste0("raw.", method, ".muscle_", size, "_set", set, "_DEDD_3.rds")))
  means1_run3 <- as.matrix(dat$means1)
  saveRDS(means1_run3, here(temp, paste0("means1_run3.", method, ".muscle_", size, "_set", set, "_DEDD.rds")))
  rm(means1_run3)
  means2_run3 <- as.matrix(dat$means2)
  saveRDS(means2_run3, here(temp, paste0("means2_run3.", method, ".muscle_", size, "_set", set, "_DEDD.rds")))
  rm(means2_run3)
  disps1_run3 <- as.matrix(dat$disps1)
  saveRDS(disps1_run3, here(temp, paste0("disps1_run3.", method, ".muscle_", size, "_set", set, "_DEDD.rds")))
  rm(disps1_run3)
  disps2_run3 <- as.matrix(dat$disps2)
  saveRDS(disps2_run3, here(temp, paste0("disps2_run3.", method, ".muscle_", size, "_set", set, "_DEDD.rds")))
  rm(disps2_run3)
  prob_run3 <- as.matrix(dat$indicators)
  saveRDS(prob_run3, here(temp, paste0("prob_run3.", method, ".muscle_", size, "_set", set, "_DEDD.rds")))
  rm(prob_run3)
  post.prop_run3 <- as.matrix(dat$proportion)
  saveRDS(post.prop_run3, here(temp, paste0("post.prop_run3.", method, ".muscle_", size, "_set", set, "_DEDD.rds")))
  rm(post.prop_run3)
  rm(dat)
  gc()
}


## Combine individual runs into overall results ####
for (set in 7:10) {
  means1_run1 <- readRDS(here(temp, paste0("means1_run1.", method, ".muscle_", size, "_set", set, "_DEDD.rds")))
  means1_run2 <- readRDS(here(temp, paste0("means1_run2.", method, ".muscle_", size, "_set", set, "_DEDD.rds")))
  means1_run3 <- readRDS(here(temp, paste0("means1_run3.", method, ".muscle_", size, "_set", set, "_DEDD.rds")))
  means1 <- rbind(means1_run1, means1_run2, means1_run3)
  rm(means1_run1, means1_run2, means1_run3)
  means2_run1 <- readRDS(here(temp, paste0("means2_run1.", method, ".muscle_", size, "_set", set, "_DEDD.rds")))
  means2_run2 <- readRDS(here(temp, paste0("means2_run2.", method, ".muscle_", size, "_set", set, "_DEDD.rds")))
  means2_run3 <- readRDS(here(temp, paste0("means2_run3.", method, ".muscle_", size, "_set", set, "_DEDD.rds")))
  means2 <- rbind(means2_run1, means2_run2, means2_run3)
  rm(means2_run1, means2_run2, means2_run3)
  lfc.mean <- log2(colMeans(means2)) - log2(colMeans(means1))
  saveRDS(lfc.mean, here(comb, paste0("lfc.mean.", method, ".muscle_", size, "_set", set, "_DEDD.rds")))
  rm(lfc.mean)
  mean.diff.untr <- means1 - means2
  saveRDS(mean.diff.untr, here(comb, paste0("mean.diff.untr.", method, ".muscle_", size, "_set", set, "_DEDD.rds")))
  rm(mean.diff.untr)
  log.means1 <- log(means1)
  rm(means1)
  log.means2 <- log(means2)
  rm(means2)
  mean.diff.log <- log.means1 - log.means2
  rm(log.means1, log.means2)
  saveRDS(mean.diff.log, here(comb, paste0("mean.diff.log.", method, ".muscle_", size, "_set", set, "_DEDD.rds")))
  rm(mean.diff.log)
  
  disps1_run1 <- readRDS(here(temp, paste0("disps1_run1.", method, ".muscle_", size, "_set", set, "_DEDD.rds")))
  disps1_run2 <- readRDS(here(temp, paste0("disps1_run2.", method, ".muscle_", size, "_set", set, "_DEDD.rds")))
  disps1_run3 <- readRDS(here(temp, paste0("disps1_run3.", method, ".muscle_", size, "_set", set, "_DEDD.rds")))
  disps1 <- rbind(disps1_run1, disps1_run2, disps1_run3)
  rm(disps1_run1, disps1_run2, disps1_run3)
  disps2_run1 <- readRDS(here(temp, paste0("disps2_run1.", method, ".muscle_", size, "_set", set, "_DEDD.rds")))
  disps2_run2 <- readRDS(here(temp, paste0("disps2_run2.", method, ".muscle_", size, "_set", set, "_DEDD.rds")))
  disps2_run3 <- readRDS(here(temp, paste0("disps2_run3.", method, ".muscle_", size, "_set", set, "_DEDD.rds")))
  disps2 <- rbind(disps2_run1, disps2_run2, disps2_run3)
  rm(disps2_run1, disps2_run2, disps2_run3)
  lfc.disp <- log2(colMeans(disps2)) - log2(colMeans(disps1))
  saveRDS(lfc.disp, here(comb, paste0("lfc.disp.", method, ".muscle_", size, "_set", set, "_DEDD.rds")))
  rm(lfc.disp)
  disp.diff.untr <- disps1 - disps2
  saveRDS(disp.diff.untr, here(comb, paste0("disp.diff.untr.", method, ".muscle_", size, "_set", set, "_DEDD.rds")))
  rm(disp.diff.untr)
  log.disps1 <- log(disps1)
  rm(disps1)
  log.disps2 <- log(disps2)
  rm(disps2)
  disp.diff.log <- log.disps1 - log.disps2
  rm(log.disps1, log.disps2)
  saveRDS(disp.diff.log, here(comb, paste0("disp.diff.log.", method, ".muscle_", size, "_set", set, "_DEDD.rds")))
  rm(disp.diff.log)
  
  prob_run1 <- readRDS(here(temp, paste0("prob_run1.", method, ".muscle_", size, "_set", set, "_DEDD.rds")))
  prob_run2 <- readRDS(here(temp, paste0("prob_run2.", method, ".muscle_", size, "_set", set, "_DEDD.rds")))
  prob_run3 <- readRDS(here(temp, paste0("prob_run3.", method, ".muscle_", size, "_set", set, "_DEDD.rds")))
  prob <- colMeans(rbind(prob_run1, prob_run2, prob_run3))
  rm(prob_run1, prob_run2, prob_run3)
  saveRDS(prob, here(comb, paste0("prob.", method, ".muscle_", size, "_set", set, "_DEDD.rds")))
  rm(prob)
  post.prop_run1 <- readRDS(here(temp, paste0("post.prop_run1.", method, ".muscle_", size, "_set", set, "_DEDD.rds")))
  post.prop_run2 <- readRDS(here(temp, paste0("post.prop_run2.", method, ".muscle_", size, "_set", set, "_DEDD.rds")))
  post.prop_run3 <- readRDS(here(temp, paste0("post.prop_run3.", method, ".muscle_", size, "_set", set, "_DEDD.rds")))
  post.prop <- mean(c(as.numeric(post.prop_run1), as.numeric(post.prop_run2), as.numeric(post.prop_run3)))
  rm(post.prop_run1, post.prop_run2, post.prop_run3)
  saveRDS(post.prop, here(comb, paste0("post.prop.", method, ".muscle_", size, "_set", set, "_DEDD.rds")))
  rm(post.prop)
  gc()
}


## Import separate overall result files and save in usual format for analysis ####
source(here("scripts", "2020-07-23_hpd_tail_prob_function.R"))
source(here("scripts", "2019-05-03_bfdr_function.R"))

for (set in 7:10) {
  dat <- readRDS(here(data.source, paste0("results.muscle_", size, "_set", set, "_DEDD.rds")))
  counts <- dat$counts
  DE <- dat$DE
  DD <- dat$DD
  DEDD <- dat$DEDD
  rm(dat)

  mean.diff.untr <- readRDS(here(comb, paste0("mean.diff.untr.", method, ".muscle_", size, "_set", set, "_DEDD.rds")))
  p.mean.untr <- apply(mean.diff.untr,2,hpd.pval)
  rm(mean.diff.untr)
  q.mean.untr <- p.adjust(p.mean.untr, method='BH')

  mean.diff.log <- readRDS(here(comb, paste0("mean.diff.log.", method, ".muscle_", size, "_set", set, "_DEDD.rds")))
  p.mean.log <- apply(mean.diff.log,2,hpd.pval)
  rm(mean.diff.log)
  q.mean.log <- p.adjust(p.mean.log, method='BH')
  
  disp.diff.untr <- readRDS(here(comb, paste0("disp.diff.untr.", method, ".muscle_", size, "_set", set, "_DEDD.rds")))
  p.disp.untr <- apply(disp.diff.untr,2,hpd.pval)
  rm(disp.diff.untr)
  q.disp.untr <- p.adjust(p.disp.untr, method='BH')
  
  disp.diff.log <- readRDS(here(comb, paste0("disp.diff.log.", method, ".muscle_", size, "_set", set, "_DEDD.rds")))
  p.disp.log <- apply(disp.diff.log,2,hpd.pval)
  rm(disp.diff.log)
  q.disp.log <- p.adjust(p.disp.log, method='BH')
  
  lfc.mean <- readRDS(here(comb, paste0("lfc.mean.", method, ".muscle_", size, "_set", set, "_DEDD.rds")))
  lfc.disp <- readRDS(here(comb, paste0("lfc.disp.", method, ".muscle_", size, "_set", set, "_DEDD.rds")))
  prob <- readRDS(here(comb, paste0("prob.", method, ".muscle_", size, "_set", set, "_DEDD.rds")))
  post.prop <- readRDS(here(comb, paste0("post.prop.", method, ".muscle_", size, "_set", set, "_DEDD.rds")))
  
  thr <- unname(sort(prob, decreasing=T)[round(length(prob) * post.prop)])
  BFDR <- bfdr(prob)
  
  res <- list(counts = counts, 
              DE = DE, 
              DD = DD, 
              DEDD = DEDD, 
              prob = prob, 
              prop = post.prop, 
              thr = thr, 
              bfdr = BFDR, 
              p.mean.untr = p.mean.untr, 
              p.mean.log = p.mean.log, 
              q.mean.untr = q.mean.untr, 
              q.mean.log = q.mean.log, 
              lfc.mean = lfc.mean, 
              p.disp.untr = p.disp.untr, 
              p.disp.log = p.disp.log, 
              q.disp.untr = q.disp.untr, 
              q.disp.log = q.disp.log, 
              lfc.disp = lfc.disp)
  
  saveRDS(res, file=here(final, paste0("results.", method, ".muscle_", size, "_set", set, "_DEDD.rds")))
  
  rm(lfc.mean, p.mean.untr, p.mean.log, q.mean.untr, q.mean.log, 
     lfc.disp, p.disp.untr, p.disp.log, q.disp.untr, q.disp.log, 
     prob, post.prop, thr, BFDR, counts, DE, DD, DEDD, res)
}


