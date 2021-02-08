library(here)
library(coda)
library(HDInterval)
library(recount)
library(edgeR)

cluster <- T

samples.per.cond <- 2
group <- factor(c(rep(1, samples.per.cond), rep(2, samples.per.cond)))

for (i in 1:10) {

  if (cluster) {
    folder <- paste0("Data/recount data/GTEx/muscle_", samples.per.cond, "_samples_per_group")
  } else {
    folder <- paste0("recount data/GTEx/muscle_", samples.per.cond, "_samples_per_group")
  }
  filename <- paste0("muscle_", samples.per.cond, "_set", i, "_DEDD")
  data <- readRDS(here(folder, paste0(filename, ".rds")))
  counts <- data$counts
  DE <- data$DEindex
  DD <- data$DDindex
  DEDD <- as.numeric(DE == 1 | DD == 1)
  FC.mean <- data$FC.mean
  FC.disp <- data$FC.disp
  rm(data)
  
  libsizes <- colSums(counts)
  nf <- calcNormFactors(counts, method="TMM")
  els <- nf * libsizes
  sf <- els / exp(mean(log(libsizes)))
  norm.counts <- t(t(counts) / sf)
  rm(counts, libsizes, nf, els, sf)

  if (cluster) {
    source(here("Data/scripts", "2020-07-17_conditional_posterior_functions_exponential_hmm.R"))
    source(here("Data/scripts", "2020-07-23_exponential_hmm_one_chain_function.R"))
    source(here("Data/scripts", "2020-07-23_exponential_hmm_three_chains_function.R"))
    # source(here("Data/scripts", "2020-07-23_exponential_hmm_adaptive_proposals_three_chains_function.R"))
    source(here("Data/scripts", "2020-07-25_exponential_hmm_adaptive_proposals_three_chains_function.R")) # mcmc.list version
    source(here("Data/scripts", "2020-09-02_run_expHMM.R"))
  } else {
    source(here("scripts", "2020-07-17_conditional_posterior_functions_exponential_hmm.R"))
    source(here("scripts", "2020-07-23_exponential_hmm_one_chain_function.R"))
    source(here("scripts", "2020-07-23_exponential_hmm_three_chains_function.R"))
    # source(here("scripts", "2020-07-23_exponential_hmm_adaptive_proposals_three_chains_function.R"))
    source(here("scripts", "2020-07-25_exponential_hmm_adaptive_proposals_three_chains_function.R")) # mcmc.list version
    source(here("scripts", "2020-09-02_run_expHMM.R"))
  }
  exp.res <- expHMM(t(norm.counts), groups=group, chain.length=3500, return.raw.results=T)
  if (cluster) {
    folder <- paste0("Data/GTEx muscle long chain results Aug 2020")
  } else {
    folder <- paste0("Results/GTEx muscle long chain results Aug 2020")
  }
  filename.exp <- paste0("raw.expHM.", filename, "_1.rds")
  saveRDS(results, file=here(folder, filename.exp))
  rm(exp.res)
  
  if (cluster) {
    source(here("Data/scripts", "2020-07-17_conditional_posterior_functions_lnonential_hmm.R"))
    source(here("Data/scripts", "2020-07-23_lnonential_hmm_one_chain_function.R"))
    source(here("Data/scripts", "2020-07-23_lnonential_hmm_three_chains_function.R"))
    source(here("Data/scripts", "2020-07-25_lnonential_hmm_adaptive_proposals_three_chains_function.R"))
    source(here("Data/scripts", "2020-07-26_run_lnHMM.R"))
  } else {
    source(here("scripts", "2020-07-17_conditional_posterior_functions_lnonential_hmm.R"))
    source(here("scripts", "2020-07-23_lnonential_hmm_one_chain_function.R"))
    source(here("scripts", "2020-07-23_lnonential_hmm_three_chains_function.R"))
    source(here("scripts", "2020-07-25_lnonential_hmm_adaptive_proposals_three_chains_function.R"))
    source(here("scripts", "2020-07-26_run_lnHMM.R"))
  }
  ln.res <- lnHMM(t(norm.counts), groups=group, chain.length=3500, return.raw.results=T)
  if (cluster) {
    folder <- paste0("Data/GTEx muscle long chain results Aug 2020")
  } else {
    folder <- paste0("Results/GTEx muscle long chain results Aug 2020")
  }
  filename.ln <- paste0("raw.lnHM.", filename, "_1.rds")
  saveRDS(results, file=here(folder, filename.ln))
  rm(ln.res)
  
  rm(norm.counts)
}

filename <- paste0("sessionInfo.GTEx_muscle_", samples.per.cond, "_set", i, "_DEDD.rds")
saveRDS(sessionInfo, file=here(folder, filename))


