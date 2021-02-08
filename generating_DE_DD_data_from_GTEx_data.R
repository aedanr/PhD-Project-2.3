library(here)
library(recount)
source(here('scripts','2019-11-29_functions_to_change_mean_disp_in_real_data.R'))

# Load base GTEx blood datasets (100, 20, 10 and 4 randomly chosen samples with RIN > 8.7 and 
# extracted using PAXgene Whole Blood)
blood_100 <- readRDS(file=here("recount data/GTEx", "blood_100.rds"))
blood_20 <- readRDS(file=here("recount data/GTEx", "blood_20.rds"))
blood_10 <- readRDS(file=here("recount data/GTEx", "blood_10.rds"))
blood_4 <- readRDS(file=here("recount data/GTEx", "blood_4.rds"))

# Initially just do one permutation for each sample size. Will take a long time to do lots of 
# permutations - as long as original simulation work. If doing multiple permutations for smaller 
# sample sizes, need to choose samples from a larger pool rather than permuting group labels 
# since there are not enough permutations (e.g. 45 permutations for 2 groups of 5).

# Create groups
group_100 <- factor(c(rep("1",50), rep("2",50)))
group_20 <- factor(c(rep("1",10), rep("2",10)))
group_10 <- factor(c(rep("1",5), rep("2",5)))
group_4 <- factor(c(rep("1",2), rep("2",2)))

# Create gene-level counts from base pair-level counts
blood_100 <- scale_counts(blood_100)
blood_20 <- scale_counts(blood_20)
blood_10 <- scale_counts(blood_10)
blood_4 <- scale_counts(blood_4)

# Apply minimum counts per million 0.5 filter
blood_100_filtered <- blood_100[which(rowMeans(apply(assay(blood_100), 2, function(x) x*1e6/sum(x))) > 0.5), ]
dim(blood_100_filtered) # 18563 100
blood_20_filtered <- blood_20[which(rowMeans(apply(assay(blood_20), 2, function(x) x*1e6/sum(x))) > 0.5), ]
dim(blood_20_filtered) # 18903 20
blood_10_filtered <- blood_10[which(rowMeans(apply(assay(blood_10), 2, function(x) x*1e6/sum(x))) > 0.5), ]
dim(blood_10_filtered) # 17937 10
blood_4_filtered <- blood_4[which(rowMeans(apply(assay(blood_4), 2, function(x) x*1e6/sum(x))) > 0.5), ]
dim(blood_4_filtered) # 17741 4

## Create initial DE/DD indices
# DE only - 10% genes DE
DEindex_DE100 <- numeric(nrow(blood_100_filtered))
DEindex_DE100[1:ceiling(nrow(blood_100_filtered)/10)] <- 1
DEindex_DE20 <- numeric(nrow(blood_20_filtered))
DEindex_DE20[1:ceiling(nrow(blood_20_filtered)/10)] <- 1
DEindex_DE10 <- numeric(nrow(blood_10_filtered))
DEindex_DE10[1:ceiling(nrow(blood_10_filtered)/10)] <- 1
DEindex_DE4 <- numeric(nrow(blood_4_filtered))
DEindex_DE4[1:ceiling(nrow(blood_4_filtered)/10)] <- 1

# DD only - 10% genes DD
DDindex_DD100 <- numeric(nrow(blood_100_filtered))
DDindex_DD100[1:ceiling(nrow(blood_100_filtered)/10)] <- 1
DDindex_DD20 <- numeric(nrow(blood_20_filtered))
DDindex_DD20[1:ceiling(nrow(blood_20_filtered)/10)] <- 1
DDindex_DD10 <- numeric(nrow(blood_10_filtered))
DDindex_DD10[1:ceiling(nrow(blood_10_filtered)/10)] <- 1
DDindex_DD4 <- numeric(nrow(blood_4_filtered))
DDindex_DD4[1:ceiling(nrow(blood_4_filtered)/10)] <- 1

# DEDD - 5% DE only, 5% DD only, 5% both
DEindex_DEDD100 <- numeric(nrow(blood_100_filtered))
DEindex_DEDD100[1:ceiling(nrow(blood_100_filtered)/10)] <- 1
DEindex_DEDD20 <- numeric(nrow(blood_20_filtered))
DEindex_DEDD20[1:ceiling(nrow(blood_20_filtered)/10)] <- 1
DEindex_DEDD10 <- numeric(nrow(blood_10_filtered))
DEindex_DEDD10[1:ceiling(nrow(blood_10_filtered)/10)] <- 1
DEindex_DEDD4 <- numeric(nrow(blood_4_filtered))
DEindex_DEDD4[1:ceiling(nrow(blood_4_filtered)/10)] <- 1
DDindex_DEDD100 <- numeric(nrow(blood_100_filtered))
DDindex_DEDD100[ceiling(nrow(blood_100_filtered)/20):(3*ceiling(nrow(blood_100_filtered)/20))] <- 1
DDindex_DEDD20 <- numeric(nrow(blood_20_filtered))
DDindex_DEDD20[ceiling(nrow(blood_20_filtered)/20):(3*ceiling(nrow(blood_20_filtered)/20))] <- 1
DDindex_DEDD10 <- numeric(nrow(blood_10_filtered))
DDindex_DEDD10[ceiling(nrow(blood_10_filtered)/20):(3*ceiling(nrow(blood_10_filtered)/20))] <- 1
DDindex_DEDD4 <- numeric(nrow(blood_4_filtered))
DDindex_DEDD4[ceiling(nrow(blood_4_filtered)/20):(3*ceiling(nrow(blood_4_filtered)/20))] <- 1
DEDDindex_DEDD100 <- as.numeric(DEindex_DEDD100==1 | DDindex_DEDD100==1)
DEDDindex_DEDD20 <- as.numeric(DEindex_DEDD20==1 | DDindex_DEDD20==1)
DEDDindex_DEDD10 <- as.numeric(DEindex_DEDD10==1 | DDindex_DEDD10==1)
DEDDindex_DEDD4 <- as.numeric(DEindex_DEDD4==1 | DDindex_DEDD4==1)
DEonlyindex_DEDD100 <- as.numeric(DEindex_DEDD100==1 & DDindex_DEDD100==0)
DEonlyindex_DEDD20 <- as.numeric(DEindex_DEDD20==1 & DDindex_DEDD20==0)
DEonlyindex_DEDD10 <- as.numeric(DEindex_DEDD10==1 & DDindex_DEDD10==0)
DEonlyindex_DEDD4 <- as.numeric(DEindex_DEDD4==1 & DDindex_DEDD4==0)
DDonlyindex_DEDD100 <- as.numeric(DEindex_DEDD100==0 & DDindex_DEDD100==1)
DDonlyindex_DEDD20 <- as.numeric(DEindex_DEDD20==0 & DDindex_DEDD20==1)
DDonlyindex_DEDD10 <- as.numeric(DEindex_DEDD10==0 & DDindex_DEDD10==1)
DDonlyindex_DEDD4 <- as.numeric(DEindex_DEDD4==0 & DDindex_DEDD4==1)
DEplusDDindex_DEDD100 <- as.numeric(DEindex_DEDD100==1 & DDindex_DEDD100==1)
DEplusDDindex_DEDD20 <- as.numeric(DEindex_DEDD20==1 & DDindex_DEDD20==1)
DEplusDDindex_DEDD10 <- as.numeric(DEindex_DEDD10==1 & DDindex_DEDD10==1)
DEplusDDindex_DEDD4 <- as.numeric(DEindex_DEDD4==1 & DDindex_DEDD4==1)

## Create DE and DD data, using integers as original GTEx data looks to be integer
## Here using only samples within each dataset to estimate true means and dispersions 
## to create new means and dispersions. This will probably give really unreliable 
## data for DE and DD for small samples. Since I have a large number of samples to 
## work from, might be better to use a bigger number of samples to generate the DE/DD 
## counts and then subsample.
# DE only
DE100 <- assay(blood_100_filtered)
DE100.DE <- diff.mean(DE100[DEindex_DE100==1, group_100=="2"], integer=T)
DE100[DEindex_DE100==1, group_100=="2"] <- DE100.DE$counts
FC_DE100 <- rep(1, nrow(DE100))
FC_DE100[DEindex_DE100==1] <- DE100.DE$FC
rm(DE100.DE)
DE20 <- assay(blood_20_filtered)
DE20.DE <- diff.mean(DE20[DEindex_DE20==1, group_20=="2"], integer=T)
DE20[DEindex_DE20==1, group_20=="2"] <- DE20.DE$counts
FC_DE20 <- rep(1, nrow(DE20))
FC_DE20[DEindex_DE20==1] <- DE20.DE$FC
rm(DE20.DE)
DE10 <- assay(blood_10_filtered)
DE10.DE <- diff.mean(DE10[DEindex_DE10==1, group_10=="2"], integer=T)
DE10[DEindex_DE10==1, group_10=="2"] <- DE10.DE$counts
FC_DE10 <- rep(1, nrow(DE10))
FC_DE10[DEindex_DE10==1] <- DE10.DE$FC
rm(DE10.DE)
DE4 <- assay(blood_4_filtered)
DE4.DE <- diff.mean(DE4[DEindex_DE4==1, group_4=="2"], integer=T)
DE4[DEindex_DE4==1, group_4=="2"] <- DE4.DE$counts
FC_DE4 <- rep(1, nrow(DE4))
FC_DE4[DEindex_DE4==1] <- DE4.DE$FC
rm(DE4.DE)

# DD only
DD100 <- assay(blood_100_filtered)
DD100.DD <- diff.disp(DD100[DDindex_DD100==1, group_100=="2"], integer=T)
DD100[DDindex_DD100==1, group_100=="2"] <- DD100.DD$counts
FC_DD100 <- rep(1, nrow(DD100))
FC_DD100[DDindex_DD100==1] <- DD100.DD$FC
rm(DD100.DD)
DD20 <- assay(blood_20_filtered)
DD20.DD <- diff.disp(DD20[DDindex_DD20==1, group_20=="2"], integer=T)
DD20[DDindex_DD20==1, group_20=="2"] <- DD20.DD$counts
FC_DD20 <- rep(1, nrow(DD20))
FC_DD20[DDindex_DD20==1] <- DD20.DD$FC
rm(DD20.DD)
DD10 <- assay(blood_10_filtered)
DD10.DD <- diff.disp(DD10[DDindex_DD10==1, group_10=="2"], integer=T)
DD10[DDindex_DD10==1, group_10=="2"] <- DD10.DD$counts
FC_DD10 <- rep(1, nrow(DD10))
FC_DD10[DDindex_DD10==1] <- DD10.DD$FC
rm(DD10.DD)
DD4 <- assay(blood_4_filtered)
DD4.DD <- diff.disp(DD4[DDindex_DD4==1, group_4=="2"], integer=T)
DD4[DDindex_DD4==1, group_4=="2"] <- DD4.DD$counts
FC_DD4 <- rep(1, nrow(DD4))
FC_DD4[DDindex_DD4==1] <- DD4.DD$FC
rm(DD4.DD)

# DEDD
DEDD100 <- assay(blood_100_filtered)
DEDD100.DE <- diff.mean(DEDD100[DEindex_DEDD100==1, group_100=="2"], integer=T)
DEDD100[DEindex_DEDD100==1, group_100=="2"] <- DEDD100.DE$counts
FC.mean_DEDD100 <- rep(1, nrow(DEDD100))
FC.mean_DEDD100[DEindex_DEDD100==1] <- DEDD100.DE$FC
DEDD100.DD <- diff.disp(DEDD100[DDindex_DEDD100==1, group_100=="2"], integer=T)
DEDD100[DDindex_DEDD100==1, group_100=="2"] <- DEDD100.DD$counts
FC.disp_DEDD100 <- rep(1, nrow(DEDD100))
FC.disp_DEDD100[DDindex_DEDD100==1] <- DEDD100.DD$FC
rm(DEDD100.DE, DEDD100.DD)
DEDD20 <- assay(blood_20_filtered)
DEDD20.DE <- diff.mean(DEDD20[DEindex_DEDD20==1, group_20=="2"], integer=T)
DEDD20[DEindex_DEDD20==1, group_20=="2"] <- DEDD20.DE$counts
FC.mean_DEDD20 <- rep(1, nrow(DEDD20))
FC.mean_DEDD20[DEindex_DEDD20==1] <- DEDD20.DE$FC
DEDD20.DD <- diff.disp(DEDD20[DDindex_DEDD20==1, group_20=="2"], integer=T)
DEDD20[DDindex_DEDD20==1, group_20=="2"] <- DEDD20.DD$counts
FC.disp_DEDD20 <- rep(1, nrow(DEDD20))
FC.disp_DEDD20[DDindex_DEDD20==1] <- DEDD20.DD$FC
rm(DEDD20.DE, DEDD20.DD)
DEDD10 <- assay(blood_10_filtered)
DEDD10.DE <- diff.mean(DEDD10[DEindex_DEDD10==1, group_10=="2"], integer=T)
DEDD10[DEindex_DEDD10==1, group_10=="2"] <- DEDD10.DE$counts
FC.mean_DEDD10 <- rep(1, nrow(DEDD10))
FC.mean_DEDD10[DEindex_DEDD10==1] <- DEDD10.DE$FC
DEDD10.DD <- diff.disp(DEDD10[DDindex_DEDD10==1, group_10=="2"], integer=T)
DEDD10[DDindex_DEDD10==1, group_10=="2"] <- DEDD10.DD$counts
FC.disp_DEDD10 <- rep(1, nrow(DEDD10))
FC.disp_DEDD10[DDindex_DEDD10==1] <- DEDD10.DD$FC
rm(DEDD10.DE, DEDD10.DD)
DEDD4 <- assay(blood_4_filtered)
DEDD4.DE <- diff.mean(DEDD4[DEindex_DEDD4==1, group_4=="2"], integer=T)
DEDD4[DEindex_DEDD4==1, group_4=="2"] <- DEDD4.DE$counts
FC.mean_DEDD4 <- rep(1, nrow(DEDD4))
FC.mean_DEDD4[DEindex_DEDD4==1] <- DEDD4.DE$FC
DEDD4.DD <- diff.disp(DEDD4[DDindex_DEDD4==1, group_4=="2"], integer=T)
DEDD4[DDindex_DEDD4==1, group_4=="2"] <- DEDD4.DD$counts
FC.disp_DEDD4 <- rep(1, nrow(DEDD4))
FC.disp_DEDD4[DDindex_DEDD4==1] <- DEDD4.DD$FC
rm(DEDD4.DE, DEDD4.DD)

# Apply 0.5 cpm filter again after changing counts
for (i in c(100, 20, 10, 4)) {
  for (j in c("DE", "DD")) {
  filter <- which(rowMeans(apply(get(paste0(j,i)), 2, function(x) x*1e6/sum(x))) > 0.5)
  assign(paste0(j,i), get(paste0(j,i))[filter, ])
  assign(paste0(j, "index_", j, i), get(paste0(j, "index_", j, i))[filter])
  assign(paste0("FC_", j, i), get(paste0("FC_", j, i))[filter])
  }
}
for (i in c("DEDD100", "DEDD20", "DEDD10", "DEDD4")) {
  filter <- which(rowMeans(apply(get(i), 2, function(x) x*1e6/sum(x))) > 0.5)
  assign(i, get(i)[filter, ])
  assign(paste0("DEindex_", i), get(paste0("DEindex_", i))[filter])
  assign(paste0("DDindex_", i), get(paste0("DDindex_", i))[filter])
  assign(paste0("DEDDindex_", i), get(paste0("DEDDindex_", i))[filter])
  assign(paste0("DEonlyindex_", i), get(paste0("DEonlyindex_", i))[filter])
  assign(paste0("DDonlyindex_", i), get(paste0("DDonlyindex_", i))[filter])
  assign(paste0("DEplusDDindex_", i), get(paste0("DEplusDDindex_", i))[filter])
  assign(paste0("FC.mean_", i), get(paste0("FC.mean_", i))[filter])
  assign(paste0("FC.disp_", i), get(paste0("FC.disp_", i))[filter])
}

# Compile and save data
DE100 <- list(counts = DE100, 
              DEindex = DEindex_DE100, 
              FC = FC_DE100)
DE20 <- list(counts = DE20, 
             DEindex = DEindex_DE20, 
             FC = FC_DE20)
DE10 <- list(counts = DE10, 
             DEindex = DEindex_DE10, 
             FC = FC_DE10)
DE4 <- list(counts = DE4, 
            DEindex = DEindex_DE4, 
            FC = FC_DE4)
DD100 <- list(counts = DD100, 
              DDindex = DDindex_DD100, 
              FC = FC_DD100)
DD20 <- list(counts = DD20, 
             DDindex = DDindex_DD20, 
             FC = FC_DD20)
DD10 <- list(counts = DD10, 
             DDindex = DDindex_DD10, 
             FC = FC_DD10)
DD4 <- list(counts = DD4, 
            DDindex = DDindex_DD4, 
            FC = FC_DD4)
DEDD100 <- list(counts = DEDD100, 
                DEindex = DEindex_DEDD100, 
                DDindex = DDindex_DEDD100, 
                DEDDindex = DEDDindex_DEDD100, 
                DEonlyindex = DEonlyindex_DEDD100, 
                DDonlyindex = DDonlyindex_DEDD100, 
                DEplusDDindex = DEplusDDindex_DEDD100, 
                FC.mean = FC.mean_DEDD100, 
                FC.disp = FC.disp_DEDD100)
DEDD20 <- list(counts = DEDD20, 
                DEindex = DEindex_DEDD20, 
                DDindex = DDindex_DEDD20, 
                DEDDindex = DEDDindex_DEDD20, 
                DEonlyindex = DEonlyindex_DEDD20, 
                DDonlyindex = DDonlyindex_DEDD20, 
                DEplusDDindex = DEplusDDindex_DEDD20, 
                FC.mean = FC.mean_DEDD20, 
                FC.disp = FC.disp_DEDD20)
DEDD10 <- list(counts = DEDD10, 
                DEindex = DEindex_DEDD10, 
                DDindex = DDindex_DEDD10, 
                DEDDindex = DEDDindex_DEDD10, 
                DEonlyindex = DEonlyindex_DEDD10, 
                DDonlyindex = DDonlyindex_DEDD10, 
                DEplusDDindex = DEplusDDindex_DEDD10, 
                FC.mean = FC.mean_DEDD10, 
                FC.disp = FC.disp_DEDD10)
DEDD4 <- list(counts = DEDD4, 
                DEindex = DEindex_DEDD4, 
                DDindex = DDindex_DEDD4, 
                DEDDindex = DEDDindex_DEDD4, 
                DEonlyindex = DEonlyindex_DEDD4, 
                DDonlyindex = DDonlyindex_DEDD4, 
                DEplusDDindex = DEplusDDindex_DEDD4, 
                FC.mean = FC.mean_DEDD4, 
                FC.disp = FC.disp_DEDD4)
saveRDS(DE100, file=here("recount data/GTEx", "blood_100_DE"))
saveRDS(DE20, file=here("recount data/GTEx", "blood_20_DE"))
saveRDS(DE10, file=here("recount data/GTEx", "blood_10_DE"))
saveRDS(DE4, file=here("recount data/GTEx", "blood_4_DE"))
saveRDS(DD100, file=here("recount data/GTEx", "blood_100_DD"))
saveRDS(DD20, file=here("recount data/GTEx", "blood_20_DD"))
saveRDS(DD10, file=here("recount data/GTEx", "blood_10_DD"))
saveRDS(DD4, file=here("recount data/GTEx", "blood_4_DD"))
saveRDS(DEDD100, file=here("recount data/GTEx", "blood_100_DEDD"))
saveRDS(DEDD20, file=here("recount data/GTEx", "blood_20_DEDD"))
saveRDS(DEDD10, file=here("recount data/GTEx", "blood_10_DEDD"))
saveRDS(DEDD4, file=here("recount data/GTEx", "blood_4_DEDD"))



