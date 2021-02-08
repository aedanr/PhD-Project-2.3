library(here)
library(recount)

## Explore GTEx data through recount ####
load(here("recount data/GTEx", "rse_gene.Rdata"))
dim(rse_gene)
# 58037 9662
dim(colData(rse_gene))
# 9662 82
info <- colData(rse_gene)
names(info)
# Data dictionaries:
attr <- read.csv(here("recount data/GTEx/Data dictionaries from GTEx", 
                      "GTEx_Analysis_v8_Annotations_SampleAttributesDD.csv"), 
                 stringsAsFactors=F)
# smrin: RINs
# smts: tissue source
# smtsd: more detailed tissue source
# smnabtcht: RNA extraction method
# smubrid: Uberon ID, anatomical location as described by the Uber Anatomy Ontology (UBERON)
# smafrze: Samples included in the GTEx Analysis Freeze
# Not much other useful info. Mostly QC. smts and smtsd will be main info of interest. 
# Maybe also smnabtch, smgebtch, smcenter - batch and centre IDs.
phen <- read.csv(here("recount data/GTEx/Data dictionaries from GTEx", 
                      "GTEx_Analysis_v8_Annotations_SubjectPhenotypesDD.csv"), 
                 stringsAsFactors=F)
# Not in colData(). Age, sex and death circumstances.

muscle <- rse_gene[, which(colData(rse_gene)$smtsd == "Muscle - Skeletal")]
dim(muscle) # 58037 475
muscle_info <- colData(muscle)
rm(rse_gene, info)

length(table(muscle_info$sampid))
# 475
length(table(muscle_info$smnabtch))
# 247
length(table(muscle_info$smgebtch))
# 120
# Not feasible to use a single NA extraction or GE batch
length(table(muscle_info$smnabtcht))
# 3
table(muscle_info$smnabtcht)
# 275 RNA Extraction from Paxgene-derived Lysate Plate Based
# 162 RNA isolation_PAXgene Tissue miRNA
# 38 RNA isolation_Trizol Manual (Cell Pellet)
# To replicate blood analysis, need 400 samples, so exclude only Trizol Manual

muscle <- muscle[, which(colData(muscle)$smnabtcht == 
                           "RNA Extraction from Paxgene-derived Lysate Plate Based" | 
                           colData(muscle)$smnabtcht == 
                           "RNA isolation_PAXgene Tissue miRNA")]
dim(muscle)
# 58037 437
muscle_info <- colData(muscle)
sum(muscle_info$smrin >= 8)
# 240
sum(muscle_info$smrin >= 7)
# 395
sum(muscle_info$smrin >= 6.9)
# 401
muscle <- muscle[, which(muscle_info$smrin >= 6.9)]
dim(muscle)
# 58037 401
muscle_info <- colData(muscle)

# Convert from read counts to gene counts
muscle <- scale_counts(muscle)

# Sample without replacement to make 10 sets of 2 groups of each of 2, 5, 10 and 20, 
# and with replacement to make 10 sets of 2 groups of 50 (not actually sampling with 
# replacement - sampling without replacement within each set of 2 groups of 50, but 
# allowing overlap between sets - so 10 independent sets of 2 groups of 50 sampled 
# without replacement)
all_indices_2 <- sample(1:ncol(muscle), 40)
set1_2 <- all_indices_2[1:4]
set2_2 <- all_indices_2[5:8]
set3_2 <- all_indices_2[9:12]
set4_2 <- all_indices_2[13:16]
set5_2 <- all_indices_2[17:20]
set6_2 <- all_indices_2[21:24]
set7_2 <- all_indices_2[25:28]
set8_2 <- all_indices_2[29:32]
set9_2 <- all_indices_2[33:36]
set10_2 <- all_indices_2[37:40]
all_indices_5 <- sample(1:ncol(muscle), 100)
set1_5 <- all_indices_5[1:10]
set2_5 <- all_indices_5[11:20]
set3_5 <- all_indices_5[21:30]
set4_5 <- all_indices_5[31:40]
set5_5 <- all_indices_5[41:50]
set6_5 <- all_indices_5[51:60]
set7_5 <- all_indices_5[61:70]
set8_5 <- all_indices_5[71:80]
set9_5 <- all_indices_5[81:90]
set10_5 <- all_indices_5[91:100]
all_indices_10 <- sample(1:ncol(muscle), 200)
set1_10 <- all_indices_10[1:20]
set2_10 <- all_indices_10[21:40]
set3_10 <- all_indices_10[41:60]
set4_10 <- all_indices_10[61:80]
set5_10 <- all_indices_10[81:100]
set6_10 <- all_indices_10[101:120]
set7_10 <- all_indices_10[121:140]
set8_10 <- all_indices_10[141:160]
set9_10 <- all_indices_10[161:180]
set10_10 <- all_indices_10[181:200]
all_indices_20 <- sample(1:ncol(muscle), 400)
set1_20 <- all_indices_20[1:40]
set2_20 <- all_indices_20[41:80]
set3_20 <- all_indices_20[81:120]
set4_20 <- all_indices_20[121:160]
set5_20 <- all_indices_20[161:200]
set6_20 <- all_indices_20[201:240]
set7_20 <- all_indices_20[241:280]
set8_20 <- all_indices_20[281:320]
set9_20 <- all_indices_20[321:360]
set10_20 <- all_indices_20[361:400]
set1_50 <- sample(1:ncol(muscle), 100)
set2_50 <- sample(1:ncol(muscle), 100)
set3_50 <- sample(1:ncol(muscle), 100)
set4_50 <- sample(1:ncol(muscle), 100)
set5_50 <- sample(1:ncol(muscle), 100)
set6_50 <- sample(1:ncol(muscle), 100)
set7_50 <- sample(1:ncol(muscle), 100)
set8_50 <- sample(1:ncol(muscle), 100)
set9_50 <- sample(1:ncol(muscle), 100)
set10_50 <- sample(1:ncol(muscle), 100)

muscle_2_set1 <- muscle[, set1_2]
saveRDS(muscle_2_set1, file=here("recount data/GTEx/muscle_2_samples_per_group", "muscle_2_set1_unfiltered.rds"))
rm(muscle_2_set1)
muscle_2_set2 <- muscle[, set2_2]
saveRDS(muscle_2_set2, file=here("recount data/GTEx/muscle_2_samples_per_group", "muscle_2_set2_unfiltered.rds"))
rm(muscle_2_set2)
muscle_2_set3 <- muscle[, set3_2]
saveRDS(muscle_2_set3, file=here("recount data/GTEx/muscle_2_samples_per_group", "muscle_2_set3_unfiltered.rds"))
rm(muscle_2_set3)
muscle_2_set4 <- muscle[, set4_2]
saveRDS(muscle_2_set4, file=here("recount data/GTEx/muscle_2_samples_per_group", "muscle_2_set4_unfiltered.rds"))
rm(muscle_2_set4)
muscle_2_set5 <- muscle[, set5_2]
saveRDS(muscle_2_set5, file=here("recount data/GTEx/muscle_2_samples_per_group", "muscle_2_set5_unfiltered.rds"))
rm(muscle_2_set5)
muscle_2_set6 <- muscle[, set6_2]
saveRDS(muscle_2_set6, file=here("recount data/GTEx/muscle_2_samples_per_group", "muscle_2_set6_unfiltered.rds"))
rm(muscle_2_set6)
muscle_2_set7 <- muscle[, set7_2]
saveRDS(muscle_2_set7, file=here("recount data/GTEx/muscle_2_samples_per_group", "muscle_2_set7_unfiltered.rds"))
rm(muscle_2_set7)
muscle_2_set8 <- muscle[, set8_2]
saveRDS(muscle_2_set8, file=here("recount data/GTEx/muscle_2_samples_per_group", "muscle_2_set8_unfiltered.rds"))
rm(muscle_2_set8)
muscle_2_set9 <- muscle[, set9_2]
saveRDS(muscle_2_set9, file=here("recount data/GTEx/muscle_2_samples_per_group", "muscle_2_set9_unfiltered.rds"))
rm(muscle_2_set9)
muscle_2_set10 <- muscle[, set10_2]
saveRDS(muscle_2_set10, file=here("recount data/GTEx/muscle_2_samples_per_group", "muscle_2_set10_unfiltered.rds"))
rm(muscle_2_set10)

muscle_5_set1 <- muscle[, set1_5]
saveRDS(muscle_5_set1, file=here("recount data/GTEx/muscle_5_samples_per_group", "muscle_5_set1_unfiltered.rds"))
rm(muscle_5_set1)
muscle_5_set2 <- muscle[, set2_5]
saveRDS(muscle_5_set2, file=here("recount data/GTEx/muscle_5_samples_per_group", "muscle_5_set2_unfiltered.rds"))
rm(muscle_5_set2)
muscle_5_set3 <- muscle[, set3_5]
saveRDS(muscle_5_set3, file=here("recount data/GTEx/muscle_5_samples_per_group", "muscle_5_set3_unfiltered.rds"))
rm(muscle_5_set3)
muscle_5_set4 <- muscle[, set4_5]
saveRDS(muscle_5_set4, file=here("recount data/GTEx/muscle_5_samples_per_group", "muscle_5_set4_unfiltered.rds"))
rm(muscle_5_set4)
muscle_5_set5 <- muscle[, set5_5]
saveRDS(muscle_5_set5, file=here("recount data/GTEx/muscle_5_samples_per_group", "muscle_5_set5_unfiltered.rds"))
rm(muscle_5_set5)
muscle_5_set6 <- muscle[, set6_5]
saveRDS(muscle_5_set6, file=here("recount data/GTEx/muscle_5_samples_per_group", "muscle_5_set6_unfiltered.rds"))
rm(muscle_5_set6)
muscle_5_set7 <- muscle[, set7_5]
saveRDS(muscle_5_set7, file=here("recount data/GTEx/muscle_5_samples_per_group", "muscle_5_set7_unfiltered.rds"))
rm(muscle_5_set7)
muscle_5_set8 <- muscle[, set8_5]
saveRDS(muscle_5_set8, file=here("recount data/GTEx/muscle_5_samples_per_group", "muscle_5_set8_unfiltered.rds"))
rm(muscle_5_set8)
muscle_5_set9 <- muscle[, set9_5]
saveRDS(muscle_5_set9, file=here("recount data/GTEx/muscle_5_samples_per_group", "muscle_5_set9_unfiltered.rds"))
rm(muscle_5_set9)
muscle_5_set10 <- muscle[, set10_5]
saveRDS(muscle_5_set10, file=here("recount data/GTEx/muscle_5_samples_per_group", "muscle_5_set10_unfiltered.rds"))
rm(muscle_5_set10)

muscle_10_set1 <- muscle[, set1_10]
saveRDS(muscle_10_set1, file=here("recount data/GTEx/muscle_10_samples_per_group", "muscle_10_set1_unfiltered.rds"))
rm(muscle_10_set1)
muscle_10_set2 <- muscle[, set2_10]
saveRDS(muscle_10_set2, file=here("recount data/GTEx/muscle_10_samples_per_group", "muscle_10_set2_unfiltered.rds"))
rm(muscle_10_set2)
muscle_10_set3 <- muscle[, set3_10]
saveRDS(muscle_10_set3, file=here("recount data/GTEx/muscle_10_samples_per_group", "muscle_10_set3_unfiltered.rds"))
rm(muscle_10_set3)
muscle_10_set4 <- muscle[, set4_10]
saveRDS(muscle_10_set4, file=here("recount data/GTEx/muscle_10_samples_per_group", "muscle_10_set4_unfiltered.rds"))
rm(muscle_10_set4)
muscle_10_set5 <- muscle[, set5_10]
saveRDS(muscle_10_set5, file=here("recount data/GTEx/muscle_10_samples_per_group", "muscle_10_set5_unfiltered.rds"))
rm(muscle_10_set5)
muscle_10_set6 <- muscle[, set6_10]
saveRDS(muscle_10_set6, file=here("recount data/GTEx/muscle_10_samples_per_group", "muscle_10_set6_unfiltered.rds"))
rm(muscle_10_set6)
muscle_10_set7 <- muscle[, set7_10]
saveRDS(muscle_10_set7, file=here("recount data/GTEx/muscle_10_samples_per_group", "muscle_10_set7_unfiltered.rds"))
rm(muscle_10_set7)
muscle_10_set8 <- muscle[, set8_10]
saveRDS(muscle_10_set8, file=here("recount data/GTEx/muscle_10_samples_per_group", "muscle_10_set8_unfiltered.rds"))
rm(muscle_10_set8)
muscle_10_set9 <- muscle[, set9_10]
saveRDS(muscle_10_set9, file=here("recount data/GTEx/muscle_10_samples_per_group", "muscle_10_set9_unfiltered.rds"))
rm(muscle_10_set9)
muscle_10_set10 <- muscle[, set10_10]
saveRDS(muscle_10_set10, file=here("recount data/GTEx/muscle_10_samples_per_group", "muscle_10_set10_unfiltered.rds"))
rm(muscle_10_set10)

muscle_20_set1 <- muscle[, set1_20]
saveRDS(muscle_20_set1, file=here("recount data/GTEx/muscle_20_samples_per_group", "muscle_20_set1_unfiltered.rds"))
rm(muscle_20_set1)
muscle_20_set2 <- muscle[, set2_20]
saveRDS(muscle_20_set2, file=here("recount data/GTEx/muscle_20_samples_per_group", "muscle_20_set2_unfiltered.rds"))
rm(muscle_20_set2)
muscle_20_set3 <- muscle[, set3_20]
saveRDS(muscle_20_set3, file=here("recount data/GTEx/muscle_20_samples_per_group", "muscle_20_set3_unfiltered.rds"))
rm(muscle_20_set3)
muscle_20_set4 <- muscle[, set4_20]
saveRDS(muscle_20_set4, file=here("recount data/GTEx/muscle_20_samples_per_group", "muscle_20_set4_unfiltered.rds"))
rm(muscle_20_set4)
muscle_20_set5 <- muscle[, set5_20]
saveRDS(muscle_20_set5, file=here("recount data/GTEx/muscle_20_samples_per_group", "muscle_20_set5_unfiltered.rds"))
rm(muscle_20_set5)
muscle_20_set6 <- muscle[, set6_20]
saveRDS(muscle_20_set6, file=here("recount data/GTEx/muscle_20_samples_per_group", "muscle_20_set6_unfiltered.rds"))
rm(muscle_20_set6)
muscle_20_set7 <- muscle[, set7_20]
saveRDS(muscle_20_set7, file=here("recount data/GTEx/muscle_20_samples_per_group", "muscle_20_set7_unfiltered.rds"))
rm(muscle_20_set7)
muscle_20_set8 <- muscle[, set8_20]
saveRDS(muscle_20_set8, file=here("recount data/GTEx/muscle_20_samples_per_group", "muscle_20_set8_unfiltered.rds"))
rm(muscle_20_set8)
muscle_20_set9 <- muscle[, set9_20]
saveRDS(muscle_20_set9, file=here("recount data/GTEx/muscle_20_samples_per_group", "muscle_20_set9_unfiltered.rds"))
rm(muscle_20_set9)
muscle_20_set10 <- muscle[, set10_20]
saveRDS(muscle_20_set10, file=here("recount data/GTEx/muscle_20_samples_per_group", "muscle_20_set10_unfiltered.rds"))
rm(muscle_20_set10)

muscle_50_set1 <- muscle[, set1_50]
saveRDS(muscle_50_set1, file=here("recount data/GTEx/muscle_50_samples_per_group", "muscle_50_set1_unfiltered.rds"))
rm(muscle_50_set1)
muscle_50_set2 <- muscle[, set2_50]
saveRDS(muscle_50_set2, file=here("recount data/GTEx/muscle_50_samples_per_group", "muscle_50_set2_unfiltered.rds"))
rm(muscle_50_set2)
muscle_50_set3 <- muscle[, set3_50]
saveRDS(muscle_50_set3, file=here("recount data/GTEx/muscle_50_samples_per_group", "muscle_50_set3_unfiltered.rds"))
rm(muscle_50_set3)
muscle_50_set4 <- muscle[, set4_50]
saveRDS(muscle_50_set4, file=here("recount data/GTEx/muscle_50_samples_per_group", "muscle_50_set4_unfiltered.rds"))
rm(muscle_50_set4)
muscle_50_set5 <- muscle[, set5_50]
saveRDS(muscle_50_set5, file=here("recount data/GTEx/muscle_50_samples_per_group", "muscle_50_set5_unfiltered.rds"))
rm(muscle_50_set5)
muscle_50_set6 <- muscle[, set6_50]
saveRDS(muscle_50_set6, file=here("recount data/GTEx/muscle_50_samples_per_group", "muscle_50_set6_unfiltered.rds"))
rm(muscle_50_set6)
muscle_50_set7 <- muscle[, set7_50]
saveRDS(muscle_50_set7, file=here("recount data/GTEx/muscle_50_samples_per_group", "muscle_50_set7_unfiltered.rds"))
rm(muscle_50_set7)
muscle_50_set8 <- muscle[, set8_50]
saveRDS(muscle_50_set8, file=here("recount data/GTEx/muscle_50_samples_per_group", "muscle_50_set8_unfiltered.rds"))
rm(muscle_50_set8)
muscle_50_set9 <- muscle[, set9_50]
saveRDS(muscle_50_set9, file=here("recount data/GTEx/muscle_50_samples_per_group", "muscle_50_set9_unfiltered.rds"))
rm(muscle_50_set9)
muscle_50_set10 <- muscle[, set10_50]
saveRDS(muscle_50_set10, file=here("recount data/GTEx/muscle_50_samples_per_group", "muscle_50_set10_unfiltered.rds"))
rm(muscle_50_set10)

indices_2 <- data.frame(set1 = set1_2, 
                        set2 = set2_2, 
                        set3 = set3_2, 
                        set4 = set4_2, 
                        set5 = set5_2, 
                        set6 = set6_2, 
                        set7 = set7_2, 
                        set8 = set8_2, 
                        set9 = set9_2, 
                        set10 = set10_2)
indices_5 <- data.frame(set1 = set1_5, 
                        set2 = set2_5, 
                        set3 = set3_5, 
                        set4 = set4_5, 
                        set5 = set5_5, 
                        set6 = set6_5, 
                        set7 = set7_5, 
                        set8 = set8_5, 
                        set9 = set9_5, 
                        set10 = set10_5)
indices_10 <- data.frame(set1 = set1_10, 
                         set2 = set2_10, 
                         set3 = set3_10, 
                         set4 = set4_10, 
                         set5 = set5_10, 
                         set6 = set6_10, 
                         set7 = set7_10, 
                         set8 = set8_10, 
                         set9 = set9_10, 
                         set10 = set10_10)
indices_20 <- data.frame(set1 = set1_20, 
                         set2 = set2_20, 
                         set3 = set3_20, 
                         set4 = set4_20, 
                         set5 = set5_20, 
                         set6 = set6_20, 
                         set7 = set7_20, 
                         set8 = set8_20, 
                         set9 = set9_20, 
                         set10 = set10_20)
indices_50 <- data.frame(set1 = set1_50, 
                         set2 = set2_50, 
                         set3 = set3_50, 
                         set4 = set4_50, 
                         set5 = set5_50, 
                         set6 = set6_50, 
                         set7 = set7_50, 
                         set8 = set8_50, 
                         set9 = set9_50, 
                         set10 = set10_50)
saveRDS(indices_2, file=here("recount data/GTEx/muscle_2_samples_per_group", 
                             "indices_2_from_whole_muscle_PAXGene_RIN_at_least_6point9.rds"))
saveRDS(indices_5, file=here("recount data/GTEx/muscle_5_samples_per_group", 
                             "indices_5_from_whole_muscle_PAXGene_RIN_at_least_6point9.rds"))
saveRDS(indices_10, file=here("recount data/GTEx/muscle_10_samples_per_group", 
                              "indices_10_from_whole_muscle_PAXGene_RIN_at_least_6point9.rds"))
saveRDS(indices_20, file=here("recount data/GTEx/muscle_20_samples_per_group", 
                              "indices_20_from_whole_muscle_PAXGene_RIN_at_least_6point9.rds"))
saveRDS(indices_50, file=here("recount data/GTEx/muscle_50_samples_per_group", 
                              "indices_50_from_whole_muscle_PAXGene_RIN_at_least_6point9.rds"))

# rm(list=ls())


#### Re-import and create DE, DD data ####
source(here('scripts','2019-11-29_functions_to_change_mean_disp_in_real_data.R'))

## Filter low-count genes and make DE, DD data
group_2 <- factor(c(rep("1",2), rep("2",2)))
group_5 <- factor(c(rep("1",5), rep("2",5)))
group_10 <- factor(c(rep("1",10), rep("2",10)))
group_20 <- factor(c(rep("1",20), rep("2",20)))
group_50 <- factor(c(rep("1",50), rep("2",50)))

for (i in 1:10) {
  for (j in c(2, 5, 10, 20, 50)) {
    # Import data
    raw.data <- readRDS(here(paste0("recount data/GTEx/muscle_", j, "_samples_per_group"), 
                             paste0("muscle_", j, "_set", i, "_unfiltered.rds")))
    # Apply 0.5 cpm filter
    raw.data <- raw.data[which(rowMeans(apply(assay(raw.data), 2, function(x) x*1e6/sum(x))) > 0.5), ]
    # Create indices for DE and DD
    rows <- nrow(raw.data)
    break1 <- ceiling(rows/20)
    break2 <- ceiling(rows/10)
    break3 <- ceiling(3*rows/20)
    DEindex_DE <- numeric(rows)
    DDindex_DD <- numeric(rows)
    DDindex_DEDD <- numeric(rows)
    DEindex_DEDD <- numeric(rows)
    DEindex_DE[1:break2] <- 1
    DDindex_DD <- DEindex_DE
    DEindex_DEDD[1:break2] <- 1
    DDindex_DEDD[(break1 + 1):break3] <- 1
    rm(rows, break1, break2, break3)
    # Create DE
    DE <- assay(raw.data)
    DE.DE_counts <- diff.mean(DE[DEindex_DE == 1, get(paste0("group_", j)) == "2"], integer=T)
    DE[DEindex_DE == 1, get(paste0("group_", j)) == "2"] <- DE.DE_counts$counts
    DE_FC <- rep(1, nrow(DE))
    DE_FC[DEindex_DE == 1] <- DE.DE_counts$FC
    rm(DE.DE_counts)
    # Create DD
    DD <- assay(raw.data)
    DD.DD_counts <- diff.disp(DD[DDindex_DD == 1, get(paste0("group_", j)) == "2"], integer=T)
    DD[DDindex_DD == 1, get(paste0("group_", j)) == "2"] <- DD.DD_counts$counts
    DD_FC <- rep(1, nrow(DD))
    DD_FC[DDindex_DD == 1] <- DD.DD_counts$FC
    rm(DD.DD_counts)
    # Create DEDD
    DEDD <- assay(raw.data)
    DEDD.DE_counts <- diff.mean(DEDD[DEindex_DEDD == 1, get(paste0("group_", j)) == "2"], integer=T)
    DEDD[DEindex_DEDD == 1, get(paste0("group_", j)) == "2"] <- DEDD.DE_counts$counts
    DEDD_FC.mean <- rep(1, nrow(DEDD))
    DEDD_FC.mean[DEindex_DEDD == 1] <- DEDD.DE_counts$FC
    rm(DEDD.DE_counts)
    DEDD.DD_counts <- diff.disp(DEDD[DDindex_DEDD == 1, get(paste0("group_", j)) == "2"], integer=T)
    DEDD[DDindex_DEDD == 1, get(paste0("group_", j)) == "2"] <- DEDD.DD_counts$counts
    DEDD_FC.disp <- rep(1, nrow(DEDD))
    DEDD_FC.disp[DDindex_DEDD == 1] <- DEDD.DD_counts$FC
    rm(DEDD.DD_counts)
    # Re-apply 0.5 cpm filter
    filter_DE <- which(rowMeans(apply(DE, 2, function(x) x*1e6/sum(x))) > 0.5)
    DE <- DE[filter_DE, ]
    DEindex_DE <- DEindex_DE[filter_DE]
    DE_FC <- DE_FC[filter_DE]
    rm(filter_DE)
    filter_DD <- which(rowMeans(apply(DD, 2, function(x) x*1e6/sum(x))) > 0.5)
    DD <- DD[filter_DD, ]
    DDindex_DD <- DDindex_DD[filter_DD]
    DD_FC <- DD_FC[filter_DD]
    rm(filter_DD)
    filter_DEDD <- which(rowMeans(apply(DEDD, 2, function(x) x*1e6/sum(x))) > 0.5)
    DEDD <- DEDD[filter_DEDD, ]
    DEindex_DEDD <- DEindex_DEDD[filter_DEDD]
    DDindex_DEDD <- DDindex_DEDD[filter_DEDD]
    DEDD_FC.mean <- DEDD_FC.mean[filter_DEDD]
    DEDD_FC.disp <- DEDD_FC.disp[filter_DEDD]
    rm(filter_DEDD)
    # Compile and save
    DE <- list(counts = DE, 
               DEindex = DEindex_DE, 
               FC = DE_FC)
    rm(DEindex_DE, DE_FC)
    DD <- list(counts = DD, 
               DDindex = DDindex_DD, 
               FC = DD_FC)
    rm(DDindex_DD, DD_FC)
    DEDD <- list(counts = DEDD, 
                 DEindex = DEindex_DEDD, 
                 DDindex = DDindex_DEDD, 
                 FC.mean = DEDD_FC.mean, 
                 FC.disp = DEDD_FC.disp)
    rm(DEindex_DEDD, DDindex_DEDD, DEDD_FC.mean, DEDD_FC.disp)
    saveRDS(raw.data, file = here(paste0("recount data/GTEx/muscle_", j, "_samples_per_group"), 
                                  paste0("muscle_", j, "_set", i, "_0point5_cpm_filter.rds")))
    saveRDS(DE, file = here(paste0("recount data/GTEx/muscle_", j, "_samples_per_group"), 
                            paste0("muscle_", j, "_set", i, "_DE.rds")))
    saveRDS(DD, file = here(paste0("recount data/GTEx/muscle_", j, "_samples_per_group"), 
                            paste0("muscle_", j, "_set", i, "_DD.rds")))
    saveRDS(DEDD, file = here(paste0("recount data/GTEx/muscle_", j, "_samples_per_group"), 
                              paste0("muscle_", j, "_set", i, "_DEDD.rds")))
    rm(DE, DD, DEDD)
  }
}





