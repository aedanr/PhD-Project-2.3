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
# Mostly undecipherable strings beginning with sm
# smrin looks like RINs, so "sm" must just be a standard prefix
# smts must be tissue source
# smtsd looks like more detailed tissue source
# smubrid must be some sort of ID; 54 values
# smnabtcht is RNA extraction method
# smafrze has 8551 "USE ME" and the rest blank; frozen tissue?
# smrdlgth is 76, 101 or 250 (most 76); must be read length
# Most are numbers, presumably a lot of QC measures, like sm350nrm
# Downloaded data dictionaries directly from GTEx:
attr <- read.csv(here("recount data/GTEx/Data dictionaries from GTEx", 
                        "GTEx_Analysis_v8_Annotations_SampleAttributesDD.csv"), 
                 stringsAsFactors=F)
# smubrid: Uberon ID, anatomical location as described by the Uber Anatomy Ontology (UBERON)
# smafrze: Samples included in the GTEx Analysis Freeze
# Not much other useful info. Mostly QC. smts and smtsd will be main info of interest. 
# Maybe also smnabtch, smgebtch, smcenter - batch and centre IDs.
phen <- read.csv(here("recount data/GTEx/Data dictionaries from GTEx", 
                      "GTEx_Analysis_v8_Annotations_SubjectPhenotypesDD.csv"), 
                 stringsAsFactors=F)
# These don't seem to be in colData(). Age, sex and death circumstances.
# Could select an age range for consistency if want to choose a smaller number of samples and 
# don't have another good way to choose. Otherwise no need to worry and should just select by 
# tissue type.

length(table(info$smtsd))
# 54 more detailed tissue sources
table(info$smtsd)
# Over 100 for quite a lot. Lots to choose from.
blood <- rse_gene[, which(colData(rse_gene)$smtsd == "Whole Blood")]
dim(blood) # 58037 456
blood_info <- colData(blood)
liver <- rse_gene[, which(colData(rse_gene)$smtsd == "Liver")]
dim(liver) # 58037 136
liver_info <- colData(liver)
cml <- rse_gene[, which(colData(rse_gene)$smtsd == "Cells - Leukemia cell line (CML)")]
dim(cml) # 58037 102
cml_info <- colData(cml)

c(length(table(blood_info$sampid)), 
  length(table(liver_info$sampid)), 
  length(table(cml_info$sampid)))
# 456 136 102; all unique
c(length(table(blood_info$smnabtch)), 
  length(table(liver_info$smnabtch)), 
  length(table(cml_info$smnabtch)))
# 196 100 1
c(length(table(blood_info$smgebtch)), 
  length(table(liver_info$smgebtch)), 
  length(table(cml_info$smgebtch)))
# 84 66 101
# Not feasible to use a single NA extraction or GE batch for blood or liver
c(length(table(blood_info$smnabtcht)), 
  length(table(liver_info$smnabtcht)), 
  length(table(cml_info$smnabtcht)))
# 2 2 1
# Blood: 432 PAXgene Blood RNA, 24 Trizol Manual (Cell Pellet)
# Liver: 88 PAXgene-derived Lysate Plate Based, 48 PaxGene Tissue miRNA
c(length(table(blood_info$smgebtcht)), 
  length(table(liver_info$smgebtcht)), 
  length(table(cml_info$smgebtcht)))
# Doesn't exist - must all be same for RNA-seq

# RNA extraction method shouldn't be something to worry about, but to be safe and 
# because there are way more than enough samples to be able to do it, will use only 
# PAXgene Blood RNA for blood. Could also do only PAXgene-derived Lysate Plate Based 
# for liver (88 samples), but for now will stick with blood.

blood <- blood[, which(colData(blood)$smnabtcht == 
                         "RNA isolation_PAXgene Blood RNA (Manual)")]
dim(blood)
# 58037 432
blood_info <- colData(blood)
# Take samples with highest RIN - 113 with RIN > 8.7
blood <- blood[, which(colData(blood)$smrin > 8.7)]
dim(blood)
# 58037 113
blood_info <- colData(blood)

# Want big sample for differential dispersion/distribution, but small for DE.
blood_100 <- blood[, sample(1:ncol(blood), 100)]
blood_20 <- blood[, sample(1:ncol(blood), 20)]
blood_10 <- blood[, sample(1:ncol(blood), 10)]
blood_4 <- blood[, sample(1:ncol(blood), 4)]
saveRDS(blood_100, file=here("recount data/GTEx", "blood_100.rds"))
saveRDS(blood_20, file=here("recount data/GTEx", "blood_20.rds"))
saveRDS(blood_10, file=here("recount data/GTEx", "blood_10.rds"))
saveRDS(blood_4, file=here("recount data/GTEx", "blood_4.rds"))




