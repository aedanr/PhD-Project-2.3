## Some code for processing GTEx data taken from MDSeq supplementary data code

# fread from data.table package allows .gct files to be read
library(data.table)
data <- fread("GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct")
dim(data) # 56200 17384
head(names(data))
length(data$Name)
head(data$Name)
head(data$Description)
head(rownames(data))
class(data)
data <- as.data.frame(data)
dim(data)
# ff1 function from MDSeq to combine ensemble gene IDs and symbols
# ff1 <- function(x) paste(x, collapse = '_')
# rownames(data) <- apply(data[1:2], 1, ff1)
rownames(data) <- data$Name
data1 <- data[,c(-1,-2)] # remove first two columns from GTEx, which are ensembl gene ID and symbol
dim(data) # 56200 17384
dim(data1) # 56200 17382
head(rownames(data1))
# plot(rowMeans(data1))
# plot(rowMeans(log(data1)))

geno <- read.delim('GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt', header = T)
dim(geno) # 22951 63
# remove sampid without 'GTEX' - taken from MDSeq
# don't know why there are sample IDs without GTEX or what they are
geno <- geno[grepl('GTEX', geno$SAMPID), ]
dim(geno) # 22734 63
# ff function from MDSeq to extract subject id
# Subject ID is GTEX-xxxxX - four numbers and a letter - which is the start of 
# the sample ID. Rest of sample IF is varying length of numbers and letters separated 
# by two hyphens
ff <- function(x) paste(unlist(strsplit(as.character(x), '-'))[1:2], collapse = '-')
geno$SUBJID <- sapply(geno$SAMPID,  ff)
# tissue
# as.data.frame(summary(as.factor(geno$SMTSD)))
head(geno$SAMPID)

load("rse_gene.Rdata")
dim(rse_gene)
# 58037 9662
library(recount)
dim(colData(rse_gene))
# 9662 82
info <- colData(rse_gene)
names(info)
head(info$sampid)
head(geno$SAMPID)
length(info$sampid) # 9662
length(geno$SAMPID) # 22951
sum(info$sampid %in% geno$SAMPID) # 9662
ids <- info$sampid[which(info$sampid %in% geno$SAMPID)]
gtex <- data1[, names(data1) %in% ids]
rcnt <- rse_gene[, info$sampid %in% ids]
info <- colData(rcnt)
dim(gtex) # 56200 8671
dim(rcnt) # 58037 9662
rcnt <- rcnt[, info$sampid %in% names(gtex)]
dim(rcnt) # 58037 8671
info <- colData(rcnt)
head(names(rcnt)) # ENSG....
length(names(rcnt)) # 58037
head(rownames(gtex)) # ENSG...
length(rownames(gtex)) # 56200
gtex <- gtex[rownames(gtex) %in% names(rcnt), ]
dim(gtex) # 54892 8671
rcnt <- rcnt[names(rcnt) %in% rownames(gtex), ]
dim(rcnt) # 54892 8671
head(names(gtex)) # sampids
length(names(gtex)) # 8671
head(rownames(rcnt)) # ENSG...
head(colnames(rcnt)) # SRR ids
head(colnames(gtex)) # sampids
class(rcnt) # RangedSummarizedExperiment
class(gtex) # data.frame
colnames(rcnt) <- info$sampid
genes <- names(rcnt)
samps <- colnames(rcnt)
gtex <- gtex[genes, samps]
dim(gtex) # 54892 8671
identical(rownames(gtex), rownames(rcnt)) # TRUE
identical(colnames(gtex), colnames(rcnt)) # TRUE
rcnt <- scale_counts(rcnt)
head(rowMeans(gtex))
head(rowMeans(assay(rcnt)))
plot(log(rowMeans(gtex)), log(rowMeans(assay(rcnt))), pch=20)

geno$SUBJID <- sapply(geno$SAMPID,  ff)
info$SUBJID <- sapply(info$sampid,  ff)
length(unique(info$SUBJID))


## Why does recount data have 9662 samples whereas GTEx has 22951, 
## with 8672 in common? Recount data might be from an older version 
## of GTEx, which might partly but not completely explain it.
## Recount uses GTEx v6 data according to f1000 article. 



