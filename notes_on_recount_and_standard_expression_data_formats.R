library(recount)
library(here)
load(here("recount data/SRP001540 (Pickrell)", "rse_gene.Rdata"))
pickrell <- rse_gene
load(here("recount data/SRP001563 (Cheung)", "rse_gene.Rdata"))
cheung <- rse_gene
rm(rse_gene)

colData(pickrell)$geo_accession
colData(cheung)$geo_accession
# GEO IDs for each sample; 160 for Pickrell, 45 for Cheung

# SummarizedExperiment/RangedSummarizedExperiment objects
# Accessors: assayNames, assays, assay, rowData, colData, dim, dimnames
# Can also find accessors by looking at variable, i.e. "> cheung"
# class, dim, metadata, assays, rownames, rowData, colnames, colData
# Also gives length of most accessors and head and tail
assayNames(pickrell)
# counts
names(assays(pickrell))
# counts
dim(assays(pickrell)$counts)
# 58037 160
dim(assays(cheung)$counts)
# 58037 45
dim(assay(pickrell))
# 58037 160
names(rowData(pickrell))
# gene_id, bp_length, symbol
names(colData(pickrell))
# 21 entries, look to be all per-sample info - sample, run, total read counts, ...
dim(pickrell)
# 58037 160
dimnames(pickrell)
# sample IDs beginning with SRR
dimnames(cheung)
# list of 2; first looks like gene IDs, beginning ENSG, second sample IDS

# From help for SummarizedExperiment:
# Rows are features; info about features is stored in a data frame accessible using rowData(). 
# DF has same number of rows as the SE object and its columns represent attributes of the 
# feature of interest, e.g. gene/transcript IDs.
# Columns are samples; info about samples is stored in a df accessible using colData(). 
# DS has same number of rows as SE object as columns, and each row provides info on the 
# sample in the corresponding column of the SE object. DF Columns represent sample attributes, 
# e.g. tissue of origin, and can be annotated using mcols(); names typically provide a short 
# identifier unique to each sample.
# SE object can also contain info about overall experiment, e.g. lab, publications, in a list 
# accessed through metadata(), which doesn't have a set format.
# SE object data accessed using assays(), which returns a SimpleList object, each element of 
# which is a matrix with dimensions equal to the SE object. Row and column names of these 
# matrices are either NULL or match those of the SE object.
# 
# Accessors:
# assays(x) - get or set assays; list or SimpleList of matrices with same dimensions as x.
# assay(x,i) - get or set ith assay element.
# assayNames(x) - get or set names of assay() elements.
# rowData(x, use.names=T) - get or set row data; DataFrame.
# colData(x) - get or set column data; DataFrame. Row names must be NULL or match column 
# names of x.
# metadata(x) - get or set experiment data; list with arbitrary content.
# dim(x) - get dimensions (features x samples) of x.
# dimnames(x) - get or set dimension names; usually a list of 2, either NULL or vectors of 
# length matching corresponding dimension.
# 
# Subsetting:
# x[i,j] <- value: create or replace a subset of x; numeric, logical, character or missing; 
# must be a SE object with dimension, dimension names, assay elements consistent with the 
# subset being replaced.
# subset(x, subset, select) - create a subset of x using an expression subset, referring to 
# columns of rowData(x), and/or select, referring to column names of colData(x).
# x$name - access or replace colData column names
# x[[i, ...]] - access or replace column i
# 
# Combining:
# cbind() - combine objects with same features but different samples (columns in assays); 
# colnames in colData(SE) must match and duplicate columns of rowData(SE) must contain the 
# same data. Data in assays are combined by name matching, or by position if all assay 
# names are null; a mixture of names and null gives an error. metadata from all objects are 
# combined into a list with no name checking.
# rbind() - combine objects with same samples bt different features (rows in assays); 
# colnames in rowData(SE) must match and duplicate columns of colData(SE) must contain the 
# same data. Data in assays are combined by name matching, or by position if all assay 
# names are null; a mixture of names and null gives an error. metadata from all objects are 
# combined into a list with no name checking.


# recount functions:
# geo_characteristics() builds a data.frame from the GEO characteristics (pheno) extracted for 
# a given sample. Names of columns correspond to field names. pheno is a DataFrame as created 
# by geo_info(), which uses GEOquery to extract information for a sample.
# 
# abstract_search() finds SRA project IDs that contain the given text in the abstract (as 
# provided by the SRAdb Bioconductor package). e.g. can search by SRA or GEO identifiers. 
# Returns a recount_abstract object - adata.frame with columns for number of samples, species, 
# abstract and project (SRA project ID).
# 
# download_study() downloads SE objects for given recount project (e.g. SRP number); e.g. from 
# project element of a recount_abstract object. Creates a folder in working directory named 
# by project name, containing a .tsv.gz file with data and a file named rse_gene.Rdata which 
# can be used to access data; load data using load(file.path(project name, 'rse_gene.Rdata)).
# 
# browse_study() takes an SRA study ID and invisibly returns the URL; by default opens in 
# browser.
# 
# geo_characteristics() creates a data.frame from the characteristics field of colData (which is 
# a CompressedCharacterList), with a column for each element, rather than a list of 
# 'element name: data'.
# 
# read_counts() computes read counts - counts provided by recount2 are base-pair counts. Not 
# sure exactly what this does but says computes read counts "using AUC".
# scale_counts() does some sort of scaling of raw counts. Seems to take into account read length 
# somehow - presumably related to fact that raw counts are base-pair counts.


project_info <- abstract_search('GSE32465')
download_study(project_info$project)
load(file.path(project_info$project, 'rse_gene.Rdata'))
browse_study(project_info$project)
colData(rse_gene)$geo_accession
geochar <- lapply(split(colData(rse_gene), seq_len(nrow(colData(rse_gene)))), geo_characteristics)
geochar <- do.call(rbind, lapply(geochar, function(x) {
  if('cells' %in% colnames(x)) {
    colnames(x)[colnames(x) == 'cells'] <- 'cell.line'
    return(x)
  } else {
    return(x)
  }
}))
sample_info <- data.frame(
  run = colData(rse_gene)$run,
  group = ifelse(grepl('uninduced', colData(rse_gene)$title), 'uninduced', 'induced'),
  gene_target = sapply(colData(rse_gene)$title, function(x) { strsplit(strsplit(x,
                                                                                'targeting ')[[1]][2], ' gene')[[1]][1] }),
  cell.line = geochar$cell.line
)
rse <- scale_counts(rse_gene)
colData(rse)$group <- sample_info$group
colData(rse)$gene_target <- sample_info$gene_target

dim(assays(rse_gene)$counts)
dim(assays(rse)$counts)
plot(density(log(rowMeans(assays(rse_gene)$counts))))
plot(density(log(rowMeans(assays(rse)$counts))))
range(rowMeans(assays(rse_gene)$counts))
range(rowMeans(assays(rse)$counts))



range(log(rowMeans(assays(pickrell)$counts)))
range(log(rowMeans(assays(read_counts(pickrell))$counts)))
range(log(rowMeans(assays(scale_counts(pickrell))$counts)))
counts <- assays(read_counts(pickrell))$counts
dim(counts)
hist(log(rowMeans(counts)), breaks=10)
counts <- counts[which(rowMeans(counts) >=1), which(colSums(counts) >= 2e6)]
dim(counts)
hist(log(rowMeans(counts)), breaks=10)
names(colData(pickrell))
geo_characteristics(colData(pickrell))
sum(colData(pickrell)$read_count_as_reported_by_sra < 2e6)
colData(pickrell)$avg_read_length



load(here('recount data/Pickrell (original recount data)', 'montpick_eset.RData'))
class(montpick.eset)
montpick.eset
phenoData(montpick.eset)
sampleNames(montpick.eset)
which(pData(montpick.eset) == "YRI")
pickold <- montpick.eset[, which(pData(montpick.eset)$population == "YRI")]
dim(pickold)
hist(log(rowMeans(exprs(pickold))))
pickold <- pickold[which(rowMeans(exprs(pickold)) >=1), which(colSums(exprs(pickold)) >= 2e6)]
dim(pickold)
hist(log(rowMeans(exprs(pickold))))

which(geo_characteristics(colData(pickrell))$cell.line %in% pData(montpick.eset)$sample.id)


