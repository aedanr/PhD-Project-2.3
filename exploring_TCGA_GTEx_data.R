library(here)
library(recount)

## Explore TCGA data through recount ####
load(here("recount data/TCGA", "rse_gene.Rdata"))
dim(rse_gene)
# 58037 11284
dim(colData(rse_gene))
# 11284 864
head(names(colData(rse_gene)), 21) # Standard 21 entries of colData for recount files
# project, sample, experiment, run, read_count_as_reported_by_sra, reads_downloaded, 
# proportion_of_reads_reported_by_sra_downloaded, paired_end, sra_misreported_paired_end, 
# mapped_read_count, auc, sharq_beta_tissue, sharq_beta_cell_type, biosample_submission_date, 
# biosample_publication_date, biosample_update_date, avg_read_length, geo_accession, 
# bigwig_file, title, characteristics
# Looks like "project" is all "TCGA", "paired_end" is all "TRUE", and all other fields except 
# "reads_downloaded", "mapped_read_count", "auc" and "bigwig_file" are NA.

# Other metadata columns:
# Columns 22-158 begin gdc_ and include info on file names and sizes, data types, platform, 
# centres, submitter, projects, demographic info, tissue sources, diagnoses, outcomes, 
# treatment, other clinical info.
# Columns 159-259 begin cgc_ and look to contain similar info.
# Columns 260-864 begin xml_ and look to contain far more detailed clinical info.
info <- colData(rse_gene)
names(info)[grep("category", names(info))]
# gdc_cases.samples.annotations.category - some reasons to possibly exclude samples
grep("sample", names(info), value=T)
# gdc_cases.samples.oct_embedded
# gdc_cases.samples.is_ffpe / cgc_sample_is_ffpe
# gdc_cases.samples.sample_type / cgc_sample_sample_type
# gdc_cases.samples.portions.annotations.notes - data for one sample only which was originally 
# wrongly identified; check if same as one with incorrect barcode from 
# gdc_cases.samples.annotations.category
# gdc_cases.samples.annotations.notes - extensive notes on a small number of samples with 
# possible issues.
grep("project", names(info), value=T)
# gdc_cases.project.name / gdc_cases.tissue_source_site.project - more or less specific cancer 
# types, identical to each other except for capitalisation and spacing
# gdc_cases.project.primary_site - site of primary tumour
# gdc_cases.project.project_id - 3/4 letter TCGA project codes
grep("type", names(info), value=T)
# cgc_file_disease_type - identical to gdc_cases.project.name except for one NA
# xml_histological_type
# xml_histological_type_other
# xml_tumor_type - only for a few samples; Primary, Type 1 or Type 2
# xml_diagnosis_subtype - only for some samples; Non-Papillary or Papillary
grep("cancer", names(info), value=T)
# xml_person_neoplasm_cancer_status
# xml_colorectal_cancer
grep("tumor", names(info), value=T)
# gdc_cases.diagnoses.tumor_stage
# cgc_case_tumor_status
# xml_primary_pathology_tumor_tissue_sites
# xml_primary_pathology_tumor_tissue_site_other
# xml_tumor_tissue_site
# xml_primary_pathology_tumor_tissue_site
# xml_tumor_location - for brain tumours only
names(info)[-unique(c(
  grep("category", names(info)), 
  grep("sample", names(info)), 
  grep("project", names(info)), 
  grep("type", names(info)), 
  grep("cancer", names(info)), 
  grep("tumor", names(info))
))]
# gdc_platform
# gdc_cases.diagnoses.primary_diagnosis / gdc_cases.diagnoses.tissue_or_organ_of_origin / 
# cgc_case_icd_o3site / cgc_case_icd10 / xml_icd_o_3_site
# cgc_file_investigation - 3/4 letter TCGA codes
# cgc_case_primary_site - mostly matches gdc_cases.project.primary_site; check
# cgc_case_histological_diagnosis - mostly matches xml_histological_type; check
# cgc_case_clinical_stage
# cgc_case_pathologic_stage - mostly matches gdc_cases.diagnoses.tumor_stage
# cgc_case_icd_o3histology / xml_icd_o_3_histology
# cgc_case_other_histological_diagnosis
# xml_other_dx
# xml_distant_metastasis_present_ind2
# xml_diagnosis - only for lung cancer, Squamous Cell Carcinoma or Adenocarcinoma
# xml_metastatic_site - mostly empty, must only be for a small subset
# xml_other_metastatic_site - mostly empty, must only be for a small subset
# xml_metastatic_site_list - mostly empty, must only be for a small subset

# Want to find a tissue source with enough normal samples without any obvious 
# inconsistencies so that I can be confident that there shouldn't be any differences in 
# mean or dispersion between randomly assigned groups. Probably need around 100 samples 
# total for differential dispersion, and probably would be good if I could use the same 
# source for real comparisons, so something with subgroups that make sense to analyse.
table(info$gdc_cases.samples.sample_type)
table(info$cgc_sample_sample_type)
# Both have 740 Solid Tissue Normal
table(info$gdc_cases.project.project_id[which(
  info$gdc_cases.samples.sample_type == "Solid Tissue Normal")])
# Most projects have a small number of normals. TCGA-BRCA most with 112, then TCGA-KIRC 
# with 72.
table(info$gdc_cases.project.name[which(
  info$gdc_cases.project.project_id == "TCGA-BRCA")])
# 1246 Breast Invasive Carcinoma
table(info$gdc_cases.project.name[which(
  info$gdc_cases.project.project_id == "TCGA-KIRC")])
# 616 Kidney Renal Clear Cell Carcinoma

## Extract data for normal samples from breast cancer project
breast.normal <- rse_gene[, which(
  colData(rse_gene)$gdc_cases.samples.sample_type == "Solid Tissue Normal" & 
  colData(rse_gene)$gdc_cases.project.project_id == "TCGA-BRCA")]
dim(breast.normal)
# 58037 112
brno.info <- colData(breast.normal)
table(brno.info$gdc_cases.samples.annotations.category) # empty
table(brno.info$gdc_cases.samples.oct_embedded) # false 45, true 67
table(brno.info$gdc_cases.samples.is_ffpe) # false 112
table(brno.info$cgc_sample_is_ffpe) # NO 112
table(brno.info$gdc_cases.samples.portions.annotations.notes) # empty
table(brno.info$gdc_cases.samples.annotations.category) # empty
table(brno.info$gdc_cases.samples.annotations.notes) # empty
table(as.character(brno.info$xml_histological_type)) # 91 of one type
table(brno.info$xml_tumor_type) # empty
table(brno.info$xml_diagnosis_subtype) # empty
table(brno.info$xml_person_neoplasm_cancer_status) # TUMOR FREE 85, WITH TUMOR 16
table(brno.info$gdc_cases.diagnoses.tumor_stage) # Very variable
table(brno.info$cgc_case_tumor_status) # TUMOR FREE 83, WITH TUMOR 16
table(brno.info$gdc_platform) # Illumina HiSeq 112
table(brno.info$gdc_cases.diagnoses.primary_diagnosis) # c50.9 106, a few others
table(brno.info$cgc_case_clinical_stage) # empty
table(brno.info$cgc_case_icd_o3histology) # 88 8500/3, a few others
table(brno.info$cgc_case_other_histological_diagnosis) # few entries
table(brno.info$xml_other_dx) # 108 No, 4 Yes
table(brno.info$xml_distant_metastasis_present_ind2) # 57 NO, 4 YES
sum(brno.info$xml_histological_type == "Infiltrating Ductal Carcinoma" & 
      brno.info$gdc_cases.diagnoses.primary_diagnosis == "c50.9" & 
      brno.info$cgc_case_icd_o3histology == "8500/3")
# 83


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
dim(blood) # 28037 456
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




