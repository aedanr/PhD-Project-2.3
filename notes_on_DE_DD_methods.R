library(MDSeq)
# genes in rows, samples in columns
data(sampleData)
dat <- sample.exprs
dim(dat)
group <- sample.pheno$group
# group is a vector, of -1 and 1 in sample data, but presumably doesn't matter

# filter.counts() by default removes genes with mean cpm < 0.05; 0.1 in vignette example.
# cpm calculated using edgeR.
# From vignette, ``highly recommend'' removing lowly expressed genes but no strict 
# recommendation on threshold.
dat.filtered <- filter.counts(dat, mean.cpm.cutoff=0.1)
dim(dat.filtered)

# In vignette, compute TMM normalisation factors, which can be incorporated as 
# normalised counts or as offsets in GLMs in MDSeq. MDSeq function normalize.counts()
# normalises using a range of methods.
# To use normalised counts:
dat.normalized <- normalize.counts(dat.filtered, group=group, method="TMM")
# Check that multiplying by factors from edgeR calcNormFactors gives same results:
libsize <- colSums(dat.filtered) # same as lib.size from DGEList(dat.filtered)
rellibsize <- libsize/exp(mean(log(libsize)))
cnf <- calcNormFactors(dat.filtered, method="TMM") # result is same as norm.factors column added if applied to a DGEList
nf <- cnf*rellibsize
rawnf <- cnf*libsize # this isn't in MDSeq vignette, not sure where I got it from

y <- DGEList(counts=dat.filtered)
y.cnf <- calcNormFactors(y)

dat.filtered[1:5,1:5]
(t(dat.filtered)*(1/y.cnf$samples$norm.factors))[1:5,1:5]
(t(t(dat.filtered)*(1/y.cnf$samples$norm.factors)))[1:5,1:5]
(t(t(dat.filtered)*(1/y.cnf$samples$norm.factors)) - 0.5)[1:5,1:5]
(ceiling(t(t(dat.filtered)*(1/y.cnf$samples$norm.factors)) - 0.5))[1:5,1:5]
(round(t(t(dat.filtered)*(1/y.cnf$samples$norm.factors))))[1:5,1:5]
dat.normalized[1:5,1:5]
# last 3 all same


library(edgeR)
?filterByExpr
y <- estimateCommonDisp(y)
names(y)
y$pseudo.counts[1:5,1:5]
dat.filtered[1:5,1:5]
(ceiling(t(t(dat.filtered)*(1/y.cnf$samples$norm.factors)) - 0.5))[1:5,1:5]
t(t(dat.filtered)*(1/y.cnf$samples$norm.factors))[1:5,1:5]
# pseudo-counts not same as normalised counts


library(baySeq)
?getLibsizes
?getPriors.NB
?getLikelihoods
groups <- list(NDE=rep(1,200), DE=rep(1,200))
replicates <- c(rep("1",200))
cD <- new("countData", data=dat.filtered, groups=groups, replicates=replicates)
libsizes(cD) <- getLibsizes(cD, estimationType="edgeR")
head(libsizes(cD)) # same as rawnf from MDSeq (except can't find rawnf in MDSeq any more)
head(libsize*y.cnf$samples$norm.factors) # also same
# so baySeq library sizes are library sizes normalised using TMM but with no library size scaling
# so counts derived from baySeq library sizes have TM applied but not library size scaling


library(DSS)
?waldTest
?estNormFactors
designs=c(rep(0,100),rep(1,100))
dat.filtered <- matrix(dat.filtered, nrow=nrow(dat.filtered))
y <- newSeqCountSet(dat.filtered, designs)
head(y@normalizationFactor)
y <- estNormFactors(y, method='lr')
head(y@normalizationFactor)
head(cnf)
y2 <- sweep(exprs(y), 2, normalizationFactor(y), FUN="/")
y2[1:5,1:5]
(t(t(exprs(y))/normalizationFactor(y)))[1:5,1:5]
# same as using sweep


library(ShrinkBayes)
DGE <- DGEList(dat.filtered)
libsize <- colSums(dat.filtered)
DGEnorm <- calcNormFactors(DGE)
normfac0 <- DGEnorm$sample[,3]
head(normfac0)
rellibsize <- libsize/exp(mean(log(libsize)))
normfac <- normfac0*rellibsize
head(normfac) # same as nf from MDSeq
pseudocounts <- round(sweep(dat.filtered,2,normfac,"/"))
pseudocounts[1:5,1:5]
# not same as dat.normalized from MDSeq or pseudocounts from edgeR
# normalises by (norm.factors * relative library size)
# so maybe these values are after TMM normalisation and library size scaling, whereas 
# MDSeq only does TMM normalisation (and maybe library size scaling is done within 
# analysis).
head(colSums(pseudocounts))

y.cnf <- calcNormFactors(y)
head(y.cnf$samples$norm.factors) # same as normfac0
cnf <- calcNormFactors(dat.filtered, method="TMM") # same as y.cnf$samples$norm.factors, normfac0
(ceiling(t(t(dat.filtered)*(1/y.cnf$samples$norm.factors)) - 0.5))[1:5,1:5]
nf <- cnf*rellibsize # same as normfac
rawnf <- cnf*libsize



# Raw lib sizes; same as libsize
head(colSums(dat.filtered))
# 41699310 35960775 53352548 39642739 39563573 37303057 
# Raw counts
head(dat.filtered[1,])
# 33      72      94     106      64      73 

# Effective lib sizes (i.e. post-TMM)
head(colSums(dat.filtered)*cnf)
# 43455775 38352446 58796310 43091710 42941054 40225176 
# TMM-normalised counts without library size scaling
head(t(t(dat.filtered)*(1/cnf))[1,])
# 31.66615 67.51006 85.29684 97.51598 58.96615 67.69699 
# Rounded TMM-normalised counts; same as dat.normalized, i.e. MDSeq normalised counts
head(round(t(t(dat.filtered)*(1/cnf)))[1,])
# 32      68      85      98      59      68 

# Lib sizes after raw lib size scaling (without TMM); geometric mean of lib sizes
head(libsize/rellibsize)
# 41293414 41293414 41293414 41293414 41293414 41293414 
# Lib-size scaled counts without TMM
head(t(t(dat.filtered)*(1/rellibsize))[1,])
# 32.67878  82.67691  72.75343 110.41371  66.79828  80.80891 
# Rounded lib-size scaled counts
head(round(t(t(dat.filtered)*(1/rellibsize))[1,]))
# 33      83      73     110      67      81 

# Lib sizes after TMM followed by lib size scaling; geom mean of post-TMM lib sizes
tmm.norm.counts <- t(t(dat.filtered)*(1/cnf))
libsize.tmm <- colSums(tmm.norm.counts)
rellibsize.tmm <- libsize.tmm/exp(mean(log(libsize.tmm)))
head(libsize.tmm/rellibsize.tmm)
# 41293414 41293414 41293414 41293414 41293414 41293414 
# same as lib size scaling without TMM
# TMM-then-lib-size scaled counts
head(t(t(tmm.norm.counts)*(1/rellibsize.tmm))[1,])
# 32.67878  82.67691  72.75343 110.41371  66.79828  80.80891 
# same as lib size scaling without TMM

# Lib sizes after lib size scaling followed by TMM
libscaled.counts <- t(t(dat.filtered)*(1/rellibsize))
cnf.ls <- calcNormFactors(libscaled.counts, method="TMM")
head(colSums(libscaled.counts)*cnf.ls)
# 43030563 44035879 45552562 44878526 44811058 44518136 
# lib-size-then-TMM-normalised counts
head(t(t(libscaled.counts)*(1/cnf.ls))[1,])
# 31.35954  77.52796  65.95101 101.59333  61.55465  74.95542 
# Rounded lib-size-then-TMM-normalised counts
head(round(t(t(libscaled.counts)*(1/cnf.ls)))[1,])
# 31      78      66     102      62      75 
# same as ShrinkBayes counts





### Following code from Soneson & Delorenzi
count.matrix <- dat.filtered
class <- c(rep(0,100),rep(1,100))
#edgeR
edgeR.dgelist = DGEList(counts = count.matrix, group = factor(class))
edgeR.dgelist = calcNormFactors(edgeR.dgelist, method = "TMM")
# DESeq
class <- factor(class)
DESeq.cds = DESeqDataSetFromMatrix(countData = count.matrix, colData = as.data.frame(class), design = ~class)
DESeq.cds = estimateSizeFactors(DESeq.cds)
# baySeq
baySeq.cd = new("countData", data = count.matrix, replicates = class, 
                groups = list(NDE = rep(1, length(class)), DE = class))
libsizes(baySeq.cd) = getLibsizes(baySeq.cd, estimationType = "edgeR")
# limma-voom
nf = calcNormFactors(count.matrix, method = "TMM")
voom.data = voom(count.matrix, design = model.matrix(~factor(class)), lib.size = colSums(count.matrix) * nf)
# ShrinkBayes
nf = calcNormFactors(count.matrix, method = "TMM") * colSums(count.matrix)/exp(mean(log(colSums(count.matrix))))
count.matrix = round(sweep(count.matrix, 2, nf, "/"))

t(t(edgeR.dgelist$counts)/edgeR.dgelist$samples$norm.factors)[1:5,1:5]
#         sample1    sample2    sample3     sample4    sample5
# gene1  31.66615  67.510057   85.29684   97.515980  58.966150
# gene3 912.56097 908.572850 1021.74727  731.369850 900.155138
# gene5   0.00000   5.625838    0.00000    9.199621   0.000000
# gene7 852.10740 951.704275 1255.85988 1169.271797 849.481102
# gene9  74.84727  41.256146  303.98342   21.159128   1.842692
t(t(counts(DESeq.cds))/DESeq.cds$sizeFactor)[1:5,1:5]
#         sample1    sample2  sample3     sample4    sample5
# gene1  31.23425   76.31772  62.9873  100.995945  61.550653
# gene3 900.11440 1027.10925 754.5075  757.469589 939.609184
# gene5   0.00000    6.35981   0.0000    9.527919   0.000000
# gene7 840.48537 1075.86779 927.3875 1210.998550 886.714092
# gene9  73.82642   46.63860 224.4760   21.914215   1.923458
baySeq.cd@data[1:5,1:5]
t(t(baySeq.cd@data)*colSums(dat.filtered)/as.numeric(libsizes(baySeq.cd)))[1:5,1:5]
# same as edgeR
2^(voom.data$E[1:5,1:5])
#           sample1    sample2      sample3    sample4     sample5
# gene1  0.77089868  1.8903618  1.607243699  2.4714730  1.50205903
# gene3 21.89582368 25.2787001 19.159365359 18.4606269 22.76376276
# gene5  0.01150595  0.1694807  0.008503935  0.2436664  0.01164387
# gene7 20.44607393 26.4781020 23.547395774 29.5068348 21.48293724
# gene9  1.80643422  1.1602910  5.706140327  0.5453485  0.05821934
count.matrix[1:5,1:5]
#       sample1 sample2 sample3 sample4 sample5
# gene1      31      78      66     102      62
# gene3     904    1043     791     762     940
# gene5       0       6       0      10       0
# gene7     844    1093     972    1218     887
# gene9      74      47     235      22       2




## Filtering
library(MDSeq)
library(edgeR)
data(sampleData)
group <- sample.pheno$group
dat <- sample.exprs
dim(dat)
# [1] 49722   200
dat.mdseq <- filter.counts(dat, mean.cpm.cutoff=0.1)
dim(dat.mdseq)
# [1] 21992   200
millionreads <- colSums(dat)/1e6
cpm <- t(t(dat)/millionreads)
cpm[1:5,1:5]
meancpm <- rowMeans(cpm)
meancpm[1:5]
sum(meancpm>0.1)
# [1] 21992
# So my method of calculating cpm is correct. But edgeR has cpm() which does the same thing.
keep.edgeR <- filterByExpr(dat, group=group)
dat.edgeR <- dat[keep.edgeR,]
dim(dat.edgeR)
# [1] 19691   200
# Default settings for filterByExpr remove slightly more genes than mean cpm > 0.1 for this data.
# But MDSeq default is 0.05, which would leave 24352.
keep.10 <- rowSums(dat)>10
dat.10 <- dat[keep.10,]
dim(dat.10)
# [1] 37357   200
# Min count across samples way less stringent for this data, but this data has a lot of samples.
keep.cpm.5.2lib <- rowSums(cpm(dat) > 0.5) >= 2
dat.cpm.5.2lib <- dat[keep.cpm.5.2lib,]
dim(dat.cpm.5.2lib)
# [1] 22019   200
# Filtering by at least 2 samples with mean cpm > 0.5 slightly more conservative than 
# mean cpm > 0.1 and edgeR default for this dataset.
keep.cpm.5.2lib.10tot <- (rowSums(cpm(dat) > 0.5) >= 2 & rowSums(dat)>10)
dat.cpm.5.2lib.10tot <- dat[keep.cpm.5.2lib.10tot,]
dim(dat.cpm.5.2lib.10tot)
# [1] 22019   200
# For this dataset with a lot of samples, all genes with mean cpm > 0.5 have total count > 10

# Try on a sample compcodeR dataset
library(compcodeR)
data2 <- generateSyntheticData(dataset='temp', n.vars=2e4, samples.per.cond=2, n.diffexp=1e3, 
                               filter.threshold.total=0, filter.threshold.mediancpm=0)
group2 <- c(0,0,1,1)
dim(data2@count.matrix)
# [1] 20000     4
data2.mdseq <- filter.counts(data2@count.matrix, mean.cpm.cutoff=0.1)
dim(data2.mdseq)
# [1] 19425     4
keep.edgeR <- filterByExpr(data2@count.matrix, group=group2)
data2.edgeR <- data2@count.matrix[keep.edgeR,]
dim(data2.edgeR)
# [1] 15579   4
keep.10 <- rowSums(data2@count.matrix)>10
data2.10 <- data2@count.matrix[keep.10,]
dim(data2.10)
# [1] 18177     4
keep.cpm.5.2lib <- rowSums(cpm(data2@count.matrix) > 0.5) >= 2
data2.cpm.5.2lib <- data2@count.matrix[keep.cpm.5.2lib,]
dim(data2.cpm.5.2lib)
# [1] 16913     4
keep.cpm.5.2lib.10tot <- (rowSums(cpm(data2@count.matrix) > 0.5) >= 2 & rowSums(data2@count.matrix)>10)
data2.cpm.5.2lib.10tot <- data2@count.matrix[keep.cpm.5.2lib.10tot,]
dim(data2.cpm.5.2lib.10tot)
# [1] 16912     4
# For dataset with very few samples, only one gene with mean cpm > 0.5 has total count <= 10
data2med <- generateSyntheticData(dataset='temp', n.vars=2e4, samples.per.cond=2, n.diffexp=1e3, 
                                  filter.threshold.total=0, filter.threshold.mediancpm=0.1)
dim(data2med@count.matrix)
# [1] 18878     4
# Can't compare directly with above results as this is a different dataset, but filtering by median 
# cpm 0.1 looks similar to row count 10 for 2 samples per group.
data2med.5 <- generateSyntheticData(dataset='temp', n.vars=2e4, samples.per.cond=2, n.diffexp=1e3, 
                                    filter.threshold.total=0, filter.threshold.mediancpm=0.5)
dim(data2med.5@count.matrix)
# [1] 16506     4
# Median cpm 0.5 similar to at least 2 samples with cpm > 0.5 (as it must be for 2 samples per group)

data10 <- generateSyntheticData(dataset='temp', n.vars=2e4, samples.per.cond=10, n.diffexp=1e3, 
                                filter.threshold.total=0, filter.threshold.mediancpm=0)
group10 <- c(rep(0,10),rep(1,10))
dim(data10@count.matrix)
# [1] 20000     20
data10.mdseq <- filter.counts(data10@count.matrix, mean.cpm.cutoff=0.1)
dim(data10.mdseq)
# [1] 19726     20
keep.edgeR <- filterByExpr(data10@count.matrix, group=group10)
data10.edgeR <- data10@count.matrix[keep.edgeR,]
dim(data10.edgeR)
# [1] 15202   20
keep.10 <- rowSums(data10@count.matrix)>10
data10.10 <- data10@count.matrix[keep.10,]
dim(data10.10)
# [1] 19960     20
keep.cpm.5.2lib <- rowSums(cpm(data10@count.matrix) > 0.5) >= 2
data10.cpm.5.2lib <- data10@count.matrix[keep.cpm.5.2lib,]
dim(data10.cpm.5.2lib)
# [1] 18838    20
keep.cpm.5.2lib.10tot <- (rowSums(cpm(data10@count.matrix) > 0.5) >= 2 & rowSums(data10@count.matrix)>10)
data10.cpm.5.2lib.10tot <- data10@count.matrix[keep.cpm.5.2lib.10tot,]
dim(data10.cpm.5.2lib.10tot)
# [1] 18838    20
data10 <- generateSyntheticData(dataset='temp', n.vars=2e4, samples.per.cond=10, n.diffexp=1e3, 
                                filter.threshold.total=0, filter.threshold.mediancpm=0.1)
dim(data10@count.matrix)
# [1] 19118    20
# median cpm 0.1 very conservative for 10 samples per group.
data10med.5 <- generateSyntheticData(dataset='temp', n.vars=2e4, samples.per.cond=10, n.diffexp=1e3, 
                                     filter.threshold.total=0, filter.threshold.mediancpm=0.5)
dim(data10med.5@count.matrix)
# [1] 16333    20
# Median cpm 0.5 a lot more strict than at least 2 samples with cpm > 0.5 with 10 samples per group, 
# as it obviously must be, and only a bit more conservative than edgeR default.

# edgeR defaults are way more strict than MDSeq, min count 10 and at least 2 samples with cpm > 0.5, 
# regardless of sample size, and are only a bit more strict than median cpm 0.5.





