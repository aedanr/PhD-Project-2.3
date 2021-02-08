library(recount)
library(here)
library(edgeR)
load(here("recount data/SRP001540 (Pickrell)", "rse_gene.Rdata"))
pickrell <- rse_gene
pickrell <- scale_counts(pickrell)
load(here("recount data/SRP001563 (Cheung)", "rse_gene.Rdata"))
cheung <- rse_gene
cheung <- scale_counts(cheung)
rm(rse_gene)


# Change mean by a factor of m without changing dispersion
# (change mean without changing variance then change variance without changing mean 
# such that dispersion matches original dispersion)
diff.mean <- function(x) {
  m <- 1.5 + rexp(1)
  new.x <- x + mean(x) * (m - 1)
  c <- (m * mean(x) + m^2 * (var(x) - mean(x))) / var(x)
  new.new.x <- sqrt(c) * (new.x - mean(new.x)) + mean(new.x)
  return(new.new.x)
}

# Change dispersion by factor of d without changing mean
# (Change variance without changing mean such that dispersion changes by factor of d)
diff.disp <- function(x) {
  d <- 1.5 + rexp(1)
  c <- (mean(x) + d * (var(x) - mean(x))) / var(x)
  new.x <- sqrt(c) * (x - mean(x)) + mean(x)
}


# Now need to make sure always returns integer counts > 0 and that var(x) - mean(x) 
# is never negative, while maintaining distributions of changes as much as possible, 
# and allow decreases as well as increases.
diff.mean <- function(x, direction, m = 1.5 + rexp(1)) {
  if (direction == "down") {
    m <- 1/m
  }
  new.x <- x + mean(x) * (m - 1)
  c <- (m * mean(x) + m^2 * max(var(x) - mean(x), 0)) / var(x)
  new.new.x <- pmax(0, round(sqrt(c) * (new.x - mean(new.x)) + mean(new.x)))
  return(new.new.x)
}
diff.disp <- function(x, direction, d=1.5 + rexp(1)) {
  if (direction == "down") {
    d <- 1/d
  }
  c <- (mean(x) + d * (var(x) - mean(x))) / var(x)
  new.x <- pmax(0, round(sqrt(c) * (x - mean(x)) + mean(x)))
  return(new.x)
}

# Re-setting negative counts to zero has potential to change mean or dispersion 
# when they're not supposed to be changed. To see whether this is enough of a 
# problem to have to find a way to avoid it, need to do some testing with 
# realistic data. Also want to make functions return factors that means and 
# dispersions are changed by, and to work efficiently on matrices rather than 
# vectors for individual genes.

diff.mean <- function(x) {
  n <- nrow(x)
  m <- 1.5 + rexp(n, 1)
  dec <- sample(seq_len(n), n/2)
  m[dec] <- 1 / m[dec]
  for (i in 1:n) {
    x[i,] <- x[i,] + mean(x[i,]) * (m[i] - 1)
    c <- (m[i] * mean(x[i,]) + m[i]^2 * max(var(x[i,]) - mean(x[i,]), 0)) / max(var(x[i,]), 1e-6)
    x[i,] <- pmax(0, sqrt(c) * (x[i,] - mean(x[i,])) + mean(x[i,]))
  }
  return(list(counts = x, 
              FC = m))
}
diff.disp <- function(x) {
  n <- nrow(x)
  d <- 1.5 + rexp(n, 1)
  dec <- sample(seq_len(n), n/2)
  d[dec] <- 1 / d[dec]
  for (i in 1:n) {
    c <- (mean(x[i,]) + d[i] * max(var(x[i,]) - mean(x[i,]), 0)) / max(var(x[i,]), 1e-6)
    x[i,] <- pmax(0, sqrt(c) * (x[i,] - mean(x[i,])) + mean(x[i,]))
  }
  return(list(counts = x, 
              FC = d))
}

dim(pickrell)
# 58037 160
group <- factor(levels=c("1","2")); group[1:80] <- "1"; group[81:160] <- "2"; group <- sample(group)
pickrell.filtered <- assay(pickrell[which(rowMeans(apply(assay(pickrell), 2, function(x) 1e6*x/sum(x))) > 0.5), ])
DE <- numeric(nrow(pickrell.filtered))
DE[1:ceiling(nrow(pickrell.filtered)/10)] <- 1
DD <- numeric(nrow(pickrell.filtered))
DD[ceiling(nrow(pickrell.filtered)/20):(3*ceiling(nrow(pickrell.filtered)/20))] <- 1
DEplusDD <- as.numeric(DE == 1 & DD == 1)
DEDD <- as.numeric(DE == 1 | DD == 1)
DEonly <- as.numeric(DE == 1 & DD == 0)
DDonly <- as.numeric(DE == 0 & DD == 1)

pickrell.altered <- pickrell.filtered
pickrell.altered.DE <- diff.mean(pickrell.altered[DE==1, group=="2"])
pickrell.altered[DE==1, group=="2"] <- pickrell.altered.DE$counts
FC.mean <- numeric(nrow(pickrell.altered))
FC.mean[DE==1] <- pickrell.altered.DE$FC
pickrell.altered.DD <- diff.disp(pickrell.altered[DD==1, group=="2"])
pickrell.altered[DD==1, group=="2"] <- pickrell.altered.DD$counts
FC.disp <- numeric(nrow(pickrell.altered))
FC.disp[DD==1] <- pickrell.altered.DD$FC
g1 <- pickrell.altered[, group=="1"]
g2 <- pickrell.altered[, group=="2"]
plot(log(rowMeans(g1)), log(rowMeans(g2)))
points(log(rowMeans(g1))[DEplusDD==1], log(rowMeans(g2))[DEplusDD==1], col='blue')
points(log(rowMeans(g1))[DEonly==1], log(rowMeans(g2))[DEonly==1], col='yellow')
points(log(rowMeans(g1))[DDonly==1], log(rowMeans(g2))[DDonly==1], col='red')
lines(c(0,11),c(0,11))
mean(rowMeans(g1)/rowMeans(g2)) # 1.04
mean((rowMeans(g1)/rowMeans(g2))[DE == 1]) # 1.46
mean((rowMeans(g1)/rowMeans(g2))[DEonly == 1]) # 1.46
mean((rowMeans(g1)/rowMeans(g2))[DD == 1]) # 1.22
mean((rowMeans(g1)/rowMeans(g2))[DDonly == 1]) # 0.98
mean((rowMeans(g1)/rowMeans(g2))[DEDD == 0]) # 0.99

filter <- which(rowMeans(apply(g1, 2, function(x) 1e6*x/sum(x))) > 0.5 & 
                  rowMeans(apply(g2, 2, function(x) 1e6*x/sum(x))) > 0.5)
DE <- DE[filter]
DD <- DD[filter]
DEDD <- DEDD[filter]
DEplusDD <- DEplusDD[filter]
DEonly <- DEonly[filter]
DDonly <- DDonly[filter]
FC.mean <- FC.mean[filter]
FC.disp <- FC.disp[filter]
pickrell.altered <- pickrell.altered[filter, ]
g1 <- pickrell.altered[, group=="1"]
g2 <- pickrell.altered[, group=="2"]
plot(log(rowMeans(g1)), log(rowMeans(g2)))
points(log(rowMeans(g1))[DEplusDD==1], log(rowMeans(g2))[DEplusDD==1], col='blue')
points(log(rowMeans(g1))[DEonly == 1], log(rowMeans(g2))[DEonly == 1], col='yellow')
points(log(rowMeans(g1))[DDonly == 1], log(rowMeans(g2))[DDonly == 1], col='red')
lines(c(0,11),c(0,11))
mean(rowMeans(g1)/rowMeans(g2)) # 1.03
mean((rowMeans(g1)/rowMeans(g2))[DE == 1]) # 1.41
mean((rowMeans(g1)/rowMeans(g2))[DEonly == 1]) # 1.39
mean((rowMeans(g1)/rowMeans(g2))[DD == 1]) # 1.19
mean((rowMeans(g1)/rowMeans(g2))[DDonly == 1]) # 0.97
mean((rowMeans(g1)/rowMeans(g2))[DDonly == 1 & FC.disp < 1]) # 0.98
mean((rowMeans(g1)/rowMeans(g2))[DDonly == 1 & FC.disp > 1]) # 0.96
mean((rowMeans(g1)/rowMeans(g2))[DEDD == 0]) # 0.99
plot(density(log(rowMeans(g1)) - log(rowMeans(g2))))
lines(density((log(rowMeans(g1)) - log(rowMeans(g2)))[DEonly == 1]))
lines(density((log(rowMeans(g1)) - log(rowMeans(g2)))[DDonly == 1]))
lines(density((log(rowMeans(g1)) - log(rowMeans(g2)))[DEplusDD == 1]))
abline(v=0)
# Changing dispersion doesn't seem to change mean

pickrell.altered.TMM <- DGEList(counts=pickrell.altered, group=group)
pickrell.altered.TMM <- calcNormFactors(pickrell.altered.TMM)
pickrell.altered.TMM <- estimateDisp(pickrell.altered.TMM, model.matrix(~group))
pickrell.altered.ql <- glmQLFit(pickrell.altered.TMM, model.matrix(~group))
pickrell.altered.ql <- glmQLFTest(pickrell.altered.ql, coef=2)
pickrell.altered.ql <- topTags(pickrell.altered.ql, n=nrow(pickrell.altered.ql), sort.by="none")$table
res.ql <- pickrell.altered.ql$FDR < 0.05
mean(res.ql) # 0.103
mean(res.ql[DE==0]) # 0.007
mean(res.ql[DDonly==1]) # 0.007
mean(res.ql[DEDD==0]) # 0.007
mean(res.ql[DE==1]) # 0.971
mean(res.ql[DEonly==1]) # 0.975
mean(res.ql[DEplusDD==1]) # 0.967
mean(res.ql[DDonly==1 & FC.disp<1]) # 0.015
mean(res.ql[DDonly==1 & FC.disp>1]) # 0

g1v <- apply(g1,1,var)
g2v <- apply(g2,1,var)
g1m <- rowMeans(g1)
g2m <- rowMeans(g2)
g1d <- (g1v-g1m)/g1m^2
g2d <- (g2v-g2m)/g2m^2
plot(density(g1d-g2d), xlim=c(-1,1))
lines(density((g1d-g2d)[DDonly==1]))
lines(density((g1d-g2d)[DEonly==1]))
lines(density((g1d-g2d)[DEplusDD==1]))
abline(v=0)
# Changing mean doesn't seem to change dispersion

