# This is the function to deal with a 'raw' vector
# Not exported.

# Returns NAs for non-numeric input or all-NA input

hdiVector <- function(object, credMass=0.95, ...) {
  result <- c(NA_real_, NA_real_)
  if(is.numeric(object)) {
    attributes(object) <- NULL
    x <- sort.int(object, method='quick')  # also removes NAs
    n <- length(x)
    if(n > 0) {
      # exclude <- ceiling(n * (1 - credMass)) # Not always the same as...
      exclude <- n - floor(n * credMass)       # Number of values to exclude
      low.poss <- x[1:exclude]             # Possible lower limits...
      upp.poss <- x[(n - exclude + 1):n]   # ... and corresponding upper limits
      best <- which.min(upp.poss - low.poss)      # Combination giving the narrowest interval
      result <- c(low.poss[best], upp.poss[best])
    }
  }
  names(result) <- c("lower", "upper")
  return(result)
}

hdiQuick <- function(x, credMass=0.95, ...) {
  n <- length(x)
  exclude <- n - floor(n * credMass)       # Number of values to exclude
  low.poss <- x[1:exclude]             # Possible lower limits...
  upp.poss <- x[(n - exclude + 1):n]   # ... and corresponding upper limits
  best <- which.min(upp.poss - low.poss)      # Combination giving the narrowest interval
  result <- c(low.poss[best], upp.poss[best])
  names(result) <- c("lower", "upper")
  return(result)
}



# x <- sort.int(rlnorm(n=10000, meanlog=1, sdlog=1), method='quick')
x <- rlnorm(n=10000, meanlog=1, sdlog=1)
plot(density(x))
hdiVector(x, 0.95)
library(HDInterval)
hdi(x)

system.time(for (i in 1:10000) {hdi(x, credMass=0.95)}) # 7.97; 0.28 if already sorted
system.time(for (i in 1:10000) {hdiVector(x, credMass=0.95)}) # 7.79; 0.15 if already sorted
system.time(for (i in 1:10000) {hdiQuick(x, credMass=0.95)}) # 0.12; 0.14 if already sorted
system.time(for (i in 1:10000) {hdi.me(x, credMass=0.95)}) # 0.11; 0.13 if already sorted



hpd.pval.old <- function(x, m=0) {
  require(HDInterval)
  y1 <- numeric(9)
  for (t in 1:9) {
    y1[t] <- sign(hdi(x, credMass = 1 - t/10)[1] - m) == sign(hdi(x, credMass = 1 - t/10)[2] + m)
    if (y1[t] == 1) {break}
  }
  if (sum(y1) == 0) {x1 <- 9} else {x1 <- min(which(y1 == 1)) - 1}
  iter <- x1/10
  y2 <- numeric(9)
  for (t in 1:9) {
    y2[t] <- sign(hdi(x, credMass = 1 - (iter + t/100))[1] - m) == sign(hdi(x, credMass = 1 - (iter + t/100))[2] + m)
    if (y2[t] == 1) {break}
  }
  if (sum(y2) == 0) {x2 <- 9} else {x2 <- min(which(y2 == 1)) - 1}
  iter <- iter + x2/100
  y3 <- numeric(9)
  for (t in 1:9) {
    y3[t] <- sign(hdi(x, credMass = 1 - (iter + t/1000))[1] - m) == sign(hdi(x, credMass = 1 - (iter + t/1000))[2] + m)
    if (y3[t] == 1) {break}
  }
  if (sum(y3) == 0) {x3 <- 9} else {x3 <- min(which(y3 == 1)) - 1}
  iter <- iter + x3/1000
  y4 <- numeric(9)
  for (t in 1:9) {
    y4[t] <- sign(hdi(x, credMass = 1 - (iter + t/10000))[1] - m) == sign(hdi(x, credMass = 1 - (iter + t/10000))[2] + m)
    if (y4[t] == 1) {break}
  }
  if (sum(y4) == 0) {x4 <- 10} else {x4 <- min(which(y4 == 1))}
  iter <- iter + x4/10000
  return(iter)
}


hdi.me <- function(x, credMass=0.95) {
  n <- length(x)
  exclude <- n - floor(n * credMass)   # Number of values to exclude
  low.poss <- x[1:exclude]             # Possible lower limits...
  upp.poss <- x[(n - exclude + 1):n]   # ... and corresponding upper limits
  best <- which.min(upp.poss - low.poss)  # Combination giving the narrowest interval
  result <- c(low.poss[best], upp.poss[best])
  return(result)
}
# Uses code from hdi() from HDInterval, but removed sorting part so can sort
# data once in hpd.pval() instead of sorting for every iteration

hpd.pval <- function(x, m=0) {
  if (max(x) < 0 | min(x) < 0) {return (1 / length(x))}
  else {
    x <- sort.int(x, method='quick')
    y1 <- numeric(9)
    for (t in 1:9) {
      y1[t] <- sign(hdi.me(x, credMass=1-t/10)[1]-m)==sign(hdi.me(x, credMass=1-t/10)[2]+m)
      if (y1[t]==1) {break}
    }
    if (sum(y1)==0) {x1 <- 9} else {x1 <- min(which(y1==1))-1}
    iter <- x1/10
    y2 <- numeric(9)
    for (t in 1:9) {
      y2[t] <- sign(hdi.me(x, credMass=1-(iter + t/100))[1]-m)==sign(hdi.me(x, credMass=1-(iter + t/100))[2]+m)
      if (y2[t]==1) {break}
    }
    if (sum(y2)==0) {x2 <- 9} else {x2 <- min(which(y2==1))-1}
    iter <- iter + x2/100
    y3 <- numeric(9)
    for (t in 1:9) {
      y3[t] <- sign(hdi.me(x, credMass=1-(iter + t/1000))[1]-m)==sign(hdi.me(x, credMass=1-(iter + t/1000))[2]+m)
      if (y3[t]==1) {break}
    }
    if (sum(y3)==0) {x3 <- 9} else {x3 <- min(which(y3==1))-1}
    iter <- iter + x3/1000
    y4 <- numeric(9)
    for (t in 1:9) {
      y4[t] <- sign(hdi.me(x, credMass=1-(iter + t/10000))[1]-m)==sign(hdi.me(x, credMass=1-(iter + t/10000))[2]+m)
      if (y4[t]==1) {break}
    }
    if (sum(y4)==0) {x4 <- 10} else {x4 <- min(which(y4==1))}
    iter <- iter + x4/10000
    return(iter)
  }
}

system.time(for (i in 1:10000) {hpd.pval.old(x)}) # 2.15 if x already sorted, 62 if not
system.time(for (i in 1:10000) {hpd.pval(x)}) # 1.14 if x already sorted, 8.39 if not


library(here)
library(compcodeR)
library(ROCR)
library(edgeR)
library(DESeq2)
source(here('scripts','2019-03-27_lognormal_hmm_adaptive_proposals_three_chains_function.R'))
counts <- generateSyntheticData(dataset='DEonly', n.vars=20000, samples.per.cond=5, 
                                n.diffexp=1000, fraction.upregulated=0.5,
                                filter.threshold.mediancpm=0.5)
DE <- counts@variable.annotations$differential.expression
group <- factor(c(rep(1,5), rep(2,5)))
nf.TMM <- calcNormFactors(counts@count.matrix)
norm.TMM <- t(t(counts@count.matrix) / nf.TMM)
dat.DESeq <- DESeqDataSetFromMatrix(countData=counts@count.matrix, colData=data.frame(group), 
                                    design=~group)
dat.DESeq <- estimateSizeFactors(dat.DESeq)
nf.DESeq <- dat.DESeq$sizeFactor
rm(dat.DESeq)
norm.DESeq <- t(t(counts@count.matrix) / nf.DESeq)
lnHM.raw <- ln_hmm_adapt_3_chains(counts=t(counts@count.matrix), groups=c(rep(1,5),rep(2,5)))
rawpost <- log(as.matrix(lnHM.raw$means1)) - log(as.matrix(lnHM.raw$means2))
lnHM.TMM <- ln_hmm_adapt_3_chains(counts=t(norm.TMM), groups=c(rep(1,5),rep(2,5)))
TMMpost <- log(as.matrix(lnHM.TMM$means1)) - log(as.matrix(lnHM.TMM$means2))
lnHM.DESeq <- ln_hmm_adapt_3_chains(counts=t(norm.DESeq), groups=c(rep(1,5),rep(2,5)))
DESeqpost <- log(as.matrix(lnHM.DESeq$means1)) - log(as.matrix(lnHM.DESeq$means2))

system.time(apply(rawpost, 2, hpd.pval.old)) # 368 s
system.time(apply(rawpost, 2, hpd.pval)) # 41 s
system.time(apply(TMMpost, 2, hpd.pval.old)) # 300 s
system.time(apply(TMMpost, 2, hpd.pval)) # 42 s
system.time(apply(DESeqpost, 2, hpd.pval.old)) # 302 s
system.time(apply(DESeqpost, 2, hpd.pval)) # 42 s
# Down from about 5 min to 3/4 min.
# For symmetric intervals, takes less than 2 s.


## Start from max density either side of zero ####
# No point trying intervals e.g. close to 1 if know that 0 is near middle - 
# not possible to have an HDI that doesn't include zero that's bigger than the 
# maximum of the density below zero and the density above zero
hpd.pval.new <- function(x, m=0) {
  if (max(x) < 0 | min(x) < 0) {return (1 / length(x))}
  else {
    x <- sort.int(x, method='quick')
    max.density <- max(which(x>0)[1] / length(x), 1 - (which(x>0)[1] / length(x)))
    y1 <- numeric(9)
    for (t in 1:(ceiling(max.density * 10))) {
      y1[t] <- sign(hdi.me(x, credMass=1-t/10)[1]-m)==sign(hdi.me(x, credMass=1-t/10)[2]+m)
      if (y1[t]==1) {break}
    }
    if (sum(y1)==0) {x1 <- 9} else {x1 <- min(which(y1==1))-1}
    iter <- x1/10
    y2 <- numeric(9)
    for (t in 1:9) {
      y2[t] <- sign(hdi.me(x, credMass=1-(iter + t/100))[1]-m)==sign(hdi.me(x, credMass=1-(iter + t/100))[2]+m)
      if (y2[t]==1) {break}
    }
    if (sum(y2)==0) {x2 <- 9} else {x2 <- min(which(y2==1))-1}
    iter <- iter + x2/100
    y3 <- numeric(9)
    for (t in 1:9) {
      y3[t] <- sign(hdi.me(x, credMass=1-(iter + t/1000))[1]-m)==sign(hdi.me(x, credMass=1-(iter + t/1000))[2]+m)
      if (y3[t]==1) {break}
    }
    if (sum(y3)==0) {x3 <- 9} else {x3 <- min(which(y3==1))-1}
    iter <- iter + x3/1000
    y4 <- numeric(9)
    for (t in 1:9) {
      y4[t] <- sign(hdi.me(x, credMass=1-(iter + t/10000))[1]-m)==sign(hdi.me(x, credMass=1-(iter + t/10000))[2]+m)
      if (y4[t]==1) {break}
    }
    if (sum(y4)==0) {x4 <- 10} else {x4 <- min(which(y4==1))}
    iter <- iter + x4/10000
    return(iter)
  }
}

x <- rlnorm(n=10000, meanlog=1, sdlog=1) - runif(1, min=0, max=25)
system.time(for (i in 1:10000) {hpd.pval(rlnorm(n=10000, meanlog=1, sdlog=1) - runif(1, min=0, max=25))}) # 34.1 s
system.time(for (i in 1:10000) {hpd.pval.new(rlnorm(n=10000, meanlog=1, sdlog=1) - runif(1, min=0, max=25))}) # 34.6 s
# Marginally longer for new version. Time taken to get max density seems to outweigh time saved by 
# not computing higher densities.


## Try binary chop ####
hpd.pval.bin <- function(x, m=0) {
if (max(x) < 0 | min(x) < 0) {return (1 / length(x))}
  else {
    x <- sort.int(x, method='quick')
    max.density <- max(which(x>0)[1] / length(x), 1 - (which(x>0)[1] / length(x)))
    jump <- max.density / 2
    density <- max.density / 2
    while (jump > 1/length(x)) {
      jump <- jump / 2
      excludes.zero <- sign(hdi.me(x, credMass = density)[1] - m) == sign(hdi.me(x, credMass = density)[2] + m)
      if (excludes.zero) {density <- density + jump}
      else {density <- density - jump}
    }
    return(1 - density)
  }
}

system.time(for (i in 1:10000) {hpd.pval(rlnorm(n=10000, meanlog=1, sdlog=1) - runif(1, min=0, max=25))}) # 14.42
system.time(for (i in 1:10000) {hpd.pval.new(rlnorm(n=10000, meanlog=1, sdlog=1) - runif(1, min=0, max=25))}) # 14.41
system.time(for (i in 1:10000) {hpd.pval.bin(rlnorm(n=10000, meanlog=1, sdlog=1) - runif(1, min=0, max=25))}) # 14.70
# No real difference between any of the methods.
system.time(for (i in 1:10000) {hpd.pval(rlnorm(n=10000, meanlog=1, sdlog=0.1) - runif(1, min=1.5, max=4))}) # 14.98
system.time(for (i in 1:10000) {hpd.pval.new(rlnorm(n=10000, meanlog=1, sdlog=0.1) - runif(1, min=1.5, max=4))}) # 15.62
system.time(for (i in 1:10000) {hpd.pval.bin(rlnorm(n=10000, meanlog=1, sdlog=0.1) - runif(1, min=1.5, max=4))}) # 16.20
# Adding first line to return 1/length if whole sample is above or below zero reduces average time by about half 
# for these examples. May not reduce as much for real data but should help.



