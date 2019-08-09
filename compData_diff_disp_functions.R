# Simulate minimal dataset with desired number of genes to get (mean, dispersion) pairs, then manually 
# create list of dispersions for second condition by multiplying appropriate dispersions by (1.5 + rexp(rate=1)), 
# then use these means and dispersions to generate data with desired differential dispersion.
# simulate.DE.DD.data gives 5% DE only, 5% DD only, 5% both, half up half down evenly distributed for each.
# simulate.DD.only.data gives 5% DD only, half up half down.

simulate.DE.DD.data <- function(dataset, n.vars, samples.per.cond, filter.threshold.mediancpm=0.5) {
  temp.data <- generateSyntheticData(dataset='temp', 
                                     n.vars=n.vars, 
                                     samples.per.cond=2, 
                                     n.diffexp=0, 
                                     filter.threshold.total=0)
  dispersion.factors <- c(rep(c(rep(1,0.025*n.vars), 1.5 + rexp(0.0125*n.vars,1), 1/(1.5 + rexp(0.0125*n.vars,1))), 2), 
                          1.5 + rexp(0.025*n.vars,1), 1/(1.5 + rexp(0.025*n.vars,1)), rep(1,0.85*n.vars))
  dispersions <- cbind(temp.data@variable.annotations$truedispersions.S1, 
                       temp.data@variable.annotations$truedispersions.S1 * dispersion.factors)
  means <- temp.data@variable.annotations$truemeans.S1
  real.data <- generateSyntheticData(dataset=dataset, 
                                     n.vars=n.vars, 
                                     samples.per.cond=samples.per.cond, 
                                     n.diffexp=0.1*n.vars, 
                                     fraction.upregulated=0.5, 
                                     relmeans=means, 
                                     dispersions=dispersions, 
                                     filter.threshold.mediancpm=filter.threshold.mediancpm)
  real.data@variable.annotations$updispersion <- as.numeric(real.data@variable.annotations$truedispersions.S2 > 
                                                              real.data@variable.annotations$truedispersions.S1)
  real.data@variable.annotations$downdispersion <- as.numeric(real.data@variable.annotations$truedispersions.S2 < 
                                                                real.data@variable.annotations$truedispersions.S1)
  real.data@variable.annotations$differential.dispersion <- as.numeric(real.data@variable.annotations$truedispersions.S2 != 
                                                                         real.data@variable.annotations$truedispersions.S1)
  real.data@variable.annotations$truelog2fcdispersion <- log2(real.data@variable.annotations$truedispersions.S2) - 
    log2(real.data@variable.annotations$truedispersions.S1)
  real.data@info.parameters$n.diffdisp <- 0.1*n.vars
  return(real.data)
}

simulate.DD.only.data <- function(dataset, n.vars, samples.per.cond, filter.threshold.mediancpm=0.5) {
  temp.data <- generateSyntheticData(dataset='temp', 
                                     n.vars=n.vars, 
                                     samples.per.cond=2, 
                                     n.diffexp=0, 
                                     filter.threshold.total=0)
  dispersion.factors <- c(1.5 + rexp(0.025*n.vars,1), 1/(1.5 + rexp(0.025*n.vars,1)), rep(1,0.95*n.vars))
  dispersions <- cbind(temp.data@variable.annotations$truedispersions.S1, 
                       temp.data@variable.annotations$truedispersions.S1 * dispersion.factors)
  means <- temp.data@variable.annotations$truemeans.S1
  real.data <- generateSyntheticData(dataset=dataset, 
                                     n.vars=n.vars, 
                                     samples.per.cond=samples.per.cond, 
                                     n.diffexp=0, 
                                     fraction.upregulated=0.5, 
                                     relmeans=means, 
                                     dispersions=dispersions, 
                                     filter.threshold.mediancpm=filter.threshold.mediancpm)
  real.data@variable.annotations$updispersion <- as.numeric(real.data@variable.annotations$truedispersions.S2 > 
                                                              real.data@variable.annotations$truedispersions.S1)
  real.data@variable.annotations$downdispersion <- as.numeric(real.data@variable.annotations$truedispersions.S2 < 
                                                                real.data@variable.annotations$truedispersions.S1)
  real.data@variable.annotations$differential.dispersion <- as.numeric(real.data@variable.annotations$truedispersions.S2 != 
                                                                         real.data@variable.annotations$truedispersions.S1)
  real.data@variable.annotations$truelog2fcdispersion <- log2(real.data@variable.annotations$truedispersions.S2) - 
    log2(real.data@variable.annotations$truedispersions.S1)
  real.data@info.parameters$n.diffdisp <- 0.05*n.vars
  return(real.data)
}

