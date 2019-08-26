library(here)
library(compcodeR)
source(here('scripts','2019-05-17_compData_diff_disp_functions.R'))

for (i in 1:50) {
  # name2 <- paste0('DEDD2.',i)
  # name5 <- paste0('DEDD5.',i)
  # name10 <- paste0('DEDD10.',i)
  # name20 <- paste0('DEDD20.',i)
  # name50 <- paste0('DEDD50.',i)
  # counts2 <- simulate.DE.DD.data(dataset=name2, n.vars=20000, samples.per.cond=2)
  # counts5 <- simulate.DE.DD.data(dataset=name5, n.vars=20000, samples.per.cond=5)
  # counts10 <- simulate.DE.DD.data(dataset=name10, n.vars=20000, samples.per.cond=10)
  # counts20 <- simulate.DE.DD.data(dataset=name20, n.vars=20000, samples.per.cond=20)
  # counts50 <- simulate.DE.DD.data(dataset=name50, n.vars=20000, samples.per.cond=50)
  # file2 <- paste0(name2,'.rds')
  # file5 <- paste0(name5,'.rds')
  # file10 <- paste0(name10,'.rds')
  # file20 <- paste0(name20,'.rds')
  # file50 <- paste0(name50,'.rds')
  # saveRDS(counts2, file=file2)
  # saveRDS(counts5, file=file5)
  # saveRDS(counts10, file=file10)
  # saveRDS(counts20, file=file20)
  # saveRDS(counts50, file=file50)
  # name2 <- paste0('DD2.',i)
  # name5 <- paste0('DD5.',i)
  # name10 <- paste0('DD10.',i)
  # name20 <- paste0('DD20.',i)
  # name50 <- paste0('DD50.',i)
  # counts2 <- simulate.DD.only.data(dataset=name2, n.vars=20000, samples.per.cond=2)
  # counts5 <- simulate.DD.only.data(dataset=name5, n.vars=20000, samples.per.cond=5)
  # counts10 <- simulate.DD.only.data(dataset=name10, n.vars=20000, samples.per.cond=10)
  # counts20 <- simulate.DD.only.data(dataset=name20, n.vars=20000, samples.per.cond=20)
  # counts50 <- simulate.DD.only.data(dataset=name50, n.vars=20000, samples.per.cond=50)
  # file2 <- paste0(name2,'.rds')
  # file5 <- paste0(name5,'.rds')
  # file10 <- paste0(name10,'.rds')
  # file20 <- paste0(name20,'.rds')
  # file50 <- paste0(name50,'.rds')
  # saveRDS(counts2, file=file2)
  # saveRDS(counts5, file=file5)
  # saveRDS(counts10, file=file10)
  # saveRDS(counts20, file=file20)
  # saveRDS(counts50, file=file50)
  # name2 <- paste0('DE2.',i)
  # name5 <- paste0('DE5.',i)
  # name10 <- paste0('DE10.',i)
  # name20 <- paste0('DE20.',i)
  # name50 <- paste0('DE50.',i)
  # counts2 <- generateSyntheticData(dataset='name2', n.vars=20000, samples.per.cond=2, n.diffexp=1000,
  #                                  fraction.upregulated=0.5, filter.threshold.mediancpm=0.5)
  # counts5 <- generateSyntheticData(dataset='name5', n.vars=20000, samples.per.cond=5, n.diffexp=1000,
  #                                  fraction.upregulated=0.5, filter.threshold.mediancpm=0.5)
  # counts10 <- generateSyntheticData(dataset='name10', n.vars=20000, samples.per.cond=10, n.diffexp=1000,
  #                                   fraction.upregulated=0.5, filter.threshold.mediancpm=0.5)
  # counts20 <- generateSyntheticData(dataset='name20', n.vars=20000, samples.per.cond=20, n.diffexp=1000,
  #                                   fraction.upregulated=0.5, filter.threshold.mediancpm=0.5)
  # counts50 <- generateSyntheticData(dataset='name50', n.vars=20000, samples.per.cond=50, n.diffexp=1000,
  #                                   fraction.upregulated=0.5, filter.threshold.mediancpm=0.5)
  # file2 <- paste0(name2,'.rds')
  # file5 <- paste0(name5,'.rds')
  # file10 <- paste0(name10,'.rds')
  # file20 <- paste0(name20,'.rds')
  # file50 <- paste0(name50,'.rds')
  # saveRDS(counts2, file=file2)
  # saveRDS(counts5, file=file5)
  # saveRDS(counts10, file=file10)
  # saveRDS(counts20, file=file20)
  # saveRDS(counts50, file=file50)
  name2 <- paste0('nodiff2.',i)
  name5 <- paste0('nodiff5.',i)
  name10 <- paste0('nodiff10.',i)
  name20 <- paste0('nodiff20.',i)
  name50 <- paste0('nodiff50.',i)
  counts2 <- generateSyntheticData(dataset='name2', n.vars=20000, samples.per.cond=2, n.diffexp=0,
                                   fraction.upregulated=0.5, filter.threshold.mediancpm=0.5)
  counts5 <- generateSyntheticData(dataset='name5', n.vars=20000, samples.per.cond=5, n.diffexp=0,
                                   fraction.upregulated=0.5, filter.threshold.mediancpm=0.5)
  counts10 <- generateSyntheticData(dataset='name10', n.vars=20000, samples.per.cond=10, n.diffexp=0,
                                    fraction.upregulated=0.5, filter.threshold.mediancpm=0.5)
  counts20 <- generateSyntheticData(dataset='name20', n.vars=20000, samples.per.cond=20, n.diffexp=0,
                                    fraction.upregulated=0.5, filter.threshold.mediancpm=0.5)
  counts50 <- generateSyntheticData(dataset='name50', n.vars=20000, samples.per.cond=50, n.diffexp=0,
                                    fraction.upregulated=0.5, filter.threshold.mediancpm=0.5)
  file2 <- paste0(name2,'.rds')
  file5 <- paste0(name5,'.rds')
  file10 <- paste0(name10,'.rds')
  file20 <- paste0(name20,'.rds')
  file50 <- paste0(name50,'.rds')
  saveRDS(counts2, file=file2)
  saveRDS(counts5, file=file5)
  saveRDS(counts10, file=file10)
  saveRDS(counts20, file=file20)
  saveRDS(counts50, file=file50)
}


# Load and re-save files using version=2 option to allow backward compatibility with R 3.4.4 on cluster
for (i in 1:50) {
  for (j in c('nodiff')) {
    for (k in c('2.', '5.', '10.', '20.','50.')) {
      dataset <- readRDS(here('Simulated data', paste0(j,k,i,'.rds')))
      saveRDS(dataset, here('Simulated data', paste0(j,k,i,'.rds')), version=2)
    }
  }
}

