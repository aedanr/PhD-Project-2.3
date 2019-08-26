library(here)
DEDD2 <- readRDS(here('Results/Dispersion estimation results Aug 2019','mse.disp.DEDD2.rds'))
DEDD5 <- readRDS(here('Results/Dispersion estimation results Aug 2019','mse.disp.DEDD5.rds'))
DEDD10 <- readRDS(here('Results/Dispersion estimation results Aug 2019','mse.disp.DEDD10.rds'))

names(DEDD2$raw.MSEs)
names(DEDD2$raw.MSEs$noDD.noDE)

par(mfrow=c(2,1), mar=c(3,2,1,1))
boxplot(DEDD2$raw.MSEs$noDD.noDE, cex.axis=0.8)
boxplot(DEDD2$raw.MSEs$noDD.DE, cex.axis=0.8)
colMeans(DEDD2$raw.MSEs$noDD.noDE)
colMeans(DEDD2$raw.MSEs$noDD.DE)

par(mfrow=c(2,1), mar=c(3,2,1,1))
boxplot(DEDD5$raw.MSEs$noDD.noDE, cex.axis=0.8)
boxplot(DEDD5$raw.MSEs$noDD.DE, cex.axis=0.8)
colMeans(DEDD5$raw.MSEs$noDD.noDE)
colMeans(DEDD5$raw.MSEs$noDD.DE)

par(mfrow=c(2,1), mar=c(3,2,1,1))
boxplot(DEDD10$raw.MSEs$noDD.noDE, cex.axis=0.8)
boxplot(DEDD10$raw.MSEs$noDD.DE, cex.axis=0.8)
colMeans(DEDD10$raw.MSEs$noDD.noDE)
colMeans(DEDD10$raw.MSEs$noDD.DE)

# Exclude edgeR trend because apart from DEDD2 it's way worse than the rest and so makes 
# boxplots more difficult to compare for other methods.



# ASMR results for 4 samples - need to remember that results for DEDD2 are actually for 
# 4 samples, etc, so HMs are best for up to (at least) 4 samples, not just 2.
# Can't compare for 2 samples because most methods won't allow as would be estimating
# for 1 sample per group. Could do 5 but seems a bit pointless since I've alread got 4. 
# Then next up is 10, where edgeR is better than HMs. Is there any point in doing, say, 8? 
# Probably not just for the sake of saying my method is better for up to x samples.
# May be best to discuss dispersion estimation completely separately from DE etc, at least 
# as the main comparison. Still should also try my dispersion estimates in edgeR/others, 
# and of course try different normalisation.
MSEdispresults <- readRDS(here('Results/Preliminary results for ASMR 2019', 
                               '2019-05-30_MSE_disp_comparison_4_samples_25_runs.rds'))
library(RColorBrewer)
par(mfrow=c(1,2), mgp=c(0.8,0.5,0), mar=c(1,3,2,1))
colours5 <- brewer.pal(5,"Accent")
colours7 <- brewer.pal(7,"Accent")
boxplot(MSEdispresults, xaxt='n', yaxt='n', main="Dispersion",cex.main=2, col=colours7, 
        ylab='MSE', cex.lab=2)
legend(x='topleft', legend=c('Hierarchical model v1', 'Hierarchical model v2', 
                             'edgeR tag', 'edgeR trend', 'DESeq2', 
                             'DSS no trend', 'DSS trend'), fill=colours7)
colMeans(MSEdispresults)
apply(MSEdispresults,2,median)

