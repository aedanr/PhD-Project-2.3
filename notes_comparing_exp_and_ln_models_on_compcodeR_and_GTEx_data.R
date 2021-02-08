########################
#### compcodeR data ####

#### Differential dispersion ####
### AUC ###
# Essentially no differences in AUC, except for small samples, higher with log and for expHM.
### pAUC ###
# Essentially no differences in pAUC, except for small samples, higher with log and for expHM.
### FDR ###
# Can only judge FDR for 20 and 50 samples per group. For 20, all are pretty bad, especially for DD only, lnHM better than expHM, 
# and no real difference between raw and log. For 50, lnHM generally closer to 0.05 than expHM, and less frequently over 0.05, and 
# untranslated generally closer to 0.05 than log, and less frequently over 0.05.
### TPR ###
# Can only judge TPR for 20 and 50 samples per group. Slightly higher for expHM than lnHM, and for log than untranslated, i.e. 
# higher TPR for situations where FDR control is poorer.
### False discovery plots ###
# FDR curves indistinguishable.

#### Differential expression ####
### AUC ###
# Essentially no differences in AUC except for small samples, for which lnHM is better than expHM and log is better than 
# untranslated.
### pAUC ###
# Essentially no difference in pAUC except for small samples, for which lnHM is better than expHM and untranslated is better than 
# log.
### FDR ###
# With differences in mean only, lnHM and untranslated give lower FDRs and closer on average to 0.05 than lnHM and log. With 
# differences in mean and dispersion, there's no real difference between expHM and lnHM, and log is generally closer to 0.05 than 
# untranslated but more often exceeds 0.05, although still with generally acceptable FDRs except for a small number of cases with 5 
# samples per group.
### TPR ###
# Can only judge TPR for at least 5 samples per group. Always higher for log than untranslated, and higher for lnHM than expHM 
# except for 20 and 50 samples per group, where expHM is very slightly higher.
### False discovery plots ###
# Virtually indistinguishable for 10 samples per group upwards. For 2 and 5, lnHM slightly better than expHM, and untranslated very 
# slightly better than log.

#### Differential distribution ####
### AUC ###
# AUC virtually always higher for lnHMM, but very small differences except for 2 to 10 samples per group and with differences in 
# dispersion only for 20 and 50 samples per group.
### pAUC ###
# Essentially no differences in pAUC, except lnHMM lower than expHMM for 50 samples per group and slightly higher for 5 samples per 
# group.
### FDR ###
# Can't really judge FDRs at all as they're all extremely low except for 2 samples per group and for 5 and 10 samples per group with 
# differences in dispersion only, where they're either undefined or close to 1. lnHMM looks to generally be better than expHMM in 
# that with 10 or more samples per group the FDRs are slightly higher, but still far below 0.05. Posterior threshold also generally 
# gives FDRs slightly closer to 0.05 than BFDR and 0.5 posterior probability threshold, but posterior threshold isn't aiming at any 
# particular FDR, so can really just say that it isn't giving problematically high FDRs. BFDR is trying to control FDR at 0.05, and 
# is very conservative.
### TPR ###
# lnHMM always gives higher TPRs than expHMM. Posterior threshold gives higher power than other methods in all situations except for 
# 50 samples per group with expHMM with differences in mean only or mean and dispersion, where BFDR is slightly higher, but lnHMM 
# with posterior threshold is always higher than expHMM with any method. BFDR is far worse than the others for 2, 50 and 10 samples 
# per group, slightly worse for 20 samples per group, and better than 0.5 posterior probability threshold for 50 samples per group.
### False discovery plots ###
# FDRs almost always lower for lnHMM than expHMM. Very little difference at low end in all cases except with 2 samples per group. At 
# higher end, lnHMM better than expHMM with differences in dispersion only or mean only, but expHMM better than lnHMM with 
# differences in mean and dispersion.


###################
#### GTEx data ####

### Differential dispersion ####
## AUC ####
# Generally very little difference between methods, but nearly always very slightly higher AUCs for log than untransformed and for 
# expHM than lnHM.
## pAUC ####
# Consistently slightly higher pAUCs for log than untransformed. No consistent pattern for expHM v lnHM - lnHM better for blood 
# data, expHM better for muscle.
## FDR ####
# Generally untransformed closer to 0.05 than log, and lnHM closer to 0.05 than expHM, but all really bad.
## TPR ####
# TPR consistently higher for log than untransformed. lnHM generally better than expHM for small samples, and expHM generally 
# slightly better than lnHM for large samples (20 or 50 per group).
## False discovery plots ####
# FDR curves virtually indistinguishable, just very slightly better for lnHM with 5 or 10 samples per group for blood, very 
# slightly better for expHM for some sample sizes for muscle, and generally very very better for log than untransformed.

### Differential expression ####
## AUC ####
# AUCs higher for untranslated than log for 2 and 5 samples per group, higher for log for 20 and 50. Consistently higher for lnHM 
# than expHM.
## pAUC ####
# Partial AUCs higher for untranslated than log, except for 50 samples per group, and consistently higher for lnHM than expHM.
## FDR ####
# FDRs always closer to 0.05 for untranslated than log, and generally closer to 0.05 for lnHM than expHM, but always really bad.
## TPR ####
# TPRs always higher for log than untransformed, and generally higher for lnHM than expHM.
## False discovery plots ####
# FDR curves lower for untransformed than log, and lower for lnHM than expHM. Difference is bigger for untransformed v log than for 
# expHM v lnHM for small samples, and all differences are nearly indistinguishable for 20 and 50 samples per group.

### Differential distribution ####
## AUC ####
# AUCs always higher for expHMM than lnHMM with differences in dispersion only, and nearly always with differences in mean and 
# dispersion.
## pAUC ####
# Partial AUCs higher for expHMM than lnHMM except  for small sample sizes with differences in mean only or mean and dispersion.
## FDR ####
# BFDR generally closer to 0.05 than posterior threshold method, although posterior threshold method doesn't aim for any particular 
# FDR. FDRs are pretty poor generally and don't improve with sample size. FDRs nearly always closer to 0.05 for expHMM than lnHMM.
## TPR ####
# Consistently highest TPRs using posterior threshold method except sometimes for 50 samples per group, and higher for lnHMM than 
# expHMM, but at expense of higher FDRs.
## False discovery plots ####
# FDR curves consistently lower for expHMM than lnHMM for differences in dispersion only, but very little difference between them 
# for differences in mean only or mean and dispersion.
