folder <- "Results/GTEx muscle simulated DE, DD, DEDD results March 2020"
missing.DE <- character(10)
missing.DEDD <- character(10)
for (i in 1:10) {
  filename <- paste0("results.muscle_50_set", i, "_DE.rds")
  res <- readRDS(here(folder, filename))
  if (length(grep("DSS.tmm", names(res))) == 0) {
    missing.DE[i] <- paste0(missing.DE[i], "tmm")
  }
  if (length(grep("DSS.rle", names(res))) == 0) {
    missing.DE[i] <- paste0(missing.DE[i], "rle")
  }
  filename <- paste0("results.muscle_50_set", i, "_DEDD.rds")
  res <- readRDS(here(folder, filename))
  if (length(grep("DSS.tmm", names(res))) == 0) {
    missing.DEDD[i] <- paste0(missing.DEDD[i], "tmm")
  }
  if (length(grep("DSS.rle", names(res))) == 0) {
    missing.DEDD[i] <- paste0(missing.DEDD[i], "rle")
  }
}
missing.DE # RLE missing for 4 and 9
missing.DEDD # Both missing for 5 and 9

## Add columns of NAs for DSS.rle for DE4 and DE9
## DSS.rle goes after DSS.tmm, everything after needs to be shifted by 3
res <- readRDS(here(folder, "results.muscle_50_set4_DE.rds"))
length(res) # 74
res$p.DSS.rle <- rep(NA, length(res$p.voom.rle))
res$q.DSS.rle <- rep(NA, length(res$p.voom.rle))
res$lfdr.DSS.rle <- rep(NA, length(res$p.voom.rle))
length(res) # 77
which(names(res) == "lfdr.DSS.tmm") # 30
res[31:77] <- res[c(75:77, 31:74)]
names(res)[31:77] <- names(res)[c(75:77, 31:74)]
saveRDS(res, here(folder, "results.muscle_50_set4_DE.rds"))

res <- readRDS(here(folder, "results.muscle_50_set9_DE.rds"))
length(res) # 74
res$p.DSS.rle <- rep(NA, length(res$p.voom.rle))
res$q.DSS.rle <- rep(NA, length(res$p.voom.rle))
res$lfdr.DSS.rle <- rep(NA, length(res$p.voom.rle))
length(res) # 77
which(names(res) == "lfdr.DSS.tmm") # 30
res[31:77] <- res[c(75:77, 31:74)]
names(res)[31:77] <- names(res)[c(75:77, 31:74)]
saveRDS(res, here(folder, "results.muscle_50_set9_DE.rds"))


## Add columns of NAs for DSS.tmm and DSS.rle for DEDD5 and DEDD9
## DSS.tmm then DSS.rle go after q.voom.rle, everything after needs to be shifted by 6
res <- readRDS(here(folder, "results.muscle_50_set5_DEDD.rds"))
length(res) # 100
res$p.DSS.tmm <- rep(NA, length(res$p.voom.tmm))
res$q.DSS.tmm <- rep(NA, length(res$p.voom.tmm))
res$lfdr.DSS.tmm <- rep(NA, length(res$p.voom.tmm))
res$p.DSS.rle <- rep(NA, length(res$p.voom.rle))
res$q.DSS.rle <- rep(NA, length(res$p.voom.rle))
res$lfdr.DSS.rle <- rep(NA, length(res$p.voom.rle))
length(res) # 106
which(names(res) == "q.voom.rle") # 28
res[29:106] <- res[c(101:106, 29:100)]
names(res)[29:106] <- names(res)[c(101:106, 29:100)]
saveRDS(res, here(folder, "results.muscle_50_set5_DEDD.rds"))

res <- readRDS(here(folder, "results.muscle_50_set9_DEDD.rds"))
length(res) # 100
res$p.DSS.tmm <- rep(NA, length(res$p.voom.tmm))
res$q.DSS.tmm <- rep(NA, length(res$p.voom.tmm))
res$lfdr.DSS.tmm <- rep(NA, length(res$p.voom.tmm))
res$p.DSS.rle <- rep(NA, length(res$p.voom.rle))
res$q.DSS.rle <- rep(NA, length(res$p.voom.rle))
res$lfdr.DSS.rle <- rep(NA, length(res$p.voom.rle))
length(res) # 106
which(names(res) == "q.voom.rle") # 28
res[29:106] <- res[c(101:106, 29:100)]
names(res)[29:106] <- names(res)[c(101:106, 29:100)]
saveRDS(res, here(folder, "results.muscle_50_set9_DEDD.rds"))



