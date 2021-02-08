library(here)
library(compcodeR)

for (i in 1:50) {
  for (j in c(2,5,10,20,50)) {
    for (k in c("DE", "DD", "DEDD")) {
      res <- readRDS(here("Results/compcodeR DE, DD, DEDD results Feb 2020", 
                          paste0("results.", k, j, ".", i, ".rds")))
      # Fix names of elements with typos
      names(res) <- gsub("H.tmmM", "HM.tmm", names(res), fixed=T)
      names(res) <- gsub("H.rleM", "HM.rle", names(res), fixed=T)
      # Add BH-adjusted FDR values that weren't saved for HMs
      p.names <- grep("HM.", names(res), value=T, fixed=T)
      q.names <- gsub("p.", "q.", p.names, fixed=T)
      q.names <- gsub("disq", "disp", q.names, fixed=T)
      for (index in (1:length(p.names))) {
        res[[q.names[index]]] <- p.adjust(get(p.names[index], res), method="BH")
      }
      # Fix incorrect DEDD computation
      if (k == "DEDD") {
        res$DEDD <- as.numeric(res$DE == 1 | res$DD == 1)
      }
      saveRDS(res, here("Results/temp_quarantine", 
                        paste0("results.", k, j, ".", i, ".rds")))
    }
  }
}

