# BFDR following Ventrucci 2011, procedure taken from cjBitSeq paper (Papasamoulis 2018) 
# and ShrinkBayes::BFDR

bfdr <- function(x) {
  ord <- sort(x, decreasing=T, index.return=T)
  bfdr <- cumsum(1-ord$x)/(1:length(x))
  bfdr[ord$ix] <- bfdr
  return(bfdr)
}
