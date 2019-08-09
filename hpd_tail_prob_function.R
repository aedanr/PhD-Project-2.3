hpd.pval <- function(x,m=0) {
  require(HDInterval)
  y1 <- numeric(9)
  for (t in 1:9) {
    y1[t] <- sign(hdi(x, credMass=1-t/10)[1]-m)==sign(hdi(x, credMass=1-t/10)[2]+m)
    if (y1[t]==1) {break}
  }
  if (sum(y1)==0) {x1 <- 9} else {x1 <- min(which(y1==1))-1}
  iter <- x1/10
  y2 <- numeric(9)
  for (t in 1:9) {
    y2[t] <- sign(hdi(x, credMass=1-(iter + t/100))[1]-m)==sign(hdi(x, credMass=1-(iter + t/100))[2]+m)
    if (y2[t]==1) {break}
  }
  if (sum(y2)==0) {x2 <- 9} else {x2 <- min(which(y2==1))-1}
  iter <- iter + x2/100
  y3 <- numeric(9)
  for (t in 1:9) {
    y3[t] <- sign(hdi(x, credMass=1-(iter + t/1000))[1]-m)==sign(hdi(x, credMass=1-(iter + t/1000))[2]+m)
    if (y3[t]==1) {break}
  }
  if (sum(y3)==0) {x3 <- 9} else {x3 <- min(which(y3==1))-1}
  iter <- iter + x3/1000
  y4 <- numeric(9)
  for (t in 1:9) {
    y4[t] <- sign(hdi(x, credMass=1-(iter + t/10000))[1]-m)==sign(hdi(x, credMass=1-(iter + t/10000))[2]+m)
    if (y4[t]==1) {break}
  }
  if (sum(y4)==0) {x4 <- 10} else {x4 <- min(which(y4==1))}
  iter <- iter + x4/10000
  return(iter)
}
