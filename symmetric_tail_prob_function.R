# Compute symmetric tail probabilities from posterior samples

tail.prob <- function(x) {
  l <- length(x)
  if (sign(min(x))==sign(max(x))) {p <- 1/l}
  else p <- min(which(sort(x)>0))/l
  if (p>0.5) {p <- 1-p}
  return(2*p)
}
