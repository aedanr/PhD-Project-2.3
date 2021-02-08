# Change mean by factor of 1.5 plus random draw from exponential(1), randomly up and down for 
# half of genes.
# For each gene, generates factor to change by, then adds mean * (factor - 1) to all counts, 
# which changes sample mean by the correct factor without changing variance. Then generates a 
# factor to change variance by without changing mean so that dispersion (estimated by moment 
# estimator) matches original and does this by adding sqrt(factor) * (distance of mean from 
# counts), but with any negative values obtained set to zero.
diff.mean <- function(x, integer=FALSE) {
  n <- nrow(x)
  m <- 1.5 + rexp(n, 1)
  dec <- sample(seq_len(n), n/2)
  m[dec] <- 1 / m[dec]
  for (i in 1:n) {
    x[i,] <- x[i,] + mean(x[i,]) * (m[i] - 1)
    c <- (m[i] * mean(x[i,]) + m[i]^2 * max(var(x[i,]) - mean(x[i,]), 0)) / max(var(x[i,]), 1e-6)
    x[i,] <- pmax(0, sqrt(c) * (x[i,] - mean(x[i,])) + mean(x[i,]))
  }
  if(integer) {
    x <- round(x)
  }
  return(list(counts = x, 
              FC = m))
}

# Changes dispersion by factor of 1.5 plus random draw from exponenial(1), randomly up and down 
# for half of genes.
# For each gene, generates factor to change by, calculates appropriate factor to change variance 
# by without changing mean and does this by adding sqrt(factor) * (distance of mean from counts), 
# but with any negative values obtained set to zero.
diff.disp <- function(x, integer=FALSE) {
  n <- nrow(x)
  d <- 1.5 + rexp(n, 1)
  dec <- sample(seq_len(n), n/2)
  d[dec] <- 1 / d[dec]
  for (i in 1:n) {
    c <- (mean(x[i,]) + d[i] * max(var(x[i,]) - mean(x[i,]), 0)) / max(var(x[i,]), 1e-6)
    x[i,] <- pmax(0, sqrt(c) * (x[i,] - mean(x[i,])) + mean(x[i,]))
  }
  if(integer) {
    x <- round(x)
  }
  return(list(counts = x, 
              FC = d))
}
