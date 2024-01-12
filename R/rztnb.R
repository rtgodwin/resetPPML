rztnb <- function(b, alpha, X) {
  X <- as.matrix(X)
  l <- exp(X %*% b)
  theta <- l / alpha
  n <- nrow(X)
  y <- rep(0, n)
  for(i in 1:n) {
    roll <- runif(1)
    probs <- alpha * ((1 / (1 + theta[i])) ^ alpha) * (theta[i] / (1 + theta[i] - (1 + theta[i]) ^ (1 - alpha)))
    k <- 1
    while(probs[k] < roll) {
      k <- k + 1
      probs <- c(probs, probs[k - 1] + ((gamma(alpha + k)) / (gamma(alpha) * gamma(k + 1))) * ((1/(1 + theta[i])) ^ alpha) * ((theta[i] / (1 + theta[i])) ^ k) * (1 / (1 - (1 + theta[i]) ^ (-alpha))))
      if(is.nan(probs[k]) == T) {break}
    }
    y[i] <- k
  }
  return(y)
}
