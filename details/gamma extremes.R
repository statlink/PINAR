# compute extremes of gamma parameter for STPNAR model
stnar_ex <- function(y, W, pl = 0.1, pu = 0.9) {
  x <- mean(W %*% y)
  c( -log(pu) / x^2, -log(pl) / x^2 )
}

# compute extremes of gamma parameter for TPNAR model

tnar_ex <- function(y, W, pl = 0.2, pu = 0.8) {
  X <- W %*% y
  qq2 <- Rfast2::rowQuantile( X, probs = c(pl, pu) )
  Rfast::colmeans(qq2)
  
  N <- dim(W)[1]
  qq <- matrix(0, 2, N)
  for (k in 1:N) {
    qq[, k] <- as.matrix( quantile(X[k, ], prob = c(pl, pu)) )
  }
  c( mean(qq[1,] ), mean(qq[2, ]) ) 
}


