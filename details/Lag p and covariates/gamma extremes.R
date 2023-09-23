
# compute extremes of gamma parameter for STPNAR model

stpnar_ex <- function(y, W, pl=0.1, pu=0.9){
  x <- mean(W%*%y)
  L=-log(pu)/x^2
  U=-log(pl)/x^2
  return(list(L=L, U=U))
}

# compute extremes of gamma parameter for TPNAR model

tpnar_ex <- function(y, W, pl=0.2, pu=0.8){
  X <- W%*%y
  N <- ncol(W)
  qq <- matrix(0, 2, N)
  for (k in 1:N){
    qq[,k] <- as.matrix(quantile(X[k,], prob=c(pl, pu)))
  }
  return(list(L=mean(qq[1,]), U=mean(qq[2,])))
}


