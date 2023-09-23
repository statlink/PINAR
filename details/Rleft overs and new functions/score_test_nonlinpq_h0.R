################################################################################
## Function for for testing linearity of Poisson NAR model with p lags, PNAR(p),
## versus nonlinear Intercept Drift (ID) alternative model
## it also includes q non time-varying covariates
################################################################################

################################################################################
################################################################################
##  scor_nonlinpq_h0() --- function for the computation of score of the nonlinear
##  Intercept Drift (ID) PNAR model, with p lag, under the null assumption
##  of linearity; it also includes q non time-varying covariates.
##  Inputs:
##     b = estimated parameters from the linear model, in the following order:
##         (intercept, p network effects, p ar effects, covariates)
##     N = number of nodes on the network
##     TT = temporal sample size
##     y = NxTT multivariate count time series
##     W = NxN row-normalized weighted adjacency matrix describing the network
##     p = number of lags in the model
##     d = lag parameter of nonlinear variable (should be between 1 and p)
##     Z = Nxq matrix of covariates (one for each column), where q is the number of 
##         covariates in the model. They must be non-negative
##  output:
##     ss = vector of quasi score
################################################################################

scor_nonlinpq_h0 <- function(b, N, TT, y, W, p, d, Z) {
  
  m <- 1 + 2 * p + max(0, ncol(Z) ) + 1
  b[m] <- 0

  z <- W %*% y 
  zf <- log1p( z )
  f <-  -b[1] * zf
  wy <- NULL
  for ( ti in (p + 1):TT )  wy <- rbind( wy, cbind(z[, (ti - 1):(ti - p)], y[, (ti - 1):(ti - p)], Z, f[, ti - d]) )
  wy <- cbind(1, wy)
  lambdat <- as.vector( wy %*% b )
  a <- wy * ( as.vector(y[, -c(1:p)]) - lambdat ) / lambdat
  matrix( Rfast::colsums(a), m, 1 )

}

scor_nonlinpq_h0_old <- function(b, N, TT, y, W, p, d, Z) {
  
  m <- 1 + 2 * p + max(0, ncol(Z) ) + 1
  ss <- matrix(0, m, 1)
  b[m] <- 0
  
  for (ti in (p + 1):TT ) {
    f <- -b[1] * log1p( W %*% y[, ti - d] )
    Xt <- cbind( 1, W %*% y[, (ti - 1):(ti - p)], y[, (ti - 1):(ti - p)], Z, f )
    lambdat <- Xt %*%b
    Dt <- diag(1/as.vector(lambdat))
    ss <- ss + t(Xt)%*%Dt%*%(y[,ti]-lambdat)
  }
  
  return(ss)
}


################################################################################
################################################################################
##  outer_nonlinpq_h0() --- function for the computation of information matrix
##  of the nonlinear Intercept Drift (ID) PNAR model, with p lag, under the null assumption
##  of linearity; it also includes q non time-varying covariates.
##  Inputs:
##     b = estimated parameters from the linear model, in the following order:
##         (intercept, p network effects, p ar effects, covariates)
##     N = number of nodes on the network
##     TT = temporal sample size
##     y = NxTT multivariate count time series
##     W = NxN row-normalized weighted adjacency matrix describing the network
##     p = number of lags in the model
##     d = lag parameter of nonlinear variable (should be between 1 and p)
##     Z = Nxq matrix of covariates (one for each column), where q is the number of 
##         covariates in the model. They must be non-negative
##  output:
##     out = information matrix
################################################################################

outer_nonlinpq_h0 <- function(b, N, TT, y, W, p, d, Z){
  
  m <- 1 + 2 * p + max(0, ncol(Z) ) + 1
  b[m] <- 0

  z <- W %*% y 
  zf <- log1p( z )
  f <-  -b[1] * zf
  wy <- NULL
  for ( ti in (p + 1):TT )  wy <- rbind( wy, cbind(z[, (ti - 1):(ti - p)], y[, (ti - 1):(ti - p)], Z, f[, ti - d]) )
  wy <- cbind(1, wy)
  lambdat <- as.vector( wy %*% b )
  a <- wy * ( as.vector(y[, -c(1:p)]) - lambdat ) / lambdat
  out1 <- 0
  k <- rep( 1:c(TT - p), each = N )
  b1 <- rowsum(a, k)
  for ( i in 1: c(TT - p) )  out1 <- out1 + tcrossprod(b1[i, ])
  out1
}

outer_nonlinpq_h0_old <- function(b, N, TT, y, W, p, d, Z){
  
  m <- 1+2*p+max(0,ncol(Z))+1
  out <- matrix(0, nrow=m, ncol=m)
  b[m] <- 0
  
  for(t in (p+1):TT){
    f <- -b[1]*log1p(W%*%y[,t-d])
    Xt <- cbind(1, W%*%y[,(t-1):(t-p)], y[,(t-1):(t-p)], Z, f)
    lambdat <- Xt%*%b
    Dt <- diag(1/as.vector(lambdat))
    ss <- t(Xt)%*%Dt%*%(y[,t]-lambdat)
    out <- out + ss%*%t(ss)
  }
  
  return(out)
}


################################################################################
################################################################################
##  hess_nonlinpq_h0() --- function for the computation of Hessian matrix
##  of the nonlinear Intercept Drift (ID) PNAR model, with p lag, 
##  under the null assumption of linearity; it also includes q non
##  time-varying covariates.
##  Inputs:
##     b = estimated parameters from the linear model, in the following order:
##         (intercept, p network effects, p ar effects, covariates)
##     N = number of nodes on the network
##     TT = temporal sample size
##     y = NxTT multivariate count time series
##     W = NxN row-normalized weighted adjacency matrix describing the network
##     p = number of lags in the model
##     d = lag parameter of nonlinear variable (should be between 1 and p)
##     Z = Nxq matrix of covariates (one for each column), where q is the number of 
##         covariates in the model. They must be non-negative
##  output:
##     hh = (-1)*Hessian matrix
################################################################################


hess_nonlinpq_h0 <- function(b, N, TT, y, W, p, d, Z){
  
  m <- 1 + 2 * p + max(0, ncol(Z) ) + 1
  b[m] <- 0
  
  z <- W %*% y 
  zf <- log1p(z)
  f <-  -b[1] * zf
  x2 <- NULL
  wy <- NULL
  for ( ti in (p + 1):TT ) {
    wy <- rbind( wy, cbind(z[, (ti - 1):(ti - p)], y[, (ti - 1):(ti - p)], Z, f[, ti - d]) )
    x2 <- c(x2, zf[, ti - d]^2)
  }
  wy <- cbind(1, wy)
  lambdat <- as.vector( wy %*% b )
  Dt <- 1 / lambdat
  #ss1 <- wy * ( as.vector(y[, -c(1:p)]) - lambdat ) / lambdat
  ss1 <- crossprod( wy, (as.vector(y[, -c(1:p)]) - lambdat)/lambdat2 )
  ct <- as.vector(y[, -c(1:p)]) / lambdat^2
  hh1 <- crossprod(wy * ct, wy)
  ss2a <- crossprod(x2 * Dt, as.vector(y[, -c(1:p)]) - lambdat2 ) * b[1]
  hh1[1, m] <- hh1[1, m] - ss1[m]/b[1]
  hh1[m, 1] <- hh1[1, m]
  hh1[m, m] <- hh1[m, m] - ss2a
  hh1

}

hess_nonlinpq_h0_old <- function(b, N, TT, y, W, p, d, Z){
  
  m <- 1+2*p+max(0,ncol(Z))+1
  hh <- matrix(0, nrow=m, ncol=m)
  b[m] <- 0
  
  for(t in (p+1):TT){
    f <- -b[1]*log1p(W%*%y[,t-d])
    Xt <- cbind(1, W%*%y[,(t-1):(t-p)], y[,(t-1):(t-p)], Z, f)
    lambdat <- Xt%*%b
    Ct <- diag(y[,t])%*%diag(1/as.vector(lambdat^2))
    hes <- t(Xt)%*%Ct%*%Xt
    Dt <- diag(1/as.vector(lambdat))
    ss <- t(Xt)%*%Dt%*%(y[,t]-lambdat)
    x2 <- log1p(W%*%y[,t-d])*log1p(W%*%y[,t-d])
    ss2 <- t(x2)%*%Dt%*%(y[,t]-lambdat)*b[1]
    hes[1,m] <- hes[1,m] - ss[m]/b[1]
    hes[m,1] <- hes[1,m]
    hes[m,m] <- hes[m,m] - ss2
    hh <- hh + hes
  }
  
  return(hh)
}


score_test_nonlinpq_h0 <- function(b, N, TT, y, W, p, d, Z){
  
  m <- 1 + 2 * p + max(0, ncol(Z) ) + 1
  b[m] <- 0

  z <- W %*% y 
  zf <- log1p( z )
  f <-  -b[1] * zf
  x2 <- NULL
  wy <- NULL
  for ( ti in (p + 1):TT ) {
    wy <- rbind( wy, cbind(z[, (ti - 1):(ti - p)], y[, (ti - 1):(ti - p)], Z, f[, ti - d]) )
    x2 <- c(x2, zf[, ti - d]^2)
  }
  wy <- cbind(1, wy)
  lambdat <- as.vector( wy %*% b )

  ## scor
  a <- wy * ( as.vector(y[, -c(1:p)]) - lambdat ) / lambdat
  S <- matrix( Rfast::colsums(a), m, 1 )

  ## hess  
  Dt <- 1 / lambdat
  #ss1 <- wy * ( as.vector(y[, -c(1:p)]) - lambdat ) / lambdat
  ct <- as.vector(y[, -c(1:p)]) / lambdat^2
  hh1 <- crossprod(wy * ct, wy)
  ss2a <- crossprod(x2 * Dt, as.vector(y[, -c(1:p)]) - lambdat ) * b[1]
  hh1[1, m] <- hh1[1, m] - scor[m]/b[1]
  hh1[m, 1] <- hh1[1, m]
  hh1[m, m] <- hh1[m, m] - ss2a
  H <- hh1

  ## out
  out1 <- 0
  k <- rep( 1:c(TT - p), each = N )
  b1 <- rowsum(a, k)
  for ( i in 1: c(TT - p) )  out1 <- out1 + tcrossprod(b1[i, ])
  B <- out1

  solveHmm <- solve(H[-m, -m])
  Sigma <- B[m, m] - H[m, -m] %*% solveHmm %*% B[-m, m]-
      B[m, -m] %*% solveHmm %*% H[-m, m] +
      H[m, -m] %*% solveHmm %*% B[-m, -m] %*% solveHmm %*% H[-m, m]

  as.numeric(S[m])^2/as.numeric(Sigma)
  
}


################################################################################
################################################################################
##  score_test_nonlinpq_h0() --- function for the computation of quasi score test
##  for testing linearity of the PNAR model versus the nonlinear Intercept Drift
##  (ID) alternative, with p lag; it also includes q non time-varying covariates. 
##  Inputs:
##     b = estimated parameters from the linear model, in the following order:
##         (intercept, p network effects, p ar effects, covariates)
##     N = number of nodes on the network
##     TT = temporal sample size
##     y = NxTT multivariate count time series
##     W = NxN row-normalized weighted adjacency matrix describing the network
##     p = number of lags in the model
##     d = lag parameter of nonlinear variable (should be between 1 and p).
##     Z = Nxq matrix of covariates (one for each column), where q is the number of 
##         covariates in the model. They must be non-negative
##  output:
##     LM = quasi score test statistic
################################################################################

score_test_nonlinpq_h0_old <- function(b, N, TT, y, W, p, d, Z){
  
  m <- 1+2*p+max(0,ncol(Z))+1
  
  S <- scor_nonlinpq_h0_old(b, N, TT, y, W, p, d, Z)
  H <- hess_nonlinpq_h0_old(b, N, TT, y, W, p, d, Z)
  B <- outer_nonlinpq_h0_old(b, N, TT, y, W, p, d, Z)
  Sigma <- B[m,m]-H[m,-m]%*%solve(H[-m,-m])%*%B[-m,m]-
      B[m,-m]%*%solve(H[-m,-m])%*%H[-m,m]+
      H[m,-m]%*%solve(H[-m,-m])%*%B[-m,-m]%*%solve(H[-m,-m])%*%H[-m,m]
  LM <- as.numeric(S[m])^2/as.numeric(Sigma)
  return(LM)
    
}


