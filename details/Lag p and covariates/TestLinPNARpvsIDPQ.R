
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

scor_nonlinpq_h0 <- function(b, N, TT, y, W, p, d, Z){
  
  m <- 1+2*p+max(0,ncol(Z))+1
  ss <- as.matrix(rep(0, m))
  b[m] <- 0
  
  for(t in (p+1):TT){
    f <- -b[1]*log1p(W%*%y[,t-d])
    Xt <- cbind(1, W%*%y[,(t-1):(t-p)], y[,(t-1):(t-p)], Z, f)
    lambdat <- Xt%*%b
    Dt <- diag(1/as.vector(lambdat))
    ss <- ss + t(Xt)%*%Dt%*%(y[,t]-lambdat)
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

score_test_nonlinpq_h0 <- function(b, N, TT, y, W, p, d, Z){
  
  m <- 1+2*p+max(0,ncol(Z))+1
  
  S <- scor_nonlinpq_h0(b, N, TT, y, W, p, d, Z)
  H <- hess_nonlinpq_h0(b, N, TT, y, W, p, d, Z)
  B <- outer_nonlinpq_h0(b, N, TT, y, W, p, d, Z)
  Sigma <- B[m,m]-H[m,-m]%*%solve(H[-m,-m])%*%B[-m,m]-
      B[m,-m]%*%solve(H[-m,-m])%*%H[-m,m]+
      H[m,-m]%*%solve(H[-m,-m])%*%B[-m,-m]%*%solve(H[-m,-m])%*%H[-m,m]
  LM <- as.numeric(S[m])^2/as.numeric(Sigma)
  return(LM)
    
}




##ignore this
#
# if(is.null(d)==TRUE){
#   LM <- vector()
#   for (d in 1:p) {
#     S <- scor_nonlinpq_h0(b, N, TT, y, W, p, d, Z)
#     H <- hess_nonlinpq_h0(b, N, TT, y, W, p, d, Z)
#     B <- outer_nonlinpq_h0(b, N, TT, y, W, p, d, Z)
#     Sigma <- B[m,m]-H[m,-m]%*%solve(H[-m,-m])%*%B[-m,m]-
#              B[m,-m]%*%solve(H[-m,-m])%*%H[-m,m]+
#              H[m,-m]%*%solve(H[-m,-m])%*%B[-m,-m]%*%solve(H[-m,-m])%*%H[-m,m]
#     LM[d] <- as.numeric(S[m])^2/as.numeric(Sigma)
#     return(max(LM))
#   }
# } else{
