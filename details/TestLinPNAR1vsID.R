
################################################################################
## Function for for testing linearity of Poisson NAR model, PNAR, with 1 lag
## versus nonlinear intercept drift alternative model
################################################################################


################################################################################
################################################################################
##  scor_nonlin_h0() --- function for the computation of score of the nonlinear
##  intercept drift (ID)PNAR model, with 1 lag, under the null assumption
##  of linearity. 
##  Inputs:
##     c = estimated parameters from the linear model
##     N = number of nodes on the network
##     TT = temporal sample size
##     y = NxTT multivariate count time series
##     W = NxN row-normalized weighted adjacency matrix describing the network
##  output:
##     ss = vector of quasi score
################################################################################

scor_nonlin_h0 <- function(c, N, TT, y, W){
  
  unon <- as.matrix(rep(1, N))
  ss <- as.matrix(rep(0, 4))
  c[4] <- 0
  
  for(t in 2:TT){
    Xt <- cbind(unon, W%*%as.matrix(y[,t-1]), as.matrix(y[,t-1]), -c[1]*log(unon+W%*%as.matrix(y[,t-1])))
    lambdat <- Xt%*%c
    Dt <- diag(1/as.vector(lambdat))
    ss <- ss + t(Xt)%*%Dt%*%(y[,t]-lambdat)
  }
  
  return(ss)
}


################################################################################
################################################################################
##  outer_nonlin_h0() --- function for the computation of information matrix
##  of the nonlinear intercept drift (ID)PNAR model, with 1 lag, 
##  under the null assumption of linearity. 
##  Inputs:
##     c = estimated parameters from the linear model
##     N = number of nodes on the network
##     TT = temporal sample size
##     y = NxTT multivariate count time series
##     W = NxN row-normalized weighted adjacency matrix describing the network
##  output:
##     out = information matrix
################################################################################

outer_nonlin_h0 <- function(c, N, TT, y, W){
  
  unon <- as.matrix(rep(1, N))
  out <- matrix(0, nrow=4, ncol=4)
  c[4] <- 0
  
  for(t in 2:TT){
    Xt <- cbind(unon, W%*%as.matrix(y[,t-1]), as.matrix(y[,t-1]), -c[1]*log(unon+W%*%as.matrix(y[,t-1])))
    lambdat <- Xt%*%c
    Dt <- diag(1/as.vector(lambdat))
    ss <- t(Xt)%*%Dt%*%(y[,t]-lambdat)
    out <- out + ss%*%t(ss)
  }
  
  return(out)
}


################################################################################
################################################################################
##  hess_nonlin_h0() --- function for the computation of Hessian matrix
##  of the nonlinear intercept drift (ID)PNAR model, with 1 lag, 
##  under the null assumption of linearity. 
##  Inputs:
##     c = estimated parameters from the linear model
##     N = number of nodes on the network
##     TT = temporal sample size
##     y = NxTT multivariate count time series
##     W = NxN row-normalized weighted adjacency matrix describing the network
##  output:
##     hh = (-1)*Hessian matrix
################################################################################

hess_nonlin_h0 <- function(c, N, TT, y, W){
  
  unon <- as.matrix(rep(1, N))
  hh <- matrix(0, nrow=4, ncol=4)
  c[4] <- 0
  
  for(t in 2:TT){
    Xt <- cbind(unon, W%*%as.matrix(y[,t-1]), as.matrix(y[,t-1]), -c[1]*log(unon+W%*%as.matrix(y[,t-1])))
    lambdat <- Xt%*%c
    Ct <- diag(y[,t])%*%diag(1/as.vector(lambdat^2))
    hes <- t(Xt)%*%Ct%*%Xt
    Dt <- diag(1/as.vector(lambdat))
    ss <- t(Xt)%*%Dt%*%(y[,t]-lambdat)
    x2 <- log(unon+W%*%y[,t-1])*log(unon+W%*%y[,t-1])
    ss2 <- t(x2)%*%Dt%*%(y[,t]-lambdat)*c[1]
    hes[1,4] <- hes[1,4] - ss[4]/c[1]
    hes[4,1] <- hes[1,4]
    hes[4,4] <- hes[4,4] - ss2
    hh <- hh + hes
  }
  
  return(hh)
}


################################################################################
################################################################################
##  score_test_nonlin_h0() --- function for the computation of quasi score test
##  for testing linearity of the PNAR model versus the nonlinear intercept drift
##  alternative, with 1 lag, 
##  Inputs:
##     c = estimated parameters from the linear model
##     N = number of nodes on the network
##     TT = temporal sample size
##     y = NxTT multivariate count time series
##     W = NxN row-normalized weighted adjacency matrix describing the network
##  output:
##     LM = quasi score test statistic
################################################################################

score_test_nonlin_h0 <- function(c, N, TT, y, W){
  
  S <- scor_nonlin_h0(c, N, TT, y, W)
  H <- hess_nonlin_h0(c, N, TT, y, W)
  B <- outer_nonlin_h0(c, N, TT, y, W)
  Sigma <- B[4,4]-H[4,1:3]%*%solve(H[1:3,1:3])%*%B[1:3,4]-
    B[4,1:3]%*%solve(H[1:3,1:3])%*%H[1:3,4]+
    H[4,1:3]%*%solve(H[1:3,1:3])%*%B[1:3,1:3]%*%solve(H[1:3,1:3])%*%H[1:3,4]
  LM <- as.numeric(S[4])^2/as.numeric(Sigma)
  
  return(LM)
}
