
################################################################################
## Function for for testing linearity of Poisson NAR model, with p lags, PNAR(p),
## versus nonlinear Smooth Transition alternative model (STPNAR)
## it also includes q non time-varying covariates
################################################################################


################################################################################
################################################################################
##  LM_gamma_stpnarpq() --- Function to optimize the score test statistic of
##  Smooth Transition model (STNAR)  with p lags,
##  under the null assumption of linearity, with respect to unknown nuisance 
##  parameter (gamma); it also includes q non time-varying covariates.
##  Input:
##    gamma = value of non identifiable nuisance parameter
##    b = estimated parameters from the linear model, in the following order:
##        (intercept, p network effects, p ar effects, covariates)
##    N = number of nodes on the network
##    TT = temporal sample size
##    y = NxTT multivariate count time series
##    W = NxN row-normalized weighted adjacency matrix describing the network
##    p = number of lags in the model
##    d = lag parameter of nonlinear variable (should be between 1 and p)
##    Z = Nxq matrix of covariates (one for each column), where q is the number of 
##        covariates in the model. They must be non-negative
##  output:
##    LM = (-1) * value of test statistic at the specified gamma
################################################################################

LM_gamma_stpnarpq <- function(gamma, b, N, TT, y, W, p, d, Z){
  
  m <- 1+3*p+max(0,ncol(Z))
  
  ss <- as.matrix(rep(0, m))
  out <- matrix(0, nrow=m, ncol=m)
  hh <- matrix(0, nrow=m, ncol=m)
  
  b[(m-p+1):m] <- 0
  
  for(t in (p+1):TT){
    xt <- as.vector(W%*%y[,t-d])
    f <- exp(-gamma*(xt*xt))
    Xp <- W%*%y[,(t-1):(t-p)]
    Xt <- cbind(1, Xp, y[,(t-1):(t-p)], Z, Xp*f)
    lambdat <- Xt[ ,1:(m-p) ] %*% b[ 1:(m-p) ]
    Dt <- diag(1/as.vector(lambdat))
    
    s <- t(Xt)%*%Dt%*%(y[,t]-lambdat)
    ss <- ss + s
    out <- out + s%*%t(s)
    
    Ct <- diag(y[,t])%*%diag(1/as.vector(lambdat^2))
    hh <- hh + t(Xt)%*%Ct%*%Xt
  }
  
  S <- ss
  H <- hh
  B <- out
  
  Sigma <- B[ (m-p+1):m , (m-p+1):m ]-
    H[ (m-p+1):m, 1:(m-p) ] %*% solve( H[ 1:(m-p), 1:(m-p) ] ) %*% B[ 1:(m-p), (m-p+1):m ]-
    B[ (m-p+1):m, 1:(m-p) ] %*% solve( H[ 1:(m-p), 1:(m-p) ] ) %*% H[ 1:(m-p), (m-p+1):m ]+
    H[ (m-p+1):m, 1:(m-p) ] %*% solve( H[ 1:(m-p), 1:(m-p) ] ) %*% B[ 1:(m-p), 1:(m-p) ] %*% 
    solve( H[ 1:(m-p), 1:(m-p) ] ) %*% H[ 1:(m-p), (m-p+1):m ]
  
  LM <- as.numeric( t( S[ (m-p+1):m ] ) %*% solve( Sigma ) %*% S[ (m-p+1):m ] )
  LM <- -1 * LM
  return( LM )
}


################################################################################
################################################################################
##  global_optimise_LM() --- Function to optimize the LM_gamma_stpnarpq() 
##  function of score test statistic under the null assumption with respect to 
##  unknown nuisance parameter (gamma). The optimization employes the Brent
##  algorithm applied on the interval [L,U]. To be sure that the a global
##  optimum is found, the optimization is performed at (I-1) consecutive
##  equidistant sub-intervals and then the minimum/maximum over them 
##  is taken as global optimum.
##  Input:
##    f=LM_gamma_stpnarpq function to be optimize. Always report it in the call.
##    L = Lower bound of search interval for gamma.
##    U = Upper bound of search interval for gamma.
##    I = number of cutting points of interval [L,U].
##    b = estimated parameters from the linear model, in the following order:
##        (intercept, p network effects, p ar effects, covariates)
##    N = number of nodes on the network
##    TT = temporal sample size
##    y = NxTT multivariate count time series
##    W = NxN row-normalized weighted adjacency matrix describing the network
##    p = number of lags in the model
##    d = lag parameter of nonlinear variable (should be between 1 and p)
##    Z = Nxq matrix of covariates (one for each column), where q is the number of 
##        covariates in the model. They must be non-negative
##    tol = tolerance of the optimization routine. Default is 1e-14.
##  output: a list with two values
##            gamma = optimum value of gamma parameter
##            supLM = value of the objective function at the optimum
##            int = list of extremes points of subintervals
################################################################################

global_optimise_LM <- function(f=LM_gamma_stpnarpq, L, U, I, b, N, TT, y, W, p, d, Z, tol= 1e-14){
  
  gami <- vector()
  supLMi <- vector()
  
  x <- seq(from=L, to=U, length.out = I)
  
  for(i in 1:(I-1)){
    
    int <- c(x[i],x[i+1])
    opt <- optimise(f=LM_gamma_stpnarpq, int, b, N, TT, y, W, p, d, Z, tol= tol)
    
    gami[i] <- opt$minimum
    supLMi[i] <- opt$objective
    
  }
  
  supLM <- min(supLMi)
  gamma <- gami[which(supLMi==supLM)[1]]
  return(list(gamma=gamma, supLM=-1 * supLM, int=x))
  
}


################################################################################
################################################################################
##  LM_gamma_stpnarpq_j() --- Function to optimize the perturbed version of
##  the score test statistic of Smooth Transition model (STNAR)  with p lags,
##  under the null assumption of linearity, with respect to unknown nuisance 
##  parameter (gamma); it also includes q non time-varying covariates.
##  Input:
##    gamma = value of non identifiable nuisance parameter
##    b = estimated parameters from the linear model, in the following order:
##        (intercept, p network effects, p ar effects, covariates)
##    N = number of nodes on the network
##    TT = temporal sample size
##    y = NxTT multivariate count time series
##    W = NxN row-normalized weighted adjacency matrix describing the network
##    p = number of lags in the model
##    d = lag parameter of nonlinear variable (should be between 1 and p)
##    Z = Nxq matrix of covariates (one for each column), where q is the number of 
##        covariates in the model. They must be non-negative
##    msn = TTx1 vector of standard normal noises
##  output:
##    LM = (-1) * value of perturbed test statistic at the specified gamma
################################################################################

LM_gamma_stpnarpq_j <- function(gamma, b, N, TT, y, W, p, d, Z, msn){
  
  m <- 1+3*p+max(0,ncol(Z))
  
  ss <- as.matrix(rep(0, m))
  out <- matrix(0, nrow=m, ncol=m)
  hh <- matrix(0, nrow=m, ncol=m)
  
  b[(m-p+1):m] <- 0
  
  for(t in (p+1):TT){
    xt <- as.vector(W%*%y[,t-d])
    f <- exp(-gamma*(xt*xt))
    Xp <- W%*%y[,(t-1):(t-p)]
    Xt <- cbind(1, Xp, y[,(t-1):(t-p)], Z, Xp*f)
    lambdat <- Xt[ ,1:(m-p) ] %*% b[ 1:(m-p) ]
    Dt <- diag(1/as.vector(lambdat))
    
    s <- t(Xt)%*%Dt%*%(y[,t]-lambdat)
    ss <- ss + s*msn[t]
    out <- out + s%*%t(s)
    
    Ct <- diag(y[,t])%*%diag(1/as.vector(lambdat^2))
    hh <- hh + t(Xt)%*%Ct%*%Xt
  }
  
  S <- ss
  H <- hh
  B <- out
  
  Sigma <- B[ (m-p+1):m , (m-p+1):m ]-
    H[ (m-p+1):m, 1:(m-p) ] %*% solve( H[ 1:(m-p), 1:(m-p) ] ) %*% B[ 1:(m-p), (m-p+1):m ]-
    B[ (m-p+1):m, 1:(m-p) ] %*% solve( H[ 1:(m-p), 1:(m-p) ] ) %*% H[ 1:(m-p), (m-p+1):m ]+
    H[ (m-p+1):m, 1:(m-p) ] %*% solve( H[ 1:(m-p), 1:(m-p) ] ) %*% B[ 1:(m-p), 1:(m-p) ] %*% 
    solve( H[ 1:(m-p), 1:(m-p) ] ) %*% H[ 1:(m-p), (m-p+1):m ]
  
  LM <- as.numeric( t( S[ (m-p+1):m ] ) %*% solve( Sigma ) %*% S[ (m-p+1):m ] )
  LM <- -1 * LM
  return( LM )
}


###########################################################################################
###########################################################################################
##  score_test_starpq_DV() --- function for the computation of Davies bound of sup-type
##  test for testing linearity of PNAR model, with p lags, versus the nonlinear 
##  Smooth Transition model (STNAR) alternative;
##  it also includes q non time-varying covariates. 
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
##     gamma_L = lower bound of grid of values for gamma parameter.
##     gamma_U =  upper bound of grid of values for gamma parameter.
##     l = number of values in the grid for gamma parameter. Default is 100.
##  output: a list with the following objects
##        DV = Davies bound of p-values for sup test
##        supLM = value of the sup test statistic in the sample y
############################################################################################

score_test_starpq_DV <- function(b, N, TT, y, W, p, d, Z, gamma_L, gamma_U, l=100){
  
  m <- 1+3*p+max(0,ncol(Z))
  LMv <- vector()
  
  gam <- seq(from=gamma_L, to=gamma_U, length.out = l)
  
  V <- 0
  
  LMv[1] <- LM_gamma_stpnarpq(gam[1], b, N, TT, y, W, p, d, Z)
  
  for(i in 2:l){
    
    gamma <- gam[i]
    LMv[i] <- LM_gamma_stpnarpq(gamma, b, N, TT, y, W, p, d, Z)
    V <- V +  abs(sqrt(-LMv[i]) - sqrt(-LMv[i-1]))
    
  }
  
  supLM <- as.numeric(max(-LMv))
  
  DV <- 1-pchisq(supLM,p) + V*supLM^(1/2*(p-1))*exp(-supLM/2)*2^(-p/2)/gamma(p/2)
  
  return(list(DV=DV, supLM=supLM))
}


###########################################################################################
###########################################################################################
##  score_test_starpq_j() --- function for the computation of bootstrap sup-type
##  test for testing linearity of PNAR model, with p lag, versus the nonlinear 
##  Smooth Transition model (STNAR) alternative;
##  it also includes q non time-varying covariates. 
##  Inputs:
##     supLM = value of the score test at the optimum gamma
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
##     J = number of bootstrap replications
##     L = Lower bound of search interval for gamma. Default 0.001.
##     U = Upper bound of search interval for gamma. Default 1000.
##     tol = tolerance of the optimization routine. Default is 1e-14.
##  output: a list with the following objects
##        pJ = bootstrap p-value sup test
##        cpJ = adjusted version of bootstrap p-value sup test
##        gammaj = optimal values of gamma parameter for score test boostrap replications
##        supLMj = values of perturbed test statistic at the optimum point gammaj
############################################################################################

score_test_starpq_j<- function(supLM, b, N, TT, y, W, p, d, Z, J, L=0.001, U=1000, tol = 1e-14){

  m <- 1+3*p+max(0,ncol(Z))
  supLMj <- vector()
  gammaj <- vector()
  pval <- vector()

  for(j in 1:J){

    snj <- rnorm(TT)

    opttj <- optimise(f=LM_gamma_stpnarpq_j, c(L,U), b, N, TT, y, W, p, d, Z, snj, tol=tol)
    gammaj[j] <- opttj$minimum
    supLMj[j] <- -1* opttj$objective

    pval[j] = ifelse( supLMj[j] >= supLM , 1, 0 )

  }

  pJ<- sum(pval)/J
  cpJ<- (sum(pval)+1)/(J+1)

  return(list(pJ=pJ, cpJ=cpJ, supLMj=supLMj, gammaj=gammaj))

}








## ignore all these, you can delete them if you prefer
################################################################################
################################################################################
##  scor_starpq_h0() --- function for the computation of score of the nonlinear
##  Smooth Transition model (STPNAR), with p lags, under the null assumption
##  of linearity; it also includes q non time-varying covariates.
##  Inputs:
##     gamma = non identifiable nuisance parameter
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
# 
# scor_starpq_h0 <- function(gamma, b, N, TT, y, W, p, d, Z){
#   
#   m <- 1+3*p+max(0,ncol(Z))
#   ss <- as.matrix(rep(0, m))
#   b[(m-p+1):m] <- 0
#   
#   for(t in (p+1):TT){
#     xt <- as.vector(W%*%y[,t-d])
#     f <- exp(-gamma*(xt*xt))
#     Xp <- W%*%y[,(t-1):(t-p)]
#     Xt <- cbind(1, Xp, y[,(t-1):(t-p)], Z, Xp*f)
#     lambdat <- Xt[ ,1:(m-p) ] %*% b[ 1:(m-p) ]
#     Dt <- diag(1/as.vector(lambdat))
#     ss <- ss + t(Xt)%*%Dt%*%(y[,t]-lambdat)
#   }
#   
#   return(ss)
# }
# 
# 
# ################################################################################
# ################################################################################
# ##  outer_starpq_h0() --- function for the computation of information matrix
# ##  of the nonlinear Smooth Transition model (STNAR), with p lags, under the 
# ##  null assumption of linearity; it also includes q non time-varying covariates.
# ##  Inputs:
# ##     gamma = non identifiable nuisance parameter
# ##     b = estimated parameters from the linear model, in the following order:
# ##         (intercept, p network effects, p ar effects, covariates)
# ##     N = number of nodes on the network
# ##     TT = temporal sample size
# ##     y = NxTT multivariate count time series
# ##     W = NxN row-normalized weighted adjacency matrix describing the network
# ##     p = number of lags in the model
# ##     d = lag parameter of nonlinear variable (should be between 1 and p)
# ##     Z = Nxq matrix of covariates (one for each column), where q is the number of 
# ##         covariates in the model. They must be non-negative
# ##  output:
# ##     out = information matrix
# ################################################################################
# 
# outer_starpq_h0 <- function(gamma, b, N, TT, y, W, p, d, Z){
#   
#   m <- 1+3*p+max(0,ncol(Z))
#   out <- matrix(0, nrow=m, ncol=m)
#   b[(m-p+1):m] <- 0
#   
#   for(t in (p+1):TT){
#     xt <- as.vector(W%*%y[,t-d])
#     f <- exp(-gamma*(xt*xt))
#     Xp <- W%*%y[,(t-1):(t-p)]
#     Xt <- cbind(1, Xp, y[,(t-1):(t-p)], Z, Xp*f)
#     lambdat <- Xt[,1:(m-p)]%*%b[1:(m-p)]
#     Dt <- diag(1/as.vector(lambdat))
#     ss <- t(Xt)%*%Dt%*%(y[,t]-lambdat)
#     out <- out + ss%*%t(ss)
#   }
#   
#   return(out)
# }
# 
# 
# ################################################################################
# ################################################################################
# ##  hess_starpq_h0() --- function for the computation of Hessian matrix
# ##  of the nonlinear Smooth Transition model (STNAR), with p lags, under the 
# ##  null assumption of linearity; it also includes q non time-varying covariates.
# ##  Inputs:
# ##     gamma = non identifiable nuisance parameter
# ##     b = estimated parameters from the linear model, in the following order:
# ##         (intercept, p network effects, p ar effects, covariates)
# ##     N = number of nodes on the network
# ##     TT = temporal sample size
# ##     y = NxTT multivariate count time series
# ##     W = NxN row-normalized weighted adjacency matrix describing the network
# ##     p = number of lags in the model
# ##     d = lag parameter of nonlinear variable (should be between 1 and p)
# ##     Z = Nxq matrix of covariates (one for each column), where q is the number of 
# ##         covariates in the model. They must be non-negative
# ##  output:
# ##     hh = (-1)*Hessian matrix
# ################################################################################
# 
# hess_starpq_h0 <- function(gamma, b, N, TT, y, W, p, d, Z){
#   
#   m <- 1+3*p+max(0,ncol(Z))
#   hh <- matrix(0, nrow=m, ncol=m)
#   b[(m-p+1):m] <- 0
#   
#   for(t in (p+1):TT){
#     xt <- as.vector(W%*%y[,t-d])
#     f <- exp(-gamma*(xt*xt))
#     Xp <- W%*%y[,(t-1):(t-p)]
#     Xt <- cbind(1, Xp, y[,(t-1):(t-p)], Z, Xp*f)
#     lambdat <- Xt[,1:(m-p)]%*%b[1:(m-p)]
#     Ct <- diag(y[,t])%*%diag(1/as.vector(lambdat^2))
#     hh <- hh + t(Xt)%*%Ct%*%Xt
#   }
#   
#   return(hh)
# }
# 
# 
# ################################################################################
# ################################################################################
# ##  scor_starpq_j() --- function for the computation of the perturbed version
# ##  of the score of the nonlinear Smooth Transition model (STNAR), with p lags,
# ##  under the null assumption of linearity; 
# ##  it also includes q non time-varying covariates.
# ##  Inputs:
# ##     gamma = non identifiable nuisance parameter
# ##     b = estimated parameters from the linear model, in the following order:
# ##         (intercept, p network effects, p ar effects, covariates)
# ##     N = number of nodes on the network
# ##     TT = temporal sample size
# ##     y = NxTT multivariate count time series
# ##     W = NxN row-normalized weighted adjacency matrix describing the network
# ##     p = number of lags in the model
# ##     d = lag parameter of nonlinear variable (should be between 1 and p)
# ##     Z = Nxq matrix of covariates (one for each column), where q is the number of 
# ##         covariates in the model. They must be non-negative
# ##     msn = TTx1 vector of standard normal noises
# ##  output:
# ##     ss = vector of perturbed quasi score
# ################################################################################
# 
# scor_starpq_j <- function(gamma, b, N, TT, y, W, p, d, Z, msn){
#   
#   m <- 1+3*p+max(0,ncol(Z))
#   ss <- as.matrix(rep(0, m))
#   b[(m-p+1):m] <- 0
#   
#   for(t in (p+1):TT){
#     
#     xt <- as.vector(W%*%y[,t-d])
#     f <- exp(-gamma*(xt*xt))
#     Xp <- W%*%y[,(t-1):(t-p)]
#     Xt <- cbind(1, Xp, y[,(t-1):(t-p)], Z, Xp*f)
#     lambdat <- Xt[,1:(m-p)]%*%b[1:(m-p)]
#     Dt <- diag(1/as.vector(lambdat))
#     sh <- t(Xt)%*%Dt%*%(y[,t]-lambdat)*msn[t]
#     ss <- ss + sh
#   }
#   
#   return(ss)
# }
# 

#
# LM_gamma_stpnarpq2 <- function(gamma, b, N, TT, y, W, p, d, Z){
#   
#   m <- 1+3*p+max(0,ncol(Z))
#   S <- scor_starpq_h0(gamma, b, N, TT, y, W, p, d, Z)
#   H <- hess_starpq_h0(gamma, b, N, TT, y, W, p, d, Z)
#   B <- outer_starpq_h0(gamma, b, N, TT, y, W, p, d, Z)
#   
#   Sigma <- B[ (m-p+1):m , (m-p+1):m ]-
#     H[ (m-p+1):m, 1:(m-p) ] %*% solve( H[ 1:(m-p), 1:(m-p) ] ) %*% B[ 1:(m-p), (m-p+1):m ]-
#     B[ (m-p+1):m, 1:(m-p) ] %*% solve( H[ 1:(m-p), 1:(m-p) ] ) %*% H[ 1:(m-p), (m-p+1):m ]+
#     H[ (m-p+1):m, 1:(m-p) ] %*% solve( H[ 1:(m-p), 1:(m-p) ] ) %*% B[ 1:(m-p), 1:(m-p) ] %*% 
#     solve( H[ 1:(m-p), 1:(m-p) ] ) %*% H[ 1:(m-p), (m-p+1):m ]
#   
#   LM <- as.numeric( t( S[ (m-p+1):m ] ) %*% solve( Sigma ) %*% S[ (m-p+1):m ] )
#   LM <- -1 * LM
#   return( LM )
# }
#
#V <- integrate(dLM_gamma_spnarpq, lower=0, upper=Inf, b=b, N=N, TT=TT, 
#               y=y, W=W, p=p, d=d, Z=Z, h=1e-8, rel.tol = 1e-8)$value
#
################################################################################
################################################################################
##  dLM_gamma_spnarpq() --- Function for the computation of the absolute
##  difference quotient of (square-rooted) the score test statistic,
##  function LM_gamma_stpnarpq(), with respect to unknown nuisance parameter
##  (gamma), for the Smooth Transition model (STNAR), with p lags,
##  under the null assumption of linearity; 
##  it also includes q non time-varying covariates.
##  Input:
##    gamma = value of non identifiable nuisance parameter
##    b = estimated parameters from the linear model, in the following order:
##        (intercept, p network effects, p ar effects, covariates)
##    N = number of nodes on the network
##    TT = temporal sample size
##    y = NxTT multivariate count time series
##    W = NxN row-normalized weighted adjacency matrix describing the network
##    p = number of lags in the model
##    d = lag parameter of nonlinear variable (should be between 1 and p)
##    Z = Nxq matrix of covariates (one for each column), where q is the number of 
##        covariates in the model. They must be non-negative
##    h = difference variable. A small positive value. Default is 1e-8.
##  output:
##    dLM = value of absolute difference quotient of the (square-rooted)
##          test statistic at the specified gamma
################################################################################
#
# dLM_gamma_spnarpq2 <- function(gamma, b, N, TT, y, W, p, d, Z, h=1e-8){
# 
#   LMh <- LM_gamma_stpnarpq(gamma+h, b, N, TT, y, W, p, d, Z)
#   LM <- LM_gamma_stpnarpq(gamma, b, N, TT, y, W, p, d, Z)
#   dLM <- abs((sqrt(-1*LMh)-sqrt(-1*LM))/h)
#   return(dLM)
# 
# }
# 
# dLM_gamma_spnarpq <- function(gamma, b, N, TT, y, W, p, d, Z, h=1e-8){
# 
# 
#   m <- 1+3*p+max(0,ncol(Z))
# 
#   ss <- as.matrix(rep(0, m))
#   out <- matrix(0, nrow=m, ncol=m)
#   hh <- matrix(0, nrow=m, ncol=m)
# 
#   ssh <- as.matrix(rep(0, m))
#   outh <- matrix(0, nrow=m, ncol=m)
#   hhh <- matrix(0, nrow=m, ncol=m)
# 
#   b[(m-p+1):m] <- 0
# 
#   for(t in (p+1):TT){
# 
#     xt <- as.vector( W %*% y[,t-d] )
# 
#     f <- exp( -gamma * ( xt * xt ) )
#     #f <- exp( -rep(gamma, N) * ( xt * xt ) )
# 
#     Xp <- W %*% y[,(t-1):(t-p)]
#     Xt <- cbind(1, Xp, y[,(t-1):(t-p)], Z, Xp*f)
#     lambdat <- Xt[,1:(m-p)]%*%b[1:(m-p)]
#     Dt <- diag(1/as.vector(lambdat))
# 
#     s <- t(Xt)%*%Dt%*%(y[,t]-lambdat)
#     ss <- ss + s
#     out <- out + s%*%t(s)
# 
#     Ct <- diag(y[,t])%*%diag(1/as.vector(lambdat^2))
#     hh <- hh + t(Xt)%*%Ct%*%Xt
# 
#     fh <- exp( -(gamma+h) * (xt * xt) )
#     #fh <- exp( -rep(gamma+h, N) * (xt*xt) )
# 
# 
#     Xth <- cbind(1, Xp, y[,(t-1):(t-p)], Z, Xp*fh)
# 
#     sh <- t(Xth)%*%Dt%*%(y[,t]-lambdat)
#     ssh <- ssh + sh
#     outh <- outh + sh%*%t(sh)
#     hhh <- hhh + t(Xth)%*%Ct%*%Xth
#   }
# 
#   S <- ss
#   H <- hh
#   B <- out
# 
#   Sigma <- B[ (m-p+1):m , (m-p+1):m ]-
#     H[ (m-p+1):m, 1:(m-p) ] %*% solve( H[ 1:(m-p), 1:(m-p) ] ) %*% B[ 1:(m-p), (m-p+1):m ]-
#     B[ (m-p+1):m, 1:(m-p) ] %*% solve( H[ 1:(m-p), 1:(m-p) ] ) %*% H[ 1:(m-p), (m-p+1):m ]+
#     H[ (m-p+1):m, 1:(m-p) ] %*% solve( H[ 1:(m-p), 1:(m-p) ] ) %*% B[ 1:(m-p), 1:(m-p) ] %*%
#     solve( H[ 1:(m-p), 1:(m-p) ] ) %*% H[ 1:(m-p), (m-p+1):m ]
# 
#   LM <- as.numeric( t( S[ (m-p+1):m ] ) %*% solve( Sigma ) %*% S[ (m-p+1):m ] )
# 
#   S <- ssh
#   H <- hhh
#   B <- outh
# 
#   Sigma <- B[ (m-p+1):m , (m-p+1):m ]-
#     H[ (m-p+1):m, 1:(m-p) ] %*% solve( H[ 1:(m-p), 1:(m-p) ] ) %*% B[ 1:(m-p), (m-p+1):m ]-
#     B[ (m-p+1):m, 1:(m-p) ] %*% solve( H[ 1:(m-p), 1:(m-p) ] ) %*% H[ 1:(m-p), (m-p+1):m ]+
#     H[ (m-p+1):m, 1:(m-p) ] %*% solve( H[ 1:(m-p), 1:(m-p) ] ) %*% B[ 1:(m-p), 1:(m-p) ] %*%
#     solve( H[ 1:(m-p), 1:(m-p) ] ) %*% H[ 1:(m-p), (m-p+1):m ]
# 
#   LMh <- as.numeric( t( S[ (m-p+1):m ] ) %*% solve( Sigma ) %*% S[ (m-p+1):m ] )
# 
#   dLM <- abs((sqrt(LMh)-sqrt(LM))/h)
#   return(-1 * dLM)
# 
# }
#
###########################################################################################
###########################################################################################
##  score_test_starpq_j() --- function for the computation of bootstrap sup-type
##  test for testing linearity of PNAR model, with p lag, versus the nonlinear 
##  Smooth Transition model (STNAR) alternative;
##  it also includes q non time-varying covariates. 
##  Inputs:
##     gamma = optimal value of nuisance parameter 
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
##     J = number of bootstrap replications
##  output: a list with the following objects
##        supLM = value of the sup test statistic in the sample y
##        DV = Davies bound of p-values for sup test
##        pJ = bootstrap p-value sup test
##        cpJ = adjusted version of bootstrap p-value sup test
##        LM = values of test statistic for different values in the grid
##        supLMj = values of perturbed sup test statistic
############################################################################################
# 
# score_test_starpq_j<- function(gamma, b, N, TT, y, W, p, d, Z, J){
#   
#   m <- 1+3*p+max(0,ncol(Z))
#   supLMj <- vector()
#   pval <- vector()
#   
#   S <- scor_starpq_h0(gamma, b, N, TT, y, W, p, d, Z)
#   H <- hess_starpq_h0(gamma, b, N, TT, y, W, p, d, Z)
#   B <- outer_starpq_h0(gamma, b, N, TT, y, W, p, d, Z)
#   
#   Sigma <- B[ (m-p+1):m , (m-p+1):m ] -
#     H[ (m-p+1):m, 1:(m-p) ] %*% solve( H[ 1:(m-p), 1:(m-p) ] ) %*% B[ 1:(m-p), (m-p+1):m ] -
#     B[ (m-p+1):m, 1:(m-p) ] %*% solve( H[ 1:(m-p), 1:(m-p) ] ) %*% H[ 1:(m-p), (m-p+1):m ] +
#     H[ (m-p+1):m, 1:(m-p) ] %*% solve( H[ 1:(m-p), 1:(m-p) ] ) %*% B[ 1:(m-p), 1:(m-p) ] %*% 
#     solve( H[ 1:(m-p), 1:(m-p) ] ) %*% H[ 1:(m-p), (m-p+1):m ]
#   
#   supLM <- as.numeric( t( S[ (m-p+1):m ] ) %*% solve( Sigma ) %*% S[ (m-p+1):m ] )
#   
#   for(j in 1:J){
#     
#     snj <- rnorm(TT)
#     Sj <- scor_starpq_j(gamma, b, N, TT, y, W, p, d, Z, snj)
#     
#     supLMj[j] <- as.numeric( t( Sj[ (m-p+1):m ] ) %*% solve( Sigma ) %*% Sj[ (m-p+1):m ] )
#     
#     pval[j]=ifelse(supLMj[j]>=supLM,1,0)
#     
#   }
#   
#   pJ<- sum(pval)/J
#   cpJ<- (sum(pval)+1)/(J+1)
#   
#   return(list(supLM=supLM, pJ=pJ, cpJ=cpJ))
# }
#
#
################################################################################
################################################################################
##  optim_gamma_stpnarpq() --- Function to optimize the LM_gamma_stpnarpq() of
##  score test statistic under the null assumption with respect to unknown
##  nuisance parameter (gamma)
##  Input:
##    b = estimated parameters from the linear model, in the following order:
##        (intercept, p network effects, p ar effects, covariates)
##    N = number of nodes on the network
##    TT = temporal sample size
##    y = NxTT multivariate count time series
##    W = NxN row-normalized weighted adjacency matrix describing the network
##    p = number of lags in the model
##    d = lag parameter of nonlinear variable (should be between 1 and p)
##    Z = Nxq matrix of covariates (one for each column), where q is the number of 
##        covariates in the model. They must be non-negative
##    method = c( "Brent", "BFGS" ) the optimization method. Default is "Brent".
##             "BFGS" uses finite different approximation gradient.
##    x0 = Starting value of optimization. Only for BFGS. Default 1.
##    L = Lower bound of search interval for gamma. Only for Brent. Default 1e-4.
##    U = Upper bound of search interval for gamma. Only for Brent. Default 10.
##  output: a list with two values
##            minimum = optimum value of gamma parameter
##            objective = (-1)*value of the objective function at the optimum
################################################################################
# 
# optim_gamma_stpnarpq <- function(b, N, TT, y, W, p, d, Z, method="Brent", x0=1, L=1e-4, U=10){
#   
#   ## please pay attention to the interval of the Brent here.
#   ## Indeed, if we set the interval to big or too small
#   ## it fails to find the optimum. So I added BFGS alternative.
#   
#   if(method=="Brent"){
#     s_gam <- optimise(f=LM_gamma_stpnarpq, c(L,U), b, N, TT, y, W, p, d, Z, tol = 1e-8)
#     return(s_gam)
#   } else{
#     
#     s_gam <- optim(par=x0,fn=LM_gamma_stpnarpq, gr=NULL, b, N, TT, y, W, p, d, Z, method="BFGS")
#     return(list(minimum=s_gam$par, objective=s_gam$value))
#   }
#   
# }
# Sj <- scor_starpq_j(gamma, b, N, TT, y, W, p, d, Z, snj)
# 
# supLMj[j] <- as.numeric( t( Sj[ (m-p+1):m ] ) %*% solve( Sigma ) %*% Sj[ (m-p+1):m ] )
#
#
# global_optimise_LM_j <- function(f=LM_gamma_stpnarpq_j, L, U, I, b, N, TT, y, W, p, d, Z, msn, tol= 1e-14){
#   
#   gami <- vector()
#   supLMi <- vector()
#   
#   x <- seq(from=L, to=U, length.out = I)
#   
#   for(i in 1:(I-1)){
#     
#     int <- c(x[i],x[i+1])
#     opt <- optimise(f=LM_gamma_stpnarpq_j, int, b, N, TT, y, W, p, d, Z, msn, tol= tol)
#     
#     gami[i] <- opt$minimum
#     supLMi[i] <- opt$objective
#     
#   }
#   
#   supLM <- min(supLMi)
#   gamma <- gami[which(supLMi==supLM)[1]]
#   return(list(gamma=gamma, supLM=-1 * supLM, int=x))
#   
# }
#
#
# score_test_starpq_j3<- function(b, N, TT, y, W, p, d, Z, J, L, U, I, tol = 1e-14){
#   
#   m <- 1+3*p+max(0,ncol(Z))
#   supLMj <- vector()
#   gammaj <- vector()
#   pval <- vector()
#   
#   optt <- global_optimise_LM(f=LM_gamma_stpnarpq, L, U, I, b, N, TT, y, W, p, d, Z, tol= tol)
#   gamma <- optt$gamma
#   supLM <- optt$supLM
#   
#   #optt <- optimise(f=LM_gamma_stpnarpq, c(1e-4,10), b, N, TT, y, W, p, d, Z, tol=tol)
#   #gamma <- optt$minimum
#   #supLM <- -1* optt$objective
#   
#   for(j in 1:J){
#     
#     snj <- rnorm(TT)
#     
#     opttj <- global_optimise_LM_j(f=LM_gamma_stpnarpq_j, L, U, I, b, N, TT, y, W, p, d, Z, snj, tol= tol)
#     gammaj[j] <- opttj$gamma
#     supLMj[j] <- opttj$supLM
#     
#     #opttj <- optimise(f=LM_gamma_stpnarpq_j, c(L,U), b, N, TT, y, W, p, d, Z, snj, tol=tol)
#     #gammaj[j] <- opttj$minimum
#     #supLMj[j] <- -1* opttj$objective
#     
#     pval[j] = ifelse( supLMj[j] >= supLM , 1, 0 )
#     
#   }
#   
#   pJ<- sum(pval)/J
#   cpJ<- (sum(pval)+1)/(J+1)
#   
#   return(list(pJ=pJ, cpJ=cpJ, supLM=supLM, gamma=gamma, supLMj=supLMj, gammaj=gammaj))
# }
#
#
# score_test_starpq_j2<- function(b, N, TT, y, W, p, d, Z, J, L=1e-4, U=10, tol = 1e-14){
#   
#   m <- 1+3*p+max(0,ncol(Z))
#   supLMj <- vector()
#   pval <- vector()
#   
#   #convj <- vector()
#   
#   optt <- optimise(f=LM_gamma_stpnarpq, c(L,U), b, N, TT, y, W, p, d, Z, tol=tol)
#   gamma <- optt$minimum
#   supLM <- -1* optt$objective
#   
#   for(j in 1:J){
#     
#     snj <- rnorm(TT)
#     supLMj[j] <- -1 * LM_gamma_stpnarpq_j(gamma, b, N, TT, y, W, p, d, Z, snj)
#     pval[j] = ifelse( supLMj[j] >= supLM , 1, 0 )
#     
#   }
#   
#   pJ<- sum(pval)/J
#   cpJ<- (sum(pval)+1)/(J+1)
#   
#   return(list(pJ=pJ, cpJ=cpJ, supLM=supLM, gamma=gamma)) #, supLMj=supLMj))
# }
