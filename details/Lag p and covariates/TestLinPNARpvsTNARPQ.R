
################################################################################
## Function for for testing linearity of Poisson NAR model, with p lags, PNAR(p),
## versus nonlinear Threshold alternative model (TPNAR)
## it also includes q non time-varying covariates
################################################################################


################################################################################
################################################################################
##  LM_gamma_tpnarpq() --- Function to optimize the score test statistic of
##  Threshold model (TNAR)  with p lags, under the null assumption of linearity,
##  with respect to unknown nuisance parameter (gamma); 
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
##  output:
##    LM = (-1) * value of test statistic at the specified gamma
################################################################################

LM_gamma_tpnarpq <- function(gamma, b, N, TT, y, W, p, d, Z){
  
  pp <- 1+2*p
  m <- pp*2+max(0,ncol(Z))
  
  ss <- as.matrix(rep(0, m))
  out <- matrix(0, nrow=m, ncol=m)
  hh <- matrix(0, nrow=m, ncol=m)
  
  b[(m-pp+1):m] <- 0
  
  for(t in (p+1):TT){
    f <- as.vector(0+(W%*%y[,t-d]<=gamma))
    Xp <- W%*%y[,(t-1):(t-p)]
    Xt <- cbind(1, Xp, y[,(t-1):(t-p)], Z, 1*f, Xp*f, y[,(t-1):(t-p)]*f)
    lambdat <- Xt[ ,1:(m-pp) ] %*% b[ 1:(m-pp) ]
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
  
  Sigma <- B[ (m-pp+1):m , (m-pp+1):m ]-
    H[ (m-pp+1):m, 1:(m-pp) ] %*% solve( H[ 1:(m-pp), 1:(m-pp) ] ) %*% B[ 1:(m-pp), (m-pp+1):m ]-
    B[ (m-pp+1):m, 1:(m-pp) ] %*% solve( H[ 1:(m-pp), 1:(m-pp) ] ) %*% H[ 1:(m-pp), (m-pp+1):m ]+
    H[ (m-pp+1):m, 1:(m-pp) ] %*% solve( H[ 1:(m-pp), 1:(m-pp) ] ) %*% B[ 1:(m-pp), 1:(m-pp) ] %*% 
    solve( H[ 1:(m-pp), 1:(m-pp) ] ) %*% H[ 1:(m-pp), (m-pp+1):m ]
  
  LM <- as.numeric( t( S[ (m-pp+1):m ] ) %*% solve( Sigma ) %*% S[ (m-pp+1):m ] )
  LM <- -1 * LM
  return( LM )
}


################################################################################
################################################################################
##  global_optimise_LM() --- Function to optimize the LM_gamma_tpnarpq() 
##  function of score test statistic under the null assumption with respect to 
##  unknown nuisance parameter (gamma). The optimization employes the Brent
##  algorithm applied on the interval [L,U]. To be sure that the a global
##  optimum is found, the optimization is performed at (I-1) consecutive
##  equidistant sub-intervals and then the minimum/maximum over them 
##  is taken as global optimum.
##  Input:
##    f=LM_gamma_tpnarpq function to be optimize. Always report it in the call.
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

global_optimise_LM <- function(f=LM_gamma_tpnarpq, L, U, I, b, N, TT, y, W, p, d, Z, tol= 1e-14){
  
  gami <- vector()
  supLMi <- vector()
  
  x <- seq(from=L, to=U, length.out = I)
  
  for(i in 1:(I-1)){
    
    int <- c(x[i],x[i+1])
    opt <- optimise(f=LM_gamma_tpnarpq, int, b, N, TT, y, W, p, d, Z, tol= tol)
    
    gami[i] <- opt$minimum
    supLMi[i] <- opt$objective
    
  }
  
  supLM <- min(supLMi)
  gamma <- gami[which(supLMi==supLM)[1]]
  return(list(gamma=gamma, supLM=-1 * supLM, int=x))
  
}


################################################################################
################################################################################
##  LM_gamma_tpnarpq_j() --- Function to optimize the perturbed version of 
##  the score test statistic of Threshold model (TNAR)  with p lags, 
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
##    LM = (-1) * value of test statistic at the specified gamma
################################################################################

LM_gamma_tpnarpq_j <- function(gamma, b, N, TT, y, W, p, d, Z, msn){
  
  pp <- 1+2*p
  m <- pp*2+max(0,ncol(Z))
  
  ss <- as.matrix(rep(0, m))
  out <- matrix(0, nrow=m, ncol=m)
  hh <- matrix(0, nrow=m, ncol=m)
  
  b[(m-pp+1):m] <- 0
  
  for(t in (p+1):TT){
    f <- as.vector(0+(W%*%y[,t-d]<=gamma))
    Xp <- W%*%y[,(t-1):(t-p)]
    Xt <- cbind(1, Xp, y[,(t-1):(t-p)], Z, 1*f, Xp*f, y[,(t-1):(t-p)]*f)
    lambdat <- Xt[ ,1:(m-pp) ] %*% b[ 1:(m-pp) ]
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
  
  Sigma <- B[ (m-pp+1):m , (m-pp+1):m ]-
    H[ (m-pp+1):m, 1:(m-pp) ] %*% solve( H[ 1:(m-pp), 1:(m-pp) ] ) %*% B[ 1:(m-pp), (m-pp+1):m ]-
    B[ (m-pp+1):m, 1:(m-pp) ] %*% solve( H[ 1:(m-pp), 1:(m-pp) ] ) %*% H[ 1:(m-pp), (m-pp+1):m ]+
    H[ (m-pp+1):m, 1:(m-pp) ] %*% solve( H[ 1:(m-pp), 1:(m-pp) ] ) %*% B[ 1:(m-pp), 1:(m-pp) ] %*% 
    solve( H[ 1:(m-pp), 1:(m-pp) ] ) %*% H[ 1:(m-pp), (m-pp+1):m ]
  
  LM <- as.numeric( t( S[ (m-pp+1):m ] ) %*% solve( Sigma ) %*% S[ (m-pp+1):m ] )
  LM <- -1 * LM
  return( LM )
}


###########################################################################################
###########################################################################################
##  score_test_tarpq_j() --- function for the computation of bootstrap sup-type
##  test for testing linearity of PNAR model, with p lag, versus the nonlinear 
##  Threshold model (TNAR) alternative; it also includes q non time-varying covariates. 
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
##     L = Lower bound of search interval for gamma.
##     U = Upper bound of search interval for gamma.
##     tol = tolerance of the optimization routine. Default is 1e-14.
##  output: a list with the following objects
##        pJ = bootstrap p-value sup test
##        cpJ = adjusted version of bootstrap p-value sup test
##        gammaj = optimal values of gamma parameter for score test boostrap replications
##        supLMj = values of perturbed test statistic at the optimum point gammaj
############################################################################################

score_test_tarpq_j<- function(supLM, b, N, TT, y, W, p, d, Z, J, L, U, tol = 1e-14){
  
  pp <- 1+2*p
  m <- pp*2+max(0,ncol(Z))
  
  supLMj <- vector()
  gammaj <- vector()
  pval <- vector()
  
  for(j in 1:J){
    
    snj <- rnorm(TT)
    
    opttj <- optimise(f=LM_gamma_tpnarpq_j, c(L,U), b, N, TT, y, W, p, d, Z, msn=snj, tol=tol)
    gammaj[j] <- opttj$minimum
    supLMj[j] <- -1* opttj$objective
    
    #supLMj[j]<- -1 * LM_gamma_tpnarpq_j(gamma, b, N, TT, y, W, p, d, Z, msn=snj)
    
    pval[j] = ifelse( supLMj[j] >= supLM , 1, 0 )
    
  }
  
  pJ<- sum(pval)/J
  cpJ<- (sum(pval)+1)/(J+1)
  
  return(list(pJ=pJ, cpJ=cpJ, supLMj=supLMj, gammaj=gammaj))
  
}

