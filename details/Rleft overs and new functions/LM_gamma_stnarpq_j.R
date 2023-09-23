################################################################################
################################################################################
##  LM_gama_stpnarpq() --- Function to optimize the score test statistic of
##  Smooth Transition model (STNAR)  with p lags,
##  under the null assumption of linearity, with respect to unknown nuisance
##  parameter (gama); it also includes q non time-varying covariates.
##  Input:
##    gama = value of non identifiable nuisance parameter
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
##    LM = (-1) * value of test statistic at the specified gama
################################################################################

LM_gama_stpnarpq <- function(gama, b, N, TT, y, W, p, d, Z) {

  dimz <- ( !is.null(Z) ) * NCOL(Z)
  m <- 1 + 3 * p + max(0, dimz)

  b[ (m - p + 1):m ] <- 0

  z <- W %*% y
  f <- exp(-gama * z^2)

  wy <- NULL
  for ( ti in (p + 1):TT )
    wy <- rbind( wy, cbind(z[, (ti - 1):(ti - p)], y[, (ti - 1):(ti - p)], Z, z[, (ti - 1):(ti - p)] * f[, ti - d] ) )
  wy <- cbind(1, wy)
  lambdat <- as.vector( wy[, 1:c(m - p)] %*% b[1:c(m - p)] )
  Dt <- 1/lambdat

  ## scor
  a <- wy * ( ( as.vector(y[, -c(1:p)]) - lambdat ) / lambdat )
  ss1 <- matrix( Rfast::colsums(a), m, 1 )

  ## hess
  ct <- as.vector(y[, -c(1:p)]) / lambdat^2
  hh1 <- crossprod(wy * ct, wy)

  ## out
  out1 <- 0
  k <- rep( 1:c(TT - p), each = N )
  b1 <- rowsum(a, k)
  for ( i in 1: c(TT - p) )  out1 <- out1 + tcrossprod(b1[i, ])

  S <- ss1
  H <- hh1
  B <- out1
  solveHmp <- solve( H[ 1:(m - p), 1:(m - p) ] )

  Sigma <- B[ (m-p+1):m , (m-p+1):m ] -
    H[ (m - p + 1):m, 1:(m - p) ] %*% solveHmp %*% B[ 1:(m - p), (m - p + 1):m ] -
    B[ (m - p + 1):m, 1:(m - p) ] %*% solveHmp %*% H[ 1:(m - p), (m - p + 1):m ] +
    H[ (m - p + 1):m, 1:(m - p) ] %*% solveHmp %*% B[ 1:(m - p), 1:(m - p) ] %*%
    solveHmp %*% H[ 1:(m - p), (m - p + 1):m ]

  - as.numeric( crossprod( S[ (m - p + 1):m ], solve( Sigma, S[ (m - p + 1):m ] ) ) )
}


################################################################################
################################################################################
##  global_optimise_LM() --- Function to optimize the LM_gama_stpnarpq()
##  function of score test statistic under the null assumption with respect to
##  unknown nuisance parameter (gama). The optimization employes the Brent
##  algorithm applied on the interval [L,U]. To be sure that the a global
##  optimum is found, the optimization is performed at (I-1) consecutive
##  equidistant sub-intervals and then the minimum/maximum over them
##  is taken as global optimum.
##  Input:
##    f=LM_gama_stpnarpq function to be optimize. Always report it in the call.
##    L = Lower bound of search interval for gama.
##    U = Upper bound of search interval for gama.
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
##            gama = optimum value of gama parameter
##            supLM = value of the objective function at the optimum
##            int = list of extremes points of subintervals
################################################################################

global_optimise_LM <- function( f = LM_gama_stpnarpq, b, y, W,
                                p, d, Z, gama_L, gama_U, len, tol = 1e-14 ) {

  if ( min(W) < 0 ) {
    stop('The adjacency matrix W contains negative values.')
  }

  W <- W / Rfast::rowsums(W)
  W[ is.na(W) ] <- 0
  dm <- dim(y)    ;    N <- dm[1]    ;    TT <- dm[2]

  lenm1 <- len - 1
  gami <- supLMi <- numeric(lenm1)

  x <- seq(from = gama_L, to = gama_U, length = len)

  for ( i in 1:lenm1 ) {
    int <- c( x[i], x[i + 1] )
    opt <- optimise( f = LM_gama_stpnarpq, int, b = b, N = N, TT = TT, 
                     y = y, W = W, p = p, d = d, Z = Z, tol = tol )
    gami[i] <- opt$minimum
    supLMi[i] <- opt$objective
  }

  supLM <- min(supLMi)
  gama <- gami[ which(supLMi == supLM)[1] ]
  list(gama = gama, supLM = -supLM, int = x)

}

################################################################################
################################################################################
##  LM_gama_stpnarpq_j() --- Function to optimize the perturbed version of
##  the score test statistic of Smooth Transition model (STNAR)  with p lags,
##  under the null assumption of linearity, with respect to unknown nuisance
##  parameter (gama); it also includes q non time-varying covariates.
##  Input:
##    gama = value of non identifiable nuisance parameter
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
##    LM = (-1) * value of perturbed test statistic at the specified gama
################################################################################

LM_gama_stpnarpq_j <- function(gama, b, N, TT, y, W, p, d, Z, msn) {

  dimz <- ( !is.null(Z) ) * NCOL(Z)
  m <- 1 + 3 * p + max(0,dimz)

  ss <- matrix(0, m, 1)
  out <- matrix(0, m, m)
  hh <- matrix(0, m, m)

  b[ (m - p + 1):m ] <- 0

  for (t in (p+1):TT) {
    xt <- as.vector(W%*%y[,t-d])
    f <- exp(-gama*(xt*xt))
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

  solveHmp  <- solve( H[ 1:(m - p), 1:(m - p) ] )
  Sigma <- B[ (m - p + 1):m , (m - p + 1):m ]-
    H[ (m - p + 1):m, 1:(m - p) ] %*% solveHmp %*% B[ 1:(m - p), (m - p + 1):m ] -
    B[ (m - p + 1):m, 1:(m - p) ] %*% solveHmp %*% H[ 1:(m - p), (m - p + 1):m ] +
    H[ (m - p + 1):m, 1:(m - p) ] %*% solveHmp %*% B[ 1:(m - p), 1:(m - p) ] %*%
    solveHmp %*% H[ 1:(m - p), (m - p + 1):m ]

  - as.numeric( t( S[ (m - p + 1):m ] ) %*% solve( Sigma, S[ (m - p + 1):m ] ) )
}


LM_gama_stpnarpq_j_old <- function(gama, b, N, TT, y, W, p, d, Z, msn){

  m <- 1+3*p+max(0,ncol(Z))

  ss <- as.matrix(rep(0, m))
  out <- matrix(0, nrow=m, ncol=m)
  hh <- matrix(0, nrow=m, ncol=m)

  b[(m-p+1):m] <- 0

  for(t in (p+1):TT){
    xt <- as.vector(W%*%y[,t-d])
    f <- exp(-gama*(xt*xt))
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

