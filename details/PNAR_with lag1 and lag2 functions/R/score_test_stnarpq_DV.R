score_test_stnarpq_DV <- function(b, y, W, p, d, Z = NULL, gama_L, gama_U, len = 100) {

  if ( min(W) < 0 ) {
    stop('The adjacency matrix W contains negative values.')
  }

  W <- W / Rfast::rowsums(W)
  W[ is.na(W) ] <- 0

  dm <- dim(y)   ;    N <- dm[1]    ;    TT <- dm[2]

  dimz <- ( !is.null(Z) ) * NCOL(Z)
  m <- 1 + 3 * p + max(0, dimz)
  b[ (m - p + 1):m ] <- 0

  LMv <- numeric(len)
  gam <- seq(from = gama_L, to = gama_U, length = len)
  V <- 0

  z <- W %*% y
  z2 <- z^2

  wy1 <- NULL
  for ( ti in (p + 1):TT )
    wy1 <- rbind( wy1, cbind(z[, (ti - 1):(ti - p)], y[, (ti - 1):(ti - p)], Z) )
  wy1 <- cbind(1, wy1)

  lambdat <- as.vector( wy1 %*% b[1:c(m - p)] )
  com <- ( as.vector( y[, -c(1:p)] ) - lambdat ) / lambdat
  ct <- as.vector( y[, -c(1:p)] ) / lambdat^2
  k <- rep( 1:c(TT - p), each = N )


  wy2 <- NULL
  f <- exp(-gam[1] * z2)
  for ( ti in (p + 1):TT )  wy2 <- rbind( wy2, z[, (ti - 1):(ti - p)] * f[, (ti - d)] )
  wy <- cbind(wy1, wy2)

  ## scor
  a <- wy * com
  S <- matrix( Rfast::colsums(a), m, 1 )

  ## hess
  H <- crossprod(wy * ct, wy)

  ## out
  #out1 <- 0
  b1 <- rowsum(a, k)
  #for ( i in 1: c(TT - p) )  out1 <- out1 + tcrossprod(b1[i, ])
  B <- crossprod(b1)
  
  solveHmp <- solve( H[ 1:(m - p), 1:(m - p) ] )
  Sigma <- B[ (m - p + 1):m , (m - p + 1):m ] -
    H[ (m - p + 1):m, 1:(m - p) ] %*% solveHmp %*% B[ 1:(m - p), (m - p + 1):m ] -
    B[ (m - p + 1):m, 1:(m - p) ] %*% solveHmp %*% H[ 1:(m - p), (m - p + 1):m ] +
    H[ (m - p + 1):m, 1:(m - p) ] %*% solveHmp %*% B[ 1:(m - p), 1:(m - p) ] %*%
    solveHmp %*% H[ 1:(m - p), (m - p + 1):m ]
  LMv[1] <-  - as.numeric( crossprod( S[ (m - p + 1):m ], solve( Sigma, S[ (m - p + 1):m ] ) ) )

  for ( j in 2:len ) {
    wy2 <- NULL
    f <- exp(-gam[j] * z2)
    for ( ti in (p + 1):TT )  wy2 <- rbind( wy2, z[, (ti - 1):(ti - p)] * f[, (ti - d)] )
    wy <- cbind(wy1, wy2)

    ## scor
    a <- wy * com
    S <- matrix( Rfast::colsums(a), m, 1 )

    ## hess
    H <- crossprod(wy * ct, wy)

    ## out
    #out1 <- 0
    b1 <- rowsum(a, k)
    #for ( i in 1: c(TT - p) )  out1 <- out1 + tcrossprod(b1[i, ])
    out1 <- crossprod(b1)
    B <- out1

    solveHmp <- solve( H[ 1:(m - p), 1:(m - p) ] )
    Sigma <- B[ (m - p + 1):m , (m - p + 1):m ] -
      H[ (m - p + 1):m, 1:(m - p) ] %*% solveHmp %*% B[ 1:(m - p), (m - p + 1):m ] -
      B[ (m - p + 1):m, 1:(m - p) ] %*% solveHmp %*% H[ 1:(m - p), (m - p + 1):m ] +
      H[ (m - p + 1):m, 1:(m - p) ] %*% solveHmp %*% B[ 1:(m - p), 1:(m - p) ] %*%
      solveHmp %*% H[ 1:(m - p), (m - p + 1):m ]
      LMv[j] <-  - as.numeric( crossprod( S[ (m - p + 1):m ], solve( Sigma, S[ (m - p + 1):m ] ) ) )
      V <- V + abs( sqrt( -LMv[j] ) - sqrt( -LMv[j - 1] ) )
  }  ##  end  for (i in 2:len) {

  supLM <- as.numeric( max(-LMv) )
  DV <- 1 - pchisq(supLM, p) + V * supLM^(0.5 * (p - 1)) * exp(-0.5 * supLM) * 2^(-0.5 * p) / gamma(p/2)
  list(DV = DV, supLM = supLM)
}




################################################################################
## Function for for testing linearity of Poisson NAR model, with p lags, PNAR(p),
## versus nonlinear Smooth Transition alternative model (STPNAR)
## it also includes q non time-varying covariates
################################################################################


################################################################################
################################################################################
##  LM_gamma_stnarpq() --- Function to optimize the score test statistic of
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

# LM_gamma_stnarpq <- function(gamma, b, N, TT, y, W, p, d, Z){
#
#   m <- 1+3*p+max(0,ncol(Z))
#
#   ss <- as.matrix(rep(0, m))
#   out <- matrix(0, nrow=m, ncol=m)
#   hh <- matrix(0, nrow=m, ncol=m)
#
#   b[(m-p+1):m] <- 0
#
#   for(t in (p+1):TT){
#     xt <- as.vector(W%*%y[,t-d])
#     f <- exp(-gamma*(xt*xt))
#     Xp <- W%*%y[,(t-1):(t-p)]
#     Xt <- cbind(1, Xp, y[,(t-1):(t-p)], Z, Xp*f)
#     lambdat <- Xt[ ,1:(m-p) ] %*% b[ 1:(m-p) ]
#     Dt <- diag(1/as.vector(lambdat))
#
#     s <- t(Xt)%*%Dt%*%(y[,t]-lambdat)
#     ss <- ss + s
#     out <- out + s%*%t(s)
#
#     Ct <- diag(y[,t])%*%diag(1/as.vector(lambdat^2))
#     hh <- hh + t(Xt)%*%Ct%*%Xt
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
#   LM <- -1 * LM
#   return( LM )
# }



###########################################################################################
###########################################################################################
##  score_test_stnarpq_DV() --- function for the computation of Davies bound of sup-type
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

# score_test_stnarpq_DV <- function(b, N, TT, y, W, p, d, Z, gamma_L, gamma_U, l=100){
#
#   m <- 1+3*p+max(0,ncol(Z))
#   LMv <- vector()
#
#   gam <- seq(from=gamma_L, to=gamma_U, length.out = l)
#
#   V <- 0
#
#   LMv[1] <- LM_gamma_stnarpq(gam[1], b, N, TT, y, W, p, d, Z)
#
#   for(i in 2:l){
#
#     gamma <- gam[i]
#     LMv[i] <- LM_gamma_stnarpq(gamma, b, N, TT, y, W, p, d, Z)
#     V <- V +  abs(sqrt(-LMv[i]) - sqrt(-LMv[i-1]))
#
#   }
#
#   supLM <- as.numeric(max(-LMv))
#
#   DV <- 1-pchisq(supLM,p) + V*supLM^(1/2*(p-1))*exp(-supLM/2)*2^(-p/2)/gamma(p/2)
#
#   return(list(DV=DV, supLM=supLM))
# }

