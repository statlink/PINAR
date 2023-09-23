################################################################################
## Function for for testing linearity of Poisson NAR model, PNAR, with 1 lag
## versus Threshold alternative model (TNAR)
################################################################################

set.seed(123)

################################################################################
################################################################################
##  scor_tar1_h0() --- function for the computation of score of the nonlinear
##  Threshold model (TNAR), with 1 lag, under the null assumption
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

scor_tar1_h0 <- function(ca, N, TT, y, W) {
  
  ca[4:6] <- 0
  Wt <- as.vector( W %*% y[, -TT] )
  It <- ( Wt <= ca[7] )
  wy <- cbind( 1, Wt, as.vector( y[, -TT] ), It, Wt * It, as.vector( y[, -TT] ) * It, 0 )
  lambdat <- wy[, 1:3] %*% ca[1:3]
  crossprod( wy, ( as.vector(y[, -1]) - lambdat ) / lambdat )

}

scor_tar1_h0_old <- function(c, N, TT, y, W) {
  
  unon <- matrix(1, N, 1)
  zeron <- matrix(0, N, 1)
  c[4:6] <- 0
  ss <- matrix(0, 7, 1)
  
  for(t in 2:TT){
    It <- 0+(W%*%y[,t-1]<=c[7])
    Xt <- cbind(unon, W%*%y[,t-1], y[,t-1], It, W%*%y[,t-1]*It, y[,t-1]*It, zeron)
    lambdat <- Xt[,1:3]%*%c[1:3]
    Dt <- diag(1/as.vector(lambdat))
    ss <- ss + t(Xt)%*%Dt%*%(y[,t]-lambdat)
  }
  
  return(ss)
}


################################################################################
################################################################################
##  outer_tar1_h0() --- function for the computation of information matrix
##  of the nonlinear Threshold model (TNAR), with 1 lag, 
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

outer_tar1_h0 <- function(ca, N, TT, y, W) {

  ca[4:6] <- 0
  Wt <- as.vector( W %*% y[, -TT] )
  It <- ( Wt <= ca[7] )
  wy <- cbind( 1, Wt, as.vector( y[, -TT] ), It, Wt * It, as.vector( y[, -TT] ) * It, 0 )
  lambdat <- as.vector( wy[, 1:3] %*% ca[1:3] )
  a <- wy  * ( as.vector(y[, -1]) - lambdat ) / lambdat   
  out <- 0
  k <- rep( 1:c(TT - 1), each = N )
  b <- rowsum(a, k)
  for (i in 1: c(TT - 1) )  out <- out + tcrossprod(b[i, ])
  out
 
}

outer_tar1_h0_old <- function(c, N, TT, y, W){
  
  unon <- as.matrix(rep(1, N))
  zeron <- as.matrix(rep(0, N))
  c[4]=c[5]=c[6]=0
  out <- matrix(0, 7, 7)
  
  for(t in 2:TT){
    It <- 0+(W%*%y[,t-1]<=c[7])
    Xt <- cbind(unon, W%*%y[,t-1], y[,t-1], It, W%*%y[,t-1]*It, y[,t-1]*It, zeron)
    lambdat <- Xt[,1:3]%*%c[1:3]
    Dt <- diag(1/as.vector(lambdat))
    ss <- t(Xt)%*%Dt%*%(y[,t]-lambdat)
    out <- out + ss%*%t(ss)
  }
  
  return(out)
}


################################################################################
################################################################################
##  hess_tar1_h0() --- function for the computation of Hessian matrix
##  of the nonlinear Threshold model (TNAR), with 1 lag, 
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

hess_tar1_h0 <- function(ca, N, TT, y, W) {
  
  ca[4:6] <- 0
  Wt <- as.vector( W %*% y[, -TT] )
  It <- ( Wt <= ca[7] )
  wy <- cbind( 1, Wt, as.vector( y[, -TT] ), It, Wt * It, as.vector( y[, -TT] ) * It, 0 )
  ct <- as.vector(y[, -1]) / as.vector( wy %*% ca )^2
  crossprod( wy * ct, wy )

}


hess_tar1_h0_old <- function(c, N, TT, y, W){
  
  unon <- as.matrix(rep(1, N))
  zeron <- as.matrix(rep(0, N))
  c[4]=c[5]=c[6]=0
  hh <- matrix(0, 7, 7)
  
  for(t in 2:TT){
    It <- 0+(W%*%y[,t-1]<=c[7])
    Xt <- cbind(unon, W%*%y[,t-1], y[,t-1], It, W%*%y[,t-1]*It, y[,t-1]*It, zeron)
    lambdat <- Xt[,1:3]%*%c[1:3]
    Ct <- diag(y[,t])%*%diag(1/as.vector(lambdat^2))
    hes <- t(Xt)%*%Ct%*%Xt
    hh <- hh + hes
  }
  
  return(hh)
}


################################################################################
################################################################################
##  scor_tar1_j() --- function for the computation of the perturbed version
##  of the score of the nonlinear Threshold model (TNAR), with 1 lag, 
##  under the null assumption of linearity. 
##  Inputs:
##     c = estimated parameters from the linear model
##     N = number of nodes on the network
##     TT = temporal sample size
##     y = NxTT multivariate count time series
##     W = NxN row-normalized weighted adjacency matrix describing the network
##     msn = TTx1 vector of standard normal noises
##  output:
##     ss = vector of perturbed quasi score
################################################################################

scor_tar1_j <- function(ca, N, TT, y, W, msn) {

  ca[4:6] <- 0
  Wt <- as.vector( W %*% y[, -TT] )
  It <- ( Wt <= ca[7] )
  wy <- cbind( 1, Wt, as.vector( y[, -TT] ), It, Wt * It, as.vector( y[, -TT] ) * It, 0 )
  lambdat <- as.vector( wy %*% ca )
  crossprod( wy, ( as.vector(y[, -1]) - lambdat ) / lambdat * rep( msn[-1], each = N ) )
  
}

scor_tar1_j_old <- function(c, N, TT, y, W, msn){
  
  unon <- as.matrix(rep(1, N))
  zeron <- as.matrix(rep(0, N))
  
  c[4]=c[5]=c[6]=0
  ss <- matrix(0, 7, 1)
  
  for(t in 2:TT){
    
    It <- 0+(W%*%y[,t-1]<=c[7])
    Xt <- cbind(unon, W%*%y[,t-1], y[,t-1], It, W%*%y[,t-1]*It, y[,t-1]*It, zeron)
    lambdat <- Xt[,1:3]%*%c[1:3]
    Dt <- diag(1/as.vector(lambdat))
    sh <- t(Xt)%*%Dt%*%(y[,t]-lambdat)*msn[t]
    ss <- ss + sh
  }
  
  return(ss)
}


scor_hess_outer_tar1_h0 <- function(ca, N, TT, y, W) {
  
  ## scor
  ca[4:6] <- 0
  Wt <- as.vector( W %*% y[, -TT] )
  It <- ( Wt <= ca[7] )
  wy <- cbind( 1, Wt, as.vector( y[, -TT] ), It, Wt * It, as.vector( y[, -TT] ) * It, 0 )
  lambdat <- as.vector( wy[, 1:3] %*% ca[1:3] )
  a <- wy  * ( as.vector(y[, -1]) - lambdat ) / lambdat   
  S <- matrix(Rfast::colsums(a), 7, 1)
  
  ## hess
  ct <- as.vector(y[, -1]) / lambdat^2
  H <- crossprod( wy * ct, wy )

  ## outer
  out <- 0
  k <- rep( 1:c(TT - 1), each = N )
  b <- rowsum(a, k)
  for (i in 1: c(TT - 1) )  out <- out + tcrossprod(b[i, ])
   
  list( S = S, H = H, B = out )    
}


################################################################################
################################################################################
##  LM_gamma_tpnar() --- Function to optimize the score test statistic of
##  Threshold model (TNAR)
##  under the null assumption with respect to unknown nuisance parameter (gamma)
##  Input:
##    gamma = value of nuisance parameter
##    par = estimated parameters from the linear model
##    N = number of nodes on the network
##    TT = temporal sample size
##    y = NxTT multivariate count time series
##    W = NxN row-normalized weighted adjacency matrix describing the network
##  output:
##    LM = (-1)*value of test statistic at the specified gamma
################################################################################

LM_gamma_tpnar <- function(gama, par, N, TT, y, W) {
  
  ca <- c(par, 0, 0, 0, gama)
  ola <- scor_hess_outer_tar1_h0(ca, N, TT, y, W)
  S <- ola$S
  H <- ola$H
  B <- ola$B
  
  solveH13 <- solve( H[1:3, 1:3] )
  Sigma <- B[4:6, 4:6] - H[4:6, 1:3] %*% solveH13 %*% B[1:3, 4:6] -
    B[4:6, 1:3] %*% solveH13 %*% H[1:3, 4:6] +
    H[4:6, 1:3] %*% solveH13 %*% B[1:3, 1:3] %*% solveH13 %*% H[1:3, 4:6]
  
  - as.numeric( crossprod(S[4:6], solve( Sigma, S[4:6] ) ) )
  
}

LM_gamma_tpnar_old <- function(gama, par, N, TT, y, W){
  
  c <- c(par, 0, 0, 0, gama)
  S <- scor_tar1_h0_old(c, N, TT, y, W)
  H <- hess_tar1_h0_old(c, N, TT, y, W)
  B <- outer_tar1_h0_old(c, N, TT, y, W)
  
  Sigma <- B[4:6,4:6]-H[4:6,1:3]%*%solve(H[1:3,1:3])%*%B[1:3,4:6]-
    B[4:6,1:3]%*%solve(H[1:3,1:3])%*%H[1:3,4:6]+
    H[4:6,1:3]%*%solve(H[1:3,1:3])%*%B[1:3,1:3]%*%solve(H[1:3,1:3])%*%H[1:3,4:6]
  
  LM <- as.numeric(t(S[4:6])%*%solve(Sigma)%*%S[4:6])
  LM <- -1*LM
  return(LM)
}


################################################################################
################################################################################
##  optim_gamma_tpnar() --- Function to optimize the LM_gamma_tpnar() of the
##  score test statistic under the null assumption with respect to unknown
##  nuisance parameter (gamma)
##  Input:
##    par = estimated parameters from the linear model
##    N = number of nodes on the network
##    TT = temporal sample size
##    y = NxTT multivariate count time series
##    W = NxN row-normalized weighted adjacency matrix describing the network
##  output: s_gam, a list with two values
##              minimum = optimum value of gamma parameter
##              objective = (-1)*value of the objective function at the optimum
################################################################################

optim_gamma_tpnar <- function(par, N, TT, y, W) { 
  Z <- W %*% y
  qq <- Rfast2::rowQuantile( Z, probs = c(0.1, 0.9) )
  optimise(f = LM_gamma_tpnar, c( min(qq[, 1]), max(qq[, 2]) ), par, N, TT, y, W)
}

optim_gamma_tpnar_old <- function(par, N, TT, y, W){
  
  #derivative-free optimization
  #golden section search algorithm looks for solution in the interval
  # min(10%-quantiles) - max(90%-quantiles)
  
  Z=W%*%y
  qq <- matrix(0, 2, nrow(Z))
  for (k in 1:nrow(Z)){
    qq[,k] <- as.matrix(quantile(Z[k,], prob=c(0.10, 0.90)))
  }
  
  s_gam <- optimise(f=LM_gamma_tpnar_old, c(min(qq[1,]),max(qq[2,])), par, N, TT, y, W)
  
 
  return(s_gam)
}


###########################################################################################
###########################################################################################
##  score_test_tar1_j() --- function for the computation of bootstrap sup-type
##  test for testing linearity of PNAR model, with 1 lag, versus the nonlinear 
##  Threshold model (TNAR) alternative. 
##  Inputs:
##     c = estimated parameters from the linear model
##     N = number of nodes on the network
##     TT = temporal sample size
##     y = NxTT multivariate count time series
##     W = NxN row-normalized weighted adjacency matrix describing the network
##     J = number of bootstrap replications
##     l = number of values in the grid of values for gamma parameter 
##     gam_L = lower bound of grid of values for gamma parameter
##     gam_U = upper bound of grid of values for gamma parameter
##  output: a list with the following objects
##        supLM = value of the sup test statistic in the sample y
##        aveLM = value of the mean test statistic in the sample y
##        pJ = bootstrap p-value sup test
##        apJ = bootstrap p-value mean test
##        cpJ = adjusted version of bootstrap p-value sup test
##        acpJ = adjusted version of bootstrap p-value mean test
##        LM = values of test statistic for different values in the grid
##        supLMj = values of perturbed sup test statistic for different values in the grid
##        aveLMj = values of perturbed mean test statistic for different values in the grid
############################################################################################

score_test_tnar1_j <- function(c, N, TT, y, W, J, len, gam_L, gam_U) {
  
  
  LMv <- vector()
  supLMj <- vector()
  pval <- vector()
  LMjv <- vector()
  aveLMj <- vector()
  apval <- vector()
  
  gam <- seq(from = gam_L, to = gam_U, length = len)
  
  for ( i in 1:len ) {
    ca[7] <- gam[i]
    ola <- scor_hess_outer_tar1_h0(ca, N, TT, y, W)
    S <- ola$S
    H <- ola$H
    B <- ola$B 
    Sigma <- B[4:6, 4:6] - H[4:6, 1:3] %*% solveH13 %*% B[1:3, 4:6] -
             B[4:6, 1:3] %*% solveH13 %*% H[1:3, 4:6] +
             H[4:6, 1:3] %*% solveH13 %*% B[1:3, 1:3] %*% solveH13 %*% H[1:3, 4:6]
    LMv[i] <- as.numeric( crossprod(S[4:6], solve( Sigma, S[4:6] ) ) )
  }
  
  supLM <- as.numeric( max(LMv) )
  aveLM <- as.numeric( mean(LMv) )
  
  for ( j in 1:J ) {
    
    snj <- rnorm(TT)
    for ( i in 1:len) {
      ca[7]<- gam[i]
      Sj <- scor_tar1_j(c, N, TT, y, W, snj)
      H <- hess_tar1_h0(c, N, TT, y, W)
      B <- outer_tar1_h0(c, N, TT, y, W)
      Sigma <- B[4:6, 4:6] - H[4:6, 1:3] %*% solveH13 %*% B[1:3, 4:6] -
               B[4:6, 1:3] %*% solveH13 %*% H[1:3, 4:6] +
               H[4:6, 1:3] %*% solveH13 %*% B[1:3, 1:3] %*% solveH13 %*% H[1:3, 4:6]
      LMjv[i] <- as.numeric( crossprod(Sj[4:6], solve( Sigma, Sj[4:6] ) ) )
    }
    
    supLMj[j] <- max(LMjv)
    pval[j] <- ( supLMj[j] >= supLM )
    aveLMj[j] <- mean(LMjv)
    apval[j] <- ( aveLMj[j] >= aveLM )
  }
  
  pJ <- sum(pval) / J
  apJ <- sum(apval) / J
  cpJ <- ( sum(pval) + 1 ) / (J + 1)
  capJ <- ( sum(apval) + 1 ) / (J + 1)
  
  list( supLM = supLM, aveLM = aveLM, pJ = pJ, apJ = apJ, cpJ = cpJ, 
        capJ = capJ, LM = LMv, supLMj = supLMj, aveLMj = aveLMj )
}


