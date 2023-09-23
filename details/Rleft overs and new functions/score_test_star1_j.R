################################################################################
################################################################################
##  scor_star1_h0() --- function for the computation of score of the nonlinear
##  Smooth Transition model (STNAR), with 1 lag, under the null assumption
##  of linearity. 
##  Inputs:
##  c = estimated parameters from the linear model
##  N = number of nodes on the network
##  TT = temporal sample size
##  y = NxTT multivariate count time series
##  W = NxN row-normalized weighted adjacency matrix describing the network
##  output:
##  ss = vector of quasi score
################################################################################

################################################################################
################################################################################
##  hess_star1_h0() --- function for the computation of Hessian matrix
##  of the nonlinear Smooth Transition model (STNAR), with 1 lag, 
##  under the null assumption of linearity.
##  output:
##  hh = (-1)*Hessian matrix
################################################################################

################################################################################
################################################################################
##  outer_star1_h0() --- function for the computation of information matrix
##  of the nonlinear Smooth Transition model (STNAR), with 1 lag, 
##  under the null assumption of linearity.
##  output:
##  out = information matrix
################################################################################

scor_hess_outer_star1_h0 <- function(ca, N, TT, y, Wt) {

  ## scor  
  ca[5] <- 0
  wy <- cbind( 1, Wt, as.vector( y[, -TT] ), 0, exp( - ca[4] * Wt^2 ) * Wt )
  lambdat <- as.vector( wy[, 1:3] %*% ca[1:3] )
  a <- wy  * ( as.vector(y[, -1]) - lambdat ) / lambdat   
  scor <- matrix( Rfast::colsums(a), 5, 1 )
  
  ## hess
  ct <- as.vector(y[, -1]) /  lambdat^2
  hess <- crossprod( wy * ct, wy )

  ## outer
  out <- 0
  k <- rep( 1:c(TT - 1), each = N )
  b <- rowsum(a, k)
  for (i in 1: c(TT - 1) )  out <- out + tcrossprod(b[i, ])
  
  list(scor = scor, hess = hess, out = out)

}



################################################################################
################################################################################
##  scor_star1_j() --- function for the computation of the perturbed version
##  of the score of the nonlinear Smooth Transition model (STNAR), with 1 lag, 
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

scor_star1_j <- function(ca, N, TT, y, Wt, msn) {

  ca[5] <- 0
  wy <- cbind( 1, Wt, as.vector( y[, -TT] ), 0, exp( - ca[4] * Wt^2 ) * Wt )
  lambdat <- as.vector( wy[, 1:3] %*% ca[1:3] )
  crossprod( wy, ( as.vector(y[, -1]) - lambdat ) / lambdat * rep( msn[-1], each = N ) ) 
  
}

################################################################################
################################################################################
##  LM_gamma_stpnar() --- Function to optimize the score test statistic of
##  Smooth Transition model (STNAR)
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

LM_gamma_stpnar <- function(gama, par, N, TT, y, Wt){
  
  ca <- c(par, gama, 0)
  ola <- scor_hess_outer_star1_h0(ca, N, TT, y, Wt)
  S <- ola$scor
  H <- ola$hess
  B <- ola$out
  
  solveH13 <- solve(H[1:3, 1:3])
  Sigma <- B[5, 5] - H[5, 1:3] %*% solveH13 %*% B[1:3, 5] -
    B[5, 1:3] %*% solveH13 %*% H[1:3, 5] +
    H[5, 1:3] %*% solveH13 %*% B[1:3, 1:3] %*% solveH13 %*% H[1:3, 5]
  
  - as.numeric(S[5])^2 / as.numeric(Sigma)
}

################################################################################
################################################################################
##  optim_gamma_stpnar() --- Function to optimize the LM_gamma_stpnar() of the
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

optim_gamma_stpnar <- function(par, N, TT, y, W) {
 optimise(f = LM_gamma_stpnar, c(0.001, 10), par, N, TT, y, W, tol = 1e-6)
}

###########################################################################################
###########################################################################################
##  score_test_star1_j() --- function for the computation of bootstrap sup-type
##  test for testing linearity of PNAR model, with 1 lag, versus the nonlinear 
##  Smooth Transition model (STNAR) alternative. 
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
##        DV = Davies bound of p-values for sup test
##        pJ = bootstrap p-value sup test
##        apJ = bootstrap p-value mean test
##        cpJ = adjusted version of bootstrap p-value sup test
##        acpJ = adjusted version of bootstrap p-value mean test
##        LM = values of test statistic for different values in the grid
##        supLMj = values of perturbed sup test statistic for different values in the grid
##        aveLMj = values of perturbed mean test statistic for different values in the grid
############################################################################################

score_test_star1_j <- function(ca, N, TT, y, W, J, len, gam_L, gam_U) {
  
  LMv <- Sigma <- numeric(len)
  supLMj <- pval <- aveLMj <- apval <- numeric(J)
  LMjv <- matrix(nrow = len, ncol = J)
  V <- 0

  gam <- seq(from = gam_L, to = gam_U, length = len)
  Wt <- as.vector( W %*% y[, -TT] )
   
  for ( i in 1:len ) {
    ca[4] <- gam[i]
    ola <- scor_hess_outer_star1_h0(ca, N, TT, y, Wt)
    S <- ola$scor
    H <- ola$hess
    B <- ola$out
    
    solveH13 <- solve( H[1:3, 1:3] )
    Sigma[i] <- B[5, 5] - H[5, 1:3] %*% solveH13 %*% B[1:3, 5] -
      B[5, 1:3] %*% solveH13 %*% H[1:3, 5] +
      H[5, 1:3] %*% solveH13 %*% B[1:3, 1:3] %*% solveH13 %*% H[1:3, 5]
    LMv[i] <- as.numeric( S[5] )^2 / as.numeric(Sigma[i])
    if ( i > 1 )  V <- V +  abs( sqrt( LMv[i] ) - sqrt( LMv[i - 1] ) )   
  }  ##  end for (i in 1:l)
  
  supLM <- as.numeric( max(LMv) )
  aveLM <- as.numeric( mean(LMv) )
  DV <- pchisq(supLM, 1, lower.tail = FALSE) + V * exp(-0.5 * supLM) / ( sqrt(2) * gamma(1/2) )
  
  ca[5] <- 0
  wy2 <- cbind( 1, Wt, as.vector( y[, -TT] ) )
  lambdat <- as.vector( wy2[, 1:3] %*% ca[1:3] )
  wy2 <- NULL
  mat <- Rfast::matrnorm(TT, J)
  rows <- rep(2:TT, each = N)
  mat <- mat[rows, ] 
  com <- mat * ( as.vector(y[, -1]) - lambdat ) / lambdat * Wt

  for ( i in 1:len ) {  
    ca[4] <- gam[i]
    wy <- exp( - ca[4] * Wt^2 )
    Sj <- Rfast::eachcol.apply(com, wy)
    LMjv[i, ] <- Sj^2 / Sigma[i]
  } ##  end for ( i in 1:1en )
    
   supLMj <- Rfast::colMaxs(LMjv, TRUE)
   pval <- ( supLMj >= supLM )
   aveLMj <- Rfast::colmeans(LMjv)
   apval <- ( aveLMj >= aveLM )
  
  pJ <- sum(pval) / J
  apJ <- sum(apval) / J
  cpJ <- ( sum(pval) + 1 ) / (J + 1)
  capJ <- ( sum(apval) + 1) / (J + 1)
  
  list( supLM = supLM, aveLM = aveLM, DV = DV, pJ = pJ, apJ = apJ, cpJ = cpJ, 
       capJ = capJ, LM = LMv, supLMj = supLMj, aveLMj = aveLMj )
}

