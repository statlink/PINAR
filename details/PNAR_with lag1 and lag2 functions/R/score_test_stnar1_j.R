################################################################################
################################################################################
##  scor_stnar1_h0() --- function for the computation of score of the nonlinear
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
##  hess_stnar1_h0() --- function for the computation of Hessian matrix
##  of the nonlinear Smooth Transition model (STNAR), with 1 lag,
##  under the null assumption of linearity.
##  output:
##  hh = (-1)*Hessian matrix
################################################################################

################################################################################
################################################################################
##  outer_stnar1_h0() --- function for the computation of information matrix
##  of the nonlinear Smooth Transition model (STNAR), with 1 lag,
##  under the null assumption of linearity.
##  output:
##  out = information matrix
################################################################################

scor_hess_outer_stnar1_h0 <- function(ca, N, TT, y, Wt) {

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
  #out1 <- 0
  k <- rep( 1:c(TT - 1), each = N )
  b1 <- rowsum(a, k)
  #for (i in 1: c(TT - 1) )  out1 <- out1 + tcrossprod(b1[i, ])
  out1 <- crossprod(b1)

  list(scor = scor, hess = hess, out = out1)

}

################################################################################
################################################################################
##  LM_gamma_stnar() --- Function to optimize the score test statistic of
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

LM_gama_stnar <- function(gama, param, N, TT, y, Wt){

  ca <- c(param, gama, 0)
  ola <- scor_hess_outer_stnar1_h0(ca, N, TT, y, Wt)
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
##  optim_gamma_stpnar() --- Function to optimize the LM_gamma_stnar() of the
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

optim_gama_stnar <- function(param, N, TT, y, W) {
  optimise(LM_gama_stnar, c(0.001, 10), param, N, TT, y, W, tol = 1e-9)
}

###########################################################################################
###########################################################################################
##  score_test_stnar1_j() --- function for the computation of bootstrap sup-type
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

score_test_stnar1_j <- function(ca, y, W, J = 499, gama_L, gama_U, len = 10) {

  if ( min(W) < 0 ) {
    stop('The adjacency matrix W contains negative values.')
  }

  W <- W / Rfast::rowsums(W)
  W[ is.na(W) ] <- 0

  dm <- dim(y)   ;    N <- dm[1]    ;    TT <- dm[2]
  LMv <- Sigma <- numeric(len)
  supLMj <- pval <- aveLMj <- apval <- numeric(J)
  LMjv <- matrix(nrow = len, ncol = J)
  V <- 0

  gam <- seq(from = gama_L, to = gama_U, length = len)
  Wt <- as.vector( W %*% y[, -TT] )

  for ( i in 1:len ) {
    ca[4] <- gam[i]
    ola <- scor_hess_outer_stnar1_h0(ca, N, TT, y, Wt)
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

  ### Bootstrap
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

