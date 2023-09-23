################################################################################
## Function for for testing linearity of Poisson NAR model, PNAR, with 1 lag
## versus Threshold alternative model (TNAR)
################################################################################


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

#Wt <- as.vector( W %*% y[, -TT] )

scor_hess_outer_tnar1_h0 <- function(ca, N, TT, y, Wt) {

  ## scor
  ca[4:6] <- 0
  It <- ( Wt <= ca[7] )
  wy <- cbind( 1, Wt, as.vector( y[, -TT] ), It, Wt * It, as.vector( y[, -TT] ) * It, 0 )
  lambdat <- as.vector( wy[, 1:3] %*% ca[1:3] )
  a <- wy  * ( as.vector(y[, -1]) - lambdat ) / lambdat
  S <- matrix(Rfast::colsums(a), 7, 1)

  ## hess
  ct <- as.vector(y[, -1]) / lambdat^2
  H <- crossprod( wy * ct, wy )

  ## outer
  #out1 <- 0
  k <- rep( 1:c(TT - 1), each = N )
  b1 <- rowsum(a, k)
  #for (i in 1: c(TT - 1) )  out1 <- out1 + tcrossprod(b1[i, ])
  out1 <- crossprod(b1)
  list( S = S, H = H, B = out1 )
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

LM_gama_tnar <- function(gama, param, N, TT, y, Wt) {

  ca <- c(param, 0, 0, 0, gama)
  ola <- scor_hess_outer_tnar1_h0(ca, N, TT, y, Wt)
  S <- ola$S
  H <- ola$H
  B <- ola$B

  solveH13 <- solve( H[1:3, 1:3] )
  Sigma <- B[4:6, 4:6] - H[4:6, 1:3] %*% solveH13 %*% B[1:3, 4:6] -
    B[4:6, 1:3] %*% solveH13 %*% H[1:3, 4:6] +
    H[4:6, 1:3] %*% solveH13 %*% B[1:3, 1:3] %*% solveH13 %*% H[1:3, 4:6]

  as.numeric( crossprod(S[4:6], solve( Sigma, S[4:6] ) ) )

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

#Wy <- W %*% y
optim_gama_tnar <- function(param, N, TT, y, Wy, Wt) {
  qq <- Rfast2::rowQuantile( Wy, probs = c(0.1, 0.9) )
  optimise( LM_gama_tnar, c( min(qq[, 1]), max(qq[, 2]) ),
            param, N, TT, y, Wt, tol = 1e-09, maximum = TRUE )
}

###########################################################################################
###########################################################################################
##  score_test_tnar1_j() --- function for the computation of bootstrap sup-type
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

score_test_tnar1_j <- function(ca, y, W, J = 499, gama_L, gama_U, len = 10) {

  if ( min(W) < 0 ) {
    stop('The adjacency matrix W contains negative values.')
  }

  W <- W / Rfast::rowsums(W)
  W[ is.na(W) ] <- 0

  dm <- dim(y)   ;    N <- dm[1]    ;    TT <- dm[2]
  solveSigma <- list()
  LMv <- supLMj <- aveLMj <- pval <- apval <- LMjv <- numeric(J)

  gam <- seq(from = gama_L, to = gama_U, length = len)
  Wt <- as.vector( W %*% y[, -TT] )

  for ( i in 1:len ) {
    ca[7] <- gam[i]
    ola <- scor_hess_outer_tnar1_h0(ca, N, TT, y, Wt)
    S <- ola$S
    H <- ola$H
    B <- ola$B
    solveH13 <- solve( H[1:3, 1:3] )
    Sigma <- B[4:6, 4:6] - H[4:6, 1:3] %*% solveH13 %*% B[1:3, 4:6] -
             B[4:6, 1:3] %*% solveH13 %*% H[1:3, 4:6] +
             H[4:6, 1:3] %*% solveH13 %*% B[1:3, 1:3] %*% solveH13 %*% H[1:3, 4:6]
    solveSigma[[ i ]] <- solve(Sigma)
    LMv[i] <- as.numeric( crossprod(S[4:6], solveSigma[[ i ]] %*% S[4:6] ) )
  }

  supLM <- as.numeric( max(LMv) )
  aveLM <- as.numeric( mean(LMv) )
  ca[4:6] <- 0

  ## Bootstrap
  mat <- Rfast::matrnorm(TT, J)
  rows <- rep(2:TT, each = N)
  mat <- mat[rows, ]
  wy1 <- cbind( 1, Wt, as.vector( y[, -TT] ) )

  for ( j in 1:J ) {
    for ( i in 1:len) {
      #ca[7] <- gam[i]
      It <- ( Wt <= gam[i] )
      #wy <- cbind( 1, Wt, as.vector( y[, -TT] ), It, Wt * It, as.vector( y[, -TT] ) * It, 0 )
      wy2 <-  cbind( It, Wt * It, as.vector( y[, -TT] ) * It )
      lambdat <- as.vector( wy1 %*% ca[1:3] )
      #Sj <- crossprod( wy2, ( as.vector(y[, -1]) - lambdat ) / lambdat * mat[, j] )
      Sj <- crossprod( wy2, mat[, j] * ( as.vector(y[, -1]) / lambdat - 1) )
      LMjv[i] <- as.numeric( crossprod(Sj, solveSigma[[ i ]] %*% Sj ) )
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



