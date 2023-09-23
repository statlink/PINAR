logl_log_lin <- function(ca, N, TT, y, W, wy) {

  nut <- as.vector( wy %*% ca)
  - sum( y[, -1] * nut - exp(nut) )

}

scor_log <- function(ca, N, TT, y, W, wy) {
  nut <- as.vector( wy %*% ca )
  a <-  - Rfast::eachcol.apply(wy, as.vector(y[, -1]) - exp(nut) )
  matrix(a, 3, 1)
}


scor_hessian_out_log <- function(ca, N, TT, y, W, wy) {

  ## scor
  Dt <- exp( as.vector( wy %*% ca ) )
  a <- wy * ( as.vector(y[, -1]) - Dt )
  scor_log <- matrix( - Rfast::colsums(a), 3, 1 )

  ## hess
  hess_log <- crossprod( wy * Dt, wy )

  ## out
  outer_log <- 0
  k <- rep( 1:c(TT - 1), each = N )
  b <- rowsum(a, k)
  for (i in 1: c(TT - 1) )  outer_log <- outer_log + tcrossprod(b[i, ])

  list(scor_log = scor_log, hess_log = hess_log, outer_log = outer_log)
}

################################################################################
################################################################################
##  log_lin_estimnar1() --- function for the constrained estimation of
##  log-linear PNAR model with 1 lag.
##  Inputs:
##    x0 = starting value of the optimization
##    N = number of nodes on the network
##    TT = temporal sample size
##    y = NxTT multivariate count time series
##    W = NxN row-normalized weighted adjacency matrix describing the network
##  output:
##    coeflin = estimated QMLE coefficients
##    selin = standard errors estimates
##    tlin = t test estimates
##    score = value of the score at the optimization point
##    aic_lin = Akaike information criterion (AIC)
##    bic_lin = Bayesian information criterion (BIC)
##    qic_lin = Quasi information criterion (QIC)
################################################################################

log_lin_estimnar1 <- function(y, W) {

  if ( min(W) < 0 ) {
    stop('The adjacency matrix W contains negative values.')
  }

  W <- W / Rfast::rowsums(W)
  W[ is.na(W) ] <- 0

  # Inequality constraints (parameters searched in the stationary region)
  # ca are the parameters to be constrained
  constr <- function(ca, N, TT, y, W, wy) {
    con <- abs(ca[2]) + abs(ca[3]) - 1
    return(con)
  }

  # Jacobian of constraints
  # ca are the parameters to be constrained
  j_constr <- function (ca, N, TT, y, W, wy) {
    j_con <- c( 0, ca[2] / abs(ca[2]), ca[3] / abs(ca[3]) )
    return(j_con)
  }

  dm <- dim(y)    ;    N <- dm[1]    ;    TT <- dm[2]
  wy <- cbind( 1, as.vector( W %*% log1p( y[, -TT] ) ), as.vector( log1p( y[, -TT] ) ) )
  z <- NULL

  ## OLS initital values
  XX <- crossprod(wy)
  Xy <- Rfast::eachcol.apply(wy, as.vector( log1p( y[, -1] ) ) )
  x0 <- solve(XX, Xy)

  # algorithm and and relative tolerance
  opts <- list("algorithm" = "NLOPT_LD_SLSQP", "xtol_rel" = 1.0e-7)

  s_qmle1 <- nloptr::nloptr(x0 = x0, eval_f = logl_log_lin, eval_grad_f = scor_log,
                    eval_g_ineq = constr, eval_jac_g_ineq = j_constr,
                    opts = opts, N = N, TT = TT, y = y, W = W, wy = wy)

  coeflin <- s_qmle1$solution

  ola <- scor_hessian_out_log(coeflin, N, TT, y, W, wy)
  S_lins <- ola$scor_log
  H_lins <- ola$hess_log
  G_lins <- ola$outer_log
  solveH_lins <- solve(H_lins)
  V_lins <- solveH_lins %*% G_lins %*% solveH_lins
  SE_lins <- sqrt( diag(V_lins) )

  tlin <- coeflin/SE_lins

  loglik <-  - s_qmle1$objective
  aic_lins <- 2 * 3 + 2 * s_qmle1$objective
  bic_lins <- log(TT) * 3 + 2 * s_qmle1$objective
  qic_lins <- 2 * sum( H_lins * V_lins ) + 2 * s_qmle1$objective

  list( coeflog = coeflin, selog = SE_lins, tlog = tlin, score = S_lins,
        loglik = loglik, aic_log = aic_lins, bic_log = bic_lins, qic_log = qic_lins )
}



