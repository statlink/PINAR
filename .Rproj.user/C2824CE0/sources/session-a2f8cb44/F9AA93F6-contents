################################################################################
## Function for estimating log-linear Poisson NAR model with p lags, PNAR(p),
## and q non time-varying covariates
################################################################################

################################################################################
################################################################################
##  logl_log_linpq() --- function for the computation of log-likelihood of the
##  log-linear PNAR model with p lags and q covariates. Inputs:
##     b = parameters of the models in the following order:
##         (intercept, p network effects, p ar effects, covariates)
##     N = number of nodes on the network
##     TT = temporal sample size
##     y = NxTT multivariate count time series
##     W = NxN row-normalized weighted adjacency matrix describing the network
##     p = number of lags in the model
##     Z = Nxq matrix of covariates (one for each column), where q is the number of
##         covariates in the model.
##  output:
##     g = (-1)* independence quasi log-likelihood
################################################################################

.logl_log_linpq <- function(b, N, TT, y, W, wy, p, Z) {
  nut <- as.vector( wy %*% b)
  - sum( as.vector( y[, -c(1:p)] ) * nut - exp(nut) )
}


################################################################################
################################################################################
##  scor_logpq() --- function for the computation of score of the
##  log-linear PNAR model with p lags and q covariates. Inputs:
##     b = parameters of the models in the following order:
##         (intercept, p network effects, p ar effects, covariates)
##     N = number of nodes on the network
##     TT = temporal sample size
##     y = NxTT multivariate count time series
##     W = NxN row-normalized weighted adjacency matrix describing the network
##     p = number of lags in the model
##     Z = Nxq matrix of covariates (one for each column), where q is the number of
##         covariates in the model.
##  output:
##     ss = (-1)* vector of quasi score
################################################################################

.scor_logpq <- function(b, N, TT, y, W, wy, p, Z) {
  nut <- as.vector( wy %*% b)
  a <-  - Rfast::eachcol.apply(wy, as.vector( y[, -c(1:p)] ) - exp(nut) )
  matrix(a, 2 * p + 1 + max( 0, ncol(Z) ), 1)
}

################################################################################
################################################################################
##  hess_logpq() --- function for the computation of Hessian matrix of the
##  log-linear PNAR model with p lags and q covariates. Inputs:
##     b = parameters of the models
##     N = number of nodes on the network
##     TT = temporal sample size
##     y = NxTT multivariate count time series
##     W = NxN row-normalized weighted adjacency matrix describing the network
##     p = number of lags in the model
##     Z = Nxq matrix of covariates (one for each column), where q is the number of
##         covariates in the model.
##  output:
##     hh = (-1)*Hessian matrix
################################################################################

################################################################################
################################################################################
##  outer_logpq() --- function for the computation of information matrix of the
##  log-linear PNAR model with p lags and q covariates. Inputs:
##     b = parameters of the models
##     N = number of nodes on the network
##     TT = temporal sample size
##     y = NxTT multivariate count time series
##     W = NxN row-normalized weighted adjacency matrix describing the network
##     p = number of lags in the model
##     Z = Nxq matrix of covariates (one for each column), where q is the number of
##         covariates in the model.
##  output:
##     out = information matrix
################################################################################

.scor_hess_outer_logpq <- function(b, N, TT, y, wy, p, Z) {

  ## scor
  nut <- as.vector( wy %*% b)
  a <- wy * ( as.vector( y[, -c(1:p)] ) - exp(nut) )
  scor <- matrix( - Rfast::colsums(a), 2 * p + 1 + max( 0, ncol(Z) ), 1)

  ## hess
  nut <- as.vector( wy %*% b )
  hh <- crossprod( wy * exp(nut), wy )

  ## outer
  #out_log <- 0
  k <- rep( 1:c(TT - p), each = N )
  b1 <- rowsum(a, k)
  #for ( i in 1: c(TT - p) )  out_log <- out_log + tcrossprod(b1[i, ])
  out_log <- crossprod(b1)

  list(scor = scor, hh = hh, out = out_log)
}

################################################################################
################################################################################
##  log_lin_estimnarpq() --- function for the constrained estimation of
##  log-linear PNAR model with p lags and q covariates. Inputs:
##    x0 = starting value of the optimization
##    N = number of nodes on the network
##    TT = temporal sample size
##    y = NxTT multivariate count time series
##    W = NxN row-normalized weighted adjacency matrix describing the network
##    p = number of lags in the model
##    Z = Nxq matrix of covariates (one for each column), where q is the number of
##         covariates in the model.
##  output:
##    coefs = estimated QMLE coefficients
##    selin = standard errors estimates
##    tlin = t test estimates
##    score = value of the score at the optimization point
##    aic_lin = Akaike information criterion (AIC)
##    bic_lin = Bayesian information criterion (BIC)
##    qic_lin = Quasi information criterion (QIC)
################################################################################

log_lin_estimnarpq <- function(y, W, p, Z = NULL, uncons = FALSE, init = NULL, xtol_rel = 1e-8, maxeval = 100) {

  if ( min(W) < 0 ) {
    stop('The adjacency matrix W contains negative values.')
  }

  if ( !is.null(Z) ) {
    Z <- model.matrix( ~., as.data.frame(Z) )
    Z <- Z[1:dim(y)[2], -1, drop = FALSE]
  }
  y <- t(y)
  W <- W / Rfast::rowsums(W)
  W[ is.na(W) ] <- 0

  dm <- dim(y)   ;    N <- dm[1]    ;    TT <- dm[2]
  ly <- log1p(y)
  z <- W %*% ly
  wy <- NULL
  for ( ti in (p + 1):TT )  wy <- rbind( wy, cbind(z[, (ti - 1):(ti - p)], ly[, (ti - 1):(ti - p)], Z) )
  wy <- cbind(1, wy)

  if( uncons ) {
    #s_qmle <- nloptr::nloptr(x0 = x0, eval_f = .logl_log_linpq, eval_grad_f = .scor_logpq,
    #                  opts = opts, N = N, TT = TT, y = y, W = W, wy = wy, p = p, Z = Z)
    colnames(wy) <- NULL
    yp <- as.vector( y[, -c(1:p)] )
    mod <- glm(yp~., data = as.data.frame(wy[, -1]), poisson)
    s_qmle <- list()
    s_qmle$solution <- mod$coefficients
    s_qmle$objective <-  - as.numeric( logLik(mod) ) - sum( lfactorial(yp) )
    m <- length(mod$coefficients)

  } else {

    if ( is.null(init) ) {
      XX <- crossprod(wy)
      Xy <- Rfast::eachcol.apply(wy, as.vector( ly[, -c(1:p)] ) )
      x0 <- solve(XX, Xy)
    } else  x0 <- init

    m <- length(x0)
    # algorithm and and relative tolerance
    opts <- list( "algorithm" = "NLOPT_LD_SLSQP", "xtol_rel" = xtol_rel, "maxeval" = maxeval )

    # Inequality constraints (parameters searched in the stationary region)
    # b are the parameters to be constrained
    constr <- function(b, N, TT, y, W, wy, p, Z) {
      con <- sum( abs( b[2:(2 * p + 1)] ) ) - 1
      return(con)
    }

    # Jacobian of constraints
    # b are the parameters to be constrained
    j_constr <- function(b, N, TT, y, W, wy, p, Z) {
      j_con <- b/abs(b)
      j_con[1] <- 0
      return(j_con)
    }

    s_qmle <- nloptr::nloptr(x0 = x0, eval_f = .logl_log_linpq, eval_grad_f = .scor_logpq,
                      eval_g_ineq = constr, eval_jac_g_ineq = j_constr, opts = opts,
                      N = N, TT = TT, y = y, W = W, wy = wy, p = p, Z = Z)
  }
  coeflin <- s_qmle$solution

  ola <- .scor_hess_outer_logpq(coeflin, N, TT, y, wy, p, Z)
  S_lins <- ola$scor
  H_lins <- ola$hh
  G_lins <- ola$out

  solveH_lins <- solve(H_lins)
  V_lins <- solveH_lins %*% G_lins %*% solveH_lins
  SE_lins <- sqrt( diag(V_lins) )

  tlin <- coeflin/SE_lins
  pval <- 2 * pnorm(abs(tlin), lower.tail = FALSE)

  loglik <-  - s_qmle$objective
  aic_lins <- 2 * m + 2 * s_qmle$objective
  bic_lins <- log(TT) * m + 2 * s_qmle$objective
  qic_lins <- 2 * sum(H_lins * V_lins ) + 2 * s_qmle$objective

  coeflog <- cbind(coeflin, SE_lins, tlin, pval)
  colnames(coeflog) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  coeflog <- as.data.frame(coeflog)
  a <- pval
  a[ which(pval > 0.1) ]   <- c("   ")
  a[ which(pval < 0.1) ]   <- c(".  ")
  a[ which(pval < 0.05) ]  <- c("*  ")
  a[ which(pval < 0.01) ]  <- c("** ")
  a[ which(pval < 0.001) ] <- c("***")
  coeflog <- cbind(coeflog, a)
  colnames(coeflog)[5] <- ""

  if ( !is.null( dim(Z) ) ) {
    rownames(coeflog) <- rownames(S_lins) <- c( "beta0", paste("beta1", 1:p, sep =""), paste("beta2", 1:p, sep =""),
                                                paste("delta", 1:dim(Z)[2], sep ="") )
  } else  rownames(coeflog) <- rownames(S_lins) <- c( "beta0", paste("beta1", 1:p, sep =""), paste("beta2", 1:p, sep ="") )
  ic <- c(aic_lins, bic_lins, qic_lins)
  names(ic) <- c("AIC", "BIC", "QIC")

 # if ( !uncons ) {
 #   if ( any( abs(S_lins) > 1e-3 ) )  {
 #     warning( paste("Optimization failed in the stationary region. Please try estimation without stationarity constraints.") )
 #   }
 # }

 # if ( uncons ) {
 #   if ( any( abs(S_lins) > 1e-3 ) )  {
 #     warning( paste("The score function is not close to zero.") )
 #   }
 # }

  result <- list( coefs = coeflog, score = S_lins, loglik = loglik, ic = ic )
  class(result) <- "PNAR"
  return(result)
}


