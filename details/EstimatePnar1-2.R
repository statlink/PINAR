
################################################################################
## Function for estimating linear Poisson NAR model, PNAR, with 1 or 2 lags
################################################################################

####Load required packages

library(moments)
library(numDeriv)
library(nloptr)
library(car)


################################################################################
################################################################################
##  logl_lin() --- function for the computation of log-likelihood of the linear
##  PNAR model with 1 lag. Inputs:
##     c = parameters of the models
##     N = number of nodes on the network
##     TT = temporal sample size
##     y = NxTT multivariate count time series
##     W = NxN row-normalized weighted adjacency matrix describing the network
##  output:
##     g = (-1)* independence quasi log-likelihood
################################################################################

logl_lin <- function(c, N, TT, y, W){
  
  unon <- as.matrix(rep(1, N))
  g <- 0
  
  for (t in 2:TT){
    Xt <- cbind(unon, W%*%y[,t-1], y[,t-1])
    lambdat <- Xt%*%c
    g <- g + t(y[,t])%*%log(lambdat)-t(unon)%*%lambdat
  }
  
  g <- -1 * g
  return(g)
}


################################################################################
################################################################################
##  scor_lin() --- function for the computation of score of the linear
##  PNAR model with 1 lag. Inputs:
##     c = parameters of the models
##     N = number of nodes on the network
##     TT = temporal sample size
##     y = NxTT multivariate count time series
##     W = NxN row-normalized weighted adjacency matrix describing the network
##  output:
##     ss = (-1)* vector of quasi score
################################################################################

scor_lin <- function(c, N, TT, y, W){
  
  unon <- as.matrix(rep(1, N))
  ss <- as.matrix(rep(0, 3))
  
  for(t in 2:TT){
    Xt <- cbind(unon, W%*%y[,t-1], y[,t-1])
    lambdat <- Xt%*%c
    Dt <- diag(1/as.vector(lambdat))
    ss <- ss + t(Xt)%*%Dt%*%(y[,t]-lambdat)
  }
  
  ss <- -1 * ss
  return(ss)
}


################################################################################
################################################################################
##  outer_lin() --- function for the computation of information matrix of the
##  linear PNAR model with 1 lag. Inputs:
##     c = parameters of the models
##     N = number of nodes on the network
##     TT = temporal sample size
##     y = NxTT multivariate count time series
##     W = NxN row-normalized weighted adjacency matrix describing the network
##  output:
##     out = information matrix
################################################################################

outer_lin <- function(c, N, TT, y, W){
  
  unon <- as.matrix(rep(1, N))
  ss <- rep(0, 3)
  out <- matrix(0, nrow=3, ncol=3)
  
  for (t in 2:TT){
    Xt <- cbind(unon, W%*%y[,t-1], y[,t-1])
    lambdat <- Xt%*%c
    Dt <- diag(1/as.vector(lambdat))
    ss <- t(Xt)%*%Dt%*%(y[,t]-lambdat)
    out <- out + ss%*%t(ss)
  }
  
  return(out)
}


################################################################################
################################################################################
##  hess_lin() --- function for the computation of Hessian matrix of the
##  linear PNAR model with 1 lag. Inputs:
##     c = parameters of the models
##     N = number of nodes on the network
##     TT = temporal sample size
##     y = NxTT multivariate count time series
##     W = NxN row-normalized weighted adjacency matrix describing the network
##  output:
##     hh = (-1)*Hessian matrix
################################################################################

hess_lin <- function(c, N, TT, y, W){
  
  unon <- as.matrix(rep(1, N))
  hh <- matrix(0, nrow=3, ncol=3)
  
  for(t in 2:TT){
    Xt <- cbind(unon, W%*%y[,t-1], y[,t-1])
    lambdat <- Xt%*%c
    Ct <- diag(y[,t])%*%diag(1/as.vector(lambdat^2))
    hh <- hh + t(Xt)%*%Ct%*%Xt
  }
  
  return(hh)
}


################################################################################
################################################################################
##  logl_lin2() --- function for the computation of log-likelihood of the linear
##  PNAR model with 2 lags. Inputs:
##     c = parameters of the models
##     N = number of nodes on the network
##     TT = temporal sample size
##     y = NxTT multivariate count time series
##     W = NxN row-normalized weighted adjacency matrix describing the network
##  output:
##     g = (-1)* independence quasi log-likelihood
################################################################################

logl_lin2 <- function(c, N, TT, y, W){
  
  unon <- as.matrix(rep(1, N))
  g <- 0
  
  for (t in 3:TT){
    Xt <- cbind(unon, W%*%y[,t-1], y[,t-1], W%*%y[,t-2], y[,t-2])
    lambdat <- Xt%*%c
    g <- g + t(y[,t])%*%log(lambdat)-t(unon)%*%lambdat
  }
  
  g <- -1 * g
  return(g)
}


################################################################################
################################################################################
##  scor_lin2() --- function for the computation of score of the linear
##  PNAR model with 2 lags. Inputs:
##     c = parameters of the models
##     N = number of nodes on the network
##     TT = temporal sample size
##     y = NxTT multivariate count time series
##     W = NxN row-normalized weighted adjacency matrix describing the network
##  output:
##     ss = (-1)* vector of quasi score
################################################################################

scor_lin2 <- function(c, N, TT, y, W){
  
  unon <- as.matrix(rep(1, N))
  ss <- as.matrix(rep(0, 5))
  
  for(t in 3:TT){
    Xt <- cbind(unon, W%*%y[,t-1], y[,t-1], W%*%y[,t-2], y[,t-2])
    lambdat <- Xt%*%c
    Dt <- diag(1/as.vector(lambdat))
    ss <- ss + t(Xt)%*%Dt%*%(y[,t]-lambdat)
  }
  
  ss <- -1 * ss
  return(ss)
}


################################################################################
################################################################################
##  outer_lin2() --- function for the computation of information matrix of the
##  linear PNAR model with 2 lags. Inputs:
##     c = parameters of the models
##     N = number of nodes on the network
##     TT = temporal sample size
##     y = NxTT multivariate count time series
##     W = NxN row-normalized weighted adjacency matrix describing the network
##  output:
##     out = information matrix
################################################################################

outer_lin2 <- function(c, N, TT, y, W){
  
  unon <- as.matrix(rep(1, N))
  ss <- rep(0, 5)
  out <- matrix(0, nrow=5, ncol=5)
  
  for (t in 3:TT){
    Xt <- cbind(unon, W%*%y[,t-1], y[,t-1], W%*%y[,t-2], y[,t-2])
    lambdat <- Xt%*%c
    Dt <- diag(1/as.vector(lambdat))
    ss <- t(Xt)%*%Dt%*%(y[,t]-lambdat)
    out <- out + ss%*%t(ss)
  }
  
  return(out)
}


################################################################################
################################################################################
##  hess_lin2() --- function for the computation of Hessian matrix of the
##  linear PNAR model with 2 lags. Inputs:
##     c = parameters of the models
##     N = number of nodes on the network
##     TT = temporal sample size
##     y = NxTT multivariate count time series
##     W = NxN row-normalized weighted adjacency matrix describing the network
##  output:
##     hh = (-1)*Hessian matrix
################################################################################

hess_lin2 <- function(c, N, TT, y, W){
  
  unon <- as.matrix(rep(1, N))
  hh <- matrix(0, nrow=5, ncol=5)
  
  for(t in 3:TT){
    Xt <- cbind(unon, W%*%y[,t-1], y[,t-1], W%*%y[,t-2], y[,t-2])
    lambdat <- Xt%*%c
    Ct <- diag(y[,t])%*%diag(1/as.vector(lambdat^2))
    hh <- hh + t(Xt)%*%Ct%*%Xt
  }
  
  return(hh)
}


################################################################################
################################################################################
##  ols.nar1() --- Compute OLS estimates of the parameters to be used 
##  as a starting value for the Quasi Maximum Likelihood optimization 
##  of linear PNAR model with 1 lag.
##  Inputs:
##    N = number of nodes on the network
##    TT = temporal sample size
##    y = NxTT multivariate count time series
##    W = NxN row-normalized weighted adjacency matrix describing the network
##  output:
##    theta = vector of OLS estimated parameters
################################################################################

ols.nar1 <- function(N, TT, y, W){
  
  XX <- matrix(0, 3, 3)
  Xy <- matrix(0, 3, 1)
  
  for (t in 2:TT){
    XXt <- cbind(rep(1,N), W%*%y[,t-1], y[,t-1])
    XX <- XX + t(XXt)%*%XXt
    Xy <- Xy + t(XXt)%*%y[,t]
  }
  
  theta <- solve(XX)%*%Xy
  return(theta)
  
}


################################################################################
################################################################################
##  ols.nar2() --- Compute OLS estimates of the parameters to be used 
##  as a starting value for the Quasi Maximum Likelihood optimization 
##  of linear PNAR model with 2 lags. 
##  Inputs:
##    N = number of nodes on the network
##    TT = temporal sample size
##    y = NxTT multivariate count time series
##    W = NxN row-normalized weighted adjacency matrix describing the network
##  output:
##    theta = vector of OLS estimated parameters
################################################################################

ols.nar2 <- function(N, TT, y, W){
  
  XX <- matrix(0, 5, 5)
  Xy <- matrix(0, 5, 1)
  
  for (t in 3:TT){
    XXt <- cbind(rep(1,N), W%*%y[,t-1], y[,t-1], W%*%y[,t-2], y[,t-2])
    XX <- XX + t(XXt)%*%XXt
    Xy <- Xy + t(XXt)%*%y[,t]
  }
  
  theta <- solve(XX)%*%Xy
  for(k in 1:length(theta)){
    if(theta[k]<0) theta[k]=0.001

  }
  return(theta)
  
}


################################################################################
################################################################################
##  lin_estimnar1() --- function for the constrained estimation of linear PNAR
##  model with 1 lag.
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

lin_estimnar1 <- function(x0, N, TT, y, W){
  
  # Lower and upper bounds (positivity constraints)
  lb <- c(0,0,0)
  ub <- c(Inf,Inf,Inf)
  
  # Inequality constraints (parameters searched in the stationary region)
  # c are the parameters to be constrained
  constr <- function (c, N, TT, y, W) {
    con <- c[2] + c[3] -1
    return (con)
  }
  
  # Jacobian of constraints
  # c are the parameters to be constrained
  j_constr <- function (c, N, TT, y, W) {
    j_con <- c(0,1,1)
    return (j_con)
  }
  
  # algorithm and relative tolerance
  opts <- list("algorithm"="NLOPT_LD_SLSQP", "xtol_rel"=1.0e-8)
  
  s_qmle1 <- nloptr(x0=x0, eval_f=logl_lin, eval_grad_f=scor_lin,
                    lb=lb, ub=ub, eval_g_ineq=constr, eval_jac_g_ineq=j_constr,
                    opts=opts, N=N, TT=TT, y=y, W=W)
  
  S_lins <- scor_lin(s_qmle1$solution, N, TT, y, W)
  
  H_lins <- hess_lin(s_qmle1$solution, N, TT, y, W)
  G_lins <- outer_lin(s_qmle1$solution, N, TT, y, W)
  V_lins <- solve(H_lins)%*%G_lins%*%solve(H_lins)
  SE_lins <- sqrt(diag(V_lins))
  
  coeflin <- s_qmle1$solution
  tlin <- coeflin/SE_lins
  
  aic_lins <- 2*3+2*s_qmle1$objective
  bic_lins <- log(ncol(y))*3+2*s_qmle1$objective
  qic_lins <- 2*sum(diag(H_lins%*%V_lins))+2*s_qmle1$objective
  
  lin_e <- list(coeflin=coeflin, selin=SE_lins, tlin=tlin, score=S_lins,
                aic_lin=aic_lins, bic_lin=bic_lins, qic_lin=qic_lins)
  
  return(lin_e)
}


################################################################################
################################################################################
##  lin_estimnar2() --- function for the constrained estimation of linear PNAR
##  model with 2 lags.
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

lin_estimnar2 <- function(x0, N, TT, y, W){
  
  # Lower and upper bounds (positivity constraints)
  lb <- rep(0,5)
  ub <- rep(Inf,5)
  
  # Inequality constraints (parameters searched in the stationary region)
  # c are the parameters to be constrained
  constr <- function (c, N, TT, y, W) {
    con <- c[2] + c[3] + c[4] + c[5] -1
    return (con)
  }
  
  # Jacobian of constraints
  # c are the parameters to be constrained
  j_constr <- function (c, N, TT, y, W) {
    j_con <- c(0,1,1,1,1)
    return (j_con)
  }
  
  # algorithm and relative tolerance
  opts <- list("algorithm"="NLOPT_LD_SLSQP", "xtol_rel"=1.0e-8)
  
  s_qmle1 <- nloptr(x0=x0, eval_f=logl_lin2, eval_grad_f=scor_lin2,
                    lb=lb, ub=ub, eval_g_ineq=constr, eval_jac_g_ineq=j_constr,
                    opts=opts, N=N, TT=TT, y=y, W=W)
  
  S_lins <- scor_lin2(s_qmle1$solution, N, TT, y, W)
  
  H_lins <- hess_lin2(s_qmle1$solution, N, TT, y, W)
  G_lins <- outer_lin2(s_qmle1$solution, N, TT, y, W)
  V_lins <- solve(H_lins)%*%G_lins%*%solve(H_lins)
  SE_lins <- sqrt(diag(V_lins))
  
  coeflin <- s_qmle1$solution
  tlin <- coeflin/SE_lins
  
  aic_lins <- 2*5+2*s_qmle1$objective
  bic_lins <- log(ncol(y))*5+2*s_qmle1$objective
  qic_lins <- 2*sum(diag(H_lins%*%V_lins))+2*s_qmle1$objective
  
  lin_e <- list(coeflin=coeflin, selin=SE_lins, tlin=tlin, score=S_lins,
                aic_lin=aic_lins, bic_lin=bic_lins, qic_lin=qic_lins)
  
  return(lin_e)
}
