
################################################################################
## Function for estimating log-linear Poisson NAR model with p lags, PNAR(p),
## and q non time-varying covariates
################################################################################

####Load required packages

#library(moments)
#library(numDeriv)
library(nloptr)
#library(car)


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


logl_log_linpq <- function(b, N, TT, y, W, p, Z){
  
  g <- 0
  ly <- log1p(y)
  
  for (t in (p+1):TT){
    
    Xt <- cbind(1, W%*%ly[,(t-1):(t-p)], ly[,(t-1):(t-p)], Z)
    nut <- Xt%*%b
    g <- g + t(y[,t])%*%nut-sum(exp(nut))
  }
  
  g <- -1 * g 
  return(g)
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

scor_logpq <- function(b, N, TT, y, W, p, Z){
  
  m <- 1+2*p+max(0,ncol(Z))
  ss <- as.matrix(rep(0, m))
  ly <- log1p(y)
  
  for(t in (p+1):TT){
    Xt <- cbind(1, W%*%ly[,(t-1):(t-p)], ly[,(t-1):(t-p)], Z)
    nut <- Xt%*%b
    ss <- ss + t(Xt)%*%(y[,t]-exp(nut))
  }
  
  ss <- -1 * ss
  return(ss)
}


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

outer_logpq <- function(b, N, TT, y, W, p, Z){
  
  m <- 1+2*p+max(0,ncol(Z))
  ss <- rep(0, m)
  out <- matrix(0, nrow=m, ncol=m)
  ly <- log1p(y)
  
  for(t in (p+1):TT){
    Xt <- cbind(1, W%*%ly[,(t-1):(t-p)], ly[,(t-1):(t-p)], Z)
    nut <- Xt%*%b
    ss <- t(Xt)%*%(y[,t]-exp(nut))
    out <- out +ss%*%t(ss)
  }
  
  return(out)
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

hess_logpq <- function(b, N, TT, y, W, p, Z){
  
  m <- 1+2*p+max(0,ncol(Z))
  hh <- matrix(0, nrow=m, ncol=m)
  ly <- log1p(y)
  
  for(t in (p+1):TT){
    Xt <- cbind(1, W%*%ly[,(t-1):(t-p)], ly[,(t-1):(t-p)], Z)
    nut <- Xt%*%b
    Dt <- diag(as.vector(exp(nut)))
    hh <- hh + t(Xt)%*%Dt%*%Xt
  }
  
  return(hh)
}



################################################################################
################################################################################
##  ols.narpq.log() --- Compute OLS estimates of the parameters to be used 
##  as a starting value for the Quasi Maximum Likelihood optimization 
##  of log-linear PNAR model with p lags and q covariates. Inputs:
##    N = number of nodes on the network
##    TT = temporal sample size
##    y = NxTT multivariate count time series
##    W = NxN row-normalized weighted adjacency matrix describing the network
##    p = number of lags in the model
##    Z = Nxq matrix of covariates (one for each column), where q is the number of 
##         covariates in the model.
##  output:
##    theta = vector of OLS estimated parameters
################################################################################

ols.narpq.log <- function(N, TT, y, W, p, Z){
  
  m <- 1+2*p+max(0,ncol(Z))
  XX <- matrix(0, m, m)
  Xy <- matrix(0, m, 1)
  ly <- log1p(y)
  
  for (t in (p+1):TT){
    XXt <- cbind(1, W%*%ly[,(t-1):(t-p)], ly[,(t-1):(t-p)], Z)
    XX <- XX + t(XXt)%*%XXt
    Xy <- Xy + t(XXt)%*%ly[,t]
  }
  
  theta <- solve(XX)%*%Xy
  return(theta)
  
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
##    uncons = logical, if TRUE an unconstrained optimization without stationarity
##             constraints is performed (default is FALSE)
##  output:
##    coeflin = estimated QMLE coefficients
##    selin = standard errors estimates
##    tlin = t test estimates
##    score = value of the score at the optimization point
##    aic_lin = Akaike information criterion (AIC)
##    bic_lin = Bayesian information criterion (BIC)
##    qic_lin = Quasi information criterion (QIC)
################################################################################

log_lin_estimnarpq <- function(x0, N, TT, y, W, p, Z, uncons=FALSE){
  
  m <- length(x0)
  
  # algorithm and and relative tolerance
  opts <- list("algorithm"="NLOPT_LD_SLSQP", "xtol_rel"=1.0e-8)
  
  if(uncons == TRUE){
    s_qmle <- nloptr(x0=x0, eval_f=logl_log_linpq, eval_grad_f=scor_logpq,
                     opts=opts, N=N, TT=TT, y=y, W=W, p=p, Z=Z)
  } else{
     
    # Inequality constraints (parameters searched in the stationary region)
    # b are the parameters to be constrained
    constr <- function (b, N, TT, y, W, p, Z) {
      con <- sum(abs(b[2:(2*p+1)])) -1
      return (con)
    }
    
    # Jacobian of constraints
    # b are the parameters to be constrained
    j_constr <- function (b, N, TT, y, W, p, Z) {
      j_con <- b/abs(b)
      j_con[1] <- 0
      return (j_con)
    }
    
    s_qmle <- nloptr(x0=x0, eval_f=logl_log_linpq, eval_grad_f=scor_logpq,
                     eval_g_ineq=constr, eval_jac_g_ineq=j_constr,
                     opts=opts, N=N, TT=TT, y=y, W=W, p=p, Z=Z)
    
    }
  
  
  
  S_lins <- scor_logpq(s_qmle$solution, N, TT, y, W, p, Z)
  
  H_lins <- hess_logpq(s_qmle$solution, N, TT, y, W, p, Z)
  G_lins <- outer_logpq(s_qmle$solution, N, TT, y, W, p, Z)
  V_lins <- solve(H_lins)%*%G_lins%*%solve(H_lins)
  SE_lins <- sqrt(diag(V_lins))
  
  coeflin <- s_qmle$solution
  tlin <- coeflin/SE_lins
  
  aic_lins <- 2*m+2*s_qmle$objective
  bic_lins <- log(ncol(y))*m+2*s_qmle$objective
  qic_lins <- 2*sum(diag(H_lins%*%V_lins))+2*s_qmle$objective
  
  lin_e <- list(coeflog=coeflin, selog=SE_lins, tlog=tlin, score=S_lins,
                aic_log=aic_lins, bic_log=bic_lins, qic_log=qic_lins)
  
  return(lin_e)
}

y = PNAR::crime
W = PNAR::crime_W
x0 = ols.narpq.log(N=nrow(y), TT=ncol(y), y=y, W=W, p=1, Z=NULL)
est1 = log_lin_estimnarpq(x0, N=nrow(y), TT=ncol(y), y=y, W=W, p=1, Z=NULL, uncons=FALSE)
est1     # fail optimization

est1.un = log_lin_estimnarpq(x0, N=nrow(y), TT=ncol(y), y=y, W=W, p=1, Z=NULL, uncons=TRUE)
est1.un                               # successful optimization
est1.un$coeflog[2]+est1.un$coeflog[3] # > 1 no stationarity 
