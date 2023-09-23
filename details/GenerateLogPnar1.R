
################################################################################
## Function for generating log-linear Poisson NAR model
################################################################################

####Load required packages

library(igraph)
library(copula)
# library(moments)
# library(numDeriv)
# library(car)

set.seed(123)


########################################################################################
########################################################################################
##  adja() --- function to generate  a network from the Stochastic Block Model (SBM).
##  To do so, for each couple of nodes it performs a Bernoulli trial with values 1
##  "draw an edge", 0 "otherwise". The probabilities of these trials are bigger 
##  if the two nodes are in the same block, lower otherwise, and they are specified 
##  based on the following arguments:
##     N = number of nodes on the network
##     K = number of blocks
##     alpha = network density
##  return values:
##     W = row-normalized weighted adjacency matrix describing the network
########################################################################################

adja <- function(N, K, alpha){
  p_in=alpha*N^(-0.3)            # probability of an edge between nodes in same block
  p_out=alpha*N^(-1)             # probability of an edge between nodes in different blocks
  pm <- matrix(p_out, nrow=K, ncol=K)             
  diag(pm) <- p_in            
  
  # generate randomly a SBM network
  gra <- sample_sbm(n=N, pref.matrix=pm, block.sizes=rep(N/K, K), directed = TRUE)
  A <- as_adjacency_matrix(gra, sparse=FALSE)  #generate adjacency matrix
  
  D <- diag(rowSums(A)^(-1))
  D[which(D==Inf)] <- 0
  W <- D%*%A                 #generate row-normalized adjacency matrix
  
  return(W)
  
}


########################################################################################
########################################################################################
##  adja_gnp() --- function to generate a network from the Erdos-Renyi model (ER).
##  To do so, for each couple of nodes it performs a Bernoulli trial with values 1
##  "draw an edge", 0 "otherwise". Each trial has the same probability of having an edge
##  specified based on the following arguments:
##     N = number of nodes on the network
##     alpha = network density
##  return values:
##     W = row-normalized weighted adjacency matrix describing the network
########################################################################################

adja_gnp <- function(N, alpha){
  
  # generate randomly a ER network
  gra <- sample_gnp(N, p=alpha*N^(-0.3), directed = TRUE, loops = FALSE)
  A <- as_adjacency_matrix(gra, sparse=FALSE)  #generate adjacency matrix
  
  D <- diag(rowSums(A)^(-1))
  D[which(D==Inf)] <- 0
  W <- D%*%A                        #generate row-normalized adjacency matrix                  #generate row-normalized adjacency matrix
  
  return(W)
  
}


################################################################################
################################################################################
##  Gmatr() --- function to generate the matrix of coefficients of the 
##  linear PNAR model, called G matrix. Inputs:
##     b0 = intercept coefficient
##     b1 = network coefficient
##     b2 = autoregressive coefficient
##     W = row-normalized weighted adjacency matrix describing the network
##     N = number of nodes on the network
##  output:
##     matr = a list of two elements:
##            b0 = a vector with the intercept coefficient
##            G = the matrix of coefficients of the linear PNAR model
################################################################################

Gmatr <- function(b0, b1, b2, W, N){
  G <- b1*W+b2*diag(N)
  b0<-b0*rep(1, length=N)
  
  matr <- list(b0, G)
  
  return(matr)
}


################################################################################
################################################################################
##  getN() --- function to count the number of events before a specified time
##  arguments:
##    x = a vector of inter-event times (assumed to be positive)
##    tt = a positive time
##  return value:
##    the number of events before time tt (possibly 0)
################################################################################

getN = function(x,tt=1) {
  if(sum(x)<tt)
    return(NA)
  if(x[1]>tt)
    return(0)
  else
    return(max(which(cumsum(x)<tt)))
}


#############################################################################################
#############################################################################################
##  poisson.MOD2.log() --- function to generate counts from Poisson log-linear 
##  PNAR(1) model with parameters:
##      d = intercept vector
##      B = matrix of coefficients
##      Time = temporal sample size
##      N = number of nodes on the network
##      copula =  from which copula generate the data:
##                "gaussian":      Gaussian copula with Toeplitz correlation matrix
##                "gaussian_rho":  Gaussian copula with constant correlation matrix
##                "student":       Student's t copula with Toeplitz correlation matrix
##                 ... we can add more
##      rho = copula parameter
##      df = degrees of freedom for Student's t copula
##  output, a list of three values:
##      y = NxTime matrix of generated counts for N time series over Time
##      log_lambda = NxTime matrix of generated Poisson log(mean) for N time series over Time
##      p2R = Toeplitz correlation matrix employed in the copula
#############################################################################################

poisson.MOD2.log <- function(d, B, Time, N, copula, rho, df)
{
  
  p2R <- rho^(1:(N-1))      ### creates column vectors of corr matrix R=(rho^|i-j|)_(i,j),
  R=toeplitz(c(1,p2R))      ### this is a symmetric Toeplitz matrix
  
  y=matrix(0, nrow=N, ncol=Time)
  nu=matrix(0, nrow=N, ncol=Time)
  nu[,1]=rep(0, N)
  
  if(copula=="gaussian") ustart=rCopula(100, normalCopula(param=p2R, dispstr = "toep", dim=N))
  
  if(copula=="gaussian_rho") ustart=rCopula(100, normalCopula(rho, dim = N))  
  
  if(copula=="student") ustart=rCopula(100, tCopula(param=p2R, dispstr = "toep", df=df, dim=N))
  
  if(copula=="clayton") ustart=rCopula(100, claytonCopula(rho, dim = N))                                            
  
  xstart <- matrix(nrow=100, ncol=N) 
  xstart  =t(t(-log(ustart))/exp(nu[,1]))
  y[,1] <-  apply(xstart,2, getN, tt=1)
  
  for (t in 2:Time){
    
    nu[,t]=d+B%*%log(1+y[,t-1])
    
    if(copula=="gaussian") u=rCopula(100, normalCopula(param=p2R, dispstr = "toep", dim=N)) 
    
    if(copula=="gaussian_rho") u=rCopula(100, normalCopula(rho, dim = N))  
    
    if(copula=="student") u=rCopula(100, tCopula(param=p2R, dispstr = "toep", df=df, dim=N))
    
    if(copula=="clayton") u=rCopula(100, claytonCopula(rho, dim = N)) 
    
    x <- matrix(nrow=100, ncol=N)
    x =t(t(-log(u))/exp(nu[ ,t]))   
    y[,t] <- apply(x,2, getN,tt=1)
    
  }
  
  res <- list(log_lambda=nu, y=y, p2R=p2R) 
  return(res)
  
}
