poisson.MOD2.tnar <- function(d, B, W, gama, a, TT, N, copula = "gaussian", corrtype = "equicorrelation", rho, dof = 1) {

  if ( min(W) < 0 ) {
    stop('The adjacency matrix W contains negative values.')
  }

  W <- W / Rfast::rowsums(W)
  W[ is.na(W) ] <- 0
  
  y <- matrix(NA, N, TT)
  ca <- matrix(NA, N, TT)
  lambda <- matrix(NA, N, TT)
  lambda[, 1] <- 1
  unon <- matrix(1, N, 1)
  zeron <- matrix(0, N, 1)

  if ( copula != "clayton" )  {
    if ( corrtype == "equicorrelation" ) {
	  p2R <- NULL
      R <- matrix(rho, nrow = N, ncol = N)
      diag(R) <- 1
    } else {
      p2R <- rho^( 1:(N - 1) )
      R <- toeplitz( c(1, p2R ) )
    }
    cholR <- chol(R)
  } else  cholR <- NULL

  u <- rcopula(n, N, rho, copula, corrtype, dof, cholR)
  x <- Rfast::eachrow( log(u), -lambda[, 1], oper = "/" )
  y[, 1] <- getN(x, tt = 1)
  x <- u <- NULL

  com <- W %*% y[, 1]
  I1a <- ( com <= gama )
  X1a <- cbind( unon, com, y[, 1], I1a, com * I1a, y[, 1] * I1a, zeron)
  ca[, 1] <- X1a[, 4:6] %*% a

  for (ti in 2:TT) {

    lambda[, ti] <- d + B %*% y[, ti - 1] + c[, ti - 1]

    u <- rcopula(n, N, rho, copula, corrtype, dof, cholR)
    x <- Rfast::eachrow( log(u), -lambda[, ti], oper = "/" )
    y[, ti] <- getN(x, tt = 1)
    x <- u <- NULL

    com <- W %*% y[, ti]
    It <- ( com <= gama )
    Xt <- cbind( unon, com, y[, ti], It, com * It, y[, ti] * It, zeron)
    ca[, ti] <- Xt[, 4:6] %*% a

  }

  list(p2R = p2R, lambda = lambda, y = y)
}


#########################################################################################
#########################################################################################
##  poisson.MOD2.tar() --- function to generate counts from Poisson non-linear
##  Threshold PNAR(1) model (TNAR) with  parameters:
##      d = intercept vector
##      B = matrix of coefficients
##      W = NxN row-normalized weighted adjacency matrix describing the network
##      a = vector of non-linear parameters
##      gamma = nuisance threshold parameter
##      TT = temporal sample size
##      N = number of nodes on the network
##      copula =  from which copula generate the data:
##                "gaussian":      Gaussian copula with Toeplitz correlation matrix
##                "gaussian_rho":  Gaussian copula with constant correlation matrix
##                "student":       Student's t copula with Toeplitz correlation matrix
##                 ... we can add more
##      rho = copula parameter
##      df = degrees of freedom for Student's t copula
##  output, a list of three values:
##      y = NxTT matrix of generated counts for N time series over TT
##      lambda = NxTT matrix of generated Poisson mean for N time series over TT
##      p2R = Toeplitz correlation matrix employed in the copula
#########################################################################################

# poisson.MOD2.tar_old <- function(d, B, W, gamma, a, TT, N, copula, rho, dof) {
#   p2R <- rho^( 1:(N - 1) )     ### creates column vectors of corr matrix R=(rho^|i-j|)_(i,j),
#   R <- toeplitz( c(1, p2R) )   ### this is a symmetric Toeplitz matrix
#
#   y <- matrix(NA, N, TT)
#   ca <- matrix(NA, N, TT)
#   lambda <- matrix(NA, N, TT)
#   lambda[, 1] <- 1
#   unon <- matrix(1, N, 1)
#   zeron <- matrix(0, N, 1)
#
#   if ( copula == "gaussian" )     u <- copula::rCopula(1000, normalCopula(param = p2R, dispstr = "toep", dim = N) )
#   if ( copula == "gaussian_rho" ) u <- copula::rCopula(100, normalCopula(rho, dim = N) )
#   if ( copula == "student" )      u <- copula::rCopula(100, tCopula(param = p2R, dispstr = "toep", df = dof, dim = N) )
#   if ( copula == "clayton" )      u <- copula::rCopula(100, claytonCopula(rho, dim = N) )
#
#   x <- Rfast::eachrow( log(u), -lambda[, 1], oper = "/" )
#   y[, 1] <- getN(x, tt = 1)
#   x <- u <- NULL
#
#   I1a <- ( W %*% y[, 1] <= gama )
#   X1a <- cbind( unon, W %*% y[, 1], y[, 1], I1a, W %*% y[, 1] * I1a, y[, 1] * I1a, zeron)
#   ca[, 1] <- X1a[, 4:6] %*% a
#
#   for (ti in 2:TT) {
#
#     lambda[, ti] <- d + B %*% y[, ti - 1] + c[, ti - 1]
#
#     if ( copula == "gaussian" )     u <- copula::rCopula(1000, normalCopula(param = p2R, dispstr = "toep", dim = N) )
#     if ( copula == "gaussian_rho" ) u <- copula::rCopula(100, normalCopula(rho, dim = N) )
#     if ( copula == "student" )      u <- copula::rCopula(100, tCopula(param = p2R, dispstr = "toep", df = dof, dim = N) )
#     if ( copula == "clayton" )      u <- copula::rCopula(100, claytonCopula(rho, dim = N) )
#
#     x <- Rfast::eachrow( log(u), -lambda[, ti], oper = "/" )
#     y[, ti] <- getN(x, tt = 1)
#     x <- u <- NULL
#
#     It <- ( W %*% y[, ti] <= gama )
#     Xt <- cbind( unon, W %*% y[, ti], y[, ti], It, W %*% y[, ti] * It, y[, ti] * It, zeron)
#     ca[, ti] <- Xt[, 4:6] %*% a
#
#   }
#
#   list(p2R = p2R, lambda = lambda, y = y)
# }
