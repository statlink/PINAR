poisson.MOD3 <- function(d, B, Time, N, copula, rho, dof) {
  
  R <- matrix(rho, nrow = N, ncol = N)
  diag(R) <- 1

  y <- matrix(0, nrow = N, ncol = Time)
  lambda <- matrix(0, nrow = N, ncol = Time)
  lambda[, 1] <- 1
  
  cholR <- chol(R)
  z <- Rfast::matrnorm(Time, N)
  z <- z %*% cholR 
  u <- pnorm(z)  
  
  x <- Rfast::eachrow( log(u), -lambda[, 1], oper = "/" )
  y[, 1] <- getN(x, tt = 1)
  x <- u <- NULL

  for ( t in 2:Time ) {
    
    lambda[, t] <- d + B %*% y[, t - 1]
    
    z <- Rfast::matrnorm(Time, N)
    z <- z %*% cholR 
    u <- pnorm(z, log.p = TRUE)
    
    x <- Rfast::eachrow( u, -lambda[, t], oper = "/" )
    y[, t] <- getN(x, tt = 1)
    x <- u <- NULL
  }
  
  list(p2R = p2R, lambda = lambda, y = y)   
}


benchmark(poisson.MOD2(d, B, Time, N, copula = "gaussian_rho", rho), poisson.MOD3(d, B, Time, N, rho = rho), times = 10 )