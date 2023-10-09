rcopula <- function(n, N, copula = "gaussian", corrtype = "equicorrelation", rho, dof = 1, cholR = NULL) {

  if ( copula != "clayton" ) {

    if ( !is.null(cholR) ) {

      if ( copula == "gaussian" ) {
        #Rfast::matrnorm(n, N, seed = seed)
        z <- matrix( rnorm(n * N), ncol = N )
        z <- z %*% cholR
        u <- pnorm(z)
      } else if ( copula == "t" ) {
        #z <- Rfast::matrnorm(n, N, seed = seed)
        z <- matrix( rnorm(n * N), ncol = N )
        w <- sqrt( dof / rchisq(n, dof) )
        z <- w * (z %*% cholR)
        u <- pt(z, dof)
      }

    } else {  ##  else cholR is not given

      if ( corrtype == "equicorrelation" ) {
        R <- matrix(rho, nrow = N, ncol = N)
        diag(R) <- 1
      } else {
        p2R <- rho^( 1:(N - 1) )
        R <- toeplitz( c(1, p2R ) )
      }

      if ( copula == "gaussian" ) {
        #z <- Rfast::rmvnorm(n, numeric(N), R, seed = seed)
        z <- matrix( rnorm(n * N), ncol = N) %*% chol(R)
        u <- pnorm(z)
      } else if ( copula == "t" ) {
        #z <- Rfast::rmvt(n, numeric(N), R, v = dof)
        z <- matrix( rnorm(n * N), ncol = N )
        w <- sqrt( dof / rchisq(n, dof))
        z <- w * ( z %*% chol(R) )
        u <- pt(z, dof)
      }  ##  end  if ( copula == "gaussian" ) {

    }  ##  end  if ( !is.null(cholR) )

  } else {  ## else copula is "clayton"
    vi <- rgamma(n, 1/rho, 1)
    z <- matrix( runif(n * N), ncol = N )
    z <- log(z) / ( - vi)
    u <- (1 + z)^(-1/rho)
  }

  u
}


