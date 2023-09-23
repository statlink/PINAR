.LM_gama_stnarpq <- function(gama, b, p, m, z2, wy1, com, ct, k, solveHmp) {

  f <- exp(-gama * z2)
  wy <- cbind(wy1, wy1[, 2:(p + 1)] * f)

  ## scor
  a <- wy * com
  S <- Rfast::colsums(a)

  ## hess
  H <- crossprod(wy * ct, wy)

  ## out
  b1 <- rowsum(a, k)
  B <- crossprod(b1)

  Sigma <- B[ (m - p + 1):m , (m - p + 1):m ] -
           H[ (m - p + 1):m, 1:(m - p) ] %*% solveHmp %*% B[ 1:(m - p), (m - p + 1):m ] -
           B[ (m - p + 1):m, 1:(m - p) ] %*% solveHmp %*% H[ 1:(m - p), (m - p + 1):m ] +
           H[ (m - p + 1):m, 1:(m - p) ] %*% solveHmp %*% B[ 1:(m - p), 1:(m - p) ] %*% solveHmp %*% H[ 1:(m - p), (m - p + 1):m ]

  - as.numeric( crossprod( S[ (m - p + 1):m ], solve( Sigma, S[ (m - p + 1):m ] ) ) )

}


global_optimise_LM_stnarpq <- function(gama_L = NULL, gama_U = NULL, len = 10, b, y, W, p, d, Z = NULL, tol = 1e-9) {

  if ( min(W) < 0 ) {
    stop('The adjacency matrix W contains negative values.')
  }

  if ( !is.null(Z) ) {
    if ( min(Z) < 0 ) {
      stop('The matrix of covariates Z contains negative values.')
    }
    Z <- model.matrix(~., as.data.frame(Z))
    Z <- Z[1:dim(y)[2], -1, drop = FALSE]
  }

  y <- t(y)
  W <- W / Rfast::rowsums(W)
  W[ is.na(W) ] <- 0
  dm <- dim(y)   ;    N <- dm[1]    ;    TT <- dm[2]

  dimz <- ( !is.null(Z) ) * NCOL(Z)
  m <- 1 + 3 * p + max(0, dimz)

  z <- W %*% y
  z2 <- z^2
  wy1 <- NULL
  for ( ti in (p + 1):TT )
    wy1 <- rbind( wy1, cbind(z[, (ti - 1):(ti - p)], y[, (ti - 1):(ti - p)], Z, z2[, ti - d] ) )
  z2 <- wy1[, dim(wy1)[2]]
  wy1 <- wy1[, -dim(wy1)[2]]
  wy1 <- cbind(1, wy1)
  lambdat <- as.vector( wy1[, 1:c(m - p)] %*% b[1:c(m - p)] )
  com <- ( as.vector( y[, -c(1:p)] ) - lambdat ) / lambdat
  ct <- as.vector( y[, -c(1:p)] ) / lambdat^2
  H <- crossprod(wy1 * ct, wy1)
  solveHmp <- solve(H)  ## solve( H[ 1:(m - p), 1:(m - p) ] )
  k <- rep( 1:c(TT - p), each = N )

  gami <- supLMi <- numeric(len)
  if ( is.null(gama_L) &  is.null(gama_U) ) {
    x <- mean(z)
    gama_L <-  -log(0.9) / x^2
    gama_U <-  -log(0.1) / x^2
  } else if ( is.null(gama_L) &  !is.null(gama_U) ) {
    x <- mean(z)
    gama_L <-  -log(0.9) / x^2
  } else if ( !is.null(gama_L) &  is.null(gama_U) ) {
    x <- mean(z)
    gama_U <-  -log(0.1) / x^2
  }
  x <- seq(gama_L, gama_U, length = len)

  for( i in 1:(len - 1) ) {
    opt <- optimise(.LM_gama_stnarpq, c( x[i], x[i + 1] ), b = b, p = p, m = m, z2 = z2, wy1 = wy1,
                         com = com, ct = ct, k = k, solveHmp = solveHmp, tol = tol)
    gami[i] <- opt$minimum
    supLMi[i] <- opt$objective
  }

  supLM <- min(supLMi)
  gama <- gami[ which.min(supLMi) ]
  list(gama = gama, supLM = -supLM, int = x)
}
