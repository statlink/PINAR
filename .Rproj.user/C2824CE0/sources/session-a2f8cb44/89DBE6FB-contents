.LM_gama_tnarpq <- function(gama, b, p, zf, wy1, com, ct, k, m, pp, solveHmpp) {

  f <- ( zf <= gama )
  wy <- cbind( wy1, f, wy1[, 2:dim(wy1)[2] ] * f )

  ## scor
  a <- wy * com
  S <- Rfast::colsums(a)

  ## hess
  H <- crossprod(wy * ct, wy)

  ## out
  b1 <- rowsum(a, k)
  B <- crossprod(b1)

  Sigma <- B[ (m - pp + 1):m , (m - pp + 1):m ] -
    H[ (m - pp + 1):m, 1:(m - pp) ] %*% solveHmpp %*% B[ 1:(m - pp), (m - pp + 1):m ] -
    B[ (m - pp + 1):m, 1:(m - pp) ] %*% solveHmpp %*% H[ 1:(m - pp), (m - pp + 1):m ] +
    H[ (m - pp + 1):m, 1:(m - pp) ] %*% solveHmpp %*% B[ 1:(m - pp), 1:(m - pp) ] %*% solveHmpp %*% H[ 1:(m - pp), (m - pp + 1):m ]

  - as.numeric( crossprod( S[ (m - pp + 1):m ], solve( Sigma, S[ (m - pp + 1):m ] ) ) )
}



global_optimise_LM_tnarpq <- function(gama_L = NULL, gama_U = NULL, len = 10, b, y, W, p, d, Z = NULL, tol = 1e-9) {

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
  pp <- 1 + 2 * p
  dimz <- ( !is.null(Z) ) * NCOL(Z)
  m <- 2 * pp + max(0, dimz)

  z <- W %*% y
  wy1 <- NULL
  for ( ti in (p + 1):TT )  wy1 <- rbind( wy1, cbind( z[, (ti - 1):(ti - p)], y[, (ti - 1):(ti - p)], Z, z[, ti - d] ) )
  zf <- wy1[, dim(wy1)[2]]
  wy1 <- cbind( 1, wy1[, - dim(wy1)[2]] )
  lambdat <- as.vector( wy1[, 1:c(m - pp)] %*% b[ 1:c(m - pp) ] )
  com <- ( as.vector( y[, -c(1:p)] ) - lambdat ) / lambdat
  ct <- as.vector( y[, -c(1:p)] ) / lambdat^2
  ## hess
  H <- crossprod(wy1 * ct, wy1)
  solveHmpp <- solve(H)
  k <- rep( 1:c(TT - p), each = N )

  nullgama_L <- is.null(gama_L)
  gami <- supLMi <- numeric(len - 1)
  if ( is.null(gama_L)  &  is.null(gama_U) ) {
    qq2 <- Rfast2::rowQuantile( z, probs = c(0.2, 0.8) )
    g <- Rfast::colmeans(qq2)
    gama_L <- max(g[1], 0.01)   ;   gama_U <- g[2]
  } else  if ( is.null(gama_L)  &  !is.null(gama_U) ) {
    qq2 <- Rfast2::rowQuantile( z, probs = c(0.2, 0.8) )
    g <- Rfast::colmeans(qq2)
    gama_L <- max(g[1], 0.01)
  } else  if ( !is.null(gama_L)  &  is.null(gama_U) ) {
    qq2 <- Rfast2::rowQuantile( z, probs = c(0.2, 0.8) )
    g <- Rfast::colmeans(qq2)
    gama_U <- g[2]
  }
  x <- seq(gama_L, gama_U, length = len)

  for( i in 1:(len - 1) ) {

    opt <- try( optimise( .LM_gama_tnarpq, c( x[i], x[i + 1] ), b = b, p = p, zf = zf, wy1 = wy1,
                     com = com, ct = ct, k = k, m = m, pp = pp, solveHmpp = solveHmpp ), silent = TRUE )
    if ( identical(class(opt), "try-error") ) {
      gami[i] <- NA
      supLMi[i] <- Inf
    } else {
      gami[i] <- opt$minimum
      supLMi[i] <- opt$objective
    }
  }

  if ( !nullgama_L ) {
    a <- min(gami, na.rm = TRUE)
    warning( paste("The minimum of the gamma values considered is", a) )
  }

  supLM <- min(supLMi)
  gama <- gami[ which.min(supLMi) ]
  list(gama = gama, supLM = -supLM, int = x)
}
