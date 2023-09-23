LM_gama_tpnarpq <- function(gama, b, N, TT, y, W, p, d, Z) {
  
  pp <- 1 + 2 * p
  dimz <- ( !is.null(Z) ) * NCOL(Z)
  m <- 2 * pp + max(0, dimz)
  b[ (m - pp + 1):m ] <- 0
  
  z <- W %*% y
  wy <- NULL
  for ( ti in (p + 1):TT ) {
    f <- ( z[, ti - d] <= gama )
    wy <- rbind( wy, cbind( z[, (ti - 1):(ti - p)], y[, (ti - 1):(ti - p)], Z, f, z[, (ti - 1):(ti - p)] * f, y[, (ti - 1):(ti - p)] * f ) )
  }
  wy <- cbind(1, wy)
  lambdat <- as.vector( wy[, 1:c(m - pp)] %*% b[ 1:c(m - pp) ] )
  
  com <- ( as.vector( y[, -c(1:p)] ) - lambdat ) / lambdat
  ct <- as.vector( y[, -c(1:p)] ) / lambdat^2
  k <- rep( 1:c(TT - p), each = N )

  ## scor
  a <- wy * com
  S <- Rfast::colsums(a)

  ## hess
  H <- crossprod(wy * ct, wy)

  ## out
  out1 <- 0
  b1 <- rowsum(a, k)
  B <- crossprod(b1)
  
  solveHmpp <- solve( H[ 1:(m - pp), 1:(m - pp) ] )
  Sigma <- B[ (m - pp + 1):m , (m - pp + 1):m ] -
    H[ (m - pp + 1):m, 1:(m - pp) ] %*% solveHmpp %*% B[ 1:(m - pp), (m - pp + 1):m ] -
    B[ (m - pp + 1):m, 1:(m - pp) ] %*% solveHmpp %*% H[ 1:(m - pp), (m - pp + 1):m ] +
    H[ (m - pp + 1):m, 1:(m - pp) ] %*% solveHmpp %*% B[ 1:(m - pp), 1:(m - pp) ] %*% solveHmpp %*% H[ 1:(m - pp), (m - pp + 1):m ]
  
  - as.numeric( crossprod( S[ (m - pp + 1):m ], solve( Sigma, S[ (m - pp + 1):m ] ) ) )
}


global_optimise_LM <- function(gama_L, gama_U, len, b, N, TT, y, W, p, d, Z, tol = 1e-9) {
  
  gami <- vector()
  supLMi <- vector()
  
  x <- seq( gama_L, gama_U, length = len )
  
  for ( i in 1:(len - 1) ) {
    int <- c( x[i], x[i + 1] )
    opt <- optimise(LM_gama_tpnarpq, c(gama_L, gama_U), b, N, TT, y, W, p, d, Z, tol = tol)
    gami[i] <- opt$minimum
    supLMi[i] <- opt$objective
  }
  
  supLM <- min(supLMi)
  gama <- gami[ which(supLMi == supLM )[1] ]
  list( gama = gama, supLM = -supLM, int = x) )
}



LM_gama_tpnarpq_j <- function(gama, b, p, zf, wy1, com, ct, k, solveHmpp, msn) {
  
  f <- ( zf <= gama )
  wy <- cbind( wy1, f, wy1[, 2:(p + 1) ] * f, wy1[, (p + 2):dim(wy1)[2] ] * f )
  
  ## scor
  a <- wy * com
  S <- Rfast::eachcol.apply(a, msn)

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




score_test_tpnarpq_j <- function(supLM, b, N, TT, y, W, p, d, Z = NULL, J = 499, 
                                gama_L = 0.01, gama_U = 5, tol = 1e-9, ncores = 1) {
  
 if ( min(W) < 0 ) {
    stop('The adjacency matrix W contains negative values.')
  }

  W <- W / Rfast::rowsums(W)
  W[ is.na(W) ] <- 0
  dm <- dim(y)    ;    N <- dm[1]    ;    TT <- dm[2]

  pp <- 1 + 2 * p
  dimz <- ( !is.null(Z) ) * NCOL(Z)
  m <- 2 * pp + max(0, dimz)

  supLMj <- gamaj <- pval <- numeric(J)
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

  msn <- Rfast::matrnorm(TT, J)
  rows <- rep( (p + 1):TT, each = N )
  msn <- msn[rows, ]
  k <- rep( 1:c(TT - p), each = N )

  for(j in 1:J){  
    opttj <- optimise( LM_gama_tpnarpq_j, c(gama_L, gama_U), b = b, p = p, zf = zf, wy1 = wy1, com = com, 
                       ct = ct, k = k, solveHmpp = solveHmpp, msn = msn[, j] )
    gamaj[j] <- opttj$minimum
    supLMj[j] <-  -opttj$objective
    pval[j] <- ( supLMj[j] >= supLM )
  }
  
  pJ <- sum(pval) / J
  cpJ <- ( sum(pval) + 1 ) / (J + 1)
  
  list( pJ = pJ, cpJ = cpJ, supLMj = supLMj, gamaj = gamaj )
  
}



