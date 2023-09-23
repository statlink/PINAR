LM_gama_stnarpq_j <- function(gama, b, p, m, z2, wy1, com, ct, msn, k, solveHmp) {

  f <- exp(-gama * z2)
  wy <- cbind(wy1, wy1[, 2:(p + 1)] * f)

  ## scor
  a <- wy * com
  S <- Rfast::eachcol.apply(a, msn)

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


score_test_stnarpq_j <- function(supLM, b, y, W, p, d, Z = NULL, J = 499, gama_L = 0.01, gama_U = 5, tol = 1e-9, ncores = 1) {

  if ( min(W) < 0 ) {
    stop('The adjacency matrix W contains negative values.')
  }

  W <- W / Rfast::rowsums(W)
  W[ is.na(W) ] <- 0
  dm <- dim(y)    ;    N <- dm[1]    ;    TT <- dm[2]

  dimz <- ( !is.null(Z) ) * NCOL(Z)
  m <- 1 + 3 * p + max(0, dimz)
  supLMj <- gamaj <- pval <- numeric(J)

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

  msn <- Rfast::matrnorm(TT, J)
  rows <- rep( (p + 1):TT, each = N )
  msn <- msn[rows, ]
  k <- rep( 1:c(TT - p), each = N )

  if ( ncores == 1 ) {
    for ( j in 1:J ) {
      opttj <- optimise( LM_gama_stnarpq_j, c( gama_L, gama_U ), b = b, p = p, m = m, z2 = z2, wy1 = wy1, 
                         com = com, ct = ct, msn = msn[, j], k = k, solveHmp = solveHmp, tol = tol )
      gamaj[j] <- opttj$minimum
      supLMj[j] <-  -opttj$objective
      pval[j] <- ( supLMj[j] >= supLM )
    }

  } else {
    oop <- options(warn = -1)
    on.exit( options(oop) )
    requireNamespace("doParallel", quietly = TRUE, warn.conflicts = FALSE)
    cl <- parallel::makePSOCKcluster(ncores)
    doParallel::registerDoParallel(cl)
    mod <- foreach(j = 1:J, .combine = rbind, .export = "LM_gama_stnarpq_j", .packages = "Rfast" ) %dopar% {
      opttj <- optimise( LM_gama_stnarpq_j, c( gama_L, gama_U ), b = b, N = N, TT = TT, y = y, p = p, d = d, 
                         m = m, z2 = z2, wy1 = wy1, com = com, ct = ct, msn = msn[, j], k = k, 
                         solveHmp = solveHmp, tol = tol )
      gamaj <- opttj$minimum
      supLMj <-  -opttj$objective
      pval <- ( supLMj >= supLM )
      return( c(gamaj, supLMj, pval) )
    }
    parallel::stopCluster(cl)
    gamaj <- mod[, 1]
    supLMj <- mod[, 2]
    pval <- mod[, 3]
  }

  pJ <- sum(pval) / J
  cpJ <- ( sum(pval) + 1 ) / (J + 1)
  list( pJ = pJ, cpJ = cpJ, supLMj = supLMj, gamaj = gamaj )
}
