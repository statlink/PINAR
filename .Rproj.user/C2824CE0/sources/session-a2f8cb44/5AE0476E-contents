log_lin_narpq_init <- function(y, W, p, Z = NULL) {

  if ( min(W) < 0 ) {
    stop('The adjacency matrix W contains negative values.')
  }

  if ( !is.null(Z) ) {
    Z <- model.matrix( ~., as.data.frame(Z) )
    Z <- Z[1:dim(y)[2], -1, drop = FALSE]
  }

  y <- t(y)
  W <- W / Rfast::rowsums(W)
  W[ is.na(W) ] <- 0

  dm <- dim(y)   ;    N <- dm[1]    ;    TT <- dm[2]
  ly <- log1p(y)
  z <- W %*% ly
  wy <- NULL
  for ( ti in (p + 1):TT )  wy <- rbind( wy, cbind(z[, (ti - 1):(ti - p)], ly[, (ti - 1):(ti - p)], Z) )
  wy <- cbind(1, wy)

  XX <- crossprod(wy)
  Xy <- Rfast::eachcol.apply(wy, as.vector( ly[, -c(1:p)] ) )
  x0 <- solve(XX, Xy)
  x0
}
