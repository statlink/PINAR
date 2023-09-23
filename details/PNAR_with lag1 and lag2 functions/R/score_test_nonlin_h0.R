score_test_nonlin_h0 <- function(ca, y, W) {

  if ( min(W) < 0 ) {
    stop('The adjacency matrix W contains negative values.')
  }

  W <- W / Rfast::rowsums(W)
  W[ is.na(W) ] <- 0

  dm <- dim(y)    ;    N <- dm[1]    ;    TT <- dm[2]

  ca[4] <- 0
  com <- as.vector( W %*% y[, -TT] )
  x2 <- log1p(com)
  wy <- cbind( 1, com, as.vector( y[, -TT] ), -ca[1] * x2 )
  lambdat <- as.vector( wy %*% ca )
  Dt <- 1 / lambdat
  ss1 <- crossprod( wy * Dt, as.vector(y[, -1]) - lambdat )
  S <- ss1

  ct <- as.vector(y[, -1]) / lambdat^2
  hh1 <- crossprod(wy * ct, wy)
  ss2a <- crossprod(x2^2 * Dt, as.vector(y[, -1]) - lambdat ) * ca[1]
  hh1[1, 4] <- hh1[1, 4] - ss1[4]/ca[1]
  hh1[4, 1] <- hh1[1, 4]
  hh1[4, 4] <- hh1[4, 4] - ss2a
  H <- hh1

  ss1 <- (wy * Dt) * ( as.vector(y[, -1]) - lambdat )
  #out1 <- 0
  k <- rep( 1:c(TT - 1), each = N )
  b1 <- rowsum(ss1, k)
  #for (i in 1: c(TT - 1) )  out1 <- out1 + tcrossprod(b1[i, ])
  out1 <- crossprod(b1)
  B <- out1

  sh13 <- solve(H[1:3, 1:3])
  Sigma <- B[4, 4] - H[4, 1:3] %*% sh13 %*% B[1:3, 4] - B[4, 1:3] %*% sh13 %*% H[1:3, 4] +
           H[4, 1:3] %*% sh13 %*% B[1:3, 1:3] %*% sh13 %*% H[1:3, 4]

  stat <- as.numeric(S[4])^2 / as.numeric(Sigma)
  pvalue <- pchisq(stat, 1, lower.tail = FALSE)
  list(stat = stat, pvalue = pvalue)
}


# score_test_nonlin_h0_old <- function(ca, N, TT, y, W) {
#
#   unon <- matrix(1, N, 1)
#   ss_der <- matrix(0, 4, 1)
#   ss <- matrix(0, 4, 4)
#   hh <- matrix(0, 4, 4)
#   out <- matrix(0, 4, 4)
#   ca[4] <- 0
#
#   for (ti in 2:TT) {
#     com <- W %*% y[, ti - 1]
#     x2 <- log1p(com)
#     Xt <- cbind( 1, com, y[, ti - 1], -ca[1] * x2 )
#     lambdat <- Xt %*% ca
#     Dt <- 1 / as.vector(lambdat)
#     ss <- crossprod(Xt * Dt, y[, ti] - lambdat)
#     ss_der <- ss_der + ss
#     ## ss_der is the S matrix
#
#     Ct <- y[, ti] / as.vector(lambdat^2)
#     hes <- crossprod(Xt * Ct, Xt)
#     ss2 <- crossprod(x2^2 * Dt, y[, ti] - lambdat ) * ca[1]
#     hes[1, 4] <- hes[1, 4] - ss[4]/ca[1]
#     hes[4, 1] <- hes[1, 4]
#     hes[4, 4] <- hes[4, 4] - ss2
#     hh <- hh + hes
#     ## hh is the H matrix
#
#     out <- out + tcrossprod(ss)
#     ## out is the B matrix
#   }
#
#   S <- ss_der
#   H <- hh
#   B <- out
#   sh13 <- solve(H[1:3, 1:3])
#   Sigma <- B[4, 4] - H[4, 1:3] %*% sh13 %*% B[1:3, 4] - B[4, 1:3] %*% sh13 %*% H[1:3, 4] +
#            H[4, 1:3] %*% sh13 %*% B[1:3, 1:3] %*% sh13 %*% H[1:3, 4]
#   as.numeric(S[4])^2/as.numeric(Sigma)
#
# }
#



