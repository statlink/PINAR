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

Gmatr <- function(b0, b1, b2, W, N) {
  G <- b1 * W + b2 * diag(N)
  b0 <- rep(b0, N)
  list(b0 = b0, G = G)
}