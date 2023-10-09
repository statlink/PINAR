adja_gnp <- function(N, alpha, directed = FALSE) {

  gra <- igraph::sample_gnp(N, p = alpha * N^(-0.3), directed = directed, loops = FALSE)
  A <- igraph::as_adjacency_matrix(gra, sparse = FALSE) 
  w <- A / Rfast::rowsums(A)
  w[ is.na(w) ] <- 0
  w

}