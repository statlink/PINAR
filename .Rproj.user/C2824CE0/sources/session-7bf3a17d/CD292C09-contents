########################################################################################
########################################################################################
##  adja() --- function to generate  a network from the Stochastic Block Model (SBM).
##  To do so, for each couple of nodes it performs a Bernoulli trial with values 1
##  "draw an edge", 0 "otherwise". The probabilities of these trials are bigger 
##  if the two nodes are in the same block, lower otherwise, and they are specified 
##  based on the following arguments:
##     N = number of nodes on the network
##     K = number of blocks
##     alpha = network density
##  return values:
##     W = row-normalized weighted adjacency matrix describing the network
########################################################################################

adja <- function(N, K, alpha, directed = FALSE) {
  p_in <- alpha * N^(-0.3)            # probability of an edge between nodes in same block
  p_out <- alpha / N             # probability of an edge between nodes in different blocks
  pm <- matrix(p_out, nrow = K, ncol = K)             
  diag(pm) <- p_in            
  gra <- igraph::sample_sbm(n = N, pref.matrix = pm, block.sizes = rep(N/K, K), directed = FALSE)
  A <- igraph::as_adjacency_matrix(gra, sparse = FALSE)  #generate adjacency matrix
  w <- A / Rfast::rowsums(A)
  w[ is.na(w) ] <- 0
  w
}