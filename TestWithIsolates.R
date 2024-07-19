#install.packages("distinctiveness")
library(distinctiveness)  #To compute distinctiveness centrality
library(igraph)

set.seed(1)  #For replication

#### Define functions ####

#Compute gamma centrality
gamma <- function(A, G = 0, normalize = FALSE) {
  O <- matrix(1, nrow = nrow(A))
  gamma <- A %*% ((A %*% O)^G)
  if (normalize) {gamma <- (nrow(A)/sum(gamma)) * gamma}
  return(gamma)
}

beta <- function(A, B = 0, normalize = FALSE) {
  I <- diag(nrow(A))
  O <- matrix(1, nrow = nrow(A))
  beta <- solve(I - (B * A)) %*% A %*% O  #Compute
  if (normalize) {beta <- beta * sqrt(nrow(A)/sum((beta)^2))}  #Normalize
  return(beta)
}

#Compute Distinctiveness D1 and D2 from adjacency matrix
distinctiveness_from_matrix_D1 <- function(adj_matrix, alpha = 1.0) {
  N <- nrow(adj_matrix)
  degrees <- colSums(adj_matrix != 0)
  distinctiveness <- adj_matrix %*% log10((N - 1) / degrees^alpha)

  return(distinctiveness)
}


distinctiveness_from_matrix_D2 <- function(adj_matrix, alpha = 1.0) {
  N <- nrow(adj_matrix)
  degrees <- colSums(adj_matrix != 0)
  distinctiveness <- (adj_matrix != 0) %*% log10((N - 1) / degrees^alpha)

  return(distinctiveness)
}


distinctiveness_from_matrix_D3 <- function(adjacency_matrix, alpha=1.0) {
  numerator <- sum(adjacency_matrix[upper.tri(adjacency_matrix)])
  denominator <- rowSums(adjacency_matrix^alpha) - adjacency_matrix^alpha + 1
  distinctiveness <- colSums(adjacency_matrix * log10(numerator / denominator))

  return(matrix(distinctiveness, ncol = 1))
}


distinctiveness_from_matrix_D4 <- function(adjacency_matrix, alpha = 1.0) {
  numerator <- adjacency_matrix * (adjacency_matrix^alpha)
  denominator <- rowSums(adjacency_matrix^alpha)
  distinctiveness <- colSums(numerator / denominator)

  return(matrix(distinctiveness, ncol = 1))
}


distinctiveness_from_matrix_D5 <- function(adj_matrix, alpha = 1.0) {
  degrees <- colSums(adj_matrix != 0)
  distinctiveness <- (adj_matrix != 0) %*% (1 / degrees^alpha)

  return(distinctiveness)
}


#Compare results from R distinctiveness package with the new functions
#Generating networks
net <- sample_pa(n = 100, m = 2, directed = FALSE)
# Add isolated nodes
num_isolated_nodes <- 5  # Number of isolated nodes to add
isolated_nodes <- make_empty_graph(n = num_isolated_nodes, directed = FALSE)
# Combine the generated network with the isolated nodes
net <- disjoint_union(net, isolated_nodes)
net_adj <- as_adjacency_matrix(net, sparse = FALSE)
wnet <- net
E(wnet)$weight <- sample(x = c(1:20), size = gsize(wnet), replace = TRUE)
wnet_adj <- as_adjacency_matrix(wnet, sparse = FALSE, attr = "weight")


#Compute metrics with on a net with isolates
alpha = 1.5
G = -alpha

#Unweighted
B = ((2/(1 + exp(1)^alpha))-1) * (1/eigen(net_adj)$values[1])
d1_from_matrix <- distinctiveness_from_matrix_D1 (net_adj, alpha = alpha)
d2_from_matrix <- distinctiveness_from_matrix_D2 (net_adj, alpha = alpha)
d3_from_matrix <- distinctiveness_from_matrix_D3 (net_adj, alpha = alpha)
d4_from_matrix <- distinctiveness_from_matrix_D4 (net_adj, alpha = alpha)
d5_from_matrix <- distinctiveness_from_matrix_D5 (net_adj, alpha = alpha)

d1 <- distinctiveness(net, measures = "D1", alpha = alpha)
d2 <- distinctiveness(net, measures = "D2", alpha = alpha)
d3 <- distinctiveness(net, measures = "D3", alpha = alpha)
d4 <- distinctiveness(net, measures = "D4", alpha = alpha)
d5 <- distinctiveness(net, measures = "D5", alpha = alpha)

gamma_results <- gamma(net_adj, G)
beta_results <- beta(net_adj, B)

# Combine all results into a single matrix
results_matrix <- cbind(d1_from_matrix, d2_from_matrix, d3_from_matrix, d4_from_matrix, d5_from_matrix,
                        d1["D1"], d2["D2"], d3["D3"], d4["D4"], d5["D5"], gamma_results, beta_results)
colnames(results_matrix) <- c("D1_from_matrix", "D2_from_matrix", "D3_from_matrix", "D4_from_matrix", "D5_from_matrix",
                              "D1", "D2", "D3", "D4", "D5", "Gamma", "Beta")
print(results_matrix)


#Weighted
B = ((2/(1 + exp(1)^alpha))-1) * (1/eigen(wnet_adj)$values[1])
d1_from_matrix <- distinctiveness_from_matrix_D1 (wnet_adj, alpha = alpha)
d2_from_matrix <- distinctiveness_from_matrix_D2 (wnet_adj, alpha = alpha)
d3_from_matrix <- distinctiveness_from_matrix_D3 (wnet_adj, alpha = alpha)
d4_from_matrix <- distinctiveness_from_matrix_D4 (wnet_adj, alpha = alpha)
d5_from_matrix <- distinctiveness_from_matrix_D5 (wnet_adj, alpha = alpha)

d1 <- distinctiveness(wnet, measures = "D1", alpha = alpha)
d2 <- distinctiveness(wnet, measures = "D2", alpha = alpha)
d3 <- distinctiveness(wnet, measures = "D3", alpha = alpha)
d4 <- distinctiveness(wnet, measures = "D4", alpha = alpha)
d5 <- distinctiveness(wnet, measures = "D5", alpha = alpha)

gamma_results <- gamma(wnet_adj, G)
beta_results <- beta(wnet_adj, B)

# Combine all results into a single matrix
results_matrix <- cbind(d1_from_matrix, d2_from_matrix, d3_from_matrix, d4_from_matrix, d5_from_matrix,
                        d1["D1"], d2["D2"], d3["D3"], d4["D4"], d5["D5"], gamma_results, beta_results)
colnames(results_matrix) <- c("D1_from_matrix", "D2_from_matrix", "D3_from_matrix", "D4_from_matrix", "D5_from_matrix",
                              "D1", "D2", "D3", "D4", "D5", "Gamma", "Beta")
print(results_matrix)