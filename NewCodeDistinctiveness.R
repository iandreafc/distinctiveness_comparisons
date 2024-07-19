#install.packages("distinctiveness")
library(distinctiveness)  #To compute distinctiveness centrality
library(igraph)

set.seed(14)  #For replication

#### Define functions ####

#Compute gamma centrality
gamma <- function(A, G = 0, normalize = FALSE) {
  O <- matrix(1, nrow = nrow(A))
  gamma <- A %*% ((A %*% O)^G)
  if (normalize) {gamma <- (nrow(A)/sum(gamma)) * gamma}
  return(gamma)
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
net <- sample_pa(n = 1000, m = 2, directed = FALSE)
net_adj <- as_adjacency_matrix(net, sparse = FALSE)
wnet <- net
E(wnet)$weight <- sample(x = c(1:20), size = gsize(wnet), replace = TRUE)
wnet_adj <- as_adjacency_matrix(wnet, sparse = FALSE, attr = "weight")


#Compare D1, D3, D4 on weighted network
d1_from_matrix <- distinctiveness_from_matrix_D1 (wnet_adj, alpha = 1.5)
d1 <- distinctiveness(wnet, measures = "D1", alpha = 1.5)
all(d1_from_matrix == d1["D1"])
d3_from_matrix <- distinctiveness_from_matrix_D3 (wnet_adj, alpha = 1.5)
d3 <- distinctiveness(wnet, measures = "D3", alpha = 1.5)
all(round(d3_from_matrix, 8) == round(d3["D3"], 8))
d4_from_matrix <- distinctiveness_from_matrix_D4 (wnet_adj, alpha = 1.5)
d4 <- distinctiveness(wnet, measures = "D4", alpha = 1.5)
all(round(d4_from_matrix, 8) == round(d4["D4"], 8))


#Compare D2, D3, D5 on UNweighted network
d2_from_matrix <- distinctiveness_from_matrix_D2 (net_adj, alpha = 1.5)
d2 <- distinctiveness(net, measures = "D2", alpha = 1.5)
all(d2_from_matrix == d2["D2"])
d3_from_matrix <- distinctiveness_from_matrix_D3 (net_adj, alpha = 1.5)
d3 <- distinctiveness(net, measures = "D3", alpha = 1.5)
all(round(d3_from_matrix, 8) == round(d3["D3"], 8))
d5_from_matrix <- distinctiveness_from_matrix_D5 (net_adj, alpha = 1.5)
d5 <- distinctiveness(net, measures = "D5", alpha = 1.5)
all(d5_from_matrix == d5["D5"])


#Check D1 and D2, D4 and D5 are the same on unweighted
d1_from_matrix <- distinctiveness_from_matrix_D1 (net_adj, alpha = 1.5)
d2_from_matrix <- distinctiveness_from_matrix_D2 (net_adj, alpha = 1.5)
all(d1_from_matrix == d2_from_matrix)
d1 <- distinctiveness(net, measures = "D1", alpha = 1.5)
d2 <- distinctiveness(net, measures = "D2", alpha = 1.5)
all(d1 == d2)

d1_from_matrix <- distinctiveness_from_matrix_D1 (net_adj, alpha = 1.5)  #unweighted
d2_from_matrix <- distinctiveness_from_matrix_D2 (wnet_adj, alpha = 1.5) #weightedd
all(d1_from_matrix == d2_from_matrix)
d1 <- distinctiveness(net, measures = "D1", alpha = 1.5)
d2 <- distinctiveness(wnet, measures = "D2", alpha = 1.5)
all(d1 == d2)


#Show that D4 and D5 are the same on unweighted networks only for alpha = 1
d4_from_matrix <- distinctiveness_from_matrix_D4 (net_adj, alpha = 1)
d5_from_matrix <- distinctiveness_from_matrix_D5 (net_adj, alpha = 1)
all(round(d4_from_matrix,8) == round(d5_from_matrix,8))
d5 <- distinctiveness(net, measures = "D5", alpha = 1)
d4 <- distinctiveness(net, measures = "D4", alpha = 1)
all(d4 == d5)

d4_from_matrix <- distinctiveness_from_matrix_D4 (net_adj, alpha = 1)  #unweighted
d5_from_matrix <- distinctiveness_from_matrix_D5 (wnet_adj, alpha = 1) #weightedd
all(round(d4_from_matrix,8) == round(d5_from_matrix,8))
d5 <- distinctiveness(wnet, measures = "D5", alpha = 1)
d4 <- distinctiveness(net, measures = "D4", alpha = 1)
all(d4 == d5)

d4_from_matrix <- distinctiveness_from_matrix_D4 (net_adj, alpha = 1.5)
d5_from_matrix <- distinctiveness_from_matrix_D5 (net_adj, alpha = 1.5)
all(round(d4_from_matrix,8) == round(d5_from_matrix,8)) #FALSE
d5 <- distinctiveness(net, measures = "D5", alpha = 1.5)
d4 <- distinctiveness(net, measures = "D4", alpha = 1.5)
all(d4 == d5) #FALSE

d4_from_matrix <- distinctiveness_from_matrix_D4 (net_adj, alpha = 1.5)  #unweighted
d5_from_matrix <- distinctiveness_from_matrix_D5 (wnet_adj, alpha = 1.5) #weightedd
all(round(d4_from_matrix,8) == round(d5_from_matrix,8)) #FALSE
d5 <- distinctiveness(wnet, measures = "D5", alpha = 1.5)
d4 <- distinctiveness(net, measures = "D4", alpha = 1.5)
all(d4 == d5) #FALSE


#Gamma D5 equivalence
alpha = 5
G = -alpha
#Unweighted TRUE
gamma_results <- gamma(net_adj, G = G, normalize = FALSE)
d5 <- distinctiveness(net, measures = "D5", alpha = alpha)
all(round(gamma_results, 8) == round(d5["D5"], 8))
#Weighted FALSE
gamma_results <- gamma(wnet_adj, G = G, normalize = FALSE)
d5 <- distinctiveness(wnet, measures = "D5", alpha = alpha)
all(round(gamma_results, 8) == round(d5["D5"], 8))