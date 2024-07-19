#setwd("C:/Users/a/Desktop/distcompare")

rm(list=ls())
library(distinctiveness)  #To compute distinctiveness centrality
library(igraph)
library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)
set.seed(1)  #For replication
increment <- 0.5  #Increment for sweeping parameter
#Set the distinctiveness alpha parameter between 0.5 and 3
lbound <- 0.5  #Alpha upper bound
ubound <- 3 #Alpha lower bound

#### Define functions ####
#Compute beta centrality
beta <- function(A, B = 0, normalize = FALSE) {
  I <- diag(nrow(A))
  O <- matrix(1, nrow = nrow(A))
  beta <- solve(I - (B * A)) %*% A %*% O  #Compute
  if (normalize) {beta <- beta * sqrt(nrow(A)/sum((beta)^2))}  #Normalize
  return(beta)
}

#Compute gamma centrality
gamma <- function(A, G = 0, normalize = FALSE) {
  O <- matrix(1, nrow = nrow(A))
  gamma <- A %*% ((A %*% O)^G)
  if (normalize) {gamma <- (nrow(A)/sum(gamma)) * gamma}
  return(gamma)
}


# Min-max normalization function
min_max_normalize <- function(x) {
  return((x - min(x)) / (max(x) - min(x)))
}

# Function to calculate the Ruzicka index
ruzicka_index <- function(A, B) {
  min_sum <- sum(pmin(A, B))
  max_sum <- sum(pmax(A, B))
  return (min_sum / max_sum)
}


# Compute the Ruzicka index matrix
compute_ruzicka_matrix <- function(df) {
  columns <- colnames(df)
  n <- length(columns)
  ruzicka_matrix <- matrix(0, n, n, dimnames = list(columns, columns))

  for (i in 1:n) {
    for (j in i:n) {
      A <- min_max_normalize(df[[i]])
      B <- min_max_normalize(df[[j]])
      ruzicka_matrix[i, j] <- ruzicka_index(A, B)
      ruzicka_matrix[j, i] <- ruzicka_matrix[i, j]
    }
  }

  ruzicka_matrix
}

#Compute distinctiveness, beta, and gamma centrality from an adjacency matrix using harmonized parameters
#'a' is the alpha parameter of distinctiveness centrality
centrality <- function(A, a = 0, normalize = FALSE, mode = "undirected") {
  centrality <- suppressMessages(distinctiveness(G = graph_from_adjacency_matrix(A, mode = mode, weighted = TRUE),
                                                 alpha = a,
                                                 normalize = normalize))[,2:6]
  centrality <- cbind(centrality, beta = beta(A = A,
                                              B = ((2/(1 + exp(1)^a))-1) * (1/eigen(A)$values[1]),
                                              normalize = normalize))
  centrality <- cbind(centrality, gamma = gamma(A = A,
                                                G = a * -1,
                                                normalize = normalize))
  rownames(centrality) <- NULL
  centrality <- round(centrality, 5)
  centrality$node <- rownames(A)
  centrality$a <- round(a,2)
  return(centrality)
}


#### Scale Free Network ####
#Unweighted
net <- sample_pa(n = 1000, m = 2, directed = FALSE)
adj_net <- as_adjacency_matrix(net, sparse = FALSE)
for (a in c(1,2,3)) {

  df <- centrality(adj_net, a = a, normalize = FALSE)

  # Create density plots
  #plot_D1 <- ggplot(df, aes(x = D1)) + geom_density(size=1.1) + labs(title = "Density Plot of D1") + theme_minimal()
  plot_D2 <- ggplot(df, aes(x = D2)) + geom_density(size=1.1) + labs(title = "Density Plot of D2") + theme_minimal()
  plot_D3 <- ggplot(df, aes(x = D3)) + geom_density(size=1.1) + labs(title = "Density Plot of D3") + theme_minimal()
  #plot_D4 <- ggplot(df, aes(x = D4)) + geom_density(size=1.1) + labs(title = "Density Plot of D4") + theme_minimal()
  plot_D5 <- ggplot(df, aes(x = D5)) + geom_density(size=1.1) + labs(title = "Density Plot of D5") + theme_minimal()
  plot_beta <- ggplot(df, aes(x = beta)) + geom_density(size=1.1) + labs(title = "Density Plot of Beta") + theme_minimal()
  plot_gamma <- ggplot(df, aes(x = gamma)) + geom_density(size=1.1) + labs(title = "Density Plot of Gamma") + theme_minimal()

  #Compare with Ruzicka index
  ruzicka_matrix <- compute_ruzicka_matrix(df[c("D2","D3","D5","beta","gamma")])
  diag(ruzicka_matrix) <- NA

  # Convert the matrix to a long format dataframe for ggplot
  melted_ruzicka_matrix <- melt(ruzicka_matrix, na.rm = TRUE)
  # Create a heatmap
  heatmap_plot <- ggplot(melted_ruzicka_matrix, aes(Var1, Var2, fill = value, label = round(value, 2))) +
    geom_tile(color = "white") +  # Add white borders to tiles
    geom_text(color = "black", size = 6, fontface = "bold") +  # Add text annotations
    scale_fill_gradient(low = "white", high = "#00e1ff", guide = "none") +
    labs(title = "Ruzicka Index Heatmap", x = "", y = "", fill = "Ruzicka Index") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  # Arrange the plots in a grid layout (3 plots per row)
  combined_plot <- grid.arrange(plot_D2, plot_D3, plot_D5, plot_beta, plot_gamma, heatmap_plot, ncol = 2,
  top = textGrob(paste0("Density Plots - Unweighted Scale Free Network - Alpha = ", as.character(a)), gp=gpar(fontsize=20,font=3)))

  # Save the plot to a PDF file with a specific name based on the value of 'a'
  pdf_filename <- paste0("ScaleFree_UNWEI_alpha", as.character(a), ".pdf")
  ggsave(pdf_filename, plot = combined_plot, width = 10, height = 14)
  cat("Plot saved as", pdf_filename, "\n")
}


#Weighted
E(net)$weight <- sample(x = c(1:80), size = gsize(net), replace = TRUE)
adj_net <- as_adjacency_matrix(net, sparse = FALSE, attr = "weight")
for (a in c(1,2,3)) {

  df <- centrality(adj_net, a = a, normalize = FALSE)

  # Create density plots
  plot_D1 <- ggplot(df, aes(x = D1)) + geom_density(size=1.1) + labs(title = "Density Plot of D1") + theme_minimal()
  #plot_D2 <- ggplot(df, aes(x = D2)) + geom_density(size=1.1) + labs(title = "Density Plot of D2") + theme_minimal()
  plot_D3 <- ggplot(df, aes(x = D3)) + geom_density(size=1.1) + labs(title = "Density Plot of D3") + theme_minimal()
  plot_D4 <- ggplot(df, aes(x = D4)) + geom_density(size=1.1) + labs(title = "Density Plot of D4") + theme_minimal()
  #plot_D5 <- ggplot(df, aes(x = D5)) + geom_density(size=1.1) + labs(title = "Density Plot of D5") + theme_minimal()
  plot_beta <- ggplot(df, aes(x = beta)) + geom_density(size=1.1) + labs(title = "Density Plot of Beta") + theme_minimal()
  plot_gamma <- ggplot(df, aes(x = gamma)) + geom_density(size=1.1) + labs(title = "Density Plot of Gamma") + theme_minimal()

  #Compare with Ruzicka index
  ruzicka_matrix <- compute_ruzicka_matrix(df[c("D1","D3","D4","beta","gamma")])
  diag(ruzicka_matrix) <- NA

  # Convert the matrix to a long format dataframe for ggplot
  melted_ruzicka_matrix <- melt(ruzicka_matrix, na.rm = TRUE)
  # Create a heatmap
  heatmap_plot <- ggplot(melted_ruzicka_matrix, aes(Var1, Var2, fill = value, label = round(value, 2))) +
    geom_tile(color = "white") +  # Add white borders to tiles
    geom_text(color = "black", size = 6, fontface = "bold") +  # Add text annotations
    scale_fill_gradient(low = "white", high = "#00e1ff", guide = "none") +
    labs(title = "Ruzicka Index Heatmap", x = "", y = "", fill = "Ruzicka Index") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  # Arrange the plots in a grid layout (3 plots per row)
  combined_plot <- grid.arrange(plot_D1, plot_D3, plot_D4, plot_beta, plot_gamma, heatmap_plot, ncol = 2,
  top = textGrob(paste0("Density Plots - Weighted Scale Free Network - Alpha = ", as.character(a)), gp=gpar(fontsize=20,font=3)))

  # Save the plot to a PDF file with a specific name based on the value of 'a'
  pdf_filename <- paste0("ScaleFree_WEI_alpha", as.character(a), ".pdf")
  ggsave(pdf_filename, plot = combined_plot, width = 10, height = 14)
  cat("Plot saved as", pdf_filename, "\n")
}