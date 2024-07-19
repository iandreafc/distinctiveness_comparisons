#setwd("/home/ubuntu/r")

rm(list=ls())
library(distinctiveness)  #To compute distinctiveness centrality
library(igraph)
library(ggplot2)
library(reshape2)
library(gridExtra)
set.seed(14)  #For replication
increment <- 0.05  #Increment for sweeping parameter
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

#Compute correlations by parameter
correl <- function(dat, lbound, ubound, increment, keep_cols) {
  cor_centrality <- data.frame(a = numeric(), d1_gamma = numeric(), d2_gamma = numeric(), d3_gamma = numeric(),
                          d4_gamma = numeric(), d5_gamma = numeric(), d1_beta = numeric(), d2_beta = numeric(),
                          d3_beta = numeric(), d4_beta = numeric(), d5_beta = numeric(), beta_gamma = numeric())
  for (a in seq(from = lbound, to = ubound, by = increment)) {
    a <- round(a,2)
    cor_centrality <- rbind(cor_centrality, data.frame(a = a,
                                D1_Gamma = cor(dat$D1[which(dat$a==a)], dat$gamma[which(dat$a==a)], method = "spearman"),
                                D2_Gamma = cor(dat$D2[which(dat$a==a)], dat$gamma[which(dat$a==a)], method = "spearman"),
                                D3_Gamma = cor(dat$D3[which(dat$a==a)], dat$gamma[which(dat$a==a)], method = "spearman"),
                                D4_Gamma = cor(dat$D4[which(dat$a==a)], dat$gamma[which(dat$a==a)], method = "spearman"),
                                D5_Gamma = cor(dat$D5[which(dat$a==a)], dat$gamma[which(dat$a==a)], method = "spearman"),
                                D1_Beta = cor(dat$D1[which(dat$a==a)], dat$beta[which(dat$a==a)], method = "spearman"),
                                D2_Beta = cor(dat$D2[which(dat$a==a)], dat$beta[which(dat$a==a)], method = "spearman"),
                                D3_Beta = cor(dat$D3[which(dat$a==a)], dat$beta[which(dat$a==a)], method = "spearman"),
                                D4_Beta = cor(dat$D4[which(dat$a==a)], dat$beta[which(dat$a==a)], method = "spearman"),
                                D5_Beta = cor(dat$D5[which(dat$a==a)], dat$beta[which(dat$a==a)], method = "spearman"),
                                Beta_Gamma = cor(dat$beta[which(dat$a==a)], dat$gamma[which(dat$a==a)], method = "spearman")
                                  ))
  }
  # Filter columns based on the list of column names to keep
  cor_centrality <- cor_centrality[, intersect(names(cor_centrality), keep_cols)]
  return(cor_centrality)
}

correlKE <- function(dat, lbound, ubound, increment, keep_cols) {
  cor_centrality <- data.frame(a = numeric(), d1_gamma = numeric(), d2_gamma = numeric(), d3_gamma = numeric(),
                          d4_gamma = numeric(), d5_gamma = numeric(), d1_beta = numeric(), d2_beta = numeric(),
                          d3_beta = numeric(), d4_beta = numeric(), d5_beta = numeric(), beta_gamma = numeric())
  for (a in seq(from = lbound, to = ubound, by = increment)) {
    a <- round(a,2)
    cor_centrality <- rbind(cor_centrality, data.frame(a = a,
                                D1_Gamma = cor(dat$D1[which(dat$a==a)], dat$gamma[which(dat$a==a)], method = "kendall"),
                                D2_Gamma = cor(dat$D2[which(dat$a==a)], dat$gamma[which(dat$a==a)], method = "kendall"),
                                D3_Gamma = cor(dat$D3[which(dat$a==a)], dat$gamma[which(dat$a==a)], method = "kendall"),
                                D4_Gamma = cor(dat$D4[which(dat$a==a)], dat$gamma[which(dat$a==a)], method = "kendall"),
                                D5_Gamma = cor(dat$D5[which(dat$a==a)], dat$gamma[which(dat$a==a)], method = "kendall"),
                                D1_Beta = cor(dat$D1[which(dat$a==a)], dat$beta[which(dat$a==a)], method = "kendall"),
                                D2_Beta = cor(dat$D2[which(dat$a==a)], dat$beta[which(dat$a==a)], method = "kendall"),
                                D3_Beta = cor(dat$D3[which(dat$a==a)], dat$beta[which(dat$a==a)], method = "kendall"),
                                D4_Beta = cor(dat$D4[which(dat$a==a)], dat$beta[which(dat$a==a)], method = "kendall"),
                                D5_Beta = cor(dat$D5[which(dat$a==a)], dat$beta[which(dat$a==a)], method = "kendall"),
                                Beta_Gamma = cor(dat$beta[which(dat$a==a)], dat$gamma[which(dat$a==a)], method = "kendall")
                                  ))
  }
  # Filter columns based on the list of column names to keep
  cor_centrality <- cor_centrality[, intersect(names(cor_centrality), keep_cols)]
  return(cor_centrality)
}


plotcor <- function(corr, savename, corrtype, title = NULL) {
  gamma_cols <- grep("_Gamma$", names(corr), value = TRUE)
  beta_cols <- grep("_Beta$", names(corr), value = TRUE)
  # Add beta_gamma to the list of beta_cols
  beta_cols <- c(beta_cols, "Beta_Gamma")

  # Melt the dataframe to long format for ggplot
  correlations_gamma_melt <- melt(corr, id.vars = "a", measure.vars = gamma_cols)
  correlations_beta_melt <- melt(corr, id.vars = "a", measure.vars = beta_cols)

  # Plot for variables ending with "_gamma"
  plot_gamma <- ggplot(correlations_gamma_melt, aes(x = a, y = value, color = variable, linetype = variable)) +
    geom_line(linewidth = 1) +  # Set line width
    scale_color_manual(values = rainbow(length(gamma_cols))) +  # Assign different colors
    scale_linetype_manual(values = rep(c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash"), length(beta_cols))) +  # Assign different line types
    labs(x = "Alpha", y = paste0(corrtype, " Correlation")) +  # Axis labels
    theme_gray() +  # Apply gray theme
    theme(
      aspect.ratio = 1,  # Make the plot square
      panel.grid.major = element_blank(),  # Remove major gridlines
      panel.grid.minor = element_blank(),  # Remove minor gridlines
      panel.border = element_blank(),  # Remove panel border
      legend.position = "bottom",  # Move legend below
      legend.direction = "horizontal"  # Make legend horizontal
    ) +  geom_vline(xintercept = 1, linetype = "solid", color = "black")  # Add vertical line at a value 1

  # Plot for variables ending with "_beta"
  plot_beta <- ggplot(correlations_beta_melt, aes(x = a, y = value, color = variable, linetype = variable)) +
    geom_line(linewidth = 1) +  # Set line width
    scale_color_manual(values = rainbow(length(beta_cols))) +  # Assign different colors
    scale_linetype_manual(values = rep(c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash"), length(beta_cols))) +  # Assign different line types
    labs(x = "Alpha", y = paste0(corrtype, " Correlation")) +  # Axis labels
    theme_gray() +  # Apply gray theme
    theme(
      aspect.ratio = 1,  # Make the plot square
      panel.grid.major = element_blank(),  # Remove major gridlines
      panel.grid.minor = element_blank(),  # Remove minor gridlines
      panel.border = element_blank(),  # Remove panel border
      legend.position = "bottom",  # Move legend below
      legend.direction = "horizontal"  # Make legend horizontal
    ) +  geom_vline(xintercept = 1, linetype = "solid", color = "black")  # Add vertical line at a value 1

  # Remove label "Variable" from the legend
  plot_gamma <- plot_gamma + labs(color = NULL, linetype = NULL)
  plot_beta <- plot_beta + labs(color = NULL, linetype = NULL)

  # Display plots side by side
  #grid.arrange(plot_gamma, plot_beta, ncol = 2, top = title)

  #Save
  pdf(paste0(savename, "_", corrtype, ".pdf"), width = 12, height = 6)
  grid.arrange(plot_gamma, plot_beta, ncol = 2, top = title)
  dev.off()
}



#### Florentine (unweighted)  ####
florentine <- matrix(c(0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,
                       0,0,0,0,0,1,1,0,1,0,0,0,0,0,0,
                       0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,
                       0,0,0,0,0,0,1,0,0,0,1,0,0,1,0,
                       0,0,1,0,0,0,0,0,0,0,1,0,0,1,0,
                       0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,
                       0,1,0,1,0,0,0,1,0,0,0,0,0,0,1,
                       0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,
                       1,1,1,0,0,0,0,0,0,0,0,1,1,0,1,
                       0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,
                       0,0,0,1,1,0,0,0,0,0,0,0,0,1,0,
                       0,0,0,0,0,0,0,0,1,0,0,0,0,1,1,
                       0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,
                       0,0,0,1,1,0,0,0,0,0,1,1,0,0,0,
                       0,0,0,0,0,0,1,0,1,0,0,1,0,0,0),
                     nrow = 15, ncol = 15,
                     dimnames = list(c("Acciaiuoli","Albizzi","Barbadori","Bischeri","Castellani","Ginori","Guadagni","Lamberteschi","Medici","Pazzi","Peruzzi","Ridolfi","Salviati","Strozzi","Tornabuoni"),
                                     c("Acciaiuoli","Albizzi","Barbadori","Bischeri","Castellani","Ginori","Guadagni","Lamberteschi","Medici","Pazzi","Peruzzi","Ridolfi","Salviati","Strozzi","Tornabuoni")))

#Compute centralities for each parameter value
print("Florentine")
dat <- data.frame(nodes = character(), a = numeric(), D1 = numeric(), D2 = numeric(), D3 = numeric(), D4 = numeric(), D5 = numeric(), beta = numeric(), gamma = numeric())
for (a in seq(from = lbound, to = ubound, by = increment)) {
  dat <- rbind(dat, centrality(florentine, a = a, normalize = FALSE))}

correlations <- correl(dat, lbound, ubound, increment,
                      c("a", "D1_Gamma", "D2_Gamma", "D3_Gamma", "D4_Gamma", "D5_Gamma",
                        "D1_Beta", "D2_Beta", "D3_Beta", "D4_Beta", "D5_Beta", "Beta_Gamma"))
plotcor (correlations, "Florentine", "Spearman", "Florentine Unweighted")

correlationsKE <- correlKE(dat, lbound, ubound, increment,
                      c("a", "D1_Gamma", "D2_Gamma", "D3_Gamma", "D4_Gamma", "D5_Gamma",
                        "D1_Beta", "D2_Beta", "D3_Beta", "D4_Beta", "D5_Beta", "Beta_Gamma"))
plotcor (correlationsKE, "Florentine", "Kendall", "Florentine Unweighted")

# Save the correlations to a CSV file
write.csv(correlations, file = "FlorentineSP.csv", row.names = FALSE)
write.csv(correlationsKE, file = "FlorentineKE.csv", row.names = FALSE)



#### Scale Free Networks ####
#Unweighted
random <- data.frame(network = numeric(), a = numeric(), D1_Gamma = numeric(), D1_Beta = numeric(),
                      D2_Gamma = numeric(), D2_Beta = numeric(), D3_Gamma = numeric(), D3_Beta = numeric(),
                      D4_Gamma = numeric(), D4_Beta = numeric(), D5_Gamma = numeric(), D5_Beta = numeric(), Gamma_Beta = numeric())

randomKE <- data.frame(network = numeric(), a = numeric(), D1_Gamma = numeric(), D1_Beta = numeric(),
                      D2_Gamma = numeric(), D2_Beta = numeric(), D3_Gamma = numeric(), D3_Beta = numeric(),
                      D4_Gamma = numeric(), D4_Beta = numeric(), D5_Gamma = numeric(), D5_Beta = numeric(), Gamma_Beta = numeric())

for (network in 1:200) {
  print(paste0("Scale-Free Network #", network))
  net <- sample_pa(n = 1000, m = 2, directed = FALSE)
  adj_net <- as_adjacency_matrix(net, sparse = FALSE)

  #Compute centralities for each parameter value
  dat <- data.frame(nodes = character(), a = numeric(), D1 = numeric(), D2 = numeric(), D3 = numeric(), D4 = numeric(), D5 = numeric(), beta = numeric(), gamma = numeric())
  for (a in seq(from = lbound, to = ubound, by = increment)) {
    dat <- rbind(dat, centrality(adj_net, a = a, normalize = FALSE))
  }

  #Compute correlations by parameter
  correlations <- correl(dat, lbound, ubound, increment,
                      c("a", "D2_Gamma", "D3_Gamma",  "D5_Gamma",
                        "D2_Beta", "D3_Beta", "D5_Beta", "Beta_Gamma"))
  correlations$network <- network
  correlations <- correlations[,c("network", "a",  "D2_Gamma", "D3_Gamma",  "D5_Gamma",
                        "D2_Beta", "D3_Beta", "D5_Beta", "Beta_Gamma")]
  random <- rbind(random, correlations)

  correlationsKE <- correlKE(dat, lbound, ubound, increment,
                      c("a", "D2_Gamma", "D3_Gamma",  "D5_Gamma",
                        "D2_Beta", "D3_Beta", "D5_Beta", "Beta_Gamma"))
  correlationsKE$network <- network
  correlationsKE <- correlationsKE[,c("network", "a",  "D2_Gamma", "D3_Gamma",  "D5_Gamma",
                        "D2_Beta", "D3_Beta", "D5_Beta", "Beta_Gamma")]
  randomKE <- rbind(randomKE, correlationsKE)
}

correlations <- aggregate(random[,c( "D2_Gamma", "D3_Gamma",  "D5_Gamma",
                                    "D2_Beta", "D3_Beta",  "D5_Beta", "Beta_Gamma")], by = list(a=random$a), FUN = mean)

plotcor (correlations, "ScaleFreeCorrUnweighted", "Spearman", "Random Scale-Free Networks - Unweighted")

# Save the correlations to a CSV file
write.csv(correlations, file = "ScaleFreeCorrUnweightedSP.csv", row.names = FALSE)

correlationsKE <- aggregate(randomKE[,c( "D2_Gamma", "D3_Gamma",  "D5_Gamma",
                                    "D2_Beta", "D3_Beta",  "D5_Beta", "Beta_Gamma")], by = list(a=randomKE$a), FUN = mean)

plotcor (correlationsKE, "ScaleFreeCorrUnweighted", "Kendall", "Random Scale-Free Networks - Unweighted")

# Save the correlations to a CSV file
write.csv(correlationsKE, file = "ScaleFreeCorrUnweightedKE.csv", row.names = FALSE)


#Weighted
random <- data.frame(network = numeric(), a = numeric(), D1_Gamma = numeric(), D1_Beta = numeric(),
                      D2_Gamma = numeric(), D2_Beta = numeric(), D3_Gamma = numeric(), D3_Beta = numeric(),
                      D4_Gamma = numeric(), D4_Beta = numeric(), D5_Gamma = numeric(), D5_Beta = numeric(), Gamma_Beta = numeric())

randomKE <- data.frame(network = numeric(), a = numeric(), D1_Gamma = numeric(), D1_Beta = numeric(),
                      D2_Gamma = numeric(), D2_Beta = numeric(), D3_Gamma = numeric(), D3_Beta = numeric(),
                      D4_Gamma = numeric(), D4_Beta = numeric(), D5_Gamma = numeric(), D5_Beta = numeric(), Gamma_Beta = numeric())

for (network in 1:200) {
  print(paste0("Scale-Free Network #", network))
  net <- sample_pa(n = 1000, m = 2, directed = FALSE)
  E(net)$weight <- sample(x = c(1:20), size = gsize(net), replace = TRUE)
  adj_net <- as_adjacency_matrix(net, sparse = FALSE, attr = "weight")

  #Compute centralities for each parameter value
  dat <- data.frame(nodes = character(), a = numeric(), D1 = numeric(), D2 = numeric(), D3 = numeric(), D4 = numeric(), D5 = numeric(), beta = numeric(), gamma = numeric())
  for (a in seq(from = lbound, to = ubound, by = increment)) {
    dat <- rbind(dat, centrality(adj_net, a = a, normalize = FALSE))
  }

  #Compute correlations by parameter
  correlations <- correl(dat, lbound, ubound, increment,
                      c("a", "D1_Gamma", "D3_Gamma",  "D4_Gamma",
                        "D1_Beta", "D3_Beta", "D4_Beta", "Beta_Gamma"))
  correlations$network <- network
  correlations <- correlations[,c("network", "a",  "D1_Gamma", "D3_Gamma",  "D4_Gamma",
                        "D1_Beta", "D3_Beta", "D4_Beta", "Beta_Gamma")]
  random <- rbind(random, correlations)

  correlationsKE <- correlKE(dat, lbound, ubound, increment,
                      c("a", "D1_Gamma", "D3_Gamma",  "D4_Gamma",
                        "D1_Beta", "D3_Beta", "D4_Beta", "Beta_Gamma"))
  correlationsKE$network <- network
  correlationsKE <- correlationsKE[,c("network", "a",  "D1_Gamma", "D3_Gamma",  "D4_Gamma",
                        "D1_Beta", "D3_Beta", "D4_Beta", "Beta_Gamma")]
  randomKE <- rbind(randomKE, correlationsKE)
}

correlations <- aggregate(random[,c( "D1_Gamma", "D3_Gamma",  "D4_Gamma",
                                    "D1_Beta", "D3_Beta",  "D4_Beta", "Beta_Gamma")], by = list(a=random$a), FUN = mean)

plotcor (correlations, "ScaleFreeCorrWeighted", "Spearman", "Random Scale-Free Networks - Weighted")

# Save the correlations to a CSV file
write.csv(correlations, file = "ScaleFreeCorrWeightedSP.csv", row.names = FALSE)

correlationsKE <- aggregate(randomKE[,c( "D1_Gamma", "D3_Gamma",  "D4_Gamma",
                                    "D1_Beta", "D3_Beta",  "D4_Beta", "Beta_Gamma")], by = list(a=randomKE$a), FUN = mean)

plotcor (correlationsKE, "ScaleFreeCorrWeighted", "Kendall", "Random Scale-Free Networks - Weighted")

# Save the correlations to a CSV file
write.csv(correlationsKE, file = "ScaleFreeCorrWeightedKE.csv", row.names = FALSE)



#### Small-World Networks ####
#Unweighted
random <- data.frame(network = numeric(), a = numeric(), D1_Gamma = numeric(), D1_Beta = numeric(),
                      D2_Gamma = numeric(), D2_Beta = numeric(), D3_Gamma = numeric(), D3_Beta = numeric(),
                      D4_Gamma = numeric(), D4_Beta = numeric(), D5_Gamma = numeric(), D5_Beta = numeric(), Gamma_Beta = numeric())

randomKE <- data.frame(network = numeric(), a = numeric(), D1_Gamma = numeric(), D1_Beta = numeric(),
                      D2_Gamma = numeric(), D2_Beta = numeric(), D3_Gamma = numeric(), D3_Beta = numeric(),
                      D4_Gamma = numeric(), D4_Beta = numeric(), D5_Gamma = numeric(), D5_Beta = numeric(), Gamma_Beta = numeric())

for (network in 1:200) {
  print(paste0("Small-World Network #", network))
  net <- sample_smallworld(dim = 1, size = 1000, nei = 2, p = 0.05)
  adj_net <- as_adjacency_matrix(net, sparse = FALSE)

  #Compute centralities for each parameter value
  dat <- data.frame(nodes = character(), a = numeric(), D1 = numeric(), D2 = numeric(), D3 = numeric(), D4 = numeric(), D5 = numeric(), beta = numeric(), gamma = numeric())
  for (a in seq(from = lbound, to = ubound, by = increment)) {
    dat <- rbind(dat, centrality(adj_net, a = a, normalize = FALSE))
  }

  #Compute correlations by parameter
  correlations <- correl(dat, lbound, ubound, increment,
                      c("a", "D2_Gamma", "D3_Gamma",  "D5_Gamma",
                        "D2_Beta", "D3_Beta", "D5_Beta", "Beta_Gamma"))
  correlations$network <- network
  correlations <- correlations[,c("network", "a",  "D2_Gamma", "D3_Gamma",  "D5_Gamma",
                        "D2_Beta", "D3_Beta", "D5_Beta", "Beta_Gamma")]
  random <- rbind(random, correlations)

  correlationsKE <- correlKE(dat, lbound, ubound, increment,
                      c("a", "D2_Gamma", "D3_Gamma",  "D5_Gamma",
                        "D2_Beta", "D3_Beta", "D5_Beta", "Beta_Gamma"))
  correlationsKE$network <- network
  correlationsKE <- correlationsKE[,c("network", "a",  "D2_Gamma", "D3_Gamma",  "D5_Gamma",
                        "D2_Beta", "D3_Beta", "D5_Beta", "Beta_Gamma")]
  randomKE <- rbind(randomKE, correlationsKE)
}

correlations <- aggregate(random[,c( "D2_Gamma", "D3_Gamma",  "D5_Gamma",
                                    "D2_Beta", "D3_Beta",  "D5_Beta", "Beta_Gamma")], by = list(a=random$a), FUN = mean)

plotcor (correlations, "SmallWorldCorrUnweighted", "Spearman", "Random Small-World Networks - Unweighted")

# Save the correlations to a CSV file
write.csv(correlations, file = "SmallWorldUnweightedSP.csv", row.names = FALSE)

correlationsKE <- aggregate(randomKE[,c( "D2_Gamma", "D3_Gamma",  "D5_Gamma",
                                    "D2_Beta", "D3_Beta",  "D5_Beta", "Beta_Gamma")], by = list(a=randomKE$a), FUN = mean)

plotcor (correlationsKE, "SmallWorldCorrUnweighted", "Kendall", "Random Small-World Networks - Unweighted")

# Save the correlations to a CSV file
write.csv(correlationsKE, file = "SmallWorldUnweightedKE.csv", row.names = FALSE)


#Weighted
set.seed(2)  #A different seed only used to avoid isolates and allow for Gamma calculation
random <- data.frame(network = numeric(), a = numeric(), D1_Gamma = numeric(), D1_Beta = numeric(),
                      D2_Gamma = numeric(), D2_Beta = numeric(), D3_Gamma = numeric(), D3_Beta = numeric(),
                      D4_Gamma = numeric(), D4_Beta = numeric(), D5_Gamma = numeric(), D5_Beta = numeric(), Gamma_Beta = numeric())

randomKE <- data.frame(network = numeric(), a = numeric(), D1_Gamma = numeric(), D1_Beta = numeric(),
                      D2_Gamma = numeric(), D2_Beta = numeric(), D3_Gamma = numeric(), D3_Beta = numeric(),
                      D4_Gamma = numeric(), D4_Beta = numeric(), D5_Gamma = numeric(), D5_Beta = numeric(), Gamma_Beta = numeric())

for (network in 1:200) {
  print(paste0("Small-World Network #", network))
  net <- sample_smallworld(dim = 1, size = 1000, nei = 2, p = 0.05)
  E(net)$weight <- sample(x = c(1:20), size = gsize(net), replace = TRUE)
  adj_net <- as_adjacency_matrix(net, sparse = FALSE, attr = "weight")

  #Compute centralities for each parameter value
  dat <- data.frame(nodes = character(), a = numeric(), D1 = numeric(), D2 = numeric(), D3 = numeric(), D4 = numeric(), D5 = numeric(), beta = numeric(), gamma = numeric())
  for (a in seq(from = lbound, to = ubound, by = increment)) {
    dat <- rbind(dat, centrality(adj_net, a = a, normalize = FALSE))
  }

  #Compute correlations by parameter
  correlations <- correl(dat, lbound, ubound, increment,
                      c("a", "D1_Gamma", "D3_Gamma",  "D4_Gamma",
                        "D1_Beta", "D3_Beta", "D4_Beta", "Beta_Gamma"))
  correlations$network <- network
  correlations <- correlations[,c("network", "a",  "D1_Gamma", "D3_Gamma",  "D4_Gamma",
                        "D1_Beta", "D3_Beta", "D4_Beta", "Beta_Gamma")]
  random <- rbind(random, correlations)

  correlationsKE <- correlKE(dat, lbound, ubound, increment,
                      c("a", "D1_Gamma", "D3_Gamma",  "D4_Gamma",
                        "D1_Beta", "D3_Beta", "D4_Beta", "Beta_Gamma"))
  correlationsKE$network <- network
  correlationsKE <- correlationsKE[,c("network", "a",  "D1_Gamma", "D3_Gamma",  "D4_Gamma",
                        "D1_Beta", "D3_Beta", "D4_Beta", "Beta_Gamma")]
  randomKE <- rbind(randomKE, correlationsKE)
}

correlations <- aggregate(random[,c( "D1_Gamma", "D3_Gamma",  "D4_Gamma",
                                    "D1_Beta", "D3_Beta",  "D4_Beta", "Beta_Gamma")], by = list(a=random$a), FUN = mean)

plotcor (correlations, "SmallWorldCorrWeighted", "Spearman", "Random Small-World Networks - Weighted")

# Save the correlations to a CSV file
write.csv(correlations, file = "SmallWorldWeightedSP.csv", row.names = FALSE)

correlationsKE <- aggregate(randomKE[,c( "D1_Gamma", "D3_Gamma",  "D4_Gamma",
                                    "D1_Beta", "D3_Beta",  "D4_Beta", "Beta_Gamma")], by = list(a=randomKE$a), FUN = mean)

plotcor (correlationsKE, "SmallWorldCorrWeighted", "Kendall", "Random Small-World Networks - Weighted")

# Save the correlations to a CSV file
write.csv(correlationsKE, file = "SmallWorldWeightedKE.csv", row.names = FALSE)